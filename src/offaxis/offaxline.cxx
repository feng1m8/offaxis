#include <memory>

#include <omp.h>

#include "offaxis/parameter.hxx"

#include "envs.hxx"
#include "histogram.hxx"
#include "kyn.hxx"
#include "raytracing.hxx"
#include "sphere.hxx"

namespace offaxis
{
    static std::valarray<double> offaxline(const std::valarray<double> &energy, const std::valarray<double> &parameter, int nside)
    {
        using namespace parameter::offaxline;

        const auto kyn(std::make_unique<const Kyn>(envs::table.at(80), parameter[a_spin], parameter[Incl]));

        const Sphere &sphere(envs::sphere.try_emplace(nside, nside).first->second);

        const auto histogram(std::make_unique<const Histogram>(energy));

        double theta2rad = utils::deg2rad(parameter[thetalp]);
        double coslp = std::cos(theta2rad);
        double sinlp = std::sin(theta2rad);
        double vlp[3] = {parameter[vr], parameter[vtheta], parameter[vphi]};
        double phi2rad = utils::deg2rad(parameter[philp]);

        Ray ray(parameter[a_spin], parameter[rlp], coslp, sinlp, vlp, parameter[Rin], parameter[Rout]);

#pragma omp parallel for firstprivate(ray)
        for (std::size_t pix = 0; pix < sphere.npix; ++pix)
        {
            auto [pr, ptheta, pphi] = sphere[pix];
            double glp = ray(pr, ptheta, pphi);

            if (glp > 0.0)
            {
                auto [gobs, cosem, lensing] = kyn->interpolate(ray->radius, ray->phi + phi2rad);
                double iobs = gobs * gobs * std::pow(glp, parameter[gamma]) * redshift(ray->radius, parameter[a_spin], ray->lambda) * cosem * lensing;
                histogram->accumulate(gobs, iobs);
            }
        }

        return histogram->get() / sphere.npix;
    }

    void offaxline(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        using namespace parameter::offaxline;

        if (envs::table.count(80) == 0)
            envs::table.try_emplace(80, envs::kydir());

        omp_set_num_threads(envs::nthreads());

        std::valarray<double> engs((1.0 + parameter[zshift]) / parameter[lineE] * energy);
        std::valarray<double> param(parameter);
        // param[lineE] = 1.0;
        // param[zshift] = 0.0;
        // param[normtype] = -1.0;
        if (parameter[Rin] < 0.0)
            param[Rin] = -parameter[Rin] * rms(parameter[a_spin]);

        flux = offaxline(engs, param, envs::nside());

        if (parameter[normtype] == 0.0)
            flux /= flux.sum();
    }
}
