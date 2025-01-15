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

        Histogram histogram(energy);

        Ray ray(parameter[rlp], utils::deg2rad(parameter[thetalp]), utils::deg2rad(parameter[philp]), &parameter[vr], parameter[a_spin], parameter[Rin], parameter[Rout]);

#pragma omp declare reduction(+ : Histogram : omp_out += omp_in) initializer(omp_priv = omp_orig)
#pragma omp parallel for firstprivate(ray) reduction(+ : histogram)
        for (std::size_t pix = 0; pix < sphere.size; ++pix)
        {
            auto [pr, ptheta, pphi] = sphere[pix];

            if (ray.tracing(pr, ptheta, pphi) == Ray::Disk)
            {
                double glp = ray.redshift();
                auto [gobs, cosem, lensing] = kyn->interpolate(ray->radius, ray->phi);
                double iobs = gobs * gobs * std::pow(glp, parameter[gamma]) * redshift(ray->radius, parameter[a_spin], ray->lambda) * cosem * lensing;
                histogram.accumulate(gobs, iobs);
            }
        }

        return histogram.get() / sphere.size;
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
