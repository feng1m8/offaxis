#include <omp.h>

#include "offaxis/parameter.hxx"

#include "envs.hxx"
#include "histogram.hxx"
#include "kyn.hxx"
#include "raytracing.hxx"
#include "sphere.hxx"

namespace offaxis
{
    static std::valarray<double> offaxline(const std::valarray<double> &energy, const std::valarray<double> &param, int nside)
    {
        using namespace parameter::offaxline;

        const Kyn kyn(envs::table.at(80), param[a_spin], param[Incl]);

        const Sphere &sphere(envs::sphere.try_emplace(nside, nside).first->second);

        Histogram histogram(energy);

        Ray ray(param[rlp], utils::deg2rad(param[thetalp]), utils::deg2rad(param[philp]), &param[vr], param[a_spin], param[Rin], param[Rout]);

#pragma omp declare reduction(+ : Histogram : omp_out += omp_in) initializer(omp_priv = omp_orig)
#pragma omp parallel for firstprivate(ray) reduction(+ : histogram)
        for (std::size_t pix = 0; pix < sphere.size; ++pix)
        {
            auto [pr, ptheta, pphi] = sphere[pix];

            if (ray.tracing(pr, ptheta, pphi) == Ray::Disk)
            {
                double glp = ray.redshift();
                auto [gobs, cosem, lensing] = kyn.interpolate(ray->radius, ray->phi);
                double iobs = gobs * gobs * std::pow(glp, param[gamma]) * redshift(ray->radius, param[a_spin], ray->lambda) * cosem * lensing;
                histogram.accumulate(gobs, iobs);
            }
        }

        return histogram.get() / sphere.size;
    }

    [[gnu::dllexport]] void offaxline(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        using namespace parameter::offaxline;
        if (parameter.size() < Nparam)
            throw std::out_of_range("RealArray index " + std::to_string(Nparam - 1) + " is out of bounds with size " + std::to_string(parameter.size()) + ".");

        if (parameter[vr] * parameter[vr] + parameter[vtheta] * parameter[vtheta] + parameter[vphi] * parameter[vphi] > 1.0)
            flux = std::valarray<double>(std::nan(""), energy.size() - 1);

        if (envs::table.count(80) == 0)
            envs::table.try_emplace(80, envs::kydir());

        omp_set_num_threads(envs::nthreads());

        std::valarray<double> engs((1.0 + parameter[zshift]) / parameter[lineE] * energy);

        std::valarray<double> param(parameter);
        if (parameter[Rin] < 0.0)
            param[Rin] = -parameter[Rin] * rms(parameter[a_spin]);

        flux = offaxline(engs, param, envs::nside());

        if (parameter[normtype] != -1.0)
            flux /= flux.sum();
    }
}
