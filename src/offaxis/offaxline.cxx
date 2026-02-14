#include <cfenv>
#include <cstdarg>

#include <omp.h>
#include <valarray>

#include "offaxis/parameter.hxx"

#include "envs.hxx"
#include "histogram.hxx"
#include "kyn.hxx"
#include "raytracing.hxx"
#include "sphere.hxx"

#include "reflection/cache.hxx"

namespace offaxis
{
    namespace utils
    {
        static std::string string_print_fotmatted(const char *format, ...)
        {
            va_list args;
            va_start(args, format);
            std::string buffer(std::vsnprintf(nullptr, 0, format, args), '\0');
            va_end(args);

            va_start(args, format);
            std::vsnprintf(buffer.data(), buffer.size() + 1, format, args);
            va_end(args);

            return buffer;
        }
    }

    static std::valarray<double> offaxline(const std::valarray<double> &energy, const std::valarray<double> &parameter, const std::tuple<long, std::filesystem::path> &environment)
    {
        using namespace parameter::offaxconv;

        auto &[nside, kydir] = environment;
        const KBHinterp kyn(envs::kbhtables.try_emplace(kydir, kydir).first->second.interp(parameter[a_spin], parameter[Incl]));
        const Sphere &sphere(envs::sphere.try_emplace(nside, nside).first->second);

        Histogram histogram(energy);

        Ray ray(parameter[rlp], utils::deg2rad(parameter[thetalp]), utils::deg2rad(parameter[philp]), &parameter[vr], parameter[a_spin], parameter[Rin], parameter[Rout]);

        omp_set_num_threads(envs::nthreads());
#pragma omp declare reduction(+ : Histogram : omp_out += omp_in) initializer(omp_priv = omp_orig)
#pragma omp parallel for firstprivate(ray) reduction(+ : histogram)
        for (std::size_t pix = 0; pix < sphere.size; ++pix)
        {
            auto palpha = sphere[pix];

            if (ray.tracing(palpha[0], palpha[1], palpha[2]) == Ray::Disk)
            {
                double glp = ray.redshift();
                auto [gobs, cosem, lensing] = kyn(ray->radius, ray->phi);
                double iobs = gobs * gobs * std::pow(glp, parameter[gamma]) * redshift(ray->radius, parameter[a_spin], 0.0) * cosem * lensing;
                histogram.accumulate(gobs, iobs);
            }
        }

        return histogram.get() / sphere.size;
    }

    [[gnu::dllexport]] void offaxline(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        using namespace parameter::offaxline;

        if (parameter.size() < Nparams)
        {
            throw std::out_of_range(utils::string_print_fotmatted("RealArray index %zu is out of bounds with size %zu.", Nparams - 1, parameter.size()));
        }

        if (parameter[vr] * parameter[vr] + parameter[vtheta] * parameter[vtheta] + parameter[vphi] * parameter[vphi] >= 1.0)
        {
            std::fetestexcept(FE_INVALID);
            flux = std::valarray<double>(std::nan(""), energy.size() - 1);
            return;
        }

        std::valarray<double> engs((1.0 + parameter[zshift]) / parameter[lineE] * energy);

        std::valarray<double> param(parameter[std::slice(rlp, gamma - rlp + 1, 1)]);
        if (parameter[Rin] < 0.0)
            param[Rin - rlp] = -parameter[Rin] * rms(parameter[a_spin]);

        // static Cache<std::valarray<double>, const std::valarray<double> &, const std::valarray<double> &, const std::tuple<long, std::filesystem::path> &> memo(offaxline, envs::cache_size());

        flux = offaxline(engs, param, {envs::nside(), envs::kydir()});

        if (parameter[normtype] != -1.0)
            flux /= flux.sum();
    }
}
