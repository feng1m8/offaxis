#include <cfenv>

#include <omp.h>

#include "offaxis/parameter.hxx"
#include "offaxline/envs.hxx"
#include "offaxline/histogram.hxx"
#include "offaxline/kbhtables.hxx"
#include "offaxline/memory.hxx"
#include "offaxline/raytracing.hxx"
#include "offaxline/sphere.hxx"

static std::valarray<double> offaxline(const std::vector<double> &energy, const std::vector<double> &parameter, std::tuple<const offaxis::Sphere *, const offaxis::KBHtables *> environment)
{
    using namespace offaxis;
    using namespace parameter::offaxconv;

    const Sphere &sphere(*std::get<0>(environment));
    const KBHinterp interp(std::get<1>(environment)->interp(parameter[a_spin], parameter[Incl]));

    Histogram histogram(energy);

    Ray ray(parameter[rlp], utils::deg2rad(parameter[thetalp]), utils::deg2rad(parameter[philp]), &parameter[vr], parameter[a_spin], parameter[Rin], parameter[Rout]);

    omp_set_num_threads(envs::nthreads());
#pragma omp declare reduction(+ : Histogram : omp_out += omp_in) initializer(omp_priv = omp_orig)
#pragma omp parallel for firstprivate(ray) reduction(+ : histogram) schedule(static, 1)
    for (std::size_t pix = 0; pix < sphere.size; ++pix)
    {
        auto palpha = sphere[pix];

        if (ray.tracing(palpha[0], palpha[1], palpha[2]) == Ray::Disk)
        {
            double glp = ray.redshift();
            auto [gobs, cosem, lensing] = interp(ray->radius, ray->phi);

            double iobs = gobs * gobs *
                          std::pow(glp, parameter[gamma]) *
                          redshift(ray->radius, parameter[a_spin], 0.0) *
                          cosem * lensing;

            histogram.accumulate(gobs, iobs);
        }
    }

    return histogram.get() / sphere.size;
}

namespace offaxis
{
    [[gnu::dllexport]] void offaxline(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        using namespace parameter::offaxline;

        if (parameter.size() < Nparams)
        {
            throw std::out_of_range(utils::fotmat("RealArray index %zu is out of bounds with size %zu.", Nparams - 1, parameter.size()));
        }

        if (parameter[vr] * parameter[vr] + parameter[vtheta] * parameter[vtheta] + parameter[vphi] * parameter[vphi] >= 1.0)
        {
            std::fetestexcept(FE_INVALID);
            flux = std::valarray(std::nan(""), energy.size() - 1);
            return;
        }

        std::valarray<double> engs((1.0 + parameter[zshift]) / parameter[lineE] * energy);

        std::vector<double> params(std::begin(parameter) + rlp, gamma + 1 + std::begin(parameter));
        if (parameter[Rin] < 0.0)
            params[Rin - rlp] = -parameter[Rin] * rms(parameter[a_spin]);

        static Memory memo(::offaxline);
        memo.max_size = envs::cache_size();
        long nside = envs::nside();
        std::filesystem::path kydir(envs::kydir());

        flux = memo(
            {std::begin(engs), std::end(engs)}, params,
            {&envs::sphere.try_emplace(nside, nside).first->second,
             &envs::kbhtables.try_emplace(kydir, kydir).first->second});

        if (parameter[normtype] != -1.0)
            flux /= flux.sum();
    }
}
