#include <cfenv>
#include <omp.h>

#include "relxill/src/Xillspec.h"

#include "offaxis/parameter.hxx"
#include "offaxline/envs.hxx"
#include "offaxline/kbhtables.hxx"
#include "offaxline/memory.hxx"
#include "offaxline/raytracing.hxx"
#include "offaxline/sphere.hxx"
#include "offaxxill/doppler.hxx"
#include "offaxxill/emission.hxx"
#include "offaxxill/spectrum.hxx"

namespace offaxis
{
    namespace offaxxillver
    {
        template <bool switch_reflfrac_boost>
        static Emission offaxis(const std::vector<double> &parameter, std::tuple<const offaxis::Sphere *, const offaxis::KBHtables *, int> environment)
        {
            using namespace parameter::offaxconv;

            const Sphere &sphere(*std::get<0>(environment));
            const KBHinterp interp(std::get<1>(environment)->interp(parameter[a_spin], parameter[Incl]));
            int n_incl = std::get<2>(environment);

            Histogram histogram(sphere.size, n_incl);

            Ray ray(parameter[rlp], utils::deg2rad(parameter[thetalp]), utils::deg2rad(parameter[philp]), &parameter[vr], parameter[a_spin], parameter[Rin], parameter[Rout]);

            omp_set_num_threads(envs::nthreads());
#pragma omp declare reduction(+ : Histogram : omp_out += omp_in) initializer(omp_priv = omp_orig)
#pragma omp parallel for firstprivate(ray) reduction(+ : histogram) schedule(static, 1)
            for (std::size_t pix = 0; pix < sphere.size; ++pix)
            {
                auto palpha = sphere[pix];

                if (ray.tracing(palpha[0], palpha[1], palpha[2]) == Ray::Disk)
                {
                    if constexpr (not switch_reflfrac_boost)
                    {
                        ++histogram.to_disk;
                    }

                    auto [gobs, cosem, lensing] = interp(ray->radius, ray->phi);

                    if (gobs >= EMIN_RELXILL_CONV and gobs < EMAX_RELXILL_CONV)
                    {
                        double glp = ray.redshift();

                        double iobs = gobs * gobs *
                                      std::pow(glp, parameter[gamma]) *
                                      redshift(ray->radius, parameter[a_spin], 0.0) *
                                      cosem * lensing;

                        histogram.accumulate(gobs, iobs, glp, cosem);
                    }
                }
                else if (ray->radius > 1000.0)
                {
                    if constexpr (not switch_reflfrac_boost)
                    {
                        ++histogram.to_infinity;
                    }
                }
            }

            return histogram.get();
        }

        static std::valarray<double> reflection(const std::vector<double> &energy, const std::vector<double> &parameter, T_PrimSpec prim_type)
        {
            using namespace parameter::offaxxillCp;

            std::vector params(std::begin(parameter) + rlp, gamma + 1 + std::begin(parameter));
            if (parameter[Rin] < 0.0)
                params[Rin - rlp] = -parameter[Rin] * rms(parameter[a_spin]);

            static Memory corona(offaxis<false>);
            corona.max_size = envs::cache_size();
            long nside = envs::nside();
            std::filesystem::path kydir(envs::kydir());

            const Emission emission = corona(
                params,
                {&envs::sphere.try_emplace(nside, nside).first->second,
                 &envs::kbhtables.try_emplace(kydir, kydir).first->second,
                 relxill::n_incl(prim_type)});

            const relxill::Spectrum spectrum(energy, parameter, prim_type);

            std::valarray flux(spectrum.convolve(emission));
            flux *= std::fabs(parameter[refl_frac]);

            if (parameter[refl_frac] >= 0.0)
                flux += spectrum.primary(parameter);

            return flux;
        }

        template <T_PrimSpec prim_type>
        static std::valarray<double> offaxxillver(const std::valarray<double> &energy, const std::valarray<double> &parameter)
        {
            using namespace parameter::offaxxillCp;

            if (parameter.size() < Nparams)
            {
                throw std::out_of_range(utils::fotmat("RealArray index %zu is out of bounds with size %zu.", Nparams - 1, parameter.size()));
            }

            if (parameter[vr] * parameter[vr] + parameter[vtheta] * parameter[vtheta] + parameter[vphi] * parameter[vphi] >= 1.0)
            {
                std::fetestexcept(FE_INVALID);
                return std::valarray(std::nan(""), energy.size() - 1);
            }

            std::valarray engs((1.0 + parameter[zshift]) * energy);

            std::valarray<double> flux(reflection(
                {std::begin(engs), std::end(engs)},
                {std::begin(parameter), std::end(parameter)},
                prim_type));

            return flux;
        };
    }

    [[gnu::dllexport]] void offaxxillCp(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        flux = offaxxillver::offaxxillver<T_PrimSpec::Nthcomp>(energy, parameter);
    }

    [[gnu::dllexport]] void offaxxill(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        using namespace parameter;

        if (parameter.size() < offaxxill::Nparams)
        {
            throw std::out_of_range(utils::fotmat("RealArray index %zu is out of bounds with size %zu.", offaxxill::Nparams - 1, parameter.size()));
        }

        std::valarray<double> param(offaxxillCp::Nparams);
        param[std::slice(offaxxillCp::rlp, offaxxillCp::logxi + 1, 1)] = parameter[std::slice(offaxxill::rlp, offaxxill::logxi + 1, 1)];
        param[offaxxillCp::logN] = 15.0;
        param[offaxxillCp::Afe] = parameter[offaxxill::Afe];

        double theta2rad = utils::deg2rad(parameter[offaxxill::thetalp]);
        double coslp = std::cos(theta2rad);
        double sinlp = std::sin(theta2rad);
        double dop = doppler(coslp, sinlp, parameter[offaxxill::philp], &parameter[offaxxill::vr], parameter[offaxxill::Incl]);
        double ginf = dop * expnu(parameter[offaxxill::rlp], coslp, sinlp, parameter[offaxxill::a_spin]);

        param[offaxxillCp::Ecut] = parameter[offaxxill::Ecut] / ginf;
        param[offaxxillCp::refl_frac] = parameter[offaxxill::refl_frac];
        param[offaxxillCp::zshift] = parameter[offaxxill::zshift];
        param[offaxxillCp::switch_reflfrac_boost] = parameter[offaxxill::switch_reflfrac_boost];

        flux = offaxxillver::offaxxillver<T_PrimSpec::CutoffPl>(energy, param);
    }
}
