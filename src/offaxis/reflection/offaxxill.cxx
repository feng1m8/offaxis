#include <omp.h>

#include "offaxis/parameter.hxx"
#include "relxill/src/Relbase.h"

#include "offaxis/envs.hxx"
#include "offaxis/kyn.hxx"
#include "offaxis/raytracing.hxx"
#include "offaxis/sphere.hxx"

#include "cache.hxx"
#include "emission.hxx"
#include "primary.hxx"
#include "spectrum.hxx"

namespace offaxis
{
    namespace offaxxillver
    {
        static Emission offaxis(const std::vector<double> &param, long nside, int n_incl)
        {
            using namespace parameter::offaxconv;

            const Kyn kyn(envs::table.at(80), param[a_spin], param[Incl]);

            const Sphere &sphere(envs::sphere.try_emplace(nside, nside).first->second);

            Histogram histogram(n_incl, sphere.size);

            Ray ray(param[rlp], utils::deg2rad(param[thetalp]), utils::deg2rad(param[philp]), &param[vr], param[a_spin], param[Rin], param[Rout]);

#pragma omp declare reduction(+ : Histogram : omp_out += omp_in) initializer(omp_priv = omp_orig)
#pragma omp parallel for firstprivate(ray) reduction(+ : histogram)
            for (std::size_t pix = 0; pix < sphere.size; ++pix)
            {
                auto palpha = sphere[pix];

                if (ray.tracing(palpha[0], palpha[1], palpha[2]) == Ray::Disk)
                {
                    ++histogram.to_disk;
                    auto [gobs, cosem, lensing] = kyn.interpolate(ray->radius, ray->phi);

                    if (gobs >= EMIN_RELXILL_CONV and gobs < EMAX_RELXILL_CONV)
                    {
                        double glp = ray.redshift();
                        double iobs = gobs * gobs * std::pow(glp, param[gamma]) * redshift(ray->radius, param[a_spin], ray->lambda) * cosem * lensing;
                        histogram.accumulate(gobs, iobs, cosem, glp);
                    }
                }
                else if (ray->radius > 1000.0)
                    ++histogram.to_infinity;
            }

            return histogram.get();
        }

        static auto reflection(const std::valarray<double> &energy, const std::valarray<double> &parameter, T_PrimSpec prim_type)
        {
            using namespace parameter;

            std::vector<double> param(offaxconv::Nparam);
            std::copy(std::begin(parameter) + offaxxillCp::rlp, std::begin(parameter) + offaxxillCp::gamma + 1, param.begin() + offaxconv::rlp);
            if (parameter[offaxxillCp::Rin] < 0.0)
                param[offaxconv::Rin] = -parameter[offaxxillCp::Rin] * rms(parameter[offaxxillCp::a_spin]);

            static Cache corona(offaxis, envs::cache_size());
            const Emission emission(corona(param, envs::nside(), relxill::n_incl(prim_type)));

            const relxill::Spectrum spectrum(parameter, prim_type);

            std::valarray<double> flux(energy.size() - 1);
            for (std::size_t i = 0; i < emission.glp.size(); ++i)
            {
                relxill::Spectrum spec = spectrum * emission.glp[i];
                relxill::Spectrum::Spec xillver(spec.xillver(energy, emission.dist[i]));

                flux += spec.norm * spec.convolve(xillver, energy, emission.hist[i]);
            }

            return std::make_tuple(flux, emission.f_refl);
        }

        template <T_PrimSpec T>
        static void offaxxillver(const std::valarray<double> &energy, const std::valarray<double> &param, std::valarray<double> &flux)
        {
            using namespace parameter::offaxxillCp;
            if (param.size() < Nparam)
                throw std::out_of_range("RealArray index " + std::to_string(Nparam - 1) + " is out of bounds with size " + std::to_string(param.size()) + ".");

            if (param[vr] * param[vr] + param[vtheta] * param[vtheta] + param[vphi] * param[vphi] > 1.0)
                flux = std::valarray<double>(std::nan(""), energy.size() - 1);

            if (envs::table.count(80) == 0)
                envs::table.try_emplace(80, envs::kydir());

            omp_set_num_threads(envs::nthreads());

            std::valarray<double> engs((1.0 + param[zshift]) * energy);

            auto [cflux, f_refl] = reflection(engs, param, T);

            flux = cflux * std::fabs(param[refl_frac]);

            if (param[switch_reflfrac_boost] == 0.0)
                flux *= f_refl;

            if (param[refl_frac] >= 0.0)
                flux += primary(engs, param, T);
        }
    }

    [[gnu::dllexport]] void offaxxillCp(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        offaxxillver::offaxxillver<T_PrimSpec::Nthcomp>(energy, parameter, flux);
    }

    [[gnu::dllexport]] void offaxxill(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        using namespace parameter;
        if (parameter.size() < offaxxill::Nparam)
            throw std::out_of_range("RealArray index " + std::to_string(offaxxill::Nparam - 1) + " is out of bounds with size " + std::to_string(parameter.size()) + ".");

        std::valarray<double> param(offaxxillCp::Nparam);
        param[std::slice(offaxxillCp::rlp, offaxxillCp::logxi + 1, 1)] = parameter[std::slice(offaxxill::rlp, offaxxill::logxi + 1, 1)];
        param[offaxxillCp::logN] = 15.0;
        param[offaxxillCp::Afe] = parameter[offaxxill::Afe];

        double theta2rad = utils::deg2rad(parameter[offaxxill::thetalp]);
        double coslp = std::cos(theta2rad);
        double sinlp = std::sin(theta2rad);
        double dop = offaxxillver::doppler(coslp, sinlp, parameter[offaxxill::philp], &parameter[offaxxill::vr], parameter[offaxxill::Incl]);
        double ginf = dop * offaxxillver::redshiftinf(parameter[offaxxill::rlp], coslp, sinlp, parameter[offaxxill::a_spin]);

        param[offaxxillCp::Ecut] = parameter[offaxxill::Ecut] / ginf;
        param[offaxxillCp::refl_frac] = parameter[offaxxill::refl_frac];
        param[offaxxillCp::zshift] = parameter[offaxxill::zshift];
        param[offaxxillCp::switch_reflfrac_boost] = parameter[offaxxill::switch_reflfrac_boost];

        offaxxillver::offaxxillver<T_PrimSpec::CutoffPl>(energy, param, flux);
    }
}
