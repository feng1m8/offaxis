#include <memory>

#include <omp.h>

#include "offaxis/parameter.hxx"
#include "relxill/src/Relbase.h"

#include "cache.hxx"
#include "offaxis/envs.hxx"
#include "offaxis/kyn.hxx"
#include "offaxis/raytracing.hxx"
#include "spectrum.hxx"
#include "offaxis/sphere.hxx"

#include "emission.hxx"
namespace offaxis
{
    namespace offaxxillver
    {
        Reflection offaxis(const std::vector<double> &parameter, int nside, int n_incl)
        static Emission offaxis(const std::vector<double> &param, int nside, int n_incl)
        {
            using namespace parameter::offaxconv;

            const auto kyn(std::make_unique<const Kyn>(envs::table.at(80), parameter[a_spin], parameter[Incl]));

            const Sphere &sphere(envs::sphere.try_emplace(nside, nside).first->second);

            Histogram histogram(n_incl, sphere.size);

            Ray ray(parameter[rlp], utils::deg2rad(parameter[thetalp]), utils::deg2rad(parameter[philp]), &parameter[vr], parameter[a_spin], parameter[Rin], parameter[Rout]);

#pragma omp declare reduction(+ : Histogram : omp_out += omp_in) initializer(omp_priv = omp_orig)
#pragma omp parallel for firstprivate(ray) reduction(+ : histogram)
            for (std::size_t pix = 0; pix < sphere.size; ++pix)
            {
                auto [pr, ptheta, pphi] = sphere[pix];

                if (ray.tracing(pr, ptheta, pphi) == Ray::Disk)
                {
                    ++histogram.to_disk;
                    auto [gobs, cosem, lensing] = kyn->interpolate(ray->radius, ray->phi);

                    if (gobs >= EMIN_RELXILL_CONV and gobs < EMAX_RELXILL_CONV)
                    {
                        double glp = ray.redshift();
                        double iobs = gobs * gobs * std::pow(glp, parameter[gamma]) * redshift(ray->radius, parameter[a_spin], ray->lambda) * cosem * lensing;
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
            // param[offaxconv::normtype] = -1.0;

            static Cache corona(offaxis, 2);
            const auto reflection(std::make_unique<const Emission>(corona(param, envs::nside(), relxill::n_incl(prim_type))));

            const relxill::Spectrum spectrum(parameter, prim_type);

            std::valarray<double> flux(energy.size() - 1);
            for (std::size_t i = 0; i < emission->glp.size(); ++i)
            {
                relxill::Spectrum spec = spectrum * emission->glp[i];
                relxill::Spectrum::Spec xillver(spec.xillver(energy, emission->dist[i]));

                flux += spec.norm * spec.convolve(xillver, energy, emission->hist[i]);
            }

            return std::make_tuple(flux, emission->f_refl);
        }

        static std::valarray<double> offaxxillver(const std::valarray<double> &energy, const std::valarray<double> &parameter, T_PrimSpec prim_type)
        {
            using namespace parameter::offaxxillCp;

            if (parameter[vr] * parameter[vr] + parameter[vtheta] * parameter[vtheta] + parameter[vphi] * parameter[vphi] > 1.0)
            {
                static int warned = 0;
                if (warned == 0)
                {
                    std::cerr << "RuntimeWarning: Cannot perform model calculation in user-defined Python model." << std::endl;
                    std::cerr << "RuntimeWarning: Velocity of corona cannot be greater than 1.0. The flux array will be filled with NAN." << std::endl;
                    ++warned;
                }
                return std::valarray<double>(std::nan(""), energy.size() - 1);
            }

            if (envs::table.count(80) == 0)
                envs::table.try_emplace(80, envs::kydir());

            omp_set_num_threads(envs::nthreads());

            std::valarray<double> engs((1.0 + parameter[zshift]) * energy);

            auto [flux, f_refl] = reflection(engs, parameter, prim_type);

            flux *= std::fabs(parameter[refl_frac]);

            if (parameter[switch_reflfrac_boost] == 0.0)
                flux *= f_refl;

            if (parameter[refl_frac] >= 0.0)
                flux += primary(engs, parameter, prim_type);

            return flux;
        }
    }

    void offaxxillCp(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        flux = offaxxillver::offaxxillver(energy, parameter, T_PrimSpec::Nthcomp);
    }

    void offaxxill(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        using namespace parameter;

        std::valarray<double> param(offaxxillCp::Nparam);
        param[std::slice(offaxxillCp::rlp, offaxxillCp::logxi + 1, 1)] = parameter[std::slice(offaxxill::rlp, offaxxill::logxi + 1, 1)];
        param[offaxxillCp::logN] = 15.0;
        param[offaxxillCp::Afe] = parameter[offaxxill::Afe];

        double theta2rad = utils::deg2rad(parameter[offaxxill::thetalp]);
        double coslp = std::cos(theta2rad);
        double sinlp = std::sin(theta2rad);
        double vlp[3] = {parameter[offaxxill::vr], parameter[offaxxill::vtheta], parameter[offaxxill::vphi]};
        double dop = offaxxillver::doppler(coslp, sinlp, parameter[offaxxill::philp], vlp, parameter[offaxxill::Incl]);
        double ginf = dop * offaxxillver::redshiftinf(parameter[offaxxill::rlp], coslp, sinlp, parameter[offaxxill::a_spin]);

        param[offaxxillCp::Ecut] = parameter[offaxxill::Ecut] / ginf;
        param[offaxxillCp::refl_frac] = parameter[offaxxill::refl_frac];
        param[offaxxillCp::zshift] = parameter[offaxxill::zshift];
        param[offaxxillCp::switch_reflfrac_boost] = parameter[offaxxill::switch_reflfrac_boost];

        flux = offaxxillver::offaxxillver(energy, param, T_PrimSpec::CutoffPl);
    }
}
