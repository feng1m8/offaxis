#include <functional>
#include <memory>

#include "relxill/src/ModelDefinition.h"

#include "offaxconv/convolve.hxx"
#include "offaxis/parameter.hxx"
#include "offaxxill/doppler.hxx"
#include "offaxxill/emission.hxx"
#include "offaxxill/spectrum.hxx"

namespace offaxis::relxill
{
    int n_incl(T_PrimSpec prim_type)
    {
        int status = EXIT_SUCCESS;

        xillTable *tab = nullptr;
        get_init_xillver_table(&tab, MOD_TYPE_RELXILLLP, convertPrimSpecType(prim_type), &status);

        return tab->n_incl;
    }

    static std::tuple<std::valarray<double>, std::valarray<double>> calc_xillver_angdep(const xillTableParam &parameter, const std::vector<double> &dist)
    {
        int status = EXIT_SUCCESS;

        std::unique_ptr<xillSpec, std::function<void(xillSpec *)>> xill_spec(get_xillver_spectra_table(&parameter, &status), free_xill_spec);
        std::valarray<double> xill_flux(xill_spec->n_ener);
        calc_xillver_angdep(std::begin(xill_flux), xill_spec.get(), dist.data(), &status);

        return {
            std::valarray<double>(xill_spec->ener, xill_spec->n_ener + 1),
            xill_flux,
        };
    }

    static std::valarray<double> get_observed_primary_spectrum(const std::vector<double> &energy, const xillTableParam &parameter, double ener_shift_source_obs)
    {
        static const std::valarray<double> energy_norm(utils::geomspace(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, N_ENER_COARSE + 1));

        int status = EXIT_SUCCESS;

        double prim_spec_source[N_ENER_COARSE];
        calc_primary_spectrum(prim_spec_source, std::begin(energy_norm), N_ENER_COARSE, &parameter, &status);
        double norm = 2.0 / calcXillverNormFromPrimarySpectrum(prim_spec_source, std::begin(energy_norm), N_ENER_COARSE, nullptr);

        std::valarray<double> flux(energy.size() - 1);
        calc_primary_spectrum(std::begin(flux), energy.data(), flux.size(), &parameter, &status, ener_shift_source_obs);

        return norm * flux;
    }

    Spectrum::Spectrum(const std::vector<double> &energy, const std::vector<double> &parameter, T_PrimSpec prim_type) : energy(energy)
    {
        using namespace parameter::offaxxillCp;

        this->parameter.gam = parameter[gamma];
        this->parameter.afe = parameter[Afe];
        this->parameter.lxi = parameter[logxi];
        this->parameter.ect = parameter[Ecut];
        this->parameter.incl = parameter[Incl];
        this->parameter.dens = parameter[logN];
        this->parameter.frac_pl_bb = 0.0;
        this->parameter.kTbb = 0.0;
        this->parameter.prim_type = convertPrimSpecType(prim_type);
        this->parameter.model_type = MOD_TYPE_RELXILLLP;
    }

    std::valarray<double> Spectrum::convolve(const offaxis::offaxxillver::Emission &emission) const
    {
        std::valarray<double> flux(this->energy.size() - 1);

        for (std::size_t i = 0; i < emission.glp.size(); ++i)
        {
            xillTableParam params(this->parameter);
            params.ect *= emission.glp[i];
            auto [engs, spec] = calc_xillver_angdep(params, emission.dist[i]);

            spec *= calc_xillver_normalization_change(1.0 / emission.glp[i], &params);

            relxill::convolveSpectrumFFTNormalized(this->energy, emission.hist[i], engs, spec);
            flux += spec;
        }

        return flux * emission.f_refl;
    }

    std::valarray<double> Spectrum::primary(const std::vector<double> &parameter) const
    {
        auto [gobs, iobs] = beaming(parameter);

        std::valarray flux(get_observed_primary_spectrum(this->energy, this->parameter, gobs));

        return iobs * flux;
    }

}
