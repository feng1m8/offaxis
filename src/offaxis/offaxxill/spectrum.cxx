#include <functional>
#include <memory>

#include "relxill/src/ModelDefinition.h"

#include "offaxconv/convolve.hxx"
#include "offaxis/parameter.hxx"
#include "spectrum.hxx"

namespace offaxis::relxill
{
    int n_incl(T_PrimSpec prim_type)
    {
        int status = EXIT_SUCCESS;

        xillTable *tab = nullptr;
        get_init_xillver_table(&tab, MOD_TYPE_RELXILLLP, convertPrimSpecType(prim_type), &status);

        return tab->n_incl;
    }

    Spectrum::Spectrum(const std::valarray<double> &parameter, T_PrimSpec prim_type)
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

        this->norm = 1.0;
    }

    Spectrum Spectrum::operator*(double ener_shift_observer_source) const
    {
        Spectrum spec(*this);
        spec.parameter.ect = this->parameter.ect * ener_shift_observer_source;
        spec.norm = calc_xillver_normalization_change(1.0 / ener_shift_observer_source, &spec.parameter);
        return spec;
    }

    Spectrum::Spec Spectrum::xillver(const std::valarray<double> &energy, const std::valarray<double> &dist) const
    {
        int status = EXIT_SUCCESS;

        std::unique_ptr<xillSpec, std::function<void(xillSpec *)>> xill_spec(get_xillver_spectra_table(&this->parameter, &status), free_xill_spec);
        std::valarray<double> xill_flux(xill_spec->n_ener);
        calc_xillver_angdep(std::begin(xill_flux), xill_spec.get(), std::begin(dist), &status);

        return {std::valarray<double>(xill_spec->ener, xill_spec->n_ener + 1), xill_flux};
    }

    std::valarray<double> Spectrum::convolve(Spectrum::Spec &xill_spec, const std::valarray<double> &energy, const std::valarray<double> &hist) const
    {
        relxill::convolveSpectrumFFTNormalized(energy, hist, xill_spec.energy, xill_spec.flux);
        return xill_spec.flux;
    }

    std::valarray<double> Spectrum::primary(const std::valarray<double> &ener, double ener_shift_source_obs) const
    {
        static const std::valarray<double> energy_norm(utils::geomspace(EMIN_XILLVER_NORMALIZATION, EMAX_XILLVER_NORMALIZATION, N_ENER_COARSE + 1));

        int status = EXIT_SUCCESS;

        double prim_spec_source[N_ENER_COARSE];
        calc_primary_spectrum(prim_spec_source, std::begin(energy_norm), N_ENER_COARSE, &this->parameter, &status);
        double norm = 2.0 / calcXillverNormFromPrimarySpectrum(prim_spec_source, std::begin(energy_norm), N_ENER_COARSE, nullptr);

        std::valarray<double> flux(ener.size() - 1);
        calc_primary_spectrum(std::begin(flux), std::begin(ener), flux.size(), &this->parameter, &status, ener_shift_source_obs);

        return norm * flux;
    }
}
