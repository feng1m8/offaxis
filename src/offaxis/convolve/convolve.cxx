#include <functional>
#include <memory>

#include "convolve.hxx"

void set_flux_outside_defined_range_to_zero(const double *ener, double *spec, int n_ener, double emin, double emax);

namespace offaxis::relxill
{
    void convolveSpectrumFFTNormalized(const std::valarray<double> &ener, const std::valarray<double> &frel, const double *ener0, std::valarray<double> &flu0)
    {
        int status = EXIT_SUCCESS;

        double fxill[N_ENER_CONV];
        rebin_spectrum(std::begin(envs::energy_conv), fxill, N_ENER_CONV, ener0, std::begin(flu0), flu0.size());

        double fout[N_ENER_CONV];
        std::unique_ptr<specCache, std::function<void(specCache *)>> spec_cache(new_specCache(), free_specCache);
        convolveSpectrumFFTNormalized(const_cast<double *>(std::begin(envs::energy_conv)), fxill, std::begin(frel), fout, N_ENER_CONV, 1, 1, 0, spec_cache.get(), &status);

        flu0.resize(ener.size() - 1);
        rebin_spectrum(std::begin(ener), std::begin(flu0), flu0.size(), std::begin(envs::energy_conv), fout, N_ENER_CONV);
    }

    void convolveSpectrumFFTNormalized(const std::valarray<double> &ener, const std::valarray<double> &frel, std::valarray<double> &flu0)
    {
        convolveSpectrumFFTNormalized(ener, frel, std::begin(ener), flu0);
        set_flux_outside_defined_range_to_zero(std::begin(ener), std::begin(flu0), flu0.size(), RELCONV_EMIN, RELCONV_EMAX);
    }
}
