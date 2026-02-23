#include <functional>
#include <memory>

#include "relxill/src/Relbase.h"

#include "offaxconv/convolve.hxx"

void set_flux_outside_defined_range_to_zero(const double *ener, double *spec, int n_ener, double emin, double emax);

namespace offaxis::relxill
{
    static specCache *new_specCache()
    {
        auto *spec = new specCache;

        spec->n_cache = 1;
        spec->nzones = 0;
        spec->n_ener = N_ENER_CONV;

        spec->conversion_factor_energyflux = nullptr;

        spec->fft_xill = new double **[1];
        spec->fft_rel = new double **[1];

        spec->fftw_xill = new fftw_complex *[1];
        spec->fftw_rel = new fftw_complex *[1];

        spec->fftw_backwards_input = fftw_alloc_complex(N_ENER_CONV);
        spec->fftw_output = new double[N_ENER_CONV];

        spec->plan_c2r = fftw_plan_dft_c2r_1d(N_ENER_CONV, spec->fftw_backwards_input, spec->fftw_output, FFTW_ESTIMATE);

        spec->xill_spec = new xillSpec *[1];

        spec->fft_xill[0] = new double *[2];
        spec->fft_rel[0] = new double *[2];

        spec->fftw_xill[0] = fftw_alloc_complex(N_ENER_CONV);
        spec->fftw_rel[0] = fftw_alloc_complex(N_ENER_CONV);

        spec->fft_xill[0][0] = new double[N_ENER_CONV];
        spec->fft_rel[0][0] = new double[N_ENER_CONV];
        spec->fft_xill[0][1] = new double[N_ENER_CONV];
        spec->fft_rel[0][1] = new double[N_ENER_CONV];

        spec->xill_spec[0] = nullptr;

        spec->out_spec = nullptr;

        return spec;
    }

    void convolveSpectrumFFTNormalized(const std::vector<double> &ener, const std::valarray<double> &frel, const std::valarray<double> &ener0, std::valarray<double> &flu0)
    {
        int status = EXIT_SUCCESS;

        double fxill[N_ENER_CONV];
        _rebin_spectrum(std::begin(envs::energy_conv), fxill, N_ENER_CONV, std::begin(ener0), std::begin(flu0), flu0.size());

        double fout[N_ENER_CONV];
        std::unique_ptr<specCache, std::function<void(specCache *)>> spec_cache(new_specCache(), free_specCache);
        convolveSpectrumFFTNormalized(std::begin(envs::energy_conv), fxill, std::begin(frel), fout, N_ENER_CONV, 1, 1, 0, spec_cache.get(), &status);

        flu0.resize(ener.size() - 1);
        _rebin_spectrum(ener.data(), std::begin(flu0), flu0.size(), std::begin(envs::energy_conv), fout, N_ENER_CONV);
    }

    void convolveSpectrumFFTNormalized(const std::valarray<double> &ener, const std::valarray<double> &frel, std::valarray<double> &flu0)
    {
        convolveSpectrumFFTNormalized({std::begin(ener), std::end(ener)}, frel, ener, flu0);
        set_flux_outside_defined_range_to_zero(std::begin(ener), std::begin(flu0), flu0.size(), RELCONV_EMIN, RELCONV_EMAX);
    }
}
