#include "relxill/src/Relbase.h"

namespace offaxis::relxill
{
    specCache *new_specCache()
    {
        auto *spec = new specCache;

        spec->n_cache = 1;
        spec->nzones = 0;
        spec->n_ener = N_ENER_CONV;

        spec->conversion_factor_energyflux = nullptr;

        spec->fft_xill = new double**[1];
        spec->fft_rel = new double**[1];

        spec->fftw_xill = new fftw_complex*[1];
        spec->fftw_rel = new fftw_complex*[1];

        spec->fftw_backwards_input = fftw_alloc_complex(N_ENER_CONV);
        spec->fftw_output = new double[N_ENER_CONV];

        spec->plan_c2r = fftw_plan_dft_c2r_1d(N_ENER_CONV, spec->fftw_backwards_input, spec->fftw_output, FFTW_ESTIMATE);

        spec->xill_spec = new xillSpec*[1];

        spec->fft_xill[0] = new double*[2];
        spec->fft_rel[0] = new double*[2];

        spec->fftw_xill[0] = fftw_alloc_complex(N_ENER_CONV);
        spec->fftw_rel[0] = fftw_alloc_complex(N_ENER_CONV);

        for (int jj = 0; jj < 2; jj++)
        {
            spec->fft_xill[0][jj] = new double[N_ENER_CONV];
            spec->fft_rel[0][jj] = new double[N_ENER_CONV];
        }
        
        spec->xill_spec[0] = nullptr;
        spec->out_spec = nullptr;

        return spec;
    }
}
