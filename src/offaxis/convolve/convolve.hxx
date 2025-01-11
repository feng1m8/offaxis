#ifndef OFFAXIS_CONVOLVE_CONVOLVE_HXX
#define OFFAXIS_CONVOLVE_CONVOLVE_HXX

#include <valarray>

#include "relxill/src/Relbase.h"

namespace offaxis
{
    namespace utils
    {
        inline std::valarray<double> geomspace(double start, double stop, std::size_t size)
        {
            std::valarray<double> geom(size);

            for (std::size_t i = 0; i < size; ++i)
                geom[i] = start * std::pow(stop / start, i / (size - 1.0));

            return geom;
        }
    }

    namespace envs
    {
        inline const std::valarray<double> energy_conv(utils::geomspace(EMIN_RELXILL_CONV, EMAX_RELXILL_CONV, N_ENER_CONV + 1));
    }

    namespace relxill
    {
        void convolveSpectrumFFTNormalized(const std::valarray<double> &ener, const std::valarray<double> &frel, const double *ener0, std::valarray<double> &flu0);
        void convolveSpectrumFFTNormalized(const std::valarray<double> &ener, const std::valarray<double> &frel, std::valarray<double> &flu0);
    }
}

#endif
