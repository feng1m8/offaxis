#ifndef OFFAXXILL_SPECTRUM_HXX
#define OFFAXXILL_SPECTRUM_HXX

#include <valarray>
#include <vector>

#include "relxill/src/ModelInfo.h"

extern "C"
{
#include "relxill/src/common.h"
}

namespace offaxis
{
    namespace offaxxillver
    {
        class Emission;
    }

    namespace relxill
    {
        int n_incl(T_PrimSpec prim_type);

        class Spectrum
        {
        public:
            Spectrum(const std::vector<double> &energy, const std::vector<double> &parameter, T_PrimSpec prim_type);

            std::valarray<double> convolve(const offaxxillver::Emission &emission) const;

            std::valarray<double> primary(const std::vector<double> &parameter) const;

        private:
            const std::vector<double> &energy;
            xillTableParam parameter;
        };
    }
}

#endif
