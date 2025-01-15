#ifndef OFFAXIS_REFLECTION_SPECTRUM_HXX
#define OFFAXIS_REFLECTION_SPECTRUM_HXX

#include <valarray>

#include "relxill/src/ModelInfo.h"

extern "C"
{
#include "relxill/src/common.h"
}

namespace offaxis::relxill
{
    int n_incl(T_PrimSpec prim_type);

    class Spectrum
    {
    public:
        struct Spec
        {
            std::valarray<double> energy;
            std::valarray<double> flux;
        };

        Spectrum(const std::valarray<double> &parameter, T_PrimSpec prim_type);

        Spectrum operator*(double ener_shift_observer_source) const;

        Spec xillver(const std::valarray<double> &energy, const std::valarray<double> &dist) const;

        std::valarray<double> convolve(Spec& xill_spec,const std::valarray<double> &energy, const std::valarray<double> &hist) const;

        double norm;

        std::valarray<double> nthcomp(const std::valarray<double> &ener, double ener_shift_source_obs) const;

    private:
        xillTableParam parameter;
    };
}

#endif
