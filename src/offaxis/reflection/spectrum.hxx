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
        Spectrum(const std::valarray<double> &parameter, T_PrimSpec prim_type);

        double norm(double ener_shift_observer_source) const;

        std::valarray<double> angdep(const std::valarray<double> &energy, const std::valarray<double> &hist, const std::valarray<double> &dist) const;

        std::valarray<double> nthcomp(const std::valarray<double> &ener, double ener_shift_source_obs) const;

        double &ect;

    private:
        xillTableParam param;
    };
}

#endif
