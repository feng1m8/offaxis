#ifndef OFFAXIS_REFLECTION_PRIMARY_HXX
#define OFFAXIS_REFLECTION_PRIMARY_HXX

#include <valarray>

#include "relxill/src/ModelInfo.h"

namespace offaxis::offaxxillver
{
    double doppler(double coslp, double sinlp, double philp, const double *vlp, double incl);
    double redshiftinf(double rlp, double coslp, double sinlp, double a_spin);

    std::valarray<double> primary(const std::valarray<double> &energy, const std::valarray<double> &parameter, T_PrimSpec prim_type);
}

#endif
