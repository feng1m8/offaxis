#ifndef OFFAXIS_REFLECTION_PRIMARY_HXX
#define OFFAXIS_REFLECTION_PRIMARY_HXX

#include <vector>

namespace offaxis
{
    double doppler(double coslp, double sinlp, double philp, const double *vlp, double incl);
    double expnu(double rlp, double coslp, double sinlp, double a_spin);

    std::tuple<double, double> beaming(const std::vector<double> &parameter);
}

#endif
