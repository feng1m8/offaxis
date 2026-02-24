#ifndef OFFAXXILL_DOPPLER_HXX
#define OFFAXXILL_DOPPLER_HXX

#include <valarray>
#include <vector>

namespace offaxis
{
    double redshift_primary(const std::valarray<double> &parameter);
    std::tuple<double, double> doppler_primary(const std::vector<double> &parameter);
}

#endif
