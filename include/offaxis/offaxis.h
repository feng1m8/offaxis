#ifndef OFFAXIS_OFFAXIS_H
#define OFFAXIS_OFFAXIS_H

#ifdef __cplusplus

#include <valarray>

namespace offaxis
{
    void offaxline(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux);
    void offaxconv(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux);
    void offaxxill(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux);
    void offaxxillCp(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux);
}

extern "C"
{
#endif

    void coffaxline(const double *energy, int Nflux, const double *parameter, int spectrum, double *flux, double *fluxVariance, const char *init);
    void coffaxconv(const double *energy, int Nflux, const double *parameter, int spectrum, double *flux, double *fluxVariance, const char *init);
    void coffaxxill(const double *energy, int Nflux, const double *parameter, int spectrum, double *flux, double *fluxVariance, const char *init);
    void coffaxxillCp(const double *energy, int Nflux, const double *parameter, int spectrum, double *flux, double *fluxVariance, const char *init);

#ifdef __cplusplus
}
#endif

#endif
