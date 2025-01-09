#include "offaxis/offaxis.h"
#include "offaxis/parameter.hxx"

template <auto F, std::size_t N>
static void xspecmodel(const double *energy, int Nflux, const double *parameter, double *flux)
{
    std::valarray<double> ener(energy, Nflux + 1);
    std::valarray<double> param(parameter, N);
    std::valarray<double> cflux(flux, Nflux);

    F(ener, param, cflux);

    std::copy(std::begin(cflux), std::end(cflux), flux);
}

extern "C"
{
    void coffaxline(const double *energy, int Nflux, const double *parameter, int spectrum, double *flux, double *fluxVariance, const char *init)
    {
        xspecmodel<offaxis::offaxline, offaxis::parameter::offaxline::Nparam>(energy, Nflux, parameter, flux);
    }

    void coffaxconv(const double *energy, int Nflux, const double *parameter, int spectrum, double *flux, double *fluxVariance, const char *init)
    {
        xspecmodel<offaxis::offaxconv, offaxis::parameter::offaxconv::Nparam>(energy, Nflux, parameter, flux);
    }

    void coffaxxill(const double *energy, int Nflux, const double *parameter, int spectrum, double *flux, double *fluxVariance, const char *init)
    {
        xspecmodel<offaxis::offaxxill, offaxis::parameter::offaxxill::Nparam>(energy, Nflux, parameter, flux);
    }

    void coffaxxillCp(const double *energy, int Nflux, const double *parameter, int spectrum, double *flux, double *fluxVariance, const char *init)
    {
        xspecmodel<offaxis::offaxxillCp, offaxis::parameter::offaxxillCp::Nparam>(energy, Nflux, parameter, flux);
    }
}
