#ifndef OFFAXIS_REFLECTION_REFLECTION_HXX
#define OFFAXIS_REFLECTION_REFLECTION_HXX

#include <valarray>
#include <vector>
#include "relxill/src/ModelInfo.h"

namespace offaxis::offaxxillver
{
    double doppler(double coslp, double sinlp, double philp, double *vlp, double incl);
    double redshiftinf(double rlp, double coslp, double sinlp, double a_spin);

    std::valarray<double> primary(const std::valarray<double> &energy, const std::valarray<double> &parameter, T_PrimSpec prim_type);

    struct Reflection
    {
        double f_refl;
        std::vector<double> glp;
        std::vector<std::valarray<double>> hist;
        std::vector<std::valarray<double>> dist;
    };

    class Emission
    {
    public:
        Emission(std::size_t npix)
        {
            this->glp.resize(npix);
            this->gobs.resize(npix);
            this->cosem.resize(npix);
            this->iobs.resize(npix);
        };

        auto at(std::size_t i)
        {
            return std::tie(glp[i], gobs[i], cosem[i], iobs[i]);
        }

        Reflection get(double f_refl, int n_incl);

    private:
        std::valarray<double> glp;
        std::valarray<double> gobs;
        std::valarray<double> cosem;
        std::valarray<double> iobs;

        std::size_t resize();
    };
}

#endif
