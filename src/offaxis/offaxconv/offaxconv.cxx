#include "offaxis/offaxis.h"
#include "offaxis/parameter.hxx"

#include "offaxline/envs.hxx"

#include "convolve.hxx"

namespace offaxis
{
    [[gnu::dllexport]] void offaxconv(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        using namespace parameter;
        using namespace parameter::offaxconv;

        if (parameter.size() < Nparams)
        {
            throw std::out_of_range(utils::fotmat("RealArray index %zu is out of bounds with size %zu.", Nparams - 1, parameter.size()));
        }

        std::valarray<double> param(offaxline::Nparams);
        param[offaxline::lineE] = 1.0;
        param[offaxline::rlp] = parameter[rlp];
        param[offaxline::thetalp] = parameter[thetalp];
        param[offaxline::philp] = parameter[philp];
        param[offaxline::vtheta] = parameter[vtheta];
        param[offaxline::vphi] = parameter[vphi];
        param[offaxline::a_spin] = parameter[a_spin];
        param[offaxline::Incl] = parameter[Incl];
        param[offaxline::Rin] = parameter[Rin];
        param[offaxline::Rout] = parameter[Rout];
        param[offaxline::gamma] = parameter[gamma];
        param[offaxline::zshift] = 0.0;
        param[offaxline::normtype] = parameter[normtype];

        std::valarray<double> line;
        offaxis::offaxline(envs::energy_conv, param, line);

        relxill::convolveSpectrumFFTNormalized(energy, line, flux);
    }
}
