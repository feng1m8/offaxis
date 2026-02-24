#include "offaxconv/convolve.hxx"
#include "offaxis/offaxis.h"
#include "offaxis/parameter.hxx"
#include "offaxline/envs.hxx"

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

        std::valarray<double> params(offaxline::Nparams);
        params[offaxline::lineE] = 1.0;
        params[offaxline::rlp] = parameter[rlp];
        params[offaxline::thetalp] = parameter[thetalp];
        params[offaxline::philp] = parameter[philp];
        params[offaxline::vr] = parameter[vr];
        params[offaxline::vtheta] = parameter[vtheta];
        params[offaxline::vphi] = parameter[vphi];
        params[offaxline::a_spin] = parameter[a_spin];
        params[offaxline::Incl] = parameter[Incl];
        params[offaxline::Rin] = parameter[Rin];
        params[offaxline::Rout] = parameter[Rout];
        params[offaxline::gamma] = parameter[gamma];
        params[offaxline::zshift] = 0.0;
        params[offaxline::normtype] = parameter[normtype];

        std::valarray<double> line;
        offaxis::offaxline(envs::energy_conv, params, line);

        relxill::convolveSpectrumFFTNormalized(energy, line, flux);
    }
}
