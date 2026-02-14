#include "offaxis/offaxis.h"
#include "offaxis/parameter.hxx"

#include "convolve.hxx"

namespace offaxis
{
    [[gnu::dllexport]] void offaxconv(const std::valarray<double> &energy, const std::valarray<double> &parameter, std::valarray<double> &flux)
    {
        using namespace parameter;
        if (parameter.size() < offaxconv::Nparams)
        {
            throw std::out_of_range("RealArray index " + std::to_string(offaxconv::Nparams - 1) + " is out of bounds with size " + std::to_string(parameter.size()) + ".");
        }

        std::valarray<double> param(offaxline::Nparams);
        param[offaxline::lineE] = 1.0;
        param[std::slice(offaxline::rlp, offaxline::gamma + 1, 1)] = parameter[std::slice(offaxconv::rlp, offaxconv::gamma + 1, 1)];
        param[offaxline::zshift] = 0.0;
        param[offaxline::normtype] = parameter[offaxconv::normtype];

        std::valarray<double> line;
        offaxis::offaxline(envs::energy_conv, param, line);

        relxill::convolveSpectrumFFTNormalized(energy, line, flux);
    }
}
