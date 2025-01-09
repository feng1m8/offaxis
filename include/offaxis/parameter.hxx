#ifndef OFFAXIS_PARAMETER_HXX
#define OFFAXIS_PARAMETER_HXX

#include <cstddef>

namespace offaxis::parameter
{
    namespace offaxline
    {
        enum : std::size_t
        {
            lineE,
            rlp,
            thetalp,
            philp,
            vr,
            vtheta,
            vphi,
            a_spin,
            Incl,
            Rin,
            Rout,
            gamma,
            zshift,
            normtype,

            Nparam,
        };
    }

    namespace offaxconv
    {
        enum : std::size_t
        {
            rlp,
            thetalp,
            philp,
            vr,
            vtheta,
            vphi,
            a_spin,
            Incl,
            Rin,
            Rout,
            gamma,
            normtype,

            Nparam,
        };
    }

    namespace offaxxill
    {
        enum : std::size_t
        {
            rlp,
            thetalp,
            philp,
            vr,
            vtheta,
            vphi,
            a_spin,
            Incl,
            Rin,
            Rout,
            gamma,
            logxi,
            Afe,
            Ecut,
            refl_frac,
            zshift,
            switch_reflfrac_boost,

            Nparam,
        };
    }

    namespace offaxxillCp
    {
        enum : std::size_t
        {
            rlp,
            thetalp,
            philp,
            vr,
            vtheta,
            vphi,
            a_spin,
            Incl,
            Rin,
            Rout,
            gamma,
            logxi,
            logN,
            Afe,
            Ecut,
            refl_frac,
            zshift,
            switch_reflfrac_boost,

            Nparam,
        };
    }
}

#endif
