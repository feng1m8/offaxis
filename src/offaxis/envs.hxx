#ifndef OFFAXIS_ENVS_HXX
#define OFFAXIS_ENVS_HXX

#include <filesystem>

#include <math.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832
#endif

namespace offaxis
{
    namespace utils
    {
        inline double deg2rad(double degree)
        {
            return M_PI / 180.0 * degree;
        }

        std::filesystem::path abspath();
    }

    namespace envs
    {
        int nthreads();
        long nside();
        std::filesystem::path kydir();
    }
}

#endif
