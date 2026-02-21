#ifndef OFFAXLINE_ENVS_HXX
#define OFFAXLINE_ENVS_HXX

#include <cmath>
#include <filesystem>

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

        inline double mod(double x1, double x2)
        {
            return x1 - x2 * std::floor(x1 / x2);
        }

        std::filesystem::path abspath();
        std::string fotmat(const char *format, ...);
    }

    namespace envs
    {
        int nthreads();
    }
}

#endif
