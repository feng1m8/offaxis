#ifndef OFFAXLINE_ENVS_HXX
#define OFFAXLINE_ENVS_HXX

#include <cmath>
#include <filesystem>

namespace offaxis
{
    namespace utils
    {
        inline constexpr double pi = 3.141592653589793238462643383279502884;

        inline double deg2rad(double degree)
        {
            return pi / 180.0 * degree;
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
