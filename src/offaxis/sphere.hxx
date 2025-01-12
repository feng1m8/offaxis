#ifndef OFFAXIS_SPHERE_HXX
#define OFFAXIS_SPHERE_HXX

#include <array>
#include <functional>
#include <map>
#include <vector>

namespace offaxis
{
    class Sphere
    {
    public:
        const std::size_t size;

        Sphere(int nside);

        std::array<std::reference_wrapper<const double>, 3> operator[](std::size_t pix) const
        {
            auto vec = this->data.cbegin() + 3 * pix;
            return {vec[0], vec[1], vec[2]};
        }

    private:
        std::vector<double> data;
    };

    namespace envs
    {
        inline std::map<int, const Sphere> sphere;
    }
}

#endif
