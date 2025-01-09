#ifndef OFFAXIS_SPHERE_HXX
#define OFFAXIS_SPHERE_HXX

#include <array>
#include <map>
#include <vector>

namespace offaxis
{
    class Sphere
    {
    public:
        Sphere(int nside);

        std::array<double, 3> operator[](std::size_t pix) const
        {
            auto vec = this->data.cbegin() + 3 * pix;
            return {vec[0], vec[1], vec[2]};
        }

        const std::size_t npix;

    private:
        std::vector<double> data;
    };

    namespace envs
    {
        inline std::map<int, const Sphere> sphere;
    }
}

#endif
