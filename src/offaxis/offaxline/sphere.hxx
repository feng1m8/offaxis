#ifndef OFFAXLINE_SPHERE_HXX
#define OFFAXLINE_SPHERE_HXX

#include <map>
#include <vector>

namespace offaxis
{
    class Sphere
    {
    public:
        Sphere(long nside);

        auto operator[](std::size_t index) const
        {
            auto i = this->data.cbegin() + 3 * index;
            return std::tie(i[0], i[1], i[2]);
        }

        const std::size_t size;

    private:
        std::vector<double> data;
    };

    namespace envs
    {
        long nside();

        inline std::map<long, const Sphere> sphere;
    }
}

#endif
