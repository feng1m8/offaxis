#ifndef OFFAXIS_SPHERE_HXX
#define OFFAXIS_SPHERE_HXX

#include <map>
#include <vector>

namespace offaxis
{
    class Sphere
    {
    public:
        Sphere(long nside);

        auto operator[](std::size_t ipix) const
        {
            return (this->data.cbegin() + 3 * ipix);
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
