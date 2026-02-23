#include <cstdlib>

#include <chealpix.h>

#include "offaxline/sphere.hxx"

namespace offaxis
{
    namespace envs
    {
        long nside()
        {
            auto env = std::getenv("OFFAXIS_NUM_SIDE");
            if (env == nullptr)
                return 64;
            else
                return std::atol(env);
        }
    }

    Sphere::Sphere(long nside) : size(nside2npix64(nside)), data(3 * this->size)
    {
        for (std::size_t i = 0; i < this->size; ++i)
            pix2vec_ring64(nside, i, 3 * i + this->data.data());
    }
}
