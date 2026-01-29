#include <chealpix.h>

#include "sphere.hxx"

namespace offaxis
{
    Sphere::Sphere(long nside) : size(nside2npix64(nside)), data(3 * this->size)
    {
        for (std::size_t i = 0; i < this->size; ++i)
            pix2vec_ring64(nside, i, 3 * i + this->data.data());
    }
}
