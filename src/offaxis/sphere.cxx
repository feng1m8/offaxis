#include <cstdint>

#include <chealpix.h>

#include "sphere.hxx"

namespace offaxis
{
    Sphere::Sphere(int nside) : size(nside2npix64(nside))
    {
        this->data.reserve(3 * this->size);

        for (std::int64_t pix = 0; pix < static_cast<std::int64_t>(this->size); ++pix)
        {
            double vec[3];
            pix2vec_ring64(nside, pix, vec);

            this->data.emplace_back(vec[0]);
            this->data.emplace_back(vec[1]);
            this->data.emplace_back(vec[2]);
        }
    }
}
