#include <cstdint>

#include <chealpix.h>

#include "sphere.hxx"

static std::vector<double> nside2vec_ring64(std::int64_t nside)
{
    std::vector<double> vecs;
    vecs.reserve(3 * nside2npix64(nside));

    for (std::int64_t pix = 0; pix < nside2npix64(nside); ++pix)
    {
        double vec[3];
        pix2vec_ring64(nside, pix, vec);
        vecs.emplace_back(vec[0]);
        vecs.emplace_back(vec[1]);
        vecs.emplace_back(vec[2]);
    }

    return vecs;
}

namespace offaxis
{
    Sphere::Sphere(int nside) : npix(nside2npix64(nside))
    {
        this->data = nside2vec_ring64(nside);
    }
}
