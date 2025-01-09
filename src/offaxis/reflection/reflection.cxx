#include <algorithm>
#include <array>
#include <numeric>

#include <gsl/gsl_statistics_double.h>

#include "relxill/src/Relbase.h"

#include "reflection.hxx"

static std::valarray<double> hist_get(const double *gobs, const double *iobs, std::size_t size)
{
    std::valarray<double> hist(N_ENER_CONV);

    static const double STEP = N_ENER_CONV / std::log(EMAX_RELXILL_CONV / EMIN_RELXILL_CONV);
    static const double LOWER = std::log(EMIN_RELXILL_CONV);

    for (std::size_t i = 0; i < size; ++i)
        hist[STEP * (std::log(gobs[i]) - LOWER)] += iobs[i];

    return hist;
}

static std::valarray<double> dist_get(const double *cosem, const double *iobs, std::size_t size, int n_incl)
{
    std::valarray<double> dist(n_incl + 1);

    for (std::size_t i = 0; i < size; ++i)
        dist[cosem[i] * n_incl] += iobs[i];

    dist[n_incl - 1] += dist[n_incl];

    dist = std::valarray<double>(std::begin(dist), n_incl);
    std::reverse(std::begin(dist), std::end(dist));

    return dist / dist.sum();
}

namespace offaxis
{
    namespace utils
    {
        static std::valarray<std::size_t> argsort(const std::valarray<double> &arr, const std::valarray<bool> &mask)
        {
            std::valarray<std::size_t> index(arr.size());
            std::iota(std::begin(index), std::end(index), 0);
            index = std::valarray<std::size_t>(index[mask]);

            std::sort(
                std::begin(index), std::end(index),
                [&arr](std::size_t i1, std::size_t i2)
                {
                    return arr[i1] < arr[i2];
                });

            return index;
        }

        static std::array<std::size_t, N_ZONES + 1> argquantile(const std::valarray<double> &weight)
        {
            std::array<std::size_t, N_ZONES + 1> index;
            index[0] = 0;

            double bin = weight.sum() / N_ZONES;

            double sum = 0.0;
            int j = 1;
            for (std::size_t i = 0; i < weight.size(); ++i)
            {
                sum += weight[i];

                if (sum > j * bin)
                {
                    index[j] = i;
                    ++j;
                }
            }

            index[N_ZONES] = weight.size();

            return index;
        }
    }

    namespace offaxxillver
    {
        std::size_t Emission::resize()
        {
            std::size_t npix = this->glp.size();

            auto index(utils::argsort(this->glp, this->glp > 0.0));

            this->glp = std::valarray<double>(this->glp[index]);
            this->gobs = std::valarray<double>(this->gobs[index]);
            this->cosem = std::valarray<double>(this->cosem[index]);
            this->iobs = std::valarray<double>(this->iobs[index]);

            return npix;
        }

        Reflection Emission::get(double f_refl, int n_incl)
        {
            Reflection reflection;
            reflection.f_refl = f_refl;

            std::size_t npix = this->resize();
            auto index(utils::argquantile(this->iobs));

            reflection.glp.reserve(N_ZONES);
            reflection.hist.reserve(N_ZONES);
            reflection.dist.reserve(N_ZONES);
            for (std::size_t i = 0; i < N_ZONES; ++i)
            {
                std::size_t j = index[i];
                std::size_t length = index[i + 1] - index[i];

                reflection.glp.emplace_back(gsl_stats_wmean(std::begin(this->iobs) + j, 1, std::begin(this->glp) + j, 1, length));
                reflection.hist.emplace_back(hist_get(std::begin(this->gobs) + j, std::begin(this->iobs) + j, length) / npix);
                reflection.dist.emplace_back(dist_get(std::begin(this->cosem) + j, std::begin(this->iobs) + j, length, n_incl));
            }

            return reflection;
        }
    }
}