#include <algorithm>

#include "relxill/src/Xillspec.h"

#include "emission.hxx"

namespace offaxis::offaxxillver
{
    Histogram::Histogram(std::size_t max_size, int n_incl) : max_size(max_size), n_incl(n_incl)
    {
        this->data.reserve(max_size);
        this->i_total = 0.0;
        this->to_disk = 0;
        this->to_infinity = 0;
    }

    void Histogram::accumulate(double gobs, double iobs, double glp, double cosem)
    {
        static const double STEP = N_ENER_CONV / std::log(EMAX_RELXILL_CONV / EMIN_RELXILL_CONV);
        static const double LOWER = std::log(EMIN_RELXILL_CONV);

        this->data.emplace_back(Data{
            .glp = glp,
            .iobs = iobs,
            .gobs = static_cast<std::size_t>(STEP * (std::log(gobs) - LOWER)),
            .cosem = static_cast<std::size_t>(cosem * this->n_incl),
        });

        this->i_total += iobs;
    }

    void Histogram::operator+=(const Histogram &rhs)
    {
        this->data.insert(this->data.cend(), rhs.data.cbegin(), rhs.data.cend());
        this->i_total += rhs.i_total;
        this->to_disk += rhs.to_disk;
        this->to_infinity += rhs.to_infinity;
    }

    Emission Histogram::get()
    {
        std::sort(this->data.begin(), this->data.end());

        Emission emission;
        emission.glp.resize(N_ZONES);
        emission.hist.resize(N_ZONES, std::valarray<double>(N_ENER_CONV));
        emission.dist.resize(N_ZONES, std::vector<double>(this->n_incl + 1));

        const double i_mean = this->i_total / N_ZONES;

        std::size_t i = 0;
        double i_sum = 0.0;

        for (std::size_t j = 0; j < N_ZONES; ++j)
        {
            double i_j = 0.0;

            double &glp = emission.glp[j];
            std::valarray<double> &hist = emission.hist[j];
            std::vector<double> &dist = emission.dist[j];

            for (int k = j + 1; i < this->data.size(); ++i)
            {
                if (i_sum >= k * i_mean and k < N_ZONES)
                    break;

                i_sum += this->data[i].iobs;
                i_j += this->data[i].iobs;
                hist[this->data[i].gobs] += this->data[i].iobs;
                dist[this->data[i].cosem] += this->data[i].iobs;
                glp += this->data[i].glp * this->data[i].iobs;
            }

            dist[this->n_incl - 1] += dist[this->n_incl];
            dist.pop_back();
            std::reverse(dist.begin(), dist.end());
            for (auto &i : dist)
                i /= i_j;

            glp /= i_j;
            hist /= this->max_size;
        }

        if (this->to_disk == 0)
            emission.f_refl = 1.0;
        else
            emission.f_refl = static_cast<double>(this->to_infinity) / static_cast<double>(this->to_disk);

        return emission;
    }
}
