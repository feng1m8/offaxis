#include <algorithm>

#include "relxill/src/Relbase.h"

#include "emission.hxx"

namespace offaxis::offaxxillver
{
    void Histogram::accumulate(double gobs, double iobs, double cosem, double glp)
    {
        static const double STEP = N_ENER_CONV / std::log(EMAX_RELXILL_CONV / EMIN_RELXILL_CONV);
        static const double LOWER = std::log(EMIN_RELXILL_CONV);

        this->data.emplace_back(Data{
            static_cast<std::size_t>(STEP * (std::log(gobs) - LOWER)),
            iobs,
            static_cast<std::size_t>(cosem * this->n_incl),
            glp,
        });

        this->itot += iobs;
    }

    void Histogram::operator+=(const Histogram &rhs)
    {
        this->data.insert(this->data.cend(), rhs.data.cbegin(), rhs.data.cend());

        if (std::is_sorted(this->data.cbegin(), this->data.cend() - rhs.data.size()) and std::is_sorted(rhs.data.cbegin(), rhs.data.cend()))
            std::inplace_merge(this->data.begin(), this->data.end() - rhs.data.size(), this->data.end());
        else
            std::sort(this->data.begin(), this->data.end());

        this->itot += rhs.itot;
        this->to_disk += rhs.to_disk;
        this->to_infinity += rhs.to_infinity;
    }

    Emission Histogram::get()
    {
        if (not std::is_sorted(this->data.cbegin(), this->data.cend()))
            std::sort(this->data.begin(), this->data.end());

        Emission emission;
        emission.glp.reserve(N_ZONES);
        emission.hist.reserve(N_ZONES);
        emission.dist.reserve(N_ZONES);

        const double i_mean = this->itot / N_ZONES;
        {
            double iobs = 0.0;
            std::size_t i = 0;
            for (int j = 1; j < N_ZONES + 1; ++j)
            {
                double glp_i = 0.0;
                double i_zone = 0.0;
                std::valarray<double> hist(N_ENER_CONV);
                std::valarray<double> dist(n_incl + 1);

                for (; i < this->data.size(); ++i)
                {
                    if (iobs >= j * i_mean and j < N_ZONES)
                        break;

                    iobs += this->data[i].iobs;
                    glp_i += this->data[i].glp * this->data[i].iobs;
                    i_zone += this->data[i].iobs;
                    hist[this->data[i].gobs] += this->data[i].iobs;
                    dist[this->data[i].cosem] += this->data[i].iobs;
                }

                emission.glp.emplace_back(glp_i / i_zone);
                emission.hist.emplace_back(hist / this->size);

                dist /= i_zone;
                dist[this->n_incl - 1] += dist[this->n_incl];
                std::reverse(std::begin(dist), std::end(dist) - 1);
                emission.dist.emplace_back(std::begin(dist), this->n_incl);
            }
        }

        emission.f_refl = static_cast<double>(this->to_infinity) / static_cast<double>(this->to_disk);

        return emission;
    }
}
