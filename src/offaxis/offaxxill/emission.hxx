#ifndef OFFAXXILL_EMISSION_HXX
#define OFFAXXILL_EMISSION_HXX

#include <valarray>
#include <vector>

namespace offaxis::offaxxillver
{
    struct Emission
    {
        std::vector<std::valarray<double>> hist;
        std::vector<std::vector<double>> dist;
        std::vector<double> glp;
        double f_refl;
    };

    class Histogram
    {
    private:
        struct Data
        {
            double glp;
            double iobs;
            std::size_t gobs;
            std::size_t cosem;

            bool operator<(const Data &rhs) const
            {
                return this->glp < rhs.glp;
            }
        };

        std::vector<Data> data;
        double i_total;

        const int n_incl;
        const std::size_t max_size;

    public:
        Histogram(std::size_t max_size, int n_incl);

        void accumulate(double gobs, double iobs, double glp, double cosem);

        void operator+=(const Histogram &rhs);

        Emission get();

        std::size_t to_disk;
        std::size_t to_infinity;
    };
}

#endif
