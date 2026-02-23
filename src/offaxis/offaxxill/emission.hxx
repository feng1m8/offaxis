#ifndef OFFAXIS_REFLECTION_EMISSION_HXX
#define OFFAXIS_REFLECTION_EMISSION_HXX

#include <valarray>
#include <vector>

namespace offaxis::offaxxillver
{
    struct Emission
    {
        double f_refl;
        std::vector<double> glp;
        std::vector<std::valarray<double>> hist;
        std::vector<std::vector<double>> dist;
    };

    class Histogram
    {
    public:
        Histogram(std::size_t max_size, int n_incl)
        {
            this->to_disk = 0;
            this->to_infinity = 0;
            this->i_total = 0.0;
            this->n_incl = n_incl;
            this->max_size = max_size;
            this->data.reserve(max_size);
        }

        void accumulate(double gobs, double iobs, double glp, double cosem);

        void operator+=(const Histogram &rhs);

        Emission get();

        std::size_t to_disk;
        std::size_t to_infinity;

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

        double i_total;
        std::vector<Data> data;
        int n_incl;
        std::size_t max_size;
    };
}

#endif
