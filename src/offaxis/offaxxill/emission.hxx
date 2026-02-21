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
        std::vector<std::valarray<double>> dist;
    };

    class Histogram
    {
    public:
        Histogram(int n_incl, std::size_t size)
        {
            this->to_disk = 0;
            this->to_infinity = 0;
            this->itot = 0.0;
            this->n_incl = n_incl;
            this->size = size;
            this->data.reserve(size);
        }

        void accumulate(double gobs, double iobs, double cosem, double glp);

        void operator+=(const Histogram &rhs);

        Emission get();

        std::size_t to_disk;
        std::size_t to_infinity;

    private:
        struct Data
        {
            std::size_t gobs;
            double iobs;
            std::size_t cosem;
            double glp;

            bool operator<(const Data &rhs) const
            {
                return this->glp < rhs.glp;
            }
        };

        double itot;
        std::vector<Data> data;
        int n_incl;
        std::size_t size;
    };
}

#endif
