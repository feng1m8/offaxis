#ifndef OFFAXIS_HISTOGRAM_HXX
#define OFFAXIS_HISTOGRAM_HXX

#include <valarray>

#include <gsl/gsl_histogram.h>

namespace offaxis
{
    class Histogram
    {
    private:
        gsl_histogram *hist;

    public:
        Histogram(const std::valarray<double> &range)
        {
            this->hist = gsl_histogram_alloc(range.size() - 1);
            gsl_histogram_set_ranges(this->hist, std::begin(range), range.size());
        }

        ~Histogram()
        {
            gsl_histogram_free(this->hist);
        }

        void accumulate(double x, double weight) const
        {
            gsl_histogram_accumulate(this->hist, x, weight);
        }

        std::valarray<double> get() const
        {
            return std::valarray<double>(this->hist->bin, this->hist->n);
        }
    };
}

#endif
