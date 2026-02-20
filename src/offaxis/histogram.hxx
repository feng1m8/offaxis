#ifndef OFFAXIS_HISTOGRAM_HXX
#define OFFAXIS_HISTOGRAM_HXX

#include <valarray>
#include <vector>

#include <gsl/gsl_histogram.h>

namespace offaxis
{
    class Histogram
    {
    private:
        gsl_histogram *hist;

    public:
        // Histogram(const std::valarray<double> &range)
        // {
        //     this->hist = gsl_histogram_alloc(range.size() - 1);
        //     gsl_histogram_set_ranges(this->hist, std::begin(range), range.size());
        // }

        Histogram(const std::vector<double> &range)
        {
            this->hist = gsl_histogram_alloc(range.size() - 1);
            gsl_histogram_set_ranges(this->hist, range.data(), range.size());
        }

        Histogram(const Histogram &other)
        {
            this->hist = gsl_histogram_clone(other.hist);
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

        void operator+=(const Histogram &rhs)
        {
            gsl_histogram_add(this->hist, rhs.hist);
        }
    };
}

#endif
