#ifndef OFFAXIS_KYN_HXX
#define OFFAXIS_KYN_HXX

#include <map>
#include <string>
#include <valarray>
#include <vector>

#define HAVE_INLINE
#include <gsl/gsl_interp2d.h>
#undef HAVE_INLINE

namespace offaxis
{
    class KBHtable
    {
        friend class Kyn;

    public:
        KBHtable(const std::string &path);

    private:
        std::vector<double> r_horizon;
        std::vector<double> inclination;
        std::vector<double> r_vector;
        std::vector<double> phi_vector;
        std::vector<std::valarray<double>> alpha;
        std::vector<std::valarray<double>> beta;
        std::vector<std::valarray<double>> lensing;
    };

    class Kyn
    {
    public:
        Kyn(const KBHtable &that, double a_spin, double Incl);
        ~Kyn();

        std::tuple<double, double, double> interpolate(double radius, double phi) const;

    private:
        const std::vector<double> &r_vector;
        const std::vector<double> &phi_vector;

        double a_spin;
        double a_spin_sq;
        double r_plus;
        double r_minus;
        double phi_a;
        double sini;
        double cosi_sq;

        std::valarray<double> alpha;
        std::valarray<double> beta;
        std::valarray<double> lensing;
        gsl_interp_accel *acc[2];
        gsl_interp2d *spline[3];
    };

    namespace envs
    {
        inline std::map<int, const KBHtable> table;
    }
}

#endif
