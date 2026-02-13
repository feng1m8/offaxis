#ifndef OFFAXIS_KYN_HXX
#define OFFAXIS_KYN_HXX

#include <filesystem>
#include <map>
#include <valarray>
#include <vector>

#define HAVE_INLINE
#include <gsl/gsl_interp2d.h>
#undef HAVE_INLINE

namespace offaxis
{
    class KBHinterp;

    class KBHtables
    {
        friend class KBHinterp;

    public:
        KBHtables(const std::filesystem::path &path);

        KBHinterp interp(double a_spin, double Incl) const;

    private:
        std::vector<double> r_horizon;
        std::vector<double> inclination;
        std::vector<double> r_vector;
        std::vector<double> phi_vector;
        std::vector<std::valarray<double>> alpha;
        std::vector<std::valarray<double>> beta;
        std::vector<std::valarray<double>> lensing;
    };

    class KBHinterp
    {
    public:
        KBHinterp(const KBHtables &that, double a_spin, double Incl);
        ~KBHinterp();

        std::tuple<double, double, double> operator()(double radius, double phi) const;

        KBHinterp(const KBHinterp &) = delete;
        KBHinterp(KBHinterp &) = delete;

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
        inline std::map<int, const KBHtables> table;
        inline std::map<std::filesystem::path, const KBHtables> kbhtables;
    }
}

#endif
