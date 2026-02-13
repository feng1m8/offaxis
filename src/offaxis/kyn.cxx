#include <CCfits/CCfits>
#include <filesystem>

#include "envs.hxx"
#include "kyn.hxx"
#include "raytracing.hxx"

namespace offaxis
{
    namespace utils
    {
        static double mod(double x1, double x2)
        {
            return x1 - x2 * std::floor(x1 / x2);
        }

        static std::valarray<double> flatten(const std::vector<std::valarray<double>> &input)
        {
            std::valarray<double> output(input.size() * input[0].size());

            for (std::size_t i = 0; i < input.size(); ++i)
                std::copy(std::begin(input[i]), std::end(input[i]), std::begin(output) + i * input[0].size());

            return output;
        }
    }

    KBHtables::KBHtables(const std::filesystem::path &path)
    {
        const CCfits::FITS fits(path.string());

        fits.extension(2).column(1).read(this->r_horizon, 0, fits.extension(2).rows());
        fits.extension(3).column(1).read(this->inclination, 0, fits.extension(3).rows());
        fits.extension(4).column(1).read(this->r_vector, 0, fits.extension(4).rows());
        fits.extension(5).column(1).read(this->phi_vector, 0, fits.extension(5).rows());

        std::size_t size = this->r_horizon.size() * this->inclination.size();

        this->alpha.reserve(size);
        this->beta.reserve(size);
        this->lensing.reserve(size);

        for (int i = 6; i < 6 + static_cast<int>(size); ++i)
        {
            std::vector<std::valarray<double>> vals;

            CCfits::Column *cols = &fits.extension(i).column("ALPHA", false);
            cols->readArrays(vals, 0, cols->rows());
            this->alpha.emplace_back(utils::flatten(vals));

            cols = &fits.extension(i).column("BETA", false);
            cols->readArrays(vals, 0, cols->rows());
            this->beta.emplace_back(utils::flatten(vals));

            cols = &fits.extension(i).column("LENSING", false);
            cols->readArrays(vals, 0, cols->rows());
            this->lensing.emplace_back(utils::flatten(vals));
        }
    }

    KBHinterp KBHtables::interp(double a_spin, double Incl) const
    {
        return KBHinterp(*this, a_spin, Incl);
    }

    KBHinterp::KBHinterp(const KBHtables &that, double a_spin, double Incl) : r_vector(that.r_vector), phi_vector(that.phi_vector)
    {
        this->a_spin = a_spin;
        this->a_spin_sq = a_spin * a_spin;

        double r_a = std::sqrt(1.0 - this->a_spin_sq);
        this->r_plus = 1.0 + r_a;
        this->r_minus = 1.0 - r_a;

        this->phi_a = a_spin / (this->r_plus - this->r_minus);

        double incl2rad = utils::deg2rad(Incl);
        this->sini = std::sin(incl2rad);
        double cosi = std::cos(incl2rad);
        this->cosi_sq = cosi * cosi;

        std::size_t i_rplus = std::lower_bound(that.r_horizon.cbegin(), that.r_horizon.cend(), this->r_plus) - that.r_horizon.cbegin();
        std::size_t i_incl = std::lower_bound(that.inclination.cbegin(), that.inclination.cend(), Incl) - that.inclination.cbegin();

        double drplus = (r_plus - that.r_horizon[i_rplus - 1]) / (that.r_horizon[i_rplus] - that.r_horizon[i_rplus - 1]);
        double dincl = (Incl - that.inclination[i_incl - 1]) / (that.inclination[i_incl] - that.inclination[i_incl - 1]);

        double weights[4] = {(1.0 - drplus) * (1.0 - dincl), (1.0 - drplus) * dincl, drplus * (1.0 - dincl), drplus * dincl};

        std::size_t index[4] = {
            (i_rplus - 1) * that.inclination.size() + (i_incl - 1),
            (i_rplus - 1) * that.inclination.size() + i_incl,
            i_rplus * that.inclination.size() + (i_incl - 1),
            i_rplus * that.inclination.size() + i_incl,
        };

        this->alpha = weights[0] * that.alpha[index[0]] + weights[1] * that.alpha[index[1]] + weights[2] * that.alpha[index[2]] + weights[3] * that.alpha[index[3]];
        this->beta = weights[0] * that.beta[index[0]] + weights[1] * that.beta[index[1]] + weights[2] * that.beta[index[2]] + weights[3] * that.beta[index[3]];
        this->lensing = weights[0] * that.lensing[index[0]] + weights[1] * that.lensing[index[1]] + weights[2] * that.lensing[index[2]] + weights[3] * that.lensing[index[3]];

        this->spline[0] = gsl_interp2d_alloc(gsl_interp2d_bilinear, this->phi_vector.size(), this->r_vector.size());
        this->spline[1] = gsl_interp2d_alloc(gsl_interp2d_bilinear, this->phi_vector.size(), this->r_vector.size());
        this->spline[2] = gsl_interp2d_alloc(gsl_interp2d_bilinear, this->phi_vector.size(), this->r_vector.size());

        gsl_interp2d_init(this->spline[0], this->phi_vector.data(), this->r_vector.data(), std::begin(this->alpha), this->phi_vector.size(), this->r_vector.size());
        gsl_interp2d_init(this->spline[1], this->phi_vector.data(), this->r_vector.data(), std::begin(this->beta), this->phi_vector.size(), this->r_vector.size());
        gsl_interp2d_init(this->spline[2], this->phi_vector.data(), this->r_vector.data(), std::begin(this->lensing), this->phi_vector.size(), this->r_vector.size());

        this->acc[0] = gsl_interp_accel_alloc();
        this->acc[1] = gsl_interp_accel_alloc();
    }

    KBHinterp::~KBHinterp()
    {
        gsl_interp2d_free(this->spline[0]);
        gsl_interp2d_free(this->spline[1]);
        gsl_interp2d_free(this->spline[2]);

        gsl_interp_accel_free(this->acc[0]);
        gsl_interp_accel_free(this->acc[1]);
    }

    std::tuple<double, double, double> KBHinterp::operator()(double radius, double phi) const
    {
        double r_kerr = radius - this->r_plus;
        double phi_kerr = this->phi_a * std::log(r_kerr / (radius - this->r_minus));
        phi_kerr = utils::mod(phi - phi_kerr + 0.5 * M_PI, 2.0 * M_PI);

        double alpha = gsl_interp2d_eval_extrap(this->spline[0], this->phi_vector.data(), this->r_vector.data(), std::begin(this->alpha), phi_kerr, r_kerr, acc[0], acc[1]);
        double beta = gsl_interp2d_eval_extrap(this->spline[1], this->phi_vector.data(), this->r_vector.data(), std::begin(this->beta), phi_kerr, r_kerr, acc[0], acc[1]);
        double lensing = gsl_interp2d_eval_extrap(this->spline[2], this->phi_vector.data(), this->r_vector.data(), std::begin(this->lensing), phi_kerr, r_kerr, acc[0], acc[1]);

        double lamda = alpha * this->sini;
        double q = std::sqrt(beta * beta + (alpha * alpha - this->a_spin_sq) * this->cosi_sq);
        double gobs = redshift(radius, this->a_spin, lamda);
        double cosem = gobs * q / radius;

        if (cosem > 1.0)
            cosem = 1.0;

        return {gobs, cosem, lensing};
    }
}
