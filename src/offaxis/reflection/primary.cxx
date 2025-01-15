#include <memory>

#include <omp.h>

#include "offaxis/parameter.hxx"
#include "ynogk_cxx/particle.hxx"

#include "cache.hxx"
#include "offaxis/envs.hxx"
#include "spectrum.hxx"
#include "offaxis/sphere.hxx"

namespace offaxis::offaxxillver
{
    double doppler(double coslp, double sinlp, double philp, double *vlp, double incl)
    {
        double phi2rad = utils::deg2rad(philp);
        double sinphi = std::sin(phi2rad);
        double cosphi = std::cos(phi2rad);

        std::valarray<double> X{sinlp * cosphi, coslp * cosphi, -sinphi};
        std::valarray<double> Z{coslp, -sinlp, 0.0};

        double incl2rad = utils::deg2rad(incl);
        std::valarray<double> obs{std::sin(incl2rad), 0.0, std::cos(incl2rad)};

        std::valarray<double> k(obs[0] * X + obs[2] * Z);

        return std::sqrt(1.0 - vlp[0] * vlp[0] - vlp[1] * vlp[1] - vlp[2] * vlp[2]) / (1.0 - vlp[0] * k[0] - vlp[1] * k[1] - vlp[2] * k[2]);
    }

    double redshiftinf(double rlp, double coslp, double sinlp, double a_spin)
    {
        double expnu, temp;
        metricgij(rlp, coslp, sinlp, a_spin, &temp, &expnu, &temp, &temp, &temp);
        return expnu;
    }

    static double to_infinity(double a_spin, double rlp, double theta, int nside)
    {
        const Sphere &sphere(envs::sphere.try_emplace(nside, nside).first->second);

        double theta2rad = utils::deg2rad(theta);
        double coslp = std::cos(theta2rad);
        double sinlp = std::sin(theta2rad);
        double vlp[3] = {0.0, 0.0, 0.0};

        ynogk::Particle ptcl(a_spin, rlp, coslp, sinlp, 1.0, vlp);

        std::size_t to_inf = 0;

#pragma omp parallel for firstprivate(ptcl) reduction(+ : to_inf)
        for (std::size_t pix = 0; pix < sphere.size; ++pix)
        {
            auto [pr, ptheta, pphi] = sphere[pix];
            ptcl.lambdaq(pr, ptheta, pphi);
            double pem = ptcl.pemdisk(0.0, 1000.0, 0.0);
            if (pem == -1.0)
                ++to_inf;
        }

        return static_cast<double>(to_inf) / static_cast<double>(sphere.size);
    }

    std::valarray<double> primary(const std::valarray<double> &energy, const std::valarray<double> &parameter, T_PrimSpec prim_type)
    {
        using namespace parameter::offaxxillCp;

        double theta2rad = utils::deg2rad(parameter[thetalp]);
        double coslp = std::cos(theta2rad);
        double sinlp = std::sin(theta2rad);
        double vlp[3] = {parameter[vr], parameter[vtheta], parameter[vphi]};

        double dinf = doppler(coslp, sinlp, parameter[philp], vlp, parameter[Incl]);
        double ginf = dinf * redshiftinf(parameter[rlp], coslp, sinlp, parameter[a_spin]);

        const auto xillver(std::make_unique<const relxill::Spectrum>(parameter, prim_type));

        static Cache to_inf(to_infinity, 2);

        double f_prim = to_inf(parameter[a_spin], parameter[rlp], parameter[thetalp], envs::nside() / 2);
        double beaming = dinf * dinf * std::pow(ginf, parameter[gamma]);

        return f_prim * beaming * xillver->nthcomp(energy, ginf);
    }
}
