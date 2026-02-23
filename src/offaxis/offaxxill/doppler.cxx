#include <tuple>
#include <valarray>

#include <omp.h>

#include "ynogk_cxx/particle.hxx"

#include "offaxis/parameter.hxx"

#include "offaxline/envs.hxx"
#include "offaxline/memory.hxx"
#include "offaxline/sphere.hxx"

namespace offaxis
{
    double doppler(double coslp, double sinlp, double philp, const double *vlp, double incl)
    {
        double phi2rad = utils::deg2rad(philp);
        double sinphi = std::sin(phi2rad);
        double cosphi = std::cos(phi2rad);

        std::valarray X{sinlp * cosphi, coslp * cosphi, -sinphi};
        std::valarray Z{coslp, -sinlp, 0.0};

        double incl2rad = utils::deg2rad(incl);
        std::valarray obs{std::sin(incl2rad), 0.0, std::cos(incl2rad)};

        std::valarray k(obs[0] * X + obs[2] * Z);

        double v_sq = vlp[0] * vlp[0] + vlp[1] * vlp[1] + vlp[2] * vlp[2];
        double v_cos = vlp[0] * k[0] + vlp[1] * k[1] + vlp[2] * k[2];

        return std::sqrt(1.0 - v_sq) / (1.0 - v_cos);
    }

    double expnu(double rlp, double coslp, double sinlp, double a_spin)
    {
        double expnu, ignore;
        metricgij(rlp, coslp, sinlp, a_spin, &ignore, &expnu, &ignore, &ignore, &ignore);
        return expnu;
    }

    static double to_infinity(double a_spin, double rlp, double theta, long nside)
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
            auto palpha = sphere[pix];
            ptcl.lambdaq(palpha[0], palpha[1], palpha[2]);
            double pem = ptcl.pemdisk(0.0, 1000.0, 0.0);
            if (pem == -1.0)
                ++to_inf;
        }

        return static_cast<double>(to_inf) / static_cast<double>(sphere.size);
    }

    std::tuple<double, double> beaming(const std::vector<double> &parameter)
    {
        using namespace parameter::offaxxillCp;

        double theta2rad = utils::deg2rad(parameter[thetalp]);
        double coslp = std::cos(theta2rad);
        double sinlp = std::sin(theta2rad);

        double dinf = doppler(coslp, sinlp, parameter[philp], &parameter[vr], parameter[Incl]);
        double ginf = dinf * expnu(parameter[rlp], coslp, sinlp, parameter[a_spin]);

        static Memory to_inf(to_infinity);
        to_inf.max_size = envs::cache_size();

        double f_prim = to_inf(parameter[a_spin], parameter[rlp], parameter[thetalp], envs::nside() / 2);
        double beaming = dinf * dinf * std::pow(ginf, parameter[gamma]);

        return {ginf, f_prim * beaming};
    }
}
