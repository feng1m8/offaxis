#ifndef OFFAXLINE_RAYTRACING_HXX
#define OFFAXLINE_RAYTRACING_HXX

#include "ynogk_cxx/particle.hxx"

namespace offaxis
{
    double redshift(double radius, double a_spin, double lamda);

    class Ray
    {
    public:
        enum Tracing : int
        {
            Disk,
            Infinity = -1,
            BlackHole = -2,
        };

        Ray(double radius, double theta2rad, double phi2rad, const double *velo, double a_spin, double Rin, double Rout);

        Tracing tracing(double pr, double ptheta, double pphi);

        double redshift();

        auto operator->()
        {
            return this->ptcl.operator->();
        }

    private:
        double phi2rad;
        double Rin;
        double Rout;
        ynogk::Particle ptcl;
    };
}

#endif
