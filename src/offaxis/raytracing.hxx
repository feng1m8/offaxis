#ifndef OFFAXIS_RAYTRACING_HXX
#define OFFAXIS_RAYTRACING_HXX

#include "ynogk_cxx/particle.hxx"

namespace offaxis
{
    class Ray
    {
    public:
        Ray(double a_spin, double rlp, double coslp, double sinlp, double *vlp, double Rin, double Rout) : ptcl(a_spin, rlp, coslp, sinlp, 1.0, vlp)
        {
            this->Rin = Rin;
            this->Rout = Rout;
            this->ptcl.metric();
        }

        double operator()(double pr, double ptheta, double pphi);

        auto operator->()
        {
            return this->ptcl.operator->();
        }

    private:
        double Rin;
        double Rout;
        ynogk::Particle ptcl;
    };

    double redshift(double radius, double a_spin, double lamda);
}

#endif
