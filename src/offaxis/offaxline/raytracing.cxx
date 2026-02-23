#include <cmath>

#include "offaxline/raytracing.hxx"

namespace offaxis
{
    double redshift(double radius, double a_spin, double lamda)
    {
        double omegas, expnu, exppsi, ignore;
        metricgij(radius, 0.0, 1.0, a_spin, &omegas, &expnu, &exppsi, &ignore, &ignore);
        double omegak = 1.0 / (a_spin + radius * std::sqrt(radius));
        double vphi = exppsi / expnu * (omegak - omegas);
        return expnu * std::sqrt(1.0 - vphi * vphi) / (1.0 - omegak * lamda);
    }

    Ray::Ray(double radius, double theta2rad, double phi2rad, const double *velo, double a_spin, double Rin, double Rout)
        : ptcl(a_spin, radius, std::cos(theta2rad), std::sin(theta2rad), 1.0, velo)
    {
        this->phi2rad = phi2rad;
        this->Rin = Rin;
        this->Rout = Rout;

        this->ptcl.metric();
    }

    Ray::Tracing Ray::tracing(double pr, double ptheta, double pphi)
    {
        this->ptcl.lambdaq(pr, ptheta, pphi);
        this->ptcl->radius = 1.0e11;

        double pem = this->ptcl.pemdisk(0.0, this->Rout, this->Rin);
        if (pem < 0.0)
            return Ray::Tracing(pem);

        this->ptcl->phi = this->phi2rad - this->ptcl.phi(pem);
        return Ray::Disk;
    }

    double Ray::redshift()
    {
        return (this->ptcl->mt.expnu * this->ptcl->f1234[4]) /
               ((1.0 - this->ptcl->mt.somiga * this->ptcl->lambda) *
                offaxis::redshift(this->ptcl->radius, this->ptcl->a_spin, this->ptcl->lambda));
    };
}
