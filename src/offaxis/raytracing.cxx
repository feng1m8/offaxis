#include <cmath>

#include "raytracing.hxx"

static double redshiftlp(ynogk::Particle &pt)
{
    return (pt->mt.expnu * pt->f1234[4]) / ((1.0 - pt->mt.somiga * pt->lambda) * offaxis::redshift(pt->radius, pt->a_spin, pt->lambda));
}

namespace offaxis
{
    double Ray::operator()(double pr, double ptheta, double pphi)
    {
        this->ptcl.lambdaq(pr, ptheta, pphi);

        this->ptcl->radius = 1.0e11;

        double pem = this->ptcl.pemdisk(0.0, this->Rout, this->Rin);
        if (pem < 0.0)
            return pem;

        this->ptcl->phi = -this->ptcl.phi(pem);

        return redshiftlp(this->ptcl);
    }

    double redshift(double radius, double a_spin, double lamda)
    {
        double omegas, expnu, exppsi, temp;
        metricgij(radius, 0.0, 1.0, a_spin, &omegas, &expnu, &exppsi, &temp, &temp);
        double omegak = 1.0 / (a_spin + radius * std::sqrt(radius));
        double vphi = exppsi / expnu * (omegak - omegas);
        return expnu * std::sqrt(1.0 - vphi * vphi) / (1.0 - omegak * lamda);
    }
}
