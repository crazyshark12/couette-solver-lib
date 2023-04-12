#include "coeffsolver.h"
#include "global.h"
#include <cmath>

double CoeffSolver::shareViscositySimple(macroParam currentPoint)
{
    double temp1 = sqrt(kB*currentPoint.temp/(M_PI*currentPoint.mixture.components[0].mass));
    double temp2 = 2*M_PI*pow(currentPoint.mixture.components[0].sigma,2);
    double omega2 = temp1*temp2;
    return (5*kB*currentPoint.temp) /(8.*omega2);
}

double CoeffSolver::lambda(macroParam currentPoint)
{
    double temp1 = sqrt(kB*currentPoint.temp/(M_PI*currentPoint.mixture.components[0].mass));
    double temp2 = 2*M_PI*pow(currentPoint.mixture.components[0].sigma,2);
    double omega2 = temp1*temp2;
    return (75*pow(kB,2)*currentPoint.temp) /(32 * currentPoint.mixture.components[0].mass * omega2 ) ;
}
