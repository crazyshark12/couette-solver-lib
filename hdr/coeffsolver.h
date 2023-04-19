#pragma once
#include "global.h"

struct CoeffSolver
{
public:

    double shareViscositySimple(macroParam currentPoint);
    double shareViscositySimple(macroParam currentPoint, double temperature);
    double lambda(macroParam currentPoint);
    double lambda(macroParam currentPoint, double temperature);

};
