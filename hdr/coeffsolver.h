#pragma once
#include "global.h"

struct CoeffSolver
{
public:

    double shareViscositySimple(macroParam currentPoint);
    double lambda(macroParam currentPoint);

};
