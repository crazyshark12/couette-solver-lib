#pragma once

#include "global.h"
#include "systemofequation.h"
struct RiemannSolver
{
    virtual void computeFlux(SystemOfEquation *system) = 0;
    virtual void computeFlux(SystemOfEquation *system, double dt, double dh){};
};

struct HLLCSolver : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system);
};


struct HLLESolver : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system);
};

struct HLLSimple : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system){};
    void computeFlux(SystemOfEquation *system, double dt, double dh);
};
