#pragma once

#include "global.h"
#include "systemofequation.h"
enum RiemannSolverType
{
    HLLCSolver,         // not correct
    HLLESolver,         // correct
    HLLSimple,          // not correct
    HLLIsentropic,      // correct
    ExacRiemanSolver,   // correct
    HLLESolverSimen     // correct
};
struct RiemannSolver
{
    RiemannSolver(){};
    virtual void computeFlux(SystemOfEquation *system){};
    virtual void computeFlux(SystemOfEquation *system, double dt, double dh){};
};

struct HLLCSolver : public RiemannSolver
{
    HLLCSolver(){};
    void computeFlux(SystemOfEquation *system);
};


struct HLLESolver : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system);
};

struct HLLESolverSimen : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system);
};

struct HLLSimple : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system, double dt, double dh);
};

struct HLLIsentropic : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system);
};

struct ExacRiemanSolver : public RiemannSolver
{
    void computeFlux(SystemOfEquation *system);

private:
     macroParam exacRiemanSolver(macroParam left, macroParam right, double Gamma);
};
