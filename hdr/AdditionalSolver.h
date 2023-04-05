#pragma once

#include <vector>
#include <mutex>

#include "global.h"
#include <utility>

typedef std::pair<vector<Matrix>, vector<Matrix>> conservative;
struct AdditionalSolver
{
public:

    void solve();
    double startValue;
    double stopValue;
    double step;
    double pressure = 0;
    double density = 0;
    double gamma = 0;
    vector <double> iterationVector;
    vector <double> rezultVector;

    macroParam ExacRiemanSolver(macroParam left, macroParam right, double Gamma);
    macroParam ExacRiemanSolverCorrect(macroParam left, macroParam right, double GammaL, double GammaR, double lambda);
    conservative SolveEvolutionExplFirstOrder(vector<Matrix> F1, Matrix F2, Matrix F3,  vector<Matrix> U1old, Matrix U2old, Matrix U3old, double dt, double delta_h);

};
