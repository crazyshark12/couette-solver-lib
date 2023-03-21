#pragma once

#include <vector>
#include <mutex>

#include "global.h"

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
    macroParam ExacRiemanSolverCorrect(macroParam left, macroParam right, double GammaL, double GammaR,double lambda = 0);

//    vector<vector<double>> SolveEvolutionExplFirstOrder(Matrix F1,   Matrix F2,    Matrix F3, Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old,
//                                                          double dt, double delta_h);
//    vector<vector<double>> SolveEvolutionExplFirstOrderForCO22(Matrix F1,   Matrix F2,    Matrix F3, Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old,
//                                                                 double dt, double delta_h, Matrix R = Matrix());
//    vector<vector<double>> SolveEvolutionExplFirstOrderForO2(Matrix F1,   Matrix F2,    Matrix F3, Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old,
//                                                                 double dt, double delta_h);
//    vector<vector<double>> SEEFOForCO2(Matrix F1,   Matrix F2,    Matrix F3, Matrix F4, Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old,
//                                                                 double dt, double delta_h, Matrix F11, Matrix F22, Matrix F33, Matrix F44, Matrix R);
//     vector<Matrix> SolveMUSCL_UL_UR(Matrix U1old, Matrix U2old, Matrix U3old, Matrix U4old, double Gamma, double dt_dx,vector<double>EnergyVibr, double energyStepTemp, double energyStartTemp,int LimType = 1, double omega = 0);

//     vector<vector<double>> SEEFOForCO23T(Matrix F1,      Matrix F2,      Matrix F3,      Matrix F4,      Matrix F5,
//                                            Matrix U1old,   Matrix U2old,   Matrix U3old,   Matrix U4old,   Matrix U5old,
//                                            double dt, double delta_h, Matrix R_1, Matrix R_2);

};
