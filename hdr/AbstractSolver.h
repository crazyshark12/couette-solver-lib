#pragma once

#include <vector>
#include <list>
#include <mutex>

#include "global.h"
#include "additionalsolver.h"

struct AbstaractSolver
{
public:
    virtual void solve() = 0;
    virtual void prepareSolving();
    void run() {solve();}
    virtual void calcRiemanPStar();
    bool breaksolve = false;
    bool pauseSolve = false;
    macroParam leftParam;
    macroParam rightParam;
    solverParams solParam;
    int typeBC = 0;
    int typeTimeStepSolution = 0;
    int typeBulkVisc = 0;
    int typeShareVisc = 0;
    int typeEnergy = 0;
    AdditionalSolver additionalSolver;

    Matrix R, P, Q_v, Q_t, T, Tv, Ent, Ent2, R_1, R_2, T12, T3, Q_v3, B_v, E_Z, PR;
    Matrix U1, U2, U3, U4,U5, pres;
protected:
    Matrix F1, F2, F3, F4, F5;
    vector <double> x;
    vector<int> vectorForParallelSolving;
    double delta_h;
    Matrix timeSolvind;

    Matrix left_density;
    Matrix right_density;
    Matrix left_velocity;
    Matrix right_velocity;
    Matrix left_pressure;
    Matrix right_pressure;
    Matrix left_Tv, right_Tv;
    Matrix Tl,Tr;
    Matrix T12L,T12R, T3L, T3R;
    vector<macroParam> rezultAfterPStart;
    mutex mutex;
    list<double> CvibrMass;
    double CVibrStartTemp;
    double CVibrStepTemp;
    vector<double> EnergyVibr;
    double energyVibrStartTemp;
    double energyVibrStepTemp;
    vector<double> EnergyVibr12;
    double energyVibrStartTemp12;
    double energyVibrStepTemp12;
    vector<double> EnergyVibr3;
    double energyVibrStartTemp3;
    double energyVibrStepTemp3;

    vector<double> Energy;
    double energyStartTemp;
    double energyStepTemp;


    void prepareVectors();
};
