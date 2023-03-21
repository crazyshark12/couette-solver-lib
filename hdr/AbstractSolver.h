#pragma once

#include <vector>
#include <list>
#include <mutex>

#include "global.h"
#include "additionalsolver.h"
#include "Mixture.h"

struct AbstaractSolver
{
public:
    void setMixture(Mixture mixture_);
    virtual void solve() = 0;
    virtual void prepareSolving();
    void run() {solve();}
    virtual void calcRiemanPStar();

    Mixture mixture;
    macroParam leftParam;
    macroParam rightParam;
    solverParams solParam;

    AdditionalSolver additionalSolver;

    Matrix U_velocity, U_energy, pres;
    vector<Matrix> U_density, R;        //U_density[i] = i-ая компонента
protected:
    Matrix F1, F2, F3;
    vector <double> x;
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



    void prepareVectors();
};
