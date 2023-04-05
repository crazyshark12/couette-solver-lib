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
    virtual void solve();
    virtual void solveFlux(vector<Matrix> U_rho_Up, Matrix U2_Up, Matrix U3_Up, vector<Matrix> U_rho_Down, Matrix U2_Down, Matrix U3_Down);
    virtual void prepareSolving();
    void run() {solve();}
    virtual void calcRiemanPStar();
    virtual void calcFliux();

    Mixture mixture;
    macroParam upParam;
    macroParam downParam;
    solverParams solParam;

    AdditionalSolver additionalSolver;

    Matrix U_rho_velocity, U_rho_energy, pres, T;
    vector<Matrix> U_rho, R;        //U_density[i] = i-ая компонента
protected:
    Matrix  F2, F3;
    vector<Matrix> F1;
    vector <double> x;
    double delta_h;                 // шаг сетки
    Matrix timeSolvind;             // вектор в котором хранятся временные шаги

    vector<Matrix> up_density;
    vector<Matrix> down_density;
    Matrix up_velocity;
    Matrix down_velocity;
    Matrix up_pressure;
    Matrix down_pressure;
    Matrix left_Tv, right_Tv;
    Matrix Tl,Tr;
    Matrix T12L,T12R, T3L, T3R;
    vector<macroParam> rezultAfterPStart;



    void prepareVectors();
};
