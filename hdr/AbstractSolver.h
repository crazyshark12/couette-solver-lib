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

    virtual void prepareSolving();
    void prepareVectors();

    virtual void solve();

    // запускает рассчет вектора потоков
    virtual void solveFlux();

    // считает поток на границе ячеек и записывает в vector<macroParam> rezultAfterPStart
    virtual void calcRiemanPStar(vector<Matrix> density_, Matrix velocity_ , Matrix pressure_);

    // записывает в переменные F1,F2 ... значения
    virtual void calcFliux(vector<Matrix> density_, Matrix velocity_ , Matrix pressure_);

    // можно скахать что это начальные данные + ещё нужно как-то описать граничные условия
    Mixture mixture;
    macroParam upParam;
    macroParam downParam;
    solverParams solParam;

protected:
    AdditionalSolver additionalSolver;

    //надо придумать название получше
    Matrix U_rho_velocity, U_rho_energy, pres, T;
    vector<Matrix> U_rho, R;        //U_density[i] = i-ая компонента

    //надо придумать название получше
    Matrix  F2, F3;
    vector<Matrix> F1;

    double delta_h;                 // шаг сетки
    Matrix timeSolvind;             // вектор в котором хранятся временные шаги

    vector<macroParam> rezultAfterPStart;


};
