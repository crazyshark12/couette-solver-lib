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


    // можно сказать что это начальные данные + ещё нужно как-то описать граничные условия
    Mixture mixture;
    macroParam upParam;
    macroParam downParam;
    solverParams solParam;

protected:

    virtual void prepareSolving();
    void prepareVectors();

    // запускает рассчет вектора потоков
    virtual void solveFlux() = 0;

    // обновляет вектор U
    virtual void updateU();

    // устанавливает временной шаг
    void setDt();

    //надо придумать название получше
    Matrix U2, U3, pres, T;
    vector<Matrix> U1, R;        //U_density[i] = i-ая компонента

    //надо придумать название получше
    Matrix  F2, F3;
    vector<Matrix> F1;

    double delta_h;                 // шаг сетки
    Matrix timeSolvind;             // вектор в котором хранятся временные шаги



};
