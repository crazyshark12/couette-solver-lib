#pragma once

#include "abstractsolver.h"
struct GodunovSolver: public AbstractSolver
{
    GodunovSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_):AbstractSolver(mixture_,startParam_,solParam_){};

    // запускает процесс решения задачи
    void solve();

protected:

    void prepareVectors();

    void computeFluxF();

    // Расчет релаксационных членов
    void computeR();

    // Расчет потоков на стыках ячеек методом HLLC
    macroParam ExacRiemanSolver(macroParam left, macroParam right, double Gamma);

    // обновляет вектор U
    void updateU();

    // обновлеяет вектор макропараметров с помощью U
    void updatePoints();

    double shareViscositySimple(double currentT);


    vector<macroParam>rezultAfterPStart;
};
