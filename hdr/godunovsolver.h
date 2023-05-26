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

    // Расчет потоков на стыках ячеек методом годунова

    macroParam ExacRiemanSolver(macroParam left, macroParam right, double Gamma);

    // bool velocity_component: 0 - касательная, 1 - нормальная
    macroParam ExacRiemanSolver(macroParam left, macroParam right, double Gamma, bool velocity_component);

    // обновляет вектор U
    void updateU();

    // обновлеяет вектор макропараметров с помощью U
    void updatePoints();


    double shareViscositySimple(double currentT);


    vector<macroParam>rezultAfterPStart;
};
