#pragma once

#include "abstractsolver.h"

struct GodunovSolver: public AbstaractSolver
{

    void solve();
    void solveFlux();

protected:
    AdditionalSolver additionalSolver;
    vector<macroParam> rezultAfterPStart;

    // считает поток на границе ячеек и записывает в vector<macroParam> rezultAfterPStart
    virtual void calcRiemanPStar(vector<Matrix> density_, Matrix velocity_ , Matrix pressure_);

    // записывает в переменные F1,F2 ... значения
    virtual void calcFliux(vector<Matrix> density_, Matrix velocity_ , Matrix pressure_);
};
