#pragma once

#include "abstractsolver.h"

struct HLLCSolver: public AbstaractSolver
{
    void solve();
    // Значения потока на границах ячеек по методу HLLC
    Matrix  hllcF2, hllcF3;
    vector<Matrix> hllcF1;

    vector<macroParam> points;
protected:

    // Расчет вектора потоков во всех ячейках
    void computeF();

    // Расчет релаксационных членов
    void computeR();

    // Расчет потоков на стыках ячеек методом HLLC
    void computeHllcF();

    // обновляет вектор U
    void updateU();

    // обновлеяет вектор макропараметров с помощью U
    void updatePoints();

    //вычисляет температуру в i-ой ячейке
    void computeT(macroParam &p, size_t i);
};
