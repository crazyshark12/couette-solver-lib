#pragma once

#include "abstractsolver.h"


struct HLLCSolver: public AbstaractSolver
{
    HLLCSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_);
    void solve();
    // Значения потока на границах ячеек по методу HLLC
    Matrix  hllcF2, hllcF3;
    vector<Matrix> hllcF1;

    void setBorderConditions(double up_velocity_, double up_temp_, double down_temp_);
    void setStartCondition(macroParam start);

protected:
    //с помощью граничных условий задаёт значения в крайних ячейках
    void useBorder();

    //подготовка размеров всех нужных векторов и инициализация начальными параметрами
    void prepareSolving();

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
