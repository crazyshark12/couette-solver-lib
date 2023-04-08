#pragma once

#include <vector>
#include <list>
#include <mutex>

#include "global.h"
#include "Mixture.h"

struct AbstaractSolver
{
public:
    void setMixture(Mixture mixture_);
    virtual void solve() = 0;

    // можно сказать что это начальные данные + ещё нужно как-то описать граничные условия
    Mixture mixture;
    macroParam startParam;
    solverParams solParam;
protected:

    virtual void prepareSolving();
    void prepareVectors();

//    // обновляет вектор U
//    virtual void updateU();

    // устанавливает временной шаг
    void setDt();

    virtual void updatePoints();

    //надо придумать название получше
    Matrix U2, U3;
    vector<Matrix> U1, R;        //U1[i] = i-ая компонента

    //надо придумать название получше
    Matrix  F2, F3;
    vector<Matrix> F1;

    // хранит значение макропараметров в каждой ячейке
    vector<macroParam> points;

//    //скорость звука в каждой ячейке
//    Matrix sound_speed;

    double delta_h;                 // шаг сетки
    Matrix timeSolvind;             // вектор в котором хранятся временные шаги
};
