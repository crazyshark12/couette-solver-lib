#pragma once

#include <vector>
#include <list>
#include <mutex>

#include "global.h"
#include "mixture.h"
#include "coeffsolver.h"
#include "bordercondition.h"
#include "datawriter.h"

struct AbstaractSolver
{
public:
    AbstaractSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_);

    // запускает процесс решения задачи
    virtual void solve() = 0;

    //устанавливает размер ячейки
    void setDelta_h(double dh);

    // устанавливает некоторые граничные условия (TODO сделать более общую структуру)
    void setBorderConditions(double up_velocity_, double up_temp_, double down_temp_);

    // устанавливает начальное распрделение температуры, плотности и скорости
    void setStartCondition(macroParam start);

    // устанавливает записыватель и поднимает флаг записи
    void setWriter(DataWriter *writer_);

    Mixture mixture;
    macroParam startParam;
    solverParams solParam;
    CoeffSolver coeffSolver;
    BorderCondition border;
    DataWriter* writer;

protected:

    //записывает текущие макропараметры points[] в папку с названием i
    void writePoints(double i);

    //заполняет начальные ячеки
    virtual void prepareSolving();

    //подготавливает размеры всех векторов
    void prepareVectors();

    // устанавливает временной шаг
    void setDt();

    virtual void updatePoints();

    // Значения потока на границах ячеек по методу HLLC
    Matrix  hllcF2, hllcF3;
    vector<Matrix> hllcF1;

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

    // шаг сетки
    double delta_h = 0;

    // вектор в котором хранятся временные шаги
    Matrix timeSolvind;

    //записывать ли данные в файл ?
    bool isWriteData = false;
};
