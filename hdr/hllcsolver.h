#pragma once

#include "abstractsolver.h"


struct HLLCSolver: public AbstaractSolver
{
    HLLCSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_);

    // запускает процесс решения задачи
    void solve();

    // устанавливает некоторые граничные условия (TODO сделать более общую структуру)
    void setBorderConditions(double up_velocity_, double up_temp_, double down_temp_);

    //устанавливает начальное распрделение температуры, плотности и скорости
    void setStartCondition(macroParam start);

    //устанавливает записыватель и поднимает флаг записи
    void setWriter(DataWriter *writer_);

    // Значения потока на границах ячеек по методу HLLC
    Matrix  hllcF2, hllcF3;
    vector<Matrix> hllcF1;
protected:
    //записывает текущие макропараметры points[] в папку с названием i
    void writePoints(double i);

    //с помощью граничных условий задаёт значения в крайних ячейках
    void useBorder();
    void UpdateBorderU();

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
    double computeT(macroParam p, size_t i);

    //записывать ли данные в файл ?
    bool isWriteData = false;
};
