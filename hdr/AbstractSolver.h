#pragma once

#include <vector>
#include <list>
#include <mutex>

#include "global.h"
#include "mixture.h"
#include "coeffsolver.h"
#include "bordercondition.h"
#include "datawriter.h"
#include "observer.h"
#include "systemofequation.h"
#include "riemannsolver.h"
struct AbstractSolver
{
public:
    AbstractSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_, SystemOfEquationType type, RiemannSolverType riemannType);

    // запускает процесс решения задачи
    virtual void solve() = 0;

    //устанавливает размер ячейки
    void setDelta_h(double dh);

    // устанавливает начальное распрделение температуры, плотности и скорости
    void setStartCondition(macroParam start);

    // устанавливает некоторые граничные условия (TODO сделать более общую структуру)
    void setBorderConditions(double up_velocity_, double up_temp_, double down_temp_);

    // устанавливает записыватель и поднимает флаг записи
    void setWriter(DataWriter *writer_);

    // устанавливает наблюдателя, который будет следить за тем, чтобы рассчёт остановился, когда будет достигнута новая точность
    void setObserver(Observer* obs);

    SystemOfEquation* system;

    RiemannSolver* riemannSolver;

    Observer* observer;

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
    virtual void prepareVectorSizes();

    // устанавливает временной шаг
    void setDt();

    void updatePoints();

    void useBorder();
    void UpdateBorderU();

    //вычисляет температуру в i-ой ячейке
    double computeT(macroParam p, size_t i);

    //функция для работы с наблюдателем, выдаёт true если нужно продолжать рассчёт
    bool observerCheck(size_t currentIteration);


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

    //смотрит ли наблюдатель за установлением равновесия?
    bool isObserverWatching = false;

    // через какое количество итераций наблюдатель должен проверять соотвествие
    //(берутся две соседних итерации, а не те, что разделены watcherIteration итерациями)
    int watcherIteration = 1;
};
