#pragma once
#include "global.h"
#include <vector>
#include <string>

struct Component
{
    double density;                 //плотность компоненты
    //... какие-то другие параметры компонент
    string name;                    //название компоненты

    double getEntalp(double T); // пример получения некоторых параметров
};
struct Mixture
{
    Mixture(std::vector<Component> components_);

    std::vector<Component> components;
    macroParam macroAverage;            //макропараметры всей смеси в целом
    int NumberOfComponents;
    double *molarFraction;              //молярная доля компонент

    //... какие-то другие параметры смеси
    //пример
    double **coeffDiffusion;            //коэффициенты диффузии
};
