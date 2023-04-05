#pragma once
#include <vector>
#include <string>

struct Component
{
    double density;                 //плотность компоненты
    //... какие-то другие параметры компонент
    std::string name;                    //название компоненты

    double getEntalp(double T); // пример получения некоторых параметров
};
struct Mixture
{
    Mixture(std::vector<Component> components_);

    std::vector<Component> components;
    int NumberOfComponents;
    //... какие-то другие параметры смеси
    //пример
    double **coeffDiffusion;            //коэффициенты диффузии
};
