#pragma once
#include <vector>
#include <string>

struct MixtureComponent
{
    double density;                 //плотность компоненты
    //... какие-то другие параметры компонент
    std::string name;                    //название компоненты

    double getEntalp(double T); // пример получения некоторых параметров
};
struct Mixture
{
    Mixture(std::vector<MixtureComponent> components_);

    std::vector<MixtureComponent> components;
    int NumberOfComponents;
    //... какие-то другие параметры смеси
    //пример
    double **coeffDiffusion;            //коэффициенты диффузии
};
