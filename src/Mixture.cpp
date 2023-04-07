#include "Mixture.h"

Mixture::Mixture(std::vector<MixtureComponent> components_)
{
    components = components_;
    NumberOfComponents = components.size();
}

double Mixture::molarMass()
{
    double sum = 0;
    for(size_t i = 0; i < NumberOfComponents; i++)
        sum += components[i].molarMass;
    return sum;
}

double Mixture::molarMass(size_t i)
{
    return components[i].molarMass;
}
