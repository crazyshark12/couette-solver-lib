#include "Mixture.h"

Mixture::Mixture(std::vector<Component> components_)
{
    components = components_;
    NumberOfComponents = components.size();
    for(size_t i = 0; i < NumberOfComponents; i++)
    {
        macroAverage.density += components[i].density;
    }
    molarFraction = new double(NumberOfComponents);
    for(size_t i = 0; i < NumberOfComponents; i++)
    {
        molarFraction[i] = components[i].density / macroAverage.density;
    }
}
