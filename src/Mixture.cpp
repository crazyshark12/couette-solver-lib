#include "Mixture.h"

Mixture::Mixture(std::vector<Component> components_)
{
    components = components_;
    NumberOfComponents = components.size();
}
