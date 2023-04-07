#include "Mixture.h"

Mixture::Mixture(std::vector<MixtureComponent> components_)
{
    components = components_;
    NumberOfComponents = components.size();
}
