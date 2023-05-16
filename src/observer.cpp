#include "observer.h"
#include <cmath>
#include <algorithm>
#include <iostream>

void Observer::setPeriodicity(int period)
{
    if(period > 0)
        periodicity = period;
    else
        cout<<" wrong period "<<endl;
}

void Observer::remember(vector<macroParam> param)
{
    previousParam = param;
}

bool Observer::checkDifference(vector<macroParam> param)
{
    double maxDifference = 0;
    auto size = param.size();
    for(size_t i = 0; i < size; i++)
    {
        double d_temp = std::abs(previousParam[i].temp - param[i].temp);
        double d_rho = std::abs(previousParam[i].density - param[i].density);
        double d_pressure = std::abs(previousParam[i].pressure - param[i].pressure);
        double d_velocity = std::abs(previousParam[i].velocity - param[i].velocity);
        double d_densityArray = 0;
        double d_fractionArray = 0;
        for(size_t j = 0; j < param[i].mixture.NumberOfComponents; j++)
        {
            double d_da = previousParam[i].densityArray[j] - param[i].densityArray[j];
            if(d_densityArray < d_da)
                d_densityArray = d_da;

            double d_fa = previousParam[i].fractionArray[j] - param[i].fractionArray[j];
            if(d_fractionArray < d_fa)
                d_fractionArray = d_fa;
        }
        double localMaxDifference = std::max({d_temp,d_rho,d_pressure, d_velocity,d_densityArray,d_fractionArray});
        if(maxDifference < localMaxDifference)
            maxDifference = localMaxDifference;
    }
    if(maxDifference > precision)
        return true;
    else
        return false;
}

int Observer::getPeriodicity()
{
    return periodicity;
}
