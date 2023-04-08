#include "godunovsolver.h"


void GodunovSolver::solve()
{
    prepareSolving();
    for(auto i  = 0; i < solParam.MaxIter; i++)
    {
        setDt();

        //считает вектор потоков
        solveFlux();

        //с помощью потоков устанавливает значение значение U_new, я эту функцию поменял, но она мне не нравится ибо возвращаемый тип неочевидный
        auto res = additionalSolver.SolveEvolutionExplFirstOrder(F1, F2,F3,U1, U2,U3,timeSolvind.last(),delta_h);

        U1 = res.first;
        U2 = res.second[0];
        U3 = res.second[1];

        //тут происходит что-то странное, видимо своего рода граничные условия
        for(size_t j = 0; j < U1.size();j++)
        {
            U1[j][0] = U1[j][1];
            size_t sz = U1[j].size();
            U1[j][sz-1] = U1[j][sz-2];
        }
        U2[0]=U2[1];
        U2[U2.size()-1] = U2[U2.size()-2];
        U3[0]=U3[1];
        U3[U3.size() - 1] =U3[U3.size() - 2];
    }
}

void GodunovSolver::solveFlux()
{
    auto density   = U1;
    auto velocity  = U2 / U1[0];
    auto pressure  = (U3 - Matrix::POW(velocity, 2)*0.5*density[0])*(solParam.Gamma-1);
    calcRiemanPStar(density,velocity,pressure);
    calcFliux(density,velocity,pressure);
}
void GodunovSolver::calcRiemanPStar(vector<Matrix> density_, Matrix velocity_ , Matrix pressure_)
{
    for(size_t i = 0; i < velocity_.size()-1; i++)
    {
        macroParam up(mixture);
        macroParam down(mixture);
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        {
            up.densityArray[j] = density_[j][i+1];
            down.densityArray[j] = density_[j][i];
        }
        up.velocity = velocity_[i+1];
        up.pressure = pressure_[i+1];
        down.velocity = velocity_[i];
        down.pressure = pressure_[i];
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolver(up,down, solParam.Gamma);
    }
    return;
}
void GodunovSolver::calcFliux(vector<Matrix> density_, Matrix velocity_ , Matrix pressure_)
{
    for(size_t i = 0; i < velocity_.size()-1; i++)
    {
        auto point = rezultAfterPStart[i];
        double tempU = pressure_[i+1]/(density_[0][i+1]*UniversalGasConstant/molMass);
        double tempD = pressure_[i]/(density_[0][i]*UniversalGasConstant/molMass);
        auto du_dx = (velocity_[i] - velocity_[i+1])/delta_h;
        double Tx = point.pressure/(point.density*UniversalGasConstant/molMass);
        //T[i] = Tx;
        double Pr = 2.0/3;
        double etta /*= additionalSolver.shareViscosity[0](upParam.temp, Tx,0,0);*/; //тут надо добавить рассчёт вязкости
        double G =  (4.0/3*etta)*du_dx;
        double k = solParam.Gamma*UniversalGasConstant/molMass*etta/(solParam.Gamma-1)/Pr;

        double dt_dx = (tempD - tempU)/delta_h;
        double q =   -k*dt_dx;

        //Тут слагаемые в целом нужно поправить под нашу задачу
        F1[0][i] = 0;
        for(size_t j = 0; j< F1.size(); j++)
            F1[j][i] = (point.density * point.velocity); // тут должна быть диффузия
        F2[i] = (/*F1.last()**/point.velocity + point.pressure -G);
        F3[i] = ((point.pressure/(solParam.Gamma - 1) + point.density* pow(point.velocity,2)/2 + point.pressure) * point.velocity -G*point.velocity + q);
    }
    return;
}
