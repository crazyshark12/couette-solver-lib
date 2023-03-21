#include <functional>

#include "abstractsolver.h"

void AbstaractSolver::setMixture(Mixture mixture_)
{
    mixture =  mixture_;
}

void AbstaractSolver::prepareSolving()
{

    U_velocity.resize(solParam.NumCell+2);
    U_energy.resize(solParam.NumCell+2);
    U_density.resize(mixture.NumberOfComponents);
    for(size_t i = 0 ; i <  U_density.size(); i++)
        U_density[i].resize(solParam.NumCell+2);

//  ???  как мы устанавливаем эти параметры ???
    leftParam.density = leftParam.pressure /(UniversalGasConstant/molMass * leftParam.temp);
    leftParam.soundSpeed = sqrt(solParam.Gamma*leftParam.pressure/leftParam.density);
    leftParam.velocity = solParam.Ma*leftParam.soundSpeed;

    rightParam.density = ((solParam.Gamma + 1)* pow(solParam.Ma,2))/(2 + (solParam.Gamma -1)* pow(solParam.Ma,2))*leftParam.density;
    rightParam.pressure = (pow(solParam.Ma,2) * 2* solParam.Gamma - (solParam.Gamma - 1))/((solParam.Gamma +1))*leftParam.pressure;
    rightParam.temp = rightParam.pressure/(rightParam.density*UniversalGasConstant/molMass);
auto g = solParam.Gamma;
auto m = solParam.Ma;
    auto temp2 = (2*g*m*m/(g+1) - (g-1)/(g+1))/((g+1)* m*m/(2 + (g-1)*m*m));

    //rightParam = additionalSolver.BoundaryCondition[typeBC](leftParam, solParam);
    rightParam.velocity = leftParam.density*leftParam.velocity/rightParam.density;

    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if(i < solParam.NumCell/2+1)
        {
            U1[i] = leftParam.density;
            U2[i] = leftParam.density*leftParam.velocity;
            U3[i] = leftParam.pressure/(solParam.Gamma-1)+0.5*pow(leftParam.velocity,2)*leftParam.density;
        }

        else
        {
            U1[i] = rightParam.density;
            U2[i] = rightParam.density*rightParam.velocity;
            U3[i] = rightParam.pressure/(solParam.Gamma-1)+0.5*pow(rightParam.velocity,2)*rightParam.density;
        }
    }
    prepareVectors();
}

void AbstaractSolver::calcRiemanPStar()
{
    const std::function<void(int&)> calcPStar = [this](int& i)
    {
        macroParam left;
        macroParam right;
        mutex.lock();
        left.density = left_density[i];
        left.velocity = left_velocity[i];
        left.pressure = left_pressure[i];
        right.density = right_density[i];
        right.velocity = right_velocity[i];
        right.pressure = right_pressure[i];
        mutex.unlock();
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolverCorrect(left,right, solParam.Gamma, solParam.Gamma);
        return;
    };
}

void AbstaractSolver::prepareVectors()
{
    double x_right =solParam.lambda*solParam.lambdaSol; //% правая граница
    delta_h = (x_right) / solParam.NumCell;
    x.clear();
    x.push_back(0+0.5*delta_h);
    for(auto i = 1; i < solParam.NumCell; i++)
        x.push_back(x[i-1] + delta_h);
    x.push_back(x_right);
    x.insert(x.begin(), 0);
    U1[0]=U1[1];
    U2[0]=solParam.typeLeftBorder*U2[1];
    U3[0]=U3[1];
    U4[0]=U4[1];
    U5[0]=U5[1];
    U1[solParam.NumCell+1]=U1[solParam.NumCell];
    U2[solParam.NumCell+1]=solParam.typeRightBorder*U2[solParam.NumCell];
    U3[solParam.NumCell+1]=U3[solParam.NumCell];
    U4[solParam.NumCell+1]=U4[solParam.NumCell];
    U5[solParam.NumCell+1]=U5[solParam.NumCell];

    timeSolvind.push_back(0);
    for(int i = 0 ; i<  solParam.NumCell+1; i++)
        vectorForParallelSolving.push_back(i);
    F1.resize(solParam.NumCell+1);
    F2.resize(solParam.NumCell+1);
    F3.resize(solParam.NumCell+1);
    F4.resize(solParam.NumCell+1);
    F5.resize(solParam.NumCell+1);
    P.resize(solParam.NumCell+2);
    Q_v.resize(solParam.NumCell+2);
    Q_v3.resize(solParam.NumCell+2);
    Q_t.resize(solParam.NumCell+2);
    R.resize(solParam.NumCell+1);
    R_1.resize(solParam.NumCell+1);
    R_2.resize(solParam.NumCell+1);
    T.resize(solParam.NumCell +2);
    Tv.resize(solParam.NumCell +2);
    T12.resize(solParam.NumCell +2);
    T3.resize(solParam.NumCell +2);
    Ent.resize(solParam.NumCell +2);
    Ent2.resize(solParam.NumCell +2);
    B_v.resize(solParam.NumCell +2);
    E_Z.resize(solParam.NumCell +2);
    PR.resize(solParam.NumCell +2);
    rezultAfterPStart.resize(solParam.NumCell+1);
}
