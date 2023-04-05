#include <functional>

#include "abstractsolver.h"

void AbstaractSolver::setMixture(Mixture mixture_)
{
    mixture =  mixture_;
}
void AbstaractSolver::prepareSolving()
{

    U_rho_velocity.resize(solParam.NumCell+2);
    U_rho_energy.resize(solParam.NumCell+2);
    U_rho.resize(mixture.NumberOfComponents);
    for(size_t i = 0 ; i <  U_rho.size(); i++)
        U_rho[i].resize(solParam.NumCell+2);

//  ???  как мы устанавливаем эти параметры ??? (взял из argonSolver)
    upParam.density = upParam.pressure /(UniversalGasConstant/molMass * upParam.temp);
    upParam.soundSpeed = sqrt(solParam.Gamma*upParam.pressure/upParam.density); // например когда и где мы задаём upParam.pressure что бы его здесь использовать
    upParam.velocity = solParam.Ma*upParam.soundSpeed;

    downParam.density = ((solParam.Gamma + 1)* pow(solParam.Ma,2))/(2 + (solParam.Gamma -1)* pow(solParam.Ma,2))*upParam.density; // так ли оно должно задаваться?
    downParam.pressure = (pow(solParam.Ma,2) * 2* solParam.Gamma - (solParam.Gamma - 1))/((solParam.Gamma +1))*upParam.pressure;
    downParam.temp = downParam.pressure/(downParam.density*UniversalGasConstant/molMass);
    auto g = solParam.Gamma;
    auto m = solParam.Ma;
    auto temp2 = (2*g*m*m/(g+1) - (g-1)/(g+1))/((g+1)* m*m/(2 + (g-1)*m*m));

    //rightParam = additionalSolver.BoundaryCondition[typeBC](leftParam, solParam);
    downParam.velocity = downParam.density*downParam.velocity/upParam.density;

    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        if(i < solParam.NumCell/2+1)
        {
            U_rho[0][i] = upParam.density;
            for(size_t j = 1; j < mixture.NumberOfComponents; j++)
                U_rho[j][i] = upParam.density * upParam.massFraction[j];
            U_rho_velocity[i] = upParam.density*upParam.velocity;
            U_rho_energy[i] = upParam.pressure/(solParam.Gamma-1)+0.5*pow(upParam.velocity,2)*upParam.density; // скорее всего иначе
        }

        else
        {
            U_rho[0][i] = downParam.density;
            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
                U_rho[j][i] = downParam.density * downParam.massFraction[j];
            U_rho_velocity[i] = downParam.density*downParam.velocity;
            U_rho_energy[i] = downParam.pressure/(solParam.Gamma-1)+0.5*pow(downParam.velocity,2)*downParam.density; // скорее всего иначе
        }
    }
    prepareVectors();
}

void AbstaractSolver::solve()
{
    prepareSolving();
    for(auto i  = 0; i < solParam.MaxIter; i++)
    {

        Matrix velocity = U_rho_velocity/U_rho[0];
        auto pressure = (U_rho_energy - Matrix::POW(velocity,2)*0.5*U_rho[0])*(solParam.Gamma - 1);
        Matrix sound_speed = Matrix::SQRT(pressure*solParam.Gamma/U_rho[0]);
        auto temp = velocity + sound_speed;                     // что это, типо абсолютная скорость?
        auto max =*std::max_element(temp.begin(), temp.end());  // видимо это нужно чтобы правильно подобрать временной шаг, чтобы соблюдался критерий КФЛ
        double dt = solParam.CFL*pow(delta_h,1)/max;

        //double dt = additionalSolver.TimeStepSolution[typeTimeStepSolution](velosity, U1, pressure, delta_h,solParam,U3);
        timeSolvind.push_back(dt);
        auto U_rho_Up = U_rho;
        for(size_t j = 0; j < U_rho.size();j++)
        {
            U_rho_Up[j].removeLast();
        }
        auto U_rho_velocity_Up = U_rho_velocity; U_rho_velocity_Up.removeLast();
        auto U_rho_energy_Up = U_rho_energy;     U_rho_energy_Up.removeLast();

        auto U_rho_Down = U_rho;
        for(size_t j = 0; j < U_rho.size();j++)
        {
            U_rho_Down[j].removeFirst();
        }
        auto U_rho_velocity_Down = U_rho_velocity; U_rho_velocity_Down.removeFirst();
        auto U_rho_energy_Down = U_rho_energy;     U_rho_energy_Down.removeFirst();

        solveFlux(U_rho_Up, U_rho_velocity_Up, U_rho_energy_Up, U_rho_Down, U_rho_velocity_Down, U_rho_energy_Down);

        auto res = additionalSolver.SolveEvolutionExplFirstOrder(F1, F2,F3,U_rho, U_rho_velocity,U_rho_energy,dt,delta_h);

        U_rho = res.first;
        U_rho_velocity = res.second[0];
        U_rho_energy = res.second[1];
        for(size_t j = 0; j < U_rho.size();j++)
        {
            U_rho[j][0] = U_rho[j][1];
            size_t sz = U_rho[j].size();
            U_rho[j][sz-1] = U_rho[j][sz-2];
        }
        U_rho_velocity[0]=U_rho_velocity[1];
        U_rho_velocity[U_rho_velocity.size()-1] = U_rho_velocity[U_rho_velocity.size()-2];
        U_rho_energy[0]=U_rho_energy[1];
        U_rho_energy[U_rho_energy.size() - 1] =U_rho_energy[U_rho_energy.size() - 2];

        //U3[U3.size() - 1] = rightParam.pressure/(solParam.Gamma-1)+0.5*pow(U2[U2.size() - 1]/U1[U1.size() - 1],2)*U1[U1.size() - 1];
    }
}
void AbstaractSolver::solveFlux(vector<Matrix> U_rho_Up, Matrix U2_Up, Matrix U3_Up, vector<Matrix> U_rho_Down, Matrix U2_Down, Matrix U3_Down)
{
    up_density   = U_rho_Up;
    down_density  = U_rho_Down;
    up_velocity  = U2_Up / up_density[0];
    down_velocity = U2_Down / down_density[0];
    up_pressure  = (U3_Up - Matrix::POW(up_velocity, 2)*0.5*up_density[0])*(solParam.Gamma-1);
    down_pressure = (U3_Down - Matrix::POW(down_velocity, 2)*0.5*down_density[0])*(solParam.Gamma-1);
    calcRiemanPStar();
    calcFliux();
}


void AbstaractSolver::calcRiemanPStar()
{
    const std::function<void(int&)> calcPStar = [this](int& i)
    {
        macroParam up(mixture);
        macroParam down(mixture);
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
            up.densityArray[j] = up_density[j][i];
        up.velocity = up_velocity[i];
        up.pressure = up_pressure[i];
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
            down.densityArray[j] = down_density[j][i];
        down.velocity = down_velocity[i];
        down.pressure = down_pressure[i];
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolver(up,down, solParam.Gamma);
        return;
    };
}
void AbstaractSolver::calcFliux()
{
    const std::function<void(int&)> calcFlux = [this](int& i)
    {
        auto point = rezultAfterPStart[i];
        double tempL = up_pressure[i]/(up_density[0][i]*UniversalGasConstant/molMass);
        double tempR = down_pressure[i]/(down_density[0][i]*UniversalGasConstant/molMass);
        auto du_dx = (down_velocity[i] - up_velocity[i])/delta_h;
        double Tx = point.pressure/(point.density*UniversalGasConstant/molMass);
        T[i] = Tx;
        double Pr = 2.0/3;
        double etta /*= additionalSolver.shareViscosity[0](upParam.temp, Tx,0,0);*/;
        double G =  (4.0/3*etta)*du_dx;
        double k = solParam.Gamma*UniversalGasConstant/molMass*etta/(solParam.Gamma-1)/Pr;

        double dt_dx = (tempR - tempL)/delta_h;
        double q =   -k*dt_dx;
        F1[0][i] = 0;
        for(size_t j = 0; j< F1.size(); j++)
            F1[j][i] = (point.density * point.velocity); // тут должна быть диффузия
        F2[i] = (/*F1.last()**/point.velocity + point.pressure -G); // тут явно какая-то ошибка
        F3[i] = ((point.pressure/(solParam.Gamma - 1) + point.density* pow(point.velocity,2)/2 + point.pressure) * point.velocity -G*point.velocity + q);
    };
}
void AbstaractSolver::  prepareVectors() // как будут задаваться условия в крайних ячейках
{
    for(size_t j = 0; j < mixture.NumberOfComponents; j++)
    {
        U_rho[j][0] = U_rho[j][1];
        U_rho[j][solParam.NumCell+1] = U_rho[j][solParam.NumCell];
    }
    U_rho_velocity[0]=U_rho_velocity[1];
    U_rho_energy[0]=U_rho_energy[1];

    U_rho_velocity[solParam.NumCell+1]=U_rho_velocity[solParam.NumCell]; // тут скорее всего будет иначе
    U_rho_energy[solParam.NumCell+1]=U_rho_energy[solParam.NumCell];

    T.resize(solParam.NumCell +2);

    timeSolvind.push_back(0);
}

