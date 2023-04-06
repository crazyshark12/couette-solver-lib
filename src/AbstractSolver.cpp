#include <functional>

#include "abstractsolver.h"

void AbstaractSolver::setMixture(Mixture mixture_)
{
    mixture =  mixture_;
}

// тут нужно задавать с помощью заранее заданных данных (давление и температуры по идее) задавать upParam и
// downParam и с помощью них заполнить две большие группы ячеек: которые примыкают к верхней пластине
// и которые примыкают к нижней ( разделение где-то посередине )
void AbstaractSolver::prepareSolving()
{

    U_rho_velocity.resize(solParam.NumCell+2);
    U_rho_energy.resize(solParam.NumCell+2);
    U_rho.resize(mixture.NumberOfComponents);
    for(size_t i = 0 ; i <  U_rho.size(); i++)
        U_rho[i].resize(solParam.NumCell+2);


//    ???  как мы устанавливаем эти параметры ??? (взял из argonSolver)
//    upParam.density = upParam.pressure /(UniversalGasConstant/molMass * upParam.temp);
//    upParam.soundSpeed = sqrt(solParam.Gamma*upParam.pressure/upParam.density);
//    upParam.velocity = solParam.Ma*upParam.soundSpeed;

//    downParam.density = ((solParam.Gamma + 1)* pow(solParam.Ma,2))/(2 + (solParam.Gamma -1)* pow(solParam.Ma,2))*upParam.density;
//    downParam.pressure = (pow(solParam.Ma,2) * 2* solParam.Gamma - (solParam.Gamma - 1))/((solParam.Gamma +1))*upParam.pressure;
//    downParam.temp = downParam.pressure/(downParam.density*UniversalGasConstant/molMass);
//    auto g = solParam.Gamma;
//    auto m = solParam.Ma;
//    auto temp2 = (2*g*m*m/(g+1) - (g-1)/(g+1))/((g+1)* m*m/(2 + (g-1)*m*m));

//    //downParam.velocity = downParam.density*downParam.velocity/upParam.density;

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
            for(size_t j = 1; j < mixture.NumberOfComponents; j++)
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
        // тут устанавливает временной шаг
        Matrix velocity = U_rho_velocity/U_rho[0];
        auto pressure = (U_rho_energy - Matrix::POW(velocity,2)*0.5*U_rho[0])*(solParam.Gamma - 1);
        Matrix sound_speed = Matrix::SQRT(pressure*solParam.Gamma/U_rho[0]);
        auto temp = velocity + sound_speed;                     // что это, типо абсолютная скорость?
        auto max =*std::max_element(temp.begin(), temp.end());  // видимо это нужно чтобы правильно подобрать временной шаг, чтобы соблюдался критерий КФЛ
        double dt = solParam.CFL*pow(delta_h,1)/max;
        timeSolvind.push_back(dt);

        //считает вектор потоков
        solveFlux();

        //где-то тут должно быть вычисление релаксационных членов


        //с помощью потоков устанавливает значение значение U_new, я эту функцию поменял, но она мне не нравится ибо возвращаемый тип неочевидный
        auto res = additionalSolver.SolveEvolutionExplFirstOrder(F1, F2,F3,U_rho, U_rho_velocity,U_rho_energy,dt,delta_h);

        U_rho = res.first;
        U_rho_velocity = res.second[0];
        U_rho_energy = res.second[1];

        //тут происходит что-то странное, видимо своего рода граничные условия
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

    }
}
void AbstaractSolver::solveFlux()
{
    auto density   = U_rho;
    auto velocity  = U_rho_velocity / U_rho[0];
    auto pressure  = (U_rho_energy - Matrix::POW(velocity, 2)*0.5*density[0])*(solParam.Gamma-1);
    calcRiemanPStar(density,velocity,pressure);
    calcFliux(density,velocity,pressure);
}


void AbstaractSolver::calcRiemanPStar(vector<Matrix> density_, Matrix velocity_ , Matrix pressure_)
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
void AbstaractSolver::calcFliux(vector<Matrix> density_, Matrix velocity_ , Matrix pressure_)
{
    for(size_t i = 0; i < velocity_.size()-1; i++)
    {
        auto point = rezultAfterPStart[i];
        double tempU = pressure_[i+1]/(density_[0][i+1]*UniversalGasConstant/molMass);
        double tempD = pressure_[i]/(density_[0][i]*UniversalGasConstant/molMass);
        auto du_dx = (velocity_[i] - velocity_[i+1])/delta_h;
        double Tx = point.pressure/(point.density*UniversalGasConstant/molMass);
        T[i] = Tx;
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
void AbstaractSolver::  prepareVectors() // подготовка размеров нужных векторов и то как будут задаваться условия в крайних ячейках
{
    for(size_t j = 0; j < mixture.NumberOfComponents; j++)
    {
        U_rho[j][0] = U_rho[j][1];
        U_rho[j][solParam.NumCell+1] = U_rho[j][solParam.NumCell];
    }
    U_rho_velocity[0]=U_rho_velocity[1];
    U_rho_energy[0]=U_rho_energy[1];

    U_rho_velocity[solParam.NumCell+1]=U_rho_velocity[solParam.NumCell];
    U_rho_energy[solParam.NumCell+1]=U_rho_energy[solParam.NumCell];

    T.resize(solParam.NumCell +2);

    timeSolvind.push_back(0);
}

