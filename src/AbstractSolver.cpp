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

    U2.resize(solParam.NumCell+2);
    U3.resize(solParam.NumCell+2);
    U1.resize(mixture.NumberOfComponents);
    for(size_t i = 0 ; i <  U1.size(); i++)
        U1[i].resize(solParam.NumCell+2);


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
            U1[0][i] = upParam.density;
            for(size_t j = 1; j < mixture.NumberOfComponents; j++)
                U1[j][i] = upParam.density * upParam.massFraction[j];
            U2[i] = upParam.density*upParam.velocity;
            U3[i] = upParam.pressure/(solParam.Gamma-1)+0.5*pow(upParam.velocity,2)*upParam.density; // скорее всего иначе
        }

        else
        {
            U1[0][i] = downParam.density;
            for(size_t j = 1; j < mixture.NumberOfComponents; j++)
                U1[j][i] = downParam.density * downParam.massFraction[j];
            U2[i] = downParam.density*downParam.velocity;
            U3[i] = downParam.pressure/(solParam.Gamma-1)+0.5*pow(downParam.velocity,2)*downParam.density; // скорее всего иначе
        }
    }
    prepareVectors();
}




void AbstaractSolver::setDt()
{
    Matrix velocity = U2/U1[0];
    auto pressure = (U3 - Matrix::POW(velocity,2)*0.5*U1[0])*(solParam.Gamma - 1);
    Matrix sound_speed = Matrix::SQRT(pressure*solParam.Gamma/U1[0]);
    auto temp = velocity + sound_speed;                     // что это, типо абсолютная скорость?
    auto max =*std::max_element(temp.begin(), temp.end());  // это нужно чтобы правильно подобрать временной шаг, чтобы соблюдался критерий КФЛ
    double dt = solParam.CFL*pow(delta_h,1)/max;
    timeSolvind.push_back(dt);
    return;
}
void AbstaractSolver::  prepareVectors() // подготовка размеров нужных векторов и то как будут задаваться условия в крайних ячейках
{
    for(size_t j = 0; j < mixture.NumberOfComponents; j++)
    {
        U1[j][0] = U1[j][1];
        U1[j][solParam.NumCell+1] = U1[j][solParam.NumCell];
    }
    U2[0]=U2[1];
    U3[0]=U3[1];

    U2[solParam.NumCell+1]=U2[solParam.NumCell];
    U3[solParam.NumCell+1]=U3[solParam.NumCell];

    T.resize(solParam.NumCell +2);

    timeSolvind.push_back(0);
}

