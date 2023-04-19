#include <functional>

#include "abstractsolver.h"

AbstaractSolver::AbstaractSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_)
{
    mixture = mixture_;
    startParam=startParam_;
    solParam =solParam_;
    delta_h = 0;
}
void AbstaractSolver::setBorderConditions(double up_velocity_, double up_temp_, double down_temp_)
{
    border.up_velocity =  up_velocity_;
    border.up_temp =  up_temp_;
    border.down_temp =  down_temp_;
    return;
}


void AbstaractSolver::setStartCondition(macroParam start)
{
    startParam = start;
    mixture = startParam.mixture;
}

void AbstaractSolver::setWriter(DataWriter *writer_)
{
    writer = writer_;
    isWriteData = true;
}

void AbstaractSolver::writePoints(double i)
{
    if(isWriteData)
        writer->writeData(points,i);
}

void AbstaractSolver::setDelta_h(double dh)
{
    delta_h = dh;
}

void AbstaractSolver::prepareSolving()
{
    prepareVectors();
    for(size_t i = 0; i < points.size(); i++)
    {
        points[i].pressure = startParam.pressure;
        points[i].temp = startParam.temp;
        points[i].density = startParam.pressure /(UniversalGasConstant/mixture.molarMass() * startParam.temp);
        points[i].soundSpeed = sqrt(solParam.Gamma*startParam.pressure/points[i].density);
        points[i].velocity = solParam.Ma*startParam.soundSpeed;
        points[i].fractionArray =  startParam.fractionArray;
        points[i].densityArray =  startParam.densityArray;
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        {
            points[i].densityArray[j] = startParam.fractionArray[j] *  points[i].density;
        }
    }

    //downParam.velocity = downParam.density*downParam.velocity/upParam.density;

    for(auto i  = 0; i < solParam.NumCell; i++) // тут возможно от 1 до solParam.NumCell-1
    {
        U1[0][i] = startParam.density;
        for(size_t j = j; j < mixture.NumberOfComponents; j++)
            U1[j][i] = startParam.densityArray[j] ;
        U2[i] = startParam.density*startParam.velocity;
        U3[i] = startParam.pressure/(solParam.Gamma-1)+0.5*pow(startParam.velocity,2)*startParam.density; // скорее всего иначе
    }

}

void AbstaractSolver::prepareVectors()
{
    U2.resize(solParam.NumCell);
    U3.resize(solParam.NumCell);
    U1.resize(mixture.NumberOfComponents);
    for(size_t i = 0 ; i <  U1.size(); i++)
        U1[i].resize(solParam.NumCell);
    points.resize(solParam.NumCell);

    F1.resize(mixture.NumberOfComponents);
    for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        F1[j].resize(solParam.NumCell);
    F2.resize(solParam.NumCell);
    F3.resize(solParam.NumCell);
    R.resize(solParam.NumCell);
    timeSolvind.push_back(0);
}




void AbstaractSolver::setDt()
{
    Matrix velocity = U2/U1[0];
    auto pressure = (U3 - Matrix::POW(velocity,2)*0.5*U1[0])*(solParam.Gamma - 1);
    auto temp = velocity;
    auto max =*std::max_element(temp.begin(), temp.end());  // это нужно чтобы правильно подобрать временной шаг, чтобы соблюдался критерий КФЛ
    double dt = solParam.CFL*pow(delta_h,1)/max;
    timeSolvind.push_back(dt);
    //timeSolvind.push_back(0.00001);
    return;
}

void AbstaractSolver::updatePoints()
{
    for(size_t i = 0; i < points.size(); i++)
    {
        points[i].velocity = U2[i]/U1[0][i];
        points[i].pressure = (U3[i] - pow(points[i].velocity,2)*0.5*U1[i][0])*(solParam.Gamma - 1);
        points[i].density = U1[0][i];
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        {
            points[i].densityArray[j] =  U1[j][i];
            points[i].fractionArray[j] = points[i].densityArray[j] / points[i].density;
        }
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
        // тут ещё должна находиться температура
    }
    return;
}


