#include <functional>

#include "abstractsolver.h"

AbstractSolver::AbstractSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_)
{
    mixture = mixture_;
    startParam=startParam_;
    solParam =solParam_;
    delta_h = 0;
}

void AbstractSolver::setStartCondition(macroParam start)
{
    startParam = start;
    mixture = startParam.mixture;
}

void AbstractSolver::setBorderConditions(double up_velocity_, double up_temp_, double down_temp_)
{
    border.up_velocity =  up_velocity_;
    border.up_temp =  up_temp_;
    border.down_temp =  down_temp_;
    return;
}

void AbstractSolver::setWriter(DataWriter *writer_)
{
    writer = writer_;
    isWriteData = true;
}

void AbstractSolver::setObserver(Observer* obs)
{
    observer = obs;
    watcherIteration = observer->getPeriodicity();
    isObserverWatching = true;
}

void AbstractSolver::writePoints(double i)
{
    if(isWriteData)
        writer->writeData(points,i);
}

void AbstractSolver::setDelta_h(double dh)
{
    delta_h = dh;
}

void AbstractSolver::prepareSolving()
{
    prepareVectors();
    for(size_t i = 1; i < points.size()-1; i++)
    {
        points[i].mixture = mixture;
        points[i].temp = startParam.temp;
        points[i].fractionArray =  startParam.fractionArray;
        points[i].pressure = startParam.pressure;
        points[i].density = startParam.pressure * mixture.molarMass()/(UniversalGasConstant * startParam.temp);
        points[i].densityArray[0] =  points[i].density;
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
        points[i].velocity_tau = solParam.Ma*points[i].soundSpeed;
        points[i].velocity_normal = 0;
        points[i].velocity = abs(points[i].velocity_tau);
    }
    // для points[0] и points[solParam.NumCell-1] (!важно что идёт после цикла!)
    useBorder();

    for(auto i  = 0; i < solParam.NumCell; i++)
    {
        U1[0][i] = points[i].density;
        for(size_t j = j; j < mixture.NumberOfComponents; j++)
            U1[j][i] = points[i].densityArray[j] ;
        U2[i] = points[i].density*points[i].velocity_tau;
        U2_normal[i] = points[i].density*points[i].velocity_normal;


        //U3[i] = points[i].pressure/(solParam.Gamma-1)+0.5*pow(points[i].velocity,2)*points[i].density;
        U3[i] = (3*points[i].pressure)/2 + 0.5*pow(points[i].velocity,2)*points[i].density;
    }
}

void AbstractSolver::prepareVectors()
{

    U1.resize(mixture.NumberOfComponents);
    U2.resize(solParam.NumCell);
    U2_normal.resize(solParam.NumCell);
    U3.resize(solParam.NumCell);
    for(size_t i = 0 ; i <  U1.size(); i++)
        U1[i].resize(solParam.NumCell);

    points.resize(solParam.NumCell);
    for(size_t i = 0; i < points.size(); i++)
        points[i].densityArray.resize(mixture.NumberOfComponents);

    F1.resize(mixture.NumberOfComponents);
    for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        F1[j].resize(solParam.NumCell);
    F2.resize(solParam.NumCell);
    F2_normal.resize(solParam.NumCell);
    F3.resize(solParam.NumCell);

    R.resize(solParam.NumCell);
    timeSolvind.push_back(0);
}


void AbstractSolver::setDt()
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

void AbstractSolver::useBorder()
{
    //0
    points[0].mixture = mixture;
    points[0].density =points[1].density;
    points[0].densityArray =points[1].densityArray;
    points[0].fractionArray =points[1].fractionArray;
    points[0].velocity_tau = -points[1].velocity_tau + 2*border.down_velocity;
    points[0].velocity_normal = 0;
    points[0].velocity = abs(points[0].velocity_tau);
    points[0].temp = -points[1].temp +  2*border.down_temp;
    // дополнительные рассчитываемые величины
    points[0].pressure = points[0].density * (UniversalGasConstant/mixture.molarMass()) * points[0].temp;
    points[0].soundSpeed = sqrt(solParam.Gamma*points[0].pressure/points[0].density);


    //solParam.NumCell-1
    points[solParam.NumCell-1].mixture = mixture;
    points[solParam.NumCell-1].density =points[solParam.NumCell-2].density;
    points[solParam.NumCell-1].densityArray =points[solParam.NumCell-2].densityArray;
    points[solParam.NumCell-1].fractionArray =points[solParam.NumCell-2].fractionArray;
    points[solParam.NumCell-1].velocity_tau = -points[solParam.NumCell-2].velocity_tau + 2*border.up_velocity;
    points[solParam.NumCell-1].velocity_normal = 0;
    points[solParam.NumCell-1].velocity = abs(points[solParam.NumCell-1].velocity_tau);
    points[solParam.NumCell-1].temp = -points[solParam.NumCell-2].temp +  2*border.up_temp;
    // дополнительные рассчитываемые величины
    points[solParam.NumCell-1].pressure = points[solParam.NumCell-1].density * (UniversalGasConstant/mixture.molarMass()) * points[solParam.NumCell-1].temp;
    points[solParam.NumCell-1].soundSpeed = sqrt(solParam.Gamma*points[solParam.NumCell-1].pressure/points[solParam.NumCell-1].density);


}
void AbstractSolver::UpdateBorderU()
{
    for(int i : { 0, solParam.NumCell-1})
    {
        U1[0][i] = points[i].density;
        for(size_t j = j; j < mixture.NumberOfComponents; j++)
            U1[j][i] = points[i].densityArray[j];
        U2[i] = points[i].density*points[i].velocity_tau;
        U2_normal[i] = points[i].density*points[i].velocity_normal;
        //U3[i] = points[i].pressure/(solParam.Gamma-1)+0.5*pow(points[i].velocity,2)*points[i].density;
        U3[i] = (3*points[i].pressure)/2 + 0.5*pow(points[i].velocity,2)*points[i].density;
    }
}
double AbstractSolver::computeT(macroParam p, size_t i) // i - номер ячейки
{
    double U = U3[i] / p.density - pow(p.velocity,2)/2;
    double n = Nav / p.mixture.molarMass() * p.density;
    return U * 2/3 * p.density / (n * kB);
//    double n = Nav / p.mixture.molarMass() * p.density;
    //    p.temp = (U3[i] - points[i].density * pow(points[i].velocity,2) / 2 )*(3/2 * n * kB ) ;
}

bool AbstractSolver::observerCheck(size_t currentIteration)
{
    if(currentIteration%watcherIteration == 0)
    {
        observer->remember(points);
        return true;
    }
    if(currentIteration%watcherIteration == 1)
    {
        return observer->checkDifference(points);
    }
    return true;
}

