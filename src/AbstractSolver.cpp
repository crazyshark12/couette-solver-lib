#include <functional>

#include "abstractsolver.h"

AbstractSolver::AbstractSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_, SystemOfEquationType type, RiemannSolverType riemannType)
{
    mixture = mixture_;
    startParam=startParam_;
    solParam =solParam_;
    delta_h = 0;
    if(type == SystemOfEquationType::couette2)
    {
        auto *tmp = new Couette2();
        system = tmp;
    }
    if(type == SystemOfEquationType::soda)
    {
        auto *tmp = new Soda();
        system = tmp;
    }
    switch(riemannType)
    {
    case RiemannSolverType::HLLCSolver:
        {
                riemannSolver = new struct HLLCSolver();
                break;
        }
    case RiemannSolverType::HLLESolver:
        {
                riemannSolver = new struct HLLESolver();
                break;
        }
    case RiemannSolverType::HLLIsentropic:
        {
                riemannSolver = new struct HLLIsentropic();
                break;
        }
    case RiemannSolverType::HLLSimple:
        {
                riemannSolver = new struct HLLSimple();
                break;
        }
    case RiemannSolverType::ExacRiemanSolver:
        {
                riemannSolver = new struct ExacRiemanSolver();
                break;
        }
    }
    system->setBorderCondition(&border);
    system->setCoeffSolver(&coeffSolver);
    system->setMixture(mixture);
    system->setNumberOfCells(solParam.NumCell);
    system->setSolverParams(solParam);
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
    prepareVectorSizes();
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
    system->prepareSolving(points);
}

void AbstractSolver::prepareVectorSizes()
{
    points.resize(solParam.NumCell);
    for(size_t i = 0; i < points.size(); i++)
        points[i].densityArray.resize(mixture.NumberOfComponents);
    system->prepareVectorSizes();
}


void AbstractSolver::setDt()
{
    double max = system->getMaxVelocity();  // это нужно чтобы правильно подобрать временной шаг, чтобы соблюдался критерий КФЛ
    double dt;
    if(max!=0)
        dt = solParam.CFL*pow(delta_h,1)/max;
    else
        dt = 0.001;
    dt = 0.001; // тут фиксированный шаг
    timeSolvind.push_back(dt);
    //timeSolvind.push_back(0.00001);
    return;
}

void AbstractSolver::updatePoints()
{
    auto size = points.size()-1;
    #pragma omp parallel for schedule(static)
    for(size_t i = 1; i < size; i++)
    {
        if(i == 50)
            double x = 50;
        points[i].velocity_tau = system->getVelocityTau(i);
        points[i].velocity_normal = system->getVelocityNormal(i);
        points[i].velocity = system->getVelocity(i);
        points[i].density = system->getDensity(i);
        points[i].pressure = system->getPressure(i);

        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        {
            points[i].densityArray[j] =  system->U[j][i]; // тут косяк TODO
            points[i].fractionArray[j] = points[i].densityArray[j] / points[i].density;
        }
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
        points[i].temp = system->getTemp(i);
    }
    useBorder();
    system->updateBorderU(points);
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

