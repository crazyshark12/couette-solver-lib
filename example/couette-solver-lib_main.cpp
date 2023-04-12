#include "hllcsolver.h"

int main()
{
    MixtureComponent argon;
    argon.name = "Ar";
    argon.density = 1.7839;
    argon.molarMass = 0.04;
    argon.mass = 6.6335209E-26;
    argon.epsilonDevK = 1.884585E-21;
    argon.sigma = 3.33E-10;
    Mixture mixture({argon});

    macroParam startParam(mixture);
    startParam.fractionArray[0] = 1;
    startParam.densityArray[0] = argon.density;
    startParam.temp = 100;

    solverParams solParam;
    solParam.NumCell     = 10;    // Число расчтеных ячеек
    solParam.Gamma    = 1.67;    // Показатель адиабаты
    solParam.CFL      = 1;    // Число Куранта
    solParam.MaxIter     = 10; // максимальное кол-во шагов по времени
    solParam.Ma       = 0;    // Число маха

    double T1wall = 300;
    double T2wall = 700;
    double velocity = 2;
    double h = 2;
    HLLCSolver solver(mixture,startParam,solParam);
    solver.setDelta_h(h / solParam.NumCell);
    solver.setBorderConditions(velocity,T2wall,T1wall);
    solver.solve();
}
