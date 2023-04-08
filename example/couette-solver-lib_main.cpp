#include "hllcsolver.h"

int main()
{
    MixtureComponent argon;
    argon.name = "Ar";
    argon.density = 1.7839;
    argon.molarMass = 0.04;
    Mixture mixture({argon});

    macroParam startParam(mixture);
    startParam.fractionArray[0] = 1;
    startParam.densityArray[0] = argon.density;
    startParam.pressure = 101325;
    startParam.temp = 273;

    solverParams solParam;
    solParam.NumCell     = 10;    // Число расчтеных ячеек
    solParam.Gamma    = 1.67;    // Показатель адиабаты
    solParam.CFL      = 1;    // Число Куранта
    solParam.MaxIter     = 10; // максимальное кол-во шагов по времени
    solParam.Ma       = 1;    // Число маха

    HLLCSolver solver(mixture,startParam,solParam);
    solver.solve();
}
