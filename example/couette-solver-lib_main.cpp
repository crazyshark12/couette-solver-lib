#include "hllcsolver.h"
#include "DataWriter.h"
#include <filesystem>

namespace fs = std::filesystem;
int main()
{
    MixtureComponent argon;
    argon.name = "Ar";
    argon.density = 1.7839;
    argon.molarMass = 0.039948;
    argon.mass = 6.633521356992E-26;
    argon.epsilonDevK = 1.8845852298E-21;
    argon.sigma = 3.33E-10;
    Mixture mixture({argon});

    macroParam startParam(mixture);
    startParam.fractionArray[0] = 1;
    startParam.pressure = 218563.81; //218563.81
    //startParam.densityArray[0] = argon.density;
    startParam.temp = 140; //140

    solverParams solParam;
    solParam.NumCell     = 20;    // Число расчтеных ячеек
    solParam.Gamma    = 1.67;    // Показатель адиабаты
    solParam.CFL      = 1;    // Число Куранта
    solParam.MaxIter     = 100; // максимальное кол-во шагов по времени
    solParam.Ma       = 0;    // Число маха

    // это меняешь под себя. Он так создаст папку data
    // если не использовать setWriter, то записи не будет, но папка создастся, ибо она в конструкторе зашита
    // он автоматически очищает папку перед новым рассчётом
    DataWriter writer("D:/couette/couette-solver-lib");

    double T1wall = 500;
    double T2wall = 500;
    double velocity = 2;
    double h = 0.2;
    HLLCSolver solver(mixture,startParam,solParam);
    solver.setWriter(&writer);
    solver.setDelta_h(h / solParam.NumCell);
    solver.setBorderConditions(velocity,T2wall,T1wall);
    solver.solve();

}
