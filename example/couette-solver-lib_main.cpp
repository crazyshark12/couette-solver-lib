//#include "hllcsolver.h"
#include "godunovsolver.h"
#include "DataWriter.h"
#include "observer.h"
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
    //startParam.pressure = 218563.81; //218563.81
    startParam.densityArray[0] = argon.density;
    startParam.density = argon.density;
    startParam.temp = 140; //140
    startParam.velocity_tau = 100;

    solverParams solParam;
    solParam.NumCell     = 1000;    // Число расчтеных ячеек с учетом двух фиктивных ячеек
    solParam.Gamma    = 1.67;    // Показатель адиабаты
    solParam.CFL      = 1;    // Число Куранта
    solParam.MaxIter     = 100000000; // максимальное кол-во итареций
    solParam.Ma       = 0.1;    // Число маха

    double precision = 0.001; // точность
    Observer watcher(precision);
    watcher.setPeriodicity(1000);

    // это меняешь под себя. Он так создаст папку data
    // если не использовать setWriter, то записи не будет, но папка создастся, ибо она в конструкторе зашита
    // он автоматически очищает папку перед новым рассчётом
    DataWriter writer("D:/couette/couette");
    DataReader reader("D:/couette/couette/prev_data");
    reader.read();
    vector<macroParam> startParameters;
    reader.getPoints(startParameters);
    double T1wall = 1000;
    double T2wall = 1000;
    double velocity = 300;
    double h = 1;
    GodunovSolver solver(mixture,startParam,solParam, SystemOfEquationType::couette2, RiemannSolverType::HLLESolver);
    writer.setDelta_h(h / (solParam.NumCell));
    solver.setWriter(&writer);
    solver.setObserver(&watcher);
    solver.setDelta_h(h / (solParam.NumCell));
    solver.setBorderConditions(velocity,T2wall,T1wall);
    solver.solve();
}
