#include <functional>

#include "abstractsolver.h"

//void AbstaractSolver::prepareSolving()
//{
//    U1.resize(solParam.NumCell+2);
//    U2.resize(solParam.NumCell+2);
//    U3.resize(solParam.NumCell+2);
//    U4.resize(solParam.NumCell+2);
//    R.resize(solParam.NumCell+2);

//    leftParam.density = leftParam.pressure /(UniversalGasConstant/molMass * leftParam.temp);
//    double Z = additionalSolver.ZCO2Vibr(leftParam.temp);
//    double Cv = 5.0/2 * kB/mass+ additionalSolver.CVibr(leftParam.temp,Z);
//    solParam.Gamma = (UniversalGasConstant/molMass + Cv)/Cv;
//    leftParam.soundSpeed = sqrt(solParam.Gamma *UniversalGasConstant/molMass * leftParam.temp);
//    leftParam.velocity = solParam.Ma*leftParam.soundSpeed;
//    leftParam.tempIntr = leftParam.temp;
//    QFile pythonFile(QDir::currentPath() + "\\Fun.py");
//    if( !pythonFile.open(QFile::ReadOnly) )
//    {
//        QMessageBox msgBox;
//        msgBox.setWindowTitle("Ошибка");
//        msgBox.setText("Нет файла с питоновским расчетом граничных значений");
//        msgBox.exec();
//        breaksolve = true;
//        return;
//    }
//    rightParam = additionalSolver.bondaryConditionPython(leftParam, solParam);
//    double leftEvibr = additionalSolver.vibrEnergy(0,leftParam.temp);
//    double leftFullEnergy = 5.0/2*kB*leftParam.temp/mass + leftEvibr;
//    double rightEVibr = additionalSolver.vibrEnergy(0,rightParam.temp);
//    double rightFullEnergy = 5.0/2*kB*rightParam.temp/mass + rightEVibr;

//    for(auto i  = 1; i < solParam.NumCell+1; i++)
//    {
//        if(i < solParam.NumCell/3 +1)
//        {
//            U1[i] = leftParam.density;
//            U2[i] = leftParam.density*leftParam.velocity;
//            U3[i] = leftParam.density*(leftFullEnergy + pow(leftParam.velocity,2)/2);
//            U4[i] = leftParam.density*leftEvibr;
//        }

//        else
//        {
//            U1[i] = rightParam.density;
//            U2[i] = rightParam.density*rightParam.velocity;
//            U3[i] = rightParam.density*(rightFullEnergy + pow(rightParam.velocity,2)/2);
//            U4[i] = rightParam.density*rightEVibr;
//        }
//    }
//    prepareVectors();
//}

void AbstaractSolver::calcRiemanPStar()
{
    const std::function<void(int&)> calcPStar = [this](int& i)
    {
        macroParam left;
        macroParam right;
        mutex.lock();
        left.density = left_density[i];
        left.velocity = left_velocity[i];
        left.pressure = left_pressure[i];
        right.density = right_density[i];
        right.velocity = right_velocity[i];
        right.pressure = right_pressure[i];
        mutex.unlock();
        rezultAfterPStart[i] = additionalSolver.ExacRiemanSolverCorrect(left,right, solParam.Gamma, solParam.Gamma);
        return;
    };
}

void AbstaractSolver::prepareVectors()
{
    double x_right =solParam.lambda*solParam.lambdaSol; //% правая граница
    delta_h = (x_right) / solParam.NumCell;
    x.clear();
    x.push_back(0+0.5*delta_h);
    for(auto i = 1; i < solParam.NumCell; i++)
        x.push_back(x[i-1] + delta_h);
    x.push_back(x_right);
    x.insert(x.begin(), 0);
    U1[0]=U1[1];
    U2[0]=solParam.typeLeftBorder*U2[1];
    U3[0]=U3[1];
    U4[0]=U4[1];
    U5[0]=U5[1];
    U1[solParam.NumCell+1]=U1[solParam.NumCell];
    U2[solParam.NumCell+1]=solParam.typeRightBorder*U2[solParam.NumCell];
    U3[solParam.NumCell+1]=U3[solParam.NumCell];
    U4[solParam.NumCell+1]=U4[solParam.NumCell];
    U5[solParam.NumCell+1]=U5[solParam.NumCell];

    timeSolvind.push_back(0);
    for(int i = 0 ; i<  solParam.NumCell+1; i++)
        vectorForParallelSolving.push_back(i);
    F1.resize(solParam.NumCell+1);
    F2.resize(solParam.NumCell+1);
    F3.resize(solParam.NumCell+1);
    F4.resize(solParam.NumCell+1);
    F5.resize(solParam.NumCell+1);
    P.resize(solParam.NumCell+2);
    Q_v.resize(solParam.NumCell+2);
    Q_v3.resize(solParam.NumCell+2);
    Q_t.resize(solParam.NumCell+2);
    R.resize(solParam.NumCell+1);
    R_1.resize(solParam.NumCell+1);
    R_2.resize(solParam.NumCell+1);
    T.resize(solParam.NumCell +2);
    Tv.resize(solParam.NumCell +2);
    T12.resize(solParam.NumCell +2);
    T3.resize(solParam.NumCell +2);
    Ent.resize(solParam.NumCell +2);
    Ent2.resize(solParam.NumCell +2);
    B_v.resize(solParam.NumCell +2);
    E_Z.resize(solParam.NumCell +2);
    PR.resize(solParam.NumCell +2);
    rezultAfterPStart.resize(solParam.NumCell+1);
}
