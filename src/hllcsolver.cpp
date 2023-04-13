#include "hllcsolver.h"
#include <algorithm>

HLLCSolver::HLLCSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_)
{
    mixture = mixture_;
    startParam=startParam_;
    solParam =solParam_;
    delta_h = 1;
}

void HLLCSolver::solve()
{
    prepareSolving();
    writePoints(-1);
    double T = 0;
    for(auto i  = 0; i < solParam.MaxIter; i++)
    {
        // Устанавливаем текущий временной шаг
        setDt();
        T += timeSolvind.last();
        // Вычисляем вектор поточных членов и релаксационных членов
        computeF();
        // HLLС
        computeHllcF();
        // Вычисляем вектор релаксационных членов
        computeR();

        // Обновляем вектор U
        updateU();

        // Обновляем вектор макропараметров
        updatePoints();

        //записать данные, если это требуется
        //writePoints(timeSolvind.last()*100000); // наносек
        writePoints(T*1000000); // микросек
    }
}

void HLLCSolver::setBorderConditions(double up_velocity_, double up_temp_, double down_temp_)
{
    border.up_velocity =  up_velocity_;
    border.up_temp =  up_temp_;
    border.down_temp =  down_temp_;
    return;
}


void HLLCSolver::setStartCondition(macroParam start)
{
    startParam = start;
    mixture = startParam.mixture;
}

void HLLCSolver::setWriter(DataWriter *writer_)
{
    writer = writer_;
    isWriteData = true;
}

void HLLCSolver::writePoints(double i)
{
    if(isWriteData)
        writer->writeData(points,i);
}

void HLLCSolver::useBorder()
{
    //0
    points[0].mixture = mixture;
    points[0].density =points[1].density;
    points[0].densityArray =points[1].densityArray;
    points[0].fractionArray =points[1].fractionArray;
    points[0].velocity = -points[1].velocity;
    points[0].temp = -points[1].temp +  2*border.down_temp;
    // дополнительные рассчитываемые величины
    points[0].pressure = points[0].density * (UniversalGasConstant/mixture.molarMass()) * points[0].temp;
    points[0].soundSpeed = sqrt(solParam.Gamma*points[0].pressure/points[0].density);


    //solParam.NumCell-1
    points[solParam.NumCell-1].mixture = mixture;
    points[solParam.NumCell-1].density =points[solParam.NumCell-2].density;
    points[solParam.NumCell-1].densityArray =points[solParam.NumCell-2].densityArray;
    points[solParam.NumCell-1].fractionArray =points[solParam.NumCell-2].fractionArray;
    points[solParam.NumCell-1].velocity = -points[solParam.NumCell-2].velocity + 2*border.up_velocity;
    points[solParam.NumCell-1].temp = -points[solParam.NumCell-2].temp +  2*border.up_temp;
    // дополнительные рассчитываемые величины
    points[solParam.NumCell-1].pressure = points[solParam.NumCell-1].density * (UniversalGasConstant/mixture.molarMass()) * points[solParam.NumCell-1].temp;
    points[solParam.NumCell-1].soundSpeed = sqrt(solParam.Gamma*points[solParam.NumCell-1].pressure/points[solParam.NumCell-1].density);


}

void HLLCSolver::prepareSolving()
{
    U1.resize(mixture.NumberOfComponents);
    U2.resize(solParam.NumCell);
    U3.resize(solParam.NumCell);
    hllcF1.resize(mixture.NumberOfComponents);
    hllcF2.resize(solParam.NumCell-1);
    hllcF3.resize(solParam.NumCell-1);

    for(size_t i = 0 ; i <  U1.size(); i++)
    {
        U1[i].resize(solParam.NumCell);
        hllcF1[i].resize(solParam.NumCell-1);
    }
    points.resize(solParam.NumCell);
    for(size_t i = 1; i < points.size(); i++)
        points[i].densityArray.resize(mixture.NumberOfComponents);


    for(size_t i = 1; i < points.size()-1; i++)
    {
        points[i].mixture = mixture;
        points[i].temp = startParam.temp;
        points[i].fractionArray =  startParam.fractionArray;
        points[i].densityArray =  startParam.densityArray;
        double density = 0;
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
            density += points[i].densityArray[j]*points[i].fractionArray[j];
        points[i].density = density;
        points[i].pressure = points[i].density * (UniversalGasConstant/mixture.molarMass()) * points[i].temp;
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
        points[i].velocity = solParam.Ma*points[i].soundSpeed;
    }
    // для points[0] и points[solParam.NumCell-1] (!важно что идёт после цикла!)
    useBorder();

    for(auto i  = 0; i < solParam.NumCell; i++)
    {
        U1[0][i] = points[i].density;
        for(size_t j = j; j < mixture.NumberOfComponents; j++)
            U1[j][i] = points[i].densityArray[j] ;
        U2[i] = points[i].density*points[i].velocity;


        U3[i] = points[i].pressure/(solParam.Gamma-1)+0.5*pow(points[i].velocity,2)*points[i].density;


//        double n = Nav / points[i].mixture.molarMass() * points[i].density;
//        U3[i] = 3/2 * n * kB * points[i].temp + points[i].density * pow(points[i].velocity,2) / 2;
    }
    prepareVectors();
}

void HLLCSolver::computeF()
{
    for(size_t i = 0 ; i < solParam.NumCell; i++)
    {
        // Временные переменные
        macroParam p0(mixture), p1(mixture), p2(mixture);

        // Забираем известные макропараметры в (.)-ах [i-1], [i], [i+1]
        p1 = points[i];
        if(i==0)
        {
            p0 = p1;
            p2 = points[i + 1];
        }
        else if(i == solParam.NumCell-1)
        {
            p0 = points[i - 1];
            p2 = p1;
        }
        else
        {
            p0 = points[i - 1];
            p2 = points[i + 1];
        }

        // Рассчитываем производные в точке i
        double dv_dy = (p2.velocity - p0.velocity) / (2.0 * delta_h);
        double dT_dy = (p2.temp - p0.temp) / (2.0 * delta_h);
        vector<double> dy_dy(mixture.NumberOfComponents);

        //учёт граничных условий
        if(i == 0 || i == solParam.NumCell-1)
            fill(dy_dy.begin(), dy_dy.end(),border.get_dyc_dy());
        else
        {
            for(size_t j = 0 ; j <mixture.NumberOfComponents; j++)
            {
                dy_dy[j] = (p2.fractionArray[j] - p0.fractionArray[j])/ (2.0 * delta_h);
            }
        }
        // Расчет поточных членов
        // .....
        // сейчас так:
        double etta = coeffSolver.shareViscositySimple(points[i]);
        double lambda = coeffSolver.lambda(points[i]);
        for(size_t j = 0 ; j <mixture.NumberOfComponents; j++)
        {
            if(j!=0)
                F1[j][i] = -p1.density * mixture.getEffDiff(j) * dy_dy[j];
            else
                F1[j][i] = 0;
        }
        F2[i] = -etta * dv_dy;
        for(size_t j = 0 ; j <mixture.NumberOfComponents; j++)
        {
            F3[i]+= - p1.density * mixture.getEffDiff(j)*dy_dy[j] * mixture.getEntalp(i);
        }
        F3[i] += -lambda*dT_dy - etta*p1.velocity*dv_dy;
    }
}

void HLLCSolver::computeR()
{
    return;
}

void HLLCSolver::computeHllcF()
{
    for(size_t i = 0 ; i < solParam.NumCell-1; i++)
    {
        // Временные переменные
        double v0, v1, rho0, rho1, H0 , H1, avg_H, S0 , S1, S_star;
        vector<double> U1_star_0(U1.size()), U1_star_1(U1.size());
        double U2_star_0,U3_star_0, U2_star_1,U3_star_1;
        double avg_a, avg_v;

        v0 = U2[i]/U1[0][i];
        v1 = U2[i+1]/U1[0][i+1];
        rho0 = sqrt(U1[0][i]);
        rho1 = sqrt(U1[0][i+1]);
        avg_v = (rho0 * v0 + rho1 * v1) / (rho0 + rho1);
        H0 = (U3[i]/U1[0][i] + points[i].pressure)/U1[0][i];
        H1 = (U3[i+1]/U1[0][i+1] + points[i+1].pressure)/U1[0][i+1];
        avg_H = (rho0 * H0 + rho1 * H1) / (rho0 + rho1);
        avg_a = sqrt((solParam.Gamma - 1)*(avg_H - 0.5 * pow(avg_v,2)));
        S0 = avg_v - avg_a;
        S1 = avg_v + avg_a;
        S_star = (points[i+1].pressure - points[i].pressure + rho0*v0*(S0 - v0) - rho1*v1*(S1 - v1)) / (rho0*(S0 - v0) - rho1*(S1 - v1));

        double coeff_0 =U1[0][i] * ((S0 - v0)/ (S0 - S_star));
        double coeff_1 =U1[0][i+1] * ((S1 - v1)/ (S1 - S_star));
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        {
            U1_star_0[j] = coeff_0;
            U1_star_1[j] = coeff_1;
        }
        U2_star_0 = coeff_0 * S_star;
        U2_star_1 = coeff_1 * S_star;
        U3_star_0 = coeff_0 * (U3[i]/(pow(U1[0][i],2)) + (S_star - v0)*(S_star + (points[i].pressure)/(U1[0][i] * (S0 - v0))));
        U3_star_1 = coeff_1 * (U3[i+1]/(pow(U1[0][i+1],2)) + (S_star - v1)*(S_star + (points[i+1].pressure)/(U1[0][i+1] * (S1 - v1))));

        vector<double> hllc_f1 (U1.size());
        double hllc_f2, hllc_f3;
        if(S0 < 0)
        {
            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
            {
                hllc_f1[j] = F1[j][i];
            }
            hllc_f2 = F2[i];
            hllc_f3 = F3[i];
        }
        else if(S1 >= 0)
        {
            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
            {
                hllc_f1[j] = F1[j][i+1];
            }
            hllc_f2 = F2[i+1];
            hllc_f3 = F3[i+1];
        }
        else if( S0 <= 0 && S_star >= 0)
        {
            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
            {
                hllc_f1[j] = F1[j][i] + S0*(U1_star_0[j] - U1[j][i]);
            }
            hllc_f2 = F2[i] + S0*(U2_star_0 - U2[i]);
            hllc_f3 = F3[i] + S0*(U3_star_0 - U3[i]);
        }
        else if( S_star <= 0 && S1 >= 0)
        {
            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
            {
                hllc_f1[j] = F1[j][i+1] + S1*(U1_star_1[j] - U1[j][i+1]);
            }
            hllc_f2 = F2[i+1] + S0*(U2_star_1 - U2[i+1]);
            hllc_f3 = F3[i+1] + S0*(U3_star_1 - U3[i+1]);
        }
        // заполнение итоговых векторов
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        {
            hllcF1[j][i] = hllc_f1[j];
        }
        hllcF2[i] = hllc_f2;
        hllcF3[i] = hllc_f3;
    }
    return;
}

void HLLCSolver::updateU()
{
    // Используем вектор потоков и релаксационные члены чтобы обновить U
    for(auto i  = 1; i < solParam.NumCell-1; i++)
    {
        for (int j = 0; j < mixture.NumberOfComponents; j++)
        {
            U1[j][i] += (R[j][i] - (hllcF1[j][i] - hllcF1[j][i - 1]) / delta_h) * timeSolvind.last();
        }
        U2[i] += (0 - (hllcF2[i] - hllcF2[i - 1]) / delta_h) * timeSolvind.last();
        U3[i] += (0 - (hllcF3[i] - hllcF3[i - 1]) / delta_h) * timeSolvind.last();
    }
}

void HLLCSolver::updatePoints()
{
    for(size_t i = 1; i < points.size()-1; i++)
    {
        points[i].velocity = U2[i]/U1[0][i];
        points[i].pressure = (U3[i] - pow(points[i].velocity,2)*0.5*U1[0][i])*(solParam.Gamma - 1);
        points[i].density = U1[0][i];
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        {
            points[i].densityArray[j] =  U1[j][i];
            points[i].fractionArray[j] = points[i].densityArray[j] / points[i].density;
        }
        points[i].soundSpeed = sqrt(solParam.Gamma*points[i].pressure/points[i].density);
        computeT(points[i],i);
    }
    useBorder();
    return;
}

void HLLCSolver::computeT(macroParam &p, size_t i) // i - номер ячейки
{
    double U = U3[i] / p.density - pow(p.velocity,2)/2;
    double n = Nav / p.mixture.molarMass() * p.density;
    p.temp = U * 2/3 * p.density / (n * kB);



//    double n = Nav / p.mixture.molarMass() * p.density;
//    p.temp = (U3[i] - points[i].density * pow(points[i].velocity,2) / 2 )*(3/2 * n * kB ) ;
}
