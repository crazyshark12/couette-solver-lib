#include "hllcsolver.h"
#include <algorithm>

//HLLCSolver::HLLCSolver(Mixture mixture_, macroParam startParam_, solverParams solParam_)
//{
//    mixture = mixture_;
//    startParam=startParam_;
//    solParam =solParam_;
//    delta_h = 1;
//}

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
        writePoints(T*1000000); // микросек
    }
}

void HLLCSolver::prepareVectors()
{
    AbstractSolver::prepareVectors();
    fluxF1.resize(mixture.NumberOfComponents);
    fluxF2.resize(solParam.NumCell-1);
    fluxF3.resize(solParam.NumCell-1);
    for(size_t i = 0 ; i <  U1.size(); i++)
    {
        fluxF1[i].resize(solParam.NumCell-1);
    }
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
//        double etta =  0;
//        double lambda = 0;
        for(size_t j = 0 ; j <mixture.NumberOfComponents; j++)
        {
            if(j!=0)
                F1[j][i] = -p1.density * mixture.getEffDiff(j) * dy_dy[j];
            else
                F1[j][i] = 0;
        }
        F2[i] = -etta * dv_dy;
        F3[i] = 0;
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
        if(i==solParam.NumCell-3)
            double x = 34;
        // Временные переменные
        double v0, v1,a0,a1, rho0, rho1, H0 , H1, avg_H, S0 , S1, S_star;
        vector<double> U1_star_0(U1.size()), U1_star_1(U1.size());
        double U2_star_0,U3_star_0, U2_star_1,U3_star_1;
        double avg_a, avg_v;

        v0 = U2[i]/U1[0][i];
        v1 = U2[i+1]/U1[0][i+1];
        rho0 = sqrt(U1[0][i]);
        rho1 = sqrt(U1[0][i+1]);

        avg_v = (rho0 * v0 + rho1 * v1) / (rho0 + rho1);
        H0 = (3*points[i].pressure)/(2*U1[0][i]);
        H1 = (3*points[i+1].pressure)/(2*U1[0][i+1]);
        avg_H = (rho0 * H0 + rho1 * H1) / (rho0 + rho1);
        avg_a = sqrt((solParam.Gamma - 1)*(avg_H - 0.5 * pow(avg_v,2)));
        S0 = (avg_v - avg_a);
        S1 = (avg_v + avg_a);
//        S0 = min(v0, v1);
//        S1 = max(v0, v1);
        S_star = (points[i+1].pressure - points[i].pressure +
                pow(rho0,2)*v0*(S0 - v0) - pow(rho1,2)*v1*(S1 - v1))
                / (pow(rho0,2)*(S0 - v0) - pow(rho1,2)*(S1 - v1));

//        if(!(S0 <= S_star && S_star <= S1))
//        {
//            std::cout <<"for i = "<< i<<std::endl;
//            std::cout <<"S0 = "<< S0<<std::endl;
//            std::cout <<"S* = "<< S_star<<std::endl;
//            std::cout <<"S1 = "<< S1<<std::endl;
//            std::cout <<"a = "<< avg_a<<std::endl;
//            std::cout <<"v1 = "<< v1<<std::endl;
//            std::cout <<"v0 = "<< v0<<std::endl;
//            std::cout <<"avg_v = "<< avg_v<<std::endl;
//            std::cout <<"H1 = "<< H1<<std::endl;
//            std::cout <<"H0 = "<< H0<<std::endl;
//            std::cout <<"rho_0 = "<< rho0<<std::endl;
//            std::cout <<"rho_1 = "<< rho0<<std::endl;
//            std::cout <<"avg_h = "<< avg_H<<std::endl;
//            std::cout <<"dp = "<< points[i+1].pressure - points[i].pressure<<std::endl;

//        }
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
        if(S0 >= 0)
        {
            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
            {
                hllc_f1[j] = F1[j][i];
            }
            hllc_f2 = F2[i];
            hllc_f3 = F3[i];
        }
        else if(S1 <= 0)
        {
            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
            {
                hllc_f1[j] = F1[j][i+1];
            }
            hllc_f2 = F2[i+1];
            hllc_f3 = F3[i+1];
        }
        else if(S0<=0 && S1>=0)
        {
            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
            {
                hllc_f1[j] = (S1 * F1[j][i] - S0 * F1[j][i+1] + S0*S1*(U1[j][i+1] - U1[j][i]))/(S1 - S0);
            }
            hllc_f2 = (S1 * F2[i] - S0 * F2[i+1] + S0*S1*(U2[i+1] - U2[i]))/(S1 - S0);
            hllc_f3 = (S1 * F3[i] - S0 * F3[i+1] + S0*S1*(U3[i+1] - U3[i]))/(S1 - S0) ;
        }
//        else if( S0 <= 0 && S_star >= 0)
//        {
//            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
//            {
//                    hllc_f1[j] = F1[j][i] + S0*(U1_star_0[j] - U1[j][i]);
//            }
//            hllc_f2 = F2[i] + S0*(U2_star_0 - U2[i]);
//            hllc_f3 = F3[i] + S0*(U3_star_0 - U3[i]);
//        }
//        else if( S_star <= 0 && S1 >= 0)
//        {
//            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
//            {
//                hllc_f1[j] = F1[j][i+1] + S1*(U1_star_1[j] - U1[j][i+1]);
//            }
//            hllc_f2 = F2[i+1] + S1*(U2_star_1 - U2[i+1]);
//            hllc_f3 = F3[i+1] + S1*(U3_star_1 - U3[i+1]);
//        }
        // заполнение итоговых векторов
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        {
            fluxF1[j][i] = hllc_f1[j];
        }
        fluxF2[i] = hllc_f2;
        fluxF3[i] = hllc_f3;
    }
//    std::cout<<"--------------"<<std::endl;
    return;
}

void HLLCSolver::updateU()
{
    // Используем вектор потоков и релаксационные члены чтобы обновить U
    for(auto i  = 1; i < solParam.NumCell-1; i++)
    {
        for (int j = 0; j < mixture.NumberOfComponents; j++)
        {
            U1[j][i] += (/*R[j][i]*/0 - (fluxF1[j][i] - fluxF1[j][i - 1]) / delta_h) * timeSolvind.last();
        }
        U2[i] += (0 - (fluxF2[i] - fluxF2[i - 1]) / delta_h) * timeSolvind.last();
        U3[i] += (0 - (fluxF3[i] - fluxF3[i - 1]) / delta_h) * timeSolvind.last();
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
        points[i].temp = computeT(points[i],i);
    }
    useBorder();
    UpdateBorderU();
    return;
}

