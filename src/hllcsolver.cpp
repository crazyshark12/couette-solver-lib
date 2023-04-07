#include "hllcsolver.h"


void HLLCSolver::solve()
{
    prepareSolving();
    for(auto i  = 0; i < solParam.MaxIter; i++)
    {
        // Устанавливаем текущий временной шаг
        setDt();
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
    }
}

void HLLCSolver::computeF()
{
//    for(size_t i = 0 ; i < solParam.NumCell+1; i++)
//    {
//        // Временные переменные
//        macroParam p0(mixture), p1(mixture), p2(mixture);

//        // Забираем известные макропараметры в (.)-ах [i-1], [i], [i+1]
//        p1 = points[i];
//        switch (i)
//        {
//            case 0:
//                p0 = p1;
//                p2 = points[i + 1];
//                break;
//            case N_CELL - 1:
//                p0 = points[i - 1];
//                p2 = p1;
//                break;
//            default:
//                p0 = points[i - 1];
//                p2 = points[i + 1];
//                break;
//        }

//        // Вспомогательные величины (кончентрации и молярные доли)
//        QVector<double> n0 = {p0.rho[0] / Mixture::mass(0),
//                              p0.rho[1] / Mixture::mass(1)};
//        QVector<double> n2 = {p2.rho[0] / Mixture::mass(0),
//                              p2.rho[1] / Mixture::mass(1)};
//        QVector<double> x0 = {n0[0] / (n0[0] + n0[1]),
//                              n0[1] / (n0[0] + n0[1])};
//        QVector<double> x2 = {n2[0] / (n2[0] + n2[1]),
//                              n2[1] / (n2[0] + n2[1])};

//        // Рассчитываем производные в точке i
//        double dv_dx = (p2.v - p0.v) / (2.0 * DX);
//        double dT_dx = (p2.t - p0.t) / (2.0 * DX);
//        double dT12_dx = (p2.t12 - p0.t12) / (2.0 * DX);
//        double dT3_dx = (p2.t3 - p0.t3) / (2.0 * DX);
//        double dp_dx = (p2.p - p0.p) / (2.0 * DX);
//        QVector<double> dx_dx = {(x2[0] - x0[0]) / (2.0 * DX),
//                                 (x2[1] - x0[1]) / (2.0 * DX)};

//        // Расчет поточных членов
//        FlowMembersDc computer;
//        computer.initialize();
//        computer.compute(p1, dx_dx, dp_dx, dT_dx, dT12_dx, dT3_dx, dv_dx);

//        // Обновляем значения вектора поточных членов в (.) [i]
//        mutex.lock();
//        for (int j = 0; j < SYSTEM_ORDER; ++j)
//        {
//            F[j][i] = computer.flow()[j];
//        }
//        mutex.unlock();
//    };
}

void HLLCSolver::computeHllcF()
{
    for(size_t i = 0 ; i < solParam.NumCell+1; i++)
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
        H0 = (U3[i]/U1[0][i] + pres[i])/U1[0][i];
        H1 = (U3[i+1]/U1[0][i+1] + pres[i+1])/U1[0][i+1];
        avg_H = (rho0 * H0 + rho1 * H1) / (rho0 + rho1);
        avg_a = sqrt((solParam.Gamma - 1)*(avg_H - 0.5 * pow(avg_v,2)));
        S0 = avg_v - avg_a;
        S1 = avg_v + avg_a;
        S_star = (pres[i+1] - pres[i] + rho0*v0*(S0 - v0) - rho1*v1*(S1 - v1)) / (rho0*(S0 - v0) - rho1*(S1 - v1));

        double coeff_0 =U1[0][i] * ((S0 - v0)/ (S0 - S_star));
        double coeff_1 =U1[0][i+1] * ((S1 - v1)/ (S1 - S_star));
        for(size_t j = 0; j < mixture.NumberOfComponents; j++)
        {
            U1_star_0[j] = coeff_0;
            U1_star_1[j] = coeff_1;
        }
        U2_star_0 = coeff_0 * S_star;
        U2_star_1 = coeff_1 * S_star;
        U3_star_0 = coeff_0 * (U3[i]/(pow(U1[0][i],2)) + (S_star - v0)*(S_star + (pres[i])/(U1[0][i] * (S0 - v0))));
        U3_star_1 = coeff_1 * (U3[i+1]/(pow(U1[0][i+1],2)) + (S_star - v1)*(S_star + (pres[i+1])/(U1[0][i+1] * (S1 - v1))));

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
    for(auto i  = 1; i < solParam.NumCell+1; i++)
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
        computeT(points[i],i);
    }
    return;
}

void HLLCSolver::computeT(macroParam &p, size_t i) // i - номер ячейки
{
    // тут должно быть совсем иначе
    double U = U3[i] / p.density - pow(p.velocity,2)/2;
    double n = Nav / p.mixture.molarMass() * p.density;
    p.temp = U * 2/3 * p.density / n;
}
