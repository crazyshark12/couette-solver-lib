#include "hllcsolver.h"


void HLLCSolver::solve()
{
    prepareSolving();
    for(auto i  = 0; i < solParam.MaxIter; i++)
    {
        // устанавливаем текущий временной шаг
        setDt();
        // Вычисляем вектор поточных членов и релаксационных членов
        computeF();
        // HLLС
        computeHllcF();
        // Вычисляем вектор релаксационных членов
        computeR();

        // Обновляем вектор U
        updateU();
    }
}

void HLLCSolver::computeF()
{
    // тут нужно будет знать вязкость и энтальпию их надо где-то рассчитывать в зависимости от смеси
}

void HLLCSolver::computeHllcF()
{
    for(size_t i = 0 ; i < solParam.NumCell; i++)
    {
        // Временные переменные
        double v0, v1, rho0, rho1, H0 , H1, avg_H, S0 , S1, S_star;
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
        //тут сделаю рассчет hllcF1,2,3
    }
}

void HLLCSolver::updateU()
{
    // Используем вектор потоков и релаксационные члены чтобы обновить U

    for(auto i  = 1; i < solParam.NumCell+1; i++)
    {
        for (int j = 1; j < U1.size(); j++)
        {
            U1[j][i] += (R[j][i] - (hllcF1[j][i] - hllcF1[j][i - 1]) / delta_h) * timeSolvind.last();
        }
        U2[i] += (0 - (hllcF2[i] - hllcF2[i - 1]) / delta_h) * timeSolvind.last();
        U3[i] += (0 - (hllcF3[i] - hllcF3[i - 1]) / delta_h) * timeSolvind.last();
    }
}
