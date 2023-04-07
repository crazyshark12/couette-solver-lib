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

void HLLCSolver::computeHllcF()
{

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
