#include "hllcsolver.h"


void HLLCSolver::solve()
{
    prepareSolving();
    for(auto i  = 0; i < solParam.MaxIter; i++)
    {
        // устанавливаем текущий временной шаг
        setDt();
        // Вычисляем вектор поточных и релаксационных членов + HLLE
        computeF();
        computeHllcF();
        computeR();

        // Обновляем вектор консервативных переменных
        updateU();
    }
}
