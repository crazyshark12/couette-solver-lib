#include "godunovsolver.h"
#include <iostream>
#include <omp.h>
void GodunovSolver::solve()
{
    prepareSolving();
    writePoints(-1);
    double T = 0;
    for(size_t i  = 0; i < solParam.MaxIter; i++)
    {
        // Устанавливаем текущий временной шаг
        setDt();
        T += timeSolvind.last();

        // FluxF1,2,3
        computeFluxF();
        // Вычисляем вектор релаксационных членов
        //computeR();

        // Обновляем вектор U
        updateU();

        // Обновляем вектор макропараметров
        updatePoints();

        //записать данные, если это требуется
        //writePoints(T*1000000); // микросек
        if(i%1 == 0)
        {
            std::cout<<i<<" iteration"<<std::endl;
            writePoints(T*1); // микросек
        }

        //проверка точности
        if(isObserverWatching)
        {
            // то есть если проверка наблюдателя не пройдена, нужно прекратить рассчёт
            if(!observerCheck(i))
                break;
        }
    }
    writePoints(T*1000000);
}

void GodunovSolver::prepareVectors()
{
    AbstractSolver::prepareVectors();
    rezultAfterPStart.resize(solParam.NumCell - 1);
}

void GodunovSolver::computeFluxF()
{
    #pragma omp parallel for schedule(static)
    for(size_t i =0; i < solParam.NumCell-1; i++)
    {
        rezultAfterPStart[i] = ExacRiemanSolver(points[i],points[i+1],solParam.Gamma);
        rezultAfterPStart[i].densityArray[0] = rezultAfterPStart[i].density;
    }

    #pragma omp parallel for schedule(static)
    for(size_t i =0; i < solParam.NumCell-1; i++)
    {
        auto point = rezultAfterPStart[i];
        double T = point.pressure/(point.density*UniversalGasConstant/point.mixture.molarMass());
        double etta = coeffSolver.shareViscositySimple(point,T);
        double lambda = coeffSolver.lambda(point,T);

        // Рассчитываем производные в точке i
        double dv_dy = (points[i+1].velocity - points[i].velocity) / (delta_h);
        double dT_dy = (points[i+1].temp - points[i].temp) / (delta_h);
        vector<double> dy_dy(mixture.NumberOfComponents);

        //учёт граничных условий
        if(i == 0 || i == solParam.NumCell-1)
            fill(dy_dy.begin(), dy_dy.end(),border.get_dyc_dy());
        else
        {
            for(size_t j = 0 ; j <mixture.NumberOfComponents; j++)
            {
                dy_dy[j] = (points[i+1].fractionArray[j] - points[i].fractionArray[j])/ (delta_h);
            }
        }
        //заполнение вектора потоков
        for(size_t j = 0 ; j <mixture.NumberOfComponents; j++)
        {
            if(j!=0)
                F1[j][i] = -point.density * mixture.getEffDiff(j) * dy_dy[j];
            else
                F1[j][i] = 0;
        }
        F2[i] = -etta * dv_dy;
        F3[i] = 0;
        for(size_t j = 0 ; j <mixture.NumberOfComponents; j++)
        {
            F3[i]+= - point.density * mixture.getEffDiff(j)*dy_dy[j] * mixture.getEntalp(i);
        }
        F3[i] += -lambda*dT_dy - etta*point.velocity*dv_dy;
    }
}

macroParam GodunovSolver::ExacRiemanSolver(macroParam left, macroParam right, double Gamma)
{
    double maxIteration = 40; // макс число итераций
    double TOL=1e-8;
    double lambda = 0; // линия на грани КО
    macroParam ret(mixture);
    ret.fractionArray[0] = 1;

    double left_soundspeed=sqrt ( Gamma*left.pressure/left.density );
    double right_soundspeed=sqrt( Gamma*right.pressure/right.density);

    double p_star= 0.5*(left.pressure+right.pressure) +
            0.125 * ( left.velocity-right.velocity ) *
            ( left.density+right.density ) *
            ( left_soundspeed+right_soundspeed );
    p_star=std::max(p_star,TOL);
    double pMin=std::min(left.pressure,right.pressure);
    double pMax=std::max(left.pressure,right.pressure);

    if ( p_star>pMax )
    {
        double temp1= sqrt ( ( 2.0/ ( Gamma+1.0 ) /left.density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left.pressure ) );
        double temp2= sqrt ( ( 2.0/ ( Gamma+1.0 ) /right.density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right.pressure ) );
        p_star= (temp1*left.pressure+temp2*right.pressure+ ( left.velocity-right.velocity ) ) / ( temp1+temp2 );
        p_star=std::max(p_star,TOL);
    }
    else if ( p_star<pMin )
    {
       double temp1= ( Gamma-1.0 ) / ( 2.0*Gamma );
       p_star= pow(( left_soundspeed+right_soundspeed+0.5*(Gamma-1.0 )*
                   ( left.velocity-right.velocity ) ) /
                   (left_soundspeed/pow(left.pressure,temp1) +
                   right_soundspeed/pow(right.pressure,temp1)), 1.0/temp1);
    }
    double f1 = 0, f2 = 0, f_d = 0 ;
    for(double iteration = 1;iteration < maxIteration; iteration++)
    {
        //LEFT
        double temp1 = sqrt ( ( 2.0/ ( Gamma+1.0 ) /left.density ) / ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left.pressure ) );

        if (p_star<=left.pressure)
            f1=2.0/ ( Gamma-1.0 ) *left_soundspeed*
                    (pow(p_star/left.pressure,(Gamma-1.0 )/(2.0*Gamma))- 1.0) ;
        else
            f1= ( p_star-left.pressure ) *temp1;
        if (p_star<=left.pressure)
            f_d= pow( p_star/left.pressure,-(Gamma+1.0 )/( 2.0*Gamma ))/
                    ( left.density*left_soundspeed );
        else
            f_d=temp1* ( 1.0-0.5* ( p_star-left.pressure ) /
                         ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *left.pressure ) );
        //RIGHT
        temp1 = sqrt ( ( 2.0/ ( Gamma+1.0 ) /right.density ) /
                       ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right.pressure ) );
        if (p_star<=right.pressure)
            f2=2.0/ ( Gamma-1.0 ) *right_soundspeed*
                    (pow(p_star/right.pressure,(Gamma-1.0 )/(2.0*Gamma))- 1.0) ;
        else
            f2= ( p_star-right.pressure ) *temp1;
        if (p_star<=right.pressure)
            f_d= f_d + pow( p_star/right.pressure,-(Gamma+1.0 )/( 2.0*Gamma ))/
                    ( right.density*right_soundspeed );
        else
            f_d=f_d + temp1* ( 1.0-0.5* ( p_star-right.pressure ) /
                         ( p_star+ ( Gamma-1.0 ) / ( Gamma+1.0 ) *right.pressure ) );
        double p_new = p_star - (f1+f2 - (left.velocity - right.velocity))/f_d;
        if(abs(p_new - p_star)/(0.5*abs(p_new + p_star)) < TOL)
            break;
        p_star = p_new;
    }
    // calculate star speed */
    double star_speed=0.5* ( left.velocity + right.velocity ) +0.5* ( f2-f1 );
    double left_star_density, left_tail_speed, left_head_speed,
            right_star_density, right_tail_speed,right_head_speed;
    //LEFT
    if ( p_star>=left.pressure ) {
            // SHOCK
        left_star_density = left.density * ( p_star / left.pressure + ( Gamma-1.0 ) / ( Gamma+1.0 ) ) /
                ( ( Gamma-1.0 ) / ( Gamma+1.0 ) * p_star / left.pressure + 1.0 );
        left_tail_speed = left.velocity -left_soundspeed * sqrt ( ( Gamma+1.0 ) / ( 2.0*Gamma ) * p_star/left.pressure +
                ( Gamma-1.0 ) / ( 2.0*Gamma ) );
        left_head_speed = left_tail_speed;
    }
    else // % left_wave_ == kRarefaction
    {
        left_star_density = left.density * pow(p_star/left.pressure,1.0/Gamma);
        left_head_speed = left.velocity - left_soundspeed;
        left_tail_speed = star_speed - sqrt ( Gamma*p_star/left_star_density );
    }
    //RIGHT
    if ( p_star>=right.pressure )
    {
        right_star_density = right.density *
                            ( p_star / right.pressure + ( Gamma-1.0 ) / ( Gamma+1.0 ) ) /
                            ( ( Gamma-1.0 ) / ( Gamma+1.0 ) * p_star / right.pressure + 1.0 );
        right_tail_speed = right.velocity +
           right_soundspeed * sqrt ( ( Gamma+1.0 ) / ( 2.0*Gamma ) * p_star/right.pressure +
           ( Gamma-1.0 ) / ( 2.0*Gamma ) );
        right_head_speed = right_tail_speed;
    }
    else // % right_wave_ == kRarefaction
    {
        right_star_density = right.density *  pow(p_star/right.pressure, 1.0/Gamma );
        right_head_speed = right.velocity + right_soundspeed;
        right_tail_speed = star_speed + sqrt ( Gamma*p_star/right_star_density );
    }

    bool is_left_of_contact = lambda  < star_speed;

    if ( is_left_of_contact )
    {// % the u is left of contact discontinuity
        if ( p_star>=left.pressure )  //the left wave is a shock
        {
            if ( lambda < left_head_speed )
            { // the u is before the shock
                ret.density  = left.density;
                ret.velocity = left.velocity;
                ret.pressure = left.pressure;
            }
            else  //% the u is behind the shock
            {
                ret.density  = left_star_density;
                ret.velocity = star_speed;
                ret.pressure = p_star;
            }
        }
        else // % the left wave is a rarefaction
        {
            if ( lambda < left_head_speed )//  % the u is before the rarefaction
            {
                ret.density  = left.density;
                ret.velocity = left.velocity;
                ret.pressure = left.pressure;
            }
            else
            {
                if ( lambda < left_tail_speed )//  % the u is inside the rarefaction
                {//% left_rarefaction (4.56)}
                    double temp1 = 2.0/ ( Gamma+1.0 ) + ( Gamma-1.0 ) / ( Gamma+1.0 )/left_soundspeed *(left.velocity - lambda);
                    ret.density = left.density *  pow(temp1, 2.0/( Gamma-1.0 ));
                    ret.pressure = left.pressure * pow(temp1, 2.0*Gamma/ ( Gamma-1.0));
                    ret.velocity = 2.0/ ( Gamma+1.0 ) * ( left_soundspeed + ( Gamma-1.0 ) /2.0*left.velocity + lambda);
                }
                else//  % the u is after the rarefaction
                {
                    ret.density  = left_star_density;
                    ret.velocity = star_speed;
                    ret.pressure = p_star;
                }
            }
        }
    }
    else// % the queried u is right of contact discontinuity
        //%------------------------------------------------------------------------
    {
        if ( p_star>=right.pressure )  //% the right wave is a shock
        {
            if ( lambda > right_head_speed )  //% the u is before the shock
            {
                ret.density  = right.density;
                ret.velocity = right.velocity;
                ret.pressure = right.pressure;
            }
            else  //% the u is behind the shock
            {
                ret.density  = right_star_density;
                ret.velocity = star_speed;
                ret.pressure = p_star;
            }
        }
        else // % the right wave is a rarefaction
        {
            if ( lambda > right_head_speed ) // % the u is before the rarefaction
            {
                ret.density  = right.density;
                ret.velocity = right.velocity;
                ret.pressure = right.pressure;
            }
            else
            {
                if ( lambda > right_tail_speed ) // % the u is inside the rarefaction
                {
                    double temp1 =2.0/ ( Gamma+1.0 ) - ( Gamma-1.0 ) / ( Gamma+1.0 ) /right_soundspeed *(right.velocity - lambda);
                    ret.density = right.density *  pow(temp1, 2.0/ ( Gamma-1.0 ) );
                    ret.pressure = right.pressure * pow(temp1, 2.0*Gamma/ ( Gamma-1.0 ) );
                    ret.velocity = 2.0/ ( Gamma+1.0 ) * ( -right_soundspeed + ( Gamma-1.0 ) /2.0*right.velocity + lambda);
                }
                else // % the u is after the rarefaction
                {
                    ret.density  = right_star_density;
                    ret.velocity = star_speed;
                    ret.pressure = p_star;
                }
            }
        }
    }
    return ret;
}

void GodunovSolver::updateU()
{
    for(auto i  = 1; i < solParam.NumCell-1; i++)
    {
        for (int j = 0; j < mixture.NumberOfComponents; j++)
        {
            U1[j][i] += (/*R[j][i]*/0 - (F1[j][i] - F1[j][i - 1]) / delta_h) * timeSolvind.last();
        }
        U2[i] += (0 - (F2[i] - F2[i - 1]) / delta_h) * timeSolvind.last();
        U3[i] += (0 - (F3[i] - F3[i - 1]) / delta_h) * timeSolvind.last();
    }
}

void GodunovSolver::updatePoints()
{
    auto size = points.size()-1;
    #pragma omp parallel for schedule(static)
    for(size_t i = 1; i < size; i++)
    {
        points[i].velocity = U2[i]/U1[0][i];
        points[i].density = U1[0][i];
        points[i].pressure = (U3[i] - pow(points[i].velocity,2)*0.5*U1[0][i])*(solParam.Gamma - 1);

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
