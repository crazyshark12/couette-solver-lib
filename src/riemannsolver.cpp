#include "riemannsolver.h"

#include <algorithm>
const double gamma = 1.4;
void HLLCSolver::computeFlux(SystemOfEquation *system)
{
    for(size_t i = 0 ; i < system->numberOfCells-1; i++)
    {
        if(i==55)
            double dg = 0;
        double u0, u1,v0,v1,a0,a1, rho0, rho1, p0, p1, E0, E1, H0 , H1, avg_H, S0 , S1, S_star;
        vector<double> U_star_0(system->systemOrder),U_star_1(system->systemOrder);
        double avg_a, avg_v;

        // тут u - нормальная составляющая, v - касательная
        u0 = system->getVelocityNormal(i);
        u1 = system->getVelocityNormal(i+1);

        v0 = system->getVelocityTau(i);
        v1 = system->getVelocityTau(i+1);

        rho0 = sqrt(system->getDensity(i));
        rho1 = sqrt(system->getDensity(i+1));

        avg_v = (rho0 * u0 + rho1 * u1) / (rho0 + rho1);
//        H0 = (3*points[i].pressure)/(2*U1[0][i]);
//        H1 = (3*points[i+1].pressure)/(2*U1[0][i+1]);
        E0 = system->getEnergy(i);
        E1 = system->getEnergy(i+1);

        p0 = system->getPressure(i);
        p1 = system->getPressure(i+1);

        H0 = E0 + p0/system->getDensity(i);
        H1 = E1 + p1/system->getDensity(i+1);
        avg_H = (rho0 * H0 + rho1 * H1) / (rho0 + rho1);
        avg_a = sqrt((gamma - 1.)*(avg_H - 0.5 * pow(avg_v,2)));
        S0 = (avg_v - avg_a);
        S1 = (avg_v + avg_a);
//        S0 = min(v0, v1);
//        S1 = max(v0, v1);
        S_star = (p1 - p0 +
                pow(rho0,2)*u0*(S0 - u0) - pow(rho1,2)*u1*(S1 - u1))
                / (pow(rho0,2)*(S0 - u0) - pow(rho1,2)*(S1 - u1));

        double coeff_0 = system->getDensity(i) * ((S0 - u0)/ (S0 - S_star));
        double coeff_1 = system->getDensity(i+1) * ((S1 - u1)/ (S1 - S_star));
        for(size_t j = 0; j < system->numberOfComponents; j++)
        {
            U_star_0[j] = coeff_0;
            U_star_1[j] = coeff_1;
        }
//        U_star_0[system->v_tau] = coeff_0 * v0;
//        U_star_1[system->v_tau] = coeff_1 * v1;

        U_star_0[system->v_normal] = system->getDensity(i)   * coeff_0 * S_star;
        U_star_1[system->v_normal] = system->getDensity(i+1) * coeff_1 * S_star;

//        U_star_0[system->energy] = coeff_0 * (E0 + (S_star - u0)*(S_star + p0/(system->getDensity(i) * (S0 - u0))));
//        U_star_1[system->energy] = coeff_1 * (E1 + (S_star - u1)*(S_star + p1/(system->getDensity(i+1) * (S1 - u1))));


        U_star_0[system->energy] = coeff_0 * (E0/system->getDensity(i) + (S_star - u0)*(S_star + p0/(system->getDensity(i) *(S0 - u0))));
        U_star_1[system->energy] = coeff_1 * (E1/system->getDensity(i+1) + (S_star - u1)*(S_star + p1/(system->getDensity(i+1)*(S1 - u1))));

        if(S0 >= 0)
        {
            for(size_t j = 0; j < system->systemOrder; j++)
            {
                system->Flux[j][i] = system->F[j][i];
            }
        }
        else if(S1 <= 0)
        {
            for(size_t j = 0; j < system->systemOrder; j++)
            {
                system->Flux[j][i] = system->F[j][i+1];
            }
        }
//        else if(S0<=0 && S1>=0)
//        {
//            for(size_t j = 0; j < mixture.NumberOfComponents; j++)
//            {
//                hllc_f1[j] = (S1 * F1[j][i] - S0 * F1[j][i+1] + S0*S1*(U1[j][i+1] - U1[j][i]))/(S1 - S0);
//            }
//            hllc_f2 = (S1 * F2[i] - S0 * F2[i+1] + S0*S1*(U2[i+1] - U2[i]))/(S1 - S0);
//            hllc_f3 = (S1 * F3[i] - S0 * F3[i+1] + S0*S1*(U3[i+1] - U3[i]))/(S1 - S0) ;
//        }
        else if( S0 <= 0 && S_star >= 0)
        {
            for(size_t j = 0; j < system->systemOrder; j++)
            {
                    system->Flux[j][i] = system->F[j][i] + S0*(U_star_0[j] - system->U[j][i]);
            }
        }
        else if( S_star <= 0 && S1 >= 0)
        {
            for(size_t j = 0; j < system->systemOrder; j++)
            {
                    system->Flux[j][i] = system->F[j][i+1] + S1*(U_star_1[j] - system->U[j][i+1]);
            }
        }
    }
    return;
}

void HLLESolver::computeFlux(SystemOfEquation *system)
{
    double gamma = 1.4;
    for(size_t i = 0 ; i < system->numberOfCells-1; i++)
    {
        double H0, H1, c0, c1, u0, u1, v0,v1,V0,V1, rho0, rho1, u_avg,v_avg, H_avg, c_avg, b0, b1, b_plus, b_minus;

        u0 = system->getVelocityNormal(i);
        u1 = system->getVelocityNormal(i+1);

        v0 = system->getVelocityTau(i);
        v1 = system->getVelocityTau(i+1);

        V0 = sqrt(pow(u0,2) + pow(v0,2));
        V1 = sqrt(pow(u1,2) + pow(v1,2));

//        H0 = (5*points[i].pressure)/(2*U1[0][i]) + pow(V0,2)/2;
//        H1 = (5*points[i+1].pressure)/(2*U1[0][i+1])+ pow(V1,2)/2;
        H0 = system->getEnergy(i) + system->getPressure(i)/system->getDensity(i);
        H1 = system->getEnergy(i+1) + system->getPressure(i+1)/system->getDensity(i+1);

        c0 = sqrt((gamma - 1)*(H0 - 0.5 * pow(V0,2))); // TODO 5/3 = gamma
        c1 = sqrt((gamma - 1)*(H1 - 0.5 * pow(V1,2)));

        rho0 = sqrt(system->getDensity(i));
        rho1 = sqrt(system->getDensity(i+1));

        u_avg = (rho0 * u0 + rho1 * u1) / (rho0 + rho1);
        v_avg = (rho0 * v0 + rho1 * v1) / (rho0 + rho1);

        H_avg = (rho0 * H0 + rho1 * H1) / (rho0 + rho1);
        c_avg = sqrt((gamma)*(H_avg - 0.5 * (pow(u_avg,2) + pow(v_avg,2)))); // TODO 5/3 = gamma

        b0 = (std::min)({u_avg - c_avg, u0 - c0});
        b1 = (std::max)({u_avg + c_avg, u1 + c1});

        b_plus = (std::max)({0., b1});
        b_minus = (std::min)({0., b0});

        for(size_t j = 0; j < system->systemOrder; j++)
        {
            system->Flux[j][i] = (b_plus * system->F[j][i] - b_minus* system->F[j][i+1])/(b_plus  - b_minus) + (b_plus * b_minus)/(b_plus - b_minus) * (system->U[j][i+1] - system->U[j][i]);
        }
    }
}


void HLLSimple::computeFlux(SystemOfEquation *system, double dt, double dh)
{
    for(size_t i = 0 ; i < system->numberOfCells - 1; i++)
    {
        double SR, SL, FL, FR, UL, UR;
        SR = dh/dt;
        SL = -dh/dt;
        if(i == 49)
            double x  = 3 ;
        for(size_t j = 0; j < system->systemOrder; j++)
        {
            FR = system->F[j][i+1];
            FL = system->F[j][i];
            UR = system->U[j][i+1];
            UL = system->U[j][i];
            system->Flux[j][i] = (SR*FL - SL*FR + SL*SR * (UR - UL))/(SR - SL);
        }
    }
}

void HLLIsentropic::computeFlux(SystemOfEquation *system)
{
    for(size_t i = 0 ; i < system->numberOfCells - 1; i++)
    {
        double u_star, a_star,ul,ur,al,ar,SR, SL, FL, FR, UL, UR;
        ul = system->getVelocity(i);
        ur = system->getVelocity(i+1);
        al = system->getSoundSpeed(i);
        ar = system->getSoundSpeed(i+1);
        u_star = 1/2*(ul-ur) + (al - ar)/(gamma - 1);
        a_star = 1/2*(al + ar) + 1/4 *(gamma - 1)*(ul-ur);
        SL = std::min({ul - al, u_star - a_star});
        SR = std::max({ur + ar, u_star + a_star});
        for(size_t j = 0; j < system->systemOrder; j++)
        {
            FR = system->F[j][i+1];
            FL = system->F[j][i];
            UR = system->U[j][i+1];
            UL = system->U[j][i];
            if(SL >= 0 )
                system->Flux[j][i] = FL;
            else if( SR <= 0)
                system->Flux[j][i] = FR;
            else
                system->Flux[j][i] = (SR*FL - SL*FR + SL*SR * (UR - UL))/(SR - SL);
        }
    }
}

void ExacRiemanSolver::computeFlux(SystemOfEquation *system)
{
    for(size_t i = 0 ; i < system->numberOfCells - 1; i++)
    {
        macroParam left,right,res;
        left.density = system->getDensity(i);
        right.density = system->getDensity(i+1);
        left.pressure = system->getPressure(i);
        right.pressure = system->getPressure(i+1);
        left.velocity = system->getVelocity(i);
        right.velocity= system->getVelocity(i+1);

        res = exacRiemanSolver(left,right, gamma);

        //тут конкретная реализация под задачу сода, ибо этот момент я видимо не продумал архитектуру,
        //надо видимо делать для точного метода какой-то особый случай, ну или здесь внутри прописывать if-ы и в системе ввести поле которое сможет сказать о типе
        //указателя базового класса, чтобы можно было определять тип объекта по указателю

        system->Flux[0][i] = res.density*res.velocity;
        system->Flux[1][i] = res.density*pow(res.velocity,2) +  res.pressure;
        double rho_e = res.pressure/(gamma - 1);
        double kinetic = res.density * pow(res.velocity,2) / 2. ;
        double E = kinetic + rho_e;
        system->Flux[2][i] =  res.velocity*(E + res.pressure);
    }
}

macroParam ExacRiemanSolver::exacRiemanSolver(macroParam left, macroParam right, double Gamma)
{
    double maxIteration = 40; // макс число итераций
    double TOL=1e-8;
    double lambda = 0; // линия на грани КО
    macroParam ret(left.mixture);;
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

void HLLESolverSimen::computeFlux(SystemOfEquation *system)
{
    for(size_t i = 0; i < system->numberOfCells-1; i++)
    {
        // Временные переменные
        double f_hlle, f0, f1, u0, u1, a0, a1, v0, v1, rho0, rho1, avg_k;
        double b0, b1, eta, avg_a, avg_v;

        // Забираем известные макропараметры в (.)-ах [i], [i+1]
        a0 = pow(system->getSoundSpeed(i),2);
        a1 = pow(system->getSoundSpeed(i+1),2);
        v0 = system->getVelocity(i);
        v1 = system->getVelocity(i+1);
        rho0 = sqrt(system->getDensity(i));
        rho1 = sqrt(system->getDensity(i+1));
        avg_v = (rho0 * v0 + rho1 * v1) / (rho0 + rho1);
        avg_k = gamma;  // у cИмена было так: 0.5 * (k[i] + k[i + 1]);

        for (size_t j = 0; j < system->systemOrder; j++)
        {
            // Забираем известные макропараметры в (.)-ах [i], [i+1]
            f0 = system->F[j][i];
            f1 = system->F[j][i + 1];
            u0 = system->U[j][i];
            u1 = system->U[j][i + 1];

            // Расчитываем сигнальные скорости
            eta = (avg_k - 1.0) / 2.0 * rho0 * rho1 / pow(rho0 + rho1, 2.0);
            avg_a = sqrt((rho0 * a0 + rho1 * a1) / (rho0 + rho1) +
                          eta * pow(v1 - v0, 2.0));
            avg_v = (rho0 * v0 + rho1 * v1) / (rho0 + rho1);

            // Расчет потока на стыке ячеек
            b0 = std::min(avg_v - avg_a, 0.0);
            b1 = std::max(avg_v + avg_a, 0.0);
            f_hlle = (b1 * f0 - b0 * f1 + b1 * b0 * (u1 - u0)) / (b1 - b0);

            // Обновляем значение [j] вектора поточных членов в (.) [i-1]
            system->Flux[j][i] = f_hlle;
        }
    }
}
