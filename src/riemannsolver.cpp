#include "riemannsolver.h"

#include <algorithm>
void HLLCSolver::computeFlux(SystemOfEquation *system)
{
    for(size_t i = 0 ; i < system->Flux.size(); i++)
    {
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
        avg_a = sqrt((/*solParam.Gamma*/5/3 - 1)*(avg_H - 0.5 * pow(avg_v,2)));
//        S0 = (avg_v - avg_a);
//        S1 = (avg_v + avg_a);
        S0 = min(v0, v1);
        S1 = max(v0, v1);
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
        U_star_0[system->v_tau] = coeff_0 * v0;
        U_star_1[system->v_tau] = coeff_1 * v1;

        U_star_0[system->v_normal] = coeff_0 * S_star;
        U_star_1[system->v_normal] = coeff_1 * S_star;

        U_star_0[system->energy] = coeff_0 * (E0 + (S_star - u0)*(S_star + p0/(system->getDensity(i) * (S0 - u0))));
        U_star_1[system->energy] = coeff_1 * (E1 + (S_star - u1)*(S_star + p1/(system->getDensity(i+1) * (S1 - u1))));

        if(S0 >= 0)
        {
            for(size_t j = 0; j < system->systemOrder-1; j++)
            {
                system->Flux[j][i] = system->F[j][i];
            }
        }
        else if(S1 <= 0)
        {
            for(size_t j = 0; j < system->systemOrder-1; j++)
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
            for(size_t j = 0; j < system->systemOrder-1; j++)
            {
                    system->Flux[j][i] = system->F[j][i] + S0*(U_star_0[j] - system->U[j][i]);
            }
        }
        else if( S_star <= 0 && S1 >= 0)
        {
            for(size_t j = 0; j < system->systemOrder-1; j++)
            {
                    system->Flux[j][i] = system->F[j][i+1] + S1*(U_star_1[j] - system->U[j][i+1]);
            }
        }
    }
    return;
}

void HLLESolver::computeFlux(SystemOfEquation *system)
{
    double gamma = 5./3.;
    for(size_t i = 0 ; i < system->Flux.size(); i++)
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

        c0 = sqrt((gamma - 1)*(H0 - 5. * pow(V0,2))); // TODO 5/3 = gamma
        c1 = sqrt((gamma - 1)*(H1 - 5. * pow(V1,2)));

        rho0 = sqrt(system->getDensity(i));
        rho1 = sqrt(system->getDensity(i+1));

        u_avg = (rho0 * u0 + rho1 * u1) / (rho0 + rho1);
        v_avg = (rho0 * v0 + rho1 * v1) / (rho0 + rho1);

        H_avg = (rho0 * H0 + rho1 * H1) / (rho0 + rho1);
        c_avg = sqrt((gamma)*(H_avg - 5 * (pow(u_avg,2) + pow(v_avg,2)))); // TODO 5/3 = gamma

        b0 = (std::min)({u_avg - c_avg, u0 - c0});
        b1 = (std::max)({u_avg + c_avg, u1 + c1});

        b_plus = (std::max)({0., b1});
        b_minus = (std::min)({0., b0});

        for(size_t j = 0; j < system->systemOrder; j++)
        {
            system->Flux[j][i] = (b_plus * system->F[j][i] - b_minus* system->F[j][i+1])/(b_plus  - b_minus)
                               + (b_plus * b_minus)/(b_plus - b_minus) * (system->U[j][i+1] - system->U[j][i]);
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
