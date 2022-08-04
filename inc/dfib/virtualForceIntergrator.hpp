#pragma once

#include"0_General.hpp"

auto virtualF_IntbySimpson(
    DfibArray& Dfib,
    divideLocal& Lo,
    grid& gA,
    size_t whichDirection
)
{

    Dfib.Val_sumz.resize(gA.nxny);

    Dfib.Val_sumz_sumy.resize(gA.nx);

    Dfib.ValSum.resize(3);


    // BoundaryCondtion for 
    // * ----------------------- z -------------------------
        
        for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
        {
            Dfib.Val_sumz[gA.icel_z(i,j)] = 0.0;
            for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
            {
                Dfib.Val_sumz[gA.icel_z(i,j)] += (
                    4.0 *Dfib.f[gA.icel(i,j,k)  + gA.iceltot*whichDirection]  * gA.Dzs[k]
                        +Dfib.f[gA.icel(i,j,k-1)+ gA.iceltot*whichDirection]  * gA.Dzs[k-1]
                        +Dfib.f[gA.icel(i,j,k+1)+ gA.iceltot*whichDirection]  * gA.Dzs[k+1]
                )/6.;
            }
        }

        // * ----------------------- y -------------------------

        for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        {
            if (Dfib.Val_sumz_sumy[i] != 0.0)
            for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
            {
                Dfib.Val_sumz[gA.icel_z(i,j)] +=(
                    4.0 *Dfib.Val_sumz[gA.icel_z(i,j)]  * gA.Dys[j]
                        +Dfib.Val_sumz[gA.icel_z(i,j-1)]* gA.Dys[j-1]
                        +Dfib.Val_sumz[gA.icel_z(i,j+1)]* gA.Dys[j+1]
                ) /6.;
            }
        }
        // * ----------------------- x -------------------------

        Dfib.ValSum[whichDirection] = 0.;
        for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
        {
            if (Dfib.Val_sumz_sumy[i] != 0.0)
            Dfib.ValSum[whichDirection] +=(
                4.0 *Dfib.Val_sumz_sumy[i]*   gA.Dxs[i] 
                    +Dfib.Val_sumz_sumy[i-1]* gA.Dxs[i-1]
                    +Dfib.Val_sumz_sumy[i+1]* gA.Dxs[i+1]
            ) /6.;
        }
}



auto virtualF_Int(
    DfibArray& Dfib,
    divideLocal& Lo,
    grid& gA,
    size_t whichDirection
)
{
    if (Dfib.ValSum.size() < 3)
        Dfib.ValSum.resize(3);

    Dfib.ValSum[whichDirection]  = 0.0;

    for (int i = Lo.i_begin ; i < Lo.i_endof ; ++i )
    for (size_t j = Lo.j_begin ; j < Lo.j_endof ; ++j )
    for (int k = Lo.k_begin ; k < Lo.k_endof ; ++k )
    {
        Dfib.ValSum[whichDirection] +=
            Dfib.f[gA.icel(i,j,k)  + gA.iceltot*whichDirection] 
                                    * (gA.Dx[i]*gA.Dy[j]*gA.Dz[k]) ;
    }

    return 0;
}


auto CD_CL(
    simpulationVariable& simu,
    divideLocal& Lo,
    DfibArray& Dfib,
    grid& gA
)
{

    virtualF_Int(Dfib, Lo, gA, 0);

    virtualF_Int(Dfib, Lo, gA, 1);

    auto area_cD    = gA.lz;
    auto area_cL    = gA.lz;

    auto cD = -2.0 * Dfib.ValSum[0] / area_cD;
 
    auto cL = -2.0 * Dfib.ValSum[1] / area_cL;

    cout << "[cD, cL] : " << cD << ", " << cL << endl; 

    std::ofstream file;

    std::string name = "Information/Time_cDcL";

    name += ".dat";

    file.open (name, std::ios::out|ios::app);
    std::string tab = " ";

    if (simu.loop == 1 ){

        std::vector<std::string> variables;
        variables.push_back("simulation time");
        variables.push_back("C<sub>D</sub>");
        variables.push_back("C<sub>L</sub>");

            file 
            << "TITLE     = \"\"\n"
            << "VARIABLES = \""
            << variables.at(0)
            << "\",\""
            << variables.at(1)
            << "\",\""
            << variables.at(2)
            << "\"\n"
            << "ZONE T=\""
            << simu.ZONE()
            << "\"";
    }

    file << "\n"<<  simu.time << tab << cD << tab  << cL ;

    file.close();

    return true;
}
