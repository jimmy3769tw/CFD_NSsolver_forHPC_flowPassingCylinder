#pragma once

#include"0_General.hpp"

void resize_variable(
    grid& gridA,
    simpulationVariable& simu,
    pressure& t1,
    velocity& T0,
    velocity& T1,
    velocity& T3,
    DfibArray& Dfib
)
{

    std::cout << "Resize Array iceltot: " << gridA.iceltot 
        << ", gridA.iceltotCal :" << gridA.iceltotCal << std::endl;

    Dfib.eta.resize(gridA.iceltot);

    Dfib.f.resize(gridA.iceltot*3);

    t1.p.resize(gridA.iceltot);

    T0.resize(gridA.iceltot);

    T1.resize(gridA.iceltot);

    T3.resize(gridA.iceltot);

    std::cout << "Resize Array finish. \n" << std::endl;
}
