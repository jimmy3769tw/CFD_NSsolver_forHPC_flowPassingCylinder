#pragma once
#include "0_General.hpp"

#ifdef OCL_ON


void UpdateT1toT3_OCL(
    OCLstruct& OCL, 
    DfibArray& Dfib,
    simpulationVariable& simu,
    pressure& t1,
    velocity& T1,
    velocity& T3,
    divideLocal& Lo,
    grid& gridA
);


#endif