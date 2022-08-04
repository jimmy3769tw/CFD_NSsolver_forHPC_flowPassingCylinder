#pragma once
#include "0_General.hpp"

#ifdef OCL_ON
void UpdateT3toT0_OCL(
    OCLstruct &OCL,
    simpulationVariable& simu,
    velocity& T3,
    velocity& T0,
    divideLocal& Lo,
    grid& gridA
);

#endif