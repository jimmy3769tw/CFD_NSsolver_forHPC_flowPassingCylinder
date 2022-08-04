#pragma once 

#include "0_General.hpp"

#ifdef OCL_ON
void CenterQuickScheme_OCL(
    OCLstruct &OCL, 
    sourceTerm& So,
    simpulationVariable& simu,
    velocity& T0,
    velocity& T1,
    divideLocal& Lo,
    grid& gridA
);

#endif
