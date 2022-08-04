#include "0_General.hpp"
#pragma once 
#ifdef OCL_ON

void SorPopeLine_OCL_pre(
    OCLstruct &OCL,
    shareMenory& ShareM,
    SORcoefficient& Sor,
    simpulationVariable& simu,
    velocity& T1,
    pressure& t1,
    divideLocal& Lo,
    grid& gridA
);


void SorPipeLine_OCL(
    OCLstruct &OCL,
    shareMenory& ShareM,
    SORcoefficient& Sor,
    simpulationVariable& simu,
    velocity& T1,
    pressure& t1,
    divideLocal& Lo,
    grid& gridA
);



#endif
