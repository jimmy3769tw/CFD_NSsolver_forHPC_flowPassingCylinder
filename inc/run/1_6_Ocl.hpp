#pragma once 
#include "0_General.hpp"

// ! for Boost Compute ----------------
#ifdef OCL_ON

bool RunOcl(
    gridArray &gridA,
    clockstruct &timer,
    simpulationVariable &simu
);


bool RunOcl_Debug(
    gridArray &gridA,
    clockstruct &timer,
    simpulationVariable &simu
);

#endif