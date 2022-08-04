#pragma once

#include"0_General.hpp"

bool runSeqOmp(
  grid& gridA,
  clockstruct& timer,
  simpulationVariable& simu,
  int argc, char **argv
);



bool debug_dwan_BC_and_quick(
    grid &gridA,
    clockstruct &timer,
    simpulationVariable &simu,
    int argc, char **argv
);



bool debug_dwan_SOR(
    grid &gridA,
    clockstruct &timer,
    simpulationVariable &simu,
    int argc, char **argv
);
