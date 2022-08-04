#pragma once

#include"0_General.hpp"

bool RunOmpSubDomain3D(
  grid& gridA,
  clockstruct& timer,
  simpulationVariable& simu
);

bool Run_OMP_Static_For(
  grid& gridA,
  clockstruct& timer,
  simpulationVariable& simu
);
