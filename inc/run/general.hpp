#pragma once

#include <random>

#include <stdio.h>      /* printf, scanf, puts, NULL */

#include <stdlib.h>     /* srand, rand */

#include <ctime>       /* time */


#include "profiling/Chronograph.hpp"

#include "resize_and_simuSetting.hpp"

#include "grid/Generator.hpp"

// --------------------------------
#include "analyses/cfl.hpp"

#include "analyses/CentralProfile.hpp"

#include "analyses/CheckL2norm.hpp"

#include "analyses/CheckSteadystate.hpp"
// --------------------------------

// --------------------------------
#include "source/EngeryEquation.hpp"

#include "source/Convection_Difussion.hpp"

#include "source/Convection_Difussion.uniformGrid.hpp"

#include "source/SmagorinskyModel.hpp"

#include "source/copyVel.hpp"
// --------------------------------

#include "InitionCondition.hpp"

#include "pressure/Transform.hpp"

// --------------------------------
#include "io/Plot3D/qfile.hpp"

#include "io/Plot3D/xfileWrite.hpp"

#include "io/Plot3D/qfileRead.hpp"
// --------------------------------

// --------------------------------
#include "matrix/PressureMatrix.hpp"

#include "matrix/ELL_sparseMatrix.hpp"

#include "matrix/SPE_sparseMatrix0.hpp"

#include "matrix/SPE_sparseMatrix.hpp"

#include "matrix/b_matrix.hpp"
// --------------------------------

// --------------------------------
#include "dfib/setEta.hpp"

#include "dfib/virtualForceIntergrator.hpp"

#include "dfib/updateUandF.hpp"
// --------------------------------

// --------------------------------
#include "solver/bicgstab.hpp"

#include "solver/bicgstab0.hpp"

#include "solver/bicgstabRe.hpp"

#include "solver/bicgstabRe2.hpp"
// --------------------------------


// --------------------------------
#include "BoundaryCondition/Boundary_Condition.hpp"

#include "BoundaryCondition/updateBC.hpp"
// --------------------------------



// --------------------------------
#include "pressure/solver/SorPipeLine.omp.hpp"

#include "pressure/solver/SorPipeLine.seq.hpp"
// --------------------------------


