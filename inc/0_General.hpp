
#pragma once

#include "BoundaryCondition/P_NEUMANN.hpp"

// ! Possion equation ------------

// #define P_SOLVER_SOR

// #define P_SOLVER_BICG_CSR

#define P_SOLVER_BICG_SPE

// #define P_SOLVER_AMGCL_BUILTIN

// #define P_SOLVER_BICG_ELL

// #define P_SOLVER_CG_ELL


// #define P_SOLVER_AMGCL_EIGEN

// #define P_SOLVER_EIGEN_CSR


// ---------------------------


#if defined (P_SOLVER_BICG_ELL)

    #define P_SOLVER_ELL

#elif defined (P_SOLVER_BICG_SPE)

    #define P_SOLVER_SPE

#elif defined (P_SOLVER_CG_ELL)

    #define P_SOLVER_ELL

#elif defined (P_SOLVER_BICG_CSR)

    #define P_SOLVER_CSR

#endif
// ---------------------------




// ! Possion equation ------------


// ! Parallel computing strategy ------------
// #define PC_SEQ
#define PC_OMP
// #define PC_HYBRID_MPI_OMP
// #define PC_OCL

/* 
    TODO :  Validation
    TODO :- [X] 0. PC_SEQ
    TODO :- [X] 2. PC_OMP
    TODO :- [X] 3. PC_HYBRID_MPI_OMP
    TODO :- [ ] 4. PC_OCL
*/

// #define MPI_DEBUG



// *-----------------

#if defined (PC_OCL)
    #define OCL_ON
    // #define OCL_DEBUG
    // #define BoostCompute_ON

#elif defined (PC_HYBRID_MPI_OMP)

    #define OMP_ON
    #define MPI_ON

#elif defined (PC_OMP)

    #define OMP_ON

#endif


// *-----------------
// ! Parallel computing strategy ------------



//! Temporal discretization ------------
#define TEMPORAL_DISCRETIZATION_1_ORDER

// #define TEMPORAL_DISCRETIZATION_2_ORDER

// #define TEMPORAL_DISCRETIZATION_3_ORDER

/* 
TODO :  Validation  
TODO :- [X] TEMPORAL_DISCRETIZATION_1_ORDER
TODO :- [ ] TEMPORAL_DISCRETIZATION_2_ORDER
TODO :- [ ] TEMPORAL_DISCRETIZATION_3_ORDER
*/
//! Temporal discretization ------------


//! convection and difussion ------------

#define CONVECTION_DIFUSSION_QUICK

// #define CONVECTION_DIFUSSION_LUD

// #define CONVECTION_DIFUSSION_UD

//* convection and difussion ------------


// ! Terbulence module ------------
// #define TERBULENCE_SMAGORINSKY
/* 
    TODO :  Validation
    TODO :- [ ] OFF
    TODO :- [ ] TERBULENCE_SMAGORINSKY
*/
// ! Terbulence module ------------




// !DFIB_Cylinder
/* 
    TODO :  Validation
    TODO :- [X] OFF
    TODO :- [X] DFIB_Cylinder-Z
    TODO :- [ ] DFIB_Cylinder-X
    TODO :- [ ] DFIB_Cylinder-Y
*/



#ifdef P_SOLVER_AMGCL_BUILTIN
    #define AMGCL_ON
#endif


#ifdef P_SOLVER_AMGCL_EIGEN
    #define AMGCL_ON
    #define EIGEN_ON
#endif


#ifdef P_SOLVER_EIGEN_CSR
    #define EIGEN_ON
#endif

#ifdef OMP_ON
    #include <omp.h>
#endif


#ifdef EIGEN_ON
    #include "import/Eigen.hpp"
#endif


#ifdef AMGCL_ON
    #include "import/Amgcl.hpp"
#endif


#ifdef OCL_ON
#include "4_OCL.hpp"
#endif

// ! HEADFILE  --------vvvvvvvv

#include "0_0_PhisicalVaribal.hpp"

#include "0_1_SimulationVariable.hpp"

#include "profiling/STL_clock.hpp"

#include "0_3_Domain.hpp"

#include "grid/structureGrid.hpp"

// --------------------------------------------------

#include "0_7_PossionEquationStruct.hpp"

#include "import/stl.hpp"

#include <dirent.h>