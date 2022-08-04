#pragma once
#include"0_General.hpp"

// *1 run
#include"run/1_1_Sequential.hpp"

#ifdef OMP_ON
#include"run/1_2_Omp.hpp"
#endif

#ifdef MPI_ON
#include"run/1_3_Mpi.hpp"
#endif

#ifdef OCL_ON
#include"run/1_6_Ocl.hpp"
#endif

// #include"1_4_OmpHybiridMpi.hpp"