#pragma once

#include"0_General.hpp"

// Pressure_transform_X_result(t1, Mx, gridA);


void Pressure_transform_X_result(
    pressure &t1, 
    MxClass &Mx ,
    divideLocal& Lo,
    grid & gridA
){
    // #pragma omp parallel firstprivate(Lo)
    {
        int count = 0;
        // #pragma omp for simd
        for (size_t i = Lo.i_begin; i < Lo.i_endof; ++i)
        for (size_t j = Lo.j_begin; j < Lo.j_endof; ++j)
        for (size_t k = Lo.k_begin; k < Lo.k_endof; ++k)
        {
            t1.p[gridA.icel(i,j,k)] = Mx.X_result[gridA.icelCal(i,j,k)];
        }
    }
}



void Pressure_transform_X_result_Dir(
    pressure &t1, 
    MxClass &Mx ,
    divideLocal& Lo,
    grid & gA
){
    // #pragma omp parallel firstprivate(Lo)
    {
        for (size_t i = Lo.i_begin; i < Lo.i_endof; ++i)
        for (size_t j = Lo.j_begin; j < Lo.j_endof; ++j)
        for (size_t k = Lo.k_begin; k < Lo.k_endof; ++k)
        {
            t1.p[gA.icel(i,j,k)] = Mx.X_result[gA.icelDir(i,j,k)];
        }
    }
}




#ifdef EIGEN_ON

void Pressure_transform_x_Eigen(
    pressure &t1, 
    MxClass &Mx ,
    divideLocal& Lo,
    gridArray & gridA
){
    const auto [nx, ny , nz] = gridA.nxyz;

    for (size_t i = Lo.i_begin; i < Lo.i_endof; ++i)
    for (size_t j = Lo.j_begin; j < Lo.j_endof; ++j)
    for (size_t k = Lo.k_begin; k < Lo.k_endof; ++k)
    {
        t1.p[gridA.icel(i,j,k)] = Mx.x_Eigen[gridA.icelCal(i,j,k)];
    }
}

#endif