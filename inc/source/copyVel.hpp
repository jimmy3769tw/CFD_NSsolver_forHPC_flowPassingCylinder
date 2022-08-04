#pragma once
#include <0_General.hpp>

void copyVel(
    const velocity& oldVel,
    velocity& curtVel
)
{

    // #if defined (PC_SEQ)

    curtVel.u = oldVel.u;
    curtVel.v = oldVel.v;    
    curtVel.w = oldVel.w;

    // #elif defined (PC_OMP)

    // #pragma omp parallel
    // {
    //     #pragma omp for simd nowait
    //     for(size_t i = 0; i< curtVel.u.size(); i++)
    //     {
    //         curtVel.u[i] = oldVel.u[i];
    //     }

    //     #pragma omp for simd nowait
    //     for(size_t i = 0; i< curtVel.v.size(); i++)
    //     {
    //         curtVel.v[i] = oldVel.v[i];
    //     }

    //     #pragma omp for simd 
    //     for(size_t i = 0; i< curtVel.w.size(); i++)
    //     {
    //         curtVel.w[i] = oldVel.w[i];
    //     }

    // }

//     #endif

}