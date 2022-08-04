// #pragma OPENCL EXTENSION cl_khr_global_int32_base_atomic:enable
// #pragma OPENCL EXTENSION cl_khr_local_int32_base_atomic:enable

__kernel
void UpdateT3toT0(
    __global double *T0_u,
    __global double *T0_v,
    __global double *T0_w,
    __global const double *T3_u,
    __global const double *T3_v,
    __global const double *T3_w,
    const uint nx,const uint ny,const uint nz,
    const uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);

    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel = i*nz*ny + j*nz + k; 

    T0_u[icel] = T3_u[icel] ;
    T0_v[icel] = T3_v[icel] ;
    T0_w[icel] = T3_w[icel] ;
}
