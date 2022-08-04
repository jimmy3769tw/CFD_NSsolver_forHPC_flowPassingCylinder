// #pragma OPENCL EXTENSION cl_khr_global_int32_base_atomic:enable
// #pragma OPENCL EXTENSION cl_khr_local_int32_base_atomic:enable

__kernel
void UpdateT1toT3(
    __global const double *T1_u,
    __global const double *T1_v,
    __global const double *T1_w,
    __global double *T3_u,
    __global double *T3_v,
    __global double *T3_w,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global const double *pressure,
    __global double *Dfib_f,
    __global const double *Dfib_eta,
    const double u_solid,
    const double v_solid,
    const double w_solid,
    const double dt,
    const double nu,
    const uint nx, const uint ny,const uint nz,const uint gCells, const uint iceltot
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);


    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;

    uint icel = i*nz*ny + j*nz + k;



    // ! ---------------  x  ---------------
    double T2u = T1_u[icel] - dt*( pressure[icel+(nz*ny)]-pressure[icel] ) / Dxs[i];

    T3_u[icel]= Dfib_eta[icel] * u_solid 
        + (1.0-0.5*(
        Dfib_eta[icel] + Dfib_eta[icel+(nz*ny)])) * T2u;

    Dfib_f[icel] = (T3_u[icel] - T2u) / dt;



    // ! ---------------  y  ---------------

    double T2v = T1_v[icel] - dt*( pressure[icel+nz]-pressure[icel] ) / Dys[j];

    T3_v[icel] = Dfib_eta[icel] * v_solid + (1.0-0.5*(Dfib_eta[icel]+Dfib_eta[icel+(1*nz)])) * T2v;

    Dfib_f[icel+iceltot]  = (T3_v[icel] - T2v) / dt;



    // ! ---------------  z  ---------------
    
    double T2w = T1_w[icel] - dt*( pressure[icel+1]-pressure[icel] ) / Dzs[k];

    T3_w[icel] = Dfib_eta[icel] * w_solid + (1.0-0.5*(Dfib_eta[icel]+Dfib_eta[icel+1])) * T2w;

    Dfib_f[icel*2*iceltot]  = (T3_w[icel]- T2w) / dt;
}

