// #pragma OPENCL EXTENSION cl_khr_global_int32_base_atomic:enable
// #pragma OPENCL EXTENSION cl_khr_local_int32_base_atomic:enable

__kernel
void BC_0_no_slip_Velocity_u(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel0 = nz*ny + j*nz + k; //i=1 
    T0_u[icel0] = 0.0;
    T0_u[icel0-nz*ny] = T0_u[icel0];

}




__kernel
void BC_0_no_slip_Velocity_v(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel0 = nz*ny + j*nz + k;//i = 1 
    T0_v[icel0] = -T0_v[icel0+nz*ny];
    T0_v[icel0-nz*ny] = T0_v[icel0];

}



__kernel
void BC_0_no_slip_Velocity_w(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel0 = nz*ny + j*nz + k;
    T0_w[icel0] = -T0_w[icel0+nz*ny];
    T0_w[icel0-nz*ny] = T0_w[icel0];

}



__kernel
void BC_1_no_slip_Velocity_u(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel1 = (nx-3)*nz*ny + j*nz + k;
    T0_u[icel1] = 0.0;
    T0_u[icel1+nz*ny] = T0_u[icel1];
    T0_u[icel1+2*nz*ny] = T0_u[icel1+nz*ny];
}


__kernel
void BC_1_no_slip_Velocity_v(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel1 = (nx-3)*nz*ny + j*nz + k;
    T0_v[icel1+nz*ny] = -T0_v[icel1]; 
    T0_v[icel1+2*nz*ny] = T0_v[icel1+nz*ny];

}


__kernel
void BC_1_no_slip_Velocity_w(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel1 = (nx-3)*nz*ny + j*nz + k;
    T0_w[icel1+nz*ny] = -T0_w[icel1];
    T0_w[icel1+2*nz*ny] = T0_w[icel1+nz*ny];

}

__kernel
void BC_2_no_slip_Velocity_u(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel2 = i*nz*ny + nz + k;
    T0_u[icel2] = -T0_u[icel2+nz];
    T0_u[icel2-nz] = T0_u[icel2];

}


__kernel
void BC_2_no_slip_Velocity_v(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel2 = i*nz*ny + nz + k;
    T0_v[icel2] = 0.0;
    T0_v[icel2-nz] = T0_v[icel2];
}




__kernel
void BC_2_no_slip_Velocity_w(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel2 = i*nz*ny + nz + k;
    T0_w[icel2] = -T0_w[icel2+1*nz];
    T0_w[icel2-nz] = T0_w[icel2];
}



__kernel
void BC_3_Dirichlet_Velocity_u(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel3 = i*nz*ny + (ny-3)*nz + k;
    T0_u[icel3+1*nz] = 2.0*1.0-T0_u[icel3];
    T0_u[icel3+2*nz] = T0_u[icel3+1*nz];
}

__kernel
void BC_3_no_slip_Velocity_v(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel3 = i*nz*ny + (ny-3)*nz + k;
    T0_v[icel3] = 0.0;
    T0_v[icel3+  nz] =  T0_v[icel3];
    T0_v[icel3+2*nz] = T0_v[icel3+nz];
}


__kernel
void BC_3_no_slip_Velocity_w(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel3 = i*nz*ny + (ny-3)*nz + k;
    T0_w[icel3+  nz] = -T0_w[icel3];
    T0_w[icel3+2*nz] = T0_w[icel3+1*nz];

}


__kernel
void BC_4_no_slip_Velocity_u(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    const int icel4 = i*nz*ny + j*nz + 1;
    T0_u[icel4] = -T0_u[icel4+1];
    T0_u[icel4-1] = T0_u[icel4];
}

__kernel
void BC_4_no_slip_Velocity_v(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    const int icel4 = i*nz*ny + j*nz + 1;
    T0_v[icel4] = -T0_v[icel4+1];
    T0_v[icel4-1] = T0_v[icel4];            
}



__kernel
void BC_4_no_slip_Velocity_w(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    const int icel4 = i*nz*ny + j*nz + 1;
    T0_w[icel4] = 0.0;
    T0_w[icel4-1] = T0_w[icel4];
}



__kernel
void BC_5_no_slip_Velocity_u(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    const int icel5 = i*nz*ny + j*nz + (nz-3);
    T0_u[icel5+1] = -T0_u[icel5];
    T0_u[icel5+2] = T0_u[icel5+1];
}

__kernel
void BC_5_no_slip_Velocity_v(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
        // # Cavity-u-no_slip
    const int icel5 = i*nz*ny + j*nz + (nz-3);
    T0_v[icel5+1] = -T0_v[icel5];
    T0_v[icel5+2] = T0_v[icel5+1];
}
__kernel
void BC_5_no_slip_Velocity_w(
    __global double * T0_u,
    __global double * T0_v,
    __global double * T0_w,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    // # Cavity-u-no_slip
    const int icel5 = i*nz*ny + j*nz + (nz-3);
    T0_w[icel5] = 0.;
    T0_w[icel5+1] =  T0_w[icel5];
    T0_w[icel5+2] = T0_w[icel5+1];
}





__kernel
void BC_0_Neumann_Pressure(
    __global double * pressure,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel0 = nz*ny + j*nz + k;
    pressure[icel0] = pressure[icel0+nz*ny];
    pressure[icel0-nz*ny] = pressure[icel0];
}






__kernel
void BC_1_Neumann_Pressure(
    __global double * pressure,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint j = get_global_id(0);
    uint k = get_global_id(1);
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel1 = (nx-3)*nz*ny + j*nz + k;
    pressure[icel1+nz*ny] = pressure[icel1];
    pressure[icel1+2*nz*ny] = pressure[icel1+nz*ny];
}




__kernel
void BC_2_Neumann_Pressure(
    __global double * pressure,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel2 = i*nz*ny + nz + k;
    pressure[icel2] = pressure[icel2+nz];
    pressure[icel2-nz] = pressure[icel2];;
}



__kernel
void BC_3_Neumann_Pressure(
    __global double * pressure,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint k = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;
    const int icel3 = i*nz*ny + (ny-3)*nz + k;
    pressure[icel3+nz] = pressure[icel3];
    pressure[icel3+2*nz] = pressure[icel3+nz];
}






__kernel
void BC_4_Neumann_Pressure(
    __global double * pressure,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    const int icel4 = i*nz*ny + j*nz + 1;
    pressure[icel4] = pressure[icel4+1]; 
    pressure[icel4-1] = pressure[icel4];  
}


__kernel
void BC_5_Neumann_Pressure(
    __global double * pressure,
    uint nx, uint ny , uint nz, uint gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    const int icel5 = i*nz*ny + j*nz + (nz-3);
    pressure[icel5+1] = pressure[ icel5];
    pressure[icel5+2] = pressure[ icel5+1];

}