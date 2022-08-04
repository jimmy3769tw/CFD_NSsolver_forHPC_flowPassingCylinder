// #pragma OPENCL EXTENSION cl_khr_global_int32_base_atomic:enable
// #pragma OPENCL EXTENSION cl_khr_local_int32_base_atomic:enable


__kernel
void QuickUnitA(
    __local double * phiP,
    __local double * phiN,
    const double Up,
    const double Un,
    const int pp,const int p,const int icel,const int n,const int nn,
    __global double * Di,
    __global double * Ds,
    const uint Idx,
    __global double * phi
){
    double phi_P, phi_N;
        if (Up > 0.0)
            phi_P  = 0.5*(phi[icel] + phi[p] )
                    -0.125L*Di[Idx+1]*Di[Idx+1]/Ds[Idx]*
                    ( (phi[p] - phi[icel] ) / Di[Idx+1] 
                    - (phi[icel] - phi[n] ) / Di[Idx]) ;
        else
            phi_P  = 0.5*(phi[icel] + phi[p] )
                    - 0.125*Di[Idx+1]*Di[Idx+1]/Ds[Idx+1]*
                    ( (phi[pp] - phi[p] )   / Di[Idx+2] 
                    - (phi[p] - phi[icel] ) / Di[Idx+1] );

        if (Un > 0.0)
            phi_N  = 0.5*(phi[n] + phi[icel] ) 
                    -0.125L*Di[Idx]*Di[Idx]/Ds[Idx-1] *
                    ( (phi[icel] - phi[n] ) / Di[Idx] 
                    - (phi[n]   - phi[nn] ) / Di[Idx-1] ) ;
        else
            phi_N  = 0.5*(phi[n] + phi[icel] ) 
                    -0.125*Di[Idx]*Di[Idx]/Ds[Idx]*
                    ( (phi[p] - phi[icel] ) / Di[Idx+1] 
                    - (phi[icel] - phi[n] ) / Di[Idx])  ;
}



__kernel
void QuickUnitB(
    __local double * phiP,
    __local double * phiN,
    const double Up,
    const double Un,
    const int pp,const int p,const int icel,const int n,const int nn,
    __global double * Di,
    __global double * Ds,
    const uint Idx,
    __global double * phi
){
    double phi_P, phi_N;
        if (Up > 0.0)
            phi_P  = 0.5*(phi[icel] + phi[p] )
                    -0.125*Ds[Idx]*Ds[Idx]/Di[Idx]*
                    ( (phi[p] - phi[icel] ) / Ds[Idx] 
                    - (phi[icel] - phi[n] ) / Ds[Idx-1]) ;
        else
            phi_P  = 0.5*(phi[icel] + phi[p] )
                    - 0.125*Ds[Idx]*Ds[Idx]/Di[Idx+1]*
                    ( (phi[pp] - phi[p] )   / Ds[Idx+1] 
                    - (phi[p] - phi[icel] ) / Ds[Idx] );

        if (Un > 0.0)
            phi_N  = 0.5*(phi[n] + phi[icel] ) 
                    -0.125*Ds[Idx-1]*Ds[Idx-1]/Di[Idx-1] *
                    ( (phi[icel] - phi[n] ) / Ds[Idx-1] 
                    - (phi[n]   - phi[nn] ) / Ds[Idx-2] ) ;
        else
            phi_N  = 0.5*(phi[n] + phi[icel] ) 
                    -0.125*Ds[Idx-1]*Ds[Idx-1]/Di[Idx] *
                    ( (phi[p] - phi[icel] ) / Ds[Idx] 
                    - (phi[icel] - phi[n] ) / Ds[Idx-1])  ;
}





__kernel
void CenterQuickScheme_u(
    __global const double *T0_u,
    __global const double *T0_v,
    __global const double *T0_w,
    __global double *T1_u,
    __global double *T1_v,
    __global double *T1_w,
    __global const double *Dx,
    __global const double *Dy,
    __global const double *Dz,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global int *NEIBcell,
    const double dt,const double nu,
    const int nx,const int ny ,const int nz,const int gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;

        const int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
        int icel = i*nz*ny + j*nz + k;
        // * -------------- convection term --------------
        const int xn  = NEIBcell[icel*12+0]; /*-x negative */ const int xnn = NEIBcell[icel*12+6]; //-2x
        const int xp  = NEIBcell[icel*12+1]; /*+x positive */ const int xpp = NEIBcell[icel*12+7]; //+2x 
        const int yn  = NEIBcell[icel*12+2]; /*-y          */ const int ynn = NEIBcell[icel*12+8]; //-2y
        const int yp  = NEIBcell[icel*12+3]; /*+y          */ const int ypp = NEIBcell[icel*12+9]; //+2y
        const int zn  = NEIBcell[icel*12+4]; /*-z          */ const int znn = NEIBcell[icel*12+10];//-2z
        const int zp  = NEIBcell[icel*12+5]; /*+z          */ const int zpp = NEIBcell[icel*12+11];//+2z
        const double Upx = 0.5 * (T0_u[xp]+T0_u[icel]); 
        const double Unx = 0.5 * (T0_u[xn]+T0_u[icel]);
        const double Upy = 0.5 * (T0_v[xp]+T0_v[icel]);
        const double Uny = 0.5 * (T0_v[yn]+T0_v[yn+xp-icel]);
        const double Upz = 0.5 * (T0_w[xp]+T0_w[icel]);
        const double Unz = 0.5 * (T0_w[zn]+T0_w[zn+xp-icel]);

        __local double 
        uPx, uNx ,
        uPy, uNy ,
        uPz, uNz ;
        QuickUnitA (&uPx, &uNx, Upx, Unx, xpp, xp, icel, xn, xnn, Dx, Dxs, i, T0_u);
        QuickUnitB (&uPy, &uNy, Upy, Uny, ypp, yp, icel, yn, ynn, Dy, Dys, j, T0_u);
        QuickUnitB (&uPz, &uNz, Upz, Unz, zpp, zp, icel, zn, znn, Dz, Dzs, k, T0_u);

        double convection =   -dt*(
                                    (uPx*Upx-uNx*Unx) / Dxs[i]
                                    +(uPy*Upy-uNy*Uny) / Dy[j]
                                    +(uPz*Upz-uNz*Unz) / Dz[k]); ///FIXME: 2021.7.12
        // * -------------- convection term --------------

        // * -------------- difussion term --------------
        double difussion =    nu*dt*(
                                    ( (T0_u[xp]-T0_u[icel]) / Dx[i+1] 
                                    - (T0_u[icel]-T0_u[xn]) / Dx[i]   ) / Dxs[i]+
                                    ( (T0_u[yp]-T0_u[icel]) / Dys[j]   
                                    - (T0_u[icel]-T0_u[yn]) / Dys[j-1] ) / Dy[j]+
                                    ( (T0_u[zp]-T0_u[icel]) / Dzs[k]   
                                    - (T0_u[icel]-T0_u[zn]) / Dzs[k-1] ) / Dz[k]);
        // // * -------------- difussion term --------------
        double temp = difussion + convection;
        T1_u[icel] = T0_u[icel] +temp ;
}




__kernel
void CenterQuickScheme_v(
    __global const double *T0_u,
    __global const double *T0_v,
    __global const double *T0_w,
    __global double *T1_u,
    __global double *T1_v,
    __global double *T1_w,
    __global const double *Dx,
    __global const double *Dy,
    __global const double *Dz,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global int *NEIBcell,
    const double dt,const double nu,
    const int nx,const int ny ,const int nz,const int gCells
)
{
    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;

        const int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
        int icel = i*nz*ny + j*nz + k;
        // * -------------- convection term --------------
        const int xn  = NEIBcell[icel*12+0]; /*-x negative */ const int xnn = NEIBcell[icel*12+6]; //-2x
        const int xp  = NEIBcell[icel*12+1]; /*+x positive */ const int xpp = NEIBcell[icel*12+7]; //+2x 
        const int yn  = NEIBcell[icel*12+2]; /*-y          */ const int ynn = NEIBcell[icel*12+8]; //-2y
        const int yp  = NEIBcell[icel*12+3]; /*+y          */ const int ypp = NEIBcell[icel*12+9]; //+2y
        const int zn  = NEIBcell[icel*12+4]; /*-z          */ const int znn = NEIBcell[icel*12+10];//-2z
        const int zp  = NEIBcell[icel*12+5]; /*+z          */ const int zpp = NEIBcell[icel*12+11];//+2z
        const double Upx = 0.5 * (T0_u[icel]+T0_u[yp]);
        const double Unx = 0.5 * (T0_u[xn]+T0_u[xn+yn-icel]); 
        const double Upy = 0.5 * (T0_v[icel]+T0_v[yp]);
        const double Uny = 0.5 * (T0_v[icel]+T0_v[yn]);
        const double Upz = 0.5 * (T0_w[icel]+T0_w[yp]);
        const double Unz = 0.5 * (T0_w[zn]+T0_w[zn+yn-icel]);
        __local double  vPx, vNx,
                        vPy, vNy,
                        vPz, vNz;
        QuickUnitB (&vPx, &vNx, Upx, Unx, xpp, xp, icel, xn, xnn, Dx, Dxs, i, T0_v);
        QuickUnitA (&vPy, &vNy, Upy, Uny, ypp, yp, icel, yn, ynn, Dy, Dys, j, T0_v);
        QuickUnitB (&vPz, &vNz, Upz, Unz, zpp, zp, icel, zn, znn, Dz, Dzs, k, T0_v);


        const double convection =   -dt*(
                                        (Upx*vPx-Unx*vNx) / Dx[i] 
                                    +(Upy*vPy-Uny*vNy) / Dys[j] 
                                    +(Upz*vPz-Unz*vNz) / Dz[k]); 
    // * -------------- convection term --------------

    // * -------------- difussion term --------------
        const double difussion = nu*dt*( 
                                ( (T0_v[xp] - T0_v[icel]) / Dxs[i]
                                - (T0_v[icel] - T0_v[xn]) / Dxs[i-1] ) / Dx[i] +
                                ( (T0_v[yp] - T0_v[icel]) / Dy[j+1] 
                                - (T0_v[icel] - T0_v[yn]) / Dy[j] ) / Dys[j] +
                                ( (T0_v[zp] - T0_v[icel]) / Dzs[k] 
                                - (T0_v[icel] - T0_v[zn]) / Dzs[k-1])/ Dz[k] );
    // * -------------- difussion term --------------
        
    
        const double temp = convection + difussion;
                    
        T1_v[icel] = T0_v[icel] + temp;

}



__kernel
void CenterQuickScheme_w(
    __global const double *T0_u,
    __global const double *T0_v,
    __global const double *T0_w,
    __global double *T1_u,
    __global double *T1_v,
    __global double *T1_w,
    __global const double *Dx,
    __global const double *Dy,
    __global const double *Dz,
    __global const double *Dxs,
    __global const double *Dys,
    __global const double *Dzs,
    __global int *NEIBcell,
    const double dt,const double nu,
    const int nx,const int ny ,const int nz,const int gCells
)
{

    uint i = get_global_id(0);
    uint j = get_global_id(1);
    uint k = get_global_id(2);
    if ( i < gCells || i >= nx-gCells ) return;
    if ( j < gCells || j >= ny-gCells ) return;
    if ( k < gCells || k >= nz-gCells ) return;

        const int icelOffSet =  (i-2)*(nz-4)*(ny-4) + (j-2)*(nz-4) + (k-2);
        const int icel = i*nz*ny + j*nz + k;
        // * -------------- convection term --------------
        const int xp  = NEIBcell[icel*12+1]; /*+x positive */ const int xpp = NEIBcell[icel*12+7]; //+2x 
        const int xn  = NEIBcell[icel*12+0]; /*-x negative */ const int xnn = NEIBcell[icel*12+6]; //-2x
        const int yn  = NEIBcell[icel*12+2]; /*-y          */ const int ynn = NEIBcell[icel*12+8]; //-2y
        const int yp  = NEIBcell[icel*12+3]; /*+y          */ const int ypp = NEIBcell[icel*12+9]; //+2y
        const int zn  = NEIBcell[icel*12+4]; /*-z          */ const int znn = NEIBcell[icel*12+10];//-2z
        const int zp  = NEIBcell[icel*12+5]; /*+z          */ const int zpp = NEIBcell[icel*12+11];//+2z
        const double Upx = 0.5 * ( T0_u[zp]+T0_u[icel] );
        const double Unx = 0.5 * ( T0_u[xn+1]+T0_u[xn] ); 
        const double Upy = 0.5 * ( T0_v[icel]+T0_v[zp] );
        const double Uny = 0.5 * ( T0_v[yn+1]+T0_v[yn] );
        const double Upz = 0.5 * ( T0_w[zp]+T0_w[icel] );
        const double Unz = 0.5 * ( T0_w[icel]+T0_w[zn] );
        __local double  wPx, wNx,
                        wPy, wNy,
                        wPz, wNz;
        QuickUnitB (&wPx, &wNx, Upx, Unx, xpp, xp, icel, xn, xnn, Dx, Dxs, i, T0_w);
        QuickUnitB (&wPy, &wNy, Upy, Uny, ypp, yp, icel, yn, ynn, Dy, Dys, j, T0_w);
        QuickUnitA (&wPz, &wNz, Upz, Unz, zpp, zp, icel, zn, znn, Dz, Dzs, k, T0_w);

        const double convection =   -dt*(
                                     (Upx*wPx-Unx*wNx) / Dx[i]
                                    +(Upy*wPy-Uny*wNy) / Dy[j]
                                    +(Upz*wPz-Unz*wNz) / Dzs[k]);
        // * -------------- convection term --------------

        // * -------------- difussion term --------------
        const double difussion = nu*dt*(
                                    ( ( T0_w[xp]-T0_w[icel] ) / Dxs[i] 
                                    - ( T0_w[icel]-T0_w[xn] ) / Dxs[i-1] ) / Dx[i]+
                                    ( ( T0_w[yp]-T0_w[icel] ) / Dys[j] 
                                    - ( T0_w[icel]-T0_w[yn] ) / Dys[j-1] ) / Dy[j]+
                                    ( ( T0_w[zp]-T0_w[icel] ) / Dz[k+1] 
                                    - ( T0_w[icel]-T0_w[zn] ) / Dz[k]   ) / Dzs[k]) ;
        // * -------------- difussion term --------------
        const double temp = convection + difussion;
                     
    T1_w[icel] = T0_w[icel] + temp;
}
