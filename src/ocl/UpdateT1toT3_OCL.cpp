#include "4_5_1_UpdateT1toT3_OCL.hpp"

#ifdef OCL_ON


void UpdateT1toT3_OCL(
    OCLstruct& OCL, 
    DfibArray& Dfib,
    simpulationVariable& simu,
    pressure& t1,
    velocity& T1,
    velocity& T3,
    divideLocal& Lo,
    grid& gridA
){
    double u_solid = 0.0; 
    double v_solid = 0.0; 
    double w_solid = 0.0; 

    using cl_double = double;
    using cl_uint = uint;
    using cl_int =int32_t;

	static cl::KernelFunctor<const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,
                             cl::Buffer&, cl::Buffer&, cl::Buffer&, 
                             const cl::Buffer&, const cl::Buffer&, const cl::Buffer&,
                             const cl::Buffer&, cl::Buffer&, const cl::Buffer&,
                             const cl_double,const cl_double, const cl_double,
                             const cl_double, const cl_double, 
                             const cl_uint, const cl_uint, const cl_uint, const cl_uint,
                             const cl_uint
                             >
                             Kernel_UpdateT1toT3(OCL.prg_m["src/4_5_1_UpdateT1toT3.cl"], "update_UandF_seq");


	auto config3D = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D3[0]-1)/OCL.D3[0]*OCL.D3[0],
												(gridA.ny+OCL.D3[1]-1)/OCL.D3[1]*OCL.D3[1],
												(gridA.nz+OCL.D3[2]-1)/OCL.D3[2]*OCL.D3[2]},
												{OCL.D3[0], OCL.D3[1], OCL.D3[2]}); 


	Kernel_UpdateT1toT3(config3D,OCL.T1_u,OCL.T1_v, OCL.T1_w,
                                 OCL.T3_u,OCL.T3_v, OCL.T3_w,
                                 OCL.Dxs, OCL.Dys, OCL.Dzs,
                                 OCL.pressure, OCL.Dfib_f, OCL.Dfib_eta,
                                 u_solid, v_solid, w_solid,
                                 simu.dt, simu.nu,
                                 gridA.nx, gridA.ny, gridA.nz,gridA.gC,
                                 gridA.iceltot);

}

#endif