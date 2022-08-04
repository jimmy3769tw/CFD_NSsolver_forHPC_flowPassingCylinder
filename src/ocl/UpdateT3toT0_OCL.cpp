#include "4_5_1_UpdateT1toT3_OCL.hpp"

#ifdef OCL_ON
void UpdateT3toT0_OCL(
    OCLstruct &OCL,
    simpulationVariable& simu,
    velocity& T3,
    velocity& T0,
    divideLocal& Lo,
    grid& gridA
){

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&, 
                            const cl::Buffer&,const cl::Buffer&,const cl::Buffer&,
                            const int,const int, const  int, const  int >
                            Kernel_UpdateT3toT0(OCL.prg_m["src/4_4_1_UpdateT3toT0.cl"], "UpdateT3toT0");

	auto config3D = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D3[0]-1)/OCL.D3[0]*OCL.D3[0],
												(gridA.ny+OCL.D3[1]-1)/OCL.D3[1]*OCL.D3[1],
												(gridA.nz+OCL.D3[2]-1)/OCL.D3[2]*OCL.D3[2]},
												{OCL.D3[0], OCL.D3[1], OCL.D3[2]}); 

	Kernel_UpdateT3toT0(config3D,OCL.T0_u, OCL.T0_v, OCL.T0_w,
                                 OCL.T3_u, OCL.T3_v, OCL.T3_w,
                                 gridA.nx, gridA.ny, gridA.nz, gridA.gC);

}
#endif