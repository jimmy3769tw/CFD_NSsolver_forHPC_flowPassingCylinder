
#include"4_0_4_BoundaryCondition_OCL.hpp"


#ifdef OCL_ON

void BoundaryCondtion_OCL(
    OCLstruct &OCL,
    simpulationVariable& simu,
    divideLocal& Lo,
    velocity& T0,
    pressure& t1,
    grid& gridA
)
{


    using cl_double = double;
    using cl_uint = uint;
    using cl_int = int32_t;


    struct Bcstruct
    {
      std::string u;
      std::string v;
      std::string w;
      std::string p;
    }static BC[6];

    if (simu.BoundaryConditionType == "Cavity-u-no_slip") {
    // x
        BC[0].u  = "BC_0_no_slip_Velocity_u";    BC[1].u  = "BC_1_no_slip_Velocity_u";
        BC[0].v  = "BC_0_no_slip_Velocity_v";    BC[1].v  = "BC_1_no_slip_Velocity_v";
        BC[0].w  = "BC_0_no_slip_Velocity_w";    BC[1].w  = "BC_1_no_slip_Velocity_w";
        BC[0].p  = "BC_0_Neumann_Pressure";      BC[1].p  = "BC_1_Neumann_Pressure"; 
    // y
        BC[2].u  = "BC_2_no_slip_Velocity_u";    BC[3].u  = "BC_3_Dirichlet_Velocity_u";
        BC[2].v  = "BC_2_no_slip_Velocity_v";    BC[3].v  = "BC_3_no_slip_Velocity_v";
        BC[2].w  = "BC_2_no_slip_Velocity_w";    BC[3].w  = "BC_3_no_slip_Velocity_w";
        BC[2].p  = "BC_2_Neumann_Pressure";      BC[3].p  = "BC_3_Neumann_Pressure"; 
    // z
        BC[4].u  = "BC_4_no_slip_Velocity_u";    BC[5].u  = "BC_4_no_slip_Velocity_u";
        BC[4].v  = "BC_4_no_slip_Velocity_v";    BC[5].v  = "BC_4_no_slip_Velocity_v";
        BC[4].w  = "BC_4_no_slip_Velocity_w";    BC[5].w  = "BC_4_no_slip_Velocity_w";
        BC[4].p  = "BC_4_Neumann_Pressure";      BC[5].p  = "BC_5_Neumann_Pressure";
    }

	const std::string sourcFileName = "src/4_0_4_BoundaryCondition.cl";

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_0_u(OCL.prg_m[sourcFileName],BC[0].u);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_0_v(OCL.prg_m[sourcFileName],BC[0].v);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_0_w(OCL.prg_m[sourcFileName],BC[0].w);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_1_u(OCL.prg_m[sourcFileName],BC[1].u);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_1_v(OCL.prg_m[sourcFileName],BC[1].v);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_1_w(OCL.prg_m[sourcFileName],BC[1].w);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_2_u(OCL.prg_m[sourcFileName],BC[2].u);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_2_v(OCL.prg_m[sourcFileName],BC[2].v);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_2_w(OCL.prg_m[sourcFileName],BC[2].w);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_3_u(OCL.prg_m[sourcFileName],BC[3].u);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_3_v(OCL.prg_m[sourcFileName],BC[3].v);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_3_w(OCL.prg_m[sourcFileName],BC[3].w);


	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_4_u(OCL.prg_m[sourcFileName],BC[4].u);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_4_v(OCL.prg_m[sourcFileName],BC[4].v);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_4_w(OCL.prg_m[sourcFileName],BC[4].w);



	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_5_u(OCL.prg_m[sourcFileName],BC[5].u);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_5_v(OCL.prg_m[sourcFileName],BC[5].v);

	static cl::KernelFunctor<cl::Buffer&,cl::Buffer&,cl::Buffer&
				,const cl_uint,const cl_uint ,const  cl_uint,const cl_uint >Kernel_BC_5_w(OCL.prg_m[sourcFileName],BC[5].w);





	auto config2D_x = cl::EnqueueArgs(OCL.queue,{ (gridA.ny+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.nz+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]}, 
												   OCL.D2[0],OCL.D2[1]); 

	auto config2D_y = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.nz+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]},
												  OCL.D2[0],OCL.D2[1]); 

	auto config2D_z = cl::EnqueueArgs(OCL.queue,{ (gridA.nx+OCL.D2[0]-1)/OCL.D2[0]*OCL.D2[0],
                                                  (gridA.ny+OCL.D2[1]-1)/OCL.D2[1]*OCL.D2[1]},
												   OCL.D2[0],OCL.D2[1]); 

        // FIXME:

    // * ----------------- BoundaryCondtion

	// FIXME: nx ny nz

	Kernel_BC_0_u(config2D_x, OCL.T0_u, OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_0_v(config2D_x, OCL.T0_u, OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_0_w(config2D_x, OCL.T0_u, OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);

	Kernel_BC_1_u(config2D_x, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_1_v(config2D_x, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_1_w(config2D_x, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);

	Kernel_BC_2_u(config2D_y, OCL.T0_u,OCL.T0_v, OCL.T0_u, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_2_v(config2D_y, OCL.T0_u,OCL.T0_v, OCL.T0_u, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_2_w(config2D_y, OCL.T0_u,OCL.T0_v, OCL.T0_u, gridA.nx, gridA.ny, gridA.nz, gridA.gC);

	Kernel_BC_3_u(config2D_y, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_3_v(config2D_y, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_3_w(config2D_y, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);

	Kernel_BC_4_u(config2D_z, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_4_v(config2D_z, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_4_w(config2D_z, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);

	Kernel_BC_5_u(config2D_z, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_5_v(config2D_z, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
	Kernel_BC_5_w(config2D_z, OCL.T0_u,OCL.T0_v, OCL.T0_w, gridA.nx, gridA.ny, gridA.nz, gridA.gC);
}

#endif