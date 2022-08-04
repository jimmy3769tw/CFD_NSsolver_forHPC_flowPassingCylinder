#include "1_1_Sequential.hpp"

#include"1_6_Ocl.hpp"

#ifdef OCL_DEBUG


bool RunOcl_Debug(
    grid &gridA,
    clockstruct &timer,
    simpulationVariable &simu
)
{
   timer.beginNew.start();
  //  * struct ---------------------------
  DfibArray Dfib;
  pressure t1;

  divide xDivide;
  divide yDivide;
  divide zDivide;
  
  velocity T0;
  velocity T1;
  velocity T3;

  sourceTerm So;
  
  SORcoefficient Sor;
  
  shareMenory ShareM;
  
  MxClass Mx;
  
  divideLocal Lo;

  OCLstruct OCL;

  //  * struct ---------------------------

  printf("pid: %d\n", getpid());

  // *         ================================  ================================

  xDivide.NumSubdomain = 1;
  yDivide.NumSubdomain = 1;
  zDivide.NumSubdomain = 1;
  resize_variable(gridA, simu, t1, So, T0, T1, T3, Dfib); // ! resize shared memory

  divideDomain(xDivide, gridA.nx - 4);
  divideDomain(yDivide, gridA.ny - 4);
  divideDomain(zDivide, gridA.nz - 4);

  // *         ============================  divid Domain ============================
  
  Lo.i_begin = xDivide.begin.at(0);
  
  Lo.j_begin = yDivide.begin.at(0);
  
  Lo.k_begin = zDivide.begin.at(0);
  
  Lo.i_endof = xDivide.endof.at(0);
  
  Lo.j_endof = yDivide.endof.at(0);
  
  Lo.k_endof = zDivide.endof.at(0);

  // *         ============================  divid Domain ============================

  // * -------------  Init the variables (First policy or data Locality) and Generate grids -------------
  
   T0.iniU(1.0, 1.0, 0.0);
   
   T1.iniU(0.0, 0.0, 0.0);
   
   T3.iniU(0.0, 0.0, 0.0);
   
   t1.init_p(0.0);

  BoundaryCondtion(simu, Lo, T0, t1, gridA); // ! dependence on generateGride

  generateGride(simu, ShareM, Lo, gridA);
  
  gridCheck_csv(simu, ShareM, Lo, gridA);

  OutputPlot3D_Xfile(simu, gridA);
  
  OutputPLOT3D_Qfile(Dfib, simu, t1, T3, Lo, gridA); ///Init
                                                        // ReadPlot3DadnCheckL2(ShareM,0,Dfib,simu,t1,T3,Lo,gridA); ///Init

  // * -------------  Init the variables (First policy or data Locality) and Generate grids -------------

  // * -------------------------- prepare for possion equation --------------------------

    timer.possion.start();


    Sor.cf.resize(gridA.iceltotCal * 8);

// CL_DEVICE_DOUBLE_FP_CONFIG

    std::fill(Dfib.eta.begin(), Dfib.eta.end(), 0.0L);

    size_t G = 1024;
    OCL.init("NVIDIA", CL_DEVICE_TYPE_GPU, 1*G); // Mesa, pocl

    std::vector<std::string> Source_File = {"src/4_0_4_BoundaryCondition.cl", //.0
                                            "src/4_4_1_UpdateT3toT0.cl",      //.1
                                            "src/4_2_1_SorPipeLine.cl",       //.2
                                            "src/4_1_1_CenterQuickscheme.cl", //.3
                                            "src/4_5_1_UpdateT1toT3.cl"       //.4
                                            };
    
    OCL.SetKernelProgram_from_SourceFile_m(Source_File);

    // filll it 
	// *create a buffer and fill it with data from an array C
	OCL.NEIBcell = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
		gridA.NEIBcell.size() * sizeof(cl_int), gridA.NEIBcell.data());

    // ? typedef int32_t cl_int
    // ? typedef double cl_double
    OCL.LocalMax = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE, 1 * sizeof(cl_double));
    OCL.iter = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE, 1 * sizeof(cl_uint));


    OCL.Dx = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
        gridA.Dx.size() * sizeof(cl_double), gridA.Dx.data());

    OCL.Dy= cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
        gridA.Dy.size() * sizeof(cl_double), gridA.Dy.data());

    OCL.Dz = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
        gridA.Dz.size()  * sizeof(cl_double), gridA.Dz.data());

    OCL.Dxs = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
        gridA.Dys.size() * sizeof(cl_double), gridA.Dxs.data());

    OCL.Dys= cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
        gridA.Dys.size() * sizeof(cl_double), gridA.Dys.data());

    OCL.Dzs = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
        gridA.Dzs.size()  * sizeof(cl_double), gridA.Dzs.data());

    // * ----------------- InitionCondition

    SorPopeLine_OCL_pre(OCL,ShareM, Sor, simu, T1, t1, Lo, gridA);

  // ! ==================  main loop begin ==================simu.loop_max
  for (size_t main_count = 1; main_count < 2; ++main_count)
  {

    simu.loop = main_count;

    // *1---------------get the T1(u strat, first velocity) --------------------

      CenterQuickScheme(So, simu, T0, T1, Lo, gridA);

    // *1--------------- get the T1(u strat, first velocity) --------------------


    OCL.T1_u = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                    T1.u.size() * sizeof(cl_double), T1.u.data());
    OCL.T1_v = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                    T1.v.size() * sizeof(cl_double), T1.v.data());
    OCL.T1_w = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                    T1.w.size() * sizeof(cl_double), T1.w.data());
    OCL.pressure = cl::Buffer(OCL.ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
                    t1.p.size() * sizeof(cl_double), t1.p.data());


    // std::vector<double> u(gridA.iceltot);
    // std::vector<double> v(gridA.iceltot);
    // std::vector<double> w(gridA.iceltot);
    // OCL.queue.enqueueReadBuffer(OCL.T1_u, CL_FALSE, 0, gridA.iceltot*sizeof(cl_double),u.data());
    // OCL.queue.enqueueReadBuffer(OCL.T1_v, CL_FALSE, 0, gridA.iceltot*sizeof(cl_double),v.data());
    // OCL.queue.enqueueReadBuffer(OCL.T1_w, CL_TRUE, 0, gridA.iceltot*sizeof(cl_double),w.data());




    SorPipeLine_OCL(OCL, ShareM, Sor, simu, T1, t1, Lo, gridA);
    // SorPipeLine_seq(ShareM, Sor, simu, T1, t1, Lo, gridA);
    SorPipeLine_seq_unpeeling(ShareM, Sor, simu, T1, t1, Lo, gridA);
    // * 2 ------------- get the pressure ------------
    std::vector<double> Check(Sor.cf.size());
    OCL.queue.enqueueReadBuffer(OCL.SorCoef, CL_TRUE, 0,Sor.cf.size() * sizeof(cl_double),Check.data());



    // for (size_t i = 0 ; i < Check.size() / 8 ; ++i){
    //     for (size_t j = 6 ; j < 7  ; j ++){
    //         // if (j != 6){
    //             if ( std::abs(Check.at(i*8+j) - Sor.cf.at(i*8+j)) >= 1.0e-4 ){
    //                 std::cout << j << ": "<< Check.at(i*8+j)  << ", " << Sor.cf.at(i*8+j) << endl;
    //             }
    //         // }
    //     }
    // }




    // OCL.queue.enqueueReadBuffer(OCL.pressure, CL_TRUE, 0,gridA.iceltot*sizeof(cl_double),t1.p.data());
    // OCL.queue.enqueueReadBuffer(OCL.T3_u, CL_FALSE, 0,gridA.iceltot*sizeof(cl_double),T3.u.data());
    // OCL.queue.enqueueReadBuffer(OCL.T3_v, CL_FALSE, 0,gridA.iceltot*sizeof(cl_double),T3.v.data());
    // OCL.queue.enqueueReadBuffer(OCL.T3_w, CL_TRUE, 0,gridA.iceltot*sizeof(cl_double),T3.w.data());
    // * 2 ------------- get the pressure ------------




    // * 3 ------------- update T1 to T3 (DFIB is included) -------------
                            timer.updateT1toT3.start();

    update_UandF_seq(Dfib, simu, t1, T1, T3, Lo, gridA); //! depend on 2.
                            
                            timer.updateT1toT3.stop();
    // * 3 ------------- update T1 to T3 (DFIB is included) -------------




    // *  ------------------ check staedy state (L2 norm) ------------------
                              timer.checkL2norm.start();

    CheckSteadyStatebyMaxVal_seq(simu, ShareM, T0, T3, Lo, gridA); //! indepent 3.
    
    CheckSteadyStatebyL2norm(simu, T0, T3, Lo, gridA);             //! indepent 3.


                              timer.checkL2norm.stop();
    // *  ------------------ check staedy state (L2 norm) ------------------




    // * .4. ------------- update T3 to T0 ------------- .4.
                    timer.updateT3toT0.start();

    UpdateT3toT0_seq(T3, T0); //! indepent 3.
                  

                      timer.updateT3toT0.stop();
    // * .4. ------------- update T3 to T0 ------------- .4.



    // * 5.------------- update the ghost cells  ------------- .5.
                          timer.BC.start();
    
    
    BoundaryCondtion(simu, Lo, T0, t1, gridA);
    
    
                          timer.BC.stop();
    // * 5. ------------- update the ghost cells  ------------- .5.

    // * 6.-------------.write plot3D formart .-------------
    if (simu.get_writefile())
    {
      OutputPLOT3D_Qfile(Dfib, simu, t1, T3, Lo, gridA);
      
      // ReadPlot3DadnCheckL2(ShareM,main_count,Dfib,simu,t1,T3,Lo,gridA);
      
      // getGhiaProfile(simu, T3, Lo, gridA);
    }
    // * 6.-------------.write plot3D formart .-------------

    timer.beginNew.stop();

    cout << "[Time]"
         << "simu :" << simu.dt * main_count << endl
         << "cal :" << timer.beginNew.elapsedTime() << endl
         << "===================================================" << endl
         << "main count :" << main_count + 1
         << ",[nx, ny, nz]=" << gridA.nx - 4 << ", " << gridA.ny - 4 << ", " << gridA.nz - 4
         << ", file :"<< simu.get_file() << endl;

    timer.beginNew.start();
    recorderTime(timer, simu, ShareM);
  }
  timer.beginNew.stop();
  // ! ==================  main loop begin ==================
  return true;
}

#endif
