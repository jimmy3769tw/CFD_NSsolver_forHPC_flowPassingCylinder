#include "0_General.hpp"

#include "run/1_1_Sequential.hpp"
#include "run/1_2_Omp.hpp"

#include "run/general.hpp"


bool runSeqOmp(
    grid &gA,
    clockstruct &timer,
    simpulationVariable &simu,
    int argc, char **argv
){

  timer.beginNew.start();
  //  * struct ---------------------------
  DfibArray Dfib;

  pressure t1;

  velocity T0 , T1, T3;

  SORcoefficient Sor;

  shareMenory ShareM;

  MxClass Mx;

  //  * struct ---------------------------

  printf("pid: %d\n", getpid());

#if defined (PC_OMP)
  // ~~~~~~~~~~~~~~~~~~~~~~~~ openMP ~~~~~~~~~~~~~~~~~~~~~~~~
  #ifdef _OPENMP
    std::cout << "OpenMP: " << _OPENMP << std::endl;
  #else
    std::cout << "Your compiler does not support OpenMP." << std::endl;
  #endif
  printf("pid: %d\n", getpid());

  int ompThreads = omp_get_max_threads() ;  // Using Inv OMP_NUM_THREADS
  omp_set_num_threads(ompThreads);
  // ~~~~~~~~~~~~~~~~~~~~~~~~ openMP ~~~~~~~~~~~~~~~~~~~~~~~~
#endif

  // ! ============================  divid Domain ============================

  divideLocal Lo;

  std::vector<int> grid_size{gA.nx, gA.ny, gA.nz};

  std::vector<int> dims{1, 1, 1};

  Lo.initTable(grid_size, dims);

  Lo.initLocal(0, 0, 0);

  // * ========================================================================

  // ! ## Init the variables (First policy or data Locality) and Generate grids 

  resize_variable(gA, simu, t1, T0, T1, T3, Dfib); // ! resize shared memory

  double ui = 1.0, vi = 0.0, wi = 0.0;

  #if defined (PC_SEQ)
    T0.iniU(ui, vi, wi);
    T1.iniU(ui, vi, wi);
    T3.iniU(ui, vi, wi);
  #elif defined (PC_OMP)
    T0.iniU_omp(ui, vi, wi);
    T1.iniU_omp(ui, vi, wi);
    T3.iniU_omp(ui, vi, wi);
  #endif

  t1.init_p(0.0);

  generateGride(simu, ShareM, Lo, gA);

  gA.io_csv();

  OutputPlot3D_Xfile(simu, gA);


  // ! ## Read the exist data .-------------
  // auto [Nblock_read, mach_read, alpha_read, reyn_read, time_read, 
  //         pPre, ucPre, vcPre, wcPre, EtaPre] = QfileRead(gA, "mx_in/in.q");

  // t1.p = pPre;
  // T1.u = T3.u = T0.u = ucPre;
  // T1.v = T3.v = T0.v = vcPre;
  // T1.w = T3.w = T0.w = wcPre;
  // PotentialFlow(simu, Dfib, T0, t1, Lo, gA);

  // simu.set_Re(reyn_read);

  // simu.set_time(time_read);

  // * --------------------------------------

  // ! ## The first BC -------------
  T1.u = T3.u = T0.u ;
  T1.v = T3.v = T0.v ;
  T1.w = T3.w = T0.w ;

  BC_staggered_main( Lo, T0, t1, gA );
  BC_staggered_copy( Lo, T0, T1, gA );

  // !## prepare for poisson equation --------------------------
    timer.p.start();


  #if  defined (P_SOLVER_ELL) \
    || defined (P_SOLVER_CSR) \
    || defined (P_SOLVER_SPE) \
    || defined (P_SOLVER_AMGCL_BUILTIN) \

    Mx.X_result.resize(gA.iceltotCal, 0.0);

  #endif




  // !### ---------------------------------------- ELL (BICG)
  #if defined (P_SOLVER_ELL)
  // ---------------------------
    timer.set_MaA.start();
    createPressureMatrix(Mx, simu, Lo, gA);
    timer.set_MaA.stop();
    std::cout << "ELL_setUptime : "<< timer.set_MaA.elapsedTime() << std::endl;
     Mx.matA.analyse();
  // ---------------------------

  // ---------------------------
    solver::bicgstab<mat::ELL_matrix<double> > pSolver(Mx.matA);
  // ---------------------------

  // ---------------------------
  // !### ---------------------------------------- CSR
  #elif defined (P_SOLVER_CSR)

    //  ----------------------------------------------------------------
    Mx.matA_csr.set( gA.createPressureMatrixCSR() );
    //  ----------------------------------------------------------------

    //  ----------------------------------------------------------------
    solver::bicgstabRe2<mat::CSR_matrix<double> > 
      pSolver(Mx.matA_csr);
    //  ----------------------------------------------------------------

  // !### ---------------------------------------- SPE
  #elif defined(P_SOLVER_SPE)

  // --------------------------------------
  Mx.matA_spe.resize(gA.nxCal, gA.nyCal, gA.nzCal);
  // --------------------------------------

  Mx.matA_spe.setupPressure( gA.gC, gA.Dx, gA.Dy, gA.Dz,gA.Dxs, gA.Dys, gA.Dzs);
    solver::bicgstabRe2<mat::SPE_matrix<double> > 
      pSolver(Mx.matA_spe);
  
  // !### ---------------------------------------- EIGEN
  #elif defined (EIGEN_ON)


    // -------------------------------------------
    timer.set_MaA.start();
    Mx.x_Eigen.resize(gA.iceltotCal);
    Mx.matA_Eigen.resize(gA.iceltotCal, gA.iceltotCal);
    // -------------------------------------------

    createPressureMatrix_Eigen(Mx, simu, Lo, gA);
    
    Mx.matA_Eigen.makeCompressed();

    timer.set_MaA.stop();

    std::cout << "[EIGEN] setUp time : "<< timer.set_MaA.elapsedTime() << std::endl;
    // -------------------------------------------

  // !### ---------------------------------------- AMGCL(buildin)
  #elif defined (P_SOLVER_AMGCL_BUILTIN)

    // -------------------------------------------
    auto [ptr, idx, values] = gA.createPressureMatrixCSR();

    auto amgcl_mat = std::tie(gA.iceltotCal, ptr, idx, values);
    // -------------------------------------------
    SolverBuiltin::params prm;

    prm.solver.tol = simu.p_criteria;

    prm.solver.maxiter = simu.p_iterMax;

    prm.solver.ns_search = true;

    prm.solver.verbose = true;
    SolverBuiltin solver_Amgcl(amgcl_mat, prm);
    // ---------------------------------
    // SolverBuiltin solver_Amgcl(amgcl_mat);
    // -------------------------------------------
    // amgcl::io::mm_write("test_io_vec.mm", rhs.data(), n, 1);

    amgcl::io::mm_write("Information/This_Matrix.mtx", amgcl_mat);
    // -------------------------------------------

  #endif




  #if  defined (P_SOLVER_ELL) \
    || defined (P_SOLVER_CSR) \
    || defined (P_SOLVER_SPE)

    pSolver.setTolerance(simu.p_criteria);

  #endif


  timer.p.stop();
  // * ----------------------------------------------------

  // ! ##  DfibGetEta(Dfib, Lo, gA) -------------
  std::fill(Dfib.eta.begin(), Dfib.eta.end(), 0.0);

  if (simu.DfibMethod == "OFF"){}
  else if (simu.DfibMethod == "DFIB_Cylinder-X")
    DFIB_CylinderX(Dfib, Lo, gA);
  else if (simu.DfibMethod == "DFIB_Cylinder-Z")
    DFIB_CylinderZ(Dfib, Lo, gA);
  else
    throw std::invalid_argument("DFIB method??");

  // * ------------------------------------------------------

    // ---------------------------------
    OutputPLOT3D_Qfile(Dfib, simu, t1, T0, Lo, gA); ///Init
    // ---------------------------------

    // ---------------------------------
    #ifdef TERBULENCE_SMAGORINSKY
    T0.Viseff.resize(gA.iceltotCal, 0.0);
    #endif
    // ---------------------------------

  // !  ------------- remove velocity (eta == 0) -------------
  // removeVelocity(Dfib, T0, Lo, gA);
  // removeVelocity(Dfib, T1, Lo, gA);
  // removeVelocity(Dfib, T3, Lo, gA);
    // * -------------------------------------------------------------
  auto perLoopWTime = timer.beginNew.elapsedTime();
  // ! ==================  main loop begin ==================
  for (simu.loop = 1; simu.get_finishloop(); simu.finishloop())
  {

    // ---------------------------------
    #if defined (TERBULENCE_SMAGORINSKY)
        SmagorinskyModel(ShareM, simu, T0, t1, Lo, gA);
    #endif
    // ---------------------------------

    // !## 1. get the T1/u* --------------------
    /*  T1(Tsstart) */ timer.convectionDifussion.start();
      ConvectionDifussion(simu, T0, T1, Lo, gA);
      // BC_updateSlid(Lo, T1, gA);
    /*  T1(Tsstart) */ timer.convectionDifussion.stop();
    //*---------------------------------------------------

    // -----------------------
    get_Cfl(T1, Lo, gA, simu);
    // simu.timestepper();
    // -----------------------



    // !## 2. get the pressure ------------ .2.
     /* get the pressure */ timer.p.start(); 


    // !### ---------------------------------------- SOR
    #if defined (P_SOLVER_SOR)

    if (simu.loop > 1)
    {
    // ---------------------------------
    #if defined (PC_OMP)
    SorPipeLine_omp(ShareM, Sor, simu, T1, t1, Lo, gA);
    #elif defined (PC_SEQ)
    SorPipeLine_seq(ShareM, Sor, simu, T1, t1, Lo, gA);
    #endif
    // ---------------------------------
    }

    // !### ---------------------------------------- CSR/SPE/ELL (BICG)(CG)
    #elif defined (P_SOLVER_BICG_CSR) \
      ||  defined (P_SOLVER_BICG_SPE) \
      ||  defined (P_SOLVER_BICG_ELL)

    if (simu.loop > 1)
    {

    // ---------------------------------
    #if defined (PC_SEQ)
    createBMatrix_seq (T1, Mx, simu, Lo, gA);
    #elif defined (PC_OMP)
    createBMatrix_omp(T1, Mx, simu, Lo, gA);
    #endif
    // ---------------------------------


    // ---------------------------------
      std::tie(simu.iters, simu.error) = pSolver.solve(Mx.matB, Mx.X_result);
    // ---------------------------------

    // ---------------------------------
      Pressure_transform_X_result(t1, Mx, Lo, gA);
    // ---------------------------------
    }


    // !### ---------------------------------------- AMGCL
    #elif defined (P_SOLVER_AMGCL_BUILTIN)
    // if (simu.loop > 1000)
    {

    // ----------------------
    #if defined (PC_SEQ)
    createBMatrix_seq(T1, Mx, simu, Lo, gA);
    #elif defined (PC_OMP)
    createBMatrix_omp(T1, Mx, simu, Lo, gA);
    #endif
    // ----------------------

    std::tie(simu.iters, simu.error) = solver_Amgcl(Mx.matB, Mx.X_result);
    Pressure_transform_X_result(t1, Mx, Lo, gA);

    }
    #endif 
    // !### ---------------------------------------- endif
    

    /* get the pressure */ timer.p.stop(); simu.printInfo();


    double min_p = 10e+9;
    auto max_p = -min_p;
    for (size_t i = Lo.i_begin; i < Lo.i_endof; ++i)
    for (size_t j = Lo.j_begin; j < Lo.j_endof; ++j)
    for (size_t k = Lo.k_begin; k < Lo.k_endof; ++k)
    {
      min_p = std::min(t1.p[gA.icel(i,j,k)], min_p);
      max_p = std::max(t1.p[gA.icel(i,j,k)], max_p);
    }

    std::cout  << "[Min_p, Max_p]:"<< min_p  << ", "<< max_p << std::endl;

    // * .2. ------------- get the pressure ------------ .2.

    // !## 3. update T1 to T3 (DFIB is included) -------------
                            timer.updateT1toT3.start();
    // ----------------------
    #if defined (PC_SEQ)
    update_UandF_seq(Dfib, simu, t1, T1, T3, Lo, gA);
    #elif defined (PC_OMP)
    update_UandF_omp(Dfib, simu, t1, T1, T3, Lo, gA); 
    #endif
    // ----------------------
    // BC_updateSlid(Lo, T3, gA); //Ahmad didn't up it !!
                            timer.updateT1toT3.stop();
    // * 3 -------------------------------------------------------------



    // !  ------------------ check staedy state (L2 norm) ------------------
                             timer.checkL2norm.start();
    #if defined (PC_SEQ)
    CheckSteadyStatebyMaxVal_seq(simu, ShareM, T0, T3, Lo, gA); 
    #elif defined (PC_OMP)
    CheckSteadyStatebyMaxVal_omp(simu, ShareM, T0, T3, Lo, gA);
    #endif
                             timer.checkL2norm.stop();
    // *  -------------------------------------------------------------------


    // ! ## 4. cpoy T3 to T0 ------------------
              timer.updateT3toT0.start();
    copyVel(T3, T0);
              timer.updateT3toT0.stop();
    // * ----------------------------------------


    // !## 5  Update the ghost cells  -------------
    timer.BC.start();
    
    BC_staggered_main( Lo, T0, t1, gA );
    BC_staggered_copy( Lo, T0, T1, gA );

    timer.BC.stop();
    // * -----------------------------------------

    // !## 6 get Cd and Cl  -------------
    if (simu.DfibMethod != "OFF")
    { CD_CL(simu, Lo, Dfib, gA);}
    // * --------------------------------

    // !## 7 write plot3D formart .-------------
    if (simu.get_writefile()){
      OutputPLOT3D_Qfile(Dfib, simu, t1, T3, Lo, gA);
    }
    // * -----------------------------------------


    // !## 8. IO loop------------------------------------------------
    timer.beginNew.stop();
    cout 
    <<     "[simu]{time, dt} ="  << simu.get_time()          
    <<     ", " << simu.dt       << endl
    <<     "[wall time] < "      << timer.beginNew.elapsedTime() - perLoopWTime 
    <<     " OF "                << timer.beginNew.elapsedTime()
    <<     " >" << endl

    // ! =============== NEXT time step ===============
    <<     "===================================================" << endl
    <<     "LOOP :"  << simu.loop + 1
    <<     gA.show() << ", file :" << simu.get_file() << endl;
    // * ----------------------------------------------------------


    timer.beginNew.start(); perLoopWTime = timer.beginNew.elapsedTime();

    recorderTime(timer, simu, ShareM);
  }
  timer.beginNew.stop();
  // ! ==================  main loop begin ==================
  return true;
}
