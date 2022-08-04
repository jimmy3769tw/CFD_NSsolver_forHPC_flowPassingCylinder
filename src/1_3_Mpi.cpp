#include"run/1_3_Mpi.hpp"
#include "0_General.hpp"
#ifdef MPI_ON
#ifdef TEMP

#include"BoundaryCondition/Boundary_Condition.hpp"


bool Run_Hy_MPI_OpenMP(
    grid &gridA,
    clockstruct &timer,
    simpulationVariable &simu,
    int argc, char **argv
)
{

  using value_type = double;

  timer.beginNew.start();
  // *         ================================  ================================


  //  * struct ---------------------------
  DfibArray Dfib;
  pressure t1;
  divide xDivide;
  divide yDivide;
  divide zDivide;
  velocity T0;
  velocity T1;
  velocity T3;
  SORcoefficient Sor;
  shareMenory ShareM;
  MxClass Mx;
  
  //  * struct ---------------------------

  printf("\npid: %d\n", getpid());

  // *         ================================  ================================




  // *         ============================  divid Domain ============================


  resize_variable(gridA, simu, t1, T0, T1, T3, Dfib); // ! resize shared memory
  
  std::vector<int> dims{10}; //!

  bool reorder = true;

  auto [mx_comm_world, 
        mpi_word_size, 
        mpi_word_rank, 
        mpi_coord, 
        mpi_neighborhood]= mpi_init(reorder, dims, argc, argv);

    //   Check *= dims;
    // }    // for (auto dims:XYZ){
    //   Check *= dims;
    // }t(reorder, dims, argc, argv);

  for (size_t i = 0 ; i < (4 - dims.size()) ; ++i){
    dims.push_back(1);
    mpi_coord.push_back(0);
  }
  std::vector<int> grid_size{gridA.nx, gridA.ny, gridA.nz};

  // *         ============================  divid Domain ============================

  divideLocal Lo; // * local domain

  Lo.initTable(grid_size, dims);

  Lo.initLocal(mpi_coord.at(0), mpi_coord.at(1), mpi_coord.at(2));

  simu.PID = mpi_word_rank;

  for(size_t i = 0; i < mpi_word_size ; ++i){
    MPI_Barrier(mx_comm_world);
    if (mpi_word_rank == i){
      cout  << "[Rank: {begin,endof}] = "<< mpi_word_rank 
            << ": [{ "  << Lo.i_begin
            << ", "     << Lo.i_endof

            << "},{ "   << Lo.j_begin
            << ", "     << Lo.j_endof

            << "},{ "   << Lo.k_begin
            << ", "     << Lo.k_endof

            << "}]"     << endl;
    }
  }

  MPI_Barrier(mx_comm_world);
  
  if (simu.PID == 0) cout << "\n<i_begin_table>,<i_endof_table>";
  
  for(size_t j = 0; j < mpi_word_size ; ++j){
    cout << std::flush;
    MPI_Barrier(mx_comm_world);
    if (mpi_word_rank == j){
      cout << "\n===== " << j << " =====" << endl;

      for(size_t i = 0; i < mpi_word_size ; ++i){
      cout << "{" << Lo.i_begin_table.at(i) << ", ";
      cout << Lo.i_endof_table.at(i) << "}, ";
    }
    }
  }

  cout << std::flush;

  MPI_Barrier(mx_comm_world);

  if (simu.PID == 0) cout << "\n<i_length_table>";

  for(size_t j = 0; j < mpi_word_size ; ++j){
    MPI_Barrier(mx_comm_world);
    if (mpi_word_rank == j){
      cout << "\n===== " << j << " =====" << endl;

      for(size_t i = 0; i < mpi_word_size ; ++i){
        cout << Lo.i_length_table.at(i) << ", ";
      }
    }
  }




  // *         ============================  divid Domain ============================

  // * -------------  Init the variables (First policy or data Locality) and Generate grids -------------
  
  
  T0.iniU(1.0, 1.0, 0.0);
  
  T1.iniU(0.0, 0.0, 0.0);
  
  T3.iniU(0.0, 0.0, 0.0);
  
  t1.init_p(0.0);


  BoundaryCondtion(simu, Lo, T0, t1, gridA); // ! dependence on generateGride


  generateGride(simu, ShareM, Lo, gridA);


  if (simu.PID == 0){
    gridCheck_csv(simu, ShareM, Lo, gridA);
  }

  if (simu.PID == 0){
    OutputPlot3D_Xfile(simu, gridA);
  }
                  // ReadPlot3DadnCheckL2(ShareM,0,Dfib,simu,t1,T3,Lo,gridA); ///Init

  // * -------------  Init the variables (First policy or data Locality) and Generate grids -------------

  // * -------------------------- prepare for possion equation --------------------------
  /* prepare for possion equation */ timer.possion.start();

  #ifdef P_SOLVER_SOR

    Sor.cf.resize(gridA.iceltotCal * 8);

  #endif

  #ifdef P_SOLVER_BICG_ELL

    createPressureMatrix(Mx, simu, Lo, gridA);
    
    Mx.X_result.resize(gridA.iceltotCal, 0.0);

    Mx.matA.mpi_init(mx_comm_world,gridA.iceltotCal);

  #endif

  #ifdef EIGEN_ON

    // Mx.matA_Eigen.resize(gridA.iceltotCal, gridA.iceltotCal);The_simple_AmgclSolver_2.cpp
    
    // createPressureMatrix_Eigen(Mx, simu, Lo, gridA);
    
    // Mx.matA_Eigen.makeCompressed();
    
    // Mx.x_Eigen.resize(gridA.iceltotCal);

  #endif


  #ifdef P_SOLVER_AMGCL_BUILTIN
    
    // createPressureMatrix(Mx, simu, Lo, gridA);
    
    // auto [ptr, indices, values] = Mx.matA.get_CSR();
    
    // Mx.X_result.resize(gridA.iceltotCal, 0.0);
    
    // size_t rows = Mx.matA.row();

    // auto amgcl_mat = std::tie(rows, ptr, indices, values);
    
    // SolverBuiltin solver_Amgcl(amgcl_mat);

    // // debug
    // cout << "\n n:            " << rows;
    // cout << "\n ptr_size:     " << ptr.size();
    // cout << "\n ptr[n]:       " << ptr.at(rows);
    // cout << "\n indices(n+1)  " << indices.size();
    // cout << "\n values :      " << values.size() << endl;

  #endif


  /* prepare for possion equation */ timer.possion.stop();
  // * -------------------------- prepare for possion equation --------------------------

  // * ------------- DfibGetEta(Dfib,Lo,gridA) -------------
  if (simu.DfibMethod == "OFF")
  {
    std::fill(Dfib.eta.begin(), Dfib.eta.end(), 0.0);
  }
  else if (simu.DfibMethod == "DFIB_Cylinder-X")
  {

    std::fill(Dfib.eta.begin(), Dfib.eta.end(), 0.0);
    DFIB_CylinderX(Dfib, Lo, gridA);

  }
  else if (simu.DfibMethod == "DFIB_Cylinder-Z")
  {
    std::fill(Dfib.eta.begin(), Dfib.eta.end(), 0.0);
    DFIB_CylinderZ(Dfib, Lo, gridA);
  }
  else
  {
    throw std::invalid_argument("DFIB method??");
  }
    // std::fill(T0.u.begin(), T0.u.end(), simu.PID);
    // std::fill(T0.v.begin(), T0.v.end(), simu.PID);
    // std::fill(T0.w.begin(), T0.w.end(), simu.PID);

    // // ! ========================  MPI no-blocking send & recv ========================
    // #ifdef MPI_DEBUG

    // mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood, T0.u, Lo, gridA);

    // mpi_iSR_double_x(2, mx_comm_world, mpi_neighborhood, T0.v, Lo, gridA );

    // mpi_iSR_double_x_Collect_to_Master(mx_comm_world, simu.PID, mpi_word_size, T0.w, Lo, gridA );

    // #endif





  // * ------------- DfibGetEta(Dfib,Lo,gridA) -------------
  if (mpi_word_rank == 0){
    // 2021.Sep 2
    OutputPLOT3D_Qfile(Dfib, simu, t1, T0, Lo, gridA); ///Init
    // OutputDataFile(0, Dfib, simu, t1, T0, Lo, gridA);
  }

  // MPI_Finalize();

  // return 0;



  #ifdef TERBULENCE_SMAGORINSKY
    T1.Viseff.resize(gridA.iceltotCal);
  #endif

  auto EachLoopTime = timer.beginNew.elapsedTime();
  // ! ==================  main loop begin ==================
  for (size_t main_count = 1; main_count < simu.loop_max; ++main_count)
  {
    simu.loop = main_count;

    // *1---------------get the T1(u strat, first velocity) --------------------
    /* get the T1(Tsstart) */ timer.convectionDifussion.start();

    // ! ========================  MPI no-blocking send & recv ========================
    #ifdef MPI_DEBUG

    mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood, T0.u, Lo, gridA);

    mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood, T0.v, Lo, gridA);

    mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood, T0.w, Lo, gridA);

    #else

    mpi_iSR_double_x(2, mx_comm_world, mpi_neighborhood, T0.u, Lo, gridA );
    
    mpi_iSR_double_x(2, mx_comm_world, mpi_neighborhood, T0.v, Lo, gridA );

    mpi_iSR_double_x(2, mx_comm_world, mpi_neighborhood, T0.w, Lo, gridA );
    #endif

    // ! ========================  MPI no-blocking send & recv ========================

  
  #ifdef TERBULENCE_SMAGORINSKY
     
      SmagorinskyModel(ShareM, simu, T0, t1, Lo, gridA);

  #endif
  
      ConvectionDifussion(simu, T0, T1, Lo, gridA);


    /* get the T1(Tsstart) */ timer.convectionDifussion.stop();
    // * .1. ------------ get the T1(u strat, first velocity) ----------------- .1.




    // * .2. ------------- get the pressure ------------ .2.
     /* get the pressure */ timer.possion.start(); 

    // ! ========================  MPI no-blocking send & recv ========================
    #ifdef MPI_DEBUG
    mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood, T1.u, Lo, gridA);
    mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood, T1.v, Lo, gridA);
    mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood, T1.w, Lo, gridA);
    #else
    mpi_iSR_double_x(1, mx_comm_world, mpi_neighborhood, T1.u, Lo, gridA );
    mpi_iSR_double_x(1, mx_comm_world, mpi_neighborhood, T1.v, Lo, gridA );
    mpi_iSR_double_x(1, mx_comm_world, mpi_neighborhood, T1.w, Lo, gridA );
    // mpi_iSR_double_x_half(-1, mx_comm_world, mpi_neighborhood, T1.w, Lo, gridA );
    #endif

    // ! ========================  MPI no-blocking send & recv ========================

    #ifdef P_SOLVER_SOR
      // SorPipeLine_seq_unpeeling(ShareM,Sor,simu,T1,t1,Lo,gridA);
        SorPipeLine(mx_comm_world, mpi_neighborhood, ShareM, Sor, simu, T1, t1, Lo, gridA);
    #endif


    #ifdef P_SOLVER_EIGEN_CSR

      // createBMatrix_Eigen(T1, Mx, simu, Lo, gridA);
    
      // Mx.solver_Eigen.compute(Mx.matA_Eigen);
    
      // Mx.solver_Eigen.setTolerance(1e-15);
    
      // Mx.x_Eigen = Mx.solver_Eigen.solve(Mx.matB_Eigen);

      // Pressure_transform_x_Eigen(t1, Mx, Lo, gridA);

      // simu.iters = Mx.solver_Eigen.iterations();
      // simu.error = Mx.solver_Eigen.error();

      // std::cout << "#iters: " << Mx.solver_Eigen.iterations() << 
      //             ", error: " << Mx.solver_Eigen.error() << std::endl
    #endif

    #ifdef P_SOLVER_AMGCL_EIGEN
      
      // createBMatrix_Eigen(T1, Mx, simu, Lo, gridA);
      
      // createBMatrix_seq(T1, Mx, simu, Lo, gridA);
      
      // SolverEigen solve(Mx.matA_Eigen);
      
      // std::tie(simu.iters, simu.error) = solve(Mx.matB_Eigen, Mx.x_Eigen);
      
      // Pressure_transform_x_Eigen(t1, Mx, Lo, gridA);

      // simu.iters = iters;

      // std::cout << "#iters: "   << iters 
      //           << ", error: "  << error 
      //           << std::endl; 

      // std::cout << solve << std::endl;  // !amgcl information
    #endif


    #ifdef P_SOLVER_BICG_ELL

      createBMatrix_seq (T1,Mx,simu,Lo,gridA);
      
      std::tie(simu.iters, simu.error) = Mx.matA.npc_bicgstab_mpi(Mx.matB, Mx.X_result);

      Pressure_transform_X_result(t1, Mx, Lo, gridA);
      if (simu.PID == 0)
      cout << "[MPI::P_SOLVER_BICG_ELL] iter: " << simu.iters << ", error " << simu.error << endl;

    #endif

    #ifdef P_SOLVER_AMGCL_BUILTIN
      // createBMatrix_seq(T1, Mx, simu, Lo, gridA);

      // std::tie(simu.iters, simu.error) = solver_Amgcl(Mx.matB, Mx.X_result);

      // solver_Amgcl(amgcl_mat,Mx.matB, Mx.X_result);

      // Pressure_transform_X_result(t1, Mx, Lo, gridA);

      // std::cout << "[P_SOLVER_AMGCL_BUILTIN] #iters " << simu.iters << ", error " << simu.error << std::endl;
    #endif

    // ! ========================  MPI no-blocking send & recv ========================
    #ifdef MPI_DEBUG
      mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood, t1.p, Lo, gridA);
    #else
      mpi_iSR_double_x(1, mx_comm_world, mpi_neighborhood, t1.p, Lo, gridA );
      // mpi_iSR_double_x_half(1, mx_comm_world, mpi_neighborhood, t1.p, Lo, gridA);
    #endif
    // ! ========================  MPI no-blocking send & recv ========================

    /* get the pressure */ timer.possion.stop();
    // * .2. ------------- get the pressure ------------ .2.


    // * 3 ------------- update T1 to T3 (DFIB is included) -------------
    /* update T1 to T3 */  timer.updateT1toT3.start();

    update_UandF_seq(Dfib, simu, t1, T1, T3, Lo, gridA); //! depend on 2.
                            
    /* update T1 to T3 */  timer.updateT1toT3.stop();
    // * 3 -------------------------------------------------------------


    // *  ------------------ check staedy state (L2 norm) ------------------
    /* check staedy state */ timer.checkL2norm.start();
    
    CheckSteadyStatebyMaxVal_mpi(mx_comm_world, simu, ShareM, T0, T3, Lo, gridA); //! indepent 3.

    // CheckSteadyStatebyL2norm(simu, T0, T3, Lo, gridA);             //! indepent 3.

    /* check staedy state */ timer.checkL2norm.stop();
    // *  --------------------------------------------------------------------




    // * .4. ------------- update T3 to T0 ------------- .4.
    /* update T3 to T0 */ timer.updateT3toT0.start();

    T0.u = T3.u;
    T0.v = T3.v;
    T0.w = T3.w;

    /* update T3 to T0 */ timer.updateT3toT0.stop();
    // * .4. ------------------------------------------- .4.





    // * 5.------------- update the ghost cells  ------------- .5.
                          timer.BC.start();
    
    BoundaryCondtion(simu, Lo, T0, t1, gridA);
    
                          timer.BC.stop();
    // * 5. -------------------------------------------------- .5.




    // * 6.-------------.write plot3D formart .-------------
    if (simu.get_writefile())
    {
      // ! ========================  MPI no-blocking send & recv ========================
      // #ifdef MPI_DEBUG
        mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood,  T3.u, Lo, gridA);
        mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood,  T3.v, Lo, gridA);
        mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood,  T3.w, Lo, gridA);
        mpi_iSR_double_x_debugger(mpi_word_size, simu.PID, mx_comm_world, mpi_neighborhood,  t1.p, Lo, gridA);
      // #else
      //   mpi_iSR_double_x_Collect_to_Master(mx_comm_world, simu.PID, mpi_word_size, T3.u, Lo, gridA );
      //   mpi_iSR_double_x_Collect_to_Master(mx_comm_world, simu.PID, mpi_word_size, T3.v, Lo, gridA );
      //   mpi_iSR_double_x_Collect_to_Master(mx_comm_world, simu.PID, mpi_word_size, T3.w, Lo, gridA );
      //   mpi_iSR_double_x_Collect_to_Master(mx_comm_world, simu.PID, mpi_word_size, t1.p, Lo, gridA );
      // #endif

      // ! ========================  MPI no-blocking send & recv ========================

      MPI_Barrier(mx_comm_world);

      if (simu.PID == 1){
        OutputPLOT3D_Qfile(Dfib, simu, t1, T3, Lo, gridA);
        getCentProfile(simu, T3, Lo, gridA);
      }
      
      // ReadPlot3DadnCheckL2(ShareM,main_count,Dfib,simu,t1,T3,Lo,gridA);

    }
    // * 6.------------------------------------------------






    MPI_Barrier(mx_comm_world);           timer.beginNew.stop();
    if (mpi_word_rank == 0){
      cout 
          << "[Time]simu :" << simu.dt * main_count << endl
          << "[Time]cal < " << timer.beginNew.elapsedTime() - EachLoopTime 
          <<" OF "  << timer.beginNew.elapsedTime() << " >" << endl

          
          << "===================================================" << endl
          << "LOOP :" << main_count + 1
          << ",[nx, ny, nz]=" << gridA.nx - 4 << ", " << gridA.ny - 4 << ", " << gridA.nz - 4
          << ", file :" << simu.get_file()
      << endl;

      EachLoopTime = timer.beginNew.elapsedTime();
      timer.beginNew.start();
      recorderTime(timer, simu, ShareM);
    }
  }
  MPI_Barrier(mx_comm_world);           timer.beginNew.stop();
  // ! ==================  main loop begin ==================
  MPI_Finalize();  //Finish program !

  return true;
}



#endif


#endif
