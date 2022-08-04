#include"run/1_2_Omp.hpp"
// BUG code
#ifdef OMP_ON

#ifdef TEMP

#include"BoundaryCondition/Boundary_Condition.hpp"

bool RunOmpSubDomain3D(
  grid& gridA,
  clockstruct& timer,
  simpulationVariable& simu
){
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
    SORcoefficient Sor;    
    MxClass Bicg;
    shareMenory ShareM;
       
    //  * struct ---------------------------


    #ifdef _OPENMP
      std::cout << "OpenMP: " << _OPENMP << std::endl;
    #else
      std::cout << "Your compiler does not support OpenMP." << std::endl;
    #endif
    printf("pid: %d\n", getpid());

    // *openMP
    int ompThreads = omp_get_max_threads() ;  // Using Inv OMP_NUM_THREADS
    omp_set_num_threads(ompThreads);

    // export OMP_SCHEDULE='STATIC'
    // export OMP_NUM_THREADS 
    // export OMP_PLACES=cores #pragma omp parallel  schedule(static,chunksize)
    // export OMP_PROC_BIND=spread #pragma omp parallel proc_bind(spread)
    // export OMP_PROC_BIND=close #pragma omp parallel proc_bind(close)

    std::cout << "ompThreads :" << ompThreads << std::endl;


    // *         ================================ SubDomain3D ================================

    xDivide.NumSubdomain    = ompThreads;
    yDivide.NumSubdomain    = 1;
    zDivide.NumSubdomain    = 1;
    resize_variable( gridA, simu, t1, T0, T1, T3, Dfib);// ! resize shared memory
    divideDomain(xDivide,gridA.nx-4);    divideDomain(yDivide,gridA.ny-4);    divideDomain(zDivide,gridA.nz-4); 
    // *         ================================ SubDomain3D ================================



    double pChangeMax =0.;

    // ! private ->  timer, simu, Dfib, Sor
    // #pragma omp parallel num_threads(ompThreads) default(none)\
    //         shared(Bicg, ShareM,pChangeMax,cout,xDivide,yDivide,zDivide,t1,T0,T1,T3, gridA)\
    //         firstprivate(timer,simu,ompThreads,Dfib,Sor)//(nonvariant && doing again and again)
    {
      divideLocal Lo;

      // *         ============================  divid Domain ============================
      simu.TID   = omp_get_thread_num();
      Lo.i_begin = xDivide.begin.at(simu.TID);        Lo.j_begin = yDivide.begin.at(0);        Lo.k_begin = zDivide.begin.at(0);
      Lo.i_endof = xDivide.endof.at(simu.TID);        Lo.j_endof = yDivide.endof.at(0);        Lo.k_endof = zDivide.endof.at(0);

      #pragma omp parallel
      {
        omp_display_affinity(NULL);
      }
      // *         ============================  divid Domain ============================


      // * -------------  Init the variables (First policy or data Locality) and Generate grids -------------
      Locality(simu,Lo,T0,T1,T3,t1,gridA); // data Locality

      BoundaryCondtion(simu, Lo, T0, t1, gridA); // ! dependence on generateGride
      
      generateGride(simu,ShareM,Lo,gridA);
  

      #pragma omp single nowait
      {
        gridCheck_csv(simu,ShareM,Lo,gridA);
      }

      #pragma omp single
      {
        OutputPlot3D_Xfile(simu,gridA);
      }

      #pragma omp barrier
      #pragma omp flush
  
      #pragma omp single
      {
        OutputPLOT3D_Qfile(Dfib, simu, t1, T3, Lo, gridA);
      }  ///Init 
      // ReadPlot3DadnCheckL2(ShareM,0,Dfib,simu,t1,T3,Lo,gridA); ///Init 
      // * -------------  Init the variables (First policy or data Locality) and Generate grids -------------

      // * -------------------------- prepare for possion equation --------------------------
                                      timer.possion.start();
      #ifdef P_SOLVER_SOR
        Sor.cf.resize(gridA.iceltotCal*8);
      #endif

                                      timer.possion.stop();
      // * -------------------------- prepare for possion equation --------------------------


      // * ------------- DfibGetEta(Dfib,Lo,gridA) -------------
      if (simu.DfibMethod == "OFF"){
        DfibInit(Dfib,Lo,gridA); //data Locality
      }
      else if (simu.DfibMethod == "DFIB_Cylinder-X"){
        DfibInit(Dfib,Lo,gridA); //data Locality
        DFIB_CylinderX(Dfib, Lo, gridA);
      }
      else if (simu.DfibMethod == "DFIB_Cylinder-Z"){
        DfibInit(Dfib,Lo,gridA); //data Locality
        DFIB_CylinderZ(Dfib, Lo, gridA);
      }
      else{
        throw std::invalid_argument("DFIB method??");
      }
      // * ------------- DfibGetEta(Dfib,Lo,gridA) -------------


      // ! ==================  main loop begin ==================
      for (size_t main_count = 1; main_count < simu.loop_max ; ++main_count){ 

        double elapTime = omp_get_wtime();
        simu.loop = main_count;
      
        // *1---------------get the T1(u strat, first velocity) --------------------
   /* get the T1(Tsstart) */ timer.convectionDifussion.start();


    #ifdef TERBULENCE_SMAGORINSKY
     
      SmagorinskyModel(ShareM, simu, T0, t1, Lo, gridA);

    #endif
  
      ConvectionDifussion(simu, T0, T1, Lo, gridA);


    /* get the T1(Tsstart) */ timer.convectionDifussion.stop();
        // *1--------------- get the T1(u strat, first velocity) --------------------





        // * 2 ------------- get the pressure ------------
                          timer.possion.start();
        #ifdef P_SOLVER_SOR
          SorPipeLineSubDomain3D_omp(ShareM,Sor,simu,T1,t1,Lo,gridA);
        #endif
                          timer.possion.stop();
        // * 2 ------------- get the pressure ------------



        // * 3 ------------- update T1 to T3 (DFIB is included) -------------
                              timer.updateT1toT3.start();
        update_UandF_seq(Dfib,simu,t1,T1,T3,Lo,gridA); //! depend on 2.
                              timer.updateT1toT3.stop();

        // * 3 ------------- update T1 to T3 (DFIB is included) -------------


        // *  ------------------ check staedy state (L2 norm) ------------------

                                timer.checkL2norm.start();
        CheckSteadyStatebyMaxVal_omp(simu, ShareM, T0, T3, Lo, gridA); //! indepent 3.
        
        CheckSteadyStatebyL2norm(simu, T0, T3, Lo, gridA);         //! indepent 3.
                                timer.checkL2norm.stop();
        // *  ------------------ check staedy state (L2 norm) ------------------

        // * 4 ------------- update T3 to T0 -------------
                          #pragma omp barrier
                      timer.updateT3toT0.start();
        T0.u = T3.u;
        T0.v = T3.v;
        T0.w = T3.w;
                      timer.updateT3toT0.stop();
        // * 4 ------------- update T3 to T0 -------------

        // * 5.------------- update the ghost cells  -------------
                              timer.BC.start();
        BoundaryCondtion(simu, Lo, T0, t1, gridA);
                              timer.BC.stop();
        // * 5.------------- update the ghost cells  -------------



        // * 6.-------------.write plot3D formart .-------------
        if ( simu.get_writefile() ){

          #pragma omp barrier
          #pragma omp flush
          #pragma omp single
          {
            OutputPLOT3D_Qfile(Dfib, simu, t1, T3, Lo, gridA);
          }
          // Output_Q_Cpp_PLOT3D(simu,T3, t1, Dfib, gridA);
          // ReadPlot3DadnCheckL2(ShareM,main_count,Dfib,simu,t1,T3,Lo,gridA);
          // getGhiaProfile(simu, T3, Lo, gridA);
        }
        // * 6.-------------.write plot3D formart .-------------
 
        #pragma omp barrier
        timer.beginNew.stop();

        #pragma omp single
        {
          cout  << "[Time]"<< "simu :" << simu.dt*main_count << endl
                << "cal :"  << timer.beginNew.elapsedTime() << endl
                << "===================================================" << endl
                << "main count :" << main_count+1
                << ",[nx, ny, nz]=" << gridA.nx - 4 << ", "<< gridA.ny - 4 << ", "<< gridA.nz - 4
          << endl;
        }
        timer.beginNew.start();

        #pragma omp critical
        {
          recorderTime(timer, simu, ShareM);
        }
      }timer.beginNew.stop();
    // ! ==================  main loop begin ==================
    }
    return true;
}


#endif
#endif
