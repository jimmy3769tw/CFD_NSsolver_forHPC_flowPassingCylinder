
    // ---------------------------------------- EIGNE (AMGCL)
    #elif defined (P_SOLVER_EIGEN_CSR)

      createBMatrix_Eigen(T1, Mx, simu, Lo, gridA);
    
      Mx.solver_Eigen.compute(Mx.matA_Eigen);
    
      Mx.solver_Eigen.setTolerance(1e-15);
    
      Mx.x_Eigen = Mx.solver_Eigen.solve(Mx.matB_Eigen);

      Pressure_transform_x_Eigen(t1, Mx, Lo, gridA);

    // ---------------------------------------- EIGEN (AMGCL)
    #elif defined (P_SOLVER_AMGCL_EIGEN)

      createBMatrix_Eigen(T1, Mx, simu, Lo, gridA);

      createBMatrix_seq(T1, Mx, simu, Lo, gridA);

      SolverEigen solve(Mx.matA_Eigen);

      std::tie(simu.iters, simu.error) = solve(Mx.matB_Eigen, Mx.x_Eigen);

      Pressure_transform_x_Eigen(t1, Mx, Lo, gridA);

      simu.iters = iters;
