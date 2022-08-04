#pragma once
#include "0_General.hpp"
// !for Amgcl ================
#include "matrix/ELL_sparseMatrix.hpp"
#include "matrix/CSR_sparseMatrix.hpp"
#include "matrix/SPE_sparseMatrix.hpp"
#include "matrix/SPE_sparseMatrix0.hpp"

// !for Amgcl ================


class MxClass
{
    public:
    mat::ELL_matrix<double> matA;
    mat::CSR_matrix<double> matA_csr;
    mat::SPE_matrix<double> matA_spe;
    std::vector <double> matB;
	std::vector <double> X_result;


#ifndef EIGEN_ON

	// Eigen::VectorXd x_Eigen;

#else

	// std::vector <double> X_result;

#endif



	//* Eigen ================================
#ifdef EIGEN_ON

    std::vector<Eigen::Triplet<double>> tripletList;

    Eigen::SparseMatrix<double, Eigen::RowMajor> matA_Eigen;
    
	Eigen::VectorXd x_Eigen;
    Eigen::VectorXd matB_Eigen;
	
	Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > solver_Eigen;

#endif
	//* Eigen ================================
};


struct SORcoefficient{
    std::vector<double> cf;
};