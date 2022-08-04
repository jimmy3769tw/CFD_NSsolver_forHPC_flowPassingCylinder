#pragma once 

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>

#include <numeric>
#include <utility> // (since C++11) std::pair
#include <algorithm> // (until C++11)
#include <tuple>
// ! row major
// #include <mpi.h>
#include <omp.h>

namespace mat
{

    template<typename T>
    class SPE_matrix0
    {

        public:

        using CSR_type = 
                // <ptr, indices, values>
                std::tuple< std::vector<int>,
                std::vector<int>,
                std::vector<T> >;
                // -----------------------

        //* ------ Constructor & Destructor ---------
        SPE_matrix0(){}

        SPE_matrix0(int nx, int ny, int nz){construct(nx, ny, nz);}

        void resize(int nx, int ny, int nz){construct(nx, ny, nz);}
        
        virtual ~SPE_matrix0(){ destruct(); };
    
        void setupPressure(
            int gC, 
            const std::vector<double> &Dx, 
            const std::vector<double> &Dy, 
            const std::vector<double> &Dz,
            const std::vector<double> &Dxs, 
            const std::vector<double> &Dys, 
            const std::vector<double> &Dzs 
        ){

            auto xp = [&](auto i){return 1.0/Dx.at(i+gC)/Dxs.at(i+gC);};
            auto yp = [&](auto j){return 1.0/Dy.at(j+gC)/Dys.at(j+gC);};
            auto zp = [&](auto k){return 1.0/Dz.at(k+gC)/Dzs.at(k+gC);};

            auto xm = [&](auto i){return 1.0/Dx.at(i+gC)/Dxs.at(i+gC-1);};
            auto ym = [&](auto j){return 1.0/Dy.at(j+gC)/Dys.at(j+gC-1);};
            auto zm = [&](auto k){return 1.0/Dz.at(k+gC)/Dzs.at(k+gC-1);};

            for (int i=0; i<nx_; ++i)
            {
                for (int j=0; j<ny_; ++j)
                {
                    for (int k=0; k<nz_; ++k)
                    {
                        int ic = k + j*nz_ + i*nynz_;
                        auto ijk =  xm(i) + ym(j) + zm(k) + 
                                    xp(i) + yp(j) + zp(k) ;
                        
                        if (i==0)     { ijk -= xm(i); }
                        if (i==nx_-1) { ijk -= xp(i); }
                        if (j==0)     { ijk -= ym(j); }
                        if (j==ny_-1) { ijk -= yp(j); }
                        if (k==0)     { ijk -= zm(k); }
                        if (k==nz_-1) { ijk -= zp(k); }

                        if (i!=0)     { val_[ic+0*ndim_] = -xm(i); }
                        if (j!=0)     { val_[ic+1*ndim_] = -yp(j); }
                        if (k!=0)     { val_[ic+2*ndim_] = -zm(k); }
                                        val_[ic+3*ndim_] =  ijk;
                        if (k!=nz_-1) { val_[ic+4*ndim_] = -zp(k); }
                        if (j!=ny_-1) { val_[ic+5*ndim_] = -ym(j); }
                        if (i!=nx_-1) { val_[ic+6*ndim_] = -xp(i); }
                    }
                }
            }
        }

        //  ! ----------------------------------- Spmv -----------------------------------

        void multiply_omp(std::vector<double> &x , std::vector<double> &r){
            #pragma omp parallel for
            for (int ic = 0; ic < ndim_ ; ++ic){
                double sum = 0.0;
                if (ic-nynz_ >= 0)       { sum +=  val_[ic+0*ndim_] * x[ic-nynz_]; }
                if (ic-nz_ >= 0)         { sum +=  val_[ic+1*ndim_] * x[ic-nz_  ]; }
                if (ic-1 >= 0)           { sum +=  val_[ic+2*ndim_] * x[ic-1    ]; }
                                           sum +=  val_[ic+3*ndim_] * x[ic      ];
                if (ic+1 <  ndim_)       { sum +=  val_[ic+4*ndim_] * x[ic+1    ]; }
                if (ic+nz_ <  ndim_)     { sum +=  val_[ic+5*ndim_] * x[ic+nz_  ]; }
                if (ic+nynz_ <  ndim_)   { sum +=  val_[ic+6*ndim_] * x[ic+nynz_]; }
                r[ic] = sum;
            }
        }

        int row() const {return ndim_;}

        private:

            int ndim_, nx_, ny_, nz_, nynz_;

            std::vector<double> val_;

            void construct(int nx, int ny, int nz){ 

                ndim_ = nx*ny*nz;
                nx_ = nx;
                ny_ = ny;
                nz_ = nz;
                nynz_ = ny*nz;
                val_.resize(ndim_*7);
            }

            void destruct(void){}

         

    };

}



