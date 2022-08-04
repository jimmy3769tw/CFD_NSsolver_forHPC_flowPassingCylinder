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

using namespace std; 

namespace mat
{
    template<typename T>
    class SPE_matrix
    {

        public:
        using CSR_type = 
                // <ptr, indices, values>
                std::tuple< std::vector<int>,
                std::vector<int>,
                std::vector<T> >;
                // -----------------------

        //* ------ Constructor & Destructor ---------
        SPE_matrix(){}

        SPE_matrix(int nx, int ny, int nz){construct(nx, ny, nz);}

        void resize(int nx, int ny, int nz){construct(nx, ny, nz);}
        
        virtual ~SPE_matrix(){ destruct(); };
    
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
                        auto ic7 = ic*7;
                        auto ijk =  xm(i) + ym(j) + zm(k) + 
                                    xp(i) + yp(j) + zp(k) ;
                        
                        if (i==0)     { ijk -= xm(i); }
                        if (i==nx_-1) { ijk -= xp(i); }
                        if (j==0)     { ijk -= ym(j); }
                        if (j==ny_-1) { ijk -= yp(j); }
                        if (k==0)     { ijk -= zm(k); }
                        if (k==nz_-1) { ijk -= zp(k); }

                        if (i!=0)     { val_[ic7+ 0] = -xm(i); }
                        if (j!=0)     { val_[ic7+ 1] = -yp(j); }
                        if (k!=0)     { val_[ic7+ 2] = -zm(k); }
                                        val_[ic7+ 3] =  ijk;
                        if (k!=nz_-1) { val_[ic7+ 4] = -zp(k); }
                        if (j!=ny_-1) { val_[ic7+ 5] = -ym(j); }
                        if (i!=nx_-1) { val_[ic7+ 6] = -xp(i); }
                    }
                }
            }
        }

        //  ! ----------------------------------- Spmv -----------------------------------

        void multiply_omp(std::vector<double> &x , std::vector<double> &r){
            #pragma omp parallel for
            for (int i = 0; i < ndim_ ; ++i){
                double sum = 0.0;
                auto i7 = i * 7;
                if (i-nynz_ >= 0)       { sum +=  val_[i7+ 0] * x[i-nynz_]; }
                if (i-nz_ >= 0)         { sum +=  val_[i7+ 1] * x[i-nz_  ]; }
                if (i-1 >= 0)           { sum +=  val_[i7+ 2] * x[i-1    ]; }
                                          sum +=  val_[i7+ 3] * x[i      ];
                if (i+1 <  ndim_)       { sum +=  val_[i7+ 4] * x[i+1    ]; }
                if (i+nz_ <  ndim_)     { sum +=  val_[i7+ 5] * x[i+nz_  ]; }
                if (i+nynz_ <  ndim_)   { sum +=  val_[i7+ 6] * x[i+nynz_]; }
                r[i] = sum;
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



