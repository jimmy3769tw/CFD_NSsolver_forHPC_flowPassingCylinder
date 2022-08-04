#pragma once 

#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <utility>  // (since C++11) std::pair
#include <algorithm> // (until C++11)
#include <tuple>
// ! row major

// #include <mpi.h>
#include <omp.h>

    // *        | col0 | col1 | col2 | 
    // *  ______|____________________|
    // * | row0 |      |      |      | 
    // * | row1 |      |      |      |
    // * | row2 |      |      |      |
    // * | row3 |      |      |      |

#define OMP_ON
#define STD_INNER_PRODUCT

namespace mat
{
    template<typename T>
    class ELL_matrix
    {

        public:
        // ---------------------
        using CSR_type = 
            // <ptr, indices, values>
            std::tuple< std::vector<int>, std::vector<int>, std::vector<T> >; 
        // ---------------------


        //! ------ Constructor & Destructor ---------
        ELL_matrix(){}
        ELL_matrix(int rows, int cols, int MaxCol) { construct(rows, cols, MaxCol); }
        ELL_matrix(int n, int MaxCol) { construct(n, n, MaxCol); }
        virtual ~ELL_matrix(){ destruct(); };

        #ifdef MPI_ON
            ELL_matrix(int rows, int cols, int MaxCol, int PID_, int begin_, int endof_);
            ELL_matrix(int n, int MaxCol, int PID_, int begin_, int endof_); 
        #endif
        //* ------------------------------------------

        // ---------------------------
        void resize(int rows, int cols, int MaxCol){ construct(rows, cols, MaxCol); }
        void resize(int n, int MaxCol) { construct(n, n, MaxCol); }
        // ---------------------------

        // ---------------------------
        int row() const { return rows_; }
        int col() const { return cols_; } 
        int MaxRow()  { return Max_colIdx_; }
        // ---------------------------

        std::vector<std::tuple<int, int, double> > get_mmt();


        // ! ----------------- set and get
        T at(int row, int col) const;
        ELL_matrix<T> & set(int row, int col, T val);
        ELL_matrix<T> & set(CSR_type c);




        void Show_ELL(void);
        // ! ------------------------------- Spmv
        std::vector<T> operator * (const std::vector<T> & x) const
        { return multiply_omp(x); }

        std::vector<T> multiply(const std::vector<T> & x) const;

        std::vector<T> multiply_omp(const std::vector<T> & x) const;
        bool multiply_omp(const std::vector<T> & x, std::vector<T> & r);

        // #if defined(MPI_ON)
        // std::vector<T> multiply_mpi(std::vector<T> & x);
        // #endif
        // * -------------------------------

        // ! ------------------------------- SpmSp
        // ELL_matrix<T> operator * (const ELL_matrix<T> & matrix) const;
        // ELL_matrix<T> multiply(const ELL_matrix<T> & matrix) const;
        // * -------------------------------


        bool analyse();

        //--------------------------------------------- // friend function
        template<typename X>
            friend std::ostream & operator << (std::ostream & os, const ELL_matrix<X> & matrix);
        //---------------------------------------------


        //---------------------------------------------
        void sort_Idx(void);
        CSR_type get_CSR(void);
        //---------------------------------------------
    
        // std::pair<int, double> npc_bicgstab(
        //                             std::vector <double> &B,
        //                             std::vector <double> &x);

        // std::pair<int, double> pc_bicgstab(
        //                             std::vector <double> &B,
        //                             std::vector <double> &x);

        // std::pair<int, double> npc_cg(
        //                             std::vector <double> &B,
        //                             std::vector <double> &x);

        // #ifdef OMP_ON
        // std::pair<int, double> npc_bicgstab_omp(
        //                             std::vector <double> &B,
        //                             std::vector <double> &x);
        // #endif
        
        // #ifdef MPI_ON
        // std::pair<int, double> npc_bicgstab_mpi(
        //                             std::vector <double> &B,
        //                             std::vector <double> &x);

        // void mpi_init(MPI_Comm comm_world, int total);

        // void allocate_vector_mpi(std::vector<T> &x);      

        // #endif


        private:
            int rows_ , cols_, Max_colIdx_;
            std::vector<std::pair<int,T>> Idx_val_;
            std::vector<int> CurrentMaxIdx_;

            // int ompBool = 1;
            // double timer[4] = {0};
            // nomber of rows
            // nomber of cols

            // bool init_pc = true;
            // std::vector<T> I_pcv;


        // #ifdef MPI_ON // !---------------- MPI_ON ------------------

        //     int m_rank_, m_size_, m_begin_, m_endof_, m_length_;

        //     MPI_Comm m_world_;

        //     std::vector<int> m_table_endof_;

        //     std::vector<int> m_table_begin_;

        //     std::vector<int> m_table_length_;

        //     // std::vector<int> m_neighborhood_;


        //     double accumulate_pow_mpi(std::vector<double> &a, std::vector<double> &b);
        //     double inner_product_mpi(std::vector<double> &a, std::vector<double> &b);

        //     std::tuple<int, int> divid(int &ID, const int q,const int r);

        // #endif      // !---------------- MPI_ON ------------------


        // #ifdef OMP_ON
        //     double accumulate_abs_omp_for(std::vector<double> &a, std::vector<double> &b);

        //     double accumulate_omp(std::vector<double> &a, std::vector<double> &b);

        //     double inner_product_omp_for(std::vector<double> &a, std::vector<double> &b);
        // #endif

        // * private function
        void construct(int rows, int cols, int Max_colIdx);


        // ----------- 
        void construct(int Maxoff, int Max_colIdx);
        void destruct(void){}
        // ----------- 


        // ----------- 
        // double  accurate_bicg_ = 1.0e-10;
        // double  accurate_cg_ = 1.0e-10;
        // int iter_max_bicg_ = 3000;
        // int iter_max_cg_ = 3000;
    };



    // ! ------------------------- construct and resize -------------------------

    template<typename T>
    inline void ELL_matrix<T>::construct(int rows, int cols, int Max_colIdx)
    {
        if (rows < 0 || cols < 0) {
			throw std::invalid_argument("Matrix dimensions cannot be negative.");
		}

        this -> rows_ = rows;
        this -> cols_ = cols;
        this -> Max_colIdx_ = Max_colIdx;
        Idx_val_.resize(Max_colIdx*rows, std::make_pair(-1, 0.0));
        CurrentMaxIdx_.resize(rows, 0);
    }

    // * ------------------------------------------------------------------------

    // #ifdef MPI_ON
    // template<typename T>
	// inline ELL_matrix<T>::ELL_matrix(int rows, int cols, int Max_colIdx, int pid, int begin, int endof)
	// {
	// 	this->construct(rows, cols, Max_colIdx);
    //     // this->mpi_init(pid, begin, endof);
	// }

    // template<typename T>
    // inline ELL_matrix<T>::ELL_matrix(int n, int Max_colIdx, int pid, int begin, int endof)
    // {
    //     this->construct(n, n, Max_colIdx);
    // //     this->mpi_init(pid, begin, endof);
    // }

    // template<typename T>
    // inline void ELL_matrix<T>::mpi_init(MPI_Comm comm_world, int total)
    // {
    //     this -> m_world_ = comm_world;

    //     MPI_Comm_rank(comm_world, &this->m_rank_);

    //     MPI_Comm_size(comm_world, &this->m_size_ );

    //     auto [q, r] = std::div(total, this->m_size_);

    //     std::tie(this->m_begin_, this->m_endof_) = this->divid(this->m_rank_, q, r);

    //     this -> m_length_ = this->m_endof_ - this->m_begin_;
        
    //     int n = this->m_size_;

    //     this -> m_table_length_.resize(n);

    //     this -> m_table_begin_.resize(n);
        
    //     this -> m_table_endof_.resize(n);

    //     for (int i = 0; i < n ;++i){
    //         std::tie(this->m_table_begin_[i], this->m_table_endof_[i]) = this->divid(i, q, r);
    //     }

    //     for (int i = 0; i < n ;++i){
    //         this->m_table_length_[i] = this->m_table_endof_[i] - this->m_table_begin_[i];
    //     }


    // }
    // #endif


    // template<typename T>
    // inline std::tuple<int, int> ELL_matrix<T>::divid(
    //         int &ID, const int q, const int r
    // ){
    //     int begin = ID*q;
        
    //     int endof = (ID+1)*q;
        
    //     if (ID < r){
    //         begin += ID;
    //         endof += ID+1;
    //     }
    //     else {
    //         begin += r;
    //         endof += r;
    //     }
        
    //     return std::make_tuple(begin, endof);
    // }






    // ------------------------------------------------------------------
    template<typename T>
    inline std::ostream & operator << (std::ostream & os, const ELL_matrix<T> & matrix){
        // ----------------
        os << std::endl;

        // ----------------

        for(int i = 0, iter; i < matrix.rows_;++i){
            for(int j = 0; j < matrix.cols_;++j){
                const int pt = i*matrix.Max_colIdx_;
                os << std::setw(4);
                for( iter=0; iter < matrix.CurrentMaxIdx_[i] ;++iter ){
                    auto [IDX,VAL] = matrix.Idx_val_.at(pt+iter);
                    if ( j == IDX ){
                        os << VAL;
                        break;
                    }
                }
                if (iter == matrix.CurrentMaxIdx_.at(i)  ){
                    os <<"x" ;
                }
            }
            os << std::endl;
        }
		return os;
    }
    // ------------------------------------------------------------------


    // * ------- set and get
    template<typename T>inline ELL_matrix<T> & ELL_matrix<T>::
    set(ELL_matrix::CSR_type c){
        // --------------
        auto [prt, idx, val] = c;

        // --------------
        rows_ = prt.size()-1;
        // --------------

        for (int row = 0; row < rows_; ++row){
            for (int j = prt[row]; j < prt[row+1] ;++j){
                set(row, idx[j], val[j]);
            }
        }
        // --------------

        cols_ = prt.size()-1;

        // !checking 
        // if (val.size() == prt.at(rows_))
        // {

        // }
        // * --------------

        return *this;
    }


    template<typename T>
    inline ELL_matrix<T> & ELL_matrix<T>::set(int row, int col, T val) {

        if (row < 0 || col < 0) {
			throw std::invalid_argument("Matrix dimensions cannot be negative.");
		}


        if (col >= cols_ ){
            std::cout << col << " >= " << cols_ << std::endl;

			throw std::invalid_argument("colIdx >= Max_colIdx_ ");
        }


        if (row >= rows_)
			throw std::invalid_argument("row >= rows_");


        if (val == T()){
            const int pt = row*this->Max_colIdx_;
            std::cout << "remove" << pt << std::endl;
            // bool NEW = true;
            for(int i = 0; i < CurrentMaxIdx_[row] ;++i ){
                int IDX;
                std::tie(IDX, std::ignore) = this->Idx_val_[pt+i];
                if (col == IDX){
                    // NEW = false; --> remove 
                    for(int j = i ; j < CurrentMaxIdx_[row] || j+1< Max_colIdx_  ;++j){
                       this->Idx_val_[pt+j] = std::make_pair(pt+j+1,pt+j+1); //FIXME?
                    }
                
                    CurrentMaxIdx_[row]--;
                    return *this;
                }
            }
            return *this;
        }
        else{
            const int pt = row*this->Max_colIdx_;
            bool NEW = true;
            for(int i = 0; i < CurrentMaxIdx_[row] ;++i ){
                
                auto [Idx_, val_] = this->Idx_val_[pt+i];
                if ( col == Idx_ ){
                    NEW = false;
                    std::cout << "Replacing" << std::endl;
                    Idx_val_[pt+i] = std::make_pair(col,val);
                    return *this;
                }
            }
        
            if (NEW){
                Idx_val_[pt+CurrentMaxIdx_[row]] = std::make_pair(col,val);
                CurrentMaxIdx_[row]++;
            }
            

            if (CurrentMaxIdx_[row] > this->Max_colIdx_){
                throw std::invalid_argument("[ELL::set] -> CurrentMaxIdx_[row] >this->Max_colIdx_");
            }

            return *this;
        }
    }

    template<typename T>
    inline bool ELL_matrix<T>::analyse(void) {
        std::vector<int> a(Max_colIdx_+1, 0);
        for (auto v:CurrentMaxIdx_){
            ++a[v];
        } 
        int i = 0;
        for (auto v:a){
            std::cout << "i=" << i++ << ", " << v << "\n"; 
        }
        return true;
    }




    template<typename T>
    inline void ELL_matrix<T>::sort_Idx(void){

        for (int row = 0; row < rows_ ; ++row){
            const int pt = row*this->Max_colIdx_;
            for (int i = 0; i < CurrentMaxIdx_[row] ; ++i){
                std::sort(Idx_val_.begin() + pt , Idx_val_.begin()+ pt +CurrentMaxIdx_[row]);
            }
        }
    }



    template<typename T>
    inline typename ELL_matrix<T>::CSR_type ELL_matrix<T>::get_CSR(void){
        sort_Idx();
        std::vector<int> CSR_idxPointer(rows_+1,0);
        std::vector<int> CSR_indices;
        std::vector<T> CSR_data;
        int iter = 0;

        for (int row = 0; row < rows_ ; ++row){
            const int pt = row*Max_colIdx_;
            for (int i = 0; i < CurrentMaxIdx_[row] ; ++i){
                auto [IDX,VAL] = Idx_val_[pt+i];
                CSR_indices.push_back(IDX);
                CSR_data.push_back(VAL);
                ++iter;
            }
            CSR_idxPointer.at(row+1) = iter;
        }

        return std::make_tuple(CSR_idxPointer, CSR_indices, CSR_data);
    }






    template<typename T>
    inline std::vector<std::tuple<int, int , double> > ELL_matrix<T>::get_mmt(){
        this->sort_Idx();

        std::vector<std::tuple<int, int , double> > mmt;


        for (int row = 0; row < rows_ ; ++row){
            const int pt = row*this->Max_colIdx_;
            for (int i = 0; i < CurrentMaxIdx_[row] ; ++i){
                auto [IDX,VAL] = this->Idx_val_[pt+i];
                mmt.push_back(std::tie(row, IDX, VAL));
            }
        }

        return mmt;
    }



    template<typename T>
    inline void ELL_matrix<T>::Show_ELL(void){
        std::cout << "   |val";
        for (int i = 0;i < this->Max_colIdx_ ;++i){ 
              std::cout<< std::setw(4) << "";
        }

        std::cout << "    |Idx\n";

        for (int i = 0;i < rows_ ;++i){   // row
            std::cout << std::setw(3) << i << "|";
            for (int j = 0;j < this->Max_colIdx_ ;++j){  //col
                T VAL;
                std::tie(std::ignore, VAL) = this->Idx_val_.at(i*this->Max_colIdx_+j);
                std::cout<< std::setw(4) << VAL << " " ;
            }
            std::cout << "|";
            for (int j = 0;j < this->Max_colIdx_ ;++j){  //col
                int IDX;
                std::tie(IDX,std::ignore) = this->Idx_val_.at(i*this->Max_colIdx_+j);
                std::cout<< std::setw(4) << IDX << " " ;
            }

            std::cout << std::endl;
        }
    }




    //  ! ---------------------------- Spmv ----------------------------
    template<typename T> 
    inline std::vector<T> ELL_matrix<T>::multiply(const std::vector<T> & x) const
    {
        std::vector<T> r(rows_, T());

        for(int i = 0; i < rows_;++i){

            const int pt = i*Max_colIdx_;
            for(int j = 0; j < CurrentMaxIdx_[i];++j){
                auto [IDX,VAL]  = Idx_val_[pt+j];
                r[i] += VAL * x[IDX];
            }
        }
        return r;
    }

    template<typename T> 
    inline std::vector<T> ELL_matrix<T>::multiply_omp(const std::vector<T> & x) const
    {
        std::vector<T> r(rows_, T());
        multiply_omp(x, r);
        return r;
    }

    template<typename T> 
    inline bool ELL_matrix<T>::multiply_omp(const std::vector<T> & x, std::vector<T> & r) 
    {
        #pragma omp parallel for
        for(int i = 0; i < rows_;++i){
            const int pt = i*Max_colIdx_;
            r[i] = T();
            for(int j = 0; j < CurrentMaxIdx_[i];++j){

                auto [IDX,VAL]  = Idx_val_[pt+j];
                r[i] += VAL * x[IDX];
            }
        }
        return true;
    }


    // template<typename T>
    // inline ELL_matrix<T> multiply(const ELL_matrix<T> & matrix) const
    // {
    //     if( cols_ != matrix.rows_ )
    //     { throw std::invalid_argument("Cannot multiply: col of left matrix != row of right matrix"); }

    //     // get Max Idx

    //     int Max_Idx_new = 0;
    //     for( int i = 0 ; i < rows_ ; ++i)
    //     {
    //         int sum = 0;
    //         for(int j=0; j<matrix.cols_; ++j){
    //             for(int k=0; k<cols_; ++k){
    //                 if ( (this->at(i,k) != T()) && (matrix.at(k,j) != T()) ) 
    //                     sum++; 
    //             }
    //         }

    //         Max_Idx_new= std::max(sum, Max_Idx_new);
    //     }

    //     ELL_matrix<T> result(rows_, matrix.cols_, Max_Idx_new);

    //     T sum;

    //     for(int i=0; i<rows_; ++i)
    //     for(int j=0; j<matrix.cols_; ++j)
    //     {
    //         sum = T();
    //         for(int k=0; k<cols_; ++k)
    //         {
    //             sum += this->at(i,k) * matrix.at(k,j);
    //         }

    //         result.set(i, j, sum);
    //     }

    //     return result;
    // }


    // #ifdef MPI_ON
    // template<typename T> 
    // inline void ELL_matrix<T>::allocate_vector_mpi(std::vector<T> &x)
    // {
    //     for(int i = 0; i < this->m_size_ ; ++i){
    //         MPI_Bcast(  (void *)&x[this->m_table_begin_[i]], 
    //                     this->m_table_length_[i], MPI_DOUBLE, i, this->m_world_);
    //     }

    //     MPI_Barrier(this->m_world_);
    // }


    // template<typename T> 
    // inline std::vector<T> ELL_matrix<T>::multiply_mpi(std::vector<T> & x_Lval_Glen)
    // {
    //     this->allocate_vector_mpi(x_Lval_Glen);

    //     if (cols_ != x_Lval_Glen.size() ){
    //         throw std::invalid_argument("col of matrix != vector.size()");
    //     }

    //     std::vector<T> result_Lval_Glen(rows_, T());

    //     // #pragma omp for
    //     for(int i = this-> m_begin_; i < this->m_endof_ ; ++i){

    //         const int pt = i*this->Max_colIdx_;
    //         // #pragma omp parallel for
    //         for(int j = 0; j < this->CurrentMaxIdx_[i];++j){
    //             auto [IDX,VAL] = this->Idx_val_[pt+j];
    //             result_Lval_Glen[i] += VAL * x_Lval_Glen[IDX];
    //         }
    //     }

    //     MPI_Barrier(this->m_world_);

    //     return result_Lval_Glen;
    // }
    // #endif

    // * ----------------------------------- multiply -----------------------------------

    // //  ! ----------------------------------- accumulate && inner_product -----------------------------------

    // template<typename T>
    // inline double ELL_matrix<T>::inner_product_mpi(std::vector<double> &a, std::vector<double> &b){
    //     #if defined(STD_INNER_PRODUCT)

    //         double temp_l = std::inner_product( a.begin()+this->m_begin_,
    //                                             a.begin()+this->m_endof_,
    //                                             b.begin()+this->m_begin_,0.0);

    //     #else

    //         double temp_l = 0.0;
    //         for(int i = this->m_begin_; i < this->m_endof_ ;++i){
    //             temp_l +=  a[i]*b[i];
    //         }

    //     #endif


    //     double temp_g;
    //     MPI_Allreduce(&temp_l, &temp_g, 1, MPI_DOUBLE, MPI_SUM, this->m_world_);
    //     return temp_g;
    // }


    // template<typename T>
    // inline double ELL_matrix<T>::accumulate_pow_mpi(std::vector<double> &a, std::vector<double> &b){
    //     #if defined(STD_INNER_PRODUCT)

    //         double temp_l = std::inner_product( a.begin()+this->m_begin_,
    //                                             a.begin()+this->m_endof_,
    //                                             b.begin()+this->m_begin_,
    //                                             0.0,std::plus<>(),
    //                                             [](auto &a, auto &b){return std::pow(a-b, 2);}
    //                                             );

    //     #else

    //         double temp_l = 0.0;
    //         for (int i = this->m_begin_ ; i < this->m_endof_ ; ++i )
    //         {
    //             temp_l += std::pow((a[i] - b[i]), 2);
    //         }

    //     #endif

    //     double temp_g;
    //     MPI_Allreduce(&temp_l, &temp_g, 1, MPI_DOUBLE, MPI_SUM, this->m_world_);
    //     return temp_g;

    // }


    // template<typename T>
    // inline double ELL_matrix<T>::inner_product_omp_for(std::vector<double> &a, std::vector<double> &b){
    //     double temp = 0.0;
    //     double time_0 = omp_get_wtime();

    //     #pragma omp parallel for schedule(static) if (ompBool) reduction(+:temp)
    //     for (int i = 0; i<a.size(); ++i ){
    //         temp += a[i] * b[i];
    //     }
    //     timer[0] +=omp_get_wtime()-time_0;

    //     return temp;
    // }

    // template<typename T>
    // inline double ELL_matrix<T>::accumulate_abs_omp_for(std::vector<double> &a, std::vector<double> &b){
    //     double temp = 0.0;

    //     #pragma omp parallel for schedule(static) if (ompBool) reduction(+:temp)
    //     for (int i = 0; i < a.size(); ++i ){
    //         temp += std::abs(a[i] - b[i]);
    //     }

    
    //     return temp;
    // }

    // template<typename T>
    // inline double ELL_matrix<T>::accumulate_omp(std::vector<double> &a, std::vector<double> &b){
    //     double temp = 0.0;

    //     #pragma omp parallel for schedule(static) if (ompBool)  reduction(+:temp)
    //     for (int i = 0; i < a.size(); ++i ){
    //         temp += std::pow((a[i] - b[i]), 2);
    //     }
    //     return temp;
    // }
    //  ! ----------------------------------- accumulate && inner_product -----------------------------------

// # include "ELL_sparseMatrix.BICG.hpp"
// # include "ELL_sparseMatrix.CG.hpp"
// # include "ELL_sparseMatrix.PreBICG.hpp"


    // template<typename T>
    // inline T ELL_matrix<T>::at(int row, int col) const {
    //     if (row < 0 || col < 0) {
	// 		throw std::invalid_argument("Matrix dimensions cannot be negative.");
	// 	}
    //     if (row >= rows_ || col >= cols_ ){
	// 		throw std::invalid_argument("row >= rows_ || colIdx >= Max_colIdx_ ");
    //     }
    //     const int pt = row*this->Max_colIdx_;

    //     for(int i=0; i < CurrentMaxIdx_[row] ;++i ){
    //         auto [IDX,VAL] = this->Idx_val_[pt+i];
    //         if (col == IDX){
    //             return VAL;
    //         }
    //     }
    //     return T();
    //     // throw std::invalid_argument("No declare!!.");
    // }








}
