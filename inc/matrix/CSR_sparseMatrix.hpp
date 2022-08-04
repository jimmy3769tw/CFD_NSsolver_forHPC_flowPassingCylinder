#pragma once 


// -------------------
#include <vector>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <numeric>
#include <utility> // (since C++11) std::pair
#include <algorithm> // (until C++11)
#include <tuple>
// -------------------

// -------------------
// #include <mpi.h>
#include <omp.h>
// -------------------


/*
! row major
*        | col0 | col1 | col2 | 
*  ______|____________________|
* | row0 |      |      |      | 
* | row1 |      |      |      |
* | row2 |      |      |      |
* | row3 |      |      |      |


! ptr_ [0, ....]
! Idx_ [ ..   ..   .]
! Val_ [ ..   ..   .]
*/


using std::vector, std::cout, std::endl, std::tuple;

namespace mat
{

    template<typename T>
    class CSR_matrix
    {


    using CSR_type = std::tuple< 
            std::vector<int>, std::vector<int>, std::vector<T> >; // <ptr, indices, values>
        public:

        //* ------ Constructor & Destructor ---------

        CSR_matrix(){ }

        CSR_matrix(int rows, int cols) { construct(rows, cols); }

        CSR_matrix(int n) { construct(n, n); }

        CSR_matrix(CSR_type c){ this->set(c); }

        virtual ~CSR_matrix(){ this->destruct(); }

        // *get infomation

        int row() const{return row_;}
        int col() const{return col_;}

        int MaxRow() const;

        void Show_CSR();

        void set(CSR_type csr);

        CSR_type get_CSR();

        // -------------------
        std::vector<T> multiply( std::vector<T> & x) ;
        // -------------------

        // -------------------
        CSR_matrix<T> multiply( CSR_matrix<T> & matrix) ;
        bool multiply( ::vector<T> & x, std::vector<T> & r);
        // -------------------

        // -------------------
        std::vector<T> multiply_omp( std::vector<T> & x) ;
        bool multiply_omp( std::vector<T> & x, std::vector<T> & r);
        // -------------------
        inline std::vector<T> operator * ( std::vector<T> & x) 
        { return this->multiply_omp(x); }

        // * resize

        void resize(int rows, int cols) { this->construct(rows, cols); }

        void resize(int n){ this->construct(n, n); }


        void init(int Idx_len, int prt_len){
            Idx_val_.clear(), ptr_.clear();
            Idx_val_.reserve(Idx_len * 7), ptr_.reserve(prt_len);
            ptr_.push_back(0);
        }

        void push_back(int idx, T val){
            Idx_val_.push_back( {idx, val} );
        }

        void finishIdx(){
            ptr_.push_back( Idx_val_.size() );  
            row_ = col_ = ptr_.size()-1; 
        }

        private:
            // ---------------------
            vector<std::pair<int, T> > Idx_val_;
            
            vector<int> ptr_;
            // ---------------------

            // ---------------------
            int row_ , col_;

            int ompBool = 1;
            // ---------------------

            // ---------------------
            vector<T> re_global_;
            // ---------------------

            // ---------------------
            bool init_pc = true;
            // ---------------------
        
            // ----------- 
            void construct(int rows, int cols);
            void destruct(){}
            // ----------- 
    };

    //* ------ Constructor & Destructor ---------


    // ! ------------------------- construct and resize -------------------------
    template<typename T> inline void CSR_matrix<T>::
    construct(int rows, int cols) {
    // ------------------
        if (rows < 0 || cols < 0) {
			throw std::invalid_argument("Matrix dimensions cannot be negative.");
		}
    // ------------------

    // ------------------
        this -> row_ = rows;

        this -> col_ = cols;

        this -> ptr_.resize(rows+1);
    // ------------------
    }
    // * ----------------------------------------------------------------------


    // ! ----------------------------------------------------------------------
    template<typename T> void CSR_matrix<T>::
    set(CSR_matrix::CSR_type csr) {

    // ------------------
        auto [prt, idx, val] = csr;
        auto len = prt.size()-1;
    // ------------------

    // ------------------
        ptr_ = prt;

        Idx_val_.resize(prt.at(len));

        for (int i = 0; i < ptr_[len] ;++i)
        { Idx_val_.at(i) = std::make_pair(idx.at(i), val.at(i)); }       

        row_ = len;

        col_ = row_;    // CSR can't not received the information regarding the number of cols!
                        // So, it assume the col_ = row_;
    // ------------------
    }
    // * ----------------------------------------------------------------------

    // ! ----------------------------------------------------------------------
    template<typename T> typename CSR_matrix<T>::CSR_type CSR_matrix<T>::
    get_CSR(){
        // -------------------------
        std::vector<int> idxout; idxout.reserve(Idx_val_.size()); 
        std::vector<T> valout;   valout.reserve(Idx_val_.size()); 
        // -------------------------

        // -------------------------
        for (auto x: Idx_val_){
            auto [A, B] = x;
            idxout.push_back(A);
            valout.push_back(B);
        }
        // -------------------------

        // -------------------------
        return make_tuple(ptr_, idxout, valout);
        // -------------------------
    }

    // * ----------------------------------------------------------------------


    //  ! ----------------------------------- Spmv -----------------------------------
    template<typename T> inline bool CSR_matrix<T>::
    multiply(std::vector<T> & x, std::vector<T> & r)
    {

        for(int i = 0; i < row_;++i){
            r[i] = T();
            for(int j = ptr_[i]; j < ptr_[i+1];++j){
                auto [IDX, VAL]  = Idx_val_[j];
                r[i] += VAL * x[IDX];
            }
        }

        return true;
    }


    template<typename T> inline bool CSR_matrix<T>::
    multiply_omp(std::vector<T> & x, std::vector<T> & r)
    {

        #pragma omp parallel for shared(r, x)
        for(int i = 0; i < row_; ++i){
            r[i] = T();
            for(int j = ptr_[i]; j < ptr_[i+1]; ++j){
                auto [IDX, VAL]  = Idx_val_[j];
                r[i] += VAL * x[IDX];
            }
        }

        return true;
    }

    template<typename T> inline std::vector<T> CSR_matrix<T>::
    multiply( std::vector<T> & x)  
    {
        re_global_.resize(row_, T());
        multiply(x, re_global_);
        return re_global_;
    }

    template<typename T> inline std::vector<T> CSR_matrix<T>::
    multiply_omp( std::vector<T> & x) 
    {
        re_global_.resize(row_);
        multiply_omp(x, re_global_ );
        return re_global_;
    }
    //  * ----------------------------------- Spmv -----------------------------------

    //  ! ----------------------------------- Information -------------------
    template<typename T>inline void CSR_matrix<T>::
    Show_CSR()
    {

        // -----------------
        std::cout << "\n|prt\n";
        for (auto v: ptr_){ std::cout << v << ", "; }
        // -----------------

        // -----------------
        std::cout << "\n|Idx\n";
        for (auto v: Idx_val_){
            auto [idx, val] = v;
            std::cout << idx << ", ";
        }
        // -----------------

        // -----------------
        std::cout << "\n|val\n";
        for (auto v: Idx_val_){
            auto [idx, val] = v;
            std::cout << val << ", ";
        }
        // -----------------

        std::cout << std::endl;
    }
    //  * ---------------------------------- -----------------------------------


}
