#pragma once 
// ! row major

/*

*        | col0 | col1 | col2 | 
*  ______|____________________|
* | row0 |      |      |      | 
* | row1 |      |      |      |
* | row2 |      |      |      |
* | row3 |      |      |      |

*/

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
#include <mpi.h>
#include <omp.h>
// -------------------
using std::vector, std::cout, std::endl, std::tuple;

namespace mat
{

    template<typename T>
    class CSR_matrix00
    {


    using CSR_type = std::tuple< 
            std::vector<int>, std::vector<int>, std::vector<T> >; // <ptr, indices, values>
        public:

        //* ------ Constructor & Destructor ---------

        CSR_matrix0(){ }

        CSR_matrix0(int rows, int cols) { construct(rows, cols); }

        CSR_matrix0(int n) { construct(n, n); }

        CSR_matrix0(CSR_type c){ this->set(c); }

        virtual ~CSR_matrix0(){ this->destruct(); }

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
        CSR_matrix0<T> multiply( CSR_matrix0<T> & matrix) ;
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


        void init(int Idx_len, int ptr_len){

            Idx_.clear(), Val_.clear(), ptr_.clear();

            Idx_.reserve(Idx_len), Val_.reserve(Idx_len), ptr_.reserve(ptr_len);

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
            vector<int> ptr_;
            vector<int> Idx_;
            vector<T> Val_;
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
        
            // ----------- construct & destruct
            void construct(int rows, int cols);
            void destruct(){}
            // ----------- 
    };

    //* ------ Constructor & Destructor ---------


    // ! ------------------------- construct and resize -------------------------
    template<typename T> inline void CSR_matrix0<T>::
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
    template<typename T> void CSR_matrix0<T>::
    set(CSR_matrix0::CSR_type csr) {

    // ------------------
        std::tie(ptr_, Idx_, Val_) = csr;
        auto len = ptr_.size()-1;
    // ------------------

    // ------------------

        row_ = len;

        col_ = row_;    // CSR can't not received the information regarding the number of cols!
                        // So, it assume the col_ = row;
    // ------------------
    }
    // * ----------------------------------------------------------------------

    // ! ----------------------------------------------------------------------
    template<typename T> typename CSR_matrix0<T>::CSR_type CSR_matrix0<T>::
    get_CSR(){ return make_tuple(ptr_, Idx_, Val_); }
    // * ----------------------------------------------------------------------


    //  ! ----------------------------------- Spmv -----------------------------------
    template<typename T> inline bool CSR_matrix0<T>::
    multiply(std::vector<T> & x, std::vector<T> & r)
    {

        for(int i = 0; i < row_;++i){

            r[i] = T();
            for(int j = ptr_[i]; j < ptr_[i+1];++j)
            {  r[i] += Val_[j] * x[Idx_[j]];}
        }

        return true;
    }


    template<typename T> inline bool CSR_matrix0<T>::
    multiply_omp(std::vector<T> & x, std::vector<T> & r)
    {

        #pragma omp parallel for shared(r, x)
        for(int i = 0; i < row_; ++i){
            r[i] = T();
            r[i] = T();
            for(int j = ptr_[i]; j < ptr_[i+1];++j)
            {  r[i] += Val_[j] * x[Idx_[j]];}
        }

        return true;
    }

    template<typename T> inline std::vector<T> CSR_matrix0<T>::
    multiply( std::vector<T> & x)  {
        re_global_.resize(row_, T());
        multiply(x, re_global_);
        return re_global_;
    }

    template<typename T> inline std::vector<T> CSR_matrix0<T>::
    multiply_omp( std::vector<T> & x) 
    {
        re_global_.resize(row_);
        multiply_omp(x, re_global_ );
        return re_global_;
    }
    //  * ----------------------------------- Spmv -----------------------------------

    //  ! ----------------------------------- Information -------------------
    template<typename T>inline void CSR_matrix0<T>::
    Show_CSR(){

        // -----------------
        std::cout << "\n|prt\n";
        for (auto v: ptr_){ std::cout << v << ", "; }
        // -----------------

        // -----------------
        std::cout << "\n|Idx\n";
        for (auto v: Idx_){
            std::cout << v << ", ";
        }
        // -----------------

        // -----------------
        std::cout << "\n|val\n";
        for (auto v: Val_){
            std::cout << v << ", ";
        }
        // -----------------

        std::cout << std::endl;
    }
    //  * ---------------------------------- -----------------------------------


}