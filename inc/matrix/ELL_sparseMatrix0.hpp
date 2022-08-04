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

#include <mpi.h>
#include <omp.h>

    // *        | col0 | col1 | col2 | 
    // *  ______|____________________|
    // * | row0 |      |      |      | 
    // * | row1 |      |      |      |
    // * | row2 |      |      |      |
    // * | row3 |      |      |      |

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

        // * -------------------------------


        bool analyse();

        //--------------------------------------------- // friend function
        template<typename X>
            friend std::ostream & operator << (std::ostream & os, const ELL_matrix<X> & matrix);
        //---------------------------------------------


        private:
            int rows_ , cols_, Max_colIdx_;
            std::vector<double> Idx_;
            std::vector<double> Val_;
            std::vector<int> CurrentMaxIdx_;

        // * private function
        void construct(int rows, int cols, int Max_colIdx);


        // ----------- 
        void construct(int Maxoff, int Max_colIdx);
        void destruct(void){}
        // ----------- 
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
        Idx_.resize(Max_colIdx*rows, -1);
        Val_.resize(Max_colIdx*rows,  0);
        CurrentMaxIdx_.resize(rows, 0);
    }




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

                    if ( j == matrix.Idx_.at(pt+iter) ){
                        os << matrix.Val_.at(pt+iter);
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

        // ! We assume that your matrix is a n by n matrix for our convienious.
        cols_ = rows_;  


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
            const int pt = row*Max_colIdx_;
            std::cout << "remove" << pt << std::endl;
            // bool NEW = true;
            for(int i = 0; i < CurrentMaxIdx_[row] ;++i ){
                int IDX;
                std::tie(IDX, std::ignore) = dx_val_[pt+i];
                if (col == IDX){
                    // NEW = false; --> remove 
                    for(int j = i ; j < CurrentMaxIdx_[row] || j+1< Max_colIdx_  ;++j){
                        Idx_[pt+j] = pt+j+1;
                        Val_[pt+j] = 0;
    
                    }
                
                    CurrentMaxIdx_[row]--;
                    return *this;
                }
            }
            return *this;
        }
        else{
            const int pt = Max_colIdx_;
            bool NEW = true;
            // -----------------------------------
            for(int i = 0; i < CurrentMaxIdx_[row] ;++i ){
                
                if ( col == Idx_[pt+i] ){
                    NEW = false;
                    std::cout << "Replacing" << std::endl;
                    Idx_[pt+i] = col;
                    Val_[pt+i] = val;
                    return *this;
                }
            }

            // -----------------------------------
            // -----------------------------------
            if (NEW){

                Idx_[pt+CurrentMaxIdx_[row]] = col;
                Val_[pt+CurrentMaxIdx_[row]] = val;
                CurrentMaxIdx_[row]++;
            }
            // -----------------------------------
            

            // -----------------------------------
            if (CurrentMaxIdx_[row] > Max_colIdx_)
            {
                throw std::invalid_argument("[ELL::set] -> CurrentMaxIdx_[row] >Max_colIdx_");
            }
            // -----------------------------------

            return *this;
        }
    }

    template<typename T>
    inline bool ELL_matrix<T>::analyse(void) {
        // -----------------------------------
        std::vector<int> a(Max_colIdx_+1, 0);
        for (auto v:CurrentMaxIdx_){
            ++a[v];
        } 
        // -----------------------------------
        int i = 0;
        // -----------------------------------
        for (auto v:a){
            std::cout << "i=" << i++ << ", " << v << "\n"; 
        }
        // -----------------------------------
        return true;
    }



    template<typename T>
    inline void ELL_matrix<T>::Show_ELL(void){
        std::cout << "   |val";
        for (int i = 0;i < Max_colIdx_ ;++i){ 
              std::cout<< std::setw(4) << "";
        }

        std::cout << "    |Idx\n";

        for (int i = 0;i < rows_ ;++i){

            // ----------------------------------------
            
            std::cout << std::setw(3) << i << "|";
            
            // ----------------------------------------
            for (int j = 0;j < Max_colIdx_ ;++j){  
                std::cout<< std::setw(4) << VAL[i*Max_colIdx_+j] << " " ;
            }
            // ----------------------------------------

            std::cout << "|";

            // ----------------------------------------
            for (int j = 0;j < Max_colIdx_ ;++j){ 
                std::cout<< std::setw(4) << Idx_[i*Max_colIdx_+j] << " " ;
            }
            // ----------------------------------------

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
                r[i] += Val_[pt+j] * x[Idx_[[pt+j]]];
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
                r[i] += Val_[pt+j] * x[Idx_[pt+j]];
            }
        }
        return true;
    }

}
