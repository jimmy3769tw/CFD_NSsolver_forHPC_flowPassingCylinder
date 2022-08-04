#pragma once
#ifndef __MPI_SPARSEMATRIXCSR_H__
#define	__MPI_SPARSEMATRIXCSR_H__

#include <vector>
#include <iostream>
#include <stdexcept>
#include <iomanip>
#include <mpi.h>
#include <omp.h>
#include <math.h> 
#include<numeric>

#include <yakutat/backend/mpi_tools.hpp>
#include <yakutat/backend/Matrix_traits.hpp>
#include <yakutat/mpi/DynamicMatrix.hpp>

#define MPI_SIZE_T MPI_UNSIGNED_LONG

namespace yakutat::mpi
{
    using namespace yakutat::backend;

    template<typename T>
    class SparseMatrixCSR
    {
        // Member types
        using self_type = SparseMatrixCSR<T>;


        public:
             // Member types
            using MatrixType = DynamicMatrix<T>;
            using size_type = typename MatrixTraits<MatrixType>::size_type;
            using value_type = T;
            using array_size_type = typename MatrixTraits<MatrixType>::array_size_type;
            using array_value_type = std::vector<value_type>;
            using data_type = std::tuple<array_size_type, array_size_type, array_value_type>;

        // --- Constructor & Destructor --- matrix partition ---
        //--------------------------------------------- 
            SparseMatrixCSR();
            SparseMatrixCSR(size_type rows, size_type cols); 
            SparseMatrixCSR(size_type n); 
            SparseMatrixCSR(const self_type & matrix); 
            virtual ~SparseMatrixCSR();
            bool partition(size_type myid, size_type nproc);
            bool partition(mpi_tools & world);
        //---------------------------------------------

        //##################################################//
        //               Basic matrix tools                 //
        //##################################################//
        //--------------------------------------------- 
        // --- Get matrix info ---
            size_type row(void) const noexcept;  
            size_type col(void) const noexcept; 
            value_type at(size_type row, size_type col) const;  

        // --- Add element ---
            self_type & set(size_type row, size_type col, value_type val); 
        // --- Accumulate element ---
            self_type & accumulate(size_type row, size_type col, value_type val); 
        // --- reserve matrix ---
            bool reserve(size_type size);
        // --- reserve matrix ---
            size_type capacity(void) const;
        // --- Clear contents ---
            bool clear(void);
        // --- fill zero ---
            bool zero(void);     
        // --- Resize matrix ---
            bool resize(size_type rows, size_type cols); 
            bool resize(size_type n); 
        // --- Output matrix data ---
            data_type data() noexcept;

        //--------------------------------------------- // matrix operators
            // --- mat1 + mat2 ---
            self_type add(const self_type & matrix) const; 
            // --- mat1 - mat2 ---
            self_type subtract(const self_type & matrix) const; 
            // --- mat * vector ---
            array_value_type  multiply(const array_value_type & x) const; 
            // --- mat1 * mat2 ---
            self_type multiply(const self_type & matrix) const;  
        //--------------------------------------------- 

        //--------------------------------------------- // operations
            self_type & operator = (const self_type & matrix);
            self_type operator + (const self_type & matrix) const;
            self_type operator - (const self_type & matrix) const;
            array_value_type operator * (const array_value_type & x) const;
            self_type operator * (const self_type & matrix) const;
            value_type operator () (const size_type & row, const size_type & col) const;
        //---------------------------------------------


        //--------------------------------------------- // friend function
            template<typename X>
				friend std::ostream & operator << (std::ostream & os, const SparseMatrixCSR<X> & matrix);
        //---------------------------------------------


        //--------------------------------------------- // MPI
            self_type & MPI_allgather(size_type itag);
        //---------------------------------------------

        private:
            size_type row_{0}, col_{0}; 

            // csr format
            array_size_type  ptr_;
            array_size_type  indexes_;
            array_value_type values_;

            //MPI variables
            mpi_tools world_;

        //--------------------------------------------- // Basic function
            bool insert(size_type pos, size_type row, size_type col, value_type val);
            bool remove(size_type pos, size_type row);
            bool copy(const self_type & matrix);
        //---------------------------------------------//
    };

    template<typename T>
    inline SparseMatrixCSR<T>::SparseMatrixCSR() {}

    template<typename T>
    inline SparseMatrixCSR<T>::SparseMatrixCSR(size_type rows, size_type cols)
    : row_(rows), col_(cols)
    {
        ptr_.resize(row_+1, 0);
    }

    template<typename T>
    inline SparseMatrixCSR<T>::SparseMatrixCSR(size_type n)
    : row_(n), col_(n)
    {
        ptr_.resize(row_+1, 0);
    }

    template<typename T>
    inline SparseMatrixCSR<T>::SparseMatrixCSR(const self_type & matrix)
    {
        copy(matrix);
    }

    template<typename T>
    inline SparseMatrixCSR<T>::~SparseMatrixCSR() {}

    template<typename T>
    inline bool SparseMatrixCSR<T>::partition(size_type myid, size_type nproc)
    {
        if(row_==0)
        {
            throw std::invalid_argument("SparseMatrixCSR partition: Your row of matrix is zero.");
            return false;
        }
        
        world_.MPI_division(row_, myid, nproc);
        return true;
    }

    template<typename T>
    inline bool SparseMatrixCSR<T>::partition(mpi_tools & world)
    {
        if(world.length==0)
        {
            throw std::invalid_argument("SparseMatrixCSR partition: Your world no division.");
            return false;
        }

        world_ = world;
        return true;
    }



    template<typename T>
    inline typename SparseMatrixCSR<T>::size_type 
    SparseMatrixCSR<T>::row(void) const noexcept
    {
        return row_;
    }

    template<typename T>
    inline typename SparseMatrixCSR<T>::size_type 
    SparseMatrixCSR<T>::col(void) const noexcept
    {
        return col_;
    }

    template<typename T>
    inline typename SparseMatrixCSR<T>::value_type 
    SparseMatrixCSR<T>::at(size_type row, size_type col) const
    {
        if (row < 0 || col < 0 || row > row_-1 || col > col_-1) 
        {
			throw std::invalid_argument("SparseMatrixCSR at: Coordination out of range.");
		}
        size_type currCol{0};
        for(size_type pos = ptr_[row]; pos<ptr_[row+1]; ++pos)
        {
            currCol = indexes_[pos];

            if(currCol == col)
            {
                return values_[pos];
            }
            else if(currCol > col)
            {
                return 0;
            }
        }
        return 0;
    }

    template<typename T>
    inline typename SparseMatrixCSR<T>::self_type & 
    SparseMatrixCSR<T>::set(size_type row, size_type col, value_type val)
    {
        if (row < 0 || col < 0 || row > row_-1 || col > col_-1) 
        {
			throw std::invalid_argument("SparseMatrixCSR set: Coordination out of range.");
		}
        size_type pos{0};
        for(pos=ptr_[row]; pos<ptr_[row+1]; ++pos)
        {
            if(indexes_[pos] >= col)
                break;
        }

        if( !(val==0) )
        {
            if( at(row, col)==0 )
            {
                insert(pos, row, col, val); 
            }
            else
            {
                values_[pos] = val; 
            }
        }
        else if( !(at(row, col)==0) && val==0 )
        {
            remove(pos, row);
        }
		return *this;
    }

    template<typename T>
    inline typename SparseMatrixCSR<T>::self_type & 
    SparseMatrixCSR<T>::accumulate(size_type row, size_type col, value_type val)
    {
        if (row < 0 || col < 0 || row > row_-1 || col > col_-1) 
        {
			throw std::invalid_argument("SparseMatrixCSR accumulate: Coordination out of range.");
		}

        set(row, col, at(row, col) + val);
        return *this;
    }

    template<typename T>
    inline bool SparseMatrixCSR<T>::reserve(size_type size)
    {
        indexes_.reserve(size);
        values_.reserve(size);
        return true;
    }


    template<typename T>
    inline typename SparseMatrixCSR<T>::size_type 
    SparseMatrixCSR<T>::capacity(void) const
    {
        return values_.capacity();
    }


    template<typename T>
    inline bool SparseMatrixCSR<T>::clear(void)
    {
        std::fill(ptr_.begin(), ptr_.end(), 0);
        std::fill(indexes_.begin(), indexes_.end(),0);
        std::fill(values_.begin(), values_.end(),0);
        
        ptr_.clear();
        indexes_.clear();
        values_.clear();
        
        row_ = 0;
        col_ = 0;
        return true;
    }


    template<typename T>
    inline bool SparseMatrixCSR<T>::resize(size_type rows, size_type cols) 
    {
        row_ = rows;
        col_ = cols;
        ptr_.resize(row_+1, 0);
        return true;
    }

    template<typename T>
    inline bool SparseMatrixCSR<T>::resize(size_type n) 
    {
        row_ = n;
        col_ = n;
        ptr_.resize(row_+1, 0);
        return true;
    }


    template<typename T>
    inline bool SparseMatrixCSR<T>::zero(void)     
    {
        std::fill(ptr_.begin(), ptr_.end(), 0);
        std::fill(indexes_.begin(), indexes_.end(),0);
        std::fill(values_.begin(), values_.end(),0);
        return true;
    }


    template<typename T>
    inline typename SparseMatrixCSR<T>::data_type 
    SparseMatrixCSR<T>::data() noexcept
    {
        return {ptr_, indexes_, values_};
    }


    template<typename T>
    inline typename SparseMatrixCSR<T>::self_type 
    SparseMatrixCSR<T>::add(const self_type & matrix) const
    {
        if(row_ != matrix.row_ || col_ != matrix.col_)
        {
            throw std::invalid_argument("SparseMatrixCSR Cannot add: szie of matrices don't match.");
        }
        self_type result(row_, col_);
        result.reserve(capacity());

        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            for(size_type j=0; j<col_; ++j)
            {
               if(at(i,j)!=0 || matrix(i,j)!=0)
                    result.set(i, j, at(i,j) + matrix(i,j) );
            }
        }
        
        return result;
    }

    template<typename T>
    inline typename SparseMatrixCSR<T>::self_type 
    SparseMatrixCSR<T>::subtract(const self_type & matrix) const
    {
        if(row_ != matrix.row_ || col_ != matrix.col_)
        {
            throw std::invalid_argument("SparseMatrixCSR Cannot subtract: szie of matrices don't match.");
        }
        self_type result(row_, col_);
        result.reserve(capacity());

        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            for(size_type j=0; j<col_; ++j)
            {
               if(at(i,j)!=0 || matrix(i,j)!=0)
                    result.set(i, j, at(i,j) - matrix(i,j) );
            }
        }
        
        return result;
    }


    template<typename T>
    inline typename SparseMatrixCSR<T>::array_value_type  
    SparseMatrixCSR<T>::multiply(const array_value_type & x) const
    {
        if( this->row_ != x.size() )
        {
			throw std::invalid_argument("SparseMatrixCSR col of matrix != vector.size()");
        }
        array_value_type result(row_, 0.0);
        value_type sum{0.0};

        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            sum = 0.0;
            for(  size_type j=ptr_[i]; j<ptr_[i+1]; ++j  )
            {
                sum += values_[j] * x[ indexes_[j] ];
            }
            result[i] = sum;
        }
        return result;
    }



    template<typename T>
    inline typename SparseMatrixCSR<T>::self_type
    SparseMatrixCSR<T>::multiply(const self_type & matrix) const
    {
        if( col_ != matrix.row_ )
        {
			throw std::invalid_argument("SparseMatrixCSR Cannot multiply: col of left matrix != row of right matrix");
        }
        self_type result(col_, matrix.row_);
        value_type sum;
        result.reserve(capacity());

        for(size_type i=world_.start; i<world_.end+1; ++i)
        {
            for(size_type j=0; j<col_; ++j)
            {
                sum = 0.0;
                for(  size_type k=ptr_[i]; k<ptr_[i+1]; ++k  )
                {
                    if(matrix.at(indexes_[k], j)!=0)
                        sum += values_[k] * matrix.at(indexes_[k], j);
                }
                result.set(i,j, sum);
            }
            
            
        }

        return result;
    }


    //--------------------------------------------- // operations
    template<typename T>
    inline typename SparseMatrixCSR<T>::self_type & 
    SparseMatrixCSR<T>::operator = (const self_type & matrix)
    {
        if(&matrix != this)
        {
            copy(matrix);
        }
        return *this;
    }

    template<typename T>
    inline typename SparseMatrixCSR<T>::value_type 
    SparseMatrixCSR<T>::operator () (const size_type & row, const size_type & col) const
    {
        return at(row, col);
    }


    template<typename T>
    inline typename SparseMatrixCSR<T>::self_type 
    SparseMatrixCSR<T>::operator + (const self_type & matrix) const
    {
        return add(matrix);
    }

    template<typename T>
    inline typename SparseMatrixCSR<T>::self_type 
    SparseMatrixCSR<T>::operator - (const self_type & matrix) const
    {
        return subtract(matrix);
    }


    template<typename T>
    inline typename SparseMatrixCSR<T>::array_value_type 
    SparseMatrixCSR<T>::operator * (const array_value_type & x) const
    {
        return multiply(x);
    }

    template<typename T>
    inline typename SparseMatrixCSR<T>::self_type 
    SparseMatrixCSR<T>::operator * (const self_type & matrix) const
    {
        return multiply(matrix);
    }
    
    //--------------------------------------------- //

    //--------------------------------------------- //  friend founction

    template<typename T> 
	inline std::ostream & operator << (std::ostream & os, const SparseMatrixCSR<T> & matrix)
	{
        os << std::endl;
        
		for (size_t i = 0; i < matrix.row_; ++i) 
        {
			for (size_t j = 0; j < matrix.col_; ++j) 
            {
				if (j != 0) 
                {
					os << "   ";
				}

				os << std::setprecision(10) << std::setw(10) << matrix.at(i, j);
			}

			
            os << std::endl;
		}
		return os;
	}
    //--------------------------------------------- //




    //--------------------------------------------- // MPI
    template<typename T>
    inline typename SparseMatrixCSR<T>::self_type & 
    SparseMatrixCSR<T>::MPI_allgather(size_type itag)
    {
        MPI_Status istat[8];
        size_type master{0};
        size_type icount{0};
        size_type num_temp, temp, sum{0};
        array_size_type num_of_elem(world_.nproc, 0), idx_temp;
        array_value_type val_temp;
        int long_size = 8;
        MPI_Type_size(MPI_LONG, &long_size);

        // =========== indexes & values =========== //
        // define num of elem vector. Its length is equal to nproc
        // and the val of elements represent num of elements of val, idx per nproc
        num_temp = indexes_.size();
        MPI_Allgather((void *)&num_temp, 1, MPI_LONG, 
        (void *)&num_of_elem[0], 1, MPI_LONG, MPI_COMM_WORLD);

        // cal total length for all nproc
        sum = std::accumulate(num_of_elem.begin(), num_of_elem.end(), 0.0);

        // have all elements for idx, val
        idx_temp.resize(sum, 0); 
        val_temp.resize(sum, 0.0);

        icount = num_of_elem[world_.myid];
        if(world_.myid > master)
        {
            MPI_Send( (void *)&indexes_[0], icount, MPI_LONG, master, itag, MPI_COMM_WORLD);
            MPI_Send( (void *)&values_[0], icount, MPI_DOUBLE, master, itag, MPI_COMM_WORLD);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);

        if(world_.myid == master)
        {
            for(size_type i = 1; i < world_.nproc; ++i)
            {
                MPI_Recv( (void *)&idx_temp[ num_of_elem[i] ], num_of_elem[i] , MPI_LONG, i, itag, MPI_COMM_WORLD, istat);
                MPI_Recv( (void *)&val_temp[ num_of_elem[i] ], num_of_elem[i] , MPI_DOUBLE, i, itag, MPI_COMM_WORLD, istat);
            }
            for(size_type i=0; i<num_of_elem[world_.myid]; ++i)
            {
                idx_temp[i] = indexes_[i];
                val_temp[i] = values_[i];
            }
        }
        indexes_ = idx_temp;
        values_  = val_temp;
        icount = sum;
        MPI_Bcast(  (void *)&indexes_[0], icount, MPI_LONG, master, MPI_COMM_WORLD);
        MPI_Bcast(  (void *)&values_[0], icount, MPI_DOUBLE, master, MPI_COMM_WORLD);
        // =========== indexes & values =========== //


        // =========== ptr =========== //
        for(size_type i=0; i<ptr_.size(); ++i)
        {
            MPI_Allreduce((void *)&ptr_[i], (void *)&temp, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
            ptr_[i] = temp;
        }
        // =========== ptr =========== //


        return *this;
    }
    //--------------------------------------------- //



    template<typename T>
    inline bool 
    SparseMatrixCSR<T>::insert(size_type pos, size_type row, size_type col, value_type val)
    {
        indexes_.insert(indexes_.begin()+pos, col);
        values_.insert(values_.begin()+pos, val);

        for(size_type i=row+1; i<row_+1; ++i)
            ptr_[i] += 1;

        return true;
    }

    template<typename T>
    inline bool 
    SparseMatrixCSR<T>::remove(size_type pos, size_type row)
    {
        indexes_.erase(indexes_.begin()+pos);
        values_.erase(values_.begin()+pos);

        for(size_type i=row+1; i<row_+1; ++i)
            ptr_[i] -= 1;
        return true;
    }

    template<typename T>
    inline bool 
    SparseMatrixCSR<T>::copy(const self_type & matrix)
    {
        row_ = matrix.row_;
        col_ = matrix.col_;
        ptr_ = matrix.ptr_;
        indexes_.reserve(matrix.indexes_.capacity());
        values_.reserve(matrix.values_.capacity());

        indexes_ = matrix.indexes_;
        values_ = matrix.values_;
        
        return true;
    }

} // end namespace yakutat
#endif