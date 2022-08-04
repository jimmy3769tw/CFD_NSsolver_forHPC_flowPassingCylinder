
#ifndef __ELL_SPARSEMATRIX_H__
#define	__ELL_SPARSEMATRIX_H__

    #include <vector>
    #include <iostream>
    #include <stdexcept>
    #include <iomanip>
    #include <mpi.h>
    #include <omp.h>
    #include <math.h>
    #include <tuple>
    
namespace mat
{

   


    template<typename T>
    class SparseMatrixELL
    {
            
        public:

            // Member types
            using size_type = size_t;
            // using size_type = int; 
            using value_type = T;
            using array_size_type = std::vector<size_type>;
            using array_value_type = std::vector<value_type>;
            using data_type = std::tuple<array_size_type&, array_size_type&, array_value_type&>; // <ptr, indices, values>

        // --- Constructor & Destructor ---
        //--------------------------------------------- 
            SparseMatrixELL();
            SparseMatrixELL(size_type rows, size_type cols, size_type max_row_nnz); 
            SparseMatrixELL(size_type n, size_type max_row_nnz); 

            // mpi
            SparseMatrixELL(size_type n, size_type max_row_nnz, int myid, int nproc); 
            SparseMatrixELL(size_type rows, size_type cols, size_type max_row_nnz, int myid, int nproc);
            virtual ~SparseMatrixELL();
        //---------------------------------------------



        //##################################################
        //               Basic matrix tools
        //##################################################
        //--------------------------------------------- 
        // --- Get matrix info ---
            size_type row(void) const;  
            size_type col(void) const; 
            size_type nnz(void) const;
            size_type max_row_nnz(void) const;
            T getValue(size_type row, size_type col) const;  

        // --- Add element ---
            SparseMatrixELL<T> & set(size_type row, size_type col, T val); 
        // --- Clear contents ---
            void clear(void);       
        // --- Reset matrix ---
            void resize(size_type rows, size_type cols, size_type max_row_nnz); 
        // --- Copy matrix ---
            SparseMatrixELL(const SparseMatrixELL<T> & matrix); 
        // --- Set row values to zero ---
            void clearRowVal(size_type row_id);
        // --- Set col values to zero ---
            void clearColVal(size_type col_id);
        // --- Get the diagonal value of sparse matrix ---
            T getDiaValue(size_type dia_id);
        // --- Set the diagonal value of sparse matrix ---
            SparseMatrixELL<T> & setDiaValue(size_type dia_id, T val); 
        // --- The element of sparse matrix += value  ---
            void accumulate(size_type row, size_type col, T value);
        // --- Get diagonal matrix of sparse matrix  ---
            SparseMatrixELL<T> diagonal(void) const;
        // --- Transpose the sparse matrix  ---
            SparseMatrixELL<T> transpose(void) const;
        // --- Get the max element of sparse matrix  ---
            T maxCoeff(void) const;
        // --- Get the min element of sparse matrix  ---
            T minCoeff(void) const;
        // --- multiply the value per row ---
            void multiplyVal2Row(size_type row, T value);
        // --- ELL transfer to CSR ---
            void ELL2CSR(void);
        // --- output CSR data ---
            data_type data();

        std::vector<size_type> data_0();
        std::vector<size_type> data_1();
        std::vector<value_type> data_2();


        SparseMatrixELL<T> add(const SparseMatrixELL<T> & matrix) const; // mat1 + mat2
        SparseMatrixELL<T> subtract(const SparseMatrixELL<T> & matrix) const; // mat1 - mat2
        std::vector<T> multiply(const std::vector<T> & x) const; // mat * vector  // mpi
        SparseMatrixELL<T> multiply(const SparseMatrixELL<T> & matrix) const; // mat1 * mat2   


        // --- Show matrix profile ---
            void show(void);
        // --- Show ell_max_idx_ ---
            void show_max_idx(void);

        //---------------------------------------------

        //--------------------------------------------- // operations
            SparseMatrixELL<T> & operator = (const SparseMatrixELL<T> & matrix);
            SparseMatrixELL<T> operator + (const SparseMatrixELL<T> & matrix) const;
            SparseMatrixELL<T> operator - (const SparseMatrixELL<T> & matrix) const;
            std::vector<T> operator * (const std::vector<T> & x) const;
            SparseMatrixELL<T> operator * (const SparseMatrixELL<T> & matrix) const;
        //---------------------------------------------

        //--------------------------------------------- // friend function
            template<typename X>
				friend std::ostream & operator << (std::ostream & os, const SparseMatrixELL<X> & matrix);
        //---------------------------------------------




        private:

        //---------------------------------------------// Local variables
            size_type m_{0}, n_{0},  max_row_nnz_{0}, nnz_{0}; // m:row, n:col
            
           // ell format
            std::vector<int> ell_col_indexes_;
            array_value_type ell_values_;
            array_size_type ell_max_idx_;

            // csr format
            array_size_type ptr_; // ia
            array_size_type indices_; // ja
            array_value_type values_; // aa


            //MPI variables
            int myid_{0}, nproc_{0};
            int start_{0}, end_{0}, count_{0};
            std::vector<int> start_list_, end_list_, count_list_;
            MPI_Status istat_[8];
        //---------------------------------------------

        //--------------------------------------------- // Basic function for SparseMatrixELL class
            void construct(size_type rows, size_type cols, size_type max_row_nnz, int myid, int nproc);
            void construct(size_type rows, size_type cols, size_type max_row_nnz);
            void insert(size_type row, size_type col, T val);
            void csr_insert(size_t index, size_t row, size_t col, T val);
            void remove(size_type row, size_type col);
            void copy(const SparseMatrixELL<T> & m);
            void destruct(void);
            void clc(void);
            void MPI_division(size_type length, int myid, int nproc);
        //---------------------------------------------

    };






    // ======================================================== // Constructor & Destructor
    template<typename T>
	inline SparseMatrixELL<T>::SparseMatrixELL()
	{
		/*empty*/
	}

    template<typename T>
	inline SparseMatrixELL<T>::SparseMatrixELL(size_type rows, size_type cols, size_type max_row_nnz)
	{
		this->construct(rows, cols, max_row_nnz);
	}

    template<typename T>
    inline SparseMatrixELL<T>::SparseMatrixELL(size_type n, size_type max_row_nnz)
    {
        this->construct(n, n, max_row_nnz);
    }

    template<typename T>
    inline SparseMatrixELL<T>::SparseMatrixELL(size_type n, size_type max_row_nnz, int myid, int nproc)
    {
        this->construct(n, n, max_row_nnz, myid, nproc);
        this->MPI_division(n, myid, nproc);
    }

    template<typename T>
	inline SparseMatrixELL<T>::SparseMatrixELL(size_type rows, size_type cols, size_type max_row_nnz, int myid, int nproc)
	{
		this->construct(rows, cols, max_row_nnz, myid, nproc);
        this->MPI_division(rows, myid, nproc);
	}

    template<typename T>
    inline SparseMatrixELL<T>::~SparseMatrixELL(void)
    {
        this->destruct();
    }
    // ======================================================== //

    
    
    



    // ======================================================== //  Basic matrix tools
    template<typename T>
	inline SparseMatrixELL<T> & SparseMatrixELL<T>::set(size_type row, size_type col, T val)
	{

        if (row < 0 || col < 0 || row > this->m_-1 || col > this->n_-1) 
        {
			throw std::invalid_argument("set: Coordination out of range.");
		}

    
        
   

        if( !(val==T()) )
        {
            if( this->getValue(row, col)==T() )
            {
                this->insert(row, col, val); 

            }
            else
            {
                size_type row_idx{0};
                row_idx = row * this->max_row_nnz_;

                for(size_type i=row_idx; i<row_idx+this->max_row_nnz_; ++i)
                {
                    if(ell_col_indexes_[i] == (int)col)
                    {
                        
                        ell_values_[i] = val;
                        break;
                    }
                }
        
            }
        }
        else if( !(this->getValue(row, col)==T()) && val==T() )
        {
            this->remove(row, col);

        }


        
		return *this;
	}



    template<typename T>
	inline T SparseMatrixELL<T>::getValue(size_type row, size_type col) const
	{
        if (row < 0 || col < 0 || row > this->m_-1 || col > this->n_-1) 
        {
			throw std::invalid_argument("getValue: Coordination out of range.");
		}

        size_type row_idx{0};
        T result{0.0};
        row_idx = row * this->max_row_nnz_;

        for(size_type i=row_idx; i<row_idx+this->max_row_nnz_; ++i)
        {
            if(ell_col_indexes_[i] == (int)col)
            {
                result = ell_values_[i];
                break;
            }
        }

       
		

		return result;
	}


    template<typename T>
	inline void SparseMatrixELL<T>::accumulate(size_type row, size_type col, T val)
    {
        if (row < 0 || col < 0 || row > this->m_-1 || col > this->n_-1) 
        {
			throw std::invalid_argument("accumulate: Coordination out of range.");
		}

        this->set(row, col, getValue(row, col) + val);

    }

    template<typename T>
    inline void SparseMatrixELL<T>::clear(void)
    {
        
       this->clc();

    }

    template<typename T>
    inline void  SparseMatrixELL<T>::resize(size_type rows, size_type cols, size_type max_row_nnz)
    {
        this->construct(rows, cols, max_row_nnz);
    }

    
    template<typename T>
	inline SparseMatrixELL<T>::SparseMatrixELL(const SparseMatrixELL<T> & matrix)
	{
		this->copy(matrix);
	}

    template<typename T>
    inline void SparseMatrixELL<T>::clearRowVal(size_type row_id)
    {

        for(size_type j=0; j<this->n_; ++j)
            this->set(row_id, j, 0);
    }

    template<typename T>
    inline void SparseMatrixELL<T>::clearColVal(size_type col_id)
    {
        for(size_type i=0; i<this->m_; ++i)
            this->set(i,col_id, 0);
    }



    
    
    template<typename T>
    inline void SparseMatrixELL<T>::show(void)
    {

        size_type row_idx{0};
        std::cout << "ELL col_indexes" << std::endl;
        for(size_type i=0; i<m_; ++i)
        {
            row_idx = this->max_row_nnz_*i;
            for(size_type j=row_idx; j<row_idx+this->max_row_nnz_; ++j)
            {
                std::cout << std::setw(4) << ell_col_indexes_[j] << " " ;
            }
            std::cout << std::endl;
        }

        row_idx = 0;
        std::cout << "ELL value" << std::endl;
        for(size_type i=0; i<m_; ++i)
        {
            row_idx = this->max_row_nnz_*i;
            for(size_type j=row_idx; j<row_idx+this->max_row_nnz_; ++j)
            {
                std::cout << std::setw(4) << ell_values_[j] << " " ;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;


        
    }
    
    template<typename T>
    inline void SparseMatrixELL<T>::show_max_idx(void)
    {
        std::cout << "show_max_idx" << std::endl;

        for(size_type i=0; i<m_; ++i)
        {   
            std::cout << ell_max_idx_[i] << std::endl;
        }
    }

    template<typename T>
    inline SparseMatrixELL<T> SparseMatrixELL<T>::diagonal(void) const
    {
        size_type m, n, max_row_nnz;
        m = this->m_;
        n = this->n_;
        max_row_nnz = this-> max_row_nnz_;
        SparseMatrixELL<T> result(m, n, 1);

        for(size_type i=0; i<this->m_; ++i)
        {
            result.set(i, i, this->getValue(i,i));
        }
        
        return result;

    }

    template<typename T>
    inline SparseMatrixELL<T> SparseMatrixELL<T>::transpose(void) const
    {
        size_type m, n, max_row_nnz;
        m = this->n_;
        n = this->m_;
        max_row_nnz = this-> max_row_nnz_;
        SparseMatrixELL<T> result(m, n, max_row_nnz);

        for(size_type i=0; i<this->m_; ++i)
        {
            for(size_type j=0; j<this->n_; ++j)
            {
                result.set(j, i, this->getValue(i,j));
            }
        }

        return result;
    }   

    template<typename T>
    inline T SparseMatrixELL<T>::maxCoeff(void) const
    {
        T result;
        for(size_t i=0; i<this->m_; ++i)
        {
            for(size_t j=0; j<this->n_; ++j)
            {
                if(result < this->getValue(i,j))
                    result = this->getValue(i,j);
            }
        }

        return result;
    }

    template<typename T>
    inline T SparseMatrixELL<T>::minCoeff(void) const
    {
        T result;
        for(size_t i=0; i<this->m_; ++i)
        {
            for(size_t j=0; j<this->n_; ++j)
            {
                if(result > this->getValue(i,j))
                    result = this->getValue(i,j);
            }
        }

        return result;
    }


    template<typename T>
    inline typename SparseMatrixELL<T>::size_type SparseMatrixELL<T>::row(void) const
    {
        return this->m_;
    } 

    template<typename T>
    inline typename SparseMatrixELL<T>::size_type SparseMatrixELL<T>::col(void) const
    {
        return this->n_;
    }

    template<typename T>
    inline typename SparseMatrixELL<T>::size_type SparseMatrixELL<T>::nnz(void) const
    {
        size_t cnt{0};
        for(size_t i=0; i<ell_values_.size(); ++i)
        {
            if(ell_values_[i]!=0) cnt++;
        }
        return cnt;
    }

    template<typename T>
    inline typename SparseMatrixELL<T>::size_type SparseMatrixELL<T>::max_row_nnz(void) const
    {
        return this->max_row_nnz_;
    }
   
    

    template<typename T>
    inline T SparseMatrixELL<T>::getDiaValue(size_type dia_id)
    {
        if (dia_id < 0 || dia_id < 0 || dia_id > this->m_-1 || dia_id > this->n_-1) 
        {
			throw std::invalid_argument("dia_id out of range.");
		}

        return this->getValue(dia_id, dia_id);

    }

    template<typename T>
    inline SparseMatrixELL<T> & SparseMatrixELL<T>::setDiaValue(size_type dia_id, T val)
    {
        return this->set(dia_id, dia_id, val);
    }


    template<typename T>
    inline SparseMatrixELL<T> SparseMatrixELL<T>::add(const SparseMatrixELL<T> & matrix) const
    {
        if(this->m_ != matrix.m_ || this->n_ != matrix.n_)
        {
            throw std::invalid_argument("Cannot add: szie of matrices don't match.");
        }

        size_t m, n, max_row_nnz;
        m = this->m_;
        n = this->n_;
        max_row_nnz = this->max_row_nnz_ + matrix.max_row_nnz;

        SparseMatrixELL<T> result(m, n, max_row_nnz);


        for(size_t i=0; i<m; ++i)
        {
            for(size_t j=0; j<n; ++j)
            {
               if(this->getValue(i,j)!=T() ||  matrix.getValue(i,j)!=T())
                    result.set(i, j, this->getValue(i,j) + matrix.getValue(i,j) );


            }
        }
        
        return result;
    }

    template<typename T>
    inline SparseMatrixELL<T> SparseMatrixELL<T>::subtract(const SparseMatrixELL<T> & matrix) const
    {
        if(this->m_ != matrix.m_ || this->n_ != matrix.n_)
        {
            throw std::invalid_argument("Cannot add: szie of matrices don't match.");
        }

        size_t m, n, max_row_nnz;
        m = this->m_;
        n = this->n_;
        max_row_nnz = this-> max_row_nnz_ + matrix.max_row_nnz;

        SparseMatrixELL<T> result(m, n, max_row_nnz);

        

        for(size_t i=0; i<m; ++i)
        {
            for(size_t j=0; j<n; ++j)
            {
               if(this->getValue(i,j)!=T() ||  matrix.getValue(i,j)!=T())
                    result.set(i, j, this->getValue(i,j) - matrix.getValue(i,j) );


            }
        }
        
        return result;
    }

    template<typename T> 
    inline std::vector<T> SparseMatrixELL<T>::multiply(const std::vector<T> & x) const
    {
        if( this->n_ != x.size() )
        {
			throw std::invalid_argument("col of matrix != vector.size()");
        }


        std::vector<T> result(this->m_, T());
        T sum;
        size_t row_idx{0};
        size_t itag{0};
        
        #pragma omp parallel for 
        for(int i=start_; i<end_+1; ++i)
        {
            sum = T();
            row_idx = i * this->max_row_nnz_;

            for(size_t j=row_idx; j<row_idx+ell_max_idx_[i]; ++j)
            {
                sum += ell_values_[j] * x[ ell_col_indexes_[j] ];
            }

            result[i] = sum;
        }

      



        return result;
    }

    template<typename T>
    inline SparseMatrixELL<T> SparseMatrixELL<T>::multiply(const SparseMatrixELL<T> & matrix) const
    {
        if( this->n_ != matrix.m_ )
        {
			throw std::invalid_argument("Cannot multiply: col of left matrix != row of right matrix");
        }
        SparseMatrixELL<T> result(this->m_, matrix.n_, this->max_row_nnz_);
        T sum;

        for(int i=0; i<this->m_; ++i)
        {
            for(int j=0; j<matrix.n_; ++j)
            {
                sum = T();
                for(int k=0; k<this->n_; ++k)
                {
                    sum += this->getValue(i,k) * matrix.getValue(k,j);
                }

                
                result.set(i, j, sum);
            }
            
        }


        return result;
    }
    

    template<typename T>
    inline void SparseMatrixELL<T>::multiplyVal2Row(size_type row, T value)
    {
        size_t row_idx{0};
        row_idx = row * this->max_row_nnz_;

        for(size_t i=row_idx; i<row_idx+this->max_row_nnz_; ++i)
        {
            ell_values_[i] = ell_values_[i] * value;
        }

    }

    template<typename T>
    inline void SparseMatrixELL<T>::ELL2CSR(void)
    {
        ptr_.resize(m_+1, 0);
        indices_.reserve(this->nnz());
        values_.reserve(this->nnz());

        size_t pos{0};
		size_t currCol{0};
        size_t row{0};
        int col{0};


        for(size_t i=0; i<ell_col_indexes_.size(); ++i)
        {
            if(ell_col_indexes_[i] >= 0)
            {
                col = ell_col_indexes_[i];
                row = (i+1)/(this->m_);
                if( (i+1)%(this->m_) == 0 ) row--;

                pos = 0;
                currCol = 0;
                for (pos = ptr_[row]; pos < ptr_[row+1]; ++pos) 
                {
                    currCol = indices_[pos];

                    if (currCol >= col) 
                    {
                        break;
                    }
                }
                csr_insert(pos, row, col, ell_values_[i]);
            }
        }
        

    }

    template<typename T>
    inline typename SparseMatrixELL<T>::data_type SparseMatrixELL<T>::data() 
    {
        return std::tie(ptr_, indices_, values_);
        // return {ptr_, indices_, values_};
    }


    template<typename T>
    inline std::vector<size_t> SparseMatrixELL<T>::data_0() 
    {
        return ptr_;
    }



    template<typename T>
    inline std::vector<size_t> SparseMatrixELL<T>::data_1() 
    {
        return indices_;
    }


    template<typename T>
    inline std::vector<T> SparseMatrixELL<T>::data_2() 
    {
        return values_;
    }


    // ======================================================== //



    // ======================================================== // operations

    template<typename T>
    inline SparseMatrixELL<T> & SparseMatrixELL<T>::operator = (const SparseMatrixELL<T> & matrix)
    {   

        if(&matrix != this)
        {
            this->copy(matrix);
        }
        return *this;
    }

    template<typename T>
	inline SparseMatrixELL<T> SparseMatrixELL<T>::operator + (const SparseMatrixELL<T> & matrix) const
	{
		return this->add(matrix);
	}

    template<typename T>
    inline SparseMatrixELL<T> SparseMatrixELL<T>::operator - (const SparseMatrixELL<T> & matrix) const
    {
        return this->subtract(matrix);
    }

    template<typename T>
    inline std::vector<T> SparseMatrixELL<T>::operator * (const std::vector<T> & x) const
    {
        return this->multiply(x);
    }
    template<typename T>
    inline SparseMatrixELL<T> SparseMatrixELL<T>::operator * (const SparseMatrixELL<T> & matrix) const
    {
        return this->multiply(matrix);
    }
    // ======================================================== //





    // ======================================================== //  friend founction

    template<typename T> 
	inline std::ostream & operator << (std::ostream & os, const SparseMatrixELL<T> & matrix)
	{
        os << std::endl;
        
		for (size_t i = 0; i < matrix.m_; ++i) 
        {
			for (size_t j = 0; j < matrix.n_; ++j) 
            {
				if (j != 0) 
                {
					os << "   ";
				}

				os << std::setw(4) << matrix.getValue(i, j);
			}

			
            os << std::endl;
		}

		return os;
	}


    // ======================================================== //


    // ======================================================== // Basic function for sparse matrix
    template<typename T>
    inline void SparseMatrixELL<T>::construct(size_type rows, size_type cols, size_type max_row_nnz, int myid, int nproc)
    {
        if (rows < 0 || cols < 0) 
        {
			throw "Matrix dimensions cannot be negative.";
		}

        this -> m_ = rows;
        this -> n_ = cols;
        this-> max_row_nnz_ = max_row_nnz;

        this->myid_  = myid;
        this->nproc_ = nproc; 


        ell_col_indexes_.resize(max_row_nnz*rows, -1);
        ell_values_.resize(max_row_nnz*rows, 0);
        ell_max_idx_.resize(rows, 0);

        

    }

    template<typename T>
    inline void SparseMatrixELL<T>::construct(size_type rows, size_type cols, size_type max_row_nnz)
    {
        if (rows < 0 || cols < 0) 
        {
			throw "Matrix dimensions cannot be negative.";
		}

        this -> m_ = rows;
        this -> n_ = cols;
        this-> max_row_nnz_ = max_row_nnz;

        ell_col_indexes_.resize(max_row_nnz*rows, -1);
        ell_values_.resize(max_row_nnz*rows, 0);
        ell_max_idx_.resize(rows, 0);
        

    }

    template<typename T>
    inline void SparseMatrixELL<T>::insert(size_type row, size_type col, T val)
    {
        size_t row_idx{0};
        row_idx = row * this->max_row_nnz_;

        for(size_t i=row_idx; i<row_idx+this->max_row_nnz_; ++i)
        {
            if(ell_col_indexes_[i] < 0)
            {
                ell_col_indexes_[i] = col;
                ell_values_[i] = val;
                ell_max_idx_[row]++;
                break;
            }
        }
    }

    template<typename T>
    inline void SparseMatrixELL<T>::csr_insert(size_t index, size_t row, size_t col, T val)
    {
     
        values_.insert(values_.begin()+index, val);
        indices_.insert(indices_.begin()+index, col); 

        for (size_t i = row+1; i < this->m_+1; ++i) 
        {
			ptr_[i] += 1;
		}
        

    }

    template<typename T>
	inline void SparseMatrixELL<T>::remove(size_type row, size_type col)
	{
        size_t row_idx{0};
        row_idx = row * this->max_row_nnz_;

        for(size_t i=row_idx; i<row_idx+this->max_row_nnz_; ++i)
        {
            if(ell_col_indexes_[i] == (int)col)
            {
                ell_col_indexes_[i] = -1;
                ell_values_[i] = T();
                ell_max_idx_[row]--;
                break;
            }
        }

		
	}

    template<typename T>
    inline void  SparseMatrixELL<T>::copy(const SparseMatrixELL<T> & matrix)
    {
        this->m_ = matrix.m_;
        this->n_ = matrix.n_;
        this->max_row_nnz_ = matrix.max_row_nnz_;



       
        ell_col_indexes_.assign(matrix.ell_col_indexes_.begin(), matrix.ell_col_indexes_.end());
        ell_values_.assign(matrix.ell_values_.begin(), matrix.ell_values_.end());
        
        
    }

    template<typename T>
    inline void SparseMatrixELL<T>::destruct(void)
    {
   
    }

    template<typename T>
    inline void SparseMatrixELL<T>::clc(void)
    {

        ell_values_.clear();
        ell_col_indexes_.clear();

    }





    template<typename T>
    void SparseMatrixELL<T>::MPI_division(size_type length, int myid, int nproc)
    {
        int Xdv = length / nproc;				    
        int Xr = length - Xdv * nproc;

        start_list_.reserve(nproc_);
        end_list_.reserve(nproc_);
        count_list_.reserve(nproc_);

        for(int i=0; i<nproc; ++i)
        {
            if(i < Xr)
            {
                start_list_.push_back( i * (Xdv + 1) + 0 );			
                end_list_.push_back( start_list_[i]+ Xdv );
            }
            else
            {
                start_list_.push_back( i * Xdv + Xr + 0 );		

                end_list_.push_back( start_list_[i] + Xdv - 1 );		

            }

            count_list_.push_back( end_list_[i] - start_list_[i] + 1 );	

        }
        start_ = start_list_[myid];
        end_ = end_list_[myid];
        count_ = count_list_[myid];
    }

    // ======================================================== //


} // end namespace mat


#endif
