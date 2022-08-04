#ifndef __HYB2_MATRIX_H__
#define	__HYB2_MATRIX_H__

    #include <vector>
    #include <iostream>
    #include <stdexcept>
    #include <iomanip>

namespace mat
{
    template<typename T>
    class SparseMatrixHYB
    {

        public:
        
        // --- Constructor & Destructor ---
        //--------------------------------------------- 
            SparseMatrixHYB(size_t n, size_t min_row_nnz); 
            SparseMatrixHYB(size_t rows, size_t cols, size_t min_row_nnz);
            virtual ~SparseMatrixHYB();
        //---------------------------------------------



        //##################################################
        //               Basic matrix tools
        //##################################################
        //--------------------------------------------- 
        // --- Get matrix info ---
            size_t row(void) const;  
            size_t col(void) const; 
            size_t nnz(void) const;
            T getValue(size_t row, size_t col) const;  

        // --- Add element ---
            SparseMatrixHYB<T> & set(size_t row, size_t col, T val); 
        // --- Clear contents ---
            void clear(void);       
        // --- Reset matrix ---
            void resize(size_t rows, size_t cols, size_t min_row_nnz); 
        // --- Copy matrix ---
            SparseMatrixHYB(const SparseMatrixHYB<T> & matrix); 
        // --- Set row values to zero ---
            void clearRowVal(size_t row_id);
        // --- Set col values to zero ---
            void clearColVal(size_t col_id);
        // --- Get the diagonal value of sparse matrix ---
            T getDiaValue(size_t dia_id);
        // --- Set the diagonal value of sparse matrix ---
            SparseMatrixHYB<T> & setDiaValue(size_t dia_id, T val); 
        // --- The element of sparse matrix += value  ---
            void accumulate(size_t row, size_t col, T value);
        // --- Get diagonal matrix of sparse matrix  ---
            SparseMatrixHYB<T> diagonal(void) const;
        // --- Transpose the sparse matrix  ---
            SparseMatrixHYB<T> transpose(void) const;
        // --- Get the max element of sparse matrix  ---
            T maxCoeff(void) const;
        // --- Get the min element of sparse matrix  ---
            T minCoeff(void) const;
        // --- Get the nnz per row of sparse matrix  ---
            std::vector<T> NNzPerRow(void) const;
        // --- Get the nnz per col of sparse matrix  ---
            std::vector<T> NNzPerCol(void) const;


        SparseMatrixHYB<T> add(const SparseMatrixHYB<T> & matrix) const; // mat1 + mat2
        SparseMatrixHYB<T> subtract(const SparseMatrixHYB<T> & matrix) const; // mat1 - mat2
        std::vector<T> multiply(const std::vector<T> & x) const; // mat * vector 
        SparseMatrixHYB<T> multiply(const SparseMatrixHYB<T> & matrix) const; // mat1 * mat2 


        // --- Show matrix profile ---
            void show(void);


        //---------------------------------------------

        //--------------------------------------------- // operations
            SparseMatrixHYB<T> & operator = (const SparseMatrixHYB<T> & matrix);
            SparseMatrixHYB<T> operator + (const SparseMatrixHYB<T> & matrix) const;
            SparseMatrixHYB<T> operator - (const SparseMatrixHYB<T> & matrix) const;
            std::vector<T> operator * (const std::vector<T> & x) const;
            SparseMatrixHYB<T> operator * (const SparseMatrixHYB<T> & matrix) const;
        //---------------------------------------------

        //--------------------------------------------- // friend function
            template<typename X>
				friend std::ostream & operator << (std::ostream & os, const SparseMatrixHYB<X> & matrix);

            template<typename X>
				friend bool operator == (const SparseMatrixHYB<X> & a, const SparseMatrixHYB<X> & b);

            template<typename X>
				friend bool operator != (const SparseMatrixHYB<X> & a, const SparseMatrixHYB<X> & b);

            
        //---------------------------------------------




        private:

        //---------------------------------------------// Local variables
            size_t m_{0}, n_{0},  min_row_nnz_{0}, nnz_{0}; // m:row, n:col


            std::vector<size_t> *ia_csr_, *ja_csr_;  //csr_rows cols list 

            std::vector<int> ell_col_indexes_;
            std::vector<T> ell_values_;
        //---------------------------------------------

        //--------------------------------------------- // Basic function for SparseMatrixHYB class
            void construct(size_t rows, size_t cols, size_t min_row_nnz);
            void insert(size_t index, size_t row, size_t col, T val);
            void remove(size_t index, size_t row, size_t col);
            void copy(const SparseMatrixHYB<T> & m);
            void destruct(void);
            void clc(void);
            size_t getCsrRowList(size_t pos) const;
            size_t getColumnList(size_t pos) const;
            T getValueList(size_t pos) const;
        //---------------------------------------------

    };






    // ======================================================== // Constructor & Destructor
    template<typename T>
    inline SparseMatrixHYB<T>::SparseMatrixHYB(size_t n, size_t min_row_nnz)
    {
        this->construct(n, n, min_row_nnz);
    }

    template<typename T>
	inline SparseMatrixHYB<T>::SparseMatrixHYB(size_t rows, size_t cols, size_t min_row_nnz)
	{
		this->construct(rows, cols, min_row_nnz);
	}

    template<typename T>
    inline SparseMatrixHYB<T>::~SparseMatrixHYB(void)
    {
        this->destruct();
    }
    // ======================================================== //

    
    
    



    // ======================================================== //  Basic matrix tools
    template<typename T>
	inline SparseMatrixHYB<T> & SparseMatrixHYB<T>::set(size_t row, size_t col, T val)
	{

        if (row < 0 || col < 0 || row > this->m_-1 || col > this->n_-1) 
        {
			throw std::invalid_argument("Coordination out of range.");
		}

       
		size_t pos = 0;
		size_t currCol = 0;
        
		for (pos = (*(this->ia_csr_))[row]; pos < (*(this->ia_csr_))[row+1] ; ++pos) 
        {
			currCol = (*(this->ja_csr_))[pos];

			if (currCol >= col) 
            {
				break;
			}
		}




        if( !(val==T()) )
        {
            if( this->getValue(row, col)==T() )
            {
                this->insert(pos, row, col, val); 

            }
            else
            {
                size_t row_idx{0};
                size_t flag{0};
                row_idx = row * this->min_row_nnz_;

                for(size_t i=row_idx; i<row_idx+this->min_row_nnz_; ++i)
                {
                    if(ell_col_indexes_[i] = (int)col)
                    {
                        ell_values_[i] = val;
                        flag = 1;
                        break;
                    }
                }
                if(flag == 0)
                    (*(this->aa_csr_))[pos] = val; 
            }
        }
        else if( !(this->getValue(row, col)==T()) && val==T() )
        {
            this->remove(pos, row, col);

        }


        
		return *this;
	}



    template<typename T>
	inline T SparseMatrixHYB<T>::getValue(size_t row, size_t col) const
	{
        if (row < 0 || col < 0 || row > this->m_-1 || col > this->n_-1) 
        {
			throw std::invalid_argument("Coordination out of range.");
		}

        size_t row_idx{0};
        size_t flag{0};

        row_idx = row * this->min_row_nnz_;

        for(size_t i=row_idx; i<row_idx+this->min_row_nnz_; ++i)
        {
            if(ell_col_indexes_[i] == (int)col)
            {
                flag = 1;
                return ell_values_[i];
            }
        }

        if(flag == 0)
        {
            size_t currCol = 0;
            // example   input=1,1
            // rows  = {0,         3,6  ,9}. 
            // cols = {0,1,2,     0,1,3  ,0,2,3,1,2,3}
            // vals = {4,-1,-1,  -1,4,-1  ,-1,4,-1,-1,-1,4}
            for (size_t pos = (*(this->ia_csr_))[row]; pos < (*(this->ia_csr_))[row+1]; ++pos)  //pos=3,4,5
            {
                currCol = (*(this->ja_csr_))[pos];  

                if (currCol == col) //when pos = 4 currCol = cols[4] = col = 1. 4 means number 4 of nnz value
                {
                    return (*(this->aa_csr_))[pos]; 

                } 
                else if (currCol > col) 
                {
                    break;
                }
            }
        }

		

		return T();
	}


    template<typename T>
	inline void SparseMatrixHYB<T>::accumulate(size_t row, size_t col, T val)
    {
        if (row < 0 || col < 0 || row > this->m_-1 || col > this->n_-1) 
        {
			throw std::invalid_argument("Coordination out of range.");
		}

        this->set(row, col, getValue(row, col) + val);

    }

    template<typename T>
    inline void SparseMatrixHYB<T>::clear(void)
    {
        
       this->clc();

    }

    template<typename T>
    inline void  SparseMatrixHYB<T>::resize(size_t rows, size_t cols, size_t min_row_nnz)
    {
        this->construct(rows, cols, min_row_nnz);
    }

    
    template<typename T>
	inline SparseMatrixHYB<T>::SparseMatrixHYB(const SparseMatrixHYB<T> & matrix)
	{
		this->copy(matrix);
	}

    template<typename T>
    inline void SparseMatrixHYB<T>::clearRowVal(size_t row_id)
    {

        for(size_t j=0; j<this->n_; ++j)
            this->set(row_id, j, 0);
    }

    template<typename T>
    inline void SparseMatrixHYB<T>::clearColVal(size_t col_id)
    {
        for(size_t i=0; i<this->m_; ++i)
            this->set(i,col_id, 0);
    }



    
    
    template<typename T>
    inline void SparseMatrixHYB<T>::show(void)
    {

        size_t row_idx{0};
        std::cout << "ELL col_indexes" << std::endl;
        for(size_t i=0; i<m_; ++i)
        {
            row_idx = this->min_row_nnz_*i;
            for(size_t j=row_idx; j<row_idx+this->min_row_nnz_; ++j)
            {
                std::cout << std::setw(4) << ell_col_indexes_[j] << " " ;
            }
            std::cout << std::endl;
        }

        row_idx = 0;
        std::cout << "ELL value" << std::endl;
        for(size_t i=0; i<m_; ++i)
        {
            row_idx = this->min_row_nnz_*i;
            for(size_t j=row_idx; j<row_idx+this->min_row_nnz_; ++j)
            {
                std::cout << std::setw(4) << ell_values_[j] << " " ;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;





        for(size_t i=0; i<this->m_+1; ++i)
        {
            std::cout << std::setw(4) << this->getCsrRowList(i) << std::endl;
        }
        std::cout <<"================="<<std::endl;
        std::cout << std::endl;
        if(this->aa_csr_ != NULL)
            for(size_t i=0; i<this->nnz(); ++i)
            {
                std::cout << std::setw(4) << this->getColumnList(i) << " "
                          << std::setw(4) << this->getValueList(i)  << std::endl;
            }
        
    }
    


    template<typename T>
    inline SparseMatrixHYB<T> SparseMatrixHYB<T>::diagonal(void) const
    {
        size_t m, n;
        m = this->m_;
        n = this->n_;
        SparseMatrixHYB<T> result(m, n);

        for(size_t i=0; i<this->m_; ++i)
        {
            for(size_t j=0; j<this->n_; ++j)
            {
                if(i == j)
                {
                    result.set(i, j, this->getValue(i,j));
                }
            }
        }
        
        return result;

    }

    template<typename T>
    inline SparseMatrixHYB<T> SparseMatrixHYB<T>::transpose(void) const
    {
        size_t m, n, z;
        m = this->n_;
        n = this->m_;
        SparseMatrixHYB<T> result(m, n);

        for(size_t i=0; i<this->m_; ++i)
        {
            for(size_t j=0; j<this->n_; ++j)
            {
                result.set(j, i, this->getValue(i,j));
            }
        }

        return result;
    }   

    template<typename T>
    inline T SparseMatrixHYB<T>::maxCoeff(void) const
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
    inline T SparseMatrixHYB<T>::minCoeff(void) const
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
    inline size_t SparseMatrixHYB<T>::row(void) const
    {
        return this->m_;
    } 

    template<typename T>
    inline size_t SparseMatrixHYB<T>::col(void) const
    {
        return this->n_;
    }

    template<typename T>
    inline size_t SparseMatrixHYB<T>::nnz(void) const
    {
        return (*(this->ia_csr_))[this->row()];
    }


   
    

    template<typename T>
    inline T SparseMatrixHYB<T>::getDiaValue(size_t dia_id)
    {
        if (dia_id < 0 || dia_id < 0 || dia_id > this->m_-1 || dia_id > this->n_-1) 
        {
			throw std::invalid_argument("dia_id out of range.");
		}

        return this->getValue(dia_id, dia_id);

    }

    template<typename T>
    inline SparseMatrixHYB<T> & SparseMatrixHYB<T>::setDiaValue(size_t dia_id, T val)
    {
        return this->set(dia_id, dia_id, val);
    }

    template<typename T>
    inline std::vector<T> SparseMatrixHYB<T>::NNzPerRow(void) const
    {
        std::vector<T> result;

        result.resize(this->m_);


        for(size_t row=0; row<this->m_; ++row)
        {
            result[row] = 0;
            for (size_t pos = (*(this->ia_csr_))[row]; pos < (*(this->ia_csr_))[row+1]; ++pos)
            {
                ++result[row];
            }
        }

        return result;
        
    }


    template<typename T>
    std::vector<T> SparseMatrixHYB<T>::NNzPerCol(void) const
    {
        std::vector<T> result;
        result.resize(this->n_);
        size_t icol=0;

        for(size_t i=0; i<this->nnz(); ++i)
        {
            icol = (*(this->ja_csr_))[i];
            ++result[icol];
        }

        return result;
    }


    template<typename T>
    inline SparseMatrixHYB<T> SparseMatrixHYB<T>::add(const SparseMatrixHYB<T> & matrix) const
    {
        if(this->m_ != matrix.m_ || this->n_ != matrix.n_)
        {
            throw std::invalid_argument("Cannot add: szie of matrices don't match.");
        }

        size_t m, n;
        m = this->m_;
        n = this->n_;

        SparseMatrixHYB<T> result(m, n);

        

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
    inline SparseMatrixHYB<T> SparseMatrixHYB<T>::subtract(const SparseMatrixHYB<T> & matrix) const
    {
        if(this->m_ != matrix.m_ || this->n_ != matrix.n_)
        {
            throw std::invalid_argument("Cannot add: szie of matrices don't match.");
        }

        size_t m, n;
        m = this->m_;
        n = this->n_;

        SparseMatrixHYB<T> result(m, n);

        

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
    inline std::vector<T> SparseMatrixHYB<T>::multiply(const std::vector<T> & x) const
    {
        if( this->n_ != x.size() )
        {
			throw std::invalid_argument("col of matrix != vector.size()");
        }


        std::vector<T> result(this->m_, T());
        T sum;
        size_t row_idx{0};
        

        for(size_t i=0; i<this->m_; ++i)
        {
            sum = T();
            row_idx = i * this->min_row_nnz_;
            for(size_t j=row_idx; j<row_idx+this->min_row_nnz_; ++j)
            {
                if(ell_col_indexes_[j] < 0)
                    break;

                sum += ell_values_[j] * x[ ell_col_indexes_[j] ];
            }

            for(size_t j=(*(this->ia_csr_))[i]; j<(*(this->ia_csr_))[i+1]; ++j  )
            {
                sum += (*(this->aa_csr_))[j] * x[ (*(this->ja_csr_))[j] ];
            }

            result[i] = sum;
        }


        return result;
    }

    template<typename T>
    inline SparseMatrixHYB<T> SparseMatrixHYB<T>::multiply(const SparseMatrixHYB<T> & matrix) const
    {
        if( this->n_ != matrix.m_ )
        {
			throw std::invalid_argument("Cannot multiply: col of left matrix != row of right matrix");
        }
        SparseMatrixHYB<T> result(this->m_, matrix.n_);
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

                
                result.set(sum, i, j);
            }
            
        }


        return result;
    }
    // ======================================================== //



    // ======================================================== // operations

    template<typename T>
    inline SparseMatrixHYB<T> & SparseMatrixHYB<T>::operator = (const SparseMatrixHYB<T> & matrix)
    {   

        if(&matrix != this)
        {
            this->copy(matrix);
        }
        return *this;
    }

    template<typename T>
	inline SparseMatrixHYB<T> SparseMatrixHYB<T>::operator + (const SparseMatrixHYB<T> & matrix) const
	{
		return this->add(matrix);
	}

    template<typename T>
    inline SparseMatrixHYB<T> SparseMatrixHYB<T>::operator - (const SparseMatrixHYB<T> & matrix) const
    {
        return this->subtract(matrix);
    }

    template<typename T>
    inline std::vector<T> SparseMatrixHYB<T>::operator * (const std::vector<T> & x) const
    {
        return this->multiply(x);
    }
    template<typename T>
    inline SparseMatrixHYB<T> SparseMatrixHYB<T>::operator * (const SparseMatrixHYB<T> & matrix) const
    {
        return this->multiply(matrix);
    }

   

    // ======================================================== //





    // ======================================================== //  friend founction

    template<typename T> 
	inline std::ostream & operator << (std::ostream & os, const SparseMatrixHYB<T> & matrix)
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

    template<typename T>
	inline bool operator == (const SparseMatrixHYB<T> & a, const SparseMatrixHYB<T> & b)
	{
        return  ( (a.vals == NULL && b.vals == NULL)
               || (a.vals != NULL && b.vals != NULL && *(a.vals) == *(b.vals)) 
                )

				&& 

                (  (a.cols == NULL && b.cols == NULL)
                || (a.cols != NULL && b.cols != NULL && *(a.cols) == *(b.cols)) 
                )

				&& 

                *(a.rows) == *(b.rows);
	}

    template<typename T>
	inline bool operator != (const SparseMatrixHYB<T> & a, const SparseMatrixHYB<T> & b)
    {
        return !( a == b );
    }


    // ======================================================== //


    // ======================================================== // Basic function for sparse matrix
    template<typename T>
    inline void SparseMatrixHYB<T>::construct(size_t rows, size_t cols, size_t min_row_nnz)
    {
        if (rows < 0 || cols < 0) 
        {
			throw "Matrix dimensions cannot be negative.";
		}

        this -> m_ = rows;
        this -> n_ = cols;
        this-> min_row_nnz_ = min_row_nnz;

        this->aa_csr_ = NULL;
		this->ja_csr_ = NULL;
		this->ia_csr_ = new std::vector<size_t>(rows+1, 0); 

        ell_col_indexes_.resize(min_row_nnz*rows, -1);
        ell_values_.resize(min_row_nnz*rows, 0);

        

    }

    template<typename T>
    inline void SparseMatrixHYB<T>::insert(size_t index, size_t row, size_t col, T val)
    {
        size_t row_idx{0};
        size_t flag{0};

        row_idx = row * this->min_row_nnz_;

        for(size_t i=row_idx; i<row_idx+this->min_row_nnz_; ++i)
        {
            if(ell_col_indexes_[i] < 0)
            {
                ell_col_indexes_[i] = col;
                ell_values_[i] = val;
                flag = 1;
                break;
            }
        }
        if(flag == 0)
        {
            if(this->aa_csr_ == NULL)
            {
            this->aa_csr_ = new std::vector<T>(1,val); 
            this->ja_csr_ = new std::vector<size_t>(1,col);
            }
            else
            {
                this->aa_csr_->insert(this->aa_csr_->begin()+index, val);
                this->ja_csr_->insert(this->ja_csr_->begin()+index, col); 
            }
            


            for (size_t i = row+1; i < this->m_+1; ++i) 
            {
                (*(this->ia_csr_))[i] += 1;
            }
        }
        
        

    }

    template<typename T>
	inline void SparseMatrixHYB<T>::remove(size_t index, size_t row, size_t col)
	{
        size_t row_idx{0};
        size_t flag{0};

        row_idx = row * this->min_row_nnz_;

        for(size_t i=row_idx; i<row_idx+this->min_row_nnz_; ++i)
        {
            if(ell_col_indexes_[i] == (int)col)
            {
                ell_col_indexes_[i] = -1;
                ell_values_[i] = T();
                flag = 1;
                break;
            }
        }

        if(flag == 0)
        {
            this->aa_csr_->erase(this->aa_csr_->begin() + index);
            this->ja_csr_->erase(this->ja_csr_->begin() + index);
            

            for (size_t i = row+1; i < this->m_+1; ++i) 
            {
                (*(this->ia_csr_))[i] -= 1;
            }
        }
		
	}

    template<typename T>
    inline void  SparseMatrixHYB<T>::copy(const SparseMatrixHYB<T> & matrix)
    {
        this->m_ = matrix.m_;
        this->n_ = matrix.n_;

        this->ia_csr_ = new std::vector<size_t>(*(matrix.ia_csr_));

        if(matrix.aa_csr_ != NULL)
        {
            this->ja_csr_ = new std::vector<size_t>(*(matrix.ja_csr_));
            this->aa_csr_ = new std::vector<T>(*(matrix.aa_csr_));
        }

        ell_col_indexes_.assign(matrix.ell_col_indexes_.begin(), matrix.ell_col_indexes_.end());
        ell_values_.assign(matrix.ell_values_.begin(), matrix.ell_values_.end());
        
        
    }

    template<typename T>
    inline void SparseMatrixHYB<T>::destruct(void)
    {
        if(this->aa_csr_ != NULL)
        {
            delete this->aa_csr_;
            delete this->ja_csr_;
        }

        
        delete this->ia_csr_;
    }

    template<typename T>
    inline void SparseMatrixHYB<T>::clc(void)
    {
        if(this->aa_csr_ != NULL)
            for(size_t pos=0; pos<this->nnz(); ++pos)
            {
                this->aa_csr_ = NULL;
                this->ja_csr_ = NULL;
            }

        for(size_t pos=0; pos<this->m_+1; ++pos)
        {
            (*(this->ia_csr_))[pos] = 0;
        }

        ell_values_.clear();
        ell_col_indexes_.clear();

    }


    template<typename T>
    inline size_t SparseMatrixHYB<T>::getCsrRowList(size_t pos) const
    {
        return (*(this->ia_csr_))[pos];
    }

    template<typename T>
    inline size_t SparseMatrixHYB<T>::getColumnList(size_t pos) const
    {
        return (*(this->ja_csr_))[pos];
    }

    template<typename T>
    inline T SparseMatrixHYB<T>::getValueList(size_t pos) const
    {
        return (*(this->aa_csr_))[pos];
    }

    // ======================================================== //


} // end namespace mat


#endif
