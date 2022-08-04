
#ifndef __CSRMATRIX_H__
#define	__CSRMATRIX_H__

    #include <algorithm>
    #include <vector>
    #include <iostream>
    #include <stdexcept>

    using std::vector;
    using std::cout;
    using std::endl;


namespace mat
{
    
    template<typename T>
    class csr_matrix
    {
        private:
        

        public:
        //--------------------------------------------- // init sparse matrix
            csr_matrix(int n); // square sparse matrix
            csr_matrix(int rows, int columns); //sparse matrix
            csr_matrix(const csr_matrix<T> & m); // copy function
            csr_matrix & set(T val, int row, int col); // add non-zero value to sparse matrix
            csr_matrix(void);  // delete memeory
        //---------------------------------------------


        //--------------------------------------------- // function of sparse matrix
            int getRowCount(void) const;  
            int getColumnCount(void) const;
            int getNnzCount(void) const;
            int getCsrRowList(int i) const;
            int getColumnList(int i) const;
            T getValueList(int i) const;
            T getValue(int row, int col) const;  // CSR core; input row, col. u can get the value in sparse matrix as well as u want
            csr_matrix<T> add(const csr_matrix<T> & matrix) const; // mat1 + mat2
            csr_matrix<T> subtract(const csr_matrix<T> & matrix) const; // mat1 - mat2
            vector<T> multiply(const vector<T> & x) const; // mat * vector 
            csr_matrix<T> multiply(const csr_matrix<T> & matrix) const; // mat1 * mat2 
        //---------------------------------------------

        //--------------------------------------------- // operations
            csr_matrix<T> & operator = (const csr_matrix<T> & matrix);
            csr_matrix<T> operator + (const csr_matrix<T> & matrix) const;
            csr_matrix<T> operator - (const csr_matrix<T> & matrix) const;
            vector<T> operator * (const vector<T> & x) const;
            csr_matrix<T> operator * (const csr_matrix<T> & matrix) const;
        //---------------------------------------------

        //--------------------------------------------- // friend founction
            template<typename X>
				friend std::ostream & operator << (std::ostream & os, const csr_matrix<X> & matrix);

            template<typename X>
				friend bool operator == (const csr_matrix<X> & a, const csr_matrix<X> & b);

            template<typename X>
				friend bool operator != (const csr_matrix<X> & a, const csr_matrix<X> & b);
        //---------------------------------------------




        private:

        //---------------------------------------------// variables
            int m, n; // m:row, n:col
            vector<T> *vals; // vals list
            vector<int> *rows, *cols;  //csr_rows cols list 
        //---------------------------------------------

        //--------------------------------------------- // Basic function for csr_matrix class
            void construct(int rows, int column);
            void insert(int index, int row, int col, T val);
            void remove(int index, int row);
            void copy(const csr_matrix<T> & m);
            void destruct(void);
        //---------------------------------------------

    };






    // ======================================================== // init sparse matrix
    template<typename T>
    inline csr_matrix<T>::csr_matrix(int n)
    {
        this->construct(n, n);
    }

    template<typename T>
	inline csr_matrix<T>::csr_matrix(int rows, int columns)
	{
		this->construct(rows, columns);
	}



    te2_Examplelate<typename T>
	inline csr_matrix<T> & csr_matrix<T>::set(T val, int row, int col)
	{

        if (row < 0 || col < 0 || row > this->m-1 || col > this->n-1) 
        {
			throw std::invalid_argument("set: Coordination out of range.");
		}

       
		int pos = 0;
		int currCol = 0;
        
		for (pos = (*(this->rows))[row]; pos < (*(this->rows))[row+1] ; ++pos) 
        {
			currCol = (*(this->cols))[pos];

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
                (*(this->vals))[pos] = val; 
            }
        }
        else if( !(this->getValue(row, col)==T()) )
        {
            this->remove(pos, row);
        }
        
		return *this;
	}

    



    template<typename T>
	inline csr_matrix<T>::csr_matrix(const csr_matrix<T> & matrix)
	{
		this->copy(matrix);
	}

    template<typename T>
    inline csr_matrix<T>::~csr_matrix(void)
    {
        this->destruct();
    }

    // ======================================================== //

    
    
    



    // ======================================================== // function of sparse matrix

    template<typename T>
    inline csr_matrix<T> csr_matrix<T>::add(const csr_matrix<T> & matrix) const
    {
        if(this->m != matrix.m || this->n != matrix.n)
        {
            throw std::invalid_argument("Cannot add: szie of matrices don't match.");
        }

        int m, n;
        m = this->m;
        n = this->n;

        csr_matrix<T> result(m, n);

        

        for(int i=0; i<this->m; ++i)
        {
            for(int j=0; j<this->n; ++j)
            {
               if(this->getValue(i,j)!=T() ||  matrix.getValue(i,j)!=T())
                    result.set( this->getValue(i,j) + matrix.getValue(i,j), i, j );


            }
        }
        
        return result;
    }

    template<typename T>
    inline csr_matrix<T> csr_matrix<T>::subtract(const csr_matrix<T> & matrix) const
    {
        if(this->m != matrix.m || this->n != matrix.n)
        {
            throw std::invalid_argument("Cannot add: szie of matrices don't match.");
        }

        int m, n;
        m = this->m;
        n = this->n;

        csr_matrix<T> result(m, n);


        for(int i=0; i<this->m; ++i)
        {
            for(int j=0; j<this->n; ++j)
            {
               if(this->getValue(i,j)!=T() ||  matrix.getValue(i,j)!=T())
                    result.set( this->getValue(i,j) - matrix.getValue(i,j), i, j );


            }
        }
        
        return result;
    }

    template<typename T> 
    inline vector<T> csr_matrix<T>::multiply(const vector<T> & x) const
    {
        if( this->n != (int)x.size() )
        {
			throw std::invalid_argument("col of matrix != vector.size()");
        }

        vector<T> result(this->m, T());
        T sum;

        for(int i=0; i<this->m; ++i)
        {
            sum = T();
            for(  int j=(*(this->rows))[i]; j<(*(this->rows))[i+1]; ++j  )
            {
                sum += (*(this->vals))[j] * x[ (*(this->cols))[j] ];
            }
            result[i] = sum;
        }
        return result;
    }

    template<typename T>
    inline csr_matrix<T> csr_matrix<T>::multiply(const csr_matrix<T> & matrix) const
    {
        if( this->n != matrix.m )
        {
			throw std::invalid_argument("Cannot multiply: col of left matrix != row of right matrix");
        }
        csr_matrix<T> result(this->m, matrix.n);
        T sum;

        for(int i=0; i<this->m; ++i)
        {
            for(int j=0; j<matrix.n; ++j)
            {
                sum = T();
                for(int k=0; k<this->n; ++k)
                {
                    sum += this->getValue(i,k) * matrix.getValue(k,j);
                }

                
                result.set(sum, i, j);
            }
            
        }


        return result;
    }



    template<typename T>
    inline int csr_matrix<T>::getRowCount(void) const
    {
        return this->m;
    }

    template<typename T>
    inline int csr_matrix<T>::getColumnCount(void) const
    {
        return this->n;
    }

    template<typename T>
    inline int csr_matrix<T>::getNnzCount(void) const
    {
        return (*(this->rows))[this->getRowCount()];
    }

    template<typename T>
    inline int csr_matrix<T>::getCsrRowList(int i) const
    {
        return (*(this->rows))[i];
    }

    template<typename T>
    inline int csr_matrix<T>::getColumnList(int i) const
    {
        return (*(this->cols))[i];
    }

    template<typename T>
    inline T csr_matrix<T>::getValueList(int i) const
    {
        return (*(this->vals))[i];
    }


   
    template<typename T>
	inline T csr_matrix<T>::getValue(int row, int col) const
	{
        if (row < 0 || col < 0 || row > this->m-1 || col > this->n-1) 
        {
			throw std::invalid_argument("getvalue: Coordination out of range.");
		}

		int currCol = 0;
        // example   input=1,1
        // rows  = {0,         3,6  ,9}. 
        // cols = {0,1,2,     0,1,3  ,0,2,3,1,2,3}
        // vals = {4,-1,-1,  -1,4,-1  ,-1,4,-1,-1,-1,4}
		for (int pos = (*(this->rows))[row]; pos < (*(this->rows))[row+1]; ++pos)  //pos=3,4,5
        {
			currCol = (*(this->cols))[pos];  

			if (currCol == col) //when pos = 4 currCol = cols[4] = col = 1. 4 means number 4 of nnz value
            {
				return (*(this->vals))[pos]; 

			} 
            else if (currCol > col) 
            {
				break;
			}
		}

		return T();
	}

    // ======================================================== //



    // ======================================================== // operations

    template<typename T>
    inline csr_matrix<T> & csr_matrix<T>::operator = (const csr_matrix<T> & matrix)
    {   
        if(&matrix != this)
        {
            this->destruct();
            this->copy(matrix);
        }
        return *this;
    }

    template<typename T>
	inline csr_matrix<T> csr_matrix<T>::operator + (const csr_matrix<T> & matrix) const
	{
		return this->add(matrix);
	}

    template<typename T>
    inline csr_matrix<T> csr_matrix<T>::operator - (const csr_matrix<T> & matrix) const
    {
        return this->subtract(matrix);
    }

    template<typename T>
    inline vector<T> csr_matrix<T>::operator * (const vector<T> & x) const
    {
        return this->multiply(x);
    }
    template<typename T>
    inline csr_matrix<T> csr_matrix<T>::operator * (const csr_matrix<T> & matrix) const
    {
        return this->multiply(matrix);
    }
    // ======================================================== //





    // ======================================================== //  friend founction

    template<typename T> 
	inline std::ostream & operator << (std::ostream & os, const csr_matrix<T> & matrix)
	{
        os << std::endl;
        
		for (int i = 0; i < matrix.m; ++i) 
        {
			for (int j = 0; j < matrix.n; ++j) 
            {
				if (j != 0) 
                {
					os << " ";
				}

				os << matrix.getValue(i, j);
			}

			
            os << std::endl;
		}

		return os;
	}

    template<typename T>
	inline bool operator == (const csr_matrix<T> & a, const csr_matrix<T> & b)
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
	inline bool operator != (const csr_matrix<T> & a, const csr_matrix<T> & b)
    {
        return !( a == b );
    }


    // ======================================================== //


    // ======================================================== // Basic function for sparse matrix
    template<typename T>
    inline void csr_matrix<T>::construct(int rows, int columns)
    {
        if (rows < 0 || columns < 0) 
        {
			throw "Matrix dimensions cannot be negative.";
		}

        this -> m = rows;
        this -> n = columns;

        this->vals = NULL;
		this->cols = NULL;
		this->rows = new vector<int>(rows+1, 0); 


    }

    template<typename T>
    inline void csr_matrix<T>::insert(int index, int row, int col, T val)
    {
        if(this -> vals == NULL)
        {
            this->vals = new vector<T>(1, val);
            this->cols = new vector<int>(1,col);
        }
        else
        {
            this->vals->insert(this->vals->begin()+index, val);
            this->cols->insert(this->cols->begin()+index, col); 
        }

        for (int i = row+1; i < this->m+1; ++i) 
        {
			(*(this->rows))[i] += 1;
		}
        

    }

    template<typename T>
	inline void csr_matrix<T>::remove(int index, int row)
	{
		this->vals->erase(this->vals->begin() + index);
		this->cols->erase(this->cols->begin() + index);


		for (int i = row+1; i < this->m+1; ++i) 
        {
			(*(this->rows))[i] -= 1;
		}
	}

    template<typename T>
    inline void  csr_matrix<T>::copy(const csr_matrix<T> & matrix)
    {
        this->m = matrix.m;
        this->n = matrix.n;

        this->rows = new vector<int>(*(matrix.rows));

        if(matrix.vals != NULL)
        {
            this->cols = new vector<int>(*(matrix.cols));
            this->vals = new vector<T>(*(matrix.vals));
        }
        
        
    }

    template<typename T>
    inline void csr_matrix<T>::destruct(void)
    {
        
        if (this->vals != NULL) 
        {
			delete this->vals;
			delete this->cols;
	    }

		delete this->rows;
        
    }
    // ======================================================== //


} // end namespace mat


#endif