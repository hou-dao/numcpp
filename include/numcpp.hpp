#ifndef NUMCPP_H_
#define NUMCPP_H_

#include <ccomplex>
#include <cstdlib>
#include <cstring>

namespace numcpp {

    template <typename T>
    class col {
    public:        
 
        int ncol;
        T  *data;
 
        col (const int _ncol): ncol(_ncol) {
            data = new T[ncol];
        }
        col (const int _ncol, T& value): ncol(_ncol) {
            data = new T[ncol];
            for (int i=0; i<ncol; ++i) {
                data[i] = value;
            }
        }
        col (const col<T>& rhs): ncol(rhs.ncol) {
            data = new T[ncol];
            memcpy(data,rhs.data,ncol*sizeof(T));
        }
        ~col () {
            delete [] data;
            ncol = 0;
        }
        col<T>& operator=(col<T>& rhs) {
            if (this != &rhs) {
                if (data) {
                    delete [] data;                    
                }
                ncol = rhs.ncol;
                data = new T[ncol];
                memcpy(data,rhs.data,ncol*sizeof(T));        
            }
            return *this;
        }
        void zeros () {
            memset(data,0,ncol*sizeof(T));
        }
        inline T& operator()(int icol) const {
            return data[icol];
        }        
    };


    template <typename T>
    class mat {
    public:        
 
        int nrow, ncol, nmat;
        T  *data;
 
        mat (const int _nrow, const int _ncol): nrow(_nrow), ncol(_ncol), nmat(_nrow*_ncol) {
            data = new T[nmat];
        }
        mat (const int _nrow, const int _ncol, T& value): nrow(_nrow), ncol(_ncol), nmat(_nrow*_ncol) {
            data = new T[nmat];
            for (int i=0; i<nmat; ++i) {
                data[i] = value;
            }
        }
        mat (const mat<T>& rhs): nrow(rhs.nrow), ncol(rhs.ncol), nmat(rhs.nmat) {
            data = new T[nmat];
            memcpy(data,rhs.data,nmat*sizeof(T));
        }
        ~mat () {
            delete [] data;
            nrow = ncol = nmat = 0;
        }
        mat<T>& operator=(mat<T>& rhs) {
            if (this != &rhs) {
                if (data) {
                    delete [] data;                    
                }
                nrow = rhs.nrow;
                ncol = rhs.ncol;
                nmat = rhs.nmat;
                data = new T[nmat];
                memcpy(data,rhs.data,nmat*sizeof(T));        
            }
            return *this;
        }
        void zeros () {
            memset(data,0,nmat*sizeof(T));
        }
        inline T& operator()(int irow, int icol) const {
            return data[irow*ncol+icol];
        }        
    };
    
    template <typename T>
    class cub {
    public:        
 
        int nslc, nrow, ncol, nmat, ncub;
        T  *data;
 
        cub (const int _nslc, const int _nrow, const int _ncol): nslc(_nslc), nrow(_nrow), ncol(_ncol), nmat(_nrow*_ncol), ncub(_nslc*_nrow*_ncol)  {
            data = new T[ncub];
        }
        cub (const int _nslc, const int _nrow, const int _ncol, T& value): nslc(_nslc), nrow(_nrow), ncol(_ncol), nmat(_nrow*_ncol), ncub(_nslc*_nrow*_ncol)  {
            data = new T[ncub];
            for (int i=0; i<ncub; ++i) {
                data[i] = value;
            }
        }
        cub (const cub<T>& rhs): nslc(_nslc), nrow(_nrow), ncol(_ncol), nmat(_nrow*_ncol), ncub(_nslc*_nrow*_ncol)  {
            data = new T[ncub];
            memcpy(data,rhs.data,ncub*sizeof(T));
        }
        ~cub () {
            delete [] data;
            nslc = nrow = ncol = nmat = ncub = 0;
        }
        cub<T>& operator=(cub<T>& rhs) {
            if (this != &rhs) {
                if (data) {
                    delete [] data;                    
                }
                nslc = rhs.nslc;
                nrow = rhs.nrow;
                ncol = rhs.ncol;
                nmat = rhs.nmat;
                ncub = rhs.ncub;
                data = new T[ncub];
                memcpy(data,rhs.data,ncub*sizeof(T));        
            }
            return *this;
        }
        void zeros () {
            memset(data,0,ncub*sizeof(T));
        }        
        inline T& operator()(int islc, int irow, int icol) const {
            return data[islc*nmat+irow*ncol+icol];
        }        
    };
}

#endif
