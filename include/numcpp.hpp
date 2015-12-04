#ifndef NUMCPP_H_
#define NUMCPP_H_

#include <ccomplex>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace numcpp {
    
    static const char _A_ = 'A';
    static const char _N_ = 'N';
    static const char _U_ = 'U';
    static const char _L_ = 'L';

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
        void copy (col<T>& rhs) {
            memcpy(data,rhs.data,ncol*sizeof(T));
        }
        // IO
        void load (const char filename[], const char type[]) {
            FILE *fp = fopen(filename,"r");
            fclose(fp);
        }
        // get Element
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
        void copy (mat<T>& rhs) {
            memcpy(data,rhs.data,nmat*sizeof(T));
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
        void copy (cub<T>& rhs) {
            memcpy(data,rhs.data,ncub*sizeof(T));
        }               
        inline T& operator()(int islc, int irow, int icol) const {
            return data[islc*nmat+irow*ncol+icol];
        }        
    };
    


    void evalst (col<double>& D, col<double>& E) {
        if (D.ncol<1 || E.ncol<0) {
            printf ("No element in D or E vectors?\n");
        } else {
            int INFO = 0;
            double *Z = 0;
            double *WORK = new double[2*D.ncol-2];
            dstev_ ((char *)&_N_,&D.ncol,D.data,E.data,Z,&D.ncol,WORK,&INFO);
            delete [] WORK;
        }
    }

    void esysst  (col<double>& D, col<double>& E, mat<double>& evec) {
        if (D.ncol<1 || E.ncol<0) {
            printf ("No element in D or E vectors?\n");
        } else {
            int INFO = 0;
            double *Z = 0;
            double *WORK = new double[2*D.ncol-2];
            dstev_ ((char *)&_V_,&D.ncol,D.data,E.data,evec.data,&D.ncol,WORK,&INFO);
            delete [] WORK;
        }
    }

    void evalsh (const mat<double>& A, col<double>& eval) {
        if (A.nrow!=A.ncol || A.nmat<=0) {
            printf ("No element in A matrix?\n");
        } else {
            int INFO = 0;
            int LWORK = 3*A.nrow-1;
            double *VECA = new double[A.nmat];
            double *WORK = new double[LWORK];
            memcpy(VECA,A.data,A.nmat*sizeof(double));
            dsyev_ ((char *)&_N_,(char *)&_L_,&A.nrow,VECA,&A.nrow,eval.data,WORK,&LWORK,&INFO);
            delete [] WORK;
            delete [] VECA;
        }
    }

    void evalsh (const mat<dcmplx>& A, col<double>& eval) {
        if (A.nrow!=A.ncol || A.nmat<=0 || A.nrow!=eval.ncol) {
            printf ("No element in A matrix?\n");
        } else {
            int INFO = 0;
            int LWORK = 2*A.nrow-1;
            dcmplx *VECA = new dcmplx[A.nmat];
            dcmplx *WORK = new dcmplx[LWORK];
            double *RWORK = new double[3*A.nrow-2];
            memcpy(VECA,A.data,A.nmat*sizeof(dcmplx));
            zheev_ ((char *)&_N_,(char *)&_L_,&A.nrow,VECA,&A.nrow,eval.data,WORK,&LWORK,RWORK,&INFO);
            delete [] RWORK;
            delete [] WORK;
            delete [] VECA;
        }
    }
    
    void esyssh (const mat<double>& A, col<double>& eval, mat<double>& evec) {
        if (A.nrow!=A.ncol || A.nmat<=0) {
            printf ("No element in A matrix?\n");
        } else {
            int INFO = 0;
            int LWORK = 3*A.nrow-1;
            double *WORK = new double[LWORK];
            memcpy(evec.data,A.data,A.nmat*sizeof(double));
            dsyev_ ((char *)&_N_,(char *)&_L_,&A.nrow,evec.data,&A.nrow,eval.data,WORK,&LWORK,&INFO);
            delete [] WORK;
        }
    }

    void esyssh (const mat<dcmplx>& A, col<double>& eval, mat<dcmplx>& evec) {
        if (A.nrow!=A.ncol || A.nmat<=0 || A.nrow!=eval.ncol) {
            printf ("No element in A matrix?\n");
        } else {
            int INFO = 0;
            int LWORK = 2*A.nrow-1;
            dcmplx *WORK = new dcmplx[LWORK];
            double *RWORK = new double[3*A.nrow-2];
            memcpy(evec.data,A.data,A.nmat*sizeof(dcmplx));       
            zheev_ ((char *)&_N_,(char *)&_L_,&A.nrow,evec.data,&A.nrow,eval.data,WORK,&LWORK,RWORK,&INFO);
            for (int i=0; i<n-1; ++i) {
                for (int j=i+1; j<n; ++j) {
                    dcmplx tmp = evec(i,j);
                    evec(i,j) = conj(evec(j,i));
                    evec(j,i) = conj(tmp);
                }
            }
            delete [] RWORK;
            delete [] WORK;
        }
    }



/*
    SEigensystem.c
    diagonalization of a complex symmetric n-by-n matrix using
    the Jacobi algorithm
    code adapted from the "Handbook" routines for complex A
    (Wilkinson, Reinsch: Handbook for Automatic Computation, p. 202)
    this file is part of the Diag library
    last modified 1 Aug 11 th
 *
    SEigensystem diagonalizes a complex symmetric n-by-n matrix.
    Input: n, A = n-by-n matrix, complex symmetric
    (only the upper triangle of A needs to be filled).
    Output: d = vector of eigenvalues, U = transformation matrix
    these fulfill diag(d) = U A U^T = U A U^-1 with U U^T = 1.
 */

    void esyszs (const mat<dcmplx>& A, col<dcmplx>& d, mat<dcmplx>& U, const int sort) {
        
        int n = A.nrow;
        double red = 0.04/(n*n*n*n);
        
        auto Sq = [](dcmplx x) {return creal(x*conj(x));}

        mat<dcmplx> &ev = mat<dcmplx>(n,2);

        for (int p=0; p<n; ++p) {
            ev(p,0) = 0;
            ev(p,1) = d(p) = A(p,p);
        }

        U.zeros();
        for (int p=0; p<n; ++p) {
            U(p,p) = 1;
        }

        for (int nsweeps=1; nsweeps<=50; ++nsweeps) {
            
            double thresh = 0;
            
            for (int q=1; q<n; ++q) {
                for (int p=0; p<q; ++p) {
                    thresh += creal(A(p,q)*conj(A(p,q)));
                }
            }
            if ( !(thresh>SYM_EPS) ) goto done;

            thresh = (nsweeps<4)?thresh*red:0;

            for (int q=1; q<n; ++q) {
                for (int p=0; p<q; ++p) {
                    dcmplx delta = A(p,q);
                    double off = creal(delta*conj(delta));
                    double sqp = creal(ev(p,1)*conj(ev(p,1)));
                    double sqq = creal(ev(q,1)*conj(ev(q,1)));
                    if (nsweeps>4 && off<SYM_EPS*(sqp+sqq)) {
                        A(p,q) = 0;
                    } else if (off > thresh) {

                        dcmplx x = 0.5*(ev(p,1)-ev(q,1));
                        dcmplx y = csqrt(x*x+delta*delta);
                        dcmplx t = x-y;
                        dcmplx s = x+y;
                        if (Sq(t) < Sq(s)) {
                            t = s;
                        }

                        t = delta/t;
                        delta *= t;
                        ev(p,0) += delta;
                        ev(q,0) -= delta;
                        ev(p,1) = d(p)+ev(p,0);
                        ev(q,1) = d(q)+ev(q,0);

                        dcmplx invc = csqrt(t*t+1);
                        s = t/invc;
                        t /= invc+1;

                        for (int j=0; j<p; ++j) {
                            dcmplx x = A(j,p);
                            dcmplx y = A(j,q);
                            A(j,p) = x+s*(y-t*x);
                            A(j,q) = y-s*(x+t*y);
                        }

                        for (int j=p+1; j<q; ++j) {
                            dcmplx x = A(p,j);
                            dcmplx y = A(j,q);
                            A(p,j) = x+s*(y-t*x);
                            A(j,q) = y-s*(x+t*y);
                        }

                        for (int j=q+1; j<n; ++j) {
                            dcmplx x = A(p,j);
                            dcmplx y = A(q,j);
                            A(p,j) = x+s*(y-t*x);
                            A(q,j) = y-s*(x+t*y);
                        }

                        A(p,q) = 0;

                        for (int j=0; j<n; ++j) {
                            dcmplx x = U(p,j);
                            dcmplx y = U(q,j);
                            U(p,j) = x+s*(y-t*x);
                            U(q,j) = y-s*(x+t*y);
                        }
                    }
                }
            }

            for (int p=0; p<n; ++p) {
                ev(p,0) = 0;
                d(p) = ev(p,1);
            }
        }

        printf("Bad convergence in SEigensystem\n");

    done:

        if (sort == 0) return;

        /* sort the eigenvalues by their real part */

        for (int p=0; p<n-1; ++p) {
            int j = p;
            dcmplx t = d(p);
            for (int q=p+1; q<n; ++q) {
                if (sort*(creal(t)-creal(d(q)) > 0) {
                    j = q;
                    t = d(q);
                }
            }
            if (j == p) continue;
            d(j) = d(p);
            d(p) = t;
            for (int q = 0; q<n; ++q) {
                dcmplx x = U(p,q);
                U(p,q) = U(j,q);
                U(j,q) = x;
            }
        }
    }
}

#endif
