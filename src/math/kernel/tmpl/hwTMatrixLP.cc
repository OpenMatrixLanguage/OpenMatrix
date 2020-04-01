/**
* @file hwTMatrixLP.cc
* @date November 2012
* Copyright (C) 2012-2018 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language ("OpenMatrix") software.
* Open Source License Information:
* OpenMatrix is free software. You can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* OpenMatrix is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
* You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
* 
* Commercial License Information: 
* For a copy of the commercial license terms and conditions, contact the Altair Legal Department at Legal@altair.com and in the subject line, use the following wording: Request for Commercial License Terms for OpenMatrix.
* Altair's dual-license business model allows companies, individuals, and organizations to create proprietary derivative works of OpenMatrix and distribute them - whether embedded or bundled with other software - under a commercial license agreement.
* Use of Altair's trademarks and logos is subject to Altair's trademark licensing policies.  To request a copy, email Legal@altair.com and in the subject line, enter: Request copy of trademark and logo usage policy.
*/

//:---------------------------------------------------------------------------
//:Description
//
//  hwTMatrix explicit <double> specialization implementation file
//  with BLAS / LAPACK functions 
//
//:---------------------------------------------------------------------------

#include <hwComplex.h>
#include <complex>
typedef std::complex<double> complexD;

//*******************************************************************
//                    BLAS / LAPACK prototypes
//*******************************************************************
#ifndef _BLAS_LAPACK_h
#define _BLAS_LAPACK_h

// y = x
extern "C" void dcopy_(int* N, double* DX, int* INCX, double* DY, int* INCY);
extern "C" void zcopy_(int* N, complexD* DX, int* INCX, complexD* DY, int* INCY);
// y = ay
extern "C" void dscal_(int* N, double* DA, double* DX, int* INCX);
extern "C" void zscal_(int* N, complexD* DA, complexD* DX, int* INCX);
// y = ax + y
extern "C" void daxpy_(int* N, double* DA, double* DX, int* INCX, double* DY, int* INCY);
extern "C" void zaxpy_(int* N, complexD* DA, complexD* DX, int* INCX, complexD* DY, int* INCY);
// C = alpha*op( A )*op( B ) + beta*C
extern "C" void dgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K,
                       double* ALPHA, double* A, int* LDA, double* B, int* LDB,
                       double* BETA, double* C, int* LDC);
extern "C" void zgemm_(char* TRANSA, char* TRANSB, int* M, int* N, int* K,
                       complexD* ALPHA, complexD* A, int* LDA, complexD* B,
                       int* LDB, complexD* BETA, complexD* C, int* LDC);
// Dot(X,Y)
extern "C" double ddot_(int* N, double* DX, int* INCX, double* DY, int* INCY);
// There are two conventions for how to return a complex value in FORTRAN.
// 1. The standard BLAS convention in Windows, which simply returns the value.
// 2. The standard BLAS convention in LINUX, which passes a pointer to the
//    return value as the first argument.
// MKL follows the second convention for both Windows and Linux.
#ifdef STANDARD_BLAS    // to be defined in OpenMatrix
  #if defined OS_WIN
    extern "C" complexD zdotc_(int* N, complexD* DX, int* INCX, complexD* DY, int* INCY);
  #else
    extern "C" void zdotc_(complexD* dotc, int* N, complexD* DX, int* INCX, complexD* DY, int* INCY);
  #endif
#else   // MKL prototype
    extern "C" void zdotc_(complexD* dotc, int* N, complexD* DX, int* INCX, complexD* DY, int* INCY);
#endif
// LU decomposition
extern "C" void dgetrf_(int* M, int* N, double* A, int* LDA, int* IPIV, int* INFO);
extern "C" void zgetrf_(int* M, int* N, complexD* A, int* LDA, int* IPIV, int* INFO);
// Matrix inversion
extern "C" void dgetri_(int* N, double* A, int* LDA, int* IPIV, double* WORK, int* LWORK, int* INFO);
extern "C" void zgetri_(int* N, complexD* A, int* LDA, int* IPIV, complexD* WORK, int* LWORK, int* INFO);
// Solve AX=B via LU
extern "C" void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV,
                       double* B, int* LDB, int* INFO);
extern "C" void zgesv_(int* N, int* NRHS, complexD* A, int* LDA, int* IPIV,
                       complexD* B, int* LDB, int* INFO);
// Solve AX=B via LU for triangular A
extern "C" void dgtsv_(int* N, int* NRHS, double* DL, double* D,
                       double* DU, double* B, int* LDB, int* INFO);
// Solve AX=B via LU for band A
extern "C" void dgbsv_(int* N, int* KL, int* KU, int* NRHS, double* AB,
                       int* LDAB, int* IPIV, double* B, int* LDB, int* INFO);
// Solve AX=B via Cholesky for SPD A
extern "C" void dposv_(char* UPLO, int* N, int* NRHS, double* A, int* LDA, double* B,
                       int* LDB, int* INFO);
// Solve AX=B via Cholesky for SPD triangular A
extern "C" void dptsv_(int* N, int* NRHS, double* D, double* E, double* B,
                       int* LDB, int* INFO);
// Solve AX=B via Cholesky for SPD BAND A
extern "C" void dpbsv_(char* UPLO, int* N, int* KD, int* NRHS, double* AB,
                       int* LDAB, double* B, int* LDB, int* INFO);
// Solve AX=B for symmetric indefinite A
extern "C" void dsysv_(char* UPLO, int* N, int* NRHS, double* A, int* LDA,
                       int* IPIV, double* B, int* LDB, double* WORK,
                       int* LWORK, int* INFO);
// Eigen decomposition
extern "C" void dsyev_(char* JOBZ, char* UPLO, int* N, double* A, int* LDA, double* W,
                       double* WORK, int* LWORK, int* INFO);
extern "C" void dgeevx_(char* balanc, char* jobvl, char* jobvr, char* sense, int* n,
                        double* a, int* lda, double* wr ,double* wi, double* vl,
                        int* ldvl, double* vr, int* ldvr, int* ilo, int* ihi,
                        double* scale, double* abnrm, double* rconde, double* rcondv,
                        double* work, int* lwork, int* iwork, int* info);
extern "C" void zheev_(char* jobz, char* uplo, int* n, complexD* a, int* lda, double* w,
                       complexD* work, int* lwork, double* rwork, int* info);
extern "C" void zgeevx_(char* balanc, char* jobvl, char* jobvr, char* sense, int* n,
                        complexD* a, int* lda, complexD* w, complexD* vl, int* ldvl,
                        complexD* vr, int* ldvr, int* ilo, int* ihi, double* scale,
                        double* abnrm, double* rconde, double* rcondv, complexD* work,
                        int* lwork, double* rwork, int* info);
extern "C" void dstev_(char* JOBZ, int* N, double* D, double* E, double* Z,
                       int* LDZ, double* WORK, int* INFO);
extern "C" void dsbev_(char* JOBZ, char* UPLO, int* N, int* KD, double* AB,
                       int* LDAB, double* W, double* Z, int* LDZ,
                       double* WORK, int* INFO);
//! Balance a matrix
extern "C" void dgebal_(char* job, int* n, double* a, int* lda, int* ilo, int* ihi,
                        double* scale, int* info);
extern "C" void zgebal_(char* job, int* n, complexD* a, int* lda, int* ilo, int* ihi,
                        double* scale, int* info);
//! Generalized Eigen decomposition
extern "C" void dsygv_(int* itype, char* jobz, char* uplo ,int* n, double* a, int* lda,
                       double* b, int* ldb, double* w, double* work, int* lwork, int* info);
extern "C" void dggev_(char* jobvl, char* jobvr ,int* n, double* a, int* lda, double* b,
                       int* ldb, double* alphar, double* alphai, double* beta, double* vl,
                       int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info);
extern "C" void zhegv_(int* itype, char* jobz, char* uplo, int* n, complexD* a, int* lda,
                       complexD* b, int* ldb, double* w, complexD* work, int* lwork,
                       double* rwork, int* info);
extern "C" void zggev_(char* jobvl, char* jobvr, int* n, complexD* a, int* lda, complexD* b,
                       int* ldb, complexD* alpha, complexD* beta, complexD* vl, int* ldvl,
                       complexD* vr, int* ldvr, complexD* work, int* lwork, double* rwork,
                       int* info);
// Cholesky decomposition of SPD matrix
extern "C" void dpotrf_(char* UPLO, int* N, double* A, int* LDA, int* INFO);
// Singular value decomposition
extern "C" void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A, int* LDA, double* S, double* U,
                        int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* INFO);
extern "C" void zgesvd_(char *jobu, char *jobvt, int *m, int *n, complexD *a, int *lda, double *s, complexD *u,
                        int *ldu, complexD *vt, int *ldvt, complexD *work, int *lwork, double *rwork, int *info);
extern "C" void dgesdd_(char* JOBZ, int* M, int* N, double* A, int* LDA, double* S, double* U,
                        int* LDU, double* VT, int* LDVT, double* WORK, int* LWORK, int* INFO);
// QR decomposition
extern "C" void dgeqrf_(int* M, int* N, double* A, int* LDA, double* TAU, double* WORK,
                        int* LWORK, int* INFO);
extern "C" void dorgqr_(int* M, int* N, int* K, double* A, int* LDA, double* TAU,
                        double* WORK, int* LWORK, int* INFO);
extern "C" void zgeqrf_(int* M, int* N, complexD* A, int* LDA, complexD* TAU, complexD* WORK,
                        int* LWORK, int* INFO);
extern "C" void zungqr_(int* M, int* N, int* K, complexD* A, int* LDA, complexD* TAU,
                        complexD* WORK, int* LWORK, int* INFO);
// Schur decomposition
extern "C" void dgees_(char* jobvs, char* sort, void* select, int* n, double* a, int* lda, int* sdim,
                       double* wr, double* wi, double* vs, int* ldvs, double* work ,int* lwork,
                       void* bwork, int* info);
extern "C" void zgees_(char* jobvs, char* sort, void* select, int* n, complexD* a, int* lda,
                       int* sdim, complexD* w, complexD* vs, int* ldvs, complexD* work,
                       int* lwork, double* rwork, void* bwork, int* info);
// Solve AX=B via QR
//extern "C" void dgels_(char* TRANS, int* M, int* N, int* NRHS, double* A, int* LDA, double* B, int* LDB,
//                       double* WORK, int* LWORK, int* INFO);
//extern "C" void zgels_(char* TRANS, int* M, int* N, int* NRHS, complexD* A, int* LDA, complexD* B, int* LDB,
//                       complexD* WORK, int* LWORK, int* INFO);
extern "C" void dgelsy_(int* M, int* N, int* NRHS, double* A, int* LDA, double* B, int* LDB, int* JPVT,
                        double* RCOND, int* RANK, double* WORK, int* LWORK, int* INFO);
extern "C" void zgelsy_(int* M, int* N, int* NRHS, complexD* A, int* LDA, complexD* B, int* LDB, int* JPVT,
                        double* RCOND, int* RANK, complexD* WORK, int* LWORK, double* RWORK, int* INFO);
// Solve AX=B via SVD
extern "C" void dgelss_(int* M, int* N, int* NRHS, double* A, int* LDA, double* B, int* LDB, double* S,
                        double* RCOND, int* RANK, double* WORK, int* LWORK, int* INFO);
extern "C" void zgelss_(int* M, int* N, int* NRHS, complexD* A, int* LDA, complexD* B, int* LDB, double* S,
                        double* RCOND, int* RANK, complexD* WORK, int* LWORK, double* RWORK, int* INFO);
// Condition number
extern "C" void dgecon_(char* NORM, int* N, double* A, int* LDA, double* ANORM, double* RCOND,
                        double* WORK, int* IWORK, int* INFO);
extern "C" void zgecon_(char* norm, int* n, complexD* a, int* lda, double* anorm, double* rcond,
                        complexD* work, double* rwork, int* info);
// Matrix norm
extern "C" double dlange_(char* norm, int* m, int* n, double* a, int* lda, double* work);
extern "C" double zlange_(char* norm, int* m, int* n, complexD* a, int* lda, double* work);
// Vector L2 norm
extern "C" double dnrm2_(int* N,double* X, int* INCX);

#endif // _BLAS_LAPACK_h

//*******************************************************************
//           hwTMatrix<double> private implementations
//*******************************************************************
// These functions are implemented before the public functions so
// that they are specialized prior to instantiation

//! Copy data
template<>
inline void hwTMatrix<double>::CopyData(void* dest, int arraySize, const void* src, int count)
{
    int inc = 1;

    if (IsReal())
        dcopy_((int*) &count, (double*) src, &inc, (double*) dest, &inc);
    else
        zcopy_((int*) &count, (complexD*) src, &inc, (complexD*) dest, &inc);
}

//! Real LU decomposition (PA = LU)
template<>
inline hwMathStatus hwTMatrix<double>::RealLU(hwTMatrix<double>& L, hwTMatrix<double>& U, hwTMatrix<int>& P ) const
{
    if (this == &L)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &U)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;
    int m = A.m_nRows;
    int n = A.m_nCols;
    int nu = _min(m, n);

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    hwTMatrix<double> B(A);

    // dimension output matrices
    status = L.Dimension(m, nu, REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(3);
        else
            status.ResetArgs();

        return status;
    }

    status = U.Dimension(nu, n, REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(4);
        else
            status.ResetArgs();

        return status;
    }

    status = P.Dimension(m, hwTMatrix<int>::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(2);
        else
            status.ResetArgs();

        return status;
    }

    if (A.Size() != 0)
    {
        // decompose the matrix
        double* a = B.m_real;
        int lda = m;
        int info;
        int* ipiv;

        try
        {
            ipiv = new int[nu];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        dgetrf_(&m, &n, a, &lda, ipiv, &info); 

        if (info != 0)
            status(HW_MATH_WARN_SINGMATRIX, 1);

        int i, j;

        for (i = 0; i < nu; ++i)
        {
            for (j = 0; j < i; ++j)
            {
                L(i, j) = B(i, j);
                U(i, j) = 0.0;
            }

            L(i, i) = 1.0;
            U(i, i) = B(i, i);

            for (j = i+1; j < nu; ++j)
            {
                L(i, j) = 0.0;
                U(i, j) = B(i, j);
            }

            for ( ; j < n; ++j)
                U(i, j) = B(i, j);
        }

        for ( ; i < m; ++i)
        {
            for (j = 0; j < n; ++j)
                L(i, j) = B(i, j);
        }

        // convert pivot vector into a permutation vector that contains the
        // permutation matrix information
        int temp;

        for (i = 0; i < m; ++i)
            P(i) = i;

        for (i = 0; i < nu; ++i)
        {
            j = ipiv[i] - 1;    // FORTRAN is 1 based

            if (i != j)
            {
                temp = P(i);
                P(i) = P(j);
                P(j) = temp;
            }
        }

        delete [] ipiv;
    }
    else
    {
        for (int i = 0; i < m; ++i)
            P(i) = i;
    }

    return status;
}

//! Complex LU decomposition (PA = LU)
template<>
inline hwMathStatus hwTMatrix<double>::ComplexLU(hwTMatrix<double>& L, hwTMatrix<double>& U, hwTMatrix<int>& P ) const
{
    if (this == &L)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &U)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;
    int m = A.m_nRows;
    int n = A.m_nCols;
    int nu = _min(m, n);

    if (A.IsReal())
        return status(HW_MATH_ERR_NEEDCOMPLEX, 1);

    hwTMatrix<double> B(A);

    // dimension output matrices
    status = L.Dimension(m, nu, COMPLEX);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(3);
        else
            status.ResetArgs();

        return status;
    }

    status = U.Dimension(nu, n, COMPLEX);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(4);
        else
            status.ResetArgs();

        return status;
    }

    status = P.Dimension(m, hwTMatrix<int>::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(2);
        else
            status.ResetArgs();

        return status;
    }

    if (A.Size() != 0)
    {
        // decompose the matrix
        complexD* a = (complexD*) B.m_complex;
        int lda = m;
        int info;
        int* ipiv;

        try
        {
            ipiv = new int[nu];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        zgetrf_(&m, &n, a, &lda, ipiv, &info); 

        if (info != 0)
            status(HW_MATH_WARN_SINGMATRIX, 1);

        int i, j;

        for (i = 0; i < nu; ++i)
        {
            for (j = 0; j < i; ++j)
            {
                L.z(i, j) = B.z(i, j);
                U.z(i, j) = hwComplex(0.0, 0.0);
            }

            L.z(i, i) = hwComplex(1.0, 0.0);
            U.z(i, i) = B.z(i, i);

            for (j = i+1; j < nu; ++j)
            {
                L.z(i, j) = hwComplex(0.0, 0.0);
                U.z(i, j) = B.z(i, j);
            }

            for ( ; j < n; ++j)
                U.z(i, j) = B.z(i, j);
        }

        for ( ; i < m; ++i)
        {
            for (j = 0; j < n; ++j)
                L.z(i, j) = B.z(i, j);
        }

        // convert pivot vector into a permutation vector that contains the
        // permutation matrix information
        int temp;

        for (i = 0; i < m; ++i)
            P(i) = i;

        for (i = 0; i < nu; ++i)
        {
            j = ipiv[i] - 1;    // FORTRAN is 1 based

            if (i != j)
            {
                temp = P(i);
                P(i) = P(j);
                P(j) = temp;
            }
        }

        delete [] ipiv;
    }
    else
    {
        for (int i = 0; i < m; ++i)
            P(i) = i;
    }

    return status;
}

//! Real asymmetric Eigen decomposition with balance option
template<>
inline hwMathStatus hwTMatrix<double>::EigenDecompReal(bool balance, hwTMatrix<double>* V, hwTMatrix<double>& D) const
{
    if (this == V)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &D)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;
    int n = A.m_nCols;

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 0);

    if (A.m_nRows != n)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    if (n == 0)
    {
        if (V)
            status = V->Dimension(0, 0, REAL);

        status = D.Dimension(0, 0, REAL);

        return status;
    }

    char BALANC = 'N';
    char JOBVL = 'N';
    char JOBVR;
    char SENSE = 'N';

    if (V)
        JOBVR = 'V'; // right Eigen vectors
    else
        JOBVR = 'N';

    int LDA = n;
    double* VL = NULL;
    int LDVL = 1;
    int LDVR = n;
    int ILO = 0;
    int IHI = 0;
    double ABNRM = 0.0;
    int LWORK = -1; // at first call of dgeevx_

    if (balance)
        BALANC = 'B';

    hwTMatrix<double> Ain(A); // copy matrix
    hwTMatrix<double> VR;
    hwTMatrix<double> WR(n, 1, REAL);
    hwTMatrix<double> WI(n, 1, REAL);
    hwTMatrix<double> SCALE(n, 1, REAL);
    hwTMatrix<double> RCONDE(n, 1, REAL);
    hwTMatrix<double> RCONDV(n, 1, REAL);
    hwTMatrix<int> IWORK(2*n-2, 1, hwTMatrix<int>::REAL);
    double WORKSZE = 0.0;
    int INFO = 0;

    if (IWORK.M() != 2*n-2)
        return status(HW_MATH_ERR_ALLOCFAILED);

    if (V)
    {
        status = VR.Dimension(n, n, REAL);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }
    }

    double* dA = Ain.m_real;
    double* dVR = VR.m_real;
    double* dWR = WR.m_real;
    double* dWI = WI.m_real;
    double* dSCALE = SCALE.m_real;
    double* dRCONDE = RCONDE.m_real;
    double* dRCONDV = RCONDV.m_real;
    int* iIWORK = IWORK.GetRealData();

    if (!dA || !dWR || !dWI || !dSCALE || !dRCONDE || !dRCONDV)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // workspace query
    dgeevx_(&BALANC, &JOBVL, &JOBVR, &SENSE, &n, dA, &LDA, dWR, dWI,
            VL, &LDVL, dVR, &LDVR, &ILO, &IHI, dSCALE, &ABNRM,
            dRCONDE, dRCONDV, &WORKSZE, &LWORK, iIWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    LWORK = static_cast<int>(WORKSZE);

    hwTMatrix<double> WORK(LWORK, 1, REAL);

    double* dWORK = WORK.m_real;

    if (!dWORK)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // decompose the matrix
    dgeevx_(&BALANC, &JOBVL, &JOBVR, &SENSE, &n, dA, &LDA, dWR, dWI,
            VL, &LDVL, dVR, &LDVR, &ILO, &IHI, dSCALE, &ABNRM,
            dRCONDE, dRCONDV, dWORK, &LWORK, iIWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    // format decomposition factors
    bool bIsComplex = false;

    for (int i = 0; i < n; i++)
    {
        if (WI(i) != 0.0)
        {
            bIsComplex = true;
            break;
        }
    }

    if (V)
    {
        if (bIsComplex)
        {
            hwTMatrix<double> VRL(n, n, COMPLEX);
            int j = -1;

            while (1) 
            {
                j++; 

                if (j >= n)
                    break;

                if (WI(j) == 0.0)
                {
                    // real eigenvalue => v(j) = vr(j)
                    for (int i = 0; i < n; ++i) 
                        VRL.z(i, j) = VR(i, j);
                }
                else
                {
                    // complex eigenvalue => conjugate pair of eigenvectors:
                    for (int i = 0; i < n; ++i) 
                    { 
                        VRL.z(i, j).Set(VR(i, j), VR(i, j+1));
                        VRL.z(i, j+1).Set(VR(i, j), -VR(i, j+1));
                    }

                    j++; 
                }
            }

            (*V) = VRL;
        }
        else
        {
            (*V) = VR;
        }
    }

    if (bIsComplex)
    {
        status = D.PackComplex(WR, &WI);

        if (!status.IsOk())
        {
            status.SetArg1(3);
            return status;
        }
    }
    else
    {
        status = D.Dimension(n, 1, REAL);

        if (!status.IsOk())
        {
            status.SetArg1(3);
            return status;
        }

        D = WR;
    }

    return status;
}

//! Complex non-Hermitian Eigen decomposition with balance option
template<>
inline hwMathStatus hwTMatrix<double>::EigenDecompComplex(bool balance, hwTMatrix<double>* V, hwTMatrix<double>& D) const
{
    if (this == V)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &D)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;
    int n = A.m_nCols;

    if (A.m_nRows != n)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    int LDA = n;
    int LDVL = 1;
    int LDVR = n;
    int ILO = 0;
    int IHI = 0;
    int LWORK = -1;
    int INFO = 0;
    double ABNRM = 0.0;
    char BALANC = 'N';
    char JOBVL = 'N';
    char JOBVR;
    char SENSE = 'N';

    if (V)
        JOBVR = 'V'; // right Eigen vectors
    else
        JOBVR = 'N';

    complexD* VL = NULL;

    hwTMatrix<double> Ain(A); // copy matrix
    hwTMatrix<double> WORKSZE(1, 1, COMPLEX);
    hwTMatrix<double> SCALE(n, 1, REAL);
    hwTMatrix<double> RCONDE(n, 1, REAL);
    hwTMatrix<double> RCONDV(n, 1, REAL);
    hwTMatrix<double> RWORK(2*n, 1, REAL);

    if (V)
    {
        status = V->Dimension(n, n, COMPLEX);

        if (!status.IsOk())
        {
            status.SetArg1(2);
            return status;
        }
    }

    status = D.Dimension(n, 1, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(4);
        return status;
    }

    complexD* dA = (complexD*) Ain.m_complex;
    complexD* dVR;
    complexD* dW = (complexD*) D.m_complex;
    complexD* dWZ = (complexD*) WORKSZE.m_complex;
    double* dSCALE = SCALE.m_real;
    double* dRCONDE = RCONDE.m_real;
    double* dRCONDV = RCONDV.m_real;
    double* dRWORK = RWORK.m_real;

    if (V)
        dVR = (complexD*) V->m_complex;
    else
        dVR = NULL;

    if (!dA || !dW || !dWZ || !dSCALE ||
        !dRCONDE || !dRCONDV || !dRWORK)
        return status(HW_MATH_ERR_ALLOCFAILED);

    if (balance)
        BALANC = 'B';

    // workspace query
    zgeevx_(&BALANC, &JOBVL, &JOBVR, &SENSE, &n, dA, &LDA, dW,
            VL, &LDVL, dVR, &LDVR, &ILO, &IHI, dSCALE, &ABNRM,
            dRCONDE, dRCONDV, dWZ, &LWORK, dRWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    LWORK = static_cast<int>(WORKSZE.z(0).Real());

    hwTMatrix<double> WORK(LWORK, 1, COMPLEX);

    complexD* dWORK = (complexD*) WORK.m_complex;

    if (!dWORK)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // decompose the matrix
    zgeevx_(&BALANC, &JOBVL, &JOBVR, &SENSE, &n, dA, &LDA, dW,
            VL, &LDVL, dVR, &LDVR, &ILO, &IHI, dSCALE, &ABNRM,
            dRCONDE, dRCONDV, dWORK, &LWORK, dRWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Balance real matrix
template<>
inline hwMathStatus hwTMatrix<double>::BalanceReal(bool noperm, hwTMatrix<double>& S,
                                       hwTMatrix<double>& P, hwTMatrix<double>& B) const
{
    if (this == &S)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &P)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 0);

    char JOB = noperm ? (false ? 'N' : 'S') : (false ? 'P' : 'B');
    int n = A.m_nCols;
    int ILO = 0;
    int IHI = 0;
    int INFO = 0;
    hwTMatrix<double> SCALE(n, 1, REAL);

    B = A;

    double* dB = B.m_real;
    double* dSCALE = SCALE.m_real;

    if (!dB || !dSCALE)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dgebal_(&JOB, &n, dB, &n, &ILO, &IHI, dSCALE, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    status = P.Dimension(n, 1, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    for (int i = 0; i < n; i++)
        P(i) = static_cast<double>(i+1);

    for (int i = n - 1; i >= IHI; i--)
    {
        int j = static_cast<int>(SCALE(i) - 1);
        double dValueI = P(i);

        P(i) = P(j);
        P(j) = dValueI;
    }

    for (int i = 0; i < ILO - 1; i++)
    {
        int j = static_cast<int>(SCALE(i) - 1);
        double dValueI = P(i);

        P(i) = P(j);
        P(j) = dValueI;
    }

    status = S.Dimension(n, 1, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    for (int i = 0; i < ILO - 1; i++)
        S(i) = 1.0;

    for (int i = ILO - 1; i < IHI; i++)
        S(i) = SCALE(i);

    for (int i = IHI; i < n; i++)
        S(i) = 1.0;

    return status;
}

//! Balance complex matrix
template<>
inline hwMathStatus hwTMatrix<double>::BalanceComplex(bool noperm, hwTMatrix<double>& S,
                                       hwTMatrix<double>& P, hwTMatrix<double>& B) const
{
    if (this == &S)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &P)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;

    if (A.IsReal())
        return status(HW_MATH_ERR_NEEDCOMPLEX, 0);

    char JOB = noperm ? (false ? 'N' : 'S') : (false ? 'P' : 'B');
    int n = A.m_nCols;
    int ILO = 0;
    int IHI = 0;
    int INFO = 0;
    hwTMatrix<double> SCALE(n, 1, REAL);

    B = A;

    complexD* dB = (complexD*) B.m_complex;
    double* dSCALE = SCALE.m_real;

    if (!dB || !dSCALE)
        return status(HW_MATH_ERR_ALLOCFAILED);

    zgebal_(&JOB, &n, dB, &n, &ILO, &IHI, dSCALE, &INFO);

    if (INFO < 0)
        return hwMathStatus(HW_MATH_ERR_DECOMPFAIL);

    status = P.Dimension(n, 1, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    for (int i = 0; i < n; i++)
        P(i) = static_cast<double>(i+1);

    for (int i = n - 1; i >= IHI; i--)
    {
        int j = static_cast<int>(SCALE(i) - 1);
        double dValueI = P(i);

        P(i) = P(j);
        P(j) = dValueI;
    }

    for (int i = 0; i < ILO - 1; i++)
    {
        int j = static_cast<int>(SCALE(i) - 1);
        double dValueI = P(i);

        P(i) = P(j);
        P(j) = dValueI;
    }

    status = S.Dimension(n, 1, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    for (int i = 0; i < ILO - 1; i++)
        S(i) = 1.0;

    for (int i = ILO - 1; i < IHI; i++)
        S(i) = SCALE(i);

    for (int i = IHI; i < n; i++)
        S(i) = 1.0;

    return status;
}

//! Generalized real asymmetric Eigen decomposition
template<>
inline hwMathStatus hwTMatrix<double>::GeneralizedEigenDecompReal(const hwTMatrix<double>& A,
                    const hwTMatrix<double>& B, hwTMatrix<double>& V, hwTMatrix<double>& D)
{
    hwMathStatus status;
    char JOBVL = 'N';
    char JOBVR = 'V';
    int INFO = 0;
    int n = A.m_nCols;
    int LDA = n;
    int LDB = n;
    int LDVL = 1;
    int LDVR = n;
    int LWORK = -1;
    double WORKSZE = 0.0;
    double* VL = NULL;

    hwTMatrix<double> A1(A); // copy matrix
    hwTMatrix<double> B1(B); // copy matrix

    hwTMatrix<double> ALPHAR(n, 1, REAL);
    hwTMatrix<double> ALPHAI(n, 1, REAL);
    hwTMatrix<double> BETA(n, 1, REAL);

    status = V.Dimension(n, n, REAL);

    if (!status.IsOk())
    {
        status.SetArg1(3);
        return status;
    }

    status = D.Dimension(n, 1, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(4);
        return status;
    }

    double* dA = A1.m_real;
    double* dB = B1.m_real;
    double* dVR = V.m_real;
    double* dALPHAR = ALPHAR.m_real;
    double* dALPHAI = ALPHAI.m_real;
    double* dBeta = BETA.m_real;

    if (!dA || !dB || !dVR || !dALPHAR || !dALPHAI || !dBeta)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dggev_(&JOBVL, &JOBVR, &n, dA, &LDA, dB, &LDB, dALPHAR, dALPHAI,
           dBeta, VL, &LDVL, dVR, &LDVR, &WORKSZE, &LWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    LWORK = static_cast<int>(WORKSZE);
    hwTMatrix<double> WORK(LWORK, 1, REAL);

    double* dW = WORK.m_real;

    if (!dW)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dggev_(&JOBVL, &JOBVR, &n, dA, &LDA, dB, &LDB, dALPHAR, dALPHAI,
           dBeta, VL, &LDVL, dVR, &LDVR, dW, &LWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    for (int i = 0; i < n; i++)
    {
        D.z(i).Real() = ALPHAR(i) / BETA(i);
        D.z(i).Imag() = ALPHAI(i) / BETA(i);
    }

    return status;
}

//! Generalized complex Eigen decomposition
template<>
inline hwMathStatus hwTMatrix<double>::GeneralizedEigenDecompComplex(const hwTMatrix<double>& A,
                    const hwTMatrix<double>& B, hwTMatrix<double>& V, hwTMatrix<double>& D)
{
    hwMathStatus status;
    char JOBVL = 'N';
    char JOBVR = 'V';
    int n = A.m_nCols;
    int LDA = n;
    int LDB = n;
    int LDVR = n;
    int LDVL = n;
    int LWORK = -1;
    int INFO;

    hwTMatrix<double> RWORK(8*n, 1, REAL);
    hwTMatrix<double> A1(A); // copy matrix
    hwTMatrix<double> B1(B); // copy matrix

    if (A1.IsReal())
    {
        status = A1.MakeComplex();

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }
    }

    if (B1.IsReal())
    {
        status = B1.MakeComplex();

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }
    }

    status = V.Dimension(n, n, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(3);
        return status;
    }

    status = D.Dimension(n, 1, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(4);
        return status;
    }

    complexD* VL = NULL;
    complexD WORKSIZE;
    hwTMatrix<double> ALPHA(n, 1, COMPLEX);
    hwTMatrix<double> BETA(n, 1, COMPLEX);

    complexD* dA = (complexD*) A1.m_complex;
    complexD* dB = (complexD*) B1.m_complex;
    complexD* dV = (complexD*) V.m_complex;
    double* dR = RWORK.m_real;
    complexD* dALPHA = (complexD*) ALPHA.m_complex;
    complexD* dBETA = (complexD*) BETA.m_complex;

    if (!dA || !dB || !dV || !dR || !dALPHA || !dBETA)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // workspace query
    zggev_(&JOBVL, &JOBVR, &n, dA, &LDA, dB, &LDB, dALPHA, dBETA,
           VL, &LDVL, dV, &LDVR, &WORKSIZE, &LWORK, dR, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    LWORK = static_cast<int>(WORKSIZE.real());

    hwTMatrix<double> WORK(LWORK, 1, COMPLEX);

    complexD* dWORK = (complexD*) WORK.m_complex;

    if (!dWORK)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // decompose the matrix
    zggev_(&JOBVL, &JOBVR, &n, dA, &LDA, dB, &LDB, dALPHA, dBETA,
           VL, &LDVL, dV, &LDVR, dWORK, &LWORK, dR, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    for (int i = 0; i < n; i++)
    {
        hwMathStatus statusSet;
        complexD r = dALPHA[i] / dBETA[i];
        D.z(i).Set(r.real(), r.imag());

        r.real();
    }

    return status;
}

//! Generalized real symmetric Eigen decomposition
template<>
inline hwMathStatus hwTMatrix<double>::GeneralizedEigenDecompRealSymmetric(const hwTMatrix<double>& A,
                                       const hwTMatrix<double>& B, hwTMatrix<double>& V, hwTMatrix<double>& D)
{
    hwMathStatus status;
    char JOBZ = 'V';
    char UPLO = 'U';
    int ITYPE = 1;
    int N = A.m_nCols;
    int LDA = N;
    int LDB = N;
    int LWORK = -1;
    int INFO = 0;
    double WORKSIZE = 0.0;

    status = D.Dimension(N, 1, REAL);

    if (!status.IsOk())
    {
        status.SetArg1(4);
        return status;
    }

    V = A;

    hwTMatrix<double> B1(B); // copy matrix

    double* dVR = V.m_real;
    double* dB = B1.m_real;
    double* dD = D.m_real;

    if (!dVR || !dB || !dD)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // workspace query
    dsygv_(&ITYPE, &JOBZ, &UPLO, &N, dVR, &LDA, dB, &LDB, dD, &WORKSIZE, &LWORK, &INFO);

    LWORK = static_cast<int>(WORKSIZE);

    hwTMatrix<double> WORK(LWORK, 1, REAL);

    double* dW = WORK.m_real;

    if (!dW)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // decompose the matrix
    dsygv_(&ITYPE, &JOBZ, &UPLO, &N, dVR, &LDA, dB, &LDB, dD, dW, &LWORK, &INFO);

    if (INFO > N)
        return status;

    return status;
}

//! Generalized complex Hermitian Eigen decomposition
template<>
inline hwMathStatus hwTMatrix<double>::GeneralizedEigenDecompComplexHermitian(const hwTMatrix<double>& A,
                    const hwTMatrix<double>& B, hwTMatrix<double>& V, hwTMatrix<double>& D)
{
    hwMathStatus status;
    char JOBZ = 'V';
    char UPLO = 'U';
    int ITYPE = 1;
    int n = A.m_nCols;
    int LDA = n;
    int LDB = n;
    int INFO;
    int LWORK;

    hwTMatrix<double> Bin(B); // copy matrix

    status = D.Dimension(n, 1, REAL);

    if (!status.IsOk())
    {
        status.SetArg1(4);
        return status;
    }

    V = A;

    if (V.IsReal())
    {
        status = V.MakeComplex();

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }
    }

    if (Bin.IsReal())
    {
        status = Bin.MakeComplex();

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }
    }

    hwTMatrix<double> RWORK(_max(1,3*A.m_nCols-2), 1, REAL);

    LWORK = _max(1, 2*A.m_nCols-1);
 
    hwTMatrix<double> WORK(LWORK, 1, COMPLEX);

    complexD* dA = (complexD*) V.m_complex;
    complexD* dB = (complexD*) Bin.m_complex;
    double* dD = D.m_real;
    double* dRWORK = RWORK.m_real;
    complexD* dWORK = (complexD*) WORK.m_complex;

    if (!dA || !dB || !dD || !dRWORK || !dWORK)
        return status(HW_MATH_ERR_ALLOCFAILED);

    zhegv_(&ITYPE, &JOBZ, &UPLO, &n, dA, &LDA, dB,
           &LDB, dD, dWORK, &LWORK, dRWORK, &INFO);

    if (INFO > n)
        status(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Real Singular Value Decomposition
template<>
inline hwMathStatus hwTMatrix<double>::RealSVD(const SVDtype& type, hwTMatrix<double>*U,
                                               hwTMatrix<double>& S, hwTMatrix<double>* VT) const
{
    if (this == U)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &S)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == VT)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;
    int m = A.m_nRows;
    int n = A.m_nCols;

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 0);

    hwTMatrix<double> Acopy(A);
    char JOBU;
    char JOBVT;
    double* u;
    double* vt;

    if (U)
    {
        if (type == SVD_STD_FLAG)
        {
            JOBU = 'A';
            status = U->Dimension(m, m, REAL);
        }
        else
        {
            JOBU = 'S';
            status = U->Dimension(m, _min(m, n), REAL);
        }

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(2);
            else
                status.ResetArgs();

            return status;
        }

        u = U->m_real;
    }
    else
    {
        JOBU = 'N';
        u = NULL;
    }

    status = S.Dimension(_min(m, n), REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(3);
        else
            status.ResetArgs();

        return status;
    }

    if (VT)
    {
        status = VT->Dimension(n, n, REAL);
        JOBVT = 'A';

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(4);
            else
                status.ResetArgs();

            return status;
        }

        vt = VT->m_real;
    }
    else
    {
        JOBVT = 'N';
        vt = NULL;
    }

    if (A.Size() != 0)
    {
        double* a = Acopy.m_real;
        double* s = S.m_real;

        // define the arguments
        int lda = m;
        int ldu = m;
        int ldvt = n;
        int info = 0;
        int lwork = -1;
        double* work = NULL;

        try
        {
            work = new double[1];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }
        work[0] = 0.;

        // workspace query
        dgesvd_(&JOBU, &JOBVT, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);

        lwork = static_cast<int>(work[0]);
        delete [] work;

        try
        {
            work = new double[lwork];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        // decompose the matrix
        dgesvd_(&JOBU, &JOBVT, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);

        delete [] work;

        if (info != 0)
            return status(HW_MATH_ERR_NOTCONVERGE);
    }
    else
    {
        if (U)
            U->Identity();

        if (VT)
            VT->Identity();
    }

    if (VT)
    {
        // transpose result for compatibility
        status = VT->Transpose();

        if (type == SVD_ECON_FLAG)
        {
            int nrows = A.m_nCols;
            int ncols = A.m_nRows;
            int mindim = (nrows < ncols) ? nrows : ncols;
            int idx = VT->m_nRows - mindim;
            if (idx > 0)
            {
                status = VT->DeleteColumns(mindim, idx);
                if (!status.IsOk())
                {
                    status.ResetArgs();
                    return status;
                }
            }
        }
    }

    return status;
}

//! Complex Singular Value Decomposition
template<>
inline hwMathStatus hwTMatrix<double>::ComplexSVD(const SVDtype& type, hwTMatrix<double>* U,
                                                  hwTMatrix<double>& S, hwTMatrix<double>* VT) const
{
    if (this == U)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &S)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == VT)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwTMatrix<double> Acopy(A);
    hwMathStatus status;
    int m = A.m_nRows;
    int n = A.m_nCols;
    char JOBU;
    char JOBVT;
    complexD* u = NULL;
    complexD* vt = NULL;

    if (U)
    {
        if (type == SVD_STD_FLAG)
        {
            JOBU = 'A';
            status = U->Dimension(m, m, COMPLEX);
        }
        else
        {
            JOBU = 'S';
            status = U->Dimension(m, _min(m, n), COMPLEX);
        }

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(2);
            else
                status.ResetArgs();

            return status;
        }

        u = (complexD*) U->m_complex;
    }
    else
    {
        JOBU = 'N';
        u = NULL;
    }

    status = S.Dimension(_min(m, n), REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(3);
        else
            status.ResetArgs();

        return status;
    }

    if (VT)
    {
        JOBVT = 'A';
        status = VT->Dimension(n, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(4);
            else
                status.ResetArgs();

            return status;
        }

        vt = (complexD*) VT->m_complex;
    }
    else
    {
        JOBVT = 'N';
        vt = NULL;
    }

    if (A.Size() != 0)
    {
        complexD* a = (complexD*) Acopy.m_complex;
        double* s = S.m_real;

        // define the arguments
        int LDA = m;
        int LDU = m;
        int LDVT = n;
        int INFO = 0;
        int LWORK = -1;
        int nrows = A.m_nCols;
        int ncols = A.m_nRows;
        int mindim = (nrows < ncols) ? nrows : ncols;
        hwTMatrix<double> RWORK(5 * mindim, 1, REAL);
        double *drwork = RWORK.m_real;
        complexD *WORKC = NULL;

        try
        {
            // 2 is minimal size
            // LWORK >=  MAX(1,2*MIN(M,N)+MAX(M,N))
            WORKC = new complexD[2];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        // workspace query
        zgesvd_(&JOBU, &JOBVT, &m, &n, a, &LDA, s, u, &LDU, vt, &LDVT, WORKC,  &LWORK, drwork, &INFO);

        LWORK = static_cast<int>(WORKC[0].real());
        delete [] WORKC;
        if (INFO < 0)
            return status(HW_MATH_ERR_NOTCONVERGE);

        try
        {
            WORKC = new complexD[LWORK * 2];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        // decompose the matrix
        zgesvd_(&JOBU, &JOBVT, &m, &n, a, &LDA, s, u, &LDU, vt, &LDVT, WORKC,  &LWORK, drwork, &INFO);
        delete [] WORKC;
        
        if (INFO != 0)
            return status(HW_MATH_ERR_NOTCONVERGE);
    }
    else
    {
        U->Identity();
        VT->Identity();
    }

    if (VT)
    {
        // hermitian result for compatibility
        VT->Conjugate();
        status = VT->Transpose();

        if (type == SVD_ECON_FLAG)
        {
            int nrows = A.m_nCols;
            int ncols = A.m_nRows;
            int mindim = (nrows < ncols) ? nrows : ncols;
            int idx = VT->m_nRows - mindim;
            if (idx > 0)
            {
                status = VT->DeleteColumns(*VT, mindim, idx);
                if (!status.IsOk())
                {
                    status.ResetArgs();
                    return status;
                }
            }
        }
    }

    return status;
}

//! Real QR decomposition (A = QR)
template<>
inline hwMathStatus hwTMatrix<double>::RealQR(hwTMatrix<double>& Q, hwTMatrix<double>& R) const
{
    if (this == &Q)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &R)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;
    int m = A.m_nRows;
    int n = A.m_nCols;

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    hwTMatrix<double> Acopy(A);

    if (A.Size() != 0)
    {
        // define the arguments
        double* q = Acopy.m_real;
        int k = _min(m, n);
        int lda = m;
        int info;
        int lwork = -1;
        double* work;
        double* tau;

        // workspace query
        try
        {
            work = new double[1];
            tau = new double[k];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        dgeqrf_(&m, &n, q, &lda, tau, work, &lwork, &info);

        lwork = static_cast<int>(work[0]);
        delete [] work;

        try
        {
            work = new double[lwork];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        // decompose the matrix
        dgeqrf_(&m, &n, q, &lda, tau, work, &lwork, &info);

        if (info != 0)
        {
            delete [] work;
            delete [] tau;

            return status(HW_MATH_ERR_NOTCONVERGE);
        }

        // construct R
        status = R.Dimension(k, n, REAL);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(3);
            else
                status.ResetArgs();

            return status;
        }

        for (int i = 0; i < k; ++i)
        {
            for (int j = 0; j < i; ++j)
                R(i, j) = 0.0;

            for (int j = i; j < n; ++j)
                R(i, j) = Acopy(i, j);
        }

        // workspace query
        delete [] work;

        lwork = -1;

        try
        {
            work = new double[1];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        dorgqr_(&m, &k, &k, q, &lda, tau, work, &lwork, &info);

        lwork = static_cast<int>(work[0]);
        delete [] work;

        try
        {
            work = new double[lwork];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        // construct Q
        dorgqr_(&m, &k, &k, q, &lda, tau, work, &lwork, &info);

        delete [] work;
        delete [] tau;

        if (info != 0)
            status(HW_MATH_ERR_NOTCONVERGE);

        status = Q.Dimension(m, k, REAL);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(2);
            else
                status.ResetArgs();

            return status;
        }

        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < k; ++j)
                Q(i, j) = Acopy(i, j);
        }
    }
    else if (m != 0)
    {
        status = Q.Dimension(m, m, REAL);
        Q.Identity();
        status = R.Dimension(m, 0, REAL);
    }
    else
    {
        status = Q.Dimension(0, 0, REAL);
        status = R.Dimension(0, n, REAL);
    }

    return status;
}

//! Complex QR decomposition (A = QR)
template<>
inline hwMathStatus hwTMatrix<double>::ComplexQR(hwTMatrix<double>& Q, hwTMatrix<double>& R) const
{
    if (this == &Q)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &R)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;
    int m = A.m_nRows;
    int n = A.m_nCols;

    if (A.IsReal())
        return status(HW_MATH_ERR_NEEDCOMPLEX, 1);

    hwTMatrix<double> Acopy(A);

    if (A.Size() != 0)
    {
        // define the arguments
        complexD* q = (complexD*) Acopy.m_complex;
        int k = _min(m, n);
        int lda = m;
        int info;
        int lwork = -1;
        complexD* work = new complexD[1];
        complexD* tau = new complexD[k];

        // workspace query
        try
        {
            work = new complexD[1];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        zgeqrf_(&m, &n, q, &lda, tau, work, &lwork, &info);

        lwork = static_cast<int>(work[0].real());
        delete [] work;

        try
        {
            work = new complexD[lwork];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        // decompose the matrix
        zgeqrf_(&m, &n, q, &lda, tau, work, &lwork, &info);

        if (info != 0)
        {
            delete [] work;
            delete [] tau;

            return status(HW_MATH_ERR_NOTCONVERGE);
        }

        // construct R
        status = R.Dimension(k, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(3);
            else
                status.ResetArgs();

            return status;
        }

        for (int i = 0; i < k; ++i)
        {
            for (int j = 0; j < i; ++j)
                R.z(i, j) = 0.0;

            for (int j = i; j < n; ++j)
                R.z(i, j) = Acopy.z(i, j);
        }

        // workspace query
        delete [] work;

        lwork = -1;

        try
        {
            work = new complexD[1];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        zungqr_(&m, &k, &k, q, &lda, tau, work, &lwork, &info);

        lwork = static_cast<int>(work[0].real());
        delete [] work;

        try
        {
            work = new complexD[lwork];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        // construct Q
        zungqr_(&m, &k, &k, q, &lda, tau, work, &lwork, &info);

        delete [] work;
        delete [] tau;

        if (info != 0)
            status(HW_MATH_ERR_NOTCONVERGE);

        status = Q.Dimension(m, k, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(2);
            else
                status.ResetArgs();

            return status;
        }

        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < k; ++j)
                Q.z(i, j) = Acopy.z(i, j);
        }
    }

    return status;
}

//! Real matrix pseudo-inversion (Moore-Penrose)
template<>
inline hwMathStatus hwTMatrix<double>::RealPinv(const hwTMatrix<double>& source)
{
    if (this == &source)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<double>& pinv = (*this);
    hwMathStatus status;

    if (source.IsEmpty())
    {
        status = pinv.Dimension(source.m_nCols, source.m_nRows, REAL);
        return status;
    }

    if (!source.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    int m = source.m_nRows;
    int n = source.m_nCols;
    int nu = _min(m, n);
    int lda = m;
    int nrhs = m;
    int ldb = _max(m, n);
    int rank;
    int info;
    double rcond = _max(m, n) * MACHEP2;
    double* s;

    try
    {
        s = new double[nu];
    }
    catch (std::bad_alloc&) 
    {
        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    hwTMatrix<double> Acopy(source);

    status = pinv.Dimension(ldb, m, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    pinv.SetElements(0.0);

    for (int i = 0; i < m; ++i)
        pinv(i, i) = 1.0;

    double* a_r = Acopy.m_real;
    double* b_r = pinv.m_real;

    // workspace query
    double* work;

    try
    {
        work = new double[1];
    }
    catch (std::bad_alloc&) 
    {
        if (s)
            delete [] s;

        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    int lwork = -1;

    dgelss_(&m, &n, &nrhs, a_r, &lda, b_r, &ldb, s, &rcond, &rank,
            work, &lwork, &info);

    lwork = static_cast<int>(work[0]);
    delete [] work;

    try
    {
        work = new double[lwork];
    }
    catch (std::bad_alloc&) 
    {
        if (s)
            delete [] s;

        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    // decompose the matrix and the solve system
    dgelss_(&m, &n, &nrhs, a_r, &lda, b_r, &ldb, s, &rcond, &rank,
            work, &lwork, &info);

    delete [] s;
    delete [] work;

    if (info != 0)
    {
        status(HW_MATH_ERR_NOTCONVERGE);
        return status;
    }

    status = pinv.Resize(n, m, false);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    return status;
}

//! Complex matrix pseudo-inversion (Moore-Penrose)
template<>
inline hwMathStatus hwTMatrix<double>::ComplexPinv(const hwTMatrix<double>& source)
{
    if (this == &source)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<double>& pinv = (*this);
    hwMathStatus status;

    if (source.IsEmpty())
    {
        status = pinv.Dimension(source.m_nCols, source.m_nRows, REAL);
        return status;
    }

    if (source.IsReal())
        return status(HW_MATH_ERR_NEEDCOMPLEX, 1);

    int m = source.m_nRows;
    int n = source.m_nCols;
    int nu = _min(m, n);
    int lda = m;
    int nrhs = m;
    int ldb = _max(m, n);
    int rank;
    int info;
    double rcond = _max(m, n) * MACHEP2;
    double* s;

    try
    {
        s = new double[nu];
    }
    catch (std::bad_alloc&) 
    {
        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    hwTMatrix<double> Acopy(source);

    status = pinv.Dimension(ldb, m, COMPLEX);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    pinv.SetElements(0.0);

    for (int i = 0; i < m; ++i)
        pinv.z(i, i) = 1.0;

    complexD* a_r = (complexD*) Acopy.m_complex;
    complexD* b_r = (complexD*) pinv.m_complex;

    // workspace query
    double* rwork;
    complexD* work;

    try
    {
        rwork = new double[5*static_cast<long int>(nu)];
        work  = new complexD[1];
    }
    catch (std::bad_alloc&) 
    {
        if (s)
            delete [] s;

        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    int lwork = -1;

    zgelss_(&m, &n, &nrhs, a_r, &lda, b_r, &ldb, s, &rcond, &rank,
            work, &lwork, rwork, &info);

    lwork = static_cast<int>(work[0].real());
    delete [] work;

    try
    {
        work = new complexD[lwork];
    }
    catch (std::bad_alloc&) 
    {
        if (s)
            delete [] s;

        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    // decompose the matrix and the solve system
    zgelss_(&m, &n, &nrhs, a_r, &lda, b_r, &ldb, s, &rcond, &rank,
            work, &lwork, rwork, &info);

    delete [] s;
    delete [] work;

    if (info != 0)
    {
        status(HW_MATH_ERR_NOTCONVERGE);
        return status;
    }

    status = pinv.Resize(n, m, false);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    return status;
}

//! Matrix pseudo-inversion (Moore-Penrose)
template<>
inline hwMathStatus hwTMatrix<double>::Pinv(const hwTMatrix<double>& source)
{
    if (source.IsReal())
        return RealPinv(source);
    else
        return ComplexPinv(source);
}

//! Real Schur Decomposition
template<>
inline hwMathStatus hwTMatrix<double>::RealSchur(hwTMatrix<double>& U,
                                       hwTMatrix<double>& T) const
{
    if (this == &U)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &T)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 0);

    if (!A.IsSquare())
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    char JOBVS = 'V';
    char SORT = 'N';
    int n = A.m_nCols;
    int LDA = _max(n, 1);
    int SDIM = 0;
    int LDVS = n;
    int LWORK = 0;
    int LWORKMIN = -1;
    double WORKOPTIMAL = 0;
    int INFO = 0;

    hwTMatrix<double> WR(n, 1, REAL);
    hwTMatrix<double> WI(n, 1, REAL);
    hwTMatrix<int> BWORK(n, 1, hwTMatrix<int>::REAL);  

    status = U.Dimension(n, n, REAL);

    if (!status.IsOk())
    {
        status.SetArg1(1);
        return status;
    }

    T = A;

    double* dT = T.m_real;
    double* dWR = WR.m_real;
    double* dWI = WI.m_real;
    double* dU = U.m_real;
    int* iBWORK = BWORK.GetRealData();
	memset(iBWORK, 0, sizeof(int) * BWORK.Size());
    
    // workspace query
    dgees_(&JOBVS, &SORT, NULL, &n, dT, &LDA, &SDIM, dWR, dWI,
           dU, &LDVS, &WORKOPTIMAL, &LWORKMIN, iBWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    LWORK = static_cast<int>(WORKOPTIMAL);

    hwTMatrix<double> WORK(LWORK, 1, REAL);
	WORK.SetElements(0.0);
	memset(iBWORK, 0, sizeof(int) * BWORK.Size());

    double* dWORK = WORK.m_real;

    // decompose matrix
    dgees_(&JOBVS, &SORT, NULL, &n, dT, &LDA, &SDIM, dWR, dWI,
           dU, &LDVS, dWORK, &LWORK, iBWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    return status;  
}

//! Complex Schur Decomposition
template<>
inline hwMathStatus hwTMatrix<double>::ComplexSchur(hwTMatrix<double>& U,
                                       hwTMatrix<double>& T) const
{
    if (this == &U)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &T)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;

    if (!A.IsSquare())
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    char JOBVS = 'V';
    char SORT = 'N';
    int n = A.m_nCols;
    int LDA = _max(n, 1);
    int SDIM = 0;
    int LDVS = n;
    int LWORKMIN = -1;
    int LWORK = 0;
    int INFO = 0;

    status = U.Dimension(n, n, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(1);
        return status;
    }

    T = A;

    if (T.IsReal())
    {
        status = T.MakeComplex();

        if (!status.IsOk())
        {
            status.SetArg1(2);
            return status;
        }
    }

    hwTMatrix<double> W(n, COMPLEX);
    hwTMatrix<double> WORK(1, COMPLEX);
    hwTMatrix<double> RWORK(n, 1, REAL);
    hwTMatrix<int> BWORK(n, 1, hwTMatrix<int>::REAL);

    complexD* dT = (complexD*) T.m_complex;
    complexD* dW = (complexD*) W.m_complex;
    complexD* dU = (complexD*) U.m_complex;
    complexD* dWORK = (complexD*) WORK.m_complex;
	memset(dWORK, 0, sizeof(complexD) * WORK.Size());

    double* dRWORK = RWORK.m_real;
    int* iBWORK = BWORK.GetRealData();
	memset(iBWORK, 0, sizeof(int) * BWORK.Size());

    // workspace query
    zgees_(&JOBVS, &SORT, NULL, &n, dT, &LDA, &SDIM, dW,
           dU, &LDVS, dWORK, &LWORKMIN, dRWORK, iBWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    if (INFO == 0)
        LWORK = _max(static_cast<int>(WORK.z(0).Real()), 2*n);
    else
        LWORK = _max(1, 2*n);

    status = WORK.Resize(LWORK);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

	memset(iBWORK, 0, sizeof(int) * BWORK.Size());

    dWORK = (complexD*) WORK.m_complex;
	memset(dWORK, 0, sizeof(complexD) * WORK.Size());

    // decompose matrix
    zgees_(&JOBVS, &SORT, NULL, &n, dT, &LDA, &SDIM, dW,
           dU, &LDVS, dWORK, &LWORK, dRWORK, iBWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    return status;  
}

//! Reciprocal condition number estimate of a real matrix
template<>
inline hwMathStatus hwTMatrix<double>::RealRCond(double& rCondNum) const
{
    hwMathStatus status;

    for (int i = 0; i < Size(); ++i)
    {
        if (IsNaN_T(m_real[i]))
        {
            rCondNum = m_real[i];
            return status;
        }
    }

    if (Size() == 1)  
    {
        double real = m_real[0];

        if (!IsNaN_T(real) && !IsInf_T(real) && !IsNegInf_T(real))
        {
            if (real == 0.0)
                rCondNum = 0.0;
            else
                rCondNum = 1.0;

            return status;
        }
    }

    hwTMatrix<double> Ain(*this); // copy matrix

    double* dA = Ain.m_real;
    double Anorm;
    char NORM = '1';
    int m = m_nRows;
    int n = m_nCols;
    int INFO = 0;

    Anorm = dlange_(&NORM, &m, &n, dA, &m, NULL);

    hwTMatrix<int> IPIV(_min(m, n), 1, hwTMatrix<int>::REAL);

    int* iPiv = IPIV.GetRealData();

    if (!dA || !iPiv)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dgetrf_(&m, &n, dA, &m, iPiv, &INFO);

    if (INFO < 0)
        return hwMathStatus(HW_MATH_ERR_DECOMPFAIL);

    INFO = 0; 
    hwTMatrix<double> WORK(4 * n, 1, REAL);
    hwTMatrix<int> IWORK(n, 1, hwTMatrix<int>::REAL);

    double* dW = WORK.m_real;
    int* iW = IWORK.GetRealData();

    if (!dW || !iW)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dgecon_(&NORM, &n, dA, &n, &Anorm, &rCondNum, dW, iW, &INFO);

    if (INFO < 0)
        return hwMathStatus(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Reciprocal condition number estimate of a complex matrix
template<>
inline hwMathStatus hwTMatrix<double>::ComplexRCond(double& rCondNum) const
{
    hwMathStatus status;
    
    for (int i = 0; i < Size(); ++i)
    {
        if (IsNaN_T(m_complex[i].Real()))
        {
            rCondNum = m_complex[i].Real();
            return status;
        }

        if (IsNaN_T(m_complex[i].Imag()))
        {
            rCondNum = m_complex[i].Imag();
            return status;
        }
    }

    if (Size() == 1)  
    {
        double real = m_complex[0].Real();
        double imag = m_complex[0].Imag();

        if ((!IsNaN_T(real) && !IsInf_T(real) && !IsNegInf_T(real)) || 
            (!IsNaN_T(imag) && !IsInf_T(imag) && !IsNegInf_T(imag)))
        {
            if (real == 0.0 && imag == 0.0)
                rCondNum = 0.0;
            else
                rCondNum = 1.0;

            return status;
        }
    }

    hwTMatrix<double> Ain(*this); // copy matrix
    double Anorm;
    char NORM = '1';
    int m = Ain.m_nRows;
    int n = Ain.m_nCols;
    int INFO = 0;

    complexD* dA = (complexD*) Ain.m_complex;

    Anorm = zlange_(&NORM, &m, &n, dA, &m, NULL);
    
    hwTMatrix<int> IPIV(_min(m, n), 1, hwTMatrix<int>::REAL);

    int* iPiv = IPIV.GetRealData();

    if (!dA || !iPiv)
        return status(HW_MATH_ERR_ALLOCFAILED);

    zgetrf_(&m, &n, dA, &m, iPiv, &INFO);

    if (INFO < 0)
        return hwMathStatus(HW_MATH_ERR_DECOMPFAIL);

    INFO = 0; 

    hwTMatrix<double> WORK(4 * n, 1, COMPLEX);
    hwTMatrix<double> RWORK(2 * n, 1, REAL);

    complexD* dW = (complexD*) WORK.m_complex;
    double* dRW = RWORK.m_real;

    if (!dW || !dRW)
        return status(HW_MATH_ERR_ALLOCFAILED);

    zgecon_(&NORM, &n, dA, &m, &Anorm, &rCondNum, dW, dRW, &INFO);

    if (INFO < 0)
        return hwMathStatus(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//*******************************************************************
//            hwTMatrix<double> public implementations
//*******************************************************************

//! Reciprocal condition number estimate of a real matrix
template<>
inline hwMathStatus hwTMatrix<double>::RCond(double& rCondNum) const
{
    hwMathStatus status;

    if (m_nRows != m_nCols)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (IsEmpty())
    {
        rCondNum = std::numeric_limits<double>::infinity();
        return status;
    }

    if (IsReal())
        status = RealRCond(rCondNum);
    else
        status = ComplexRCond(rCondNum);

    return status;
}

//! Solve linear system AX=B with QR decomposition, where X = *this
template<>
inline hwMathStatus hwTMatrix<double>::QRSolve(const hwTMatrix<double>& A, const hwTMatrix<double>& B)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<double>& X = (*this);
    hwMathStatus status;
    int m = A.m_nRows;
    int n = A.m_nCols;
    int bn = B.m_nCols;
    int lda = m;
    int ldb = _max(m, n);
    int info;
    int lwork = -1;
    int rank;
    int* jpvt;
    double rcond = 1.0e-12;
    // char TRANS = 'N';

    if (B.m_nRows != m)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (A.IsReal() && B.IsReal())
    {
        // create copies of A and B;
        hwTMatrix<double> Acopy(A);
        hwTMatrix<double> Bcopy(ldb, bn, REAL);

        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < bn; ++j)
                Bcopy(i, j) = B(i, j);
        }

        status = X.Dimension(n, bn, REAL);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        if (A.IsEmpty())
        {
            X.SetElements(0.0);
            return status;
        }

        double* a = Acopy.m_real;
        double* b = Bcopy.m_real;

        // workspace query
        double* work;

        try
        {
            jpvt = new int[n];
            work = new double[1];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        memset(jpvt, 0, n*sizeof(int));

        // dgels_(&TRANS, &m, &n, &bn, a, &lda, b, &ldb, work, &lwork, &info);
        dgelsy_(&m, &n, &bn, a, &lda, b, &ldb, jpvt, &rcond, &rank,
                work, &lwork, &info);

        lwork = static_cast<int>(work[0]);
        delete [] work;

        try
        {
            work = new double[lwork];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        // decompose the matrix and solve
        // dgels_(&TRANS, &m, &n, &bn, a, &lda, b, &ldb, work, &lwork, &info);
        dgelsy_(&m, &n, &bn, a, &lda, b, &ldb, jpvt, &rcond, &rank,
                work, &lwork, &info);

        if (info != 0 || rank < _min(m, n))
            status(HW_MATH_WARN_MTXNOTFULLRANK, 1);

        delete [] jpvt;
        delete [] work;

        // copy solution back into X
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < bn; ++j)
                X(i, j) = Bcopy(i, j);
        }
    }
    else
    {
        status = X.Dimension(n, bn, COMPLEX);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        // create copies of A and B;
        hwTMatrix<double> Acopy(m, n, COMPLEX);
        hwTMatrix<double> Bcopy(ldb, bn, COMPLEX);

        if (!A.IsReal() && !B.IsReal())
        {
            Acopy = A;

            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < bn; ++j)
                    Bcopy.z(i, j) = B.z(i, j);
            }
        }
        else if (A.IsReal() && !B.IsReal())
        {
            status = Acopy.PackComplex(A);

            if (!status.IsOk())
            {
                status.ResetArgs();
                return status;
            }

            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < bn; ++j)
                    Bcopy.z(i, j) = B.z(i, j);
            }
        }
        else // if (!A.IsReal() && B.IsReal())
        {
            Acopy = A;

            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < bn; ++j)
                    Bcopy.z(i, j) = B(i, j);
            }
        }

        complexD* a = (complexD*) Acopy.m_complex;
        complexD* b = (complexD*) Bcopy.m_complex;

        // workspace query
        double* rwork;
        complexD* work;

        try
        {
            jpvt  = new int[n];
            rwork = new double[2 * static_cast<long int>(n)];
            work  = new complexD[1];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        memset(jpvt, 0, n*sizeof(int));

        // zgels_(&TRANS, &m, &n, &bn, a, &lda, b, &ldb, work, &lwork, &info);
        zgelsy_(&m, &n, &bn, a, &lda, b, &ldb, jpvt, &rcond, &rank,
                work, &lwork, rwork, &info);

        lwork = static_cast<int>(work[0].real());
        delete [] work;

        try
        {
            work = new complexD[lwork];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        // decompose the matrix and solve
        // zgels_(&TRANS, &m, &n, &bn, a, &lda, b, &ldb, work, &lwork, &info);
        zgelsy_(&m, &n, &bn, a, &lda, b, &ldb, jpvt, &rcond, &rank,
                work, &lwork, rwork, &info);

        if (info != 0 || rank < _min(m, n))
            status(HW_MATH_WARN_MTXNOTFULLRANK, 1);

        delete [] jpvt;
        delete [] rwork;
        delete [] work;

        // copy solution back into X
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < bn; ++j)
                X.z(i, j) = Bcopy.z(i, j);
        }
    }

    return status;
}

//! Solve linear system AX=B, where X = *this
//! uses LU decomposition when square, and QR otherwise
template<>
inline hwMathStatus hwTMatrix<double>::LSolve(const hwTMatrix<double>& A, const hwTMatrix<double>& B)
{
    hwMathStatus status;
    hwMathStatus status2;

    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<double>& X = (*this);

    // check dimensions
    int m = A.m_nRows;
    int n = A.m_nCols;
    int bn = B.m_nCols;

    if (m != n)
        return QRSolve(A, B);

    if (B.m_nRows != m)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    // check matrix condition
    double rcond;

    status = A.RCond(rcond);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (rcond < MACHEP2)
    {
        status2(HW_MATH_WARN_SINGMATRIXDIV, 1);
    }

    int lda = m;
    int ldb = m;
    int* ipiv;
    int info;
    hwTMatrix<double> Acopy(A);

    try
    {
        ipiv = new int[m];
    }
    catch (std::bad_alloc&) 
    {
        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    if (A.IsReal() && B.IsReal())
    {
        status = X.Dimension(n, bn, REAL);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        if (X.IsEmpty())
            return status;

        X = B;

        double* x_r = X.m_real;
        double* a_r = Acopy.m_real;

        dgesv_(&n, &bn, a_r, &lda, ipiv, x_r, &ldb, &info);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        status = X.Dimension(n, bn, COMPLEX);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        X = B;

        complexD* a_c = (complexD*) Acopy.m_complex;
        complexD* x_c = (complexD*) X.m_complex;

        zgesv_(&n, &bn, a_c, &lda, ipiv, x_c, &ldb, &info);
    }
    else if (A.IsReal() && !B.IsReal())
    {
        status = X.Dimension(n, bn, COMPLEX);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        X = B;

        status = Acopy.MakeComplex();

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        complexD* a_c = (complexD*) Acopy.m_complex;
        complexD* x_c = (complexD*) X.m_complex;

        zgesv_(&n, &bn, a_c, &lda, ipiv, x_c, &ldb, &info);
    }
    else // if (!A.IsReal() && B.IsReal())
    {
        status = X.Dimension(n, bn, COMPLEX);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        status = X.PackComplex(B);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        complexD* a_c = (complexD*) Acopy.m_complex;
        complexD* x_c = (complexD*) X.m_complex;

        zgesv_(&n, &bn, a_c, &lda, ipiv, x_c, &ldb, &info);
    }

    delete [] ipiv;

    if (info != 0)
    {
        status(HW_MATH_WARN_SINGMATRIXDIV, 1);
        X.SetElements(std::numeric_limits<double>::infinity());
    }
    else if (!status2.IsOk())   // preserve earlier warning
        status = status2;

    return status;
}

//! Add two matrices so that (*this) = A + B
template<>
inline hwMathStatus hwTMatrix<double>::Add(const hwTMatrix<double>& A, const hwTMatrix<double>& B)
{
    hwMathStatus status;

    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    // check dimensions
    int m = A.m_nRows;
    int n = A.m_nCols;
    int size = A.Size();

    if (B.m_nRows != m || B.m_nCols != n)
    {
        if (size == 0 && B.Size() == 1)
            status = Dimension(A.m_nRows, A.m_nCols, REAL);
        else if (size == 1 && B.Size() == 0)
            status = Dimension(B.m_nRows, B.m_nCols, REAL);
        else
            status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        return status;
    }

    // prepare for LAPACK function call
    int inc = 1;
    double a[2] = {1.0, 0.0};   // real or complex

    if (A.IsReal() && B.IsReal())
        status = Dimension(m, n, REAL);
    else
        status = Dimension(m, n, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(0);
        return status;
    }

    if (A.IsReal() && B.IsReal())
    {
        if (size != 0)
        {
            double* t_r = m_real;
            double* a_r = (double*) A.m_real;
            double* b_r = (double*) B.m_real;

            dcopy_(&size, a_r, &inc, t_r, &inc);
            daxpy_(&size, a, b_r, &inc, t_r, &inc);
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        complexD* t_c = (complexD*) m_complex;
        complexD* a_c = (complexD*) A.m_complex;
        complexD* b_c = (complexD*) B.m_complex;

        zcopy_(&size, a_c, &inc, t_c, &inc);
        zaxpy_(&size, (complexD*) a, b_c, &inc, t_c, &inc);
    }
    else if (A.IsReal() && !B.IsReal())
    {
        double* t_c = (double*) m_complex;
        double* a_r = A.m_real;
        double* b_c = (double*) B.m_complex;

        int inc2 = 2;
        int size2 = size<<1;
        dcopy_(&size2, b_c, &inc, t_c, &inc);
        daxpy_(&size, a, a_r, &inc, t_c, &inc2);
    }
    else // if (!A.IsReal() && B.IsReal())
    {
        double* t_c = (double*) m_complex;
        double* a_c = (double*) A.m_complex;
        double* b_r = B.m_real;

        int inc2 = 2;
        int size2 = size<<1;
        dcopy_(&size2, a_c, &inc, t_c, &inc);
        daxpy_(&size, a, b_r, &inc, t_c, &inc2);
    }

    return status;
}

//! Add a matrix to the calling object so that (this) += A
template<>
inline hwMathStatus hwTMatrix<double>::AddEquals(const hwTMatrix<double>& A)
{
    hwMathStatus status;

    if (this == &A)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    // check dimensions
    int m = m_nRows;
    int n = m_nCols;
    int size = Size();

    if (A.m_nRows != m || A.m_nCols != n)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (IsReal() && !A.IsReal())
    {
        status = MakeComplex();

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }
    }

    // prepare for LAPACK function call
    int inc = 1;
    double a[2] = {1.0, 0.0};   // real or complex

    if (IsReal() && A.IsReal())
    {
        if (size != 0)
        {
            double* t_r = m_real;
            double* a_r = (double*) A.m_real;

            daxpy_(&size, a, a_r, &inc, t_r, &inc);
        }
    }
    else if (!IsReal() && !A.IsReal())
    {
        complexD* t_c = (complexD*) m_complex;
        complexD* a_c = (complexD*) A.m_complex;

        zaxpy_(&size, (complexD*) a, a_c, &inc, t_c, &inc);
    }
    else // if (!IsReal() && A.IsReal())
    {
        double* t_c = (double*) m_complex;
        double* a_r = (double*) A.m_real;

        int inc2 = 2;
        daxpy_(&size, a, a_r, &inc, t_c, &inc2);
    }

    return status;
}

//! Subtract one matrix from another so that (*this) = A - B
template<>
inline hwMathStatus hwTMatrix<double>::Subtr(const hwTMatrix<double>& A, const hwTMatrix<double>& B)
{
    hwMathStatus status;

    if (this == &A)
        return status(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    // check dimensions
    int m = A.m_nRows;
    int n = A.m_nCols;
    int size = A.Size();

    if (B.m_nRows != m || B.m_nCols != n)
    {
        if (size == 0 && B.Size() == 1)
            status = Dimension(A.m_nRows, A.m_nCols, REAL);
        else if (size == 1 && B.Size() == 0)
            status = Dimension(B.m_nRows, B.m_nCols, REAL);
        else
            status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        return status;
    }

    // prepare for LAPACK function call
    int inc = 1;
    double a[2] = {-1.0, 0.0};   // real or complex

    if (A.IsReal() && B.IsReal())
        status = Dimension(m, n, REAL);
    else
        status = Dimension(m, n, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(0);
        return status;
    }

    if (A.IsReal() && B.IsReal())
    {
        if (size != 0)
        {
            double* t_r = m_real;
            double* a_r = A.m_real;
            double* b_r = B.m_real;

            dcopy_(&size, a_r, &inc, t_r, &inc);
            daxpy_(&size, a, b_r, &inc, t_r, &inc);
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        complexD* t_c = (complexD*) m_complex;
        complexD* a_c = (complexD*) A.m_complex;
        complexD* b_c = (complexD*) B.m_complex;

        zcopy_(&size, a_c, &inc, t_c, &inc);
        zaxpy_(&size, (complexD*) a, b_c, &inc, t_c, &inc);
    }
    else if (A.IsReal() && !B.IsReal())
    {
        complexD* t_c = (complexD*) m_complex;
        double* a_r = A.m_real;
        complexD* b_c = (complexD*) B.m_complex;

        SetElements(0.0);

        int inc2 = 2;
        dcopy_(&size, a_r, &inc, (double*) t_c, &inc2);
        zaxpy_(&size, (complexD*) a, b_c, &inc, t_c, &inc);
    }
    else // if (!A.IsReal() && B.IsReal())
    {
        complexD* t_c = (complexD*) m_complex;
        complexD* a_c = (complexD*) A.m_complex;
        double* b_r = B.m_real;

        int inc2 = 2;
        zcopy_(&size, a_c, &inc, t_c, &inc);
        daxpy_(&size, a, b_r, &inc, (double*) t_c, &inc2);
    }

    return status;
}

//! Subtract a matrix from the calling object so that (this) -= A
template<>
inline hwMathStatus hwTMatrix<double>::SubtrEquals(const hwTMatrix<double>& A)
{
    hwMathStatus status;

    if (this == &A)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    // check dimensions
    int m = m_nRows;
    int n = m_nCols;
    int size = Size();

    if (A.m_nRows != m || A.m_nCols != n)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    // prepare for LAPACK function call
    int inc = 1;
    double a[2] = {-1.0, 0.0};   // real or complex

    if (IsReal() && A.IsReal())
    {
        if (size != 0)
        {
            double* t_r = m_real;
            double* a_r = A.m_real;

            daxpy_(&size, a, a_r, &inc, t_r, &inc);
        }
    }
    else if (!IsReal() && !A.IsReal())
    {
        complexD* t_c = (complexD*) m_complex;
        complexD* a_c = (complexD*) A.m_complex;

        zaxpy_(&size, (complexD*) &a, a_c, &inc, t_c, &inc);
    }
    else if (IsReal() && !A.IsReal())
    {
        status = MakeComplex();

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }

        complexD* t_c = (complexD*) m_complex;
        complexD* a_c = (complexD*) A.m_complex;

        zaxpy_(&size, (complexD*) &a, a_c, &inc, t_c, &inc);
    }
    else // if (!IsReal() && A.IsReal())
    {
        double* t_c = (double*) m_complex;
        double* a_r = A.m_real;

        int inc2 = 2;
        daxpy_(&size, a, a_r, &inc, t_c, &inc2);
    }

    return status;
}

//! Multiply two matrices so that (*this) = A * B
template<>
inline hwMathStatus hwTMatrix<double>::Mult(const hwTMatrix<double>& A, const hwTMatrix<double>& B)
{
    hwMathStatus status;

    if (this == &A)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    if (this == &B)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    // get dimensions info
    int m = A.m_nRows;
    int n = B.m_nCols;
    int k = A.m_nCols;

    if (B.m_nRows != k)
    {
        if (A.Size() == 0 && B.Size() == 1)
            status = Dimension(A.m_nRows, A.m_nCols, REAL);
        else if (A.Size() == 1 && B.Size() == 0)
            status = Dimension(B.m_nRows, B.m_nCols, REAL);
        else
            status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        return status;
    }

    // prepare for LAPACK function call
    char TRANSA = 'N';
    char TRANSB = 'N';
    int LDA = A.m_nRows;
    int LDB = B.m_nRows;
    int LDC = A.m_nRows;    // = C.m_nRows
    double ALPHA[2] = {1.0, 0.0};   // real or complex
    double BETA[2] = {0.0, 0.0};    // real or complex

    if (A.IsReal() && B.IsReal())
    {
        status = Dimension(m, n, REAL);

        if (!status.IsOk())
        {
            status.SetArg1(3);
            return status;
        }

        if (A.Size() != 0 && B.Size() != 0)
        {
            double* t_r = m_real;
            double* a_r = A.m_real;
            double* b_r = B.m_real;

            dgemm_(&TRANSA, &TRANSB, &m, &n, &k, ALPHA, a_r, &LDA, b_r, &LDB, BETA, t_r, &LDC);
        }
        else if (Size() != 0)
        {
            SetElements(0.0);
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        status = Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            status.SetArg1(3);
            return status;
        }

        complexD* t_c = (complexD*) m_complex;
        complexD* a_c = (complexD*) A.m_complex;
        complexD* b_c = (complexD*) B.m_complex;

        zgemm_(&TRANSA, &TRANSB, &m, &n, &k, (complexD*) &ALPHA,
               a_c, &LDA, b_c, &LDB, (complexD*) &BETA, t_c, &LDC);
    }
    else if (A.IsReal() && !B.IsReal())
    {
        status = Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            status.SetArg1(3);
            return status;
        }

        if (A.Size() != 0)
        {
            hwTMatrix<double> AC(A);
            status = AC.MakeComplex();

            if (!status.IsOk())
            {
                status.ResetArgs();
                return status;
            }

            complexD* t_c = (complexD*) m_complex;
            complexD* a_c = (complexD*) AC.m_complex;
            complexD* b_c = (complexD*) B.m_complex;

            zgemm_(&TRANSA, &TRANSB, &m, &n, &k, (complexD*) &ALPHA,
                   a_c, &LDA, b_c, &LDB, (complexD*) &BETA, t_c, &LDC);
        }
        else if (Size() != 0)
        {
            Dimension(m, n, REAL);
            SetElements(0.0);
        }
    }
    else // if (!A.IsReal() && B.IsReal())
    {
        status = Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            status.SetArg1(3);
            return status;
        }

        if (B.Size() != 0)
        {
            hwTMatrix<double> BC(B);
            status = BC.MakeComplex();

            if (!status.IsOk())
            {
                status.ResetArgs();
                return status;
            }

            complexD* t_c = (complexD*) m_complex;
            complexD* a_c = (complexD*) A.m_complex;
            complexD* b_c = (complexD*) BC.m_complex;

            zgemm_(&TRANSA, &TRANSB, &m, &n, &k, (complexD*) &ALPHA,
                   a_c, &LDA, b_c, &LDB, (complexD*) &BETA, t_c, &LDC);
        }
        else if (Size() != 0)
        {
            status = Dimension(m, n, REAL);
            SetElements(0.0);
        }
    }

    return status;
}

//! Multiply a matrix by a real number so that (*this) = A * x
template<>
inline hwMathStatus hwTMatrix<double>::Mult(const hwTMatrix<double>& A, double x)
{
    hwMathStatus status;

    if (this == &A)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    // get dimensions info
    int m = A.m_nRows;
    int n = A.m_nCols;
    int size = A.Size();

    // prepare for LAPACK function call
    int inc = 1;

    if (A.IsReal())
    {
        status = Dimension(m, n, REAL);

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }

        if (size != 0)
        {
            double* t_r = m_real;
            double* a_r = A.m_real;

            dcopy_(&size, a_r, &inc, t_r, &inc);
            dscal_(&size, (double*) &x, t_r, &inc);
        }
    }
    else
    {
        status = Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }

        complexD* t_c = (complexD*) m_complex;
        complexD* a_c = (complexD*) A.m_complex;

        int size2 = size<<1;
        zcopy_(&size, a_c, &inc, t_c, &inc);
        dscal_(&size2, (double*) &x, (double*) t_c, &inc);
    }

    return status;
}

//! Multiply a matrix by a complex number so that (*this) = A * z
template<>
inline hwMathStatus hwTMatrix<double>::Mult(const hwTMatrix<double>& A, const hwTComplex<double>& z)
{
    hwMathStatus status;    
    
    if (this == &A)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    // get dimensions info
    int m = A.m_nRows;
    int n = A.m_nCols;
    int size = A.Size();
    int inc = 1;

    status = Dimension(m, n, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(0);
        return status;
    }

    complexD* t_c = (complexD*) m_complex;

    if (A.IsReal())
    {
        if (size != 0)
        {
            hwTMatrix<double> AC(A);

            status = AC.MakeComplex();

            if (!status.IsOk())
            {
                status.ResetArgs();
                return status;
            }

            complexD* a_c = (complexD*) AC.m_complex;
            zcopy_(&size, a_c, &inc, t_c, &inc);
            zscal_(&size, (complexD*) &z, t_c, &inc);
        }
    }
    else
    {
        complexD* a_c = (complexD*) A.m_complex;
        zcopy_(&size, a_c, &inc, t_c, &inc);
        zscal_(&size, (complexD*) &z, t_c, &inc);
    }

    return status;
}

//! Multiply a matrix by a real number so that (*this)* = x
template<>
inline void hwTMatrix<double>::MultEquals(double x)
{
    int size = Size();
    int inc = 1;

    if (IsReal())
    {
        if (size != 0)
        {
            double* t_r = m_real;
            dscal_(&size, (double*) &x, t_r, &inc);
        }
    }
    else
    {
        double* t_c = (double*) m_complex;

        int size2 = size<<1;
        dscal_(&size2, (double*) &x, t_c, &inc);
    }
}

//! Multiply a matrix by a complex number so that (*this)* = z
template<>
inline hwMathStatus hwTMatrix<double>::MultEquals(const hwTComplex<double>& z)
{
    int size = Size();
    int inc = 1;
    hwMathStatus status;

    if (IsReal())
    {
        status = MakeComplex();

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }
    }

    if (size != 0)
    {
        complexD* t_c = (complexD*) m_complex;
        zscal_(&size, (complexD*) &z, t_c, &inc);
    }

    return status;
}

//! Divide a matrix by a matrix so that (*this) = A \ B, or A(*this) = B
template<>
inline hwMathStatus hwTMatrix<double>::DivideLeft(const hwTMatrix<double>& A, const hwTMatrix<double>& B)
{
    hwMathStatus status;

    if (A.m_nRows == A.m_nCols)
    {
        if (B.IsEmpty() && A.Size() == 1)   // treat A as a scalar
            return Dimension(B.m_nRows, B.m_nCols, hwTMatrix<double>::REAL);

        status = LSolve(A, B);
    }
    else
    {
        status = QRSolve(A, B);
    }

    if (status == HW_MATH_WARN_MTXNOTFULLRANK || status == HW_MATH_WARN_SINGMATRIXDIV)
    {
        // status = SVDSolve(A, B);     - this will fail if B contains NaN
        hwTMatrix<double> P;
        hwMathStatus status2 = P.Pinv(A);

        if (!status2.IsOk())
            return status2;

        (*this) = P * B;
    }

    return status;
}

//! Divide a matrix by a matrix so that (*this) = B / A, or (*this)A = B
template<>
inline hwMathStatus hwTMatrix<double>::DivideRight(const hwTMatrix<double>& B, const hwTMatrix<double>& A)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<double>& X = (*this);
    hwMathStatus status;
    hwTMatrix<double> AT;
    hwTMatrix<double> BT;
    hwTMatrix<double> XT;

    status = AT.Transpose(A);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = BT.Transpose(B);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = XT.Transpose(X);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = XT.DivideLeft(AT, BT);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 1)
            status.SetArg1(2);
        else if (status.GetArg1() == 2)
            status.SetArg1(1);

        if (status.GetArg2() == 2)
            status.SetArg2(1);
    }

    hwMathStatus status2 = X.Transpose(XT);

    if (!status2.IsOk())
    {
        status2.ResetArgs();
        return status2;
    }

    return status;
}

//! Implement the /= operator with a matrix argument
template<>
inline hwMathStatus hwTMatrix<double>::DivideEquals(const hwTMatrix<double>& A)
{
    hwTMatrix<double> temp(*this);
    hwMathStatus status;

    status = DivideRight(A, temp);

    if (!status.IsOk())
        MakeEmpty();

    return status;
}

//! Real determinant
template<>
inline hwMathStatus hwTMatrix<double>::Determinant(double& det) const
{
    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;
    int n = A.m_nCols;

    if (A.IsEmpty())
    {
        det = 1.0;
        return status;
    }

    if (A.m_nRows != n)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 0);

    hwTMatrix<double> B(A);
    double* a = B.m_real;

    // decompose the matrix
    int lda = n;
    int info;
    int* ipiv;

    try
    {
        ipiv = new int[n];
    }
    catch (std::bad_alloc&) 
    {
        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    dgetrf_(&n, &n, a, &lda, ipiv, &info); 

    if (info < 0)
        return status(HW_MATH_ERR_NOTCONVERGE);

    // compute determinant from trace
    int pivSign = 1;

    for (int i = 0; i < n-1; ++i)
    {
        if (i != ipiv[i]-1)
            pivSign = -pivSign;
    }

    delete [] ipiv;

    det = (double) pivSign * B(0, 0);

    for (int i = 1; i < n; ++i)
        det *= B(i, i);

    return status;
}

//! Complex determinant
template<>
inline hwMathStatus hwTMatrix<double>::Determinant(hwTComplex<double>& det) const
{
    const hwTMatrix<double>& A = (*this);
    hwMathStatus status;
    int n = A.m_nCols;

    if (A.IsEmpty())
    {
        det.Set(1.0, 0.0);
        return status;
    }

    if (A.m_nRows != n)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    if (A.IsReal())
    {
        double det_r;
        status = Determinant(det_r);

        if (!status.IsOk())
            return status;

        det.Set(det_r, 0.0);
    }

    hwTMatrix<double> B(A);
    complexD* a_c = (complexD*) B.m_complex;

    // decompose the matrix
    int lda = n;
    int info;
    int* ipiv;

    try
    {
        ipiv = new int[n];
    }
    catch (std::bad_alloc&) 
    {
        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    zgetrf_(&n, &n, a_c, &lda, ipiv, &info); 

    if (info < 0)
        return status(HW_MATH_ERR_NOTCONVERGE);

    // compute determinant from trace
    int i;
    int pivSign = 1;

    for (i = 0; i < n-1; ++i)
    {
        if (i != ipiv[i]-1)
            pivSign = -pivSign;
    }

    delete [] ipiv;

    det = B.z(0, 0) * (double) pivSign;

    for (i = 1; i < n; ++i)
        det *= B.z(i, i);

    return status;
}

//! Invert a matrix
template<>
inline hwMathStatus hwTMatrix<double>::Inverse(const hwTMatrix<double>& source)
{
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<double>& inv = (*this);

    // check dimensions
    int n = source.m_nCols;

    if (source.m_nRows != n)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (n == 0)
        return inv.Dimension(0, 0, REAL);

    // prepare for LAPACK function call
    inv = source;

    int lda = n;
    int info;
    int lwork = -1;
    int* ipiv = nullptr;

    try
    {
        ipiv = new int[n];
    }
    catch (std::bad_alloc&) 
    {
        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    if (source.IsReal())
    {
        // decompose the matrix
        double* a_r = inv.m_real;
        double* work = NULL;

        dgetrf_(&n, &n, a_r, &lda, ipiv, &info); 

        if (info)
        {
            if (info > 0)
            {
                inv.SetElements(std::numeric_limits<double>::infinity());
                return status(HW_MATH_WARN_SINGMATRIX, 1);
            }
            else
            {
                return status(HW_MATH_ERR_INVALIDINPUT, 1);
            }
        }
            
        // workspace query
        try
        {
            work = new double[1];
        }
        catch (std::bad_alloc&) 
        {
            if (ipiv)
                delete [] ipiv;

            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        dgetri_(&n, a_r, &lda, ipiv, work, &lwork, &info);

        lwork = static_cast<int>(work[0]);
        delete [] work;

        // back substitute
        try
        {
            work = new double[lwork];
        }
        catch (std::bad_alloc&) 
        {
            if (ipiv)
                delete [] ipiv;

            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        dgetri_(&n, a_r, &lda, ipiv, work, &lwork, &info);

        delete [] work;
    }
    else
    {
        // decompose the matrix
        complexD* a_c = (complexD*) inv.m_complex;
        complexD* work = NULL;

        zgetrf_(&n, &n, a_c, &lda, ipiv, &info); 

        if (info != 0)
        {
            inv.SetElements(std::numeric_limits<double>::infinity());
            return status(HW_MATH_WARN_SINGMATRIX, 1);
        }

        // workspace query
        try
        {
            work = new complexD[1];
        }
        catch (std::bad_alloc&) 
        {
            if (ipiv)
                delete [] ipiv;

            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        zgetri_(&n, a_c, &lda, ipiv, work, &lwork, &info);

        lwork = static_cast<int>(work[0].real());
        delete [] work;

        try
        {
            work = new complexD[lwork];
        }
        catch (std::bad_alloc&) 
        {
            if (ipiv)
                delete [] ipiv;

            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        // back substitute
        zgetri_(&n, a_c, &lda, ipiv, work, &lwork, &info);

        delete [] work;
    }

    delete [] ipiv;

    return status;
}

//! Matrix exponential (not a matrix of exponentials)
template<>
inline hwMathStatus hwTMatrix<double>::MatExp(const hwTMatrix<double>& power)
{
    hwMathStatus status;

    if (this == &power)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<double>& E = (*this);

    if (!power.IsSquare())
        return status(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (power.IsReal())
        status = E.Dimension(power.m_nRows, power.m_nCols, REAL);
    else
        status = E.Dimension(power.m_nRows, power.m_nCols, COMPLEX);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(2);
        else
            status.ResetArgs();

        return status;
    }

    //// scale A
    // get mag of largest row
    int n = power.m_nRows;
    double max = 0.0;
    double value;

    for (int i = 0; i < n; ++i)
    {
        value = 0.0;

        if (power.IsReal())
        {
            for (int j = 0; j < n; ++j)
                value += fabs(power(i, j));
        }
        else
        {
            for (int j = 0; j < n; ++j)
                value += power.z(i, j).Mag();
        }

        if (value > max)
            max = value;
    }

    // get scale and divide
    int exponent;
    double mantissa;
    hwTMatrix<double> B;

    mantissa = frexp(max, &exponent);

    if (exponent + 1 > 0)
        ++exponent;
    else
        exponent = 0;

    value = pow(2.0, (double) exponent);
    status = B.Divide(power, value);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    // Pade approximation for exp(A)
    hwTMatrix<double> X(B);
    hwTMatrix<double> I(n, n, REAL);
    hwTMatrix<double> tempM;
    hwTMatrix<double> D;
    double c = 0.5;

    I.Identity();
    tempM = B * c;
    E = I + tempM;
    D = I - tempM;

    int q = 6;
    bool p = true;

    for (int k = 2; k <= q; ++k)
    {
        c = c * (static_cast<long int>(q-k)+1.0) / (k*(2.0*q-k+1.0));
        X = B * X;
        tempM = X * c;
        E = E + tempM;

        if (p)
            D = D + tempM;
        else
            D = D - tempM;

        p = !p;
    }

    tempM = E;
    status = E.DivideLeft(D, tempM);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 2) { }
        else
            status.ResetArgs();

        return status;
    }

    // Undo scaling by repeated squaring
    for (int k = 1; k <= exponent; ++k)
        E = E * E;

    return status;
}

//! Solve tridigonal linear system AX=B with LU decomposition, where X = *this
template<>
inline hwMathStatus hwTMatrix<double>::LSolveT(const hwTMatrix<double>& D, const hwTMatrix<double>& DL,
                                               const hwTMatrix<double>& DU, const hwTMatrix<double>& B)
{
    hwMathStatus status;

    if (!D.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (!D.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!DL.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 2);

    if (!DL.IsVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    if (!DU.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 3);

    if (!DU.IsVector())
        return status(HW_MATH_ERR_VECTOR, 3);

    if (!B.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 4);

    if (DL.Size() != D.Size() - 1)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (DU.Size() != D.Size() - 1)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 3);

    int n = D.Size();
    int nrhs = B.m_nCols;
    int info;

    if (B.m_nRows != n)
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 4);

    *this = B;

    if (B.IsEmpty())
        return status;

    if (IsEmpty())
        return status(HW_MATH_ERR_ALLOCFAILED);

    hwTMatrix<double> DC(D);
    hwTMatrix<double> DLC(DL);
    hwTMatrix<double> DUC(DU);
    double* x = m_real;
    double* d = DC.m_real;
    double* dl = DLC.m_real;
    double* du = DUC.m_real;

    if (!x)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dgtsv_(&n, &nrhs, dl, d, du, x, &n, &info);

    if (info != 0)
        status(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Solve band linear system AX=B with LU decomposition, where X = *this
template<>
inline hwMathStatus hwTMatrix<double>::LSolveB(const hwTMatrix<double>& A, int kl, int ku,
                                               const hwTMatrix<double>& B)
{
    // To format A see SizeBandMatrix() and SetBandMatrixElem()
    hwMathStatus status;

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (kl < 0 || kl > A.m_nCols - 1)
        return status(HW_MATH_ERR_KLUD, 2);

    if (ku < 0 || ku > A.m_nCols - 1)
        return status(HW_MATH_ERR_KLUD, 3);

    int ldab = 2*kl+ku+1;

    if (A.m_nRows != ldab)
        return status(HW_MATH_ERR_ARRAYDIM, 2, 3);

    if (!B.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 4);

    if (B.m_nRows != A.m_nCols)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 4);

    int info;
    int n = A.m_nCols;
    int nrhs = B.m_nCols;
    hwTMatrix<int> ipiv(n, hwTMatrix<int>::REAL);     // vector of pivots
    hwTMatrix<double> decomp(A);

    *this = B;

    if (B.IsEmpty())
        return status;

    if (IsEmpty())
        return status(HW_MATH_ERR_ALLOCFAILED);

    int* p = ipiv.GetRealData();
    double* d = decomp.m_real;
    double* x = m_real;

    if (!d || !x)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dgbsv_(&n, &kl, &ku, &nrhs, d, &ldab, p, x, &n, &info);

    if (info != 0)
        status(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Solve symmetric positive definite linear system AX=B with Cholesky decomposition,
//! where X = *this
template<>
inline hwMathStatus hwTMatrix<double>::LSolveSPD(const hwTMatrix<double>& A, const hwTMatrix<double>& B)
{
    hwMathStatus status;

    if (this == &A)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    if (this == &B)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<double>& X = (*this);
    int m = A.m_nRows;
    int n = A.m_nCols;

    if (m < n)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    if (B.m_nRows != m)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (!B.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 2);

    hwTMatrix<double> AA(A);

    X = B;

    if (B.IsEmpty())
        return status;

    if (X.IsEmpty())
        return status(HW_MATH_ERR_ALLOCFAILED);

    double* a = AA.m_real;
    double* x = X.m_real;

    if (!a || !x)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // decompose the matrix and solve
    int lda = m;
    int ldb = m;
    int info;
    int bn = B.m_nCols;
    char UPLO = 'L';

    dposv_(&UPLO, &n, &bn, a, &lda, x, &ldb, &info);

    if (info != 0)
        status(HW_MATH_ERR_MTXNOTSPD, 1);

    return status;
}

//! Solve SPD tridigonal linear system AX=B with Cholesky decomposition,
//! where X = *this
template<>
inline hwMathStatus hwTMatrix<double>::LSolveSPDT(const hwTMatrix<double>& D, const hwTMatrix<double>& E,
                                                  const hwTMatrix<double>& B)
{
    hwMathStatus status;

    if (!D.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (!D.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!E.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 2);

    if (!E.IsVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    if (!B.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 3);

    int n = D.Size();

    if (E.Size() != n - 1)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (B.m_nRows != n)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 3);

    *this = B;

    if (B.IsEmpty())
        return status;

    if (IsEmpty())
        return status(HW_MATH_ERR_ALLOCFAILED);

    hwTMatrix<double> DD(D);
    hwTMatrix<double> EE(E);
    double* d = DD.m_real;
    double* e = EE.m_real;
    double* x = m_real;
    int info;
    int nrhs = B.m_nCols;

    if (!d || !e || !x)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dptsv_(&n, &nrhs, d, e, x, &n, &info);

    if (info != 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Solve SPD band linear system AX=B with Cholesky decomposition,
//! where X = *this
template<>
inline hwMathStatus hwTMatrix<double>::LSolveSPDB(const hwTMatrix<double>& A, int kd,
                                       const hwTMatrix<double>& B)
{
    // To format A see SizeSymBandMatrix() and SetSymBandMatrixElem()
    hwMathStatus status;

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (kd < 0 || kd > A.m_nCols - 1)
        return hwMathStatus(HW_MATH_ERR_KLUD, 2);

    int ldab = kd+1;

    if (A.m_nRows != ldab)
        return status(HW_MATH_ERR_ARRAYSIZE, 1);

    if (!B.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 3);

    int n = A.m_nCols;

    if (B.m_nRows != n)
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 3);

    char uplo = 'L';
    int info;
    int nrhs = B.m_nCols;
    hwTMatrix<double> decomp(A);

    *this = B;

    if (B.IsEmpty())
        return status;

    if (IsEmpty())
        return status(HW_MATH_ERR_ALLOCFAILED);

    double* d = decomp.m_real;
    double* x = m_real;

    if (!d || !x)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dpbsv_(&uplo, &n, &kd, &nrhs, d, &ldab, x, &n, &info);

    if (info != 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Solve symmetric indefinite linear system AX=B, where X = *this
template<>
inline hwMathStatus hwTMatrix<double>::LSolveSI(const hwTMatrix<double>& A, const hwTMatrix<double>& B)
{
    hwMathStatus status;

    // check dimensions
    int m = A.m_nRows;
    int n = A.m_nCols;

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (!B.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 2);

    if (A.m_nRows != n)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (B.m_nRows != m)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    // prepare for LAPACK function call
    char uplo = 'L';
    int info;
    int lwork = -1;
    int nrhs = B.m_nCols;
    hwTMatrix<int> ipiv(n, hwTMatrix<int>::REAL);     // vector of pivots
    hwTMatrix<double> decomp(A);

    *this = B;

    if (B.IsEmpty())
        return status;

    if (IsEmpty())
        return status(HW_MATH_ERR_ALLOCFAILED);

    double* d = decomp.m_real;
    int* p = ipiv.GetRealData();
    double* x = m_real;
    double* work;
    
    if (!d || !p || !x)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // workspace query
    try
    {
        work = new double[1];
    }
    catch (std::bad_alloc&) 
    {
        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    dsysv_(&uplo, &n, &nrhs, d, &n, p, x, &n, work, &lwork, &info);

    if (info != 0)
        return hwMathStatus(HW_MATH_ERR_DECOMPFAIL, 0);

    lwork = static_cast<int>(work[0]);
    delete [] work;

    try
    {
        work = new double[lwork];
    }
    catch (std::bad_alloc&) 
    {
        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    // solve the system
    dsysv_(&uplo, &n, &nrhs, d, &n, p, x, &n, work, &lwork, &info);

    delete [] work;

    return status;
}

//! Solve real linear system AX=B with SVD, where X = *this
template<>
inline hwMathStatus hwTMatrix<double>::SVDSolve(const hwTMatrix<double>& A, const hwTMatrix<double>& B)
{
    hwMathStatus status;

    if (this == &A)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    if (this == &B)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<double>& X = (*this);

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (!B.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 2);

    int m = A.m_nRows;
    int n = A.m_nCols;

    if (B.m_nRows != m)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    int nu = _min(m, n);
    int lda = m;
    int nrhs = B.m_nCols;
    int ldb = _max(m, n);
    int rank;
    int info;
    double rcond = 1.0e-12;
    double* s;

    try
    {
        s = new double[nu];
    }
    catch (std::bad_alloc&) 
    {
        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    hwTMatrix<double> Acopy(A);
    hwTMatrix<double> Bcopy(B);

    status = Bcopy.Resize(ldb, nrhs, true);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = X.Dimension(n, nrhs, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (A.IsEmpty() || X.IsEmpty())
        return status;

    double* a_r = Acopy.m_real;
    double* b_r = Bcopy.m_real;

    // workspace query
    double* work;
    int lwork = -1;

    try
    {
        work = new double[1];
    }
    catch (std::bad_alloc&) 
    {
        if (s)
            delete [] s;

        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    dgelss_(&m, &n, &nrhs, a_r, &lda, b_r, &ldb, s, &rcond, &rank,
            work, &lwork, &info);

    lwork = static_cast<int>(work[0]);
    delete [] work;

    try
    {
        work = new double[lwork];
    }
    catch (std::bad_alloc&) 
    {
        if (s)
            delete [] s;

        return status(HW_MATH_ERR_ALLOCFAILED);
    }

    // decompose the matrix and the solve system
    dgelss_(&m, &n, &nrhs, a_r, &lda, b_r, &ldb, s, &rcond, &rank,
            work, &lwork, &info);

    if (info != 0)
    {
        delete [] s;
        delete [] work;

        status(HW_MATH_ERR_NOTCONVERGE);
        return status;
    }

    int i, j;
    double thresh = 1.0e-9;
    double min = thresh * s[0];

    for (i = 0; i < nu; ++i)
    {
        if (s[i] < min)
            s[i] = 0.0;
    }

    if (s[nu-1] == 0.0)
    {
        if (m > n)
            status(HW_MATH_WARN_MATRIXDEPENDCOL);
        else
            status(HW_MATH_WARN_MTXNOTFULLRANK);
    }

    for (i = 0; i < n; ++i)
    {
        for (j = 0; j < nrhs; ++j)
            X(i, j) = Bcopy(i, j);
    }

    delete [] s;
    delete [] work;

    return status;
}

//! LU decomposition (PA=LU, where A = *this and P is a permutation matrix)
template<>
inline hwMathStatus hwTMatrix<double>::LU(hwTMatrix<int>& P, hwTMatrix<double>& L, hwTMatrix<double>& U) const
{
    if (IsReal())
        return RealLU(L, U, P);
    else
        return ComplexLU(L, U, P);
}

//! Cholesky decomposition
template<>
inline hwMathStatus hwTMatrix<double>::Csky(hwTMatrix<double>& T, bool upper) const
{
    hwMathStatus status;

    if (this == &T)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    int n = A.m_nCols;

    if (A.m_nRows != n)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (A.Size() != 0)
    {
        T = A;

        double* t = T.m_real;

        if (!t)
            return status(HW_MATH_ERR_ALLOCFAILED);

        // decompose the matrix
        int i, j;
        int lda = n;
        int info;
        char UPLO;

        if (upper)
            UPLO = 'U';
        else
            UPLO = 'L';

        dpotrf_(&UPLO, &n, t, &lda, &info);

        if (info != 0)
        {
            status(HW_MATH_ERR_MTXNOTSPD, 1);

            if (info > 0)
            {
                n = info - 1;
                hwMathStatus status2 = T.Resize(n, n);

                if (!status2.IsOk())
                {
                    status2.ResetArgs();
                    return status2;
                }
            }
            else
            {
                return status;
            }
        }

        if (upper)
        {
            for (i = 1; i < n; ++i)
            {
                for (j = 0; j < i; ++j)
                    T(i, j) = 0.0;
            }
        }
        else
        {
            for (j = 1; j < n; ++j)
            {
                for (i = 0; i < j; ++i)
                    T(i, j) = 0.0;
            }
        }
    }
    else
        status = T.Dimension(0, 0, REAL);

    return status;
}

//! Eigen decomposition
template<>
inline hwMathStatus hwTMatrix<double>::Eigen(bool balance, hwTMatrix<double>* V, hwTMatrix<double>& D) const
{
    hwMathStatus status;

    if (m_nRows != m_nCols)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (IsEmpty())
    {
        if (V)
            status = V->Dimension(0, 0, REAL);

        status = D.Dimension(0, 0, REAL);

        return status;
    }

    if (!IsFinite())
        return status(HW_MATH_ERR_NONFINITEDATA, 1);

    if (IsReal())
        status = EigenDecompReal(balance, V, D);
    else
        status = EigenDecompComplex(balance, V, D);

    return status;
}

//! Real symmetric Eigen decomposition
template<>
inline hwMathStatus hwTMatrix<double>::EigenDecompRealSymmetric(hwTMatrix<double>& V, hwTMatrix<double>& D) const
{
    hwMathStatus status;

    if (this == &V)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    if (this == &D)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    int n = A.m_nCols;

    if (!A.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 0);

    if (A.m_nRows != n)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    char JOBZ = 'V';
    char UPLO = 'U';

    int LDA = n;
    int INFO = 0;
    double WORKSZE = 0.0;
    int LWORK = -1;

    V = A;

    status = D.Dimension(n, REAL);

    if (!status.IsOk())
    {
        status.SetArg1(2);
        return status;
    }

    double* dVR = V.m_real;
    double* dD = D.m_real;

    if (!dVR)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // workspace query
    dsyev_(&JOBZ, &UPLO, &n, dVR, &LDA, dD, &WORKSZE, &LWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    LWORK = static_cast<int>(WORKSZE);
    hwTMatrix<double> WORK(LWORK, 1, REAL);

    double* dWORK = (double*) WORK.m_real;

    if (!dWORK)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // decompose the matrix
    dsyev_(&JOBZ, &UPLO, &n, dVR, &LDA, dD, dWORK, &LWORK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Complex Hermitian Eigen decomposition
template<>
inline hwMathStatus hwTMatrix<double>::EigenDecompComplexHermitian(hwTMatrix<double>& V, hwTMatrix<double>& D) const
{
    hwMathStatus status;

    if (this == &V)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    if (this == &D)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    const hwTMatrix<double>& A = (*this);
    int n = A.m_nCols;

    if (A.m_nRows != n)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    char JOBZ = 'V';
    char UPLO = 'U';
    int LDA = n;
    int LWORK = -1;
    int INFO = 0;

    hwTMatrix<double> WORKSZE(1, 1, COMPLEX);
    hwTMatrix<double> RWORK(_max(1,3*A.m_nCols-2), 1, REAL);

    V = A;

    if (V.IsReal())
        status = V.MakeComplex();

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = D.Dimension(n, 1, REAL);

    if (!status.IsOk())
    {
        status.SetArg1(3);
        return status;
    }

    complexD* dV = (complexD*) V.m_complex;
    double* dD = D.m_real;
    complexD* dWS = (complexD*) WORKSZE.m_complex;
    double* dRWK = RWORK.m_real;

    if (!dV || !dD || !dWS || !dRWK)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // workspace query
    zheev_(&JOBZ, &UPLO, &n, dV, &LDA, dD, dWS, &LWORK, dRWK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    LWORK = static_cast<int>(WORKSZE.z(0).Real());

    hwTMatrix<double> WORK(LWORK, 1, COMPLEX);

    complexD* dW = (complexD*) WORK.m_complex;

    if (!dW)
        return status(HW_MATH_ERR_ALLOCFAILED);

    // decompose the matrix
    zheev_(&JOBZ, &UPLO, &n, dV, &LDA, dD, dW, &LWORK, dRWK, &INFO);

    if (INFO < 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Balance matrix
template<>
inline hwMathStatus hwTMatrix<double>::Balance(bool noperm, hwTMatrix<double>& S,
                                       hwTMatrix<double>& P, hwTMatrix<double>& B) const
{
    hwMathStatus status;

    if (IsEmpty())
    {
        status = S.Dimension(0, 0, REAL);
        status = P.Dimension(0, 0, REAL);
        status = B.Dimension(0, 0, REAL);

        return status;
    }

    // Square matrix required
    if (m_nRows != m_nCols)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    // Finite matrix required
    if (!IsFinite())
        return status(HW_MATH_ERR_NONFINITEDATA, 0);

    if (IsReal())
        status = BalanceReal(noperm, S, P, B);
    else
        status = BalanceComplex(noperm, S, P, B);

    return status;
}

//! Eigen decomposition - Symmetric/Hermitian
template<>
inline hwMathStatus hwTMatrix<double>::EigenSH(hwTMatrix<double>* V, hwTMatrix<double>& D) const
{
    hwMathStatus status;

    if (m_nRows != m_nCols)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    if (IsEmpty())
    {
        if (V)
            status = V->Dimension(0, 0, REAL);

        status = D.Dimension(0, 0, REAL);

        return status;
    }

    if (!IsFinite())
        return status(HW_MATH_ERR_NONFINITEDATA, 0);

    if (IsReal())
        status = EigenDecompRealSymmetric(*V, D);
    else
        status = EigenDecompComplexHermitian(*V, D);

    return status;
}

//! Eigen decomposition of a symmetric tridigonal matrix
template<>
inline hwMathStatus hwTMatrix<double>::EigenST(const hwTMatrix<double>& D, const hwTMatrix<double>& E,
                                               hwTMatrix<double>& W, hwTMatrix<double>& Z)
{
    hwMathStatus status;

    if (!D.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (!D.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!E.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 2);

    if (!E.IsVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    int n = D.Size();

    if (E.Size() != n - 1)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    status = Z.Dimension(n, n, REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(4);
        else
            status.ResetArgs();

        return status;
    }

    hwTMatrix<double> EE(E);
    W = D;

    if (D.IsEmpty())
        return status;

    char jobz = 'V';
    int info;
    hwTMatrix<double> work(2 * D.Size() - 2, REAL);

    double* w = W.m_real;
    double* z = Z.m_real;
    double* wrk = work.m_real;
    double* e = EE.m_real;

    if (!w || !z || !wrk || !e)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dstev_(&jobz, &n, w, e, z, &n, wrk, &info);

    if (info != 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Eigen decomposition of a symmetric band matrix
template<>
inline hwMathStatus hwTMatrix<double>::EigenSB(int kd, hwTMatrix<double>& W, hwTMatrix<double>& Z) const
{
    // To format A see SizeSymBandMatrix() and SetSymBandMatrixElem()
    hwMathStatus status;

    if (!IsReal())
        return status(HW_MATH_ERR_COMPLEX, 0);

    if (kd < 0 || kd > m_nCols - 1)
        return status(HW_MATH_ERR_KLUD, 1);

    int n = m_nCols;
    int ldab = kd + 1;

    if (m_nRows != ldab)
        return status(HW_MATH_ERR_ARRAYDIM, 0, 1);

    status = W.Dimension(n, REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(2);
        else
            status.ResetArgs();

        return status;
    }

    status = Z.Dimension(n, n, REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(3);
        else
            status.ResetArgs();

        return status;
    }

    char jobz = 'V';
    char uplo = 'L';
    int info;
    hwTMatrix<double> work(3*n-2, REAL);
    hwTMatrix<double> AA(*this);
    double* a = AA.m_real;
    double* w = W.m_real;
    double* z = Z.m_real;
    double* wrk = work.m_real;

    if (!a || !w || !z || !wrk)
        return status(HW_MATH_ERR_ALLOCFAILED);

    dsbev_(&jobz, &uplo, &n, &kd, a, &ldab, w, z, &n, wrk, &info);

    if (info != 0)
        return status(HW_MATH_ERR_DECOMPFAIL);

    return status;
}

//! Generalized Eigen decomposition
template<>
inline hwMathStatus hwTMatrix<double>::Eigen(const hwTMatrix<double>& A, const hwTMatrix<double>& B,
                    hwTMatrix<double>& V, hwTMatrix<double>& D) const
{
    hwMathStatus status;

    if (A.m_nRows != A.m_nCols)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (B.m_nRows != B.m_nCols)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 2);

    if (A.m_nRows != B.m_nRows)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (A.IsEmpty())
    {
        status = V.Dimension(0, 0, REAL);
        status = D.Dimension(0, 0, REAL);

        return status;
    }

    if (!A.IsFinite())
        return status(HW_MATH_ERR_NONFINITEDATA, 1);

    if (!B.IsFinite())
        return status(HW_MATH_ERR_NONFINITEDATA, 2);

    if (A.IsReal() && B.IsReal())
    {
        if (A.IsSymmetric() && B.IsSymmetric())
            status = GeneralizedEigenDecompRealSymmetric(A, B, V, D);
        else
            status = GeneralizedEigenDecompReal(A, B, V, D);
    }
    else
    {
        if (A.IsHermitian() && B.IsHermitian())
            status = GeneralizedEigenDecompComplexHermitian(A, B, V, D);
        else
            status = GeneralizedEigenDecompComplex(A, B, V, D);
    }

    return status;
}

//! QR decomposition (A=QR, where A = *this)
template<>
inline hwMathStatus hwTMatrix<double>::QR(hwTMatrix<double>& Q, hwTMatrix<double>& R) const
{
    if (IsReal())
        return RealQR(Q, R);
    else
        return ComplexQR(Q, R);
}

enum SVDtype
{
    SVD_STD_FLAG,   // full SVD
    SVD_ECON_FLAG,  // economy size with S returned as a matrix
    SVD_ECON_0_FLAG // economy size with S returned as a vector
};

//! Singular Value Decomposition
template<>
inline hwMathStatus hwTMatrix<double>::SVD(int flagSvd, hwTMatrix<double>*U, hwTMatrix<double>& S, hwTMatrix<double>* V) const
{
    SVDtype type;
    hwMathStatus status;

    switch (flagSvd)
    {
    case 0:
        type = SVD_STD_FLAG;
        break;
    case 1:
        type = SVD_ECON_FLAG;
        break;
    case 2:
        type = SVD_ECON_0_FLAG;
        break;
    default:
        return status(HW_MATH_ERR_INVALIDINPUT, 1);
        break;
    }

    if (!IsFinite())
        return hwMathStatus(HW_MATH_ERR_NONFINITEDATA, 0);

    if (IsReal())
        status = RealSVD(type, U, S, V);
    else
        status = ComplexSVD(type, U, S, V);

    return status;
}

// Schur decomposition
template<>
inline hwMathStatus hwTMatrix<double>::Schur(bool real, hwTMatrix<double>& U, hwTMatrix<double>& T) const
{
    hwMathStatus status;

    if (!IsSquare())
        return status(HW_MATH_ERR_MTXNOTSQUARE, 0);

    if (IsEmpty())
    {
        // empty matrix
        // nothing to do
        // U = []
        // T = []
        return status;
    }   

    if (IsReal())
    {
        if (real)
            status = RealSchur(U, T);
        else
            status = ComplexSchur(U, T);
    }
    else
        status = ComplexSchur(U, T);        

    return status;
}

//! Condition number
template<>
inline hwMathStatus hwTMatrix<double>::Cond(double& condNum) const
{
    hwMathStatus status;

    if (IsEmpty())
    {
        condNum = 0.0;
        return status;
    }

    hwTMatrix<double> S;

    status = SVD(0, NULL, S, NULL);

    if (!status.IsOk())
        return status;

    int nu = S.Size();

    if (nu)
    {
        double largest = S(0);
        double smallest = S(nu-1);

        if ((largest == 0.0) || (smallest == 0.0))
            condNum = std::numeric_limits<double>::infinity();
        else
            condNum = largest / smallest;
    }
    else
    {
        condNum = 0.0;
    }

    return status;
}

//! Rank
template<>
inline hwMathStatus hwTMatrix<double>::Rank(int& rank) const
{
    hwMathStatus status;

    if (IsEmpty())
    {
        rank = 0;
        return status;
    }

    hwTMatrix<double> S;

    status = SVD(0, NULL, S, NULL);

    if (!status.IsOk())
        return status;

    int m = m_nRows;
    int n = m_nCols;
    int nu = _min(m, n);
    double tol = _max(m, n) * MachPrecision(S(0));

    rank = 0;

    for (int i = 0; i < nu; ++i)
    {
        if (S(i) > tol)
            ++rank;
        else
            break;
    }

    return status;
}

//! Rank with a tolerance interval
template<>
inline hwMathStatus hwTMatrix<double>::Rank(double tol, int& rank) const
{
    hwMathStatus status;
    hwTMatrix<double> S;

    status = SVD(0, NULL, S, NULL);

    if (!status.IsOk())
        return status;

    int m = m_nRows;
    int n = m_nCols;
    int nu = _min(m, n);

    rank = 0;

    for (int i = 0; i < nu; ++i)
    {
        if (S(i) > tol)
            ++rank;
        else
            break;
    }

    return status;
}

//! L2 vector norm
template<>
inline hwMathStatus hwTMatrix<double>::L2Norm(double& norm) const
{
    hwMathStatus status;

    if (IsReal())
    {
        if (!m_real)
        {
            norm = 0.0;
            return status;
        }

        int n = Size();
        int inc = 1;

        norm = dnrm2_(&n, m_real, &inc);
    }
    else
    {
        status = L2NormSq(norm);

        if (status.IsOk())
            norm = sqrt(norm);
    }

    return status;
}

//! Matrix p norm
template<>
inline hwMathStatus hwTMatrix<double>::Norm(double& norm, int p) const
{
    hwMathStatus status;

    if (p == 1)
    {
        char NORM;
        double* Work = NULL;
        int m, n;

        if (IsEmpty())
        {
            norm = 0.0;
            return status;
        }

        if (IsVector())
        {
            NORM = '1';
            m = Size();
            n = 1;
        }
        else    //matrix
        {
            NORM = '1';
            m = m_nRows;
            n = m_nCols;
        }

        if (IsReal())
            norm = dlange_(&NORM, &m, &n, m_real, &m, Work);
        else
            norm = zlange_(&NORM, &m, &n, (complexD*) m_complex, &m, Work);
    }
    else if (p == 2)
    {
        if (IsVector())
        {
            status = L2Norm(norm);
        }
        else    //matrix
        {
            if (IsEmpty())
            {
                norm = 0.0;
                return status;
            }

            hwTMatrix<double> S;

            status = SVD(0, NULL, S, NULL);

            if (!status.IsOk())
                return status;

             norm = S(0);
        }
    }
    else if (p > 2)
    {
        if (IsEmpty())
        {
            norm = 0.0;
            return status;
        }

        if (!IsVector())
        {
            status(HW_MATH_ERR_INVALIDINPUT, 2);
            return status;
        }
            
        int count = Size();
        norm = 0.0;

        if (IsReal())
        {
            double* a_r = m_real;

            while (count--)
                norm += CustomPow(fabs(*a_r++), (double) p);
        }
        else
        {
            hwComplex* a_c = m_complex;

            while (count--)
                norm += CustomPow((a_c++)->Mag(), (double) p);
        }

        norm = CustomPow(norm, 1.0 / p);
    }
    else
    {
        status(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    return status;
}

//! Matrix norm - Frobenius or inf
template<>
inline hwMathStatus hwTMatrix<double>::Norm(double& norm, const char* type) const
{
    hwMathStatus status;

    if (!type)
        return status(HW_MATH_ERR_NOTSTRING, 2);

    if (!strcmp(type, "-inf"))
    {
        if (!IsVector())
            return status(HW_MATH_ERR_INVALIDINPUT, 2);

        int count = Size();

        if (IsReal())
        {
            if (!m_real)
            {
                norm = 0.0;
                return status;
            }

            double* a_r = m_real;

            norm = fabs(*a_r++);
            count--;

            while (count--)
            {
                if (fabs(*a_r) < norm)
                    norm = fabs(*a_r);

                ++a_r;
            }
        }
        else
        {
            hwComplex* a_c = m_complex;

            norm = (a_c++)->Mag();
            count--;

            while (count--)
            {
                double mag = a_c->Mag();

                if (mag < norm)
                    norm = mag;

                ++a_c;
            }
        }

        return status;
    }

    int m, n;
    char NORM;
    double* Work = NULL;

    if (IsVector())
    {
        m = Size();
        n = 1;
    }
    else
    {
        m = m_nRows;
        n = m_nCols;
    }

    if (!strcmp(type, "inf"))
    {
        NORM = 'I';

        try
        {
            Work = new double[m];
        }
        catch (std::bad_alloc&) 
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }
    }
    else if (!strcmp(type, "fro"))
        NORM = 'F';
    else
        return status(HW_MATH_ERR_INVALIDINPUT, 2);

    if (IsEmpty())
    {
        norm = 0.0;
        return status;
    }

    if (IsReal())
        norm = dlange_(&NORM, &m, &n, m_real, &m, Work);
    else    // complex
        norm = zlange_(&NORM, &m, &n, (complexD*) m_complex, &m, Work);

    if (Work)
        delete [] Work;

    return status;
}

//! Normalize the object vector
template<>
inline hwMathStatus hwTMatrix<double>::Normalize()
{
    int count = Size();
    double norm;
    hwMathStatus status;

    if (!IsVector())
        return status(HW_MATH_ERR_VECTOR, 0);

    status = Norm(norm, 2);

    if (!status.IsOk())
        return status;

    double* v = m_real;

    while (count--)
        v[count] /= norm;

    return hwMathStatus();
}

//! Dot product of two real vectors
template<>
inline hwMathStatus hwTMatrix<double>::Dot(const hwTMatrix<double>& A, const hwTMatrix<double>& B, double& dot)
{
    if (A.m_complex)
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);

    if (!A.IsEmptyOrVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

    if (B.m_complex)
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);

    if (!B.IsEmptyOrVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);

    int count = A.Size();

    if (B.Size() != count)
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    int inc = 1;
    double* a_r = A.m_real;
    double* b_r = B.m_real;

    if (!a_r)
    {
        dot = 0.0;
        return hwMathStatus();
    }

    dot = ddot_(&count, a_r, &inc, b_r, &inc);

    return hwMathStatus();
}

//! Dot product of two vectors, either real or complex
template<>
inline hwMathStatus hwTMatrix<double>::Dot(const hwTMatrix<double>& A, const hwTMatrix<double>& B, hwTComplex<double>& dot)
{
    hwMathStatus status;

    if (!A.IsEmptyOrVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!B.IsEmptyOrVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    if (A.m_real && B.m_real)
    {
        status = Dot(A, B, dot.Real());
        dot.Imag() = 0.0;
        return status;
    }

    int count = A.Size();

    if (B.Size() != count)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (A.m_real && B.m_complex)
    {
        hwTMatrix<double> temp;
        
        status = temp.PackComplex(A);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        status = Dot(temp, B, dot);
    }
    else if (A.m_complex && B.m_real)
    {
        hwTMatrix<double> temp;
        
        status = temp.PackComplex(B);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        status = Dot(A, temp, dot);
    }
    else // (A.m_complex && B.m_complex)
    {
        int inc = 1;
        complexD* a_c = (complexD*) A.m_complex;
        complexD* b_c = (complexD*) B.m_complex;

        if (!a_c)
        {
            dot = 0.0;
            return status;
        }

        zdotc_((complexD*) (&dot), &count, a_c, &inc, b_c, &inc);
    }

    return status;
}

//! Raise the object matrix to an integer power
template<>
inline hwMathStatus hwTMatrix<double>::Power(const hwTMatrix<double>& base, int power)
{
    if (this == &base)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    // This algorithm uses a divide and conquer stategy, but then performs
    // the calculations in reverse order. To compute M^p find n such that
    //      M^p = M^(2n)     if p is even
    //      M^p = M^(2n+1)   if p is odd
    // Now set p=n and repeat. The result is a sequence of values of n that
    // ends with 1. The sequence is used in reverse order to accumulate
    // multiples of M using as many squaring options as possible to minimize
    // the number of matrix multiplies.
    hwTMatrix<double>& result = (*this);
    hwMathStatus status;

    if (base.m_nRows != base.m_nCols)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (power < 0)
        return status(HW_MATH_WARN_NONNONNEGINT, 2);

    if (power == 0)
    {
        status = result.Dimension(base.m_nRows, base.m_nCols, base.Type());

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        result.Identity();
        return status;
    }

    if (power == 1)
    {
        result = base;
        return status;
    }

    // accumulate the total
    int pow_cur2;                   // 2 * the accumulated power
    int pow_tgt;                    // target power for each iteration
    hwTMatrix<double> A(base.m_nRows, base.m_nRows, base.Type()); // accumulator

    A = base;
    pow_cur2 = 2;

    while (true)
    {
        // find pow_tgt
        pow_tgt = power;

        while (pow_tgt > pow_cur2 + 1)
            pow_tgt >>= 1;

        // update result
        result = A * A;

        if (pow_tgt != pow_cur2)    // if pow_tgt is odd
        {
            A = result;
            result = A * base;
        }

        // update accumulator unless done
        if (pow_tgt < power)
        {
            A = result;
            pow_cur2 = pow_tgt << 1;
        }
        else
            break;
    }

    return status;
}

//! Raise the object matrix to a possibly non integer power
//! A symmetric matrix is assumed
template<>
inline hwMathStatus hwTMatrix<double>::Power(const hwTMatrix<double>& base, double power)
{
    if (this == &base)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;
    hwTMatrix<double>& result = (*this);

    if (IsInteger(power, 1.0e-10) != HW_MATH_ERR_NONINTEGER)
    {
        if (power >= 0 && power < 9)
        {
            return result.Power(base, static_cast<int>(power));
        }
        else if (power > -9 && power < 0)
        {
            hwTMatrix<double> C(base.m_nRows, base.m_nRows, REAL); // accumulator
            status = C.Inverse(base);

            if (!status.IsOk())
            {
                if (status.GetArg1() != 1)
                    status.ResetArgs();

                if (!status.IsWarning())
                    return status;
            }

            hwMathStatus status2 = result.Power(C, static_cast<int>(-power));

            if (!status2.IsOk())
                return status2;
            else
                return status;
        }
    }

    hwTMatrix<double> W;
    hwTMatrix<double> WP;
    hwTMatrix<double> V;
    hwTMatrix<double> D;

    if (base.IsHermitian())
        status = base.EigenSH(&V, W);
    else
        status = base.Eigen(true, &V, W);

    if (!status.IsOk())
    {
        if (IsInteger(power, 1.0e-10) != HW_MATH_ERR_NONINTEGER)     // decomp method failed, see if int algorithm applies
            return result.Power(base, static_cast<int>(power));

        else
        {
            status.ResetArgs();
            return status;
        }
    }

    status = WP.PowerByElems(W, power);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = D.Diag(WP, 0);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    // construct result matrix from Eigen decomposition
    status = result.DivideRight(V*D, V);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    return status;
}

//! Raise the object matrix (symmetric) to a complex power
template<>
inline hwMathStatus hwTMatrix<double>::Power(const hwTMatrix<double>& base, const hwTComplex<double>& power)
{
    if (this == &base)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    if (power.Imag() == 0.0)
        return Power(base, power.Real());

    hwMathStatus status;
    hwTMatrix<double>& result = (*this);
    hwTMatrix<double> W;
    hwTMatrix<double> WP;
    hwTMatrix<double> V;
    hwTMatrix<double> D;

    if (base.IsHermitian())
        status = base.EigenSH(&V, W);
    else
        status = base.Eigen(true, &V, W);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = WP.PowerByElems(W, power);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = D.Diag(WP, 0);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    // construct result matrix from Eigen decomposition
    status = result.DivideRight(V*D, V);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    return status;
}

//! 2D convolution of two matrices
template<>
inline hwMathStatus hwTMatrix<double>::Conv2D(const hwTMatrix<double>& X,
                                              const hwTMatrix<double>& Y,
                                              int row1, int col1, int row2, int col2)
{
    if (this == &X)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &Y)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    // 2D linear convolution
    // X is reversed and passed over Y
    hwMathStatus status;

    if (X.IsEmpty())
        return Dimension(0, 0, REAL);

    if (X.m_complex)
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 1);

    if (Y.IsEmpty())
        return Dimension(0, 0, REAL);

    if (Y.m_complex)
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 2);

    int xm = X.m_nRows;
    int xn = X.m_nCols;
    int ym = Y.m_nRows;
    int yn = Y.m_nCols;
    int cm = xm + ym - 1;       // # rows in conv
    int cn = xn + yn - 1;       // # cols in conv

    if (row1 < 0)
        return status(HW_MATH_ERR_INVALIDINDEX, 3);

    if (row2 > cm - 1)
        return status(HW_MATH_ERR_INVALIDINDEX, 5);

    if (row2 == -1)
        row2 = cm - 1;

    if (col1 < 0)
        return status(HW_MATH_ERR_INVALIDINDEX, 4);

    if (col2 > cn - 1)
        return status(HW_MATH_ERR_INVALIDINDEX, 6);

    if (col2 == -1)
        col2 = cn - 1;

    status = Dimension(row2 - row1 + 1, col2 - col1 + 1, REAL);

    if (!status.IsOk())
    {
        if (cm < 1 || cn < 1)
            return hwMathStatus();

        if (status.GetArg1() != 0)
            status.ResetArgs();

        return status;
    }

    // flip the smaller of X or Y in both dimensions
    int incS = 1;
    int incD = -1;
    int count = xm * xn;
    int count2 = ym * yn;
    hwTMatrix<double> flip;
    bool flipX;

    if (count2 > count) // flip X
    {
        flipX = true;
        status = flip.Dimension(xm, xn, REAL);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        dcopy_(&count, const_cast<double*> (X.m_real), &incS, flip.m_real, &incD);
    }
    else    // flip Y
    {
        flipX = false;
        count = count2;
        status = flip.Dimension(ym, yn, REAL);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        dcopy_(&count, const_cast<double*> (Y.m_real), &incS, flip.m_real, &incD);
    }

    // find each conv(i, j)
    int start_xi;                       // X column indexing
    int start_xj;                       // X row indexing
    int start_yi, stop_yi;              // Y column window
    int start_yj, stop_yj;              // Y row window
    int partial_i = _min(xm, ym) - 1;   // size of partial col overlaps
    int partial_j = _min(xn, yn) - 1;   // size of partial row overlaps
    double value;

    for (int j = col1; j <= col2; ++j)
    {
        if (j < partial_j)
        {
            start_yj = 0;
            stop_yj = j + 1;
            start_xj = j;
        }
        else if (j < cn - partial_j)
        {
            if (xn < yn)
            {
                start_yj = j - xn + 1;
                stop_yj = j + 1;
                start_xj = xn - 1;
            }
            else
            {
                start_yj = 0;
                stop_yj = yn;
                start_xj = j;
            }
        }
        else
        {
            start_yj = j - xn + 1;
            stop_yj = yn;
            start_xj = xn - 1;
        }

        for (int i = row1; i <= row2; ++i)
        {
            if (i < partial_i)
            {
                start_yi = 0;
                stop_yi = i + 1;
                start_xi = i;
            }
            else if (i < cm - partial_i)
            {
                if (xm < ym)
                {
                    start_yi = i - xm + 1;
                    stop_yi = i + 1;
                    start_xi = xm - 1;
                }
                else
                {
                    start_yi = 0;
                    stop_yi = ym;
                    start_xi = i;
                }
            }
            else
            {
                start_yi = i - xm + 1;
                stop_yi = ym;
                start_xi = xm - 1;
            }

            // multiply the 2D overlap of Y and the reversed, shifted X
            int incx  = 1;
            int incy  = 1;
            value     = 0.0;
            count     = stop_yi - start_yi;

            if (flipX)
            {
                start_xi = xm - 1 - start_xi;
                double* xp = const_cast<double*> (&(flip(start_xi, xn - 1 - start_xj)));
                double* yp = const_cast<double*> (&(Y(start_yi, start_yj)));

                for (int jj = start_yj; jj < stop_yj; ++jj)
                {
                    value += ddot_(&count, xp, &incx, yp, &incy);
                    xp    += xm;
                    yp    += ym;
                }
            }
            else
            {
                start_xi = start_xi - (count - 1);
                start_yi = ym - 1 - start_yi - (count - 1);
                double* xp = const_cast<double*> (&(X(start_xi, start_xj)));
                double* yp = const_cast<double*> (&(flip(start_yi, yn - 1 - start_yj)));

                for (int jj = start_yj; jj < stop_yj; ++jj)
                {
                    value += ddot_(&count, xp, &incx, yp, &incy);
                    xp    -= xm;
                    yp    -= ym;
                }
            }

            (*this)(i - row1, j - col1) = value;
        }
    }

    return status;
}

//*******************************************************************
//                  hwTMatrix<double> band matrix
//*******************************************************************
//! Dimension a band matrix with CDS storage
template<>
inline hwMathStatus hwTMatrix<double>::DimensionBandMatrix(int n, int kl, int ku, bool initZero)
{
    // prepare A to store a CDS matrix by column
    int ldab = 2*kl+ku+1;
    hwMathStatus status;

    status = Dimension(ldab, n, REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(4);
        else
            status.ResetArgs();
    }

    if (initZero)
        SetElements(0.0);

    // check last for API reasons
    if (kl < 0 || kl > n - 1)
        return status(HW_MATH_ERR_KLUD, 2);

    if (ku < 0 || ku > n - 1)
        return status(HW_MATH_ERR_KLUD, 3);

    return status;
}

//! Dimension a symmetric band matrix
template<>
inline hwMathStatus hwTMatrix<double>::DimensionSymBandMatrix(int n, int kd, bool initZero)
{
    //! Dimension a symmetric band matrix with CDS storage
    int ldab = kd + 1;
    hwMathStatus status;

    status = Dimension(ldab, n, REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(3);
        else
            status.ResetArgs();
    }

    if (initZero)
        SetElements(0.0);

    // check last for API reasons
    if (kd < 0 || kd > n - 1)
        return status(HW_MATH_ERR_KLUD, 2);

    return status;}

//! Set a band matrix element
template<>
inline hwMathStatus hwTMatrix<double>::SetBandMatrixElem(int i, int j, double value, int kl, int ku)
{
    int row = kl+ku+i-j;
    int col = j;

    if (kl < 0 || kl > m_nRows - 1)
        return hwMathStatus(HW_MATH_ERR_KLUD, 5);

    if (ku < 0 || ku > m_nRows - 1)
        return hwMathStatus(HW_MATH_ERR_KLUD, 6);

    this->operator()(row, col) = value;

    return hwMathStatus();
}

//! Set a symmetric band matrix element
template<>
inline hwMathStatus hwTMatrix<double>::SetSymBandMatrixElem(int i, int j, double value, int kd)
{
    int row = i-j;
    int col = j;

    if (kd < 0 || kd > m_nRows - 1)
        return hwMathStatus(HW_MATH_ERR_KLUD, 5);

    this->operator()(row, col) = value;

    return hwMathStatus();
}
// ********** END SPECIALIZED MATRIX SECTION **********
