/**
* @file hwTMatrix.cc
* @date March 2012
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
//  hwTMatrix template function implementation file
//
//:---------------------------------------------------------------------------

#include <string.h>
#include <limits>
#include <iostream>
#include <sstream>
#include <string>

#include <Globals.h>
#include <hwMathStatus.h>
#include <tmpl/hwTComplex.h>

#ifdef max
   #undef max
#endif

// ****************************************************
//       Construction / Destruction / Assignment
// ****************************************************

//! Construct empty matrix
template<typename T1, typename T2>
hwTMatrix<T1, T2>::hwTMatrix()
    : m_refCount(1)
{
    m_bits.ownData = 0;
    MakeEmpty();
}

//! Construct vector
template<typename T1, typename T2>
hwTMatrix<T1, T2>::hwTMatrix(int size, DataType dataType)
                             : m_refCount(1), m_capacity(0)
{
    if (size)
    {
        m_nCols = 1;
        m_nRows = size;
    }
    else
    {
        m_nCols = 0;
        m_nRows = 0;
    }

    m_bits.ownData = 1;
    m_real_memory = nullptr;
    m_complex_memory = nullptr;
    m_real = nullptr;
    m_complex = nullptr;

    if (dataType == REAL)
        m_bits.realData = 1;
    else    // dataType == COMPLEX
        m_bits.realData = 0;

    Allocate(dataType);
}

//! Construct matrix
template<typename T1, typename T2>
hwTMatrix<T1, T2>::hwTMatrix(int m, int n, DataType dataType)
                             : m_refCount(1), m_capacity(0)
{
    m_nCols = n;
    m_nRows = m;
    m_bits.ownData = 1;
    m_real_memory = nullptr;
    m_complex_memory = nullptr;
    m_real = nullptr;
    m_complex = nullptr;

    if (dataType == REAL)
        m_bits.realData = 1;
    else    // dataType == COMPLEX
        m_bits.realData = 0;

    Allocate(dataType);
}

//! Construct vector with external data
template<typename T1, typename T2>
hwTMatrix<T1, T2>::hwTMatrix(int size, void* data, DataType dataType)
                             : m_refCount(1), m_capacity(0)
{
    if (size)
    {
        m_nCols = 1;
        m_nRows = size;
    }
    else
    {
        m_nCols = 0;
        m_nRows = 0;
    }

    m_real_memory = nullptr;
    m_complex_memory = nullptr;
    m_bits.ownData = 0;

    if (dataType == REAL)
    {
        m_bits.realData = 1;
        m_real = reinterpret_cast<T1*> (data);
        m_complex = nullptr;
    }
    else    // dataType == COMPLEX
    {
        m_bits.realData = 0;
        m_real = nullptr;
        m_complex = reinterpret_cast<T2*> (data);
    }
}

//! Construct matrix with external data
template<typename T1, typename T2>
hwTMatrix<T1, T2>::hwTMatrix(int m, int n, void* data,
                             DataType dataType)
                             : m_refCount(1), m_capacity(0)
{
    m_nCols = n;
    m_nRows = m;
    m_real_memory = nullptr;
    m_complex_memory = nullptr;
    m_bits.ownData = 0;

    if (dataType == REAL)
    {
        m_bits.realData = 1;
        m_real = reinterpret_cast<T1*> (data);
        m_complex = nullptr;
    }
    else    // dataType == COMPLEX
    {
        m_bits.realData = 0;
        m_real = nullptr;
        m_complex = reinterpret_cast<T2*> (data);
    }
}

//! Copy constructor
template<typename T1, typename T2>
hwTMatrix<T1, T2>::hwTMatrix(const hwTMatrix<T1, T2>& source)
    : m_refCount(1)
{
    m_nCols = source.m_nCols;
    m_nRows = source.m_nRows;
    m_bits.ownData = 1;
    m_real_memory = nullptr;
    m_complex_memory = nullptr;
    m_real = nullptr;
    m_complex = nullptr;

    hwMathStatus status = Copy(source);

    if (!status.IsOk())
    {
        MakeEmpty();
        return;
    }
}

//! Destructor
template<typename T1, typename T2>
hwTMatrix<T1, T2>::~hwTMatrix()
{
    if (m_refCount > 1)
        return;

    Deallocate();
}

//! Implement the = operator
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator=(const hwTMatrix<T1, T2>& rhs)
{
    if (this == &rhs)
        return *this;

    hwMathStatus status = Copy(rhs);

    if (!status.IsOk())
        MakeEmpty();

    return *this;
}

// ****************************************************
//             Data Type and Ownership
// ****************************************************

//! Determine if the matrix contains real or complex data
//! (i.e. no non-zero imaginary components)
template<typename T1, typename T2>
bool hwTMatrix<T1, T2>::IsRealData() const
{
    if (m_complex)
    {
        for (int i = 0; i < Size(); ++i)
        {
            if (m_complex[i].Imag() != 0.0)
                return false;
        }
    }

    return true;
}

//! Determine if the matrix is diagonal
template<typename T1, typename T2>
bool hwTMatrix<T1, T2>::IsDiag() const
{
    if (IsReal())
    {
        for (int jj = 0; jj < m_nCols; ++jj)
        {
            for (int ii = 0; ii < m_nCols; ++ii)
            {
                if (ii == jj)
                    continue;

                if ((*this)(ii, jj) != static_cast<T1> (0))
                    return false;
            }
        }
    }
    else
    {
        for (int jj = 0; jj < m_nCols; ++jj)
        {
            for (int ii = 0; ii < m_nCols; ++ii)
            {
                if (ii == jj)
                    continue;

                if (z(ii, jj) != static_cast<T1> (0))
                    return false;
            }
        }
    }

    return true;
}

// ****************************************************
//                Set Matrix Dimensions
// ****************************************************

//! Dimension matrix as a vector, specifying real or complex
//! Applies to empty or populated matrix, does not preserve existing data
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Dimension(int size, DataType dataType)
{
    if (size)
        return Dimension(size, 1, dataType);
    else
        return Dimension(0, 0, dataType);
}

//! Dimension matrix for double indexing, specifying real or complex
//! Applies to empty or populated matrix, does not preserve existing data
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Dimension(int m, int n, DataType dataType)
{
    if (m_real)
    {
        if (dataType == REAL)
        {
            if (m_nRows == m && m_nCols == n)
                return hwMathStatus();
            else
                MakeEmpty();
        }
        else
            MakeEmpty();
    }
    else if (m_complex)
    {
        if (dataType == COMPLEX)
        {
            if (m_nRows == m && m_nCols == n)
                return hwMathStatus();
            else
                MakeEmpty();
        }
        else
            MakeEmpty();
    }

    m_nRows = m;
    m_nCols = n;

    try
    {
        Allocate(dataType);
    }
    catch (std::bad_alloc&)
    {
        return hwMathStatus(HW_MATH_ERR_ALLOCFAILED, 0);
    }

    if (!m_complex)
        m_bits.realData = 1;
    else
        m_bits.realData = 0;

    return hwMathStatus();
}

//! Resize the matrix as a vector, leaving the data type unchanged
//! Applies to populated matrix, can only resize an empty matrix if it remains empty
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Resize(int size, bool initZero)
{
    if (size)
        return Resize(size, 1, initZero);
    else
        return Resize(0, 0);
}

//! Resize the matrix for double indexing, leaving the data type unchanged
//! Applies to populated matrix, can only resize an empty matrix if it remains empty
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Resize(int m, int n, bool initZero)
{
    if (m < 0)
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 1);

    if (n < 0)
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 2);

    int m_old = m_nRows;
    int n_old = m_nCols;

	if ((m_old == m) && (n_old == n))
		return hwMathStatus();

    if (m == 0 || n == 0)
    {
        MakeEmpty();
        m_nRows = m;
        m_nCols = n;
        return hwMathStatus();
    }

    m_nRows = m;
    m_nCols = n;

    // this section should be in SetCapacity(), but getting it there
    // requires rethinking the allocation/copy functionality.
    // simply change dimensions if capacity is sufficient, and
    // object is a vector or only higher dimensions are changing
    if (m <= m_capacity / n)
    {
        bool capacityIsOK;

        if (n_old == 1 && n == 1)
        {
            capacityIsOK = true;
        }
        else if (m_old == m)
        {
            capacityIsOK = true;
        }
        else
        {
            // move columns to new locations, starting with the last
            // column in memory. Do later if needed. Reallocate for now.
            capacityIsOK = false;
        }

        if (capacityIsOK)
        {
            if (initZero)
            {
                ZeroBlock(0, m_old-1, n_old, n-1);
                ZeroBlock(m_old, m-1, 0, n-1);
            }

            return hwMathStatus();
        }
    }

    if (m_real)
    {
        T1* m_real_old = m_real;
        char* m_real_memory_old = m_real_memory;
        m_real = nullptr;
        m_real_memory = nullptr;

        try
        {
            Allocate(REAL);
        }
        catch (std::bad_alloc&)
        {
            return hwMathStatus(HW_MATH_ERR_ALLOCFAILED, 0);
        }

        CopyBlock(m_real_old, m_old, n_old, 0, m_old-1, 0, n_old-1, 0, 0);

        if (initZero)
        {
            ZeroBlock(0, m_old-1, n_old, n-1);
            ZeroBlock(m_old, m-1, 0, n-1);
        }

        if (m_bits.ownData)
            FreeMemory(m_real_memory_old, m_real_old);
        else
            m_bits.ownData = 1;
    }
    else if (m_complex)
    {
        T2* m_complex_old = m_complex;
        char* m_complex_memory_old = m_complex_memory;
        m_complex = nullptr;
        m_complex_memory = nullptr;

        try
        {
            Allocate(COMPLEX);
        }
        catch (std::bad_alloc&)
        {
            return hwMathStatus(HW_MATH_ERR_ALLOCFAILED, 0);
        }

        CopyBlock(m_complex_old, m_old, n_old, 0, m_old-1, 0, n_old-1, 0, 0);

        if (initZero)
        {
            ZeroBlock(0, m_old-1, n_old, n-1);
            ZeroBlock(m_old, m-1, 0, n-1);
        }

        if (m_bits.ownData)
            FreeMemory(m_complex_memory_old, m_complex_old);
        else
            m_bits.ownData = 1;
    }
    else if (m != 0 && n != 0)    // cannot allocate an empty matrix with Resize
	    return hwMathStatus(HW_MATH_ERR_NOTALLOWED);

    return hwMathStatus();
}

//! Change the dimensions of a matrix while maintaining the same number of elements
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Reshape(int m, int n)
{
    if (m < 0 && m != -1)
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 1);

    if (n < 0 && n != -1)
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 2);

    int size = Size();

    if (m == -1)
    {
        if (n == -1)
            return hwMathStatus(HW_MATH_ERR_MATRIXRESHAPE1, 1, 2);

        m = size / n;
    }
    else if (n == -1)
    {
        n = size / m;
    }

    if (m * n != Size())
        return hwMathStatus(HW_MATH_ERR_MATRIXRESHAPE2, 1, 2);

    m_nRows = m;
    m_nCols = n;

    return hwMathStatus();
}

// ****************************************************
//                  Matrix Properties
// ****************************************************

//! Determine if the matrix is symmetric
template<typename T1, typename T2>
bool hwTMatrix<T1, T2>::IsSymmetric(T1 tol) const
{
    int i, j;
    int m = m_nRows;
    int n = m_nCols;

    if (m != n)
        return false;

    if (tol == (T1) 0)
    {
        if (m_real)
        {
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < i; ++j)
                {
                    if ((*this)(i, j) != (*this)(j, i))
                        return false;
                }
            }
        }

        if (m_complex)
        {
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < i; ++j)
                {
                    if (z(i, j) != z(j, i))
                        return false;
                }
            }
        }
    }
    else
    {
        T1 norm;
        hwMathStatus status = Norm(norm, "inf");

        if (m_real)
        {
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < i; ++j)
                {
                    if (fabs((*this)(i, j) - (*this)(j, i)) > norm * tol)
                        return false;
                }
            }
        }

        if (m_complex)
        {
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < i; ++j)
                {
                    if ((z(i, j) - z(j, i)).Mag() > norm * tol)
                        return false;
                }
            }
        }
    }

    return true;
}

//! Determine if the matrix is Hermitian
template<typename T1, typename T2>
bool hwTMatrix<T1, T2>::IsHermitian(T1 tol) const
{
    int i, j;
    int m = m_nRows;
    int n = m_nCols;

    if (m != n)
        return false;

    if (tol == (T1) 0)
    {
        if (m_real)
        {
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < i; ++j)
                {
                    if ((*this)(i, j) != (*this)(j, i))
                        return false;
                }
            }
        }

        if (m_complex)
        {
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < i; ++j)
                {
                    if (z(i, j) != z(j, i).Conjugate())
                        return false;
                }

                if (z(i, i).Imag() != 0.0)  // real diagonal
                    return false;
            }
        }
    }
    else
    {
        T1 norm;
        hwMathStatus status = Norm(norm, "inf");

        if (m_real)
        {
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < i; ++j)
                {
                    if (fabs((*this)(i, j) - (*this)(j, i)) > norm * tol)
                        return false;
                }
            }
        }

        if (m_complex)
        {
            for (i = 0; i < m; ++i)
            {
                for (j = 0; j < i; ++j)
                {
                    if ((z(i, j) - z(j, i).Conjugate()).Mag() > norm * tol)
                        return false;
                }
            }

            if (fabs(z(i, i).Imag()) > norm * tol)  // real diagonal
                return false;
        }
    }

    return true;
}

#include <GeneralFuncs.h>

//! Determine if the matrix contains non-finite elements
template<typename T1, typename T2>
bool hwTMatrix<T1, T2>::IsFinite() const
{
    int i = 0;
    int size = Size();

    if (IsReal())
    {
        for (i = 0; i < size; ++i) 
        {
            if (IsNaN_T(m_real[i]) || IsInf_T(m_real[i]) || IsNegInf_T(m_real[i]))
                return false;
        }
    }
    else
    {
        for (i = 0; i < size; ++i) 
        {
            T1 real = m_complex[i].Real();
            T1 imag = m_complex[i].Imag();

            if (IsNaN_T(real) || IsInf_T(real) || IsNegInf_T(real))
                return false;
            if (IsNaN_T(imag) || IsInf_T(imag) || IsNegInf_T(imag))
                return false;
        }
    }

    return true;
}

//! Determine if the matrix contains non-finite elements
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::IsFinite(hwTMatrix<bool>& R) const
{
    int i = 0;
    int size = Size();
    hwMathStatus status;

    status = R.Dimension(m_nRows, m_nCols, hwTMatrix<bool>::REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (IsReal())
    {
        for (i = 0; i < size; ++i) 
        {
            if (IsNaN_T(m_real[i]) || IsInf_T(m_real[i]) || IsNegInf_T(m_real[i])) 
                R(i) = false;
            else
                R(i) = true;
        }
    }
    else
    {
        for (i = 0; i < size; ++i) 
        {
            T1 real = m_complex[i].Real();
            T1 imag = m_complex[i].Imag();

            if (IsNaN_T(real) || IsInf_T(real) || IsNegInf_T(real) || IsNaN_T(imag) || IsInf_T(imag) || IsNegInf_T(imag))
			{
                R(i) = false;
			}
			else
			{
				R(i) = true;
			}
        }
    }

    return status;
}

//! Determine if the matrix contains non-finite elements
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::IsFinite(hwTMatrix<double>& R) const
{
    int i = 0;
    int size = Size();
    hwMathStatus status;

    status = R.Dimension(m_nRows, m_nCols, hwTMatrix<double>::REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (IsReal())
    {
        for (i = 0; i < size; ++i) 
        {
            if (IsNaN_T(m_real[i]) || IsInf_T(m_real[i]) || IsNegInf_T(m_real[i])) 
                R(i) = 0.0;
            else
                R(i) = 1.0;
        }
    }
    else
    {
        for (i = 0; i < size; ++i) 
        {
            T1 real = m_complex[i].Real();
            T1 imag = m_complex[i].Imag();

            if (IsNaN_T(real) || IsInf_T(real) || IsNegInf_T(real) || IsNaN_T(imag) || IsInf_T(imag) || IsNegInf_T(imag))
			{
                R(i) = 0.0;
			}
			else
			{
				R(i) = 1.0;
			}
        }
    }

    return status;
}

// ****************************************************
//         Access Functions for Real Elements
// ****************************************************

//! Return a reference to the real matrix element at the specified single index
template<typename T1, typename T2>
T1& hwTMatrix<T1, T2>::operator()(int index, IndexDir dir)
{
    // m = m_nRows
    // n = m_nCols
    // index = col*m + row
    if (dir == BY_COL)
    {
        // row = index % m;
        // col = index / m;
        return m_real[index];
    }
    else    // dir == BY_ROW
    {
        // row = index / n;
        // col = index % n;
        return m_real[(index%m_nCols)*m_nRows + (index/m_nCols)];
    }
}

//! Return a const reference to the real data element at the specified single index
template<typename T1, typename T2>
const T1& hwTMatrix<T1, T2>::operator()(int index, IndexDir dir) const
{
    if (dir == BY_COL)
        return m_real[index];
    else    // dir == BY_ROW
        return m_real[(index%m_nCols)*m_nRows + (index/m_nCols)];
}

//! Return a reference to the real matrix element at the specified index pair
template<typename T1, typename T2>
T1& hwTMatrix<T1, T2>::operator()(int i, int j)
{
    // m = m_nRows
    // n = m_nCols
    // index = j*m + i
    return m_real[j*m_nRows + i];
}

//! Return a reference to the real data element at the specified index pair
template<typename T1, typename T2>
const T1& hwTMatrix<T1, T2>::operator()(int i, int j) const
{
    return m_real[j*m_nRows + i];
}

//! Set the element at the specified single index to a real value
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::SetElement(int index, T1 real, IndexDir dir)
{
    if (m_real)
        (*this)(index, dir) = real;
    else if (m_complex)
        z(index, dir) = real;
}

//! Set the element at the specified index pair to a real value
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::SetElement(int i, int j, T1 real)
{
    if (m_real)
        (*this)(i, j) = real;
    else if (m_complex)
        z(i, j) = real;
}

//! Set every matrix element to a specified real value
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::SetElements(T1 value)
{
    int size = Size();

    if (m_real)
    {
        for (int i = 0; i < size; ++i)
            m_real[i] = value;
    }
    else if (m_complex)
    {
        for (int i = 0; i < size; ++i)
            m_complex[i] = value;
    }
}

// ****************************************************
//        Access Functions for Complex Elements
// ****************************************************

//! Return a reference to the complex matrix element at the specified single index
template<typename T1, typename T2>
T2& hwTMatrix<T1, T2>::z(int index, IndexDir dir)
{
    if (dir == BY_COL)
        return m_complex[index];
    else    // dir == BY_ROW
        return m_complex[(index%m_nCols)*m_nRows + (index/m_nCols)];
}

//! Return a reference to the complex matrix element at the specified single index
template<typename T1, typename T2>
const T2& hwTMatrix<T1, T2>::z(int index, IndexDir dir) const
{
    if (dir == BY_COL)
        return m_complex[index];
    else    // dir == BY_ROW
        return m_complex[(index%m_nCols)*m_nRows + (index/m_nCols)];
}

//! Return a reference to the complex matrix element at the specified index pair
template<typename T1, typename T2>
T2& hwTMatrix<T1, T2>::z(int i, int j)
{
    return m_complex[j*m_nRows + i];
}

//! Return a reference to the complex matrix element at the specified index pair
template<typename T1, typename T2>
const T2& hwTMatrix<T1, T2>::z(int i, int j) const
{
    return m_complex[j*m_nRows + i];
}

//! Set the element at the specified single index to a complex value
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::SetElement(int index, T2 cplx, IndexDir dir)
{
    if (m_real)
        MakeComplex();

    z(index, dir) = cplx;
}

//! Set the element at the specified index pair to a complex value
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::SetElement(int i, int j, T2 cplx)
{
    if (m_real)
        MakeComplex();

    z(i, j) = cplx;
}

//! Set every matrix element to a specified complex value
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::SetElements(T2 value)
{
    if (IsReal())
        return hwMathStatus(HW_MATH_ERR_ARRAYTYPE);

    int size = Size();

    for (int i = 0; i < size; ++i)
        m_complex[i] = value;

    return hwMathStatus();
}

// ****************************************************
//             Real / Complex Conversions
// ****************************************************

//! Convert the matrix to complex
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::MakeComplex()
{
    if (m_complex || !m_real)      // already complex, or empty
        return hwMathStatus();

    // Note: an empty matrix is real by convention, but there is
    // no reason to trigger an error. The client will still have
    // to dimension the matrix somewhere.

    try
    {
        Allocate(COMPLEX);
    }
    catch (std::bad_alloc&)
    {
        return hwMathStatus(HW_MATH_ERR_ALLOCFAILED, 0);
    }

    int size = Size();

    for (int i = 0; i < size; ++i)
    {
        m_complex[i].Real() = m_real[i];
        m_complex[i].Imag() = (T1) 0;
    }

    if (m_bits.ownData)
        FreeMemory(m_real_memory, m_real);
    else
        m_bits.ownData = 1;

    m_real = nullptr;
    m_bits.realData = 0;

    return hwMathStatus();
}

//! Pack real and imaginary components into a complex matrix
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::PackComplex(const hwTMatrix<T1, T2>& real, const hwTMatrix<T1, T2>* imag)
{
    hwMathStatus status;

    if (!real.IsReal())
        return status(HW_MATH_ERR_ARRAYTYPE, 1);

    if (imag)
    {
        if (real.m_nRows != imag->m_nRows)
            return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        if (real.m_nCols != imag->m_nCols)
            return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        if (!imag->IsReal())
            return status(HW_MATH_ERR_ARRAYTYPE, 2);
    }

    status = Dimension(real.m_nRows, real.m_nCols, COMPLEX);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    int size = Size();
    T1* pReal = real.m_real;
    T1* pComplex = reinterpret_cast<T1*> (m_complex);

    if (imag)
    {
        T1* pImag = imag->m_real;

        while (size--)
        {
            *pComplex++ = *pReal++;
            *pComplex++ = *pImag++;
        }
    }
    else
    {
        while (size--)
        {
            *pComplex++ = *pReal++;
            *pComplex++ = (T1) 0;
        }
    }

    return status;
}

//! Unpack real and imaginary components from a complex matrix
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::UnpackComplex(hwTMatrix<T1, T2>* real, hwTMatrix<T1, T2>* imag) const
{
    hwMathStatus status;

    if (IsReal())
        return status(HW_MATH_ERR_ARRAYTYPE, 0);

    if (real)
    {
        status = real->Dimension(m_nRows, m_nCols, REAL);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }
    }

    if (imag)
    {
        status = imag->Dimension(m_nRows, m_nCols, REAL);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }
    }

    int size = Size();
    T1* pComplex = reinterpret_cast<T1*> (m_complex);

    if (real && imag)
    {
        T1* pReal = real->m_real;
        T1* pImag = imag->m_real;

        while (size--)
        {
            *pReal++ = *pComplex++;
            *pImag++ = *pComplex++;
        }
    }
    else if (real)
    {
        T1* pReal = real->m_real;

        while (size--)
        {
            *pReal++ = *pComplex++;
            pComplex++;             // skip imag part
        }
    }
    else if (imag)
    {
        T1* pImag = imag->m_real;

        while (size--)
        {
            pComplex++;             // skip real part
            *pImag++ = *pComplex++;
        }
    }

    return status;
}

// ****************************************************
//                  Vector Operations
// ****************************************************

//! Insert a vector at the specified location of a source vector
//! and store the result in the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::InsertElements(const hwTMatrix<T1, T2>& source, int startElem, const hwTMatrix<T1, T2>& elems)
{
    // client must make sure source and elems have same data type
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    if (!source.IsVector())
        return status(HW_MATH_ERR_NOVECTOR, 1);

    if (!elems.IsVector())
        return status(HW_MATH_ERR_NOVECTOR, 3);

    int sourceSize = source.Size();
    int elemsSize = elems.Size();

    if (startElem < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 2);

    if (startElem > sourceSize)
        return status(HW_MATH_ERR_INVALIDINDEX, 2);

    if (source.IsReal())
    {
        if (!elems.IsReal())
            return status(HW_MATH_ERR_ARRAYTYPE, 1, 3);
    }
    else if (elems.IsReal())
    {
        return status(HW_MATH_ERR_ARRAYTYPE, 1, 3);
    }

    status = Dimension(sourceSize+elemsSize, 1, elems.Type());

    if (!status.IsOk())
        return status;

    if (m_real)
    {
        CopyBlock(source.m_real, sourceSize, 1, 0, startElem-1, 0, 0, 0, 0);
        CopyBlock(source.m_real, sourceSize, 1, startElem, sourceSize-1, 0, 0, startElem+elemsSize, 0);
        CopyBlock(elems.m_real, elemsSize, 1, 0, elemsSize-1, 0, 0, startElem, 0);
    }
    else if (m_complex)
    {
        CopyBlock(source.m_complex, sourceSize, 1, 0, startElem-1, 0, 0, 0, 0);
        CopyBlock(source.m_complex, sourceSize, 1, startElem, sourceSize-1, 0, 0, startElem+elemsSize, 0);
        CopyBlock(elems.m_complex, elemsSize, 1, 0, elemsSize-1, 0, 0, startElem, 0);
    }

    if (source.m_nCols != 1)
        Transpose();

    return status;
}

//! Insert a set zero elements at the specified location of a source vector
//! and store the result in the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::InsertElements(const hwTMatrix<T1, T2>& source, int startElem, int numElems)
{
    // client must make sure source and elems have same data type
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    if (!source.IsVector())
        return status(HW_MATH_ERR_NOVECTOR, 1);

    int sourceSize = source.Size();

    if (startElem < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 2);

    if (startElem > sourceSize)
        return status(HW_MATH_ERR_INVALIDINDEX, 2);

    if (numElems < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 3);

    status = Dimension(sourceSize+numElems, 1, source.Type());

    if (!status.IsOk())
        return status;

    if (m_real)
    {
        CopyBlock(source.m_real, sourceSize, 1, 0, startElem-1, 0, 0, 0, 0);
        CopyBlock(source.m_real, sourceSize, 1, startElem, sourceSize-1, 0, 0, startElem+numElems, 0);
    }
    else if (m_complex)
    {
        CopyBlock(source.m_complex, sourceSize, 1, 0, startElem-1, 0, 0, 0, 0);
        CopyBlock(source.m_complex, sourceSize, 1, startElem, sourceSize-1, 0, 0, startElem+numElems, 0);
    }

    ZeroBlock(startElem, startElem+numElems-1, 0, 0);

    if (source.m_nCols != 1)
        Transpose();

    return status;
}

//! Delete a number of elements from a source vector at the specified location
//! and store the result in the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::DeleteElements(const hwTMatrix<T1, T2>& source, int startElem, int numElems)
{
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    if (!source.IsVector())
        return status(HW_MATH_ERR_NOVECTOR, 1);

    int sourceSize = source.Size();

    if (startElem < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 2);

    if (numElems <= 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 3);

    if (startElem >= sourceSize + 1 - numElems)
        return status(HW_MATH_ERR_INVALIDINDEX, 2, 3);

    status = Dimension(sourceSize-numElems, 1, source.Type());

    if (!status.IsOk())
        return status;

    if (source.m_real)
    {
        CopyBlock(source.m_real, sourceSize, 1, 0, startElem-1, 0, 0, 0, 0);
        CopyBlock(source.m_real, sourceSize, 1, startElem+numElems, sourceSize-1, 0, 0, startElem, 0);
    }
    else if (source.m_complex)
    {
        CopyBlock(source.m_complex, sourceSize, 1, 0, startElem-1, 0, 0, 0, 0);
        CopyBlock(source.m_complex, sourceSize, 1, startElem+numElems, sourceSize-1, 0, 0, startElem, 0);
    }

    if (source.m_nRows == 1)
        Transpose();

    return status;
}

//! Delete a number of elements from the calling object at the specified location
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::DeleteElements(int startElem, int numElems)
{
    if (m_real)
    {
        hwTMatrix<T1, T2> temp(m_nRows, m_nCols, (void*) m_real, REAL);
        temp.OwnData(true);
        temp.m_real_memory = m_real_memory;
        m_real = nullptr;
        m_real_memory = nullptr;
        MakeEmpty();
        return DeleteElements(temp, startElem, numElems);
    }
    else if (m_complex)
    {
        hwTMatrix<T1, T2> temp(m_nRows, m_nCols, (void*) m_complex, COMPLEX);
        temp.OwnData(true);
        temp.m_complex_memory = m_complex_memory;
        m_complex = nullptr;
        m_complex_memory = nullptr;
        MakeEmpty();
        return DeleteElements(temp, startElem, numElems);
    }

    return hwMathStatus();
}

// ****************************************************
//                    Row Operations
// ****************************************************

//! Create a single vector from the rows of a matrix and store
//! the result in the object on which the function is called
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::ConcatRows(const hwTMatrix<T1, T2>& source)
{
    int size = source.Size();
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    status = Dimension(size, source.Type());

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (m_real)
    {
        for (int i = 0; i < size; ++i)
            (*this)(i) = source(i, BY_ROW);
    }

    if (m_complex)
    {
        for (int i = 0; i < size; ++i)
            z(i) = source.z(i, BY_ROW);
    }

    return status;
}

//! Insert a set of row vectors at the specified location of a source matrix
//! and store the result in the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::InsertRows(const hwTMatrix<T1, T2>& source, int startRow, const hwTMatrix<T1, T2>& rows)
{
    // client must make sure source and rows have same data type
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    if (startRow < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 2);

    if (startRow > source.m_nRows)
        return status(HW_MATH_ERR_INVALIDINDEX, 2);

    if (source.m_nCols != rows.m_nCols)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 3);

    if (source.Type() != rows.Type())
        return status(HW_MATH_ERR_ARRAYTYPE, 1, 3);

    status = Dimension(source.m_nRows+rows.m_nRows, source.m_nCols, source.Type());

    if (!status.IsOk())
        return status;

    if (m_real)
    {
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, 0, startRow-1, 0, source.m_nCols-1, 0, 0);
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, startRow, source.m_nRows-1, 0, source.m_nCols-1, startRow+rows.m_nRows, 0);
        CopyBlock(rows.m_real, rows.m_nRows, rows.m_nCols, 0, rows.m_nRows-1, 0, rows.m_nCols-1, startRow, 0);
    }
    else if (m_complex)
    {
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, 0, startRow-1, 0, source.m_nCols-1, 0, 0);
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, startRow, source.m_nRows-1, 0, source.m_nCols-1, startRow+rows.m_nRows, 0);
        CopyBlock(rows.m_complex, rows.m_nRows, rows.m_nCols, 0, rows.m_nRows-1, 0, rows.m_nCols-1, startRow, 0);
    }

    return status;
}

//! Insert a set of zeroed rows at the specified location of a source matrix
//! and store the result in the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::InsertRows(const hwTMatrix<T1, T2>& source, int startRow, int numRows)
{
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    if (startRow < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 2);

    if (startRow > source.m_nRows)
        return status(HW_MATH_ERR_INVALIDINDEX, 2);

    if (numRows < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 3);

    status = Dimension(source.m_nRows+numRows, source.m_nCols, source.Type());

    if (!status.IsOk())
        return status;

    if (m_real)
    {
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, 0, startRow-1, 0, source.m_nCols-1, 0, 0);
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, startRow, source.m_nRows-1, 0, source.m_nCols-1, startRow+numRows, 0);
    }
    else if (m_complex)
    {
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, 0, startRow-1, 0, source.m_nCols-1, 0, 0);
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, startRow, source.m_nRows-1, 0, source.m_nCols-1, startRow+numRows, 0);
    }

    ZeroBlock(startRow, startRow+numRows-1, 0, source.m_nCols-1);

    return status;
}

//! Delete a number of rows from a source matrix at the specified location
//! and store the result in the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::DeleteRows(const hwTMatrix<T1, T2>& source,
                                           int startRow, int numRows)
{
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    if (startRow < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 2);

    if (numRows <= 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 3);

    if (startRow >= source.m_nRows + 1 - numRows)
        return status(HW_MATH_ERR_INVALIDINDEX, 2, 3);

    status = Dimension(source.m_nRows-numRows, source.m_nCols, source.Type());

    if (!status.IsOk())
        return status;

    if (source.m_real)
    {
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, 0, startRow-1, 0, source.m_nCols-1, 0, 0);
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, startRow+numRows, source.m_nRows-1, 0, source.m_nCols-1, startRow, 0);
    }
    else if (source.m_complex)
    {
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, 0, startRow-1, 0, source.m_nCols-1, 0, 0);
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, startRow+numRows, source.m_nRows-1, 0, source.m_nCols-1, startRow, 0);
    }

    return status;
}

//! Delete a number of rows from the calling object at the specified location
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::DeleteRows(int startRow, int numRows)
{
    if (m_real)
    {
        hwTMatrix<T1, T2> temp(m_nRows, m_nCols, (void*) m_real, REAL);
        temp.OwnData(true);
        temp.m_real_memory = m_real_memory;
        m_real = nullptr;
        m_real_memory = nullptr;
        MakeEmpty();
        return DeleteRows(temp, startRow, numRows);
    }
    else if (m_complex)
    {
        hwTMatrix<T1, T2> temp(m_nRows, m_nCols, (void*) m_complex, COMPLEX);
        temp.OwnData(true);
        temp.m_complex_memory = m_complex_memory;
        m_complex = nullptr;
        m_complex_memory = nullptr;        MakeEmpty();
        return DeleteRows(temp, startRow, numRows);
    }

    return hwMathStatus();
}

//! Read the row at the specified location
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::ReadRow(int rowNum, hwTMatrix<T1, T2>& row) const
{
    hwMathStatus status;

    if (rowNum < 0 || rowNum > m_nRows-1)
        return status(HW_MATH_ERR_INVALIDINDEX, 1);

    status = row.Dimension(1, m_nCols, Type());

    if (!status.IsOk())
        return status;

    if (IsReal())
    {
        row.CopyBlock(m_real, m_nRows, m_nCols, rowNum, rowNum, 0, m_nCols-1, 0, 0);
    }
    else
    {
        row.CopyBlock(m_complex, m_nRows, m_nCols, rowNum, rowNum, 0, m_nCols-1, 0, 0);
    }

    return status;
}

//! Write the row at the specified location
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::WriteRow(int rowNum, const hwTMatrix<T1, T2>& row)
{
    hwMathStatus status;

    if (rowNum < 0 || rowNum > m_nRows-1)
        return status(HW_MATH_ERR_INVALIDINDEX, 1);

    if (row.m_nRows != 1)
        return status(HW_MATH_ERR_ARRAYSIZE, 2);

    if (row.m_nCols != m_nCols)
        return status(HW_MATH_ERR_ARRAYSIZE, 2);

    if (IsReal())
    {
        if (row.IsReal())
            CopyBlock(row.m_real, 1, m_nCols, 0, 0, 0, m_nCols-1, rowNum, 0);
        else
        {
            MakeComplex();
            CopyBlock(row.m_complex, 1, m_nCols, 0, 0, 0, m_nCols-1, rowNum, 0);
        }
    }
    else
    {
        if (row.IsReal())
        {
            for (int j = 0; j < m_nCols; ++j)
                z(rowNum, j) = row(j);
        }
        else
            CopyBlock(row.m_complex, 1, m_nCols, 0, 0, 0, m_nCols-1, rowNum, 0);
    }

    return status;
}

// ****************************************************
//                   Column Operations
// ****************************************************

//! Create a single vector from the columns of a matrix and store
//! the result in the object on which the function is called
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::ConcatColumns(const hwTMatrix<T1, T2>& source)
{
    int size = source.Size();

    if (this == &source)
        return hwMathStatus(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    Copy(source);
    m_nRows = size;
    m_nCols = 1;

    return hwMathStatus();
}

//! Insert a set of column vectors at the specified location of a source matrix
//! and store the result in the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::InsertColumns(const hwTMatrix<T1, T2>& source, int startCol, const hwTMatrix<T1, T2>& cols)
{
    // client must make sure source and cols have same data type
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    if (startCol < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 2);

    if (startCol > source.m_nCols)
        return status(HW_MATH_ERR_INVALIDINDEX, 2);

    if (source.m_nRows != cols.m_nRows)
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 3);

    if (source.Type() != cols.Type())
        return status(HW_MATH_ERR_ARRAYTYPE, 1, 3);

    status = Dimension(source.m_nRows, source.m_nCols+cols.m_nCols, source.Type());

    if (!status.IsOk())
        return status;

    if (m_real)
    {
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, 0, startCol-1, 0, 0);
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, startCol, source.m_nCols-1, 0, startCol+cols.m_nCols);
        CopyBlock(cols.m_real, cols.m_nRows, cols.m_nCols, 0, cols.m_nRows-1, 0, cols.m_nCols-1, 0, startCol);
    }

    if (m_complex)
    {
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, 0, startCol-1, 0, 0);
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, startCol, source.m_nCols-1, 0, startCol+cols.m_nCols);
        CopyBlock(cols.m_complex, cols.m_nRows, cols.m_nCols, 0, cols.m_nRows-1, 0, cols.m_nCols-1, 0, startCol);
    }

    return status;
}

//! Insert a set of zeroed columns at the specified location of a source matrix
//! and store the result in the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::InsertColumns(const hwTMatrix<T1, T2>& source, int startCol, int numCols)
{
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    if (startCol < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 2);

    if (startCol > source.m_nCols)
        return status(HW_MATH_ERR_INVALIDINDEX, 2);

    if (numCols < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 3);

    status = Dimension(source.m_nRows, source.m_nCols+numCols, source.Type());

    if (!status.IsOk())
        return status;

    if (m_real)
    {
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, 0, startCol-1, 0, 0);
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, startCol, source.m_nCols-1, 0, startCol+numCols);
    }

    if (m_complex)
    {
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, 0, startCol-1, 0, 0);
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, startCol, source.m_nCols-1, 0, startCol+numCols);
    }

    ZeroBlock(0, source.m_nRows-1, startCol, startCol+numCols-1);

    return status;
}

//! Delete a number of columns from a source matrix at the specified location
//! and store the result in the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::DeleteColumns(const hwTMatrix<T1, T2>& source,
                                              int startCol, int numCols)
{
    hwMathStatus status;

    if (this == &source)
        return status(HW_MATH_ERR_NOTALLOWED);      // disallow in-place operation

    if (startCol < 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 2);

    if (numCols <= 0)
        return status(HW_MATH_ERR_NONNONNEGINT, 3);

    if (startCol >= source.m_nCols + 1 - numCols)
        return status(HW_MATH_ERR_INVALIDINDEX, 2, 3);

    status = Dimension(source.m_nRows, source.m_nCols-numCols, source.Type());

    if (!status.IsOk())
        return status;

    if (source.m_real)
    {
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, 0, startCol-1, 0, 0);
        CopyBlock(source.m_real, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, startCol+numCols, source.m_nCols-1, 0, startCol);
    }

    if (source.m_complex)
    {
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, 0, startCol-1, 0, 0);
        CopyBlock(source.m_complex, source.m_nRows, source.m_nCols, 0, source.m_nRows-1, startCol+numCols, source.m_nCols-1, 0, startCol);
    }

    return status;
}

//! Delete a number of columns from the calling object at the specified location
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::DeleteColumns(int startCol, int numCols)
{
    if (m_real)
    {
        hwTMatrix<T1, T2> temp(m_nRows, m_nCols, (void*) m_real, REAL);
        temp.OwnData(true);
        temp.m_real_memory = m_real_memory;
        m_real = nullptr;
        m_real_memory = nullptr;
        MakeEmpty();
        return DeleteColumns(temp, startCol, numCols);
    }
    else if (m_complex)
    {
        hwTMatrix<T1, T2> temp(m_nRows, m_nCols, (void*) m_complex, COMPLEX);
        temp.OwnData(true);
        temp.m_complex_memory = m_complex_memory;
        m_complex = nullptr;
        m_complex_memory = nullptr;
        MakeEmpty();
        return DeleteColumns(temp, startCol, numCols);
    }

    return hwMathStatus();
}

//! Read the column at the specified location
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::ReadColumn(int colNum, hwTMatrix<T1, T2>& col) const
{
    hwMathStatus status;

    if (colNum < 0 || colNum > m_nCols-1)
        return status(HW_MATH_ERR_INVALIDINDEX, 1);

    status = col.Dimension(m_nRows, 1, Type());

    if (!status.IsOk())
        return status;

    if (IsReal())
    {
        col.CopyBlock(m_real, m_nRows, m_nCols, 0, m_nRows-1, colNum, colNum, 0, 0);
    }
    else
    {
        col.CopyBlock(m_complex, m_nRows, m_nCols, 0, m_nRows-1, colNum, colNum, 0, 0);
    }

    return status;
}

//! Write the column at the specified location
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::WriteColumn(int colNum, const hwTMatrix<T1, T2>& col)
{
    if (colNum < 0 || colNum > m_nCols-1)
        return hwMathStatus(HW_MATH_ERR_INVALIDINDEX, 1);

    if (col.m_nCols != 1)
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 2);

    if (col.m_nRows != m_nRows)
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 2);

    if (IsReal())
    {
        if (col.IsReal())
            CopyBlock(col.m_real, m_nRows, 1, 0, m_nRows-1, 0, 0, 0, colNum);
        else
        {
            MakeComplex();
            CopyBlock(col.m_complex, m_nRows, 1, 0, m_nRows-1, 0, 0, 0, colNum);
        }
    }
    else
    {
        if (col.IsReal())
        {
            for (int i = 0; i < m_nRows; ++i)
                z(i, colNum) = col(i);
        }
        else
            CopyBlock(col.m_complex, m_nRows, 1, 0, m_nRows-1, 0, 0, 0, colNum);
    }

    return hwMathStatus();
}

// ****************************************************
//                   Submatrix Operations
// ****************************************************

//! Read a submatrix of a source, starting at the specified location, writing to the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::ReadSubmatrix(int startRow, int startCol, int numRows, int numCols, const hwTMatrix<T1, T2>& source)
{
    hwMathStatus status;

    int sourceM = source.m_nRows;
    int sourceN = source.m_nCols;

    if (startRow < 0 || startRow > sourceM - numRows)
        return status(HW_MATH_ERR_INVALIDINDEX, 1);

    if (startCol < 0 || startCol > sourceN - numCols)
        return status(HW_MATH_ERR_INVALIDINDEX, 2);

    status = Dimension(numRows, numCols, source.Type());

    if (!status.IsOk())
        return status;

    if (IsReal())
        CopyBlock(source.m_real, sourceM, sourceN,
            startRow, startRow + numRows, startCol, startCol + numCols, 0, 0);
    else
        CopyBlock(source.m_complex, sourceM, sourceN,
            startRow, startRow + numRows, startCol, startCol + numCols, 0, 0);

    return status;
}

//! Write a submatrix source to the calling object, starting at the specified location
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::WriteSubmatrix(int startRow, int startCol, const hwTMatrix<T1, T2>& source)
{
    hwMathStatus status;

    if (startRow < 0)
        return status(HW_MATH_ERR_INVALIDINDEX, 1);

    if (startCol < 0)
        return status(HW_MATH_ERR_INVALIDINDEX, 2);

    if (source.IsReal())
    {
        if (!IsReal())
        {
            hwTMatrix<T1, T2> temp;
            temp.PackComplex(source);

            return WriteSubmatrix(startRow, startCol, temp);
        }
    }
    else    // source is complex
    {
        if (IsReal())
            MakeComplex();
    }

    int sourceM = source.m_nRows;
    int sourceN = source.m_nCols;
    int targetM = m_nRows;
    int targetN = m_nCols;

    if (startRow + sourceM > targetM)
        targetM = startRow + sourceM;

    if (startCol + sourceN > targetN)
        targetN = startCol + sourceN;

    if (IsEmpty())
    {
        status = Dimension(targetM, targetN, source.Type());
        SetElements((T1) 0);
    }
    else
    {
        status = Resize(targetM, targetN, true);
    }

    if (!status.IsOk())
        return status;

    if (IsReal())
        CopyBlock(source.m_real, sourceM, sourceN,
                  0, sourceM-1, 0, sourceN-1, startRow, startCol);
    else
        CopyBlock(source.m_complex, sourceM, sourceN,
                  0, sourceM-1, 0, sourceN-1, startRow, startCol);

    return status;
}

// ****************************************************
//           Misc Linear Algebra Operations
// ****************************************************

//! Set the diagonal to ones, with zeros elsewhere
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::Identity()
{
    ZeroBlock(0, m_nRows-1, 0, m_nCols-1);

    int max = _min(m_nRows, m_nCols);

    if (m_real)
    {
        for (int i = 0; i < max; ++i)
            (*this)(i, i) = (T1) 1;
    }

    if (m_complex)
    {
        for (int i = 0; i < max; ++i)
            z(i, i).Set((T1) 1, (T1) 0);
    }
}

//! Generate a diagonal matrix or extract a diagonal from a matrix
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Diag(const hwTMatrix<T1, T2>& source, int k)
{
    hwMathStatus status;

    if (source.IsVector())
    {
        int size = source.Size();
        int n = size + abs(k);

        status = Dimension(n, n, source.Type());

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        SetElements(0.0);

        if (source.IsReal())   // real
        {

            if (k >= 0)
            {
                for (int i = 0; i < size; i++)
                    (*this)(i, i+k) = source(i);
            }
            else
            {
                for (int i = 0; i < size; i++)
                    (*this)(i-k, i) = source(i);
            }
        }
        else                // complex
        {
            if (k >= 0)
            {
                for (int i = 0; i < size; i++)
                    z(i, i+k) = source.z(i);
            }
            else
            {
                for (int i = 0; i < size; i++)
                    z(i-k, i) = source.z(i);
            }
        }
    }
    else    // matrix
    {
        int m = source.m_nRows;
        int n = source.m_nCols;
        int length;

        if (m >= n)
        {
            if (k >= 0)
                length = n - k;
            else if (k < n - m)
                length = m + k;
            else
                length = n;
        }
        else
        {
            if (k <= 0)
                length = m + k;
            else if (k > n - m)
                length = n - k;
            else
                length = m;
        }

        if (length > 0)
        {
            status = Dimension(length, 1, source.Type());

            if (!status.IsOk())
            {
                status.ResetArgs();
                return status;
            }

            if (source.IsReal())   // real
            {
                if (k >= 0)
                {
                    for (int i = 0; i < length; i++)
                        (*this)(i) = source(i, i+k);
                }
                else
                {
                    for (int i = 0; i < length; i++)
                         (*this)(i) = source(i-k, i);
                }
            }
            else                // complex
            {
                if (k >= 0)
                {
                    for (int i = 0; i < length; i++)
                        z(i) = source.z(i, i+k);
                }
                else
                {
                    for (int i = 0; i < length; i++)
                        z(i) = source.z(i-k, i);
                }
            }
        }
        else if (m == 0 || n == 0)
        {
            status = Dimension(0, 0, REAL);
        }
        else
        {
            status = Dimension(0, 1, REAL);
        }
    }

    return status;
}

//! Transpose the matrix in place
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Transpose()
{
    if (!m_bits.ownData)
    {
        // move data to another object and transpose out of place
        void* data;

        if (IsReal())
            data = (void*) m_real;
        else
            data = (void*) m_complex;

        hwTMatrix<T1, T2> temp(m_nRows, m_nCols, data, Type());
        return Transpose(temp);
    }

    if (IsVector())
    {
        int temp = m_nCols;
        m_nCols = m_nRows;
        m_nRows = temp;
    }
    else if (IsSquare())
    {
        if (m_real)
        {
            T1 tempT1;

            for (int i = 0; i < m_nCols; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    tempT1 = (*this)(i, j);
                    (*this)(i, j) = (*this)(j, i);
                    (*this)(j, i) = tempT1;
                }
            }
        }

        if (m_complex)
        {
            T2 tempT2;

            for (int i = 0; i < m_nCols; ++i)
            {
                for (int j = 0; j < i; ++j)
                {
                    tempT2 = z(i, j);
                    z(i, j) = z(j, i);
                    z(j, i) = tempT2;
                }
            }
        }
    }
    else
    {
        // eventually replace this with an in place transpose
        int m = m_nRows;
        int n = m_nCols;

        int temp = m_nCols;
        m_nCols = m_nRows;
        m_nRows = temp;

        if (m_real)
        {
            T1* temp_real = m_real;
            char* temp_real_memory = m_real_memory;
            hwTMatrix<T1, T2> tempM(m, n, temp_real, REAL);

            m_real = nullptr;

            try
            {
                Allocate(REAL);
            }
            catch (std::bad_alloc&)
            {
                return hwMathStatus(HW_MATH_ERR_ALLOCFAILED, 0);
            }

            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                    (*this)(j, i) = tempM(i, j);
            }

            FreeMemory(temp_real_memory, temp_real);
        }

        if (m_complex)
        {
            T2* temp_complex = m_complex;
            char* temp_complex_memory = m_complex_memory;
            hwTMatrix<T1, T2> tempM(m, n, temp_complex, COMPLEX);

            m_complex = nullptr;

            try
            {
                Allocate(COMPLEX);
            }
            catch (std::bad_alloc&)
            {
                return hwMathStatus(HW_MATH_ERR_ALLOCFAILED, 0);
            }

            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                    z(j, i) = tempM.z(i, j);
            }

            FreeMemory(temp_complex_memory, temp_complex);
        }
    }

    return hwMathStatus();
}

//! Transpose the matrix
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Transpose(const hwTMatrix<T1, T2>& source)
{
    hwMathStatus status = Dimension(source.m_nCols, source.m_nRows, source.Type());

    if (source.IsVector())
    {
        int size = Size();

        if (m_real)
        {
            for (int i = 0; i < size; ++i)
                m_real[i] = source.m_real[i];
        }

        if (m_complex)
        {
            for (int i = 0; i < size; ++i)
                m_complex[i] = source.m_complex[i];
        }
    }
    else
    {
        int m = m_nRows;
        int n = m_nCols;

        if (m_real)
        {
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                    (*this)(i, j) = source(j, i);
            }
        }

        if (m_complex)
        {
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                    z(i, j) = source.z(j, i);
            }
        }
    }

    return status;
}

//! Conjugate the matrix in place
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::Conjugate()
{
    if (IsReal())
        return;

    int size = Size();

    for (int i = 0; i < size; ++i)
        m_complex[i].Imag() = -m_complex[i].Imag();
}

//! Conjugate the matrix
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::Conjugate(const hwTMatrix<T1, T2>& source)
{
    if (source.IsReal())
    {
        *this = source;
    }
    else
    {
        hwMathStatus status = Dimension(source.m_nRows, source.m_nCols, COMPLEX);

        int size = Size();

        for (int i = 0; i < size; ++i)
        {
            m_complex[i].Real() =  source.m_complex[i].Real();
            m_complex[i].Imag() = -source.m_complex[i].Imag();
        }
    }
}

//! Transpose and conjugate the matrix in place
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::Hermitian()
{
    Transpose();
    Conjugate();
}

//! Transpose and conjugate the matrix
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::Hermitian(const hwTMatrix<T1, T2>& source)
{
    Transpose(source);
    Conjugate();
}

// ****************************************************
//               Arithmetic Operations
// ****************************************************

//! Add two matrices so that (*this) = A + B
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Add(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    // check dimensions
    hwMathStatus status;
    int m = A.m_nRows;
    int n = A.m_nCols;
    int size = A.Size();

    if (A.m_nRows != B.m_nRows || A.m_nCols != B.m_nCols)
    {
        if (size == 0 && B.Size() == 1)
            status = Dimension(0, 0, REAL);
        else if (size == 1 && B.Size() == 0)
            status = Dimension(0, 0, REAL);
        else
            status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        return status;
    }

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
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                (*this)(i, j) = A(i, j) + B(i, j);
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                z(i, j) = A.z(i, j) + B.z(i, j);
        }
    }
    else if (A.IsReal() && !B.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                z(i, j) = A(i, j) + B.z(i, j);
        }
    }
    else // if (!A.IsReal() && B.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                z(i, j) = A.z(i, j) + B(i, j);
        }
    }

    return status;
}

//! Add a real number to a matrix so that (*this) = A + real
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Add(const hwTMatrix<T1, T2>& A, T1 real)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;
    int count = A.Size();

    status = Dimension(A.m_nRows, A.m_nCols, A.Type());

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (A.IsReal())
    {
        T1* t_r = m_real;
        T1* a_r = A.m_real;

        while (count--)
            *t_r++ = (*a_r++) + real;
    }
    else
    {
        T2* t_c = m_complex;
        T2* a_c = A.m_complex;

        while (count--)
            *t_c++ = (*a_c++) + real;
    }

    return status;
}

//! Add a complex number to a matrix so that (*this) = A + cmplx
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Add(const hwTMatrix<T1, T2>& A, const T2& cmplx)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    int m = A.m_nRows;
    int n = A.m_nCols;
    int count = A.Size();
    hwMathStatus status = Dimension(m, n, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(0);
        return status;
    }

    T2* t_c = m_complex;

    if (A.IsReal())
    {
        T1* a_r = A.m_real;

        while (count--)
            *t_c++ = (*a_r++) + cmplx;
    }
    else
    {
        T2* a_c = A.m_complex;

        while (count--)
            *t_c++ = (*a_c++) + cmplx;
    }

    return status;
}

//! Add a matrix to the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::AddEquals(const hwTMatrix<T1, T2>& A)
{
    // check dimensions
    hwMathStatus status;
    int m = m_nRows;
    int n = m_nCols;
    int size = Size();

    if (A.m_nRows != m || A.m_nCols != n)
    {
        if (size == 0 && A.Size() == 1)
            status = Dimension(0, 0, REAL);
        else if (size == 1 && A.Size() == 0)
            status = Dimension(0, 0, REAL);
        else
            status(HW_MATH_ERR_ARRAYSIZE, 0, 1);

        return status;
    }

    if (IsReal() && !A.IsReal())
    {
        status = MakeComplex();

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }
    }

    if (IsReal() && A.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                (*this)(i, j) += A(i, j);
        }
    }
    else if (!IsReal() && !A.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                z(i, j) += A.z(i, j);
        }
    }
    else // if (!IsReal() && A.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                z(i, j) += A(i, j);
        }
    }

    return status;
}

//! Add a real number to the calling object
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::AddEquals(T1 real)
{
    int count = Size();

    if (IsReal())
    {
        T1* t_r = m_real;

        while (count--)
            *t_r++ += real;
    }
    else
    {
        T2* t_c = m_complex;

        while (count--)
            *t_c++ += real;
    }
}

//! Add a complex number to the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::AddEquals(const T2& cmplx)
{
    hwMathStatus status;
    int count = Size();

    if (IsReal())
    {
        status = MakeComplex();

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }
    }

    T2* t_c = m_complex;

    while (count--)
        *t_c++ += cmplx;

    return status;
}

//! Subtract two matrices so that (*this) = A - B
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Subtr(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B)
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

    if (A.m_nRows != B.m_nRows || A.m_nCols != B.m_nCols)
    {
        if (size == 0 && B.Size() == 1)
            status = Dimension(0, 0, REAL);
        else if (size == 1 && B.Size() == 0)
            status = Dimension(0, 0, REAL);
        else
            status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        return status;
    }

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
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                (*this)(i, j) = A(i, j) - B(i, j);
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                z(i, j) = A.z(i, j) - B.z(i, j);
        }
    }
    else if (A.IsReal() && !B.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                z(i, j) = A(i, j) - B.z(i, j);
        }
    }
    else // if (!A.IsReal() && B.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                z(i, j) = A.z(i, j) - B(i, j);
        }
    }

    return status;
}

//! Add a real number to a matrix so that (*this) = A - real
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Subtr(const hwTMatrix<T1, T2>& A, T1 real)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;
    int count = A.Size();

    status = Dimension(A.m_nRows, A.m_nCols, A.Type());

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (A.IsReal())
    {
        T1* t_r = m_real;
        T1* a_r = A.m_real;

        while (count--)
            *t_r++ = (*a_r++) - real;
    }
    else
    {
        T2* t_c = m_complex;
        T2* a_c = A.m_complex;

        while (count--)
            *t_c++ = (*a_c++) - real;
    }

    return status;
}

//! Add a complex number to a matrix so that (*this) = A - cmplx
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Subtr(const hwTMatrix<T1, T2>& A, const T2& cmplx)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    int m = A.m_nRows;
    int n = A.m_nCols;
    int count = A.Size();
    hwMathStatus status = Dimension(m, n, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(0);
        return status;
    }

    T2* t_c = m_complex;

    if (A.IsReal())
    {
        T1* a_r = A.m_real;

        while (count--)
            *t_c++ = (*a_r++) - cmplx;
    }
    else
    {
        T2* a_c = A.m_complex;

        while (count--)
            *t_c++ = (*a_c++) - cmplx;
    }

    return status;
}

//! Subtract a matrix from a real number so that (*this) = real - A
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Subtr(T1 real, const hwTMatrix<T1, T2>& A)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    int count = A.Size();
    hwMathStatus status = Dimension(A.m_nRows, A.m_nCols, A.Type());

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (A.IsReal())
    {
        T1* t_r = m_real;
        T1* a_r = A.m_real;

        while (count--)
            *t_r++ = real - (*a_r++);
    }
    else
    {
        T2* t_c = m_complex;
        T2* a_c = A.m_complex;

        while (count--)
            *t_c++ = real - (*a_c++);
    }

    return status;
}

//! Subtract a matrix from a complex number so that (*this) = cmplx - A
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Subtr(const T2& cmplx, const hwTMatrix<T1, T2>& A)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    int count = A.Size();
    hwMathStatus status = Dimension(A.m_nRows, A.m_nCols, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(0);
        return status;
    }

    T2* t_c = m_complex;

    if (A.IsReal())
    {
        T1* a_r = A.m_real;

        while (count--)
            *t_c++ = cmplx - (*a_r++);
    }
    else
    {
        T2* a_c = A.m_complex;

        while (count--)
            *t_c++ = cmplx - (*a_c++);
    }

    return status;
}

//! Negate a matrix
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Negate(const hwTMatrix<T1, T2>& A)
{
    int count = A.Size();
    hwMathStatus status = Dimension(A.m_nRows, A.m_nCols, A.Type());

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (A.IsReal())
    {
        T1* t_r = m_real;
        T1* a_r = A.m_real;

        while (count--)
            *t_r++ = -(*a_r++);
    }
    else
    {
        T2* t_c = m_complex;
        T2* a_c = A.m_complex;

        while (count--)
            *t_c++ = -(*a_c++);
    }

    return status;
}

//! Subtract a matrix from the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::SubtrEquals(const hwTMatrix<T1, T2>& A)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;
    int m = m_nRows;
    int n = m_nCols;
    int size = Size();

    if (A.m_nRows != m || A.m_nCols != n)
    {
        if (size == 0 && A.Size() == 1)
            status = Dimension(0, 0, REAL);
        else if (size == 1 && A.Size() == 0)
            status = Dimension(0, 0, REAL);
        else
            status(HW_MATH_ERR_ARRAYSIZE, 0, 1);

        return status;
    }

    if (IsReal() && !A.IsReal())
    {
        status = MakeComplex();

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }
    }

    if (IsReal() && A.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                (*this)(i, j) -= A(i, j);
        }
    }
    else if (!IsReal() && !A.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                z(i, j) -= A.z(i, j);
        }
    }
    else // if (!IsReal() && A.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
                z(i, j) -= A(i, j);
        }
    }

    return status;
}

//! Subtract a real number from the calling object
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::SubtrEquals(T1 real)
{
    int count = Size();

    if (IsReal())
    {
        T1* t_r = m_real;

        while (count--)
            *t_r++ -= real;
    }
    else
    {
        T2* t_c = m_complex;

        while (count--)
            *t_c++ -= real;
    }
}

//! Subtract a complex number from the calling object
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::SubtrEquals(const T2& cmplx)
{
    hwMathStatus status;
    int count = Size();

    if (IsReal())
    {
        status = MakeComplex();

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }
    }

    T2* t_c = m_complex;

    while (count--)
        *t_c++ -= cmplx;

    return status;
}

//! Multiply two matrices so that (*this) = A * B
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Mult(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;
    int m = A.m_nRows;
    int n = B.m_nCols;
    int bm = B.m_nRows;

    if (A.m_nCols != bm)
    {
        if (A.Size() == 1)
        {
            if (A.IsReal())
                return Mult(B, A(0));
            else
                return Mult(B, A.z(0));
        }
        else if (B.Size() == 1)
        {
            if (B.IsReal())
                return Mult(A, B(0));
            else
                return Mult(A, B.z(0));
        }
        else
            status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        return status;
    }

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
        T1 temp_r;
        SetElements((T1) 0);

        // optimized for column storage
        for (int j = 0; j < n; j++)
        {
            for (int k = 0; k < bm; k++)
            {
                temp_r = B(k, j);

                for (int i = 0; i < m; i++)
                    (*this)(i, j) += A(i, k) * temp_r;
            }
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                z(i, j).Set((T1) 0, (T1) 0);

                for (int k = 0; k < bm; ++k)
                    z(i, j) += A.z(i, k) * B.z(k, j);
            }
        }
    }
    else if (A.IsReal() && !B.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                z(i, j).Set((T1) 0, (T1) 0);

                for (int k = 0; k < bm; ++k)
                    z(i, j) += B.z(k, j) * A(i, k);
            }
        }
    }
    else // if (!A.IsReal() && B.IsReal())
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                z(i, j).Set((T1) 0, (T1) 0);

                for (int k = 0; k < bm; ++k)
                    z(i, j) += A.z(i, k) * B(k, j);
            }
        }
    }

    return status;
}       

//! Multiply a matrix and a real number so that (*this) = A * real
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Mult(const hwTMatrix<T1, T2>& A, T1 real)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    int size = A.Size();
    hwMathStatus status = Dimension(A.m_nRows, A.m_nCols, A.Type());

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (A.IsReal())
    {
        for (int i = 0; i < size; ++i)
            (*this)(i) = A(i) * real;
    }
    else
    {
        for (int i = 0; i < size; ++i)
            z(i) = A.z(i) * real;
    }

    return status;
}

//! Multiply a matrix and a complex number so that (*this) = A * cmplx
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Mult(const hwTMatrix<T1, T2>& A, const T2& cmplx)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    int size = A.Size();
    hwMathStatus status = Dimension(A.m_nRows, A.m_nCols, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(0);
        return status;
    }

    if (A.IsReal())
    {
        for (int i = 0; i < size; ++i)
            z(i) = A(i) * cmplx;
    }
    else
    {
        for (int i = 0; i < size; ++i)
            z(i) = A.z(i) * cmplx;
    }

    return status;
}

//! Multiply the calling object by a matrix
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::MultEquals(const hwTMatrix<T1, T2>& A)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;

    if (A.m_nRows != A.m_nCols)
        return status(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (A.m_nRows != m_nCols)
        return status(HW_MATH_ERR_ARRAYSIZE, 0, 1);

    hwTMatrix<T1, T2> C(*this);

    status = Mult(C, A);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.ResetArgs();

        return status;
    }

    return status;
}

//! Multiply the calling object by a real number
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::MultEquals(T1 real)
{
    int count = Size();

    if (IsReal())
    {
        T1* t_r = m_real;

        while (count--)
            *t_r++ *= real;
    }
    else
    {
        T2* t_c = m_complex;

        while (count--)
            *t_c++ *= real;
    }
}

//! Multiply the calling object by a complex number
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::MultEquals(const T2& cmplx)
{
    hwMathStatus status;
    int count = Size();

    if (IsReal())
    {
        status = MakeComplex();

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }
    }

    T2* t_c = m_complex;

    while (count--)
        *t_c++ *= cmplx;

    return status;
}

//! Divide a matrix by a real number so that (*this) = A / real
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Divide(const hwTMatrix<T1, T2>& A, T1 real)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    int size = A.Size();
    hwMathStatus status = Dimension(A.m_nRows, A.m_nCols, A.Type());

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (A.IsReal())
    {
        for (int i = 0; i < size; ++i)
            (*this)(i) = A(i) / real;
    }
    else
    {
        for (int i = 0; i < size; ++i)
            z(i) = A.z(i) / real;
    }

    return status;
}

//! Divide a matrix by a complex number so that (*this) = A / cmplx
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Divide(const hwTMatrix<T1, T2>& A, const T2& cmplx)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    int size = A.Size();
    hwMathStatus status = Dimension(A.m_nRows, A.m_nCols, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(0);
        return status;
    }

    if (A.IsReal())
    {
        for (int i = 0; i < size; ++i)
            z(i) = A(i) / cmplx;
    }
    else
    {
        for (int i = 0; i < size; ++i)
            z(i) = A.z(i) / cmplx;
    }

    return status;
}

//! Divide a real number by each element of a matrix
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Divide(T1 real, const hwTMatrix<T1, T2>& A)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;
    int size = A.Size();

    status = Dimension(A.m_nRows, A.m_nCols, A.Type());

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (A.IsReal())
    {
        for (int i = 0; i < size; ++i)
            (*this)(i) = real / A(i);
    }
    else
    {
        T2 c(real, (T1) 0);

        for (int i = 0; i < size; ++i)
            z(i) = c / A.z(i);
    }

    return status;
}

//! Divide a complex number by each element of a matrix
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Divide(const T2& cmplx, const hwTMatrix<T1, T2>& A)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;
    int size = A.Size();

    status = Dimension(A.m_nRows, A.m_nCols, COMPLEX);

    if (!status.IsOk())
    {
        status.SetArg1(0);
        return status;
    }

    if (A.IsReal())
    {
        for (int i = 0; i < size; ++i)
            z(i) = cmplx / A(i);
    }
    else
    {
        for (int i = 0; i < size; ++i)
            z(i) = cmplx / A.z(i);
    }

    return status;
}

//! Divide the calling object by a real number
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::DivideEquals(T1 real)
{
    int count = Size();

    if (IsReal())
    {
        T1* t_r = m_real;

        while (count--)
            *t_r++ /= real;
    }
    else
    {
        T2* t_c = m_complex;

        while (count--)
            *t_c++ /= real;
    }
}

//! Divide the calling object by a complex number
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::DivideEquals(const T2& cmplx)
{
    hwMathStatus status;
    int count = Size();

    if (IsReal())
    {
        status = MakeComplex();

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }
    }

    T2* t_c = m_complex;

    while (count--)
        *t_c++ /= cmplx;

    return status;
}

//! Multiply two matrices entry by entry so that (*this) = A .* B
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::MultByElems(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<T1, T2>& C = (*this);
    hwMathStatus status;
    int count = A.Size();
    int m = A.m_nRows;
    int n = A.m_nCols;

    if (B.m_nRows != m || B.m_nCols != n)
    {
        if (count == 0 && B.Size() == 1)
            status = C.Dimension(0, 0, REAL);
        else if (count == 1 && B.Size() == 0)
            status = C.Dimension(0, 0, REAL);
        else
            status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        return status;
    }

    if (A.IsReal() && B.IsReal())
    {
        status = C.Dimension(m, n, REAL);

        if (!status.IsOk())
            return status;

        T1* a_r = A.m_real;
        T1* b_r = B.m_real;
        T1* c_r = C.m_real;

        // note: This could be implemented with DSBMV, which may be faster with parallelization
        while (count--)
            *c_r++ = (*a_r++) * (*b_r++);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
            return status;

        T2* a_c = A.m_complex;
        T2* b_c = B.m_complex;
        T2* c_c = C.m_complex;

        while (count--)
            *c_c++ = (*a_c++) * (*b_c++);
    }
    else if (A.IsReal() && !B.IsReal())
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(2);
            else
                status.ResetArgs();

            return status;
        }

        T1* a_r = A.m_real;
        T2* b_c = B.m_complex;
        T2* c_c = C.m_complex;

        while (count--)
            *c_c++ = (*b_c++) * (*a_r++);
    }
    else // if (!A.IsReal() && B.IsReal())
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(2);
            else
                status.ResetArgs();

            return status;
        }

        T2* a_c = A.m_complex;
        T1* b_r = B.m_real;
        T2* c_c = C.m_complex;

        while (count--)
            *c_c++ = (*a_c++) * (*b_r++);
    }

    return status;
}

//! Divide two matrices entry by entry so that (*this) = A ./ B
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::DivideByElems(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<T1, T2>& C = (*this);
    hwMathStatus status;
    int count = A.Size();
    int m = A.m_nRows;
    int n = A.m_nCols;

    if (B.m_nRows != m || B.m_nCols != n)
    {
        if (A.Size() == 1)
        {
            if (A.IsReal())
                return Divide(A(0), B);
            else
                return Divide(A.z(0), B);
        }
        else if (B.Size() == 1)
        {
            if (B.IsReal())
                return Divide(A, B(0));
            else
                return Divide(A, B.z(0));
        }
        else
            return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    if (A.IsReal() && B.IsReal())
    {
        status = C.Dimension(m, n, REAL);

        if (!status.IsOk())
            return status;

        T1* a_r = A.m_real;
        T1* b_r = B.m_real;
        T1* c_r = C.m_real;

        while (count--)
            *c_r++ = (*a_r++) / (*b_r++);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
            return status;

        T2* a_c = A.m_complex;
        T2* b_c = B.m_complex;
        T2* c_c = C.m_complex;

        while (count--)
            *c_c++ = (*a_c++) / (*b_c++);
    }
    else if (A.IsReal() && !B.IsReal())
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(2);
            else
                status.ResetArgs();

            return status;
        }

        T1* a_r = A.m_real;
        T2* b_c = B.m_complex;
        T2* c_c = C.m_complex;

        while (count--)
            *c_c++ = (*a_r++) / (*b_c++);
    }
    else // if (!A.IsReal() && B.IsReal())
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(2);
            else
                status.ResetArgs();

            return status;
        }

        T2* a_c = A.m_complex;
        T1* b_r = B.m_real;
        T2* c_c = C.m_complex;

        while (count--)
            *c_c++ = (*a_c++) / (*b_r++);
    }

    return status;
}

//! Raise each object matrix entry to the specified real power so that (*this) = B.^power
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::PowerByElems(const hwTMatrix<T1, T2>& Base, T1 power)
{
    if (this == &Base)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<T1, T2>& C = (*this);
    hwMathStatus status;
    int m = Base.m_nRows;
    int n = Base.m_nCols;
    int count = -1;
    hwTMatrix<T1, T2>* B = (hwTMatrix<T1, T2>*) &Base;

    // check for negative real bases to non-integer exponents
    if (Base.IsReal() && !IsInteger(power).IsOk())
    {
        T1* b_r = Base.m_real;
        count = Base.Size();

        while (count--)
        {
            if (*b_r < 0.0)
                break;

            ++b_r;
        }
    }

    if (count != -1)
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
            return status;

        int index = 0;
        T1* b_r = Base.m_real;
        count = Size();

        while (index < count)
        {
            C.z(index) = hwTComplex<T1>::pow_c(*b_r++, power);
            ++index;
        }
    }
    else if (Base.IsReal())
    {
        status = C.Dimension(m, n, REAL);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(3);
            else
                status.ResetArgs();

            return status;
        }

        T1* b_r = B->m_real;
        T1* c_r = C.m_real;
        count = Size();

        while (count--)
        {
            *c_r = CustomPow(*b_r, power);
            ++c_r;
            ++b_r;
        }
    }
    else
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(3);
            else
                status.ResetArgs();

            return status;
        }

        count = Size();

        while (count--)
            C.z(count) = hwTComplex<T1>::pow(Base.z(count), power);
    }

    return status;
}

//! Raise each object matrix entry to the specified complex power so that (*this) = B.^power
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::PowerByElems(const hwTMatrix<T1, T2>& Base, T2 power)
{
    if (this == &Base)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<T1, T2>& C = (*this);
    int m = Base.m_nRows;
    int n = Base.m_nCols;
    hwTMatrix<T1, T2>* B = (hwTMatrix<T1, T2>*) &Base;

    hwMathStatus status = C.Dimension(m, n, COMPLEX);

    int count = Size();

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (Base.IsReal())
    {
        T1* b_r = B->m_real;
        T2* c_c = C.m_complex;

        while (count--)
        {
            *c_c = hwTComplex<T1>::pow(*b_r, power);
            ++c_c;
            ++b_r;
        }
    }
    else
    {
        while (count--)
            C.z(count) = hwTComplex<T1>::pow(Base.z(count), power);
    }

    return status;
}

//! Raise each object matrix entry to the power of the argument
//! matrix entry so that (*this) = B.^Power
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::PowerByElems(const hwTMatrix<T1, T2>& Base,
                                             const hwTMatrix<T1, T2>& Pow)
{
    if (this == &Base)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &Pow)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    if (Base.m_nRows != Pow.m_nRows || Base.m_nCols != Pow.m_nCols)
    {
        if (Base.Size() == 1)
        {
            if (Base.IsReal())
                return PowerByElems(Base(0), Pow);
            else
                return PowerByElems(Base.z(0), Pow);
        }
        else if (Pow.Size() == 1)
        {
            if (Pow.IsReal())
                return PowerByElems(Base, Pow(0));
            else
                return PowerByElems(Base, Pow.z(0));
        }
        else
            return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    hwTMatrix<T1, T2>& C = (*this);
    hwMathStatus status;
    int m = Base.m_nRows;
    int n = Base.m_nCols;
    int count = -1;

    // check for negative bases to real non-integer exponents
    if (Base.IsReal() && Pow.IsReal())
    {
        const T1* b_r = Base.GetRealData();
        const T1* p_r = Pow.GetRealData();
        count = Base.Size();

        while (count--)
        {
            if (*b_r < 0.0 && !IsInteger(*p_r).IsOk())
                break;

            ++b_r;
            ++p_r;
        }
    }

    if (count != -1)
    {
        status = Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(2);
            else
                status.ResetArgs();

            return status;
        }

        int index = 0;
        const T1* b_r = Base.GetRealData();
        const T1* p_r = Pow.GetRealData();
        count = Size();

        while (index < count)
        {
            z(index) = hwTComplex<T1>::pow_c(*b_r++, *p_r++);
            ++index;
        }
    }
    else if (Base.IsReal() && Pow.IsReal())
    {
        status = C.Dimension(m, n, REAL);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(3);
            else
                status.ResetArgs();

            return status;
        }

        T1* b_r = Base.m_real;
        T1* p_r = Pow.m_real;
        T1* c_r = C.m_real;
        count = Size();

        while (count--)
            *c_r++ = CustomPow((*b_r++), (*p_r++));
    }
    else if (!Base.IsReal() && !Pow.IsReal())
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(3);
            else
                status.ResetArgs();

            return status;
        }

        int index = 0;
        count = Size();

        while (index < count)
        {
            C.z(index) = hwTComplex<T1>::pow(Base.z(index), Pow.z(index));
            ++index;
        }
    }
    else if (Base.IsReal() && !Pow.IsReal())
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(3);
            else
                status.ResetArgs();

            return status;
        }

        int index = 0;
        T1* b_r = Base.m_real;
        count = Size();

        while (index < count)
        {
            C.z(index) = hwTComplex<T1>::pow(*b_r++, Pow.z(index));
            ++index;
        }
    }
    else // if (!Base.IsReal() && Pow.IsReal())
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(3);
            else
                status.ResetArgs();

            return status;
        }

        int index = 0;
        T1* p_r = Pow.m_real;
        count = Size();

        while (index < count)
        {
            C.z(index) = hwTComplex<T1>::pow(Base.z(index), *p_r++);
            ++index;
        }
    }

    return status;
}

//! Raise a real value to the exponent specified by each entry in a
//! matrix so that (*this) = base.^Pow
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::PowerByElems(T1 base, const hwTMatrix<T1, T2>& Pow)
{
    if (this == &Pow)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<T1, T2>& C = (*this);
    hwMathStatus status;
    int m = Pow.m_nRows;
    int n = Pow.m_nCols;
    int count = -1;

    // check for negative base to real non-integer exponents
    if (base < 0.0 && Pow.IsReal())
    {
        T1* p_r = Pow.m_real;
        count = Pow.Size();

        while (count--)
        {
            if (!IsInteger(*p_r++).IsOk())
                break;
        }
    }

    if (count != -1)
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
            return status;

        int index = 0;
        T1* p_r = Pow.m_real;
        count = Size();

        while (index < count)
        {
            C.z(index) = hwTComplex<T1>::pow_c(base, *p_r++);
            ++index;
        }
    }
    else if (Pow.IsReal())
    {
        status = C.Dimension(m, n, REAL);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(3);
            else
                status.ResetArgs();

            return status;
        }

        T1* p_r = Pow.m_real;
        T1* c_r = C.m_real;
        count = Size();

        while (count--)
            *c_r++ = CustomPow(base, (*p_r++));
    }
    else // !Power.IsReal()
    {
        status = C.Dimension(m, n, COMPLEX);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 0)
                status.SetArg1(3);
            else
                status.ResetArgs();

            return status;
        }

        int index = 0;
        count = Size();

        while (index < count)
        {
            C.z(index) = hwTComplex<T1>::pow(base, Pow.z(index));
            ++index;
        }
    }

    return status;
}

//! Raise a complex value to the exponent specified by each entry in
//! the object matrix so that (*this) = base.^A
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::PowerByElems(T2 base, const hwTMatrix<T1, T2>& Pow)
{
    if (this == &Pow)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<T1, T2>& C = (*this);
    int m = Pow.m_nRows;
    int n = Pow.m_nCols;
    int count = Pow.Size();
    int index = 0;

    hwMathStatus status = C.Dimension(m, n, COMPLEX);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (Pow.IsReal())
    {
        T1* p_r = Pow.m_real;

        while (index < count)
        {
            C.z(index) = hwTComplex<T1>::pow(base, *p_r++);
            ++index;
        }
    }
    else // !Power.IsReal()
    {
        while (index < count)
        {
            C.z(index) = hwTComplex<T1>::pow(base, Pow.z(index));
            ++index;
        }
    }

    return status;
}

// ****************************************************
//               Arithmetic Operators
// ****************************************************

//! Implement the == operator
template<typename T1, typename T2>
bool hwTMatrix<T1, T2>::operator==(const hwTMatrix<T1, T2>& A) const
{
    return IsEqual(A, (T1) 0);
}

//! Determine if the matrix is equal to another
template<typename T1, typename T2>
bool hwTMatrix<T1, T2>::IsEqual(const hwTMatrix<T1, T2>& A, T1 tol) const
{
    int m = m_nRows;
    int n = m_nCols;

    if (A.m_nRows != m || A.m_nCols != n)
        return false;

    int size = Size();

    if (m_real)
    {
        if (A.m_real)
        {
            for (int i = 0; i < size; ++i)
            {
                if (!AreEqual(m_real[i], A.m_real[i], tol))
                    return false;
            }
        }
        else if (A.m_complex)
        {
            for (int i = 0; i < size; ++i)
            {
                if (!A.m_complex[i].IsReal(tol))
                    return false;

                if (!A.m_complex[i].IsEqual(m_real[i], tol))
                    return false;
            }
        }
    }
    else if (A.m_complex)
    {
        for (int i = 0; i < size; ++i)
        {
            if (!m_complex[i].IsEqual(A.m_complex[i], tol))
                return false;
        }
    }
    else if (A.m_real)
    {
        for (int i = 0; i < size; ++i)
        {
            if (!m_complex[i].IsReal(tol))
                return false;

            if (!m_complex[i].IsEqual(A.m_real[i], tol))
                return false;
        }
    }

    return true;
}

//! Implement the != operator
template<typename T1, typename T2>
bool hwTMatrix<T1, T2>::operator!=(const hwTMatrix<T1, T2>& A) const
{
    return (!(*this == A));
}

//! Implement the + operator with a matrix argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator+(const hwTMatrix<T1, T2>& A) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Add(*this, A);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the + operator with a real argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator+(T1 real) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Add(*this, real);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the + operator with a complex argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator+(const T2& cmplx) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Add(*this, cmplx);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the += operator with a matrix argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator+=(const hwTMatrix<T1, T2>& A)
{
    hwMathStatus status = AddEquals(A);

    if (!status.IsOk())
        MakeEmpty();

    return *this;
}

//! Implement the += operator with a real argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator+=(T1 real)
{
    AddEquals(real);
    return *this;
}

//! Implement the += operator with a complex argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator+=(const T2& cmplx)
{
    hwMathStatus status = AddEquals(cmplx);

    if (!status.IsOk())
        MakeEmpty();

    return *this;
}

//! Implement the - operator with a matrix argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator-(const hwTMatrix<T1, T2>& A) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Subtr(*this, A);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the - operator with a real argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator-(T1 real) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Subtr(*this, real);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the - operator with a complex argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator-(const T2& cmplx) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Subtr(*this, cmplx);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the - operator as a urnary prefix
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator-() const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Negate(*this);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the -= operator with a matrix argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator-=(const hwTMatrix<T1, T2>& A)
{
    hwMathStatus status = SubtrEquals(A);

    if (!status.IsOk())
        MakeEmpty();

    return *this;
}

//! Implement the -= operator with a real argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator-=(T1 real)
{
    SubtrEquals(real);
    return *this;
}

//! Implement the -= operator with a complex argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator-=(const T2& cmplx)
{
    hwMathStatus status = SubtrEquals(cmplx);

    if (!status.IsOk())
        MakeEmpty();

    return *this;
}

//! Implement the * operator with a matrix argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator*(const hwTMatrix<T1, T2>& A) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Mult(*this, A);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the * operator with a real argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator*(T1 real) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Mult(*this, real);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the * operator with a complex argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator*(const T2& cmplx) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Mult(*this, cmplx);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the *= operator with a matrix argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator*=(const hwTMatrix<T1, T2>& B)
{
    hwMathStatus status = MultEquals(B);

    if (!status.IsOk())
        MakeEmpty();

    return *this;
}

//! Implement the *= operator with a real argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator*=(T1 real)
{
    MultEquals(real);
    return *this;
}

//! Implement the *= operator with a complex argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator*=(const T2& cmplx)
{
    hwMathStatus status = MultEquals(cmplx);

    if (!status.IsOk())
        MakeEmpty();

    return *this;
}

//! Implement the / operator with a matrix argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator/(const hwTMatrix<T1, T2>& B) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.DivideRight(*this, B);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the / operator with a real argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator/(T1 real) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Divide(*this, real);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the / operator with a complex argument
template<typename T1, typename T2>
hwTMatrix<T1, T2> hwTMatrix<T1, T2>::operator/(const T2& cmplx) const
{
    hwTMatrix<T1, T2> result;
    hwMathStatus status = result.Divide(*this, cmplx);

    if (!status.IsOk())
        result.MakeEmpty();

    return result;
}

//! Implement the /= operator with a matrix argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator/=(const hwTMatrix<T1, T2>& A)
{
    hwMathStatus status = DivideEquals(A);

    if (!status.IsOk())
        MakeEmpty();

    return *this;
}

//! Implement the /= operator with a real argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator/=(T1 real)
{
    DivideEquals(real);
    return *this;
}

//! Implement the /= operator with a complex argument
template<typename T1, typename T2>
hwTMatrix<T1, T2>& hwTMatrix<T1, T2>::operator/=(const T2& cmplx)
{
    hwMathStatus status = DivideEquals(cmplx);

    if (!status.IsOk())
        MakeEmpty();

    return *this;
}

// ****************************************************
//                 Vector Operations
// ****************************************************

//! Square of L2 vector norm
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::L2NormSq(T1& normSq) const
{
    if (!IsVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 0);

    int count = Size();

    normSq = 0.0;

    if (IsReal())
    {
        T1* a_r = m_real;

        while (count--)
            normSq += (*a_r) * (*a_r++);
    }
    else    // complex
    {
        T2* a_c = m_complex;

        while (count--)
            normSq += (a_c++)->MagSq();
    }

    return hwMathStatus();
}

//! Dot product of two vectors
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Dot(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B, T1& dot)
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
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 0, 1);

    T1* a_r = A.m_real;
    T1* b_r = B.m_real;

    dot = (T1) 0;

    while (count--)
    {
        dot += (*a_r) * (*b_r);
        ++a_r;
        ++b_r;
    }

    return hwMathStatus();
}

//! Kronecker product of two matrices
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Kronecker(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B)
{
    hwMathStatus status;
    int am = A.m_nRows;
    int an = A.m_nCols;
    int bm = B.m_nRows;
    int bn = B.m_nCols;
    hwTMatrix<T1, T2> block;

    if (A.IsReal())
    {
        status = Dimension(am*bm, an*bn, B.Type());

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }

        for (int j = 0; j < an; ++j)
        {
            for (int i = 0; i < am; ++i)
            {
                block = A(i,j) * B;
                status = WriteSubmatrix(i*bm, j*bn, block);
            }
        }
    }
    else
    {
        status = Dimension(am*bm, an*bn, COMPLEX);

        if (!status.IsOk())
        {
            status.SetArg1(0);
            return status;
        }

        for (int j = 0; j < an; ++j)
        {
            for (int i = 0; i < am; ++i)
            {
                block = A.z(i,j) * B;
                status = WriteSubmatrix(i*bm, j*bn, block);
            }
        }
    }

    return status;
}

//! Cross product of two vectors
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Cross(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwTMatrix<T1, T2>& cross = (*this);
    hwMathStatus status;

    if (A.m_complex)
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);

    if (A.Size() != 3)
        return status(HW_MATH_ERR_ARRAYSIZE, 1);

    if (B.m_complex)
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);

    if (B.Size() != 3)
        return status(HW_MATH_ERR_ARRAYSIZE, 2);

    cross.Dimension(3, REAL);

    if (!status.IsOk())
        return status;

    cross(0) = A(1)*B(2) - B(1)*A(2);
    cross(1) = B(0)*A(2) - A(0)*B(2);
    cross(2) = A(0)*B(1) - B(0)*A(1);

    return status;
}

//! Linear convolution of two vectors
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::ConvLin(const hwTMatrix<T1, T2>& X, const hwTMatrix<T1, T2>& Y)
{
    if (this == &X)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &Y)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    // X is reversed and passed over Y
    hwMathStatus status;

    if (X.m_nRows != 1 && X.m_nCols != 1)  // allow empty vector
        return status(HW_MATH_ERR_VECTOR, 1);

    if (Y.m_nRows != 1 && Y.m_nCols != 1)  // allow empty vector
        return status(HW_MATH_ERR_VECTOR, 2);

    int xn = X.Size();
    int yn = Y.Size();
    int cn = xn + yn - 1;			// size of conv
    int partial = _min(xn, yn) - 1; // size of partial overlaps
    int count;

    if (X.IsReal() && Y.IsReal())
    {
        if (X.m_nRows == 1)
            status = Dimension(1, cn, REAL);
        else
            status = Dimension(cn, 1, REAL);

        if (!cn)
            return status;

        if (!status.IsOk())
        {
            if (status.GetArg1() != 0)
                status.ResetArgs();

            return status;
        }

        T1 value;
        T1* x_start = X.m_real;
        T1* y_start = Y.m_real;
        T1* x;
        T1* y;
        T1* c = m_real;

        // leading partial overlap
        for (int i = 0; i < partial; ++i)
        {
            x = x_start + i;
            y = y_start;
            count = i + 1;
            value = (T1) 0;

            while (count--)
                value += (*x--) * (*y++);

            c[i] = value;
        }

        // complete overlap
        if (xn < yn)
        {
            for (int i = partial; i < yn; ++i)
            {
                x = x_start + partial;
                y = y_start - partial + i;
                count = xn;
                value = (T1) 0;

                while (count--)
                    value += (*x--) * (*y++);

                c[i] = value;
            }
        }
        else
        {
            for (int i = partial; i < xn; ++i)
            {
                x = x_start + i;
                y = y_start;
                count = yn;
                value = (T1) 0;

                while (count--)
                    value += (*x--) * (*y++);

                c[i] = value;
            }
        }

        // trailing partial overlap
        for (int i = cn - partial; i < cn; ++i)
        {
            x = x_start + xn - 1;
            y = y_start - xn + 1 + i;
            count = cn - i;
            value = (T1) 0;

            while (count--)
                value += (*x--) * (*y++);

            c[i] = value;
        }
    }
    else if (X.m_complex && Y.m_complex)
    {
        if (X.m_nRows == 1)
            status = Dimension(1, cn, COMPLEX);
        else
            status = Dimension(cn, 1, COMPLEX);

        if (!cn)
            return status;

        if (!status.IsOk())
        {
            if (status.GetArg1() != 0)
                status.ResetArgs();

            return status;
        }

        T2 value;
        T2* x_start = X.m_complex;
        T2* y_start = Y.m_complex;
        T2* x;
        T2* y;
        T2* c = m_complex;

        // leading partial overlap
        for (int i = 0; i < partial; ++i)
        {
            x = x_start + i;
            y = y_start;
            count = i + 1;
            value = hwTComplex<T1>(0.0, 0.0);

            while (count--)
                value += (*x--) * (*y++);

            c[i] = value;
        }

        // complete overlap
        if (xn < yn)
        {
            for (int i = partial; i < yn; ++i)
            {
                x = x_start + partial;
                y = y_start - partial + i;
                count = xn;
                value = hwTComplex<T1>(0.0, 0.0);

                while (count--)
                    value += (*x--) * (*y++);

                c[i] = value;
            }
        }
        else
        {
            for (int i = partial; i < xn; ++i)
            {
                x = x_start + i;
                y = y_start;
                count = yn;
                value = hwTComplex<T1>(0.0, 0.0);

                while (count--)
                    value += (*x--) * (*y++);

                c[i] = value;
            }
        }

        // trailing partial overlap
        for (int i = cn - partial; i < cn; ++i)
        {
            x = x_start + xn - 1;
            y = y_start - xn + 1 + i;
            count = cn - i;
            value = hwTComplex<T1>(0.0, 0.0);

            while (count--)
                value += (*x--) * (*y++);

            c[i] = value;
        }
    }
    else if (X.m_real)
    {
        hwTMatrix<T1, T2> temp;
        temp.PackComplex(X, nullptr);

        return ConvLin(temp, Y);
    }
    else if (Y.m_real)
    {
        hwTMatrix<T1, T2> temp;
        temp.PackComplex(Y, nullptr);

        return ConvLin(X, temp);
    }

    return status;
}

//! Linear correlation of two vectors
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::CorrLin(const hwTMatrix<T1, T2>& X, const hwTMatrix<T1, T2>& Y)
{
    if (this == &X)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &Y)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    // linear correlation - time domain
    // X is passed over Y
    hwMathStatus status;

    if (X.m_complex)
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 1);

    if (!X.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (Y.m_complex)
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 2);

    if (!Y.IsVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    int xn = X.Size();
    int yn = Y.Size();
    int cn = xn + yn - 1;       // size of corr

    if (X.m_nRows == 1)
        status = Dimension(1, cn, REAL);
    else
        status = Dimension(cn, 1, REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() != 0)
            status.ResetArgs();

        return status;
    }

    int partial = _min(xn, yn) - 1; // size of partial overlaps
    int count;
    T1 value;
    T1* x_start = X.m_real;
    T1* y_start = Y.m_real;
    T1* x;
    T1* y;
    T1* c = m_real;

    // leading partial overlap
    for (int i = 0; i < partial; ++i)
    {
        x = x_start;
        y = y_start + yn - 1 - i;
        count = i + 1;
        value = (T1) 0;

        while (count--)
            value += (*x++) * (*y++);

        c[i] = value;
    }

    // complete overlap
    if (xn < yn)
    {
        for (int i = partial; i < yn; ++i)
        {
            x = x_start;
            y = y_start + yn - 1 - i;
            count = xn;
            value = (T1) 0;

            while (count--)
                value += (*x++) * (*y++);

            c[i] = value;
        }
    }
    else
    {
        for (int i = partial; i < xn; ++i)
        {
            x = x_start - yn + 1 + i;
            y = y_start;
            count = yn;
            value = (T1) 0;

            while (count--)
                value += (*x++) * (*y++);

            c[i] = value;
        }
    }

    // trailing partial overlap
    for (int i = cn - partial; i < cn; ++i)
    {
        x = x_start - yn + 1 + i;
        y = y_start;
        count = cn - i;
        value = (T1) 0;

        while (count--)
            value += (*x++) * (*y++);

        c[i] = value;
    }

    return status;
}

//! Linear correlation of matrix columns in the time domain
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::CorrLin(const hwTMatrix<T1, T2>& X)
{
    hwMathStatus status;

    if (X.m_complex)
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 1);

    if (X.IsVector())
        return status(HW_MATH_ERR_NOVECTOR, 1);

    int m = X.m_nRows;
    int n = X.m_nCols;
    const double* col1;
    const double* col2;
    const double* colr;

    if (m == 0)
    {
        status = Dimension(0, n*n, REAL);
        return status;
    }

    status = Dimension(2*m-1, n*n, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    for (int i = 0; i < n; ++i)
    {
        col1 = &X(0, i);
        hwTMatrix<T1, T2> vec1(m, 1, (void*) col1, REAL);

        for (int j = 0; j < n; ++j)
        {
            col2 = &X(0, j);
            colr = &(*this)(0, n*i+j);
            hwTMatrix<T1, T2> vec2(m, 1, (void*) col2, REAL);
            hwTMatrix<T1, T2> vecr(2*m-1, 1, (void*) colr, REAL);

            status = vecr.CorrLin(vec1, vec2);

            if (!status.IsOk())
            {
                if (status.GetArg1() != 0)
                    status.ResetArgs();

                return status;
            }
        }
    }

    return status;
}

//! 2D convolution of two matrices
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Conv2D(const hwTMatrix<T1, T2>& X,
                                       const hwTMatrix<T1, T2>& Y,
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

    // find each conv(i, j)
    int start_xi;                       // X column indexing
    int start_xj, index_xj;             // X row indexing
    int start_yi, stop_yi;              // Y column window
    int start_yj, stop_yj;              // Y row window
    int partial_i = _min(xm, ym) - 1;   // size of partial col overlaps
    int partial_j = _min(xn, yn) - 1;   // size of partial row overlaps
    T1 value;

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
            value = (T1) 0;
            index_xj = start_xj;

            for (int jj = start_yj; jj < stop_yj; ++jj)
            {
                const T1* xp = &(X(start_xi, index_xj));
                const T1* yp = &(Y(start_yi, jj));

                for (int ii = start_yi; ii < stop_yi; ++ii)
                {
                    value += (*xp--) * (*yp++);
                }

                --index_xj;
            }

            (*this)(i - row1, j - col1) = value;
        }
    }

    return status;
}

//! 2D convolution of a matrix with a column vector and a row vector
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Conv2D(const hwTMatrix<T1, T2>& col, const hwTMatrix<T1, T2>& row,
                                       const hwTMatrix<T1, T2>& X,
                                       int row1, int col1, int row2, int col2)
{
    if (this == &col)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &row)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &X)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    // convolve each column of X with vector 'col', and then
    // convolve each row of the result with the  vector 'row'.
    // if 'col' is a column vector and 'row' is a row vector,
    // this is the same as conv2(col*row,X).

    // 'col' is treated as a column vector, regardless of orientation
    // 'row' is treated as a row vector, regardless of orientation
    hwMathStatus status;

    if (!col.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (col.m_complex)
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 1);

    if (!row.IsVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    if (row.m_complex)
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 2);

    if (X.IsEmpty())
        return Dimension(0, 0, REAL);

    if (X.m_complex)
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 3);

    int xm = X.m_nRows;
    int xn = X.m_nCols;
    int colSize = col.Size();
    int rowSize = row.Size();
    int cm = xm + colSize - 1;       // # rows in conv
    int cn = xn + rowSize - 1;       // # cols in conv

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

        return status;
    }

    // find each conv(i, j)
    int start_xi, stop_xi;                      // matrix column window
    int start_xj, stop_xj;                      // matrix row window
    int start_col;                              // column vector indexing
    int start_row, index_row;                   // row vector indexing
    int partial_col = _min(xm, colSize) - 1;    // size of partial col overlaps
    int partial_row = _min(xn, rowSize) - 1;    // size of partial row overlaps
    T1 value;
    hwTMatrix<T1, T2> col_temp(xm, 1, REAL);

    for (int j = col1; j <= col2; ++j)
    {
        if (j < partial_row)
        {
            start_xj = 0;
            stop_xj = j + 1;
            start_row = j;
        }
        else if (j < cn - partial_row)
        {
            if (xn < rowSize)
            {
                start_xj = 0;
                stop_xj = xn;
                start_row = j;
            }
            else
            {
                start_xj = j - rowSize + 1;
                stop_xj = j + 1;
                start_row = rowSize - 1;
            }
        }
        else
        {
            start_xj = j - rowSize + 1;
            stop_xj = xn;
            start_row = rowSize - 1;
        }

        // multiply the overlap of each X row and the reversed, shifted 'row'
        // put output in a temporary col
        col_temp.SetElements((T1) 0);
        index_row = start_row;

        for (int jj = start_xj; jj < stop_xj; ++jj)
        {
            T1* ct = col_temp.m_real;
            const T1* xp = &(X(0, jj));

            value = row(index_row);

            for (int i = 0; i < xm; ++i)
                (*ct++) += value * (*xp++);

            --index_row;
        }

        // convolve 'col' with the temporary col
        for (int i = row1; i <= row2; ++i)
        {
            if (i < partial_col)
            {
                start_xi = 0;
                stop_xi = i + 1;
                start_col = i;
            }
            else if (i < cm - partial_col)
            {
                if (xm < colSize)
                {
                    start_xi = 0;
                    stop_xi = xm;
                    start_col = i;
                }
                else
                {
                    start_xi = i - colSize + 1;
                    stop_xi = i + 1;
                    start_col = colSize - 1;
                }
            }
            else
            {
                start_xi = i - colSize + 1;
                stop_xi = xm;
                start_col = colSize - 1;
            }

            // multiply the overlap of col_temp and the reversed, shifted 'col'
            const T1* cp = &(col(start_col));
            T1* ct = &(col_temp(start_xi));
            value = (T1)0;

            for (int ii = start_xi; ii < stop_xi; ++ii)
            {
                value += (*cp--) * (*ct++);
            }

            (*this)(i - row1, j - col1) = value;
        }
    }

    return status;
}

// ****************************************************
//            Magnitude / Phase Operations
// ****************************************************

//! Compute the magnitudes of a matrix of values
template<>
inline hwMathStatus hwTMatrix<double>::Abs(const hwTMatrix<double>& A)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;

    status = Dimension(A.m_nRows, A.m_nCols, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    int size = Size();

    if (A.IsReal())
    {
        for (int i = 0; i < size; ++i)
            (*this)(i) = fabs(A(i));
    }
    else
    {
        for (int i = 0; i < size; ++i)
            (*this)(i) = A.z(i).Mag();
    }

    return status;
}

//! Compute the squared magnitudes of a matrix of values
template<>
inline hwMathStatus hwTMatrix<double>::AbsSq(const hwTMatrix<double>& A)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;

    status = Dimension(A.m_nRows, A.m_nCols, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    int size = Size();

    if (A.IsReal())
    {
        for (int i = 0; i < size; ++i)
            (*this)(i) = A(i) * A(i);
    }
    else
    {
        for (int i = 0; i < size; ++i)
            (*this)(i) = A.z(i).MagSq();
    }

    return status;
}

//! Compute the phases of a matrix of values
template<>
inline hwMathStatus hwTMatrix<double>::Phase(const hwTMatrix<double>& A)
{
    if (this == &A)
        return hwMathStatus(HW_MATH_ERR_NOTIMPLEMENT);

    hwMathStatus status;

    status = Dimension(A.m_nRows, A.m_nCols, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (A.IsReal())
    {
        SetElements(0.0);
    }
    else
    {
        int size = Size();

        for (int i = 0; i < size; ++i)
            (*this)(i) = A.z(i).Arg();
    }

    return status;
}

//! Compute the hypot on corresponding elements of two matrices
template<>
inline hwMathStatus hwTMatrix<double>::Hypot(const hwTMatrix<double>& A, const hwTMatrix<double>& B)
{
	hwMathStatus status;

	if (A.IsEmpty() && B.IsEmpty())
	{
		if ( (A.m_nRows == B.m_nRows) && (A.m_nCols == B.m_nCols) )
		{
			// hypot(zeros(1,0), zeros(1, 0))
			// hypot([], [])
			status = Dimension(A.m_nRows, A.m_nCols, REAL);
			if (!status.IsOk())
			{
				return status;
			}
		}
		else
		{
			// error: hypot(zeros(1,0), zeros(0, 3))
			return status(HW_MATH_ERR_ARRAYDIM, 1, 2);
		}
	}
	else if ( A.IsEmpty() || B.IsEmpty() )
	{
		if ( A.IsEmpty() && (B.Size() == 1) )
		{
			// hypot(ones(1, 0), 1)
			// hypot(ones(1, 0), i)
			status = Dimension(A.m_nRows, A.m_nCols, REAL);
			if (!status.IsOk())
			{
				return status;
			}
		}
		else if ( (A.Size() == 1) && B.IsEmpty() )
		{
			// hypot(1, ones(1, 0))
			// hypot(1+i, ones(1, 0))
			status = Dimension(B.m_nRows, B.m_nCols, REAL);
			if (!status.IsOk())
			{
				return status;
			}
		}
		else
		{
			// error: hypot(zeros(1, 0), [2 4])
			// error: hypot(zeros(1, 0), [2i 4])
			return status(HW_MATH_ERR_ARRAYDIM, 1, 2);
		}
	}
	else if ( (A.m_nRows == B.m_nRows && A.m_nCols == B.m_nCols) ||
		(A.m_nRows == 1 && A.m_nCols == 1) ||
		(B.m_nRows == 1 && B.m_nCols == 1) )
	{
		int m = _max(A.m_nRows, B.m_nRows);
		int n = _max(A.m_nCols, B.m_nCols);

		status = Dimension(m, n, REAL);

		if (!status.IsOk())
		{
			return status;
		}

		if (A.IsReal() && B.IsReal())
		{
			if (A.m_nRows == B.m_nRows && A.m_nCols == B.m_nCols)
			{
				for (int i = 0; i < A.m_nRows * A.m_nCols; i++)
				{
					(*this)(i) = _hypot(A(i), B(i));
				}
			}
			else if (A.m_nRows == 1 && A.m_nCols == 1)
			{
				for (int i = 0; i < B.m_nRows * B.m_nCols; i++)
				{
					(*this)(i) = _hypot(A(0), B(i));
				}
			}
			else // (B.m_nRows == 1 && B.m_nCols == 1)
			{
				for (int i = 0; i < A.m_nRows * A.m_nCols; i++)
				{
					(*this)(i) = _hypot(A(i), B(0));
				}
			}
		}
		else 
		{
			if (A.IsReal())
			{
				if (A.m_nRows == B.m_nRows && A.m_nCols == B.m_nCols)
				{
					for (int i = 0; i < A.m_nRows * A.m_nCols; i++)
					{
						(*this)(i) = sqrt(A(i)*A(i) + B.z(i).MagSq());
					}
				}
				else if (A.m_nRows == 1 || A.m_nCols == 1)
				{	
					double value = A(0) * A(0);
					for (int i = 0; i < B.m_nRows * B.m_nCols; i++)
					{
						(*this)(i) = sqrt(value + B.z(i).MagSq());
					}
				}
				else // (B.m_nRows == 1 && B.m_nCols == 1)
				{
					double value = B.z(0).MagSq();
					for (int i = 0; i < A.m_nRows * A.m_nCols; i++)
					{
						(*this)(i) = sqrt(A(i)*A(i) + value);
					}
				}
			}

			else if (B.IsReal())
			{
				if (A.m_nRows == B.m_nRows && A.m_nCols == B.m_nCols)
				{
					for (int i = 0; i < A.m_nRows * A.m_nCols; i++)
					{
						(*this)(i) = sqrt(A.z(i).MagSq() + B(i)*B(i));
					}
				}
				else if (A.m_nRows == 1 && A.m_nCols == 1)
				{
					double value = A.z(0).MagSq();

					for (int i = 0; i < B.m_nRows * B.m_nCols; i++)
					{
						(*this)(i) = sqrt(value + B(i)*B(i));
					}
				}
				else // (B.m_nRows == 1 && B.m_nCols == 1)
				{
					double value = B(0) * B(0);

					for (int i = 0; i < A.m_nRows * A.m_nCols; i++)
					{
						(*this)(i) = sqrt(A.z(i).MagSq() + value);
					}
				}
			}
			else // (A.m_complex && B.m_complex)
			{
				if (A.m_nRows == B.m_nRows && A.m_nCols == B.m_nCols)
				{
					for (int i = 0; i < A.m_nRows * A.m_nCols; i++)
					{
						(*this)(i) = sqrt(A.z(i).MagSq() + B.z(i).MagSq());
					}
				}
				else if (A.m_nRows == 1 && A.m_nCols == 1)
				{
					double value = A.z(0).MagSq();

					for (int i = 0; i < B.m_nRows * B.m_nCols; i++)
					{
						(*this)(i) = sqrt(value + B.z(i).MagSq());
					}
				}
				else // (B.m_nRows == 1 && B.m_nCols == 1)
				{
					double value = B.z(0).MagSq();

					for (int i = 0; i < A.m_nRows * A.m_nCols; i++)
					{
						(*this)(i) = sqrt(A.z(i).MagSq() + value);
					}
				}
			}
		}
	}
	else
	{
		return status(HW_MATH_ERR_ARRAYDIM, 1, 2);
	}
	return status;
}

//! Unwrap a matrix of phases
template<>
inline hwMathStatus hwTMatrix<double>::UnwrapMat(const hwTMatrix<double>& phase,
                                                 int dim, double tol)
{
    hwMathStatus status;

    if (this == &phase)
        return status(HW_MATH_ERR_NOTIMPLEMENT);

    if (phase.m_complex)
        return status(HW_MATH_ERR_COMPLEX, 1);

    status = Dimension(phase.m_nRows, phase.m_nCols, REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (IsEmpty())
        return status;

    if (tol <= 0.0)
        return status(HW_MATH_ERR_NONPOSITIVE, 3);

    if (tol <= PI)
        tol = PI;

    int i, j;
    int num_jumps;
    int m = phase.m_nRows;
    int n = phase.m_nCols;
    double delta;
    double phase1;
    double phase2;
    double jump = 2.0 * PI;

    if (dim == 1)       // by column
    {
        for (j = 0; j < n; ++j)
        {
            phase1 = phase(0, j);
            (*this)(0, j) = phase1;

            for (i = 1; i < m; ++i)
            {
                phase2 = phase(i, j);
                delta = phase2 - phase1;

                if (delta >= tol)
                {
                    num_jumps = static_cast<int>(ceil((delta - PI) / jump));
                    phase2 -= num_jumps * jump;
                }
                else if (-delta >= tol)
                {
                    num_jumps = static_cast<int>(ceil(-(delta + PI) / jump));
                    phase2 += num_jumps * jump;
                }

                (*this)(i, j) = phase2;
                phase1 = phase2;
            }
        }
    }
    else if (dim == 2)  // by row
    {
        for (i = 0; i < m; ++i)
        {
            phase1 = phase(i, 0);
            (*this)(i, 0) = phase1;

            for (j = 1; j < n; ++j)
            {
                phase2 = phase(i, j);
                delta = phase2 - phase1;

                if (delta >= tol)
                {
                    num_jumps = static_cast<int>(ceil((delta - PI) / jump));
                    phase2 -= num_jumps * jump;
                }
                else if (-delta >= tol)
                {
                    num_jumps = static_cast<int>(ceil(-(delta + PI) / jump));
                    phase2 += num_jumps * jump;
                }

                (*this)(i, j) = phase2;
                phase1 = phase2;
            }
        }
    }
    else
        return status(HW_MATH_ERR_ARRAYDIM, 2);
    
    return status;
}

//! Unwrap a vector of phases
template<>
inline hwMathStatus hwTMatrix<double>::UnwrapVec(const hwTMatrix<double>& phase, double tol)
{
    int m = phase.m_nRows;
    int n = phase.m_nCols;
    hwMathStatus status;

    if (m == 1)
        status = UnwrapMat(phase, 2, tol);
    else if (n == 1)
        status = UnwrapMat(phase, 1, tol);
    else
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 3)
            status.SetArg1(2);
    }

    return hwMathStatus();
}

// ****************************************************
//              Special Matrix Support
// ****************************************************

//! Generate Hankel Matrix
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Hankel(const hwTMatrix<T1, T2>& C, const hwTMatrix<T1, T2>* R)
{
    if (!C.IsVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

    int m = C.Size();
    hwTMatrix<T1, T2> MasterCol(m, 1, C.Type());
    hwMathStatus status;
    hwMathStatus statusW;

    if (R)
    {
        if (!R->IsVector())
            return status(HW_MATH_ERR_VECTOR, 2);

        int n = R->Size();

        // check R(0) == C(end)
        if (C.IsReal())
        {
            if (R->IsReal())
            {
                if ((*R)(0) != C(m - 1))
                    statusW(HW_MATH_WARN_HANKEL, 1, 2);
            }
            else
            {
                if (R->z(0) != C(m - 1))
                    statusW(HW_MATH_WARN_HANKEL, 1, 2);
            }
        }
        else
        {
            if (R->IsReal())
            {
                if (C.z(m - 1) != (*R)(0))
                    statusW(HW_MATH_WARN_HANKEL, 1, 2);
            }
            else
            {
                if (C.z(m - 1) != R->z(0))
                    statusW(HW_MATH_WARN_HANKEL, 1, 2);
            }
        }

        hwTMatrix<T1, T2> R1;
        R1.DeleteElements(*R, 0, 1);
        status = MasterCol.InsertElements(C, m, R1);

        if (!status.IsOk())
            return status;

        status = Dimension(m, n, MasterCol.Type());
    }
    else
    {
        status = MasterCol.InsertElements(C, m, m - 1);

        if (!status.IsOk())
            return status;

        status = Dimension(m, m, MasterCol.Type());
    }

    if (!status.IsOk())
        return status;

    if (IsReal())
    {
        double* start = MasterCol.GetRealData();

        for (int j = 0; j < N(); ++j)
        {
            hwTMatrix<T1, T2> Col(m, 1, (void*) start, hwTMatrix<T1, T2>::REAL);
            WriteColumn(j, Col);
            start++;
        }
    }
    else
    {
        hwTComplex<T1>* start = MasterCol.GetComplexData();

        for (int j = 0; j < N(); ++j)
        {
            hwTMatrix<T1, T2> Col(m, 1, (void*) start, hwTMatrix<T1, T2>::COMPLEX);
            WriteColumn(j, Col);
            start++;
        }
    }

    return statusW;
}

//! Generate Toeplitz Matrix
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Toeplitz(const hwTMatrix<T1, T2>& C, const hwTMatrix<T1, T2>* R)
{
    if (!C.IsVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

    int m = C.Size();
    hwTMatrix<T1, T2> MasterCol;
    const hwTMatrix<T1, T2>* RC;
    int n;
    hwMathStatus status;
    hwMathStatus statusW;

    if (R)
    {
        if (!R->IsVector())
            return status(HW_MATH_ERR_VECTOR, 2);

        n = R->Size();

        // check R(0) == C(0)
        if (C.IsReal())
        {
            if (R->IsReal())
            {
                if ((*R)(0) != C(0))
                    statusW(HW_MATH_WARN_TOEPLITZ, 1, 2);
            }
            else
            {
                if (R->z(0) != C(0))
                    statusW(HW_MATH_WARN_TOEPLITZ, 1, 2);
            }
        }
        else
        {
            if (R->IsReal())
            {
                if (C.z(0) != (*R)(0))
                    statusW(HW_MATH_WARN_TOEPLITZ, 1, 2);
            }
            else
            {
                if (C.z(0) != R->z(0))
                    statusW(HW_MATH_WARN_TOEPLITZ, 1, 2);
            }
        }

        RC = R;
    }
    else
    {
        n = m;
        RC = &C;
    }
    
    hwTMatrix<T1, T2> RV(n - 1, 1, RC->Type());

    if (RC->IsReal())
    {
        for (int i = 0; i < n - 1; ++i)
            RV(i) = (*RC)(n - 1 - i);
    }
    else
    {
        for (int i = 0; i < n - 1; ++i)
            RV.z(i) = RC->z(n - 1 - i);
    }

    status = MasterCol.InsertElements(RV, n - 1, C);

    if (!status.IsOk())
        return status;

    status = Dimension(m, n, MasterCol.Type());

    if (!status.IsOk())
        return status;

    if (IsReal())
    {
        double* start = MasterCol.GetRealData() + (n - 1);

        for (int j = 0; j < n; ++j)
        {
            hwTMatrix<T1, T2> Col(m, 1, (void*)start, hwTMatrix<T1, T2>::REAL);
            WriteColumn(j, Col);
            start--;
        }
    }
    else
    {
        hwTComplex<T1>* start = MasterCol.GetComplexData() + (n - 1);

        for (int j = 0; j < n; ++j)
        {
            hwTMatrix<T1, T2> Col(m, 1, (void*)start, hwTMatrix<T1, T2>::COMPLEX);
            WriteColumn(j, Col);
            start--;
        }
    }

    return statusW;
}

// ****************************************************
//                 Private Utilities
// ****************************************************

//! Compute memory capacity for the matrix data
template<>
inline void hwTMatrix<double>::SetCapacity(DataType dataType)
{
    if (m_nCols < 0 || m_nRows < 0)
    {
        m_capacity = -1;    // force allocaction failure
        return;
    }

    if (!m_nCols || !m_nRows)
    {
        m_capacity = 0;
        return;
    }

    constexpr size_t maxSize = std::numeric_limits<int>::max();

    if ((size_t) m_nCols > maxSize / (size_t) m_nRows)
    {
        m_capacity = -1;    // force allocaction failure
        return;
    }

    int size = Size();

    if (size > m_capacity)
    {
        if (m_capacity > maxSize - m_capacity)
        {
            m_capacity = static_cast<int>(maxSize);
        }
        else
        {
            m_capacity = _max(2 * m_capacity, size);

            if (dataType == REAL)   // extend for 16-byte alignment
                m_capacity += 2;
            else
                m_capacity++;
        }
    }
}

//! Compute memory capacity for the matrix data
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::SetCapacity(DataType dataType)
{
    if (m_nCols < 0 || m_nRows < 0)
    {
        m_capacity = -1;    // force allocaction failure
        return;
    }

    if (!m_nCols || !m_nRows)
    {
        m_capacity = 0;
        return;
    }

    constexpr size_t maxSize = std::numeric_limits<int>::max();
    
    if ((size_t) m_nCols > maxSize / (size_t) m_nRows)
    {
        m_capacity = -1;    // force allocaction failure
        return;
    }

    int size = Size();

    if (size > m_capacity)
    {
        if (m_capacity > maxSize - m_capacity)
            m_capacity = static_cast<int>(maxSize);
        else
            m_capacity = _max(2 * m_capacity, size);
    }
}

//! Delete the matrix data
template<>
inline void hwTMatrix<double>::Deallocate()
{
    if (m_bits.ownData)
    {
        if (m_real_memory)
            delete [] m_real_memory;
        else if (m_real)
            delete [] m_real;	    // ownership of external data was assumed

        if (m_complex_memory)
            delete [] m_complex_memory;
        else if (m_complex)
            delete [] m_complex;    // ownership of external data was assumed
    }
}

//! Delete the matrix data
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::Deallocate()
{
    if (m_bits.ownData)
    {
        if (m_real)
            delete [] m_real;

        if (m_complex)
            delete [] m_complex;
    }
}

//! Allocate memory for the matrix data
template<>
inline void hwTMatrix<double>::Allocate(DataType dataType)
{
    SetCapacity(dataType);

    if (m_capacity == 0)
    {
        return;
    }

    // allocate contiguous block
    if (dataType == REAL)
    {
        if (m_real)
        {
            MakeEmpty();
            m_capacity = -1;    // force allocaction failure
        }

        try
        {
            m_real_memory = new char[m_capacity * sizeof(double)];
            m_real = reinterpret_cast<double*>
                     ((reinterpret_cast<std::ptrdiff_t> (m_real_memory) + 15) & ~ 0xF);   // 16-byte alignment
        }
        catch (std::bad_alloc&)
        {
            Deallocate();
            m_real_memory = nullptr;
            throw;
        }
    }
    else
    {
        if (m_complex)
        {
            MakeEmpty();
            m_capacity = -1;    // force allocaction failure
        }

        try
        {
            m_complex_memory = new char[m_capacity * sizeof(hwTComplex<double>)];
            m_complex = reinterpret_cast<hwTComplex<double>*>
                        ((reinterpret_cast<std::ptrdiff_t> (m_complex_memory) + 15) & ~ 0xF);   // 16-byte alignment
        }
        catch (std::bad_alloc&)
        {
            Deallocate();
            m_complex_memory = nullptr;
            throw;
        }
    }
}

//! Allocate memory for the matrix data
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::Allocate(DataType dataType)
{
    SetCapacity(dataType);

    if (m_capacity == 0)
    {
        return;
    }

    // allocate contiguous block
    if (dataType == REAL)
    {
        if (m_real)
        {
            MakeEmpty();
            m_capacity = -1;    // force allocaction failure
        }

        try
        {
            m_real = new T1[m_capacity];
        }
        catch (std::bad_alloc&)
        {
            Deallocate();
            m_real = nullptr;
            throw;
        }
    }
    else
    {
        if (m_complex)
        {
            MakeEmpty();
            m_capacity = -1;    // force allocaction failure
        }

        try
        {
            m_complex = new T2[m_capacity];
        }
        catch (std::bad_alloc&)
        {
            Deallocate();
            m_complex = nullptr;
            throw;
        }
    }
}

//! Release real memory
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::FreeMemory(char*& pMemory, T1*& pReal)
{
    if (pMemory)     // aligned memory
    {
        delete [] pMemory;
        pMemory = nullptr;
        pReal = nullptr;
    }
    else if (pReal)
    {
        delete [] (T1*) pReal;
        pReal = nullptr;
    }
}

//! Release complex memory
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::FreeMemory(char*& pMemory, T2*& pComplex)
{
    if (pMemory)     // aligned memory
    {
        delete [] pMemory;
        pMemory = nullptr;
        pComplex = nullptr;
    }
    else if (pComplex)
    {
        delete [] (T2*) pComplex;
        pComplex = nullptr;
    }
}

//! Set matrix to empty condition
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::MakeEmpty()
{
    Deallocate();
    m_capacity = 0;
    m_real_memory = nullptr;
    m_complex_memory = nullptr;
    m_real = nullptr;
    m_complex = nullptr;
    m_nCols = 0;
    m_nRows = 0;
    m_bits.ownData = 1;
    m_bits.realData = 1;
}

//! Copy function
template<typename T1, typename T2>
hwMathStatus hwTMatrix<T1, T2>::Copy(const hwTMatrix<T1, T2>& source)
{
    bool sizeChange;

    if (m_nCols != source.m_nCols || m_nRows != source.m_nRows)
    {
        if (Size() != source.Size())
            sizeChange = true;
        else
            sizeChange = false;

        m_nCols = source.m_nCols;
        m_nRows = source.m_nRows;
    }
    else
        sizeChange = false;

    m_capacity = source.m_capacity;
    m_bits.realData = source.m_bits.realData;

    if (m_bits.ownData)
    {
        if (source.m_real)
        {
            if (!m_real)
            {
                try
                {
                    Allocate(REAL);
                }
                catch (std::bad_alloc&)
                {
                    return hwMathStatus(HW_MATH_ERR_ALLOCFAILED, 0);
                }
            }
            else if (sizeChange)
            {
                FreeMemory(m_real_memory, m_real);

                try
                {
                    Allocate(REAL);
                }
                catch (std::bad_alloc&)
                {
                    return hwMathStatus(HW_MATH_ERR_ALLOCFAILED, 0);
                }
            }
        }
        else if (m_real)
        {
            FreeMemory(m_real_memory, m_real);
        }

        if (source.m_complex)
        {
            if (!m_complex)
            {
                try
                {
                    Allocate(COMPLEX);
                }
                catch (std::bad_alloc&)
                {
                    return hwMathStatus(HW_MATH_ERR_ALLOCFAILED, 0);
                }
            }
            else if (sizeChange)
            {
                FreeMemory(m_complex_memory, m_complex);

                try
                {
                    Allocate(COMPLEX);
                }
                catch (std::bad_alloc&)
                {
                    return hwMathStatus(HW_MATH_ERR_ALLOCFAILED, 0);
                }
            }
        }
        else if (m_complex)
        {
            FreeMemory(m_complex_memory, m_complex);
        }
    }
    else if (sizeChange)
    {
        // copying data is not allowed when the target matrix does not
        // own its data and there is a change of dimension
        MakeEmpty();
        return hwMathStatus(HW_MATH_ERR_NOTALLOWED);
    }

    int newSize = Size();

    if (newSize > 0)
    {
        if (m_real)
            CopyData(m_real, newSize, source.m_real, newSize);
        else if (m_complex)
            CopyData(m_complex, newSize, source.m_complex, newSize);
    }

    return hwMathStatus();
}

//! Copy a real submatrix from another matrix to *this
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::CopyBlock(const T1* real, int m, int n, int row1, int row2,
                                  int col1, int col2, int ii, int jj)
{
    // This private function performs no bounds checks, but will not write beyond the *this bounds
    // T1 is an m by n matrix that contains the submatrix to be copied
    // row1 and row2 are the first and last rows the submatrix
    // col1 and col2 are the first and last columns the submatrix
    // ii and jj are the indices of the first element to which the submatrix will be copied
    int m_max = _min(row1+m_nRows-1-ii, row2);
    int n_max = _min(col1+m_nCols-1-jj, col2);
    int offset_t;   // target offset
    int offset_s;   // source offset
    int numElems;
    int numRows = m_max - row1 + 1;

    if (numRows == m && numRows == m_nRows)  // contiguous memory for both source and target
    {
        numElems = numRows * (n_max-col1+1);

        if (numElems)
        {
            offset_t = jj * m_nRows;
            offset_s = col1 * m;
            CopyData(m_real + offset_t, numElems, real + offset_s, numElems);
        }
    }
    else
    {
        numElems = numRows;

        if (numElems)
        {
            offset_t = jj * m_nRows + ii;
            offset_s = col1 * m + row1;

            for (int i = col1; i <= n_max; ++i)
            {
                CopyData(m_real + offset_t, numElems, real + offset_s, numElems);
                offset_t += m_nRows;
                offset_s += m;
            }
        }
    }
}

//! Copy a complex submatrix from another matrix to *this
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::CopyBlock(const T2* cmplx, int m, int n, int row1, int row2,
                                  int col1, int col2, int ii, int jj)
{
    // This private function performs no bounds checks, but will not write beyond the *this bounds
    // T1 is an m by n matrix that contains the submatrix to be copied
    // row1 and row2 are the first and last rows the submatrix
    // col1 and col2 are the first and last columns the submatrix
    // ii and jj are the indices of the first element to which the submatrix will be copied
    int m_max = _min(row1+m_nRows-1-ii, row2);
    int n_max = _min(col1+m_nCols-1-jj, col2);
    int offset_t;   // target offset
    int offset_s;   // source offset
    int numElems;
    int numRows = m_max - row1 + 1;

    if (numRows == m && numRows == m_nRows)  // contiguous memory for both source and target
    {
        numElems = numRows * (n_max-col1+1);

        if (numElems)
        {
            offset_t = jj * m_nRows;
            offset_s = col1 * m;
            CopyData(m_complex + offset_t, numElems, cmplx + offset_s, numElems);
        }
    }
    else
    {
        numElems = numRows;

        if (numElems)
        {
            offset_t = jj * m_nRows + ii;
            offset_s = col1 * m + row1;

            T2* start_t = m_complex + offset_t;
            T2* start_s = (T2*) cmplx + offset_s;

            for (int i = col1; i <= n_max; ++i)
            {
                CopyData(start_t, numElems, start_s, numElems);
                start_t += m_nRows;
                start_s += m;
            }
        }
    }
}

//! Copy data
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::CopyData(void* dest, int arraySize, const void* src, int count)
{
    // memcpy cannot be used here
    if (m_real)
    {
        T1* destTemp = reinterpret_cast<T1*>(dest);
        const T1* srcTemp = reinterpret_cast<T1*>(const_cast<void *>(src));
        
        for (int i = 0; i < count; ++i)
            destTemp[i] = srcTemp[i];
    }
    else if (m_complex)
    {
        T2* destTemp = reinterpret_cast<T2*>(dest);
        const T2* srcTemp = reinterpret_cast<T2*>(const_cast<void *>(src));
        
        for (int i = 0; i < count; ++i)
            destTemp[i] = srcTemp[i];
    }
}

//! Set a submatrix of *this to zeros
template<typename T1, typename T2>
void hwTMatrix<T1, T2>::ZeroBlock(int row1, int row2, int col1, int col2)
{
    int numBytes;
    int numRows = row2 - row1 + 1;

    if (m_real)
    {
        T1* start = m_real + (static_cast<long int>(col1) * m_nRows + row1);

        if (numRows == m_nRows)  // contiguous memory
        {
            numBytes = numRows * (col2-col1+1) * sizeof(T1);

            if (numBytes > 0)
                memset(start, 0, numBytes);
        }
        else
        {
            numBytes = numRows * sizeof(T1);

            if (numBytes > 0)
            {
                for (int i = col1; i <= col2; ++i)
                {
                    memset(start, 0, numBytes);
                    start += m_nRows;
                }
            }
        }
    }

    if (m_complex)
    {
        T2* start = m_complex + (static_cast<long int>(col1) * m_nRows + row1);

        if (numRows == m_nRows)  // contiguous memory
        {
            numBytes = numRows * (col2-col1+1) * sizeof(T2);

            if (numBytes > 0)
                memset(start, 0, numBytes);
        }
        else
        {
            numBytes = numRows * sizeof(T2);

            if (numBytes > 0)
            {
                for (int i = col1; i <= col2; ++i)
                {
                    memset(start, 0, numBytes);
                    start += m_nRows;
                }
            }
        }
    }
}
