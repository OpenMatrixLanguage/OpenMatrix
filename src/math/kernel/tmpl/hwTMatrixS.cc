/**
* @file hwTMatrixS.cc
* @date May 2019
* Copyright (C) 2014-2019 Altair Engineering, Inc.
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
//  hwTMatrixS template function implementation file
//
//:---------------------------------------------------------------------------

#include <algorithm>
#include <hwMathException.h>
#include <hwSliceArg.h>
#include <tmpl/hwTComplex.h>
// #include "mkl_dss.h"
// #include "mkl_spblas.h"

// ****************************************************
//                  Error handling
// ****************************************************
//  Clients of hwTMatrixS should use try blocks
//  Here is an example of an apprporiate catch block
//
//  catch (hwMathException& except)
//  {
//      const char* msg = except.what();
//  }
// ****************************************************

// ****************************************************
//       Construction / Destruction / Assignment
// ****************************************************

//! Construct empty matrix
template<typename T1, typename T2>
hwTMatrixS<T1, T2>::hwTMatrixS()
    : m_nRows(0), m_nCols(0), m_refCount(1)
{

}

//! Construct sparse from dense
template<typename T1, typename T2>
hwTMatrixS<T1, T2>::hwTMatrixS(const hwTMatrix<T1, T2>& A)
    : m_nRows(A.M()), m_nCols(A.N()), m_refCount(1)
{
    int nnz = 0;
    m_pointerB.resize(m_nCols);
    m_pointerE.resize(m_nCols);

    if (m_nCols)
        m_pointerB[0] = 0;

    // populate m_rows
    if (A.IsReal())
    {
        m_values.Dimension(1, hwTMatrix<T1, T2>::REAL);

        for (int j = 0; j < m_nCols; ++j)
        {
            for (int i = 0; i < m_nRows; ++i)
            {
                if (A(i, j) != 0.0)
                {
                    m_values.Resize(nnz+1);
                    m_values(nnz) = A(i, j);
                    m_rows.push_back(i);
                    ++nnz;
                }
            }

            m_pointerE[j] = nnz;

            if (j != m_nCols - 1)
                m_pointerB[j + 1] = nnz;
        }
    }
    else
    {
        m_values.Dimension(1, hwTMatrix<T1, T2>::COMPLEX);

        for (int j = 0; j < m_nCols; ++j)
        {
            for (int i = 0; i < m_nRows; ++i)
            {
                if (A.z(i, j) != 0.0)
                {
                    m_values.Resize(nnz + 1);
                    m_values.z(nnz) = A.z(i, j);
                    m_rows.push_back(i);
                    ++nnz;
                }
            }

            m_pointerE[j] = nnz;

            if (j != m_nCols - 1)
                m_pointerB[j + 1] = nnz;
        }
    }

    if (!nnz)
    {
        m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);
    }
}

// Comparison function to stable_sort pairs using only first elements
template<typename T>
bool compare(const std::pair<long long int, T>& a, const std::pair<long long int, T>& b)
{
    return a.first < b.first;
}

//! Construct sparse from vector representation
template<typename T1, typename T2>
hwTMatrixS<T1, T2>::hwTMatrixS(const std::vector<int>&  ivec,
                               const std::vector<int>&  jvec,
                               const hwTMatrix<T1, T2>& V,
                               int                      m,
                               int                      n,
                               const char*              option)
    : m_nRows(m), m_nCols(n), m_refCount(1)
{
    int nnz = V.Size();

    if (ivec.size() != nnz)
    {
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 3);
    }

    if (jvec.size() != nnz)
    {
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    if (m == -1 && n == -1)
    {
        m_nRows = 0;
        m_nCols = 0;

        for (int i = 0; i < nnz; ++i)
        {
            int ii = ivec[i];
            int jj = jvec[i];

            if (ii < 0)
                throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);

            if (jj < 0)
                throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 2);

            if (ii > m_nRows - 1)
                m_nRows = ii + 1;

            if (jj > m_nCols - 1)
                m_nCols = jj + 1;
        }
    }
    else if (m < 0 || n < 0)
    {
        throw hwMathException(HW_MATH_ERR_ARRAYDIM);
    }

    // check for zero values
    int nz = 0;

    if (V.IsReal())
    {
        for (int i = 0; i < nnz; ++i)
        {
            if (V(i) == 0.0)
                ++nz;
        }
    }
    else
    {
        for (int i = 0; i < nnz; ++i)
        {
            if (V.z(i) == 0.0)
                ++nz;
        }
    }

    nnz -= nz;

    // populate members
    hwMathStatus status = m_values.Dimension(nnz, V.Type());
    m_rows.resize(nnz);
    m_pointerB.resize(m_nCols);
    m_pointerE.resize(m_nCols);

    if (m_nRows == 0 || m_nCols == 0)
        return;

    m_pointerB[0] = 0;
    m_pointerE[0] = 0;
    int col;
    long long int index_prev = -1;
    int col_prev   = -1;
    int nDuplicates = 0;

    nz = 0;

    if (V.IsReal())
    {
        std::vector<std::pair<long long int, T1>> pairs(nnz);

        for (int i = 0; i < nnz + nz; ++i)
        {
            if (V(i) != 0.0)
            {
                if (ivec[i] >= m_nRows)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);

                if (jvec[i] >= m_nCols)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 2);

                pairs[i-nz].first = jvec[i] * static_cast<long long int> (m_nRows) + ivec[i];
                pairs[i-nz].second = V(i);
            }
            else
            {
                ++nz;
            }
        }

        stable_sort(pairs.begin(), pairs.end(), compare<T1>);

        for (int i = 0; i < nnz; ++i)
        {
            long long int index = pairs[i].first;

            if (index == index_prev)
            {
                ++nDuplicates;
                if (option && !strcmp(option, "unique"))
                {
                    m_values(i - nDuplicates) = pairs[i].second;
                }
                else
                {
                    m_values(i - nDuplicates) += pairs[i].second;
                }
            }
            else
            {
                m_values(i - nDuplicates) = pairs[i].second;
                m_rows[i - nDuplicates] = static_cast<int> (index % m_nRows);
                index_prev = index;
                col = static_cast<int> (index / m_nRows);

                if (col != col_prev)
                {
                    if (col_prev > -1)
                    {
                        for (int k = col_prev + 1; k <= col; ++k)
                            m_pointerB[k] = m_pointerE[k] = m_pointerE[k - 1];
                    }

                    col_prev = col;
                }

                ++m_pointerE[col];
            }
        }
    }
    else   // complex V
    {
        std::vector<std::pair<long long int, T2>> pairs(nnz);

        for (int i = 0; i < nnz; ++i)
        {
            if (V.z(i) != 0.0)
            {
                pairs[i-nz].first = jvec[i] * static_cast<long long int> (m_nRows) + ivec[i];
                pairs[i-nz].second = V.z(i);
            }
            else
            {
                ++nz;
            }
        }

        stable_sort(pairs.begin(), pairs.end(), compare<T2>);

        for (int i = 0; i < nnz; ++i)
        {
            long long int index = pairs[i].first;

            if (index == index_prev)
            {
                ++nDuplicates;
                if (option && !strcmp(option, "unique"))
                {
                    m_values.z(i - nDuplicates) = pairs[i].second;
                }
                else
                {
                    m_values.z(i - nDuplicates) += pairs[i].second;
                }
            }
            else
            {
                m_values.z(i - nDuplicates) = pairs[i].second;
                m_rows[i - nDuplicates] = static_cast<int> (index % m_nRows);
                index_prev = index;
                col = static_cast<int> (index / m_nRows);

                if (col != col_prev)
                {
                    if (col_prev > -1)
                    {
                        for (int k = col_prev + 1; k <= col; ++k)
                            m_pointerB[k] = m_pointerE[k] = m_pointerE[k - 1];
                    }

                    col_prev = col;
                }

                ++m_pointerE[col];
            }
        }
    }

    if (nDuplicates)
    {
        nnz   -= nDuplicates;
        status = m_values.Resize(nnz);
        m_rows.resize(nnz);

        if (V.IsReal())
        {
            for (int ii = 0; ii < nnz; ++ii)
            {
                int row;
                int col;
                T1 value;
    
                NZinfo(ii, row, col, value);
    
                if (value == 0.0)
                {
                    ZeroElement(row, col);
                    --nnz;
                    --ii;
                }
            }
        }
        else
        {
            for (int ii = 0; ii < nnz; ++ii)
            {
                int row;
                int col;
                T2 value;

                NZinfo(ii, row, col, value);

                if (value == 0.0)
                {
                    ZeroElement(row, col);
                    --nnz;
                    --ii;
                }
            }
        }
    }
}

//! Copy constructor
template<typename T1, typename T2>
hwTMatrixS<T1, T2>::hwTMatrixS(const hwTMatrixS<T1, T2>& source)
    : m_nRows(source.m_nRows), m_nCols(source.m_nCols), m_values(source.m_values),
      m_rows(source.m_rows), m_pointerB(source.m_pointerB),
      m_pointerE(source.m_pointerE), m_refCount(1)
{
}

//! Destructor
template<typename T1, typename T2>
hwTMatrixS<T1, T2>::~hwTMatrixS()
{
    if (m_refCount > 1)
        return;
}

//! Implement the = operator
template<typename T1, typename T2>
hwTMatrixS<T1, T2>& hwTMatrixS<T1, T2>::operator=(const hwTMatrixS<T1, T2>& rhs)
{
    if (this == &rhs)
        return *this;

    m_refCount = 1;
    m_nRows    = rhs.m_nRows;
    m_nCols    = rhs.m_nCols;
    m_values   = rhs.m_values;
    m_rows     = rhs.m_rows;
    m_pointerB = rhs.m_pointerB;
    m_pointerE = rhs.m_pointerE;

    return *this;
}

// ****************************************************
//                 Dense Functions
// ****************************************************
//! Create dense from sparse
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Full(hwTMatrix<T1, T2>& A) const
{
    hwMathStatus status;
    int i = 0;

    if (IsReal())
    {
        status = A.Dimension(m_nRows, m_nCols, hwTMatrix<T1, T2>::REAL);

        if (!status.IsOk())
            throw hwMathException(HW_MATH_ERR_OUTOFMEMORY);

        A.SetElements(static_cast<T1> (0));

        for (int j = 0; j < m_nCols; ++j)
        {
            for (; i < m_pointerE[j]; ++i)
            {
                A(m_rows[i], j) = m_values(i);
            }
        }
    }
    else
    {
        status = A.Dimension(m_nRows, m_nCols, hwTMatrix<T1, T2>::COMPLEX);

        if (!status.IsOk())
            throw hwMathException(HW_MATH_ERR_OUTOFMEMORY);

        A.SetElements(static_cast<T1> (0));

        for (int j = 0; j < m_nCols; ++j)
        {
            for (; i < m_pointerE[j]; ++i)
            {
                A.z(m_rows[i], j) = m_values.z(i);
            }
        }
    }
}

//! Find the real data for a storage index
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::NZinfo(int index, int& row, int& col, T1& value) const
{
    value = m_values(index);
    row   = m_rows[index];

    for (col = 0; col < m_nCols; ++col)
    {
        if (index < m_pointerE[col])
            break;
    }
}

//! Find the complex data for a storage index
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::NZinfo(int index, int& row, int& col, T2& value) const
{
    value = m_values.z(index);
    row = m_rows[index];

    for (col = 0; col < m_nCols; ++col)
    {
        if (index < m_pointerE[col])
            break;
    }
}

//! Find the non-zero data
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::NZinfo(int first, int last, std::vector<int>& row,
                                std::vector<int>& column, hwTMatrix<T1, T2>& value) const
{
    if (first < 0)
        throw hwMathException(HW_MATH_ERR_BADRANGE);

    if (last < first || last > NNZ() - 1)
        throw hwMathException(HW_MATH_ERR_BADRANGE);

    row.resize(last - first + 1);
    column.resize(last - first + 1);
    int col = 0;

    if (IsReal())
    {
        hwMathStatus status = value.Dimension(last - first + 1, 1, hwTMatrix<T1, T2>::REAL);

        for (int i = first; i <= last; ++i)
        {
            row[i - first] = m_rows[i];
            value(i - first) = m_values(i);

            for (; col < m_nCols; ++col)
            {
                if (i < m_pointerE[col])
                    break;
            }

            column[i - first] = col;
        }
    }
    else
    {
        hwMathStatus status = value.Dimension(last - first + 1, hwTMatrix<T1, T2>::COMPLEX);

        for (int i = first; i <= last; ++i)
        {
            row[i - first] = m_rows[i];
            value.z(i - first) = m_values.z(i);

            for (; col < m_nCols; ++col)
            {
                if (i < m_pointerE[col])
                    break;
            }

            column[i - first] = col;
        }
    }
}

//! Find the single indicies of the non-zero data
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::NZinfo(int first, int last, std::vector<int>& index) const
{
    if (first < 0)
        throw hwMathException(HW_MATH_ERR_BADRANGE);

    if (last < first || last > NNZ() - 1)
        throw hwMathException(HW_MATH_ERR_BADRANGE);

    index.resize(last - first + 1);

    int col = 0;

    for (int i = first; i <= last; ++i)
    {
        for (; col < m_nCols; ++col)
        {
            if (i < m_pointerE[col])
                break;
        }

        index[i - first] = m_nRows * col + m_rows[i];
    }
}

// ****************************************************
//         Access Functions for Real Elements
// ****************************************************

//! Return a reference to the real data element at the specified index
//! Client is responsible for bound checking
template<typename T1, typename T2>
T1& hwTMatrixS<T1, T2>::operator()(int index)
{
    int row = index % m_nRows;
    int col = index / m_nRows;
    int k;

    for (k = m_pointerB[col]; k < m_pointerE[col]; ++k)
    {
        if (m_rows[k] < row)
            continue;

        if (m_rows[k] > row)
            break;

        return m_values(k);
    }

    // insert a zero to be reassigned
    if (m_values.IsEmpty())
    {
        hwMathStatus status = m_values.Dimension(k + 1, hwTMatrix<T1, T2>::REAL);
        m_rows.resize(1);
        m_rows[0] = row;
    }
    else
    {
        hwTMatrix<T1, T2> copy(m_values);
        hwMathStatus status = m_values.InsertElements(copy, k);
        std::vector<int>::iterator it = m_rows.begin() + k;
        m_rows.insert(it, row);
    }

    for (int j = col; j < m_nCols - 1; ++j)
    {
        ++m_pointerB[j + 1];
        ++m_pointerE[j];
    }

    ++m_pointerE[m_nCols - 1];

    return m_values(k);
}

//! Return the real data element at the specified index
//! Client is responsible for bound checking
template<typename T1, typename T2>
T1 hwTMatrixS<T1, T2>::operator()(int index) const
{
    int row = index % m_nRows;
    int col = index / m_nRows;

    for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
    {
        if (m_rows[k] < row)
            continue;

        if (m_rows[k] > row)
            break;  // return a 0

        return m_values(k);
    }

    return static_cast<T1> (0);
}

//! Return a reference to the real data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
T1& hwTMatrixS<T1, T2>::operator()(int i, int j)
{
    int k;

    for (k = m_pointerB[j]; k < m_pointerE[j]; ++k)
    {
        if (m_rows[k] < i)
            continue;

        if (m_rows[k] > i)
            break;

        return m_values(k);
    }

    // insert a zero to be reassigned
    if (m_values.IsEmpty())
    {
        hwMathStatus status = m_values.Dimension(k + 1, hwTMatrix<T1, T2>::REAL);
        m_rows.resize(1);
        m_rows[0] = i;
    }
    else
    {
        hwTMatrix<T1, T2> copy(m_values);
        hwMathStatus status = m_values.InsertElements(copy, k);
        std::vector<int>::iterator it = m_rows.begin() + k;
        m_rows.insert(it, i);
    }
    
    for ( ; j < m_nCols - 1; ++j)
    {
        ++m_pointerB[j + 1];
        ++m_pointerE[j];
    }

    ++m_pointerE[m_nCols - 1];

    return m_values(k);
}

//! Return the real data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
T1 hwTMatrixS<T1, T2>::operator()(int i, int j) const
{
    for (int k = m_pointerB[j]; k < m_pointerE[j]; ++k)
    {
        if (m_rows[k] < i)
            continue;

        if (m_rows[k] > i)
            break;  // return a 0

        return m_values(k);
    }

    return static_cast<T1> (0);
}

//! Zero the element at the specified single index
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::ZeroElement(int index)
{
    int row = index % m_nRows;
    int col = index / m_nRows;

    ZeroElement(row, col);
}

//! Zero the element at the specified specified indices
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::ZeroElement(int i, int j)
{
    if (i < 0 || i >= m_nRows)
        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

    if (j < 0 || j >= m_nCols)
        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

    for (int k = m_pointerB[j]; k < m_pointerE[j]; ++k)
    {
        if (m_rows[k] < i)
            continue;

        if (m_rows[k] > i)
            break;

        m_values.DeleteElements(k);
        m_rows.erase(m_rows.begin() + k);

        for (int jj = j + 1; jj < m_nCols; ++jj)
            --m_pointerB[jj];

        for (int jj = j; jj < m_nCols; ++jj)
            --m_pointerE[jj];

        break;
    }
}

// ****************************************************
//         Access Functions for Complex Elements
// ****************************************************

//! Return a reference to the complex data element at the specified index
//! Client is responsible for bound checking
template<typename T1, typename T2>
T2& hwTMatrixS<T1, T2>::z(int index)
{
    int row = index % m_nRows;
    int col = index / m_nRows;
    int k;

    for (k = m_pointerB[col]; k < m_pointerE[col]; ++k)
    {
        if (m_rows[k] < row)
            continue;

        if (m_rows[k] > row)
            break;

        return m_values.z(k);
    }

    // insert a zero to be reassigned
    if (m_values.IsEmpty())
    {
        hwMathStatus status = m_values.Dimension(k + 1, hwTMatrix<T1, T2>::COMPLEX);
        m_rows.resize(1);
        m_rows[0] = row;
    }
    else
    {
        hwTMatrix<T1, T2> copy(m_values);
        hwMathStatus status = m_values.InsertElements(copy, k);
        std::vector<int>::iterator it = m_rows.begin() + k;
        m_rows.insert(it, row);
    }

    for (int j = col; j < m_nCols - 1; ++j)
    {
        ++m_pointerB[j + 1];
        ++m_pointerE[j];
    }

    ++m_pointerE[m_nCols - 1];

    return m_values.z(k);
}

//! Return the complex data element at the specified index
//! Client is responsible for bound checking
template<typename T1, typename T2>
T2 hwTMatrixS<T1, T2>::z(int index) const
{
    int row = index % m_nRows;
    int col = index / m_nRows;

    for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
    {
        if (m_rows[k] < row)
            continue;

        if (m_rows[k] > row)
            break;  // return a 0

        return m_values.z(k);
    }

    return T2(static_cast<T1> (0), static_cast<T1> (0));
}

//! Return a reference to the complex data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
T2& hwTMatrixS<T1, T2>::z(int i, int j)
{
    int k;

    for (k = m_pointerB[j]; k < m_pointerE[j]; ++k)
    {
        if (m_rows[k] < i)
            continue;

        if (m_rows[k] > i)
            break;

        return m_values.z(k);
    }

    // insert a zero to be reassigned
    if (m_values.IsEmpty())
    {
        hwMathStatus status = m_values.Dimension(k + 1, hwTMatrix<T1, T2>::COMPLEX);
        m_rows.resize(1);
        m_rows[0] = i;
    }
    else
    {
        hwTMatrix<T1, T2> copy(m_values);
        hwMathStatus status = m_values.InsertElements(copy, k);
        std::vector<int>::iterator it = m_rows.begin() + k;
        m_rows.insert(it, i);
    }

    for ( ; j < m_nCols - 1; ++j)
    {
        ++m_pointerB[j + 1];
        ++m_pointerE[j];
    }

    ++m_pointerE[m_nCols - 1];

    return m_values.z(k);
}

//! Return the complex data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
T2 hwTMatrixS<T1, T2>::z(int i, int j) const
{
    for (int k = m_pointerB[j]; k < m_pointerE[j]; ++k)
    {
        if (m_rows[k] < i)
            continue;

        if (m_rows[k] > i)
            break;  // return a 0

        return m_values.z(k);
    }

    return T2(static_cast<T1> (0), static_cast<T1> (0));
}

// ****************************************************
//                Set Matrix Dimensions
// ****************************************************

//! Change the dimensions of a matrix while maintaining the same number of elements
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Reshape(int m, int n)
{
    if (m < 0 && m != -1)
        throw hwMathException(HW_MATH_ERR_ARRAYDIM, 2);

    int size = Size();

    if (m == -1)
    {
        if (n == -1)
            throw hwMathException(HW_MATH_ERR_MATRIXRESHAPE1, 1, 2);

        m = size / n;
    }
    else if (n == -1)
    {
        n = size / m;
    }

    if (m * n != size)
        throw hwMathException(HW_MATH_ERR_MATRIXRESHAPE2, 1, 2);

    // reassign element locations
    int col = 0;
    std::vector<int> accumulator(n);

    for (int ii = 0; ii < m_rows.size(); ++ii)
    {
        // get current (row, col)
        int row = m_rows[ii];

        for (col = 0; col < m_nCols; ++col)
        {
            if (ii < m_pointerE[col])
                break;
        }

        // get current (row, col)
        int index = m_nRows * col + row;
        m_rows[ii] = index % m;
        ++accumulator[index / m];
    }

    m_pointerB.resize(n);
    m_pointerE.resize(n);
    m_pointerE[0] = accumulator[0];

    for (int ii = 1; ii < n; ++ii)
    {
        m_pointerE[ii] = accumulator[ii] + m_pointerE[ii - 1];
        m_pointerB[ii] = m_pointerE[ii - 1];
    }

    m_nRows = m;
    m_nCols = n;
}

// ****************************************************
//             Real / Complex Conversions
// ****************************************************

//! Pack real and imaginary components into a complex matrix
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::PackComplex(const hwTMatrixS<T1, T2>& real, const hwTMatrixS<T1, T2>* imag)
{
    m_nRows = real.m_nRows;
    m_nCols = real.m_nCols;
    m_rows = real.m_rows;
    m_pointerB = real.m_pointerB;
    m_pointerE = real.m_pointerE;

    if (imag)
    {
        if (m_rows != imag->m_rows)
            throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

        if (m_pointerB != imag->m_pointerB)
            throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

        if (m_pointerE != imag->m_pointerE)
            throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

        hwMathStatus status = m_values.PackComplex(real.m_values, &imag->m_values);
    }
    else
    {
        hwMathStatus status = m_values.PackComplex(real.m_values, nullptr);
    }
}

//! Unpack real and imaginary components from a complex matrix
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::UnpackComplex(hwTMatrixS<T1, T2>* real, hwTMatrixS<T1, T2>* imag) const
{
    hwMathStatus status;

    if (real)
    {
        hwTMatrix<T1, T2> realD;

        if (imag)
        {
            hwTMatrix<T1, T2> imagD;
            hwMathStatus status = m_values.UnpackComplex(&realD, &imagD);
            imag->m_nRows    = m_nRows;
            imag->m_nCols    = m_nCols;
            imag->m_values   = imagD;
            imag->m_rows     = m_rows;
            imag->m_pointerB = m_pointerB;
            imag->m_pointerE = m_pointerE;
        }
        else
        {
            m_values.UnpackComplex(&realD, nullptr);
        }

        real->m_nRows    = m_nRows;
        real->m_nCols    = m_nCols;
        real->m_values   = realD;
        real->m_rows     = m_rows;
        real->m_pointerB = m_pointerB;
        real->m_pointerE = m_pointerE;
    }
    else if (imag)
    {
        hwTMatrix<T1, T2> imagD;
        hwMathStatus status = m_values.UnpackComplex(nullptr, &imagD);
        imag->m_nRows    = m_nRows;
        imag->m_nCols    = m_nCols;
        imag->m_values   = imagD;
        imag->m_rows     = m_rows;
        imag->m_pointerB = m_pointerB;
        imag->m_pointerE = m_pointerE;
    }
}

// ****************************************************
//                 Slice Operations
// ****************************************************

//! Read a matrix slice from the calling object, as if the calling object is being
//! sliced on the right hand side of an equals sign
//! lhs_matrix = matrix(slice args)
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::SliceRHS(const std::vector<hwSliceArg>& sliceArg,
                                  hwTMatrixS<T1, T2>& lhsMatrix) const
{
    if (sliceArg.size() == 1)
    {
        if (sliceArg[0].IsScalar())
        {
            throw hwMathException(HW_MATH_ERR_NOTIMPLEMENT);
        }
        else if (sliceArg[0].IsColon())
        {
            lhsMatrix.m_nRows = Size();
            lhsMatrix.m_nCols = 1;
            lhsMatrix.m_values = m_values;
            lhsMatrix.m_rows.resize(m_rows.size());
            lhsMatrix.m_pointerB.resize(1);
            lhsMatrix.m_pointerE.resize(1);
            lhsMatrix.m_pointerB[0] = 0;
            lhsMatrix.m_pointerE[0] = static_cast<int> (m_rows.size());

            for (int ii = 0; ii < m_rows.size(); ++ii)
            {
                int col;

                for (col = 0; col < m_nCols; ++col)
                {
                    if (ii < m_pointerE[col])
                        break;
                }

                lhsMatrix.m_rows[ii] = m_nRows * col + m_rows[ii];
            }
        }
        else if (sliceArg[0].IsVector())
        {
            if (!IsVector())
                throw hwMathException(HW_MATH_ERR_NOTALLOWED);  // orientation of argument is unknown

            if (m_nCols > 1)
            {
                lhsMatrix.m_nRows = 1;
                lhsMatrix.m_nCols = static_cast<int> (sliceArg[0].Vector().size());
            }
            else
            {
                lhsMatrix.m_nRows = static_cast<int> (sliceArg[0].Vector().size());
                lhsMatrix.m_nCols = 1;
            }

            hwMathStatus status = lhsMatrix.m_values.Dimension(1, m_values.Type());
            lhsMatrix.m_rows.empty();
            lhsMatrix.m_pointerB.resize(lhsMatrix.m_nCols);
            lhsMatrix.m_pointerE.resize(lhsMatrix.m_nCols);
            lhsMatrix.m_pointerB[0] = 0;
            int nnz = 0;

            for (int ii = 0; ii < sliceArg[0].Vector().size(); ++ii)
            {
                int index = sliceArg[0].Vector()[ii];
                int row = index % m_nRows;
                int col = index / m_nRows;
                int col_nnz = 0;

                if (col < 0 || col >= lhsMatrix.m_nCols)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                {
                    if (m_rows[k] < row)
                        continue;

                    if (m_rows[k] > row)
                        break;

                    lhsMatrix.m_rows.push_back(0);
                    ++nnz;
                    ++col_nnz;
                    status = lhsMatrix.m_values.Resize(nnz);

                    if (m_values.IsReal())
                    {
                        lhsMatrix.m_values(nnz - 1) = m_values(k);
                    }
                    else
                    {
                        lhsMatrix.m_values.z(nnz - 1) = m_values.z(k);
                    }
                }

                lhsMatrix.m_pointerE[ii] = col_nnz;    // not cummulative
            }

            for (int ii = 1; ii < lhsMatrix.m_nCols; ++ii)
            {
                lhsMatrix.m_pointerE[ii] += lhsMatrix.m_pointerE[ii - 1];
                lhsMatrix.m_pointerB[ii]  = lhsMatrix.m_pointerE[ii - 1];
            }

            if (!nnz)
            {
                status = lhsMatrix.m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);     // lhs is empty
            }
        }
    }
    else if (sliceArg.size() == 2)
    {
        if (sliceArg[0].IsScalar())
        {
            int row = sliceArg[0].Scalar();

            if (row < 0 || row >= m_nRows)
                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();
                throw hwMathException(HW_MATH_ERR_NOTIMPLEMENT);
            }
            else if (sliceArg[1].IsColon())
            {
                lhsMatrix.m_nRows = 1;
                lhsMatrix.m_nCols = m_nCols;
                hwMathStatus status = lhsMatrix.m_values.Dimension(1, m_values.Type());
                lhsMatrix.m_rows.empty();
                lhsMatrix.m_pointerB.resize(m_nCols);
                lhsMatrix.m_pointerE.resize(m_nCols);
                lhsMatrix.m_pointerB[0] = 0;
                int nnz = 0;

                for (int col = 0; col < m_nCols; ++col)
                {
                    int col_nnz = 0;

                    for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                    {
                        if (m_rows[k] < row)
                            continue;

                        if (m_rows[k] > row)
                            break;

                        lhsMatrix.m_rows.push_back(0);
                        ++nnz;
                        ++col_nnz;
                        status = lhsMatrix.m_values.Resize(nnz);

                        if (m_values.IsReal())
                        {
                            lhsMatrix.m_values(nnz - 1) = m_values(k);
                        }
                        else
                        {
                            lhsMatrix.m_values.z(nnz - 1) = m_values.z(k);
                        }

                        lhsMatrix.m_pointerE[col] = col_nnz;    // not cummulative
                    }
                }

                if (!nnz)
                {
                    status = lhsMatrix.m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);     // lhs is empty
                }
            }
            else if (sliceArg[1].IsVector())
            {
                lhsMatrix.m_nRows = 1;
                lhsMatrix.m_nCols = static_cast<int> (sliceArg[1].Vector().size());
                hwMathStatus status = lhsMatrix.m_values.Dimension(1, m_values.Type());
                lhsMatrix.m_rows.empty();
                lhsMatrix.m_pointerB.resize(lhsMatrix.m_nCols);
                lhsMatrix.m_pointerE.resize(lhsMatrix.m_nCols);
                lhsMatrix.m_pointerB[0] = 0;
                int nnz = 0;

                for (int index = 0; index < lhsMatrix.m_nCols; ++index)
                {
                    int col     = sliceArg[1].Vector()[index];
                    int col_nnz = 0;

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                    {
                        if (m_rows[k] < row)
                            continue;

                        if (m_rows[k] > row)
                            break;

                        lhsMatrix.m_rows.push_back(0);
                        ++nnz;
                        ++col_nnz;
                        status = lhsMatrix.m_values.Resize(nnz);

                        if (m_values.IsReal())
                        {
                            lhsMatrix.m_values(nnz - 1) = m_values(k);
                        }
                        else
                        {
                            lhsMatrix.m_values.z(nnz - 1) = m_values.z(k);
                        }
                    }

                    lhsMatrix.m_pointerE[index] = col_nnz;    // not cummulative
                }

                if (!nnz)
                {
                    status = lhsMatrix.m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);     // lhs is empty
                }
            }

            for (int ii = 1; ii < lhsMatrix.m_nCols; ++ii)
            {
                lhsMatrix.m_pointerE[ii] += lhsMatrix.m_pointerE[ii - 1];
                lhsMatrix.m_pointerB[ii]  = lhsMatrix.m_pointerE[ii - 1];
            }
        }
        else if (sliceArg[0].IsColon())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();

                if (col < 0 || col >= m_nCols)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                int start = m_pointerB[col];
                int nnz = m_pointerE[col] - start;
                lhsMatrix.m_nRows = m_nRows;
                lhsMatrix.m_nCols = 1;
                hwMathStatus status = lhsMatrix.m_values.Dimension(nnz, m_values.Type());
                lhsMatrix.m_rows.resize(nnz);
                lhsMatrix.m_pointerB.resize(1);
                lhsMatrix.m_pointerE.resize(1);
                lhsMatrix.m_pointerB[0] = 0;
                lhsMatrix.m_pointerE[0] = nnz;

                for (int index = 0; index < nnz; ++index)
                {
                    if (m_values.IsReal())
                    {
                        lhsMatrix.m_values(index) = m_values(index + start);
                    }
                    else
                    {
                        lhsMatrix.m_values.z(index) = m_values.z(index + start);
                    }

                    lhsMatrix.m_rows[index] = m_rows[index + start];
                }
            }
            else if (sliceArg[1].IsColon())
            {
                lhsMatrix.m_nRows    = m_nRows;
                lhsMatrix.m_nCols    = m_nCols;
                lhsMatrix.m_values   = m_values;
                lhsMatrix.m_rows     = m_rows;
                lhsMatrix.m_pointerB = m_pointerB;
                lhsMatrix.m_pointerE = m_pointerE;
            }
            else if (sliceArg[1].IsVector())
            {
                lhsMatrix.m_nRows = m_nRows;
                lhsMatrix.m_nCols = static_cast<int> (sliceArg[1].Vector().size());
                hwMathStatus status = lhsMatrix.m_values.Dimension(1, m_values.Type());
                lhsMatrix.m_rows.empty();
                lhsMatrix.m_pointerB.resize(lhsMatrix.m_nCols);
                lhsMatrix.m_pointerE.resize(lhsMatrix.m_nCols);
                lhsMatrix.m_pointerB[0] = 0;
                int nnz = 0;

                for (int index = 0; index < lhsMatrix.m_nCols; ++index)
                {
                    int col = sliceArg[1].Vector()[index];

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                    {
                        lhsMatrix.m_rows.push_back(m_rows[k]);
                        ++nnz;
                        status = lhsMatrix.m_values.Resize(nnz);

                        if (m_values.IsReal())
                        {
                            lhsMatrix.m_values(nnz - 1) = m_values(k);
                        }
                        else
                        {
                            lhsMatrix.m_values.z(nnz - 1) = m_values.z(k);
                        }
                    }

                    lhsMatrix.m_pointerE[index] = m_pointerE[col] - m_pointerB[col];    // not cummulative
                }

                for (int ii = 1; ii < lhsMatrix.m_nCols; ++ii)
                {
                    lhsMatrix.m_pointerE[ii] += lhsMatrix.m_pointerE[ii - 1];
                    lhsMatrix.m_pointerB[ii]  = lhsMatrix.m_pointerE[ii - 1];
                }

                if (!nnz)
                {
                    status = lhsMatrix.m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);     // lhs is empty
                }
            }
        }
        else if (sliceArg[0].IsVector())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();
                lhsMatrix.m_nRows = static_cast<int> (sliceArg[0].Vector().size());
                lhsMatrix.m_nCols = 1;
                hwMathStatus status = lhsMatrix.m_values.Dimension(1, m_values.Type());
                lhsMatrix.m_rows.empty();
                lhsMatrix.m_pointerB.resize(1);
                lhsMatrix.m_pointerE.resize(1);
                lhsMatrix.m_pointerB[0] = 0;
                int nnz = 0;

                if (col < 0 || col >= m_nCols)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                for (int index = 0; index < lhsMatrix.m_nRows; ++index)
                {
                    int row = sliceArg[0].Vector()[index];

                    if (row < 0 || row >= m_nRows)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                    {
                        if (m_rows[k] < row)
                            continue;

                        if (m_rows[k] > row)
                            break;

                        lhsMatrix.m_rows.push_back(nnz);
                        ++nnz;
                        status = lhsMatrix.m_values.Resize(nnz);

                        if (m_values.IsReal())
                        {
                            lhsMatrix.m_values(nnz - 1) = m_values(k);
                        }
                        else
                        {
                            lhsMatrix.m_values.z(nnz - 1) = m_values.z(k);
                        }
                    }
                }

                lhsMatrix.m_pointerE[0] = nnz;

                if (!nnz)
                {
                    status = lhsMatrix.m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);     // lhs is empty
                }
            }
            else if (sliceArg[1].IsColon())
            {
                lhsMatrix.m_nRows = static_cast<int> (sliceArg[0].Vector().size());
                lhsMatrix.m_nCols = m_nCols;
                hwMathStatus status = lhsMatrix.m_values.Dimension(1, m_values.Type());
                lhsMatrix.m_rows.empty();
                lhsMatrix.m_pointerB.resize(m_nCols);
                lhsMatrix.m_pointerE.resize(m_nCols);
                lhsMatrix.m_pointerB[0] = 0;
                int nnz = 0;

                for (int col = 0; col < lhsMatrix.m_nCols; ++col)
                {
                    int col_nnz = 0;

                    for (int index = 0; index < lhsMatrix.m_nRows; ++index)
                    {
                        int row = sliceArg[0].Vector()[index];

                        for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                        {
                            if (m_rows[k] < row)
                                continue;

                            if (m_rows[k] > row)
                                break;

                            lhsMatrix.m_rows.push_back(index);
                            ++nnz;
                            ++col_nnz;
                            status = lhsMatrix.m_values.Resize(nnz);

                            if (m_values.IsReal())
                            {
                                lhsMatrix.m_values(nnz - 1) = m_values(k);
                            }
                            else
                            {
                                lhsMatrix.m_values.z(nnz - 1) = m_values.z(k);
                            }
                        }
                    }

                    lhsMatrix.m_pointerE[col] = col_nnz;    // not cummulative
                }

                for (int ii = 1; ii < m_nCols; ++ii)
                {
                    lhsMatrix.m_pointerE[ii] += lhsMatrix.m_pointerE[ii - 1];
                    lhsMatrix.m_pointerB[ii]  = lhsMatrix.m_pointerE[ii - 1];
                }

                if (!nnz)
                {
                    status = lhsMatrix.m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);     // lhs is empty
                }
            }
            else if (sliceArg[1].IsVector())
            {
                lhsMatrix.m_nRows = static_cast<int> (sliceArg[0].Vector().size());
                lhsMatrix.m_nCols = static_cast<int> (sliceArg[1].Vector().size());
                hwMathStatus status = lhsMatrix.m_values.Dimension(1, m_values.Type());
                lhsMatrix.m_rows.empty();
                lhsMatrix.m_pointerB.resize(lhsMatrix.m_nCols);
                lhsMatrix.m_pointerE.resize(lhsMatrix.m_nCols);
                lhsMatrix.m_pointerB[0] = 0;
                int nnz = 0;

                for (int colIndex = 0; colIndex < lhsMatrix.m_nCols; ++colIndex)
                {
                    int col = sliceArg[1].Vector()[colIndex];
                    int col_nnz = 0;

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    for (int rowIndex = 0; rowIndex < lhsMatrix.m_nRows; ++rowIndex)
                    {
                        int row = sliceArg[0].Vector()[rowIndex];

                        if (row < 0 || row >= m_nRows)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                        {
                            if (m_rows[k] < row)
                                continue;

                            if (m_rows[k] > row)
                                break;

                            lhsMatrix.m_rows.push_back(rowIndex);
                            ++nnz;
                            ++col_nnz;
                            status = lhsMatrix.m_values.Resize(nnz);

                            if (m_values.IsReal())
                            {
                                lhsMatrix.m_values(nnz - 1) = m_values(k);
                            }
                            else
                            {
                                lhsMatrix.m_values.z(nnz - 1) = m_values.z(k);
                            }
                        }
                    }

                    lhsMatrix.m_pointerE[colIndex] = col_nnz;    // not cummulative
                }

                for (int ii = 1; ii < lhsMatrix.m_nCols; ++ii)
                {
                    lhsMatrix.m_pointerE[ii] += lhsMatrix.m_pointerE[ii - 1];
                    lhsMatrix.m_pointerB[ii] = lhsMatrix.m_pointerE[ii - 1];
                }

                if (!nnz)
                {
                    status = lhsMatrix.m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);     // lhs is empty
                }
            }
        }
    }
    else
    {
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
    }
}

//! Write to a matrix slice of the calling object, as if the calling object is being
//! sliced on the left hand side of an equals sign
//! matrix(slice args) = rhs_matrix
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::SliceLHS(const std::vector<hwSliceArg>& sliceArg,
                                  const hwTMatrix<T1, T2>& rhsMatrix)
{
    if (!rhsMatrix.IsReal())
        MakeComplex();

    if (sliceArg.size() == 1)
    {
        if (sliceArg[0].IsScalar())
        {
            // int index = sliceArg[0].Scalar();
            throw hwMathException(HW_MATH_ERR_NOTIMPLEMENT);
        }
        else if (sliceArg[0].IsColon())
        {
            if (rhsMatrix.IsVector())
            {
                int size = Size();

                if (rhsMatrix.Size() == size)
                {
                    m_values = rhsMatrix;
                    m_rows.resize(size);
                    
                    for (int ii = 0; ii < size; ++ii)
                        m_rows[ii] = ii % m_nRows;

                    for (int ii = 0; ii < m_nCols; ++ii)
                    {
                        m_pointerE[ii] = m_nRows * (ii + 1);

                        if (ii < m_nCols - 1)
                            m_pointerB[ii + 1] = m_pointerE[ii];
                    }

                    if (IsReal())
                    {
                        for (int ii = 0; ii < size; ++ii)
                        {
                            if (m_values(ii) == 0.0)
                                ZeroElement(ii);
                        }
                    }
                    else
                    {
                        for (int ii = 0; ii < size; ++ii)
                        {
                            if (m_values.z(ii) == 0.0)
                                ZeroElement(ii);
                        }
                    }
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                }
            }
            else if (rhsMatrix.Is0x0())
            {
                m_nRows = 0;
                m_nCols = 0;
                m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);
                m_rows.empty();
                m_pointerB.resize(0);
                m_pointerE.resize(0);
            }
            else
            {
                throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
            }
        }
        else if (sliceArg[0].IsVector())
        {
            int vecSize = static_cast<int> (sliceArg[0].Vector().size());

            if (rhsMatrix.IsVector() && rhsMatrix.Size() == vecSize)
            {
                int size = Size();

                for (int ii = 0; ii < vecSize; ++ii)
                {
                    int index = sliceArg[0].Vector()[ii];

                    if (index < 0 || index >= size)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    if (rhsMatrix.IsReal())
                    {
                        if (rhsMatrix(ii) == 0.0)
                            ZeroElement(index);
                        else if (IsReal())
                            operator()(index) = rhsMatrix(ii);
                        else
                            z(index) = rhsMatrix(ii);
                    }
                    else
                    {
                        if (rhsMatrix.z(ii) == 0.0)
                            ZeroElement(index);
                        else
                            z(index) = rhsMatrix.z(ii);
                    }
                }
            }
            else
            {
                throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
            }
        }
    }
    else if (sliceArg.size() == 2)
    {
        if (sliceArg[0].IsScalar())
        {
            int row = sliceArg[0].Scalar();

            if (row < 0 || row >= m_nRows)
                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

            if (sliceArg[1].IsScalar())
            {
                // int col = sliceArg[1].Scalar();
                throw hwMathException(HW_MATH_ERR_NOTIMPLEMENT);
            }
            else if (sliceArg[1].IsColon())
            {
                if (rhsMatrix.IsVector())
                {
                    if (m_nCols == rhsMatrix.M() || m_nCols == rhsMatrix.N())
                    {
                        for (int ii = 0; ii < m_nCols; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii) == 0.0)
                                    ZeroElement(row, ii);
                                else if (IsReal())
                                    operator()(row, ii) = rhsMatrix(ii);
                                else
                                    z(row, ii) = rhsMatrix(ii);
                            }
                            else
                            {
                                if (rhsMatrix.z(ii) == 0.0)
                                    ZeroElement(row, ii);
                                else
                                    z(row, ii) = rhsMatrix.z(ii);
                            }
                        }
                    }
                    else
                    {
                        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                    }
                }
                else if (rhsMatrix.Is0x0())
                {
                    for (int col = 0; col < m_nCols; ++col)
                    {
                        for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                        {
                            if (m_rows[k] < row)
                                continue;

                            if (m_rows[k] > row)
                            {
                                --m_rows[k];
                                continue;
                            }

                            m_values.DeleteElements(k);
                            m_rows.erase(m_rows.begin() + k);

                            for (int jj = col+1; jj < m_nCols; ++jj)
                                --m_pointerB[jj];

                            for (int jj = col; jj < m_nCols; ++jj)
                                --m_pointerE[jj];

                            --k;
                        }
                    }

                    --m_nRows;
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                }
            }
            else if (sliceArg[1].IsVector())
            {
                if (rhsMatrix.IsVector())
                {
                    if (sliceArg[1].Vector().size() == rhsMatrix.Size())
                    {
                        for (int colIndex = 0; colIndex < rhsMatrix.Size(); ++colIndex)
                        {
                            int col = sliceArg[1].Vector()[colIndex];

                            if (col < 0 || col >= m_nCols)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(colIndex) == 0.0)
                                    ZeroElement(row, col);
                                else if (IsReal())
                                    operator()(row, col) = rhsMatrix(colIndex);
                                else
                                    z(row, col) = rhsMatrix(colIndex);
                            }
                            else
                            {
                                if (rhsMatrix.z(colIndex) == 0.0)
                                    ZeroElement(row, col);
                                else
                                    z(row, col) = rhsMatrix.z(colIndex);
                            }
                        }
                    }
                    else
                    {
                        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                    }
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                }
            }
        }
        else if (sliceArg[0].IsColon())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();

                if (col < 0 || col >= m_nCols)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                if (rhsMatrix.IsVector())
                {
                    if (m_nRows == rhsMatrix.M() || m_nRows == rhsMatrix.N())
                    {
                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii) == 0.0)
                                    ZeroElement(ii, col);
                                else if (IsReal())
                                    operator()(ii, col) = rhsMatrix(ii);
                                else
                                    z(ii, col) = rhsMatrix(ii);
                            }
                            else
                            {
                                if (rhsMatrix.z(ii) == 0.0)
                                    ZeroElement(ii, col);
                                else
                                    z(ii, col) = rhsMatrix.z(ii);
                            }
                        }
                    }
                }
                else if (rhsMatrix.Is0x0())
                {
                    int nnzCol = m_pointerE[col] - m_pointerB[col];

                    for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                    {
                        m_values.DeleteElements(k);
                        m_rows.erase(m_rows.begin() + k);
                        --m_pointerE[col];
                        --k;
                    }

                    m_pointerB.erase(m_pointerB.begin() + col);
                    m_pointerE.erase(m_pointerE.begin() + col);
                    --m_nCols;

                    for (int k = col; k < m_nCols; ++k)
                    {
                        m_pointerB[col] -= nnzCol;
                        m_pointerE[col] -= nnzCol;
                    }
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                }
            }
            else if (sliceArg[1].IsColon())
            {
                if (rhsMatrix.M() == m_nRows && rhsMatrix.N() == m_nCols)
                {
                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii, jj) == 0.0)
                                    ZeroElement(ii, jj);
                                else if (IsReal())
                                    operator()(ii, jj) = rhsMatrix(ii, jj);
                                else
                                    z(ii, jj) = rhsMatrix(ii, jj);
                            }
                            else
                            {
                                if (rhsMatrix.z(ii, jj) == 0.0)
                                    ZeroElement(ii, jj);
                                else
                                    z(ii, jj) = rhsMatrix.z(ii, jj);
                            }
                        }
                    }
                }
                else if (rhsMatrix.Is0x0())
                {
                    m_nRows = 0;
                    m_nCols = 0;
                    m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);
                    m_rows.empty();
                    m_pointerB.resize(0);
                    m_pointerE.resize(0);
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                }
            }
            else if (sliceArg[1].IsVector())
            {
                if (rhsMatrix.M() == m_nRows && rhsMatrix.N() == sliceArg[1].Vector().size())
                {
                    for (int jj = 0; jj < rhsMatrix.N(); ++jj)
                    {
                        int col = sliceArg[1].Vector()[jj];

                        if (col < 0 || col >= m_nCols)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii, jj) == 0.0)
                                    ZeroElement(ii, col);
                                else if (IsReal())
                                    operator()(ii, col) = rhsMatrix(ii, jj);
                                else
                                    z(ii, col) = rhsMatrix(ii, jj);
                            }
                            else
                            {
                                if (rhsMatrix.z(ii, jj) == 0.0)
                                    ZeroElement(ii, col);
                                else
                                    z(ii, col) = rhsMatrix.z(ii, jj);
                            }
                        }
                    }
                }
                else if (rhsMatrix.Is0x0())
                {
                    std::vector<int> cols = sliceArg[1].Vector();
                    std::sort(cols.begin(), cols.end());

                    for (int jj = 0; jj < cols.size(); ++jj)
                    {
                        int col = cols[jj] - jj;

                        if (col < 0 || col >= m_nCols)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        int nnzCol = m_pointerE[col] - m_pointerB[col];

                        for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                        {
                            m_values.DeleteElements(k);
                            m_rows.erase(m_rows.begin() + k);
                            --m_pointerE[col];
                            --k;
                        }

                        m_pointerB.erase(m_pointerB.begin() + col);
                        m_pointerE.erase(m_pointerE.begin() + col);
                        --m_nCols;

                        for (int k = col; k < m_nCols; ++k)
                        {
                            m_pointerB[k] -= nnzCol;
                            m_pointerE[k] -= nnzCol;
                        }
                    }
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                }
            }
        }
        else if (sliceArg[0].IsVector())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();

                if (rhsMatrix.IsVector())
                {
                    if (sliceArg[0].Vector().size() == rhsMatrix.Size())
                    {
                        for (int rowIndex = 0; rowIndex < rhsMatrix.Size(); ++rowIndex)
                        {
                            int row = sliceArg[0].Vector()[rowIndex];

                            if (row < 0 || row >= m_nRows)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(rowIndex) == 0.0)
                                    ZeroElement(row, col);
                                else if (IsReal())
                                    operator()(row, col) = rhsMatrix(rowIndex);
                                else
                                    z(row, col) = rhsMatrix(rowIndex);
                            }
                            else
                            {
                                if (rhsMatrix.z(rowIndex) == 0.0)
                                    ZeroElement(row, col);
                                else
                                    z(row, col) = rhsMatrix.z(rowIndex);
                            }
                        }
                    }
                    else
                    {
                        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                    }
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                }
            }
            else if (sliceArg[1].IsColon())
            {
                if (rhsMatrix.M() == sliceArg[0].Vector().size() && rhsMatrix.N() == m_nCols)
                {
                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        for (int ii = 0; ii < rhsMatrix.M(); ++ii)
                        {
                            int row = sliceArg[0].Vector()[ii];

                            if (row < 0 || row >= m_nRows)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii, jj) == 0.0)
                                    ZeroElement(row, jj);
                                else if (IsReal())
                                    operator()(row, jj) = rhsMatrix(ii, jj);
                                else
                                    z(row, jj) = rhsMatrix(ii, jj);
                            }
                            else
                            {
                                if (rhsMatrix.z(ii, jj) == 0.0)
                                    ZeroElement(row, jj);
                                else
                                    z(row, jj) = rhsMatrix.z(ii, jj);
                            }
                        }
                    }
                }
                else if (rhsMatrix.Is0x0())
                {
                    std::vector<int> rows = sliceArg[0].Vector();
                    std::sort(rows.begin(), rows.end());

                    for (int ii = 0; ii < rows.size(); ++ii)
                    {
                        int row = rows[ii] - ii;

                        if (row < 0 || row >= m_nRows)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        for (int col = 0; col < m_nCols; ++col)
                        {
                            for (int k = m_pointerB[col]; k < m_pointerE[col]; ++k)
                            {
                                if (m_rows[k] < row)
                                    continue;

                                if (m_rows[k] > row)
                                {
                                    --m_rows[k];
                                    continue;
                                }

                                m_values.DeleteElements(k);
                                m_rows.erase(m_rows.begin() + k);

                                for (int jj = col + 1; jj < m_nCols; ++jj)
                                    --m_pointerB[jj];

                                for (int jj = col; jj < m_nCols; ++jj)
                                    --m_pointerE[jj];

                                --k;
                            }
                        }

                        --m_nRows;
                    }
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                }
            }
            else if (sliceArg[1].IsVector())
            {
                if (rhsMatrix.M() == sliceArg[0].Vector().size() && rhsMatrix.N() == sliceArg[1].Vector().size())
                {
                    for (int jj = 0; jj < rhsMatrix.N(); ++jj)
                    {
                        int col = sliceArg[1].Vector()[jj];

                        if (col < 0 || col >= m_nCols)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        for (int ii = 0; ii < rhsMatrix.M(); ++ii)
                        {
                            int row = sliceArg[0].Vector()[ii];

                            if (row < 0 || row >= m_nRows)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii, jj) == 0.0)
                                    ZeroElement(row, col);
                                else if (IsReal())
                                    operator()(row, col) = rhsMatrix(ii, jj);
                                else
                                    z(row, col) = rhsMatrix(ii, jj);
                            }
                            else
                            {
                                if (rhsMatrix.z(ii, jj) == 0.0)
                                    ZeroElement(row, col);
                                else
                                    z(row, col) = rhsMatrix.z(ii, jj);
                            }
                        }
                    }
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_ARRAYSIZE);
                }
            }
        }
    }
    else
    {
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
    }
}

//! Write to a matrix slice of the calling object, as if the calling object is being
//! sliced on the left hand side of an equals sign
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::SliceLHS(const std::vector<hwSliceArg>& sliceArg, T1 real)
{
    if (!IsReal())
    {
        return SliceLHS(sliceArg, hwTComplex<T1>(real, static_cast<T1> (0)));
    }

    if (sliceArg.size() == 1)
    {
        if (sliceArg[0].IsScalar())
        {
            int index = sliceArg[0].Scalar();

            if (real != static_cast<T1> (0))
                operator()(index) = real;
            else
                ZeroElement(index);
        }
        else if (sliceArg[0].IsColon())
        {
            if (real != static_cast<T1> (0))
            {
                int size = Size();
                hwMathStatus status = m_values.Dimension(size, 1, hwTMatrix<T1, T2>::REAL);
                m_rows.resize(size);

                m_values.SetElements(real);

                for (int ii = 0; ii < size; ++ii)
                    m_rows[ii] = ii % m_nRows;

                for (int ii = 0; ii < m_nCols; ++ii)
                {
                    m_pointerE[ii] = m_nRows * (ii + 1);

                    if (ii < m_nCols - 1)
                        m_pointerB[ii + 1] = m_pointerE[ii];
                }
            }
            else
            {
                hwMathStatus status = m_values.Dimension(0, 0, hwTMatrix<T1, T2>::REAL);
                m_rows.clear();

                for (int ii = 0; ii < m_nCols; ++ii)
                {
                    m_pointerE[ii] = 0;

                    if (ii < m_nCols - 1)
                        m_pointerB[ii + 1] = m_pointerE[ii];
                }
            }
        }
        else if (sliceArg[0].IsVector())
        {
            int vecSize = static_cast<int> (sliceArg[0].Vector().size());
            int size = Size();

            for (int ii = 0; ii < vecSize; ++ii)
            {
                int index = static_cast<int> (sliceArg[0].Vector()[ii]);

                if (index < 0 || index >= size)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                if (real != static_cast<T1> (0))
                    operator()(index) = real;
                else
                    ZeroElement(index);
            }
        }
    }
    else if (sliceArg.size() == 2)
    {
        if (sliceArg[0].IsScalar())
        {
            int row = sliceArg[0].Scalar();

            if (row < 0 || row >= m_nRows)
                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();

                if (real != static_cast<T1> (0))
                    operator()(row, col) = real;
                else
                    ZeroElement(row, col);
            }
            else if (sliceArg[1].IsColon())
            {
                for (int ii = 0; ii < m_nCols; ++ii)
                {
                    if (real != static_cast<T1> (0))
                        operator()(row, ii) = real;
                    else
                        ZeroElement(row, ii);
                }
            }
            else if (sliceArg[1].IsVector())
            {
                for (int colIndex = 0; colIndex < sliceArg[1].Vector().size(); ++colIndex)
                {
                    int col = sliceArg[1].Vector()[colIndex];

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    if (real != static_cast<T1> (0))
                        operator()(row, col) = real;
                    else
                        ZeroElement(row, col);
                }
            }
        }
        else if (sliceArg[0].IsColon())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();

                if (col < 0 || col >= m_nCols)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                for (int ii = 0; ii < m_nCols; ++ii)
                {
                    if (real != static_cast<T1> (0))
                        operator()(ii, col) = real;
                    else
                        ZeroElement(ii, col);
                }
            }
            else if (sliceArg[1].IsColon())
            {
                for (int jj = 0; jj < m_nCols; ++jj)
                {
                    if (real != static_cast<T1> (0))
                    {
                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            operator()(ii, jj) = real;
                        }
                    }
                    else
                    {
                        hwMathStatus status = m_values.Dimension(0, 0, hwTMatrix<T1, T2>::REAL);
                        m_rows.clear();

                        for (int ii = 0; ii < m_nCols; ++ii)
                        {
                            m_pointerE[ii] = 0;

                            if (ii < m_nCols - 1)
                                m_pointerB[ii + 1] = m_pointerE[ii];
                        }
                    }
                }
            }
            else if (sliceArg[1].IsVector())
            {
                for (int jj = 0; jj < sliceArg[1].Vector().size(); ++jj)
                {
                    int col = sliceArg[1].Vector()[jj];

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    for (int ii = 0; ii < m_nRows; ++ii)
                    {
                        if (real != static_cast<T1> (0))
                            operator()(ii, col) = real;
                        else
                            ZeroElement(ii, col);
                    }
                }
            }
        }
        else if (sliceArg[0].IsVector())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();

                for (int rowIndex = 0; rowIndex < sliceArg[0].Vector().size(); ++rowIndex)
                {
                    int row = sliceArg[0].Vector()[rowIndex];

                    if (row < 0 || row >= m_nRows)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    if (real != static_cast<T1> (0))
                        operator()(row, col) = real;
                    else
                        ZeroElement(row, col);
                }
            }
            else if (sliceArg[1].IsColon())
            {
                for (int jj = 0; jj < m_nCols; ++jj)
                {
                    for (int ii = 0; ii < sliceArg[0].Vector().size(); ++ii)
                    {
                        int row = sliceArg[0].Vector()[ii];

                        if (row < 0 || row >= m_nRows)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        if (real != static_cast<T1> (0))
                            operator()(row, jj) = real;
                        else
                            ZeroElement(row, jj);
                    }
                }
            }
            else if (sliceArg[1].IsVector())
            {
                for (int jj = 0; jj < sliceArg[1].Vector().size(); ++jj)
                {
                    int col = sliceArg[1].Vector()[jj];

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    for (int ii = 0; ii < sliceArg[0].Vector().size(); ++ii)
                    {
                        int row = sliceArg[0].Vector()[ii];

                        if (row < 0 || row >= m_nRows)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        if (real != static_cast<T1> (0))
                            operator()(row, col) = real;
                        else
                            ZeroElement(row, col);
                    }
                }
            }
        }
    }
    else
    {
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
    }
}

//! Write to a matrix slice of the calling object, as if the calling object is being
//! sliced on the left hand side of an equals sign
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::SliceLHS(const std::vector<hwSliceArg>& sliceArg, const T2& cmplx)
{
    if (cmplx == static_cast<T1> (0))
    {
        return SliceLHS(sliceArg, static_cast<T1> (0));
    }

    if (IsReal())
    {
        MakeComplex();
    }

    if (sliceArg.size() == 1)
    {
        if (sliceArg[0].IsScalar())
        {
            int index = sliceArg[0].Scalar();
            z(index) = cmplx;
        }
        else if (sliceArg[0].IsColon())
        {
            int size = Size();

            hwMathStatus status = m_values.Dimension(size, 1, hwTMatrix<T1, T2>::COMPLEX);
            m_rows.resize(size);

            m_values.SetElements(cmplx);

            for (int ii = 0; ii < size; ++ii)
                m_rows[ii] = ii % m_nRows;

            for (int ii = 0; ii < m_nCols; ++ii)
            {
                m_pointerE[ii] = m_nRows * (ii + 1);

                if (ii < m_nCols - 1)
                    m_pointerB[ii + 1] = m_pointerE[ii];
            }
        }
        else if (sliceArg[0].IsVector())
        {
            int vecSize = static_cast<int> (sliceArg[0].Vector().size());
            int size = Size();

            for (int ii = 0; ii < vecSize; ++ii)
            {
                int index = static_cast<int> (sliceArg[0].Vector()[ii]);

                if (index < 0 || index >= size)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                z(index) = cmplx;
            }
        }
    }
    else if (sliceArg.size() == 2)
    {
        if (sliceArg[0].IsScalar())
        {
            int row = sliceArg[0].Scalar();

            if (row < 0 || row >= m_nRows)
                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();
                z(row, col) = cmplx;
            }
            else if (sliceArg[1].IsColon())
            {
                for (int ii = 0; ii < m_nCols; ++ii)
                {
                    z(row, ii) = cmplx;
                }
            }
            else if (sliceArg[1].IsVector())
            {
                for (int colIndex = 0; colIndex < sliceArg[1].Vector().size(); ++colIndex)
                {
                    int col = sliceArg[1].Vector()[colIndex];

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    z(row, col) = cmplx;
                }
            }
        }
        else if (sliceArg[0].IsColon())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();

                if (col < 0 || col >= m_nCols)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                for (int ii = 0; ii < m_nCols; ++ii)
                {
                    z(col, ii) = cmplx;
                }
            }
            else if (sliceArg[1].IsColon())
            {
                for (int jj = 0; jj < m_nCols; ++jj)
                {
                    for (int ii = 0; ii < m_nRows; ++ii)
                    {
                        z(ii, jj) = cmplx;
                    }
                }
            }
            else if (sliceArg[1].IsVector())
            {
                for (int jj = 0; jj < sliceArg[1].Vector().size(); ++jj)
                {
                    int col = sliceArg[1].Vector()[jj];

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    for (int ii = 0; ii < m_nRows; ++ii)
                    {
                        z(ii, col) = cmplx;
                    }
                }
            }
        }
        else if (sliceArg[0].IsVector())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();

                for (int rowIndex = 0; rowIndex < sliceArg[0].Vector().size(); ++rowIndex)
                {
                    int row = sliceArg[0].Vector()[rowIndex];

                    if (row < 0 || row >= m_nRows)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    z(row, col) = cmplx;
                }
            }
            else if (sliceArg[1].IsColon())
            {
                for (int jj = 0; jj < m_nCols; ++jj)
                {
                    for (int ii = 0; ii < sliceArg[0].Vector().size(); ++ii)
                    {
                        int row = sliceArg[0].Vector()[ii];

                        if (row < 0 || row >= m_nRows)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        z(row, jj) = cmplx;
                    }
                }
            }
            else if (sliceArg[1].IsVector())
            {
                for (int jj = 0; jj < sliceArg[1].Vector().size(); ++jj)
                {
                    int col = sliceArg[1].Vector()[jj];

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    for (int ii = 0; ii < sliceArg[0].Vector().size(); ++ii)
                    {
                        int row = sliceArg[0].Vector()[ii];

                        if (row < 0 || row >= m_nRows)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        z(row, col) = cmplx;
                    }
                }
            }
        }
    }
    else
    {
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
    }
}

// ****************************************************
//                   Submatrix Operations
// ****************************************************

//! Write a submatrix source to the calling object, starting at the specified location
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::WriteSubmatrix(int startRow, int startCol, const hwTMatrixS<T1, T2>& source)
{
    if (startRow < 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);

    if (startCol < 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 2);

    if (source.IsReal())
    {
        if (!IsReal())
        {
            hwTMatrixS<T1, T2> temp;
            temp.PackComplex(source);

            return WriteSubmatrix(startRow, startCol, temp);
        }
    }
    else    // source is complex
    {
        if (IsReal())
            MakeComplex();
    }

    int sourceM = source.M();
    int sourceN = source.N();
    int targetM = m_nRows;
    int targetN = m_nCols;

    if (startRow + sourceM > targetM)
        targetM = startRow + sourceM;

    if (startCol + sourceN > targetN)
        targetN = startCol + sourceN;

    int old_cols = m_nCols;
    int nnz = static_cast<int> (m_rows.size());

    m_nRows = targetM;
    m_nCols = targetN;
    m_pointerB.resize(m_nCols);
    m_pointerE.resize(m_nCols);

    for (int i = old_cols; i < m_nCols; ++i)
    {
        m_pointerB[i] = nnz;
        m_pointerE[i] = nnz;
    }

    int col = 0;
    int k;

    for (int index = 0; index < source.m_rows.size(); ++index)
    {
        for ( ; col < source.m_nCols; ++col)
        {
            if (index < source.m_pointerE[col])
                break;
        }

        int ii = startRow + source.m_rows[index];
        int jj = startCol + col;
        bool loc_exists = false;
    
        for (k = m_pointerB[jj]; k < m_pointerE[jj]; ++k)
        {
            if (m_rows[k] < ii)
                continue;
    
            if (m_rows[k] > ii)
                break;
    
            loc_exists = true;
            break;
        }
    
        if (!loc_exists)
        {
            // insert a location to be reassigned
            hwMathStatus status;

            if (m_values.IsEmpty())
            {
                if (source.IsReal())
                    status = m_values.Dimension(k + 1, hwTMatrix<T1, T2>::REAL);
                else
                    status = m_values.Dimension(k + 1, hwTMatrix<T1, T2>::COMPLEX);

                m_rows.resize(1);
                m_rows[0] = ii;
            }
            else
            {
                hwTMatrix<T1, T2> copy(m_values);
                status = m_values.InsertElements(copy, k);
                std::vector<int>::iterator it = m_rows.begin() + k;
                m_rows.insert(it, ii);
            }
    
            for (; jj < m_nCols - 1; ++jj)
            {
                ++m_pointerB[jj + 1];
                ++m_pointerE[jj];
            }
    
            ++m_pointerE[m_nCols - 1];
        }
    
        if (IsReal())
            m_values(k) = source.m_values(index);
        else
            m_values.z(k) = source.m_values.z(index);
    }
}

//! Write a submatrix source to the calling object, starting at the specified location
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::WriteSubmatrix(int startRow, int startCol, const hwTMatrix<T1, T2>& source)
{
    if (startRow < 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);

    if (startCol < 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 2);

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

    int sourceM = source.M();
    int sourceN = source.N();
    int targetM = m_nRows;
    int targetN = m_nCols;

    if (startRow + sourceM > targetM)
        targetM = startRow + sourceM;

    if (startCol + sourceN > targetN)
        targetN = startCol + sourceN;

    int old_cols = m_nCols;
    int nnz = static_cast<int> (m_rows.size());

    m_nRows = targetM;
    m_nCols = targetN;
    m_pointerB.resize(m_nCols);
    m_pointerE.resize(m_nCols);

    for (int i = old_cols; i < m_nCols; ++i)
    {
        m_pointerB[i] = nnz;
        m_pointerE[i] = nnz;
    }

    for (int j = 0; j < source.N(); ++j)
    {
        int k = 0;

        for (int i = 0; i < source.M(); ++i)
        {
            if (IsReal())
            {
                if (source(i, j) == static_cast<T1> (0))
                    continue;
            }
            else if (source.z(i, j) == static_cast<T1> (0))
            {
                continue;
            }

            int ii = startRow + i;
            int jj = startCol + j;
            bool loc_exists = false;

            for (k = m_pointerB[jj]; k < m_pointerE[jj]; ++k)
            {
                if (m_rows[k] < ii)
                    continue;

                if (m_rows[k] > ii)
                    break;

                loc_exists = true;
                break;
            }

            if (!loc_exists)
            {
                // insert a location to be reassigned
                hwMathStatus status;

                if (m_values.IsEmpty())
                {
                    status = m_values.Dimension(k + 1, source.Type());
                    m_rows.resize(1);
                    m_rows[0] = ii;
                }
                else
                {
                    hwTMatrix<T1, T2> copy(m_values);
                    status = m_values.InsertElements(copy, k);
                    std::vector<int>::iterator it = m_rows.begin() + k;
                    m_rows.insert(it, ii);
                }

                for (; jj < m_nCols - 1; ++jj)
                {
                    ++m_pointerB[jj + 1];
                    ++m_pointerE[jj];
                }

                ++m_pointerE[m_nCols - 1];
            }

            if (IsReal())
                m_values(k) = source(i, j);
            else
                m_values.z(k) = source.z(i, j);
        }
    }
}

//! Concatenate two matrices vertically
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::ConcatVertical(const hwTMatrixS<T1, T2>& A, const hwTMatrixS<T1, T2>& B)
{
    if (A.IsReal() && !B.IsReal())
    {
        hwMatrixS AC;
        AC.PackComplex(A);
        return ConcatVertical(AC, B);
    }
    else if (!A.IsReal() && B.IsReal())
    {
        hwMatrixS BC;
        BC.PackComplex(B);
        return ConcatVertical(A, BC);
    }

    m_nCols = A.N();

    if (B.N() != m_nCols)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    m_nRows = A.M() + B.M();
    m_pointerB.resize(m_nCols);
    m_pointerE.resize(m_nCols);

    for (int i = 0; i < m_nCols; ++i)
    {
        m_pointerB[i] = A.m_pointerB[i] + B.m_pointerB[i];
        m_pointerE[i] = A.m_pointerE[i] + B.m_pointerE[i];
    }

    m_rows.reserve(m_pointerE.back());
    std::vector<int>::const_iterator itA = A.m_rows.begin();
    int am = A.M();

    if (A.IsReal())
    {
        m_values.Dimension(m_pointerE.back(), 1, hwTMatrix<T1, T2>::REAL);
        const T1* realA = A.m_values.GetRealData();
        const T1* realB = B.m_values.GetRealData();
        T1* realC = m_values.GetRealData();
        int nnz = 0;
        int nnzB = 0;

        for (int i = 0; i < m_nCols; ++i)
        {
            int nnzCol = A.m_pointerE[i] - A.m_pointerB[i];
            int numBytes = nnzCol * sizeof(T1);
            nnz += nnzCol;
            m_rows.insert(m_rows.end(), itA, itA + nnzCol);
            itA += nnzCol;

            memcpy_s(realC, numBytes, realA, numBytes);
            realC += nnzCol;
            realA += nnzCol;

            nnzCol = B.m_pointerE[i] - B.m_pointerB[i];
            numBytes = nnzCol * sizeof(T1);
            m_rows.resize(nnz + nnzCol);

            for (int k = 0; k < nnzCol; ++k)
            {
                m_rows[nnz + k] = B.m_rows[nnzB + k] + am;
            }

            nnz += nnzCol;
            nnzB += nnzCol;

            memcpy_s(realC, numBytes, realB, numBytes);
            realC += nnzCol;
            realB += nnzCol;
        }
    }
    else    // complex
    {
        m_values.Dimension(m_pointerE.back(), 1, hwTMatrix<T1, T2>::COMPLEX);
        const T2* cmplxA = A.m_values.GetComplexData();
        const T2* cmplxB = B.m_values.GetComplexData();
        T2* cmplxC = m_values.GetComplexData();
        int nnz = 0;
        int nnzB = 0;

        for (int i = 0; i < m_nCols; ++i)
        {
            int nnzCol = A.m_pointerE[i] - A.m_pointerB[i];
            int numBytes = nnzCol * sizeof(T2);
            nnz += nnzCol;
            m_rows.insert(m_rows.end(), itA, itA + nnzCol);
            itA += nnzCol;

            memcpy_s(cmplxC, numBytes, cmplxA, numBytes);
            cmplxC += nnzCol;
            cmplxA += nnzCol;

            nnzCol = B.m_pointerE[i] - B.m_pointerB[i];
            numBytes = nnzCol * sizeof(T2);
            m_rows.resize(nnz + nnzCol);

            for (int k = 0; k < nnzCol; ++k)
            {
                m_rows[nnz + k] = B.m_rows[nnzB + k] + am;
            }

            nnz += nnzCol;
            nnzB += nnzCol;

            memcpy_s(cmplxC, numBytes, cmplxB, numBytes);
            cmplxC += nnzCol;
            cmplxB += nnzCol;
        }
    }
}

//! Concatenate two matrices horizontally
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::ConcatHorizontal(const hwTMatrixS<T1, T2>& A, const hwTMatrixS<T1, T2>& B)
{
    if (A.IsReal() && !B.IsReal())
    {
        hwMatrixS AC;
        AC.PackComplex(A);
        return ConcatHorizontal(AC, B);
    }
    else if (!A.IsReal() && B.IsReal())
    {
        hwMatrixS BC;
        BC.PackComplex(B);
        return ConcatHorizontal(A, BC);
    }

    m_nRows = A.M();

    if (B.M() != m_nRows)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    int an = A.N();
    m_nCols = an + B.N();

    m_pointerB.reserve(m_nCols);
    m_pointerB = A.m_pointerB;
    m_pointerB.resize(m_nCols);

    m_pointerE.reserve(m_nCols);
    m_pointerE = A.m_pointerE;
    m_pointerE.resize(m_nCols);

    int nnzA = A.NNZ();

    for (int i = an; i < m_nCols; ++i)
    {
        m_pointerB[i] = B.m_pointerB[i - an] + nnzA;
        m_pointerE[i] = B.m_pointerE[i - an] + nnzA;
    }

    m_rows.reserve(m_pointerE.back());
    m_rows = A.m_rows;
    m_rows.insert(m_rows.end(), B.m_rows.begin(), B.m_rows.end());

    if (A.IsReal())
    {
        m_values.Dimension(m_pointerE.back(), 1, hwTMatrix<T1, T2>::REAL);
        m_values.WriteSubmatrix(0, 0, A.m_values);
        m_values.WriteSubmatrix(nnzA, 0, B.m_values);
    }
    else    // complex
    {
        m_values.Dimension(m_pointerE.back(), 1, hwTMatrix<T1, T2>::COMPLEX);
        m_values.WriteSubmatrix(0, 0, A.m_values);
        m_values.WriteSubmatrix(nnzA, 0, B.m_values);
    }
}

//! Expand and assign memory of a dimensioned sparse matrix, column-wise
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::FillEmptyDimensionedRows(int startRow, const hwTMatrixS<T1, T2>& B)
{
    // This function is similar to ConcatVertical except that
    // *this already has the appropriate output dimensions upon input
    // Warning: startRow must be one more that the last populated row of A
    hwTMatrixS<T1, T2> A(*this);
    
    if (startRow < 0 || startRow + B.M() > m_nRows)
        throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);
    
    int nRows_save = m_nRows;
    A.m_nRows = startRow;
    MakeEmpty();
    ConcatVertical(A, B);
    m_nRows = nRows_save;
}

//! Expand and assign memory of a dimensioned sparse matrix, row-wise
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::FillEmptyDimensionedColumns(int startCol, const hwTMatrixS<T1, T2>& B)
{
    // This function is similar to ConcatHorizontal, except that
    // *this already has the appropriate output dimensions upon input
    // Warning: startCol must be one more that the last populated column of A
    hwTMatrixS<T1, T2> A(*this);
    
    if (startCol < 0 || startCol + B.N() > m_nCols)
        throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);
    
    int nCols_save = m_nCols;
    A.m_nCols = startCol;
    A.m_pointerB.resize(startCol);
    A.m_pointerE.resize(startCol);
    MakeEmpty();
    ConcatHorizontal(A, B);
    m_nCols = nCols_save;
    m_pointerB.resize(m_nCols);
    m_pointerE.resize(m_nCols);
}

// ****************************************************
//           Misc Linear Algebra Operations
// ****************************************************

//! Generate a diagonal matrix or extract a diagonal from a matrix
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Diag(const hwTMatrixS<T1, T2>& source, int k)
{
    typename hwTMatrix<T1, T2>::DataType type;

    if (source.IsReal())
        type = hwTMatrix<T1, T2>::REAL;
    else
        type = hwTMatrix<T1, T2>::COMPLEX;

    if (source.IsVector())
    {
        int length = source.Size();
        int m = length + abs(k);
        int n = m;
        int nnz = source.NNZ();
        std::vector<int> ivec(nnz);
        std::vector<int> jvec(nnz);
        hwTMatrix<T1, T2> diag(nnz, 1, type);
        int row;
        int col;

        if (source.IsReal())   // real
        {
            T1 real;

            if (k >= 0)
            {
                for (int i = 0; i < nnz; ++i)
                {
                    source.NZinfo(i, row, col, real);
                    ivec[i] = col;
                    jvec[i] = col + k;
                    diag(i) = real;
                }
            }
            else
            {
                for (int i = 0; i < nnz; ++i)
                {
                    source.NZinfo(i, row, col, real);
                    ivec[i] = col - k;
                    jvec[i] = col;
                    diag(i) = real;
                }
            }
        }
        else                   // complex
        {
            T2 cmplx;

            if (k >= 0)
            {
                for (int i = 0; i < nnz; ++i)
                {
                    source.NZinfo(i, row, col, cmplx);
                    ivec[i] = col;
                    jvec[i] = col + k;
                    diag.z(i) = cmplx;
                }
            }
            else
            {
                for (int i = 0; i < nnz; ++i)
                {
                    source.NZinfo(i, row, col, cmplx);
                    ivec[i] = col - k;
                    jvec[i] = col;
                    diag.z(i) = cmplx;
                }
            }

        }

        hwTMatrixS temp(ivec, jvec, diag, m, n);
        (*this) = temp;
    }
    else    // matrix
    {
        hwMathStatus status;
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

        std::vector<int> ivec;
        std::vector<int> jvec;
        hwTMatrix<T1, T2> diag;
        int nnz = 0;

        if (length > 0)
        {
            status = diag.Dimension(1, 1, type);

            if (!status.IsOk())
            {
                throw hwMathException(HW_MATH_ERR_OUTOFMEMORY);
            }

            if (source.IsReal())   // real
            {
                if (k >= 0)
                {
                    for (int i = 0; i < length; i++)
                    {
                        T1 real = source(i, i + k);

                        if (real != static_cast<T1> (0))
                        {
                            ivec.push_back(i);
                            jvec.push_back(0);
                            ++nnz;
                            diag.Resize(nnz, 1);
                            diag(nnz - 1) = real;
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < length; i++)
                    {
                        T1 real = source(i - k, i);

                        if (real != static_cast<T1> (0))
                        {
                            ivec.push_back(i);
                            jvec.push_back(0);
                            ++nnz;
                            diag.Resize(nnz, 1);
                            diag(nnz - 1) = real;
                        }
                    }
                }
            }
            else                // complex
            {
                if (k >= 0)
                {
                    for (int i = 0; i < length; i++)
                    {
                        const T2& cplx = source.z(i, i + k);

                        if (cplx != static_cast<T1> (0))
                        {
                            ivec.push_back(i);
                            jvec.push_back(0);
                            ++nnz;
                            diag.Resize(nnz, 1);
                            diag.z(nnz - 1) = cplx;
                        }
                    }
                }
                else
                {
                    for (int i = 0; i < length; i++)
                    {
                        const T2& cplx = source.z(i - k, i);

                        if (cplx != static_cast<T1> (0))
                        {
                            ivec.push_back(i);
                            jvec.push_back(0);
                            ++nnz;
                            diag.Resize(nnz, 1);
                            diag.z(nnz - 1) = cplx;
                        }
                    }
                }
            }

            if (!nnz)
                status = diag.Dimension(0, hwTMatrix<T1, T2>::REAL);

            hwTMatrixS temp(ivec, jvec, diag, length, 1);
            (*this) = temp;
        }
        else if (m == 0 || n == 0)
        {
            MakeEmpty();
        }
        else
        {
            MakeEmpty();
            m_nCols = 1;
        }
    }
}

//! Generate a diagonal matrix from a vector
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Diag(const hwTMatrixS<T1, T2>& source, int m, int n)
{
    typename hwTMatrix<T1, T2>::DataType type;

    if (source.IsReal())
        type = hwTMatrix<T1, T2>::REAL;
    else
        type = hwTMatrix<T1, T2>::COMPLEX;

    if (!source.IsVector())
        throw hwMathException(HW_MATH_ERR_VECTOR, 1);

    if (m < 0)
        throw hwMathException(HW_MATH_ERR_ARRAYDIM, 2);

    if (n < 0)
        throw hwMathException(HW_MATH_ERR_ARRAYDIM, 3);

    int nnz = source.NNZ();
    std::vector<int> ivec;
    std::vector<int> jvec;
    hwTMatrix<T1, T2> diag;
    int row;
    int col;
    hwMathStatus status;

    if (source.IsReal())   // real
    {
        T1 real;

        for (int i = 0; i < nnz; ++i)
        {
            source.NZinfo(i, row, col, real);

            if (col >= m || col >= n)
            {
                break;
            }

            ivec.push_back(col);
            jvec.push_back(col);

            if (i == 0)
                status = diag.Dimension(1, 1, type);
            else
                status = diag.Resize(i + 1, 1);

            diag(i) = real;
        }
    }
    else                   // complex
    {
        T2 cmplx;

        for (int i = 0; i < nnz; ++i)
        {
            source.NZinfo(i, row, col, cmplx);

            if (col >= m || col >= n)
                break;

            ivec.push_back(col);
            jvec.push_back(col);

            if (i == 0)
                status = diag.Dimension(1, 1, type);
            else
                status = diag.Resize(i + 1, 1);

            diag.z(i) = cmplx;
        }
    }

    hwTMatrixS temp(ivec, jvec, diag, m, n);
    (*this) = temp;
}

//! Transpose the matrix in place
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Transpose()
{
    // This can probably be done with mkl_dcsrcsc and mkl_zcsrcsc
    int nnz = static_cast<int>(m_rows.size());
    std::vector<int> jvec(nnz);
    int col = 0;

    for (int col = 0; col < m_nCols; ++col)
    {
        for (int ii = m_pointerB[col]; ii < m_pointerE[col]; ++ii)
        {
            jvec[ii] = col;
        }
    }

    hwTMatrixS<T1, T2> trans(jvec, m_rows, m_values, m_nCols, m_nRows);
    *this = trans;
}

//! Transpose the matrix argument
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Transpose(const hwTMatrixS<T1, T2>& source)
{
    // This can probably be done with mkl_dcsrcsc and mkl_zcsrcsc
    int nnz = static_cast<int> (source.m_rows.size());
    std::vector<int> jvec(nnz);
    int col = 0;

    for (int col = 0; col < source.m_nCols; ++col)
    {
        for (int ii = source.m_pointerB[col]; ii < source.m_pointerE[col]; ++ii)
        {
            jvec[ii] = col;
        }
    }

    hwTMatrixS<T1, T2> trans(jvec, source.m_rows, source.m_values, source.m_nCols, source.m_nRows);
    *this = trans;
}

//! Conjugate the matrix in place
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Conjugate()
{
    m_values.Conjugate();
}

//! Conjugate the matrix argument
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Conjugate(const hwTMatrixS<T1, T2>& source)
{
    *this = source;
    m_values.Conjugate();
}

//! Transpose and conjugate the matrix in place
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Hermitian()
{
    Transpose();
    Conjugate();
}

//! Transpose and conjugate the matrix argument
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Hermitian(const hwTMatrixS<T1, T2>& source)
{
    Transpose(source);
    Conjugate();
}

//! Sum along a direction
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Sum(const hwTMatrixS<T1, T2>& source, bool cols)
{
    hwMathStatus status;
    hwTMatrix<T1, T2> sum;

    if (cols)   // sum each column
    {
        int index = 0;

        if (source.IsReal())
        {
            status = sum.Dimension(1, source.m_nCols, hwTMatrix<T1, T2>::REAL);
            sum.SetElements(0.0);

            for (int col = 0; col < source.m_nCols; ++col)
            {
                for (int row = source.m_pointerB[col]; row < source.m_pointerE[col]; ++row)
                    sum(0, col) += source(source.m_rows[index++], col);
            }
        }
        else
        {
            status = sum.Dimension(1, source.m_nCols, hwTMatrix<T1, T2>::COMPLEX);
            sum.SetElements(0.0);

            for (int col = 0; col < source.m_nCols; ++col)
            {
                for (int row = source.m_pointerB[col]; row < source.m_pointerE[col]; ++row)
                    sum.z(0, col) += source.z(source.m_rows[index++], col);
            }
        }
    }
    else   // sum each row
    {
        int nnz = source.NNZ();

        if (source.IsReal())
        {
            status = sum.Dimension(source.m_nRows, 1, hwTMatrix<T1, T2>::REAL);
            sum.SetElements(0.0);

            for (int index = 0; index < nnz; ++index)
                sum(source.m_rows[index], 0) += source.m_values(index);
        }
        else
        {
            status = sum.Dimension(source.m_nRows, 1, hwTMatrix<T1, T2>::COMPLEX);
            sum.SetElements(0.0);

            for (int index = 0; index < nnz; ++index)
                sum.z(source.m_rows[index], 0) += source.m_values.z(index);
        }
    }

    hwTMatrixS<T1, T2> temp(sum);
    (*this) = temp;
}

// ****************************************************
//               Arithmetic Operations
// ****************************************************

//! Add two matrices, sum.Add(A,B)
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Add(const hwTMatrixS<T1, T2>& A, const hwTMatrixS<T1, T2>& B)
{
    if (A.NNZ() < B.NNZ())
        return Add(B, A);

    if (this == &A)
        throw hwMathException(HW_MATH_ERR_NOTIMPLEMENT);
    if (this == &B)
        throw hwMathException(HW_MATH_ERR_NOTIMPLEMENT);

    // check dimensions
    int m_nRows = A.m_nRows;
    int m_nCols = A.m_nCols;

    if (A.m_nRows != B.m_nRows || A.m_nCols != B.m_nCols)
    {
        if (A.Size() == 0 && B.Size() == 1)
        {
        }
        else if (A.Size() == 1 && B.Size() == 0)
        {
            m_nRows = B.m_nRows;
            m_nCols = B.m_nCols;
        }
        else
        {
            throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);
        }

        return;
    }

    (*this) = A;

    if (IsReal() && !B.IsReal())
        MakeComplex();

    int col = 0;
    int k;

    for (int index = 0; index < B.m_rows.size(); ++index)
    {
        for (; col < B.m_nCols; ++col)
        {
            if (index < B.m_pointerE[col])
                break;
        }

        int ii = B.m_rows[index];
        int jj = col;
        bool loc_exists = false;

        for (k = m_pointerB[jj]; k < m_pointerE[jj]; ++k)
        {
            if (m_rows[k] < ii)
                continue;

            if (m_rows[k] > ii)
                break;

            if (IsReal())
            {
                m_values(k) += B.m_values(index);

                if (m_values(k) == 0.0)
                    ZeroElement(m_rows[k], jj);
            }
            else if (B.IsReal())
            {
                m_values.z(k) += B.m_values(index);

                if (m_values.z(k) == 0.0)
                    ZeroElement(m_rows[k], jj);
            }
            else
            {
                m_values.z(k) += B.m_values.z(index);

                if (m_values.z(k) == 0.0)
                    ZeroElement(m_rows[k], jj);
            }

            loc_exists = true;
            break;
        }

        if (!loc_exists)
        {
            // insert a location to be reassigned
            hwTMatrix<T1, T2> copy(m_values);
            hwMathStatus status = m_values.InsertElements(copy, k);
            std::vector<int>::iterator it = m_rows.begin() + k;
            m_rows.insert(it, ii);

            for (; jj < m_nCols - 1; ++jj)
            {
                ++m_pointerB[jj + 1];
                ++m_pointerE[jj];
            }

            ++m_pointerE[m_nCols - 1];

            if (IsReal())
                m_values(k)   = B.m_values(index);
            else if (B.IsReal())
                m_values.z(k) = B.m_values(index);
            else
                m_values.z(k) = B.m_values.z(index);
        }
    }
}
//! Subtract two sparse matrices, diff.Subtr(A,B)
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Subtr(const hwTMatrixS<T1, T2>& A, const hwTMatrixS<T1, T2>& B)
{
    hwTMatrixS<T1, T2> C;
    C.Negate(B);
    return Add(A, C);
}

/*
// This function temporarily resides in /oml/BuiltInFuncs.* to avoid propogating
// the MKL header #includes into all client projects.
// Multiply a sparse matrix by a full matrix, prod = A * B
template<>
inline void hwTMatrixS<double>::Mult(const hwTMatrix<double>& B, hwTMatrix<double>& prod) const
{
    hwTMatrixS<double>& A = const_cast<hwTMatrixS<double>&> (*this);

    // get dimensions info
    int m = A.m_nRows;
    int n = B.N();
    int k = A.m_nCols;
    hwMathStatus status;

    if (B.M() != k)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 0, 1);

    int* pr = const_cast<int*> (A.rows());
    int* pb = const_cast<int*> (A.pointerB());
    int* pe = const_cast<int*> (A.pointerE());

    // MKL matrix and description
    sparse_matrix_t A_MKL;
    struct matrix_descr DSC;
    DSC.type = SPARSE_MATRIX_TYPE_GENERAL;
    sparse_status_t mkl_status;

    if (A.IsReal() && B.IsReal())
    {
        status = prod.Dimension(A.m_nRows, B.N(), hwTMatrix<double>::REAL);

        double* pV = const_cast<double*> (A.GetRealData());

        // populate MKL matrix
        mkl_status = mkl_sparse_d_create_csc(&A_MKL, SPARSE_INDEX_BASE_ZERO,
            m_nRows, m_nCols, pb, pe, pr, pV);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        // perform MKL operation
        if (n == 1)
        {
            mkl_status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,
                1.0, A_MKL, DSC, B.GetRealData(), 0.0, prod.GetRealData());
        }
        else
        {
            mkl_status = mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE,
                1.0, A_MKL, DSC, SPARSE_LAYOUT_COLUMN_MAJOR,
                B.GetRealData(), n, m, 0.0, prod.GetRealData(), m);
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        status = prod.Dimension(A.m_nRows, B.N(), hwTMatrix<double>::COMPLEX);

        hwTComplex<double>* pV = const_cast<hwTComplex<double>*> (A.GetComplexData());

        // populate MKL matrix
        mkl_status = mkl_sparse_z_create_csc(&A_MKL, SPARSE_INDEX_BASE_ZERO,
            m_nRows, m_nCols, pb, pe, pr, reinterpret_cast<MKL_Complex16*> (pV));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        // perform MKL operation
        hwTComplex<double>* pB = const_cast<hwTComplex<double>*> (B.GetComplexData());
        MKL_Complex16*      pP = reinterpret_cast<MKL_Complex16*> (prod.GetComplexData());
        struct _MKL_Complex16 one;
        struct _MKL_Complex16 zero;

        one.real = 1.0;
        one.imag = 0.0;
        zero.real = 0.0;
        zero.imag = 0.0;

        if (n == 1)
        {
            mkl_status = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE,
                one, A_MKL, DSC, reinterpret_cast<MKL_Complex16*> (pB), zero, pP);
        }
        else
        {
            mkl_status = mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE,
                one, A_MKL, DSC, SPARSE_LAYOUT_COLUMN_MAJOR,
                reinterpret_cast<MKL_Complex16*> (pB), n, m, zero, pP, m);
        }
    }
    else if (!A.IsReal())
    {
        hwTMatrix<double> C;
        C.PackComplex(B);
        return Mult(C, prod);
    }
    else    // !B.IsReal()
    {
        hwTMatrixS<double> C;
        C.PackComplex(A);
        return C.Mult(B, prod);
    }

    if (mkl_status != SPARSE_STATUS_SUCCESS)
        throw hwMathException(HW_MATH_ERR_UNKNOWN);
}
*/
//! Multiply a matrix and a real number
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Mult(const hwTMatrixS<T1, T2>& A, T1 real)
{
    hwMathStatus status;

    if (real == static_cast<T1> (0))
    {
        status = m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);
        m_rows.clear();
        m_pointerB.clear();
        m_pointerB.resize(m_nCols);
        m_pointerE.clear();
        m_pointerE.resize(m_nCols);
    }
    else
    {
        m_nRows    = A.m_nRows;
        m_nCols    = A.m_nCols;
        m_rows     = A.m_rows;
        m_pointerB = A.m_pointerB;
        m_pointerE = A.m_pointerE;
        status = m_values.Mult(A.m_values, real);
    }

    if (!status.IsOk())
        throw hwMathException(status.GetMsgCode());
}

//! Multiply a matrix and a complex number
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Mult(const hwTMatrixS<T1, T2>& A, const T2& cmplx)
{
    hwMathStatus status;

    if (cmplx == static_cast<T1> (0))
    {
        status = m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);
        m_rows.clear();
        m_pointerB.clear();
        m_pointerB.resize(m_nCols);
        m_pointerE.clear();
        m_pointerE.resize(m_nCols);
    }
    else
    {
        m_nRows    = A.m_nRows;
        m_nCols    = A.m_nCols;
        m_rows     = A.m_rows;
        m_pointerB = A.m_pointerB;
        m_pointerE = A.m_pointerE;
        status = m_values.Mult(A.m_values, cmplx);
    }

    if (!status.IsOk())
        throw hwMathException(status.GetMsgCode());
}
/*
// This function temporarily resides in /oml/BuiltInFuncs.* to avoid propogating
// the MKL header #includes into all client projects.
// Divide a full matrix by a sparse matrix on the left side, Q = (*this) \ B
template<>
inline void hwTMatrixS<double>::DivideLeft(const hwTMatrix<double>& B, hwTMatrix<double>& Q) const
{
    const hwTMatrixS<double>& A = (*this);
    hwMathStatus status;

    if (A.m_nRows != A.m_nCols)
        throw hwMathException(HW_MATH_ERR_MTXNOTSQUARE, 0);

    if (A.m_nRows != B.M())
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 0, 1);

    // get dimensions info
    // MKL matrix and description
    _MKL_DSS_HANDLE_t handle;
    _INTEGER_t error;
    MKL_INT opt = MKL_DSS_DEFAULTS;
//    MKL_INT opt = MKL_DSS_ZERO_BASED_INDEXING;
    MKL_INT sym = MKL_DSS_NON_SYMMETRIC;
    MKL_INT type = MKL_DSS_INDEFINITE;

    // initialize the solver
    error = dss_create(handle, opt);

    if (error != MKL_DSS_SUCCESS)
        throw hwMathException(HW_MATH_ERR_UNKNOWN);

    // define the non-zero structure of the matrix
    hwTMatrixS<double> AT;
    AT.Transpose(A);

    int* pr = const_cast<int*> (AT.rows());
    int* pb = const_cast<int*> (AT.pointerB());
    int* pe = const_cast<int*> (AT.pointerE());
    int nnz = NNZ();

    std::vector<int> colcount;
    std::vector<int> newrows;
    int base = 1;
    colcount.push_back(base);

    for (int ii = 0; ii < m_nCols; ++ii)
        colcount.push_back(pe[ii] + base);

    for (int ii = 0; ii < nnz; ++ii)
        newrows.push_back(pr[ii] + base);

    int* pcc = colcount.data();
    int* prr = newrows.data();

    error = dss_define_structure(handle, sym, pcc, m_nRows, m_nCols, prr, nnz);

    if (error != MKL_DSS_SUCCESS)
        throw hwMathException(HW_MATH_ERR_UNKNOWN);

    // reorder the matrix
    error = dss_reorder(handle, opt, 0);

    if (error != MKL_DSS_SUCCESS)
        throw hwMathException(HW_MATH_ERR_UNKNOWN);

    if (A.IsReal() && B.IsReal())
    {
        status = Q.Dimension(A.m_nCols, B.N(), hwTMatrix<double>::REAL);

        // factor the matrix
        error = dss_factor_real(handle, type, AT.GetRealData());

        if (error != MKL_DSS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        // solve system
        int n = B.N();

        error = dss_solve_real(handle, opt, B.GetRealData(), n, Q.GetRealData());

        if (error != MKL_DSS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        status = Q.Dimension(A.m_nCols, B.N(), hwTMatrix<double>::COMPLEX);

        // factor the matrix
        error = dss_factor_complex(handle, type, AT.GetComplexData());

        if (error != MKL_DSS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        // solve system
        int n = B.N();

        error = dss_solve_complex(handle, opt, B.GetComplexData(), n,
            Q.GetComplexData());

        if (error != MKL_DSS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);
    }
    else if (!A.IsReal())
    {
        hwTMatrix<double> C;
        C.PackComplex(B);
        return DivideLeft(C, Q);
    }
    else    // !B.IsReal()
    {
        hwTMatrixS<double> C;
        C.PackComplex(A);
        return C.DivideLeft(B, Q);
    }

    // deallocate solver storage
    error = dss_delete(handle, opt);

    if (error != MKL_DSS_SUCCESS)
        throw hwMathException(HW_MATH_ERR_UNKNOWN);
}
*/
//! Divide a matrix by a real number
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Divide(const hwTMatrixS<T1, T2>& A, T1 real)
{
    m_nRows    = A.m_nRows;
    m_nCols    = A.m_nCols;
    m_rows     = A.m_rows;
    m_pointerB = A.m_pointerB;
    m_pointerE = A.m_pointerE;

    hwMathStatus status = m_values.Divide(A.m_values, real);

    if (!status.IsOk())
        throw hwMathException(status.GetMsgCode());
}

//! Divide a matrix by a complex number
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Divide(const hwTMatrixS<T1, T2>& A, const T2& cmplx)
{
    m_nRows    = A.m_nRows;
    m_nCols    = A.m_nCols;
    m_rows     = A.m_rows;
    m_pointerB = A.m_pointerB;
    m_pointerE = A.m_pointerE;

    hwMathStatus status = m_values.Divide(A.m_values, cmplx);

    if (!status.IsOk())
        throw hwMathException(status.GetMsgCode());
}

//! Negate a sparse matrix
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Negate(const hwTMatrixS<T1, T2>& A)
{
    m_nRows    = A.m_nRows;
    m_nCols    = A.m_nCols;
    m_rows     = A.m_rows;
    m_pointerB = A.m_pointerB;
    m_pointerE = A.m_pointerE;

    hwMathStatus status = m_values.Negate(A.m_values);

    if (!status.IsOk())
        throw hwMathException(status.GetMsgCode());
}
/*
//! Add a sparse matrix and a full matrix, A.Add(B,sum)
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Add(const hwTMatrix<T1, T2>& B, hwTMatrix<T1, T2>& sum) const
{
    return hwMathStatus();
}

//! Subtract two sparse matrices, sum.Subtr(A,B)
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::Subtr(const hwTMatrixS<T1, T2>& A, const hwTMatrixS<T1, T2>& B)
{
    return hwMathStatus();
}
*/
// ****************************************************
//               Arithmetic Operators
// ****************************************************
//! Implement the == operator
template<typename T1, typename T2>
bool hwTMatrixS<T1, T2>::operator==(const hwTMatrixS<T1, T2>& A) const
{
    // TODO:
    return false;
}
    
//! Implement the != operator
template<typename T1, typename T2>
bool hwTMatrixS<T1, T2>::operator!=(const hwTMatrixS<T1, T2>& A) const
{
    return !((*this) == A);
}

// ****************************************************
//                 Private Utilities
// ****************************************************

//! Set matrix to empty condition
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::MakeEmpty()
{
    m_nRows = 0;
    m_nCols = 0;
    m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);
    m_rows.clear();
    m_pointerB.resize(1);
    m_pointerE.resize(1);
    m_pointerB[0] = 0;
    m_pointerE[0] = 0;
}
