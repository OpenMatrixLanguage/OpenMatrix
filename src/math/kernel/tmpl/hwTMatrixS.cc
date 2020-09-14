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

// Functions that use the headers below temporarily reside
// in /oml/BuiltInFuncs.*. The may be relocated to /math/sparse
// in order to avoid propagating the MKL header dependencies
// to all clients.
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

// Comparison function to stable_sort pairs using first elements
template<typename T1, typename T2>
bool compare(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b)
{
    return a.first < b.first;
}

// Comparison function to sort pairs using fabs(first elements)
template<typename T1, typename T2>
bool compare2(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b)
{
	return fabs(a.first) > fabs(b.first);
}

// Comparison function to sort pairs using second elements
template<typename T1, typename T2>
bool compare3(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b)
{
	return a.second < b.second;
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
    long long int index_prev = -1;
    int col_prev = -1;
	int nDistinctIndices = 0;
	std::vector<int> colNum(nnz);

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

        stable_sort(pairs.begin(), pairs.end(), compare<long long int, T1>);

        for (int i = 0; i < nnz; ++i)
        {
            long long int index = pairs[i].first;

            if (index == index_prev)
            {
                if (option && !strcmp(option, "unique"))
                {
					m_values(nDistinctIndices - 1) = pairs[i].second;
				}
                else
                {
					m_values(nDistinctIndices - 1) += pairs[i].second;
				}
            }
            else
            {
				m_values(nDistinctIndices) = pairs[i].second;
				m_rows[nDistinctIndices] = static_cast<MKL_INT> (index% m_nRows);
                index_prev = index;
                int col = static_cast<int> (index / m_nRows);
				colNum[nDistinctIndices] = col;

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
				++nDistinctIndices;
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

        stable_sort(pairs.begin(), pairs.end(), compare<long long int, T2>);

        for (int i = 0; i < nnz; ++i)
        {
            long long int index = pairs[i].first;

            if (index == index_prev)
            {
                if (option && !strcmp(option, "unique"))
                {
					m_values.z(nDistinctIndices - 1) = pairs[i].second;
				}
				else
				{
					m_values.z(nDistinctIndices - 1) += pairs[i].second;
				}
            }
            else
            {
                m_values.z(nDistinctIndices) = pairs[i].second;
                m_rows[nDistinctIndices] = static_cast<int> (index % m_nRows);
                index_prev = index;
                int col = static_cast<int> (index / m_nRows);
				colNum[nDistinctIndices] = col;

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
				++nDistinctIndices;
			}
        }
    }

    if (nDistinctIndices != nnz && (option == nullptr || strcmp(option, "unique")))
    {
		// rebuild matrix, removing zero values
		std::vector<MKL_INT> nzCol(m_nCols);
		hwMatrix copyV(m_values);
		std::vector<MKL_INT> copyR = m_rows;
		m_rows.clear();
		m_rows.resize(nDistinctIndices);
		MKL_INT* srcR  = copyR.data();
		MKL_INT* destR = m_rows.data();
		int prevZero   = -1;

		nz = 0;

		if (m_values.IsReal())
		{
			T1* srcV = copyV.GetRealData();
			hwMathStatus status = m_values.Dimension(nDistinctIndices, 1, hwMatrix::REAL);
			T1* destV = m_values.GetRealData();

			for (int ii = 0; ii < nDistinctIndices; ++ii)
			{
				if (srcV[ii] == 0.0)
				{
					// copy data values between stored zeros
					std::size_t count = ii - prevZero - 1;
					std::size_t numBytes = count * sizeof(T1);

					if (count)
					{
						memcpy_s(destV, numBytes, srcV + prevZero + 1, numBytes);
						destV += count;
					}

					// copy row numbers between stored zeros
					numBytes = count * sizeof(MKL_INT);

					if (count)
					{
						memcpy_s(destR, numBytes, srcR + prevZero + 1, numBytes);
						destR += count;
					}

					// count stored zeros per column
					++nzCol[colNum[ii]];

					prevZero = ii;
					++nz;
				}
			}

			// copy remaining data values
			status = m_values.Resize(nDistinctIndices - nz, 1);
			std::size_t count = nDistinctIndices - 1 - prevZero;
			std::size_t numBytes = count * sizeof(T1);

			if (count)
				memcpy_s(destV, numBytes, srcV + prevZero + 1, numBytes);

			// copy remaining row numbers
			m_rows.resize(nDistinctIndices - nz);
			numBytes = count * sizeof(MKL_INT);

			if (count)
				memcpy_s(destR, numBytes, srcR + prevZero + 1, numBytes);
		}
		else
		{
			T2* srcV = copyV.GetComplexData();
			hwMathStatus status = m_values.Dimension(nDistinctIndices, 1, hwMatrix::COMPLEX);
			T2* destV = m_values.GetComplexData();

			for (int ii = 0; ii < nDistinctIndices; ++ii)
			{
				if (srcV[ii] == 0.0)
				{
					// copy data values between stored zeros
					std::size_t count = ii - prevZero - 1;
					std::size_t numBytes = count * sizeof(T2);

					if (count)
					{
						memcpy_s(destV, numBytes, srcV + prevZero + 1, numBytes);
						destV += count;
					}

					// copy row numbers between stored zeros
					numBytes = count * sizeof(MKL_INT);

					if (count)
					{
						memcpy_s(destR, numBytes, srcR + prevZero + 1, numBytes);
						destR += count;
					}

					// count stored zeros per column
					++nzCol[colNum[ii]];

					prevZero = ii;
					++nz;
				}
			}

			// copy remaining data values
			status = m_values.Resize(nDistinctIndices - nz, 1);
			std::size_t count = nDistinctIndices - 1 - prevZero;
			std::size_t numBytes = count * sizeof(T2);

			if (count)
				memcpy_s(destV, numBytes, srcV + prevZero + 1, numBytes);

			// copy remaining row numbers
			m_rows.resize(nDistinctIndices - nz);
			numBytes = count * sizeof(MKL_INT);

			if (count)
				memcpy_s(destR, numBytes, srcR + prevZero + 1, numBytes);
		}

		// adjust column counts
		for (int ii = 1; ii < m_nCols; ++ii)	// make zero count cumulative
			nzCol[ii] += nzCol[ii - 1];

		for (int ii = 1; ii < m_nCols; ++ii)
			m_pointerB[ii] -= nzCol[ii - 1];

		for (int ii = 0; ii < m_nCols; ++ii)
			m_pointerE[ii] -= nzCol[ii];

/*		// sort, prune, unsort (alternate method, leave for now)
		std::vector<std::pair<T1, long long int>> pairs(nDistinctIndices);

		for (int i = 0; i < nDistinctIndices; ++i)
		{
			pairs[i].first = m_values(i);
			pairs[i].second = colNum[i] * static_cast<long long int> (m_nRows) + m_rows[i];
		}

		sort(pairs.begin(), pairs.end(), compare2<T1, long long int>);

		int firstZero = -1;

		for (int i = nDistinctIndices - 1; i > -1; --i)
		{
			if (pairs[i].first != 0.0)
			{
				firstZero = i + 1;
				break;
			}
		}

		std::vector<int> nzCol(m_nCols);

		for (int i = firstZero; i < nDistinctIndices; ++i)
		{
			long long int index = pairs[i].second;
			int col = static_cast<int> (index / m_nRows);
			++nzCol[col];
		}

		pairs.resize(firstZero);
		sort(pairs.begin(), pairs.end(), compare3<T1, long long int>);

		hwMathStatus status = m_values.Dimension(firstZero, 1, hwMatrix::REAL);
		m_rows.clear();
		m_rows.resize(firstZero);

		for (int i = 0; i < firstZero; ++i)
		{
			long long int index = pairs[i].second;
			m_values(i) = pairs[i].first;
			m_rows[i] = static_cast<int> (index % m_nRows);
		}

		for (int ii = 1; ii < m_nCols; ++ii)	// make zero count cumulative
			nzCol[ii] += nzCol[ii - 1];

		// adjust column counts
		for (int ii = 1; ii < m_nCols; ++ii)
			m_pointerB[ii] -= nzCol[ii];

		for (int ii = 0; ii < m_nCols; ++ii)
			m_pointerE[ii] -= nzCol[ii];
*/
    }
}

//! Construct real sparse from MKL representation
template<typename T1, typename T2>
hwTMatrixS<T1, T2>::hwTMatrixS(int            nRows,
                               int            nCols,
                               const MKL_INT* pBeginCount,
                               const MKL_INT* pEndCount,
                               const MKL_INT* pRowNum,
                               const T1*      pValues)
{
    m_nRows = nRows;
    m_nCols = nCols;
    m_pointerB.resize(nCols);
    m_pointerE.resize(nCols);
    int nnz = *(pEndCount + nCols - 1);		// number of stored values

	// Note: MKL allows computed zero values to remain in the
	// output, but we exclude them for consistency with Octave.

	// Note: MKL does not guarantee that row numbers are sorted
	// within each column, but we do. It is possibile that Octave
	// only sorts row numbers when displayed. This issue needs to
	// be revisited.

	// count stored zeros
	std::vector<MKL_INT> nzCol(m_nCols);	// number of stored zeros in each column
	int nz = 0;								// total number of stored zeros
	int col = 0;

	for (int i = 0; i < nnz; ++i)
	{
		if (pValues[i] == 0.0)
		{
			for ( ; col < m_nCols; ++col)
			{
				if (i < pEndCount[col])
				{
					++nzCol[col];
					break;
				}
			}

			++nz;
		}
	}

	nnz -= nz;
	nz = 0;

	// populate data and sorted row vectors, removing stored zeros
    m_rows.resize(nnz);
    hwMathStatus status = m_values.Dimension(nnz, hwTMatrix<T1, T2>::REAL);

    for (int i = 0; i < nCols; ++i)
    {
        int start = *(pBeginCount + i);
        int nnzCol = *(pEndCount + i) - start;	// stored values per column
		int nz_col = 0;		// counter for stored zeros per column

        std::vector<std::pair<int, T1>> pairs(nnzCol - nzCol[i]);

        for (int j = 0; j < nnzCol; ++j)
        {
			if (pValues[start + j] != 0.0)
			{
				pairs[j - nz_col].first = pRowNum[start + j];
				pairs[j - nz_col].second = pValues[start + j];
			}
			else
			{
				++nz_col;
			}
        }

        sort(pairs.begin(), pairs.end(), compare<int, T1>);
		nnzCol -= nz_col;

        for (int j = 0; j < nnzCol; ++j)
        {
            m_rows[start + j - nz] = pairs[j].first;
            m_values(start + j - nz) = pairs[j].second;
        }

		nz += nz_col;
    }

	// populate column counter vectors
	std::size_t count = m_nCols * sizeof(MKL_INT);
	memcpy_s(m_pointerB.data(), count, pBeginCount, count);
	memcpy_s(m_pointerE.data(), count, pEndCount, count);

	for (int ii = 1; ii < m_nCols; ++ii)	// make zero count cumulative
		nzCol[ii] += nzCol[ii - 1];

	for (int ii = 1; ii < m_nCols; ++ii)
		m_pointerB[ii] -= nzCol[ii - 1];

	for (int ii = 0; ii < m_nCols; ++ii)
		m_pointerE[ii] -= nzCol[ii];
}

//! Construct complex sparse from MKL representation
template<typename T1, typename T2>
hwTMatrixS<T1, T2>::hwTMatrixS(int            nRows,
                               int            nCols,
                               const MKL_INT* pBeginCount,
                               const MKL_INT* pEndCount,
                               const MKL_INT* pRowNum,
                               const T2*      pValues)
{
	m_nRows = nRows;
	m_nCols = nCols;
	m_pointerB.resize(nCols);
	m_pointerE.resize(nCols);
	int nnz = *(pEndCount + nCols - 1);		// number of stored values

	// Note: MKL allows computed zero values to remain in the output, but
	// we exclude them for consistency with Octave.

	// Note: MKL does not guarantee that row numbers are sorted
	// within each column, but we do. It is possibile that Octave
	// only sorts row numbers when displayed. This issue needs to
	// be revisited.

	// count stored zeros
	std::vector<MKL_INT> nzCol(m_nCols);	// number of stored zeros in each column
	int nz = 0;								// total number of stored zeros
	int col = 0;

	for (int i = 0; i < nnz; ++i)
	{
		if (pValues[i] == 0.0)
		{
			for (; col < m_nCols; ++col)
			{
				if (i < pEndCount[col])
				{
					++nzCol[col];
					break;
				}
			}

			++nz;
		}
	}

	nnz -= nz;
	nz = 0;

	// populate data and sorted row vectors, removing stored zeros
	m_rows.resize(nnz);
	hwMathStatus status = m_values.Dimension(nnz, hwTMatrix<T1, T2>::COMPLEX);

	for (int i = 0; i < nCols; ++i)
	{
		int start = *(pBeginCount + i);
		int nnzCol = *(pEndCount + i) - start;	// stored values per column
		int nz_col = 0;		// counter for stored zeros per column

		std::vector<std::pair<int, T2>> pairs(nnzCol - nzCol[i]);

		for (int j = 0; j < nnzCol; ++j)
		{
			if (pValues[start + j] != 0.0)
			{
				pairs[j - nz_col].first = pRowNum[start + j];
				pairs[j - nz_col].second = pValues[start + j];
			}
			else
			{
				++nz_col;
			}
		}

		sort(pairs.begin(), pairs.end(), compare<int, T2>);
		nnzCol -= nz_col;

		for (int j = 0; j < nnzCol; ++j)
		{
			m_rows[start + j - nz] = pairs[j].first;
			m_values.z(start + j - nz) = pairs[j].second;
		}

		nz += nz_col;
	}

	// populate column counter vectors
	std::size_t count = m_nCols * sizeof(MKL_INT);
	memcpy_s(m_pointerB.data(), count, pBeginCount, count);
	memcpy_s(m_pointerE.data(), count, pEndCount, count);

	for (int ii = 1; ii < m_nCols; ++ii)	// make zero count cumulative
		nzCol[ii] += nzCol[ii - 1];

	for (int ii = 1; ii < m_nCols; ++ii)
		m_pointerB[ii] -= nzCol[ii - 1];

	for (int ii = 0; ii < m_nCols; ++ii)
		m_pointerE[ii] -= nzCol[ii];
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
            lhsMatrix.m_rows.clear();
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

                lhsMatrix.m_pointerE[ii] = col_nnz;    // not cumulative
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
                lhsMatrix.m_rows.clear();
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

                        lhsMatrix.m_pointerE[col] = col_nnz;    // not cumulative
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
                lhsMatrix.m_rows.clear();
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

                    lhsMatrix.m_pointerE[index] = col_nnz;    // not cumulative
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
                lhsMatrix.m_rows.clear();
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

                    lhsMatrix.m_pointerE[index] = m_pointerE[col] - m_pointerB[col];    // not cumulative
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
                lhsMatrix.m_rows.clear();
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
                lhsMatrix.m_rows.clear();
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

                    lhsMatrix.m_pointerE[col] = col_nnz;    // not cumulative
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
                lhsMatrix.m_rows.clear();
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

                    lhsMatrix.m_pointerE[colIndex] = col_nnz;    // not cumulative
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

                    RemoveStoredZeros();
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
                m_rows.clear();
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
                bool removeZeros = false;
                const hwTMatrixS<T1, T2>* const_this = this;

                for (int ii = 0; ii < vecSize; ++ii)
                {
                    int index = sliceArg[0].Vector()[ii];

                    if (index < 0 || index >= size)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    if (rhsMatrix.IsReal())
                    {
                        if (rhsMatrix(ii) != static_cast<T1> (0))
                        {
                            if (IsReal())
                                operator()(index) = rhsMatrix(ii);
                            else
                                z(index) = rhsMatrix(ii);
                        }
                        else if (IsReal())
                        {
                            if (const_this->operator()(index) != static_cast<T1> (0))
                            {
                                operator()(index) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                        else if (const_this->z(index) != static_cast<T1> (0))
                        {
                            z(index) = static_cast<T1> (0);
                            removeZeros = true;
                        }
                    }
                    else if (rhsMatrix.z(ii) != static_cast<T1> (0))
                    {
                        z(index) = rhsMatrix.z(ii);
                    }
                    else if (const_this->z(index) != static_cast<T1> (0))
                    {
                        z(index) = static_cast<T1> (0);
                        removeZeros = true;
                    }
                }

                if (removeZeros)
                    RemoveStoredZeros();
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
                        bool removeZeros = false;
                        const hwTMatrixS<T1, T2>* const_this = this;

                        for (int ii = 0; ii < m_nCols; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(row, ii) = rhsMatrix(ii);
                                    else
                                        z(row, ii) = rhsMatrix(ii);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(row, ii) != static_cast<T1> (0))
                                    {
                                        operator()(row, ii) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(row, ii) != static_cast<T1> (0))
                                {
                                    z(row, ii) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(ii) != static_cast<T1> (0))
                            {
                                z(row, ii) = rhsMatrix.z(ii);
                            }
                            else if (const_this->z(row, ii) != static_cast<T1> (0))
                            {
                                z(row, ii) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }

                        if (removeZeros)
                            RemoveStoredZeros();
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
                        bool removeZeros = false;
                        const hwTMatrixS<T1, T2>* const_this = this;

                        for (int colIndex = 0; colIndex < rhsMatrix.Size(); ++colIndex)
                        {
                            int col = sliceArg[1].Vector()[colIndex];

                            if (col < 0 || col >= m_nCols)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(colIndex) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(row, col) = rhsMatrix(colIndex);
                                    else
                                        z(row, col) = rhsMatrix(colIndex);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(row, col) != static_cast<T1> (0))
                                    {
                                        operator()(row, col) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(row, col) != static_cast<T1> (0))
                                {
                                    z(row, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else
                            {
                                if (rhsMatrix.z(colIndex) != static_cast<T1> (0))
                                {
                                    z(row, col) = rhsMatrix.z(colIndex);
                                }
                                else if (const_this->z(row, col) != static_cast<T1> (0))
                                {
                                    z(row, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                        }

                        if (removeZeros)
                            RemoveStoredZeros();
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
                        bool removeZeros = false;
                        const hwTMatrixS<T1, T2>* const_this = this;

                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(ii, col) = rhsMatrix(ii);
                                    else
                                        z(ii, col) = rhsMatrix(ii);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(ii, col) != static_cast<T1> (0))
                                    {
                                        operator()(ii, col) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(ii, col) != static_cast<T1> (0))
                                {
                                    z(ii, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else
                            {
                                if (rhsMatrix.z(ii) != static_cast<T1> (0))
                                {
                                    z(ii, col) = rhsMatrix(ii);
                                }
                                else if (const_this->z(ii, col) != static_cast<T1> (0))
                                {
                                    z(ii, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }

                            if (removeZeros)
                                RemoveStoredZeros();
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
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii, jj) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(ii, jj) = rhsMatrix(ii, jj);
                                    else
                                        z(ii, jj) = rhsMatrix(ii, jj);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(ii, jj) != static_cast<T1> (0))
                                    {
                                        operator()(ii, jj) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(ii, jj) != static_cast<T1> (0))
                                {
                                    z(ii, jj) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(ii, jj) != static_cast<T1> (0))
                            {
                                z(ii, jj) = rhsMatrix.z(ii, jj);
                            }
                            else if (const_this->z(ii, jj) != static_cast<T1> (0))
                            {
                                z(ii, jj) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
                }
                else if (rhsMatrix.Is0x0())
                {
                    m_nRows = 0;
                    m_nCols = 0;
                    m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);
                    m_rows.clear();
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
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < rhsMatrix.N(); ++jj)
                    {
                        int col = sliceArg[1].Vector()[jj];

                        if (col < 0 || col >= m_nCols)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii, jj) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(ii, col) = rhsMatrix(ii, jj);
                                    else
                                        z(ii, col) = rhsMatrix(ii, jj);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(ii, col) != static_cast<T1> (0))
                                    {
                                        operator()(ii, col) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(ii, col) != static_cast<T1> (0))
                                {
                                    z(ii, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(ii, jj) != static_cast<T1> (0))
                            {
                                z(ii, col) = rhsMatrix.z(ii, jj);
                            }
                            else if (const_this->z(ii, col) != static_cast<T1> (0))
                            {
                                z(ii, col) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
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
                        bool removeZeros = false;
                        const hwTMatrixS<T1, T2>* const_this = this;

                        for (int rowIndex = 0; rowIndex < rhsMatrix.Size(); ++rowIndex)
                        {
                            int row = sliceArg[0].Vector()[rowIndex];

                            if (row < 0 || row >= m_nRows)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(rowIndex) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(row, col) = rhsMatrix(rowIndex);
                                    else
                                        z(row, col) = rhsMatrix(rowIndex);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(row, col) != static_cast<T1> (0))
                                    {
                                        operator()(row, col) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(row, col) != static_cast<T1> (0))
                                {
                                    z(row, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(rowIndex) != static_cast<T1> (0))
                            {
                                z(row, col) = rhsMatrix.z(rowIndex);
                            }
                            else if (const_this->z(row, col) != static_cast<T1> (0))
                            {
                                z(row, col) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }

                        if (removeZeros)
                            RemoveStoredZeros();
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
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        for (int ii = 0; ii < rhsMatrix.M(); ++ii)
                        {
                            int row = sliceArg[0].Vector()[ii];

                            if (row < 0 || row >= m_nRows)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii, jj) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(row, jj) = rhsMatrix(ii, jj);
                                    else
                                        z(row, jj) = rhsMatrix(ii, jj);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(row, jj) != static_cast<T1> (0))
                                    {
                                        operator()(row, jj) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(row, jj) != static_cast<T1> (0))
                                {
                                    z(row, jj) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(ii, jj) != static_cast<T1> (0))
                            {
                                z(row, jj) = rhsMatrix.z(ii, jj);
                            }
                            else if (const_this->z(row, jj) != static_cast<T1> (0))
                            {
                                z(row, jj) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
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
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

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
                                if (rhsMatrix(ii, jj) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(row, col) = rhsMatrix(ii, jj);
                                    else
                                        z(row, col) = rhsMatrix(ii, jj);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(row, col) != static_cast<T1> (0))
                                    {
                                        operator()(row, col) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(row, col) != static_cast<T1> (0))
                                {
                                    z(row, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(ii, jj) != static_cast<T1> (0))
                            {
                                z(row, col) = rhsMatrix.z(ii, jj);
                            }
                            else if (const_this->z(row, col) != static_cast<T1> (0))
                            {
                                z(row, col) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
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
//! matrix(slice args) = rhs_matrix
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::SliceLHS(const std::vector<hwSliceArg>& sliceArg,
                                  const hwTMatrixS<T1, T2>& rhsMatrix)
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
                    m_values = rhsMatrix.m_values;
                    m_rows.resize(NNZ());
                    m_pointerB.clear();
                    m_pointerE.clear();
                    m_pointerB.resize(m_nCols);
                    m_pointerE.resize(m_nCols);

                    int skip = 0;

                    for (int ii = 0; ii < size; ++ii)
                    {
                        if (IsReal())
                        {
                            if (rhsMatrix(ii) != static_cast<T1> (0))
                            {
                                m_rows[ii - skip] = ii % m_nRows;
                                ++m_pointerE[ii / m_nRows];
                            }
                            else
                            {
                                ++skip;
                            }
                        }
                        else if (rhsMatrix.z(ii) != static_cast<T1> (0))
                        {
                            m_rows[ii - skip] = ii % m_nRows;
                            ++m_pointerE[ii / m_nRows];
                        }
                        else
                        {
                            ++skip;
                        }
                    }

                    for (int ii = 1; ii < m_nCols; ++ii)
                    {
                        m_pointerB[ii]  = m_pointerE[ii - 1];
                        m_pointerE[ii] += m_pointerE[ii - 1];
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
                m_rows.clear();
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
                bool removeZeros = false;
                const hwTMatrixS<T1, T2>* const_this = this;

                for (int ii = 0; ii < vecSize; ++ii)
                {
                    int index = sliceArg[0].Vector()[ii];

                    if (index < 0 || index >= size)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    if (rhsMatrix.IsReal())
                    {
                        if (rhsMatrix(ii) != static_cast<T1> (0))
                        {
                            if (IsReal())
                                operator()(index) = rhsMatrix(ii);
                            else
                                z(index) = rhsMatrix(ii);
                        }
                        else if (IsReal())
                        {
                            if (const_this->operator()(index) != static_cast<T1> (0))
                            {
                                operator()(index) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                        else if (const_this->z(index) != static_cast<T1> (0))
                        {
                            z(index) = static_cast<T1> (0);
                            removeZeros = true;
                        }
                    }
                    else if (rhsMatrix.z(ii) != static_cast<T1> (0))
                    {
                        z(index) = rhsMatrix.z(ii);
                    }
                    else if (const_this->z(index) != static_cast<T1> (0))
                    {
                        z(index) = static_cast<T1> (0);
                        removeZeros = true;
                    }
                }

                if (removeZeros)
                    RemoveStoredZeros();
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
                        bool removeZeros = false;
                        const hwTMatrixS<T1, T2>* const_this = this;

                        for (int ii = 0; ii < m_nCols; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(row, ii) = rhsMatrix(ii);
                                    else
                                        z(row, ii) = rhsMatrix(ii);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(row, ii) != static_cast<T1> (0))
                                    {
                                        operator()(row, ii) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(row, ii) != static_cast<T1> (0))
                                {
                                    z(row, ii) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(ii) != static_cast<T1> (0))
                            {
                                z(row, ii) = rhsMatrix.z(ii);
                            }
                            else if (const_this->z(row, ii) != static_cast<T1> (0))
                            {
                                z(row, ii) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }

                        if (removeZeros)
                            RemoveStoredZeros();
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

                            for (int jj = col + 1; jj < m_nCols; ++jj)
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
                        bool removeZeros = false;
                        const hwTMatrixS<T1, T2>* const_this = this;

                        for (int colIndex = 0; colIndex < rhsMatrix.Size(); ++colIndex)
                        {
                            int col = sliceArg[1].Vector()[colIndex];

                            if (col < 0 || col >= m_nCols)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(colIndex) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(row, col) = rhsMatrix(colIndex);
                                    else
                                        z(row, col) = rhsMatrix(colIndex);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(row, col) != static_cast<T1> (0))
                                    {
                                        operator()(row, col) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(row, col) != static_cast<T1> (0))
                                {
                                    z(row, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else
                            {
                                if (rhsMatrix.z(colIndex) != static_cast<T1> (0))
                                {
                                    z(row, col) = rhsMatrix.z(colIndex);
                                }
                                else if (const_this->z(row, col) != static_cast<T1> (0))
                                {
                                    z(row, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                        }

                        if (removeZeros)
                            RemoveStoredZeros();
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
                        bool removeZeros = false;
                        const hwTMatrixS<T1, T2>* const_this = this;

                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(ii, col) = rhsMatrix(ii);
                                    else
                                        z(ii, col) = rhsMatrix(ii);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(ii, col) != static_cast<T1> (0))
                                    {
                                        operator()(ii, col) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(ii, col) != static_cast<T1> (0))
                                {
                                    z(ii, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else
                            {
                                if (rhsMatrix.z(ii) != static_cast<T1> (0))
                                {
                                    z(ii, col) = rhsMatrix(ii);
                                }
                                else if (const_this->z(ii, col) != static_cast<T1> (0))
                                {
                                    z(ii, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }

                            if (removeZeros)
                                RemoveStoredZeros();
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
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii, jj) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(ii, jj) = rhsMatrix(ii, jj);
                                    else
                                        z(ii, jj) = rhsMatrix(ii, jj);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(ii, jj) != static_cast<T1> (0))
                                    {
                                        operator()(ii, jj) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(ii, jj) != static_cast<T1> (0))
                                {
                                    z(ii, jj) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(ii, jj) != static_cast<T1> (0))
                            {
                                z(ii, jj) = rhsMatrix.z(ii, jj);
                            }
                            else if (const_this->z(ii, jj) != static_cast<T1> (0))
                            {
                                z(ii, jj) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
                }
                else if (rhsMatrix.Is0x0())
                {
                    m_nRows = 0;
                    m_nCols = 0;
                    m_values.Dimension(0, hwTMatrix<T1, T2>::REAL);
                    m_rows.clear();
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
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < rhsMatrix.N(); ++jj)
                    {
                        int col = sliceArg[1].Vector()[jj];

                        if (col < 0 || col >= m_nCols)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii, jj) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(ii, col) = rhsMatrix(ii, jj);
                                    else
                                        z(ii, col) = rhsMatrix(ii, jj);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(ii, col) != static_cast<T1> (0))
                                    {
                                        operator()(ii, col) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(ii, col) != static_cast<T1> (0))
                                {
                                    z(ii, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(ii, jj) != static_cast<T1> (0))
                            {
                                z(ii, col) = rhsMatrix.z(ii, jj);
                            }
                            else if (const_this->z(ii, col) != static_cast<T1> (0))
                            {
                                z(ii, col) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
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
                        bool removeZeros = false;
                        const hwTMatrixS<T1, T2>* const_this = this;

                        for (int rowIndex = 0; rowIndex < rhsMatrix.Size(); ++rowIndex)
                        {
                            int row = sliceArg[0].Vector()[rowIndex];

                            if (row < 0 || row >= m_nRows)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(rowIndex) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(row, col) = rhsMatrix(rowIndex);
                                    else
                                        z(row, col) = rhsMatrix(rowIndex);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(row, col) != static_cast<T1> (0))
                                    {
                                        operator()(row, col) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(row, col) != static_cast<T1> (0))
                                {
                                    z(row, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(rowIndex) != static_cast<T1> (0))
                            {
                                z(row, col) = rhsMatrix.z(rowIndex);
                            }
                            else if (const_this->z(row, col) != static_cast<T1> (0))
                            {
                                z(row, col) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }

                        if (removeZeros)
                            RemoveStoredZeros();
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
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        for (int ii = 0; ii < rhsMatrix.M(); ++ii)
                        {
                            int row = sliceArg[0].Vector()[ii];

                            if (row < 0 || row >= m_nRows)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (rhsMatrix.IsReal())
                            {
                                if (rhsMatrix(ii, jj) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(row, jj) = rhsMatrix(ii, jj);
                                    else
                                        z(row, jj) = rhsMatrix(ii, jj);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(row, jj) != static_cast<T1> (0))
                                    {
                                        operator()(row, jj) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(row, jj) != static_cast<T1> (0))
                                {
                                    z(row, jj) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(ii, jj) != static_cast<T1> (0))
                            {
                                z(row, jj) = rhsMatrix.z(ii, jj);
                            }
                            else if (const_this->z(row, jj) != static_cast<T1> (0))
                            {
                                z(row, jj) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
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
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

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
                                if (rhsMatrix(ii, jj) != static_cast<T1> (0))
                                {
                                    if (IsReal())
                                        operator()(row, col) = rhsMatrix(ii, jj);
                                    else
                                        z(row, col) = rhsMatrix(ii, jj);
                                }
                                else if (IsReal())
                                {
                                    if (const_this->operator()(row, col) != static_cast<T1> (0))
                                    {
                                        operator()(row, col) = static_cast<T1> (0);
                                        removeZeros = true;
                                    }
                                }
                                else if (const_this->z(row, col) != static_cast<T1> (0))
                                {
                                    z(row, col) = static_cast<T1> (0);
                                    removeZeros = true;
                                }
                            }
                            else if (rhsMatrix.z(ii, jj) != static_cast<T1> (0))
                            {
                                z(row, col) = rhsMatrix.z(ii, jj);
                            }
                            else if (const_this->z(row, col) != static_cast<T1> (0))
                            {
                                z(row, col) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
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
            bool removeZeros = false;
            const hwTMatrixS<T1, T2>* const_this = this;

            for (int ii = 0; ii < vecSize; ++ii)
            {
                int index = static_cast<int> (sliceArg[0].Vector()[ii]);

                if (index < 0 || index >= size)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                if (real != static_cast<T1> (0))
                {
                    operator()(index) = real;
                }
                else if (const_this->operator()(index) != static_cast<T1> (0))
                {
                    operator()(index) = static_cast<T1> (0);
                    removeZeros = true;
                }
            }

            if (removeZeros)
                RemoveStoredZeros();
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
                if (real != static_cast<T1> (0))
                {
                    for (int jj = 0; jj < m_nCols; ++jj)
                        operator()(row, jj) = real;
                }
                else
                {
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        if (const_this->operator()(row, jj) != static_cast<T1> (0))
                        {
                            operator()(row, jj) = static_cast<T1> (0);
                            removeZeros = true;
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
                }
            }
            else if (sliceArg[1].IsVector())
            {
                bool removeZeros = false;
                const hwTMatrixS<T1, T2>* const_this = this;

                for (int colIndex = 0; colIndex < sliceArg[1].Vector().size(); ++colIndex)
                {
                    int col = sliceArg[1].Vector()[colIndex];

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    if (real != static_cast<T1> (0))
                    {
                        operator()(row, col) = real;
                    }
                    else if (const_this->operator()(row, col) != static_cast<T1> (0))
                    {
                        operator()(row, col) = static_cast<T1> (0);
                        removeZeros = true;
                    }
                }

                if (removeZeros)
                    RemoveStoredZeros();
            }
        }
        else if (sliceArg[0].IsColon())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();

                if (col < 0 || col >= m_nCols)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                if (real != static_cast<T1> (0))
                {
                    for (int ii = 0; ii < m_nRows; ++ii)
                        operator()(ii, col) = real;
                }
                else
                {
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int ii = 0; ii < m_nRows; ++ii)
                    {
                        if (const_this->operator()(ii, col) != static_cast<T1> (0))
                        {
                            operator()(ii, col) = static_cast<T1> (0);
                            removeZeros = true;
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
                }
            }
            else if (sliceArg[1].IsColon())
            {
                for (int jj = 0; jj < m_nCols; ++jj)
                {
                    if (real != static_cast<T1> (0))
                    {
                        for (int ii = 0; ii < m_nRows; ++ii)
                            operator()(ii, jj) = real;
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
                if (real != static_cast<T1> (0))
                {
                    for (int jj = 0; jj < sliceArg[1].Vector().size(); ++jj)
                    {
                        int col = sliceArg[1].Vector()[jj];

                        if (col < 0 || col >= m_nCols)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        for (int ii = 0; ii < m_nRows; ++ii)
                            operator()(ii, col) = real;
                    }
                }
                else
                {
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < sliceArg[1].Vector().size(); ++jj)
                    {
                        int col = sliceArg[1].Vector()[jj];

                        if (col < 0 || col >= m_nCols)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (const_this->operator()(ii, col) != static_cast<T1> (0))
                            {
                                operator()(ii, col) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
                }
            }
        }
        else if (sliceArg[0].IsVector())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();
                bool removeZeros = false;
                const hwTMatrixS<T1, T2>* const_this = this;

                for (int rowIndex = 0; rowIndex < sliceArg[0].Vector().size(); ++rowIndex)
                {
                    int row = sliceArg[0].Vector()[rowIndex];

                    if (row < 0 || row >= m_nRows)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    if (real != static_cast<T1> (0))
                    {
                        operator()(row, col) = real;
                    }
                    else if (const_this->operator()(row, col) != static_cast<T1> (0))
                    {
                        operator()(row, col) = static_cast<T1> (0);
                        removeZeros = true;
                    }
                }

                if (removeZeros)
                    RemoveStoredZeros();
            }
            else if (sliceArg[1].IsColon())
            {
                if (real != static_cast<T1> (0))
                {
                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        for (int ii = 0; ii < sliceArg[0].Vector().size(); ++ii)
                        {
                            int row = sliceArg[0].Vector()[ii];

                            if (row < 0 || row >= m_nRows)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            operator()(row, jj) = real;
                        }
                    }
                }
                else
                {
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        for (int ii = 0; ii < sliceArg[0].Vector().size(); ++ii)
                        {
                            int row = sliceArg[0].Vector()[ii];

                            if (row < 0 || row >= m_nRows)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (const_this->operator()(row, jj) != static_cast<T1> (0))
                            {
                                operator()(row, jj) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
                }
            }
            else if (sliceArg[1].IsVector())
            {
                bool removeZeros = false;
                const hwTMatrixS<T1, T2>* const_this = this;

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
                        {
                            operator()(row, col) = real;
                        }
                        else if (const_this->operator()(row, col) != static_cast<T1> (0))
                        {
                            operator()(row, col) = static_cast<T1> (0);
                            removeZeros = true;
                        }
                    }
                    if (removeZeros)
                        RemoveStoredZeros();
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

            if (cmplx != static_cast<T1> (0))
                z(index) = cmplx;
            else
                ZeroElement(index);
        }
        else if (sliceArg[0].IsColon())
        {
            if (cmplx != static_cast<T1> (0))
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
            bool removeZeros = false;
            const hwTMatrixS<T1, T2>* const_this = this;

            for (int ii = 0; ii < vecSize; ++ii)
            {
                int index = static_cast<int> (sliceArg[0].Vector()[ii]);

                if (index < 0 || index >= size)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                if (cmplx != static_cast<T1> (0))
                {
                    z(index) = cmplx;
                }
                else if (const_this->z(index) != static_cast<T1> (0))
                {
                    z(index) = static_cast<T1> (0);
                    removeZeros = true;
                }
            }

            if (removeZeros)
                RemoveStoredZeros();
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

                if (cmplx != static_cast<T1> (0))
                    z(row, col) = cmplx;
                else
                    ZeroElement(row, col);
            }
            else if (sliceArg[1].IsColon())
            {
                if (cmplx != static_cast<T1> (0))
                {
                    for (int jj = 0; jj < m_nCols; ++jj)
                        z(row, jj) = cmplx;
                }
                else
                {
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        if (const_this->z(row, jj) != static_cast<T1> (0))
                        {
                            z(row, jj) = static_cast<T1> (0);
                            removeZeros = true;
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
                }
            }
            else if (sliceArg[1].IsVector())
            {
                bool removeZeros = false;
                const hwTMatrixS<T1, T2>* const_this = this;

                for (int colIndex = 0; colIndex < sliceArg[1].Vector().size(); ++colIndex)
                {
                    int col = sliceArg[1].Vector()[colIndex];

                    if (col < 0 || col >= m_nCols)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    if (cmplx != static_cast<T1> (0))
                    {
                        z(row, col) = cmplx;
                    }
                    else if (const_this->z(row, col) != static_cast<T1> (0))
                    {
                        z(row, col) = static_cast<T1> (0);
                        removeZeros = true;
                    }
                }

                if (removeZeros)
                    RemoveStoredZeros();
            }
        }
        else if (sliceArg[0].IsColon())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();

                if (col < 0 || col >= m_nCols)
                    throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                if (cmplx != static_cast<T1> (0))
                {
                    for (int ii = 0; ii < m_nRows; ++ii)
                        z(ii, col) = cmplx;
                }
                else
                {
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int ii = 0; ii < m_nRows; ++ii)
                    {
                        if (const_this->z(ii, col) != static_cast<T1> (0))
                        {
                            z(ii, col) = static_cast<T1> (0);
                            removeZeros = true;
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
                }
            }
            else if (sliceArg[1].IsColon())
            {
                for (int jj = 0; jj < m_nCols; ++jj)
                {
                    if (cmplx != static_cast<T1> (0))
                    {
                        for (int ii = 0; ii < m_nRows; ++ii)
                            z(ii, jj) = cmplx;
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
                if (cmplx != static_cast<T1> (0))
                {
                    for (int jj = 0; jj < sliceArg[1].Vector().size(); ++jj)
                    {
                        int col = sliceArg[1].Vector()[jj];

                        if (col < 0 || col >= m_nCols)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        for (int ii = 0; ii < m_nRows; ++ii)
                            z(ii, col) = cmplx;
                    }
                }
                else
                {
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < sliceArg[1].Vector().size(); ++jj)
                    {
                        int col = sliceArg[1].Vector()[jj];

                        if (col < 0 || col >= m_nCols)
                            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                        for (int ii = 0; ii < m_nRows; ++ii)
                        {
                            if (const_this->z(ii, col) != static_cast<T1> (0))
                            {
                                z(ii, col) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
                }
            }
        }
        else if (sliceArg[0].IsVector())
        {
            if (sliceArg[1].IsScalar())
            {
                int col = sliceArg[1].Scalar();
                bool removeZeros = false;
                const hwTMatrixS<T1, T2>* const_this = this;

                for (int rowIndex = 0; rowIndex < sliceArg[0].Vector().size(); ++rowIndex)
                {
                    int row = sliceArg[0].Vector()[rowIndex];

                    if (row < 0 || row >= m_nRows)
                        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                    if (cmplx != static_cast<T1> (0))
                    {
                        z(row, col) = cmplx;
                    }
                    else if (const_this->z(row, col) != static_cast<T1> (0))
                    {
                        z(row, col) = static_cast<T1> (0);
                        removeZeros = true;
                    }
                }

                if (removeZeros)
                    RemoveStoredZeros();
            }
            else if (sliceArg[1].IsColon())
            {
                if (cmplx != static_cast<T1> (0))
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
                else
                {
                    bool removeZeros = false;
                    const hwTMatrixS<T1, T2>* const_this = this;

                    for (int jj = 0; jj < m_nCols; ++jj)
                    {
                        for (int ii = 0; ii < sliceArg[0].Vector().size(); ++ii)
                        {
                            int row = sliceArg[0].Vector()[ii];

                            if (row < 0 || row >= m_nRows)
                                throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

                            if (const_this->z(row, jj) != static_cast<T1> (0))
                            {
                                z(row, jj) = static_cast<T1> (0);
                                removeZeros = true;
                            }
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
                }
            }
            else if (sliceArg[1].IsVector())
            {
                bool removeZeros = false;
                const hwTMatrixS<T1, T2>* const_this = this;

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

                        if (cmplx != static_cast<T1> (0))
                        {
                            z(row, col) = cmplx;
                        }
                        else if (const_this->z(row, col) != static_cast<T1> (0))
                        {
                            z(row, col) = static_cast<T1> (0);
                            removeZeros = true;
                        }
                    }

                    if (removeZeros)
                        RemoveStoredZeros();
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

//! Remove stored zeros
template<typename T1, typename T2>
void hwTMatrixS<T1, T2>::RemoveStoredZeros()
{
    if (IsReal())
    {
        hwTMatrixS<T1, T2> copy(m_nRows, m_nCols, pointerB(), pointerE(), rows(),
                                GetRealData());
        MakeEmpty();
        *this = copy;
    }
    else
    {
        hwTMatrixS<T1, T2> copy(m_nRows, m_nCols, pointerB(), pointerE(), rows(),
                                GetComplexData());
        MakeEmpty();
        *this = copy;
    }
}
