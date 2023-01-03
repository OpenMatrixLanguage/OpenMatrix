/**
* @file hwTMatrixN_MKL.cc
* @date August 2022
* Copyright (C) 2022 Altair Engineering, Inc.
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
//  hwTMatrixN explicit <double> specialization implementation file
//  with BLAS / LAPACK or other MKL functions
//
//:---------------------------------------------------------------------------

//*******************************************************************
//                    BLAS / LAPACK prototypes
//*******************************************************************

#include <complex>
typedef std::complex<double> complexD;

// y = x
extern "C" void dcopy_(int* N, double* DX, int* INCX, double* DY, int* INCY);
extern "C" void zcopy_(int* N, complexD* DX, int* INCX, complexD* DY, int* INCY);

//*******************************************************************
//            hwTMatrixN<double> public implementations
//*******************************************************************

//! Copy data
template<>
inline void hwTMatrixN<double>::CopyData(void* dest, const void* src,
                                         int count) const
{
    int inc = 1;

    if (IsReal())
        dcopy_((int*)&count, (double*)src, &inc, (double*)dest, &inc);
    else
        zcopy_((int*)&count, (complexD*)src, &inc, (complexD*)dest, &inc);
}

//! Copy data
template<>
inline void hwTMatrixN<double>::CopyData(void* dest, int stride_dest,
                                         const void* src, int stride_src,
                                         int count) const
{
    if (IsReal())
        dcopy_((int*)&count, (double*)src, &stride_src, (double*)dest, &stride_dest);
    else
        zcopy_((int*)&count, (complexD*)src, &stride_src, (complexD*)dest, &stride_dest);
}

//! Reorder matrix dimensions, a generalized transpose
template<>
inline void hwTMatrixN<double>::Permute(const hwTMatrixN<double>& source, const std::vector<int>& permuteVec)
{
    PermuteCheck(permuteVec);

    // dimension permuted matrix
    MakeEmpty();
    int numDims = static_cast<int> (permuteVec.size());
    int sourceDims = static_cast<int> (source.m_dim.size());
    int maxRHSdim = -1;
    int maxLHSdim = -1;

    if (sourceDims > numDims)
        throw hwMathException(HW_MATH_ERR_PERMVEC1);

    m_dim.resize(numDims);

    for (int i = 0; i < numDims; ++i)
    {
        int k = permuteVec[i];

        if (k < sourceDims)
        {
            m_dim[i] = source.m_dim[k];

            if (maxRHSdim == -1 || source.m_dim[k] > source.m_dim[maxRHSdim])
            {
                maxLHSdim = i;
                maxRHSdim = k;
            }
        }
        else
        {
            m_dim[i] = 1;
        }
    }

    Allocate(source.Type());
    m_bits.realData = source.m_bits.realData;

    // copy data
    m_rhsMatrixIndex.clear();
    m_rhsMatrixIndex.resize(numDims);
    m_lhsMatrixIndex.clear();
    m_lhsMatrixIndex.resize(numDims);
    int strideRHS = source.Stride(maxRHSdim);
    int strideLHS = Stride(maxLHSdim);
    int numVecs = m_size / m_dim[maxLHSdim];

    if (m_real)
    {
        for (int i = 0; i < numVecs; ++i)
        {
            int startRHS = source.Index(m_rhsMatrixIndex);
            int startLHS = Index(m_lhsMatrixIndex);
            const double* realRHS = source.GetRealData() + startRHS;
            double* realLHS = m_real + startLHS;

            CopyData(realLHS, strideLHS, realRHS, strideRHS, m_dim[maxLHSdim]);

            // advance rhs and lhs matrix indices
            for (int j = 0; j < numDims; ++j)
            {
                if (j == maxLHSdim)
                {
                    continue;
                }

                int k = permuteVec[j];

                // increment index k if possible
                if (k < sourceDims)
                {
                    if (m_lhsMatrixIndex[j] < m_dim[j] - 1)
                    {
                        ++m_rhsMatrixIndex[k];
                        ++m_lhsMatrixIndex[j];
                        break;
                    }

                    // index k is maxed out, so reset and continue to k+1
                    m_rhsMatrixIndex[k] = 0;
                    m_lhsMatrixIndex[j] = 0;
                }
            }
        }
    }
    else if (m_complex)
    {
        for (int i = 0; i < numVecs; ++i)
        {
            int startRHS = source.Index(m_rhsMatrixIndex);
            int startLHS = Index(m_lhsMatrixIndex);
            const hwTComplex<double>* realRHS = source.GetComplexData() + startRHS;
            hwTComplex<double>* realLHS = m_complex + startLHS;

            CopyData(realLHS, strideLHS, realRHS, strideRHS, m_dim[maxLHSdim]);

            // advance rhs matrix indices
            for (int j = 0; j < numDims; ++j)
            {
                if (j == maxLHSdim)
                {
                    continue;
                }

                int k = permuteVec[j];

                // increment index k if possible
                if (k < sourceDims)
                {
                    if (m_lhsMatrixIndex[j] < m_dim[j] - 1)
                    {
                        ++m_rhsMatrixIndex[k];
                        ++m_lhsMatrixIndex[j];
                        break;
                    }

                    // index k is maxed out, so reset and continue to k+1
                    m_rhsMatrixIndex[k] = 0;
                    m_lhsMatrixIndex[j] = 0;
                }
            }
        }
    }

    // discard any singleton dimensions that may have been created
    while (m_dim.size() > 2 && m_dim.back() == 1)
    {
        m_dim.pop_back();
    }
}
