/**
* @file hwTMatrixN.cc
* @date June 2014
* Copyright (C) 2014-2018 Altair Engineering, Inc.  
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
//  hwTMatrixN template function implementation file
//
//:---------------------------------------------------------------------------

#include <algorithm>
#include <limits>
#include "../hwMathException.h"
#include "../hwSliceArg.h"
#include "hwTComplex.h"

// ****************************************************
//                  Error handling
// ****************************************************
//  Clients of hwTMatrixN should use try blocks
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
hwTMatrixN<T1, T2>::hwTMatrixN()
    : m_refCount(1)
{
    m_bits.ownData = 0;
    MakeEmpty();
}

//! Construct dimensioned matrix
template<typename T1, typename T2>
hwTMatrixN<T1, T2>::hwTMatrixN(const std::vector<int>& dim, DataType dataType)
    :  m_refCount(1), m_capacity(0)
{
    SetDimensions(dim);
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

//! Construct dimensioned matrix with external data
template<typename T1, typename T2>
hwTMatrixN<T1, T2>::hwTMatrixN(const std::vector<int>& dim, void* data, DataType dataType)
    : m_dim(dim), m_refCount(1), m_capacity(0)
{
    m_bits.ownData = 0;
    m_real_memory = nullptr;
    m_complex_memory = nullptr;

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

    ComputeSize();
}

//! Copy constructor
template<typename T1, typename T2>
hwTMatrixN<T1, T2>::hwTMatrixN(const hwTMatrixN<T1, T2>& source)
    : m_real(nullptr), m_complex(nullptr),
      m_real_memory(nullptr), m_complex_memory(nullptr),
      m_capacity(0), m_refCount(1)
{
    Copy(source);
}

//! Destructor
template<typename T1, typename T2>
hwTMatrixN<T1, T2>::~hwTMatrixN()
{
    if (m_refCount > 1)
        return;

    Deallocate();
}

//! Implement the = operator
template<typename T1, typename T2>
hwTMatrixN<T1, T2>& hwTMatrixN<T1, T2>::operator=(const hwTMatrixN<T1, T2>& rhs)
{
    if (this == &rhs)
        return *this;

    try
    {
        Copy(rhs);
    }
    catch (hwMathException& except)
    {
        MakeEmpty();
        except.Status().ResetArgs();
        throw;
    }

    return *this;
}

// ****************************************************
//                 Private Utilities
// These function are placed high in the file to ensure
// specialization prior to instantiation
// ****************************************************

//! Allocate memory for the matrix data
template<>
inline void hwTMatrixN<double>::Allocate(DataType dataType)
{
    ComputeSize();
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
            throw hwMathException(HW_MATH_ERR_ALLOCFAILED);
        }

        try
        {
            m_real_memory = new char[(m_capacity + 8) * sizeof(double)];
            m_real = reinterpret_cast<double*>
                ((reinterpret_cast<std::ptrdiff_t> (m_real_memory) + 63) & ~0x3F);   // 64-byte alignment
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
            throw hwMathException(HW_MATH_ERR_ALLOCFAILED);
        }

        try
        {
            m_complex_memory = new char[(m_capacity + 4) * sizeof(hwTComplex<double>)];
            m_complex = reinterpret_cast<hwTComplex<double>*>
                ((reinterpret_cast<std::ptrdiff_t> (m_complex_memory) + 63) & ~0x3F);   // 64-byte alignment
        }
        catch (std::bad_alloc&)
        {
            Deallocate();
            m_complex_memory = nullptr;
            throw;
        }
    }
}
/*
// This would be nice to have, but needs to be in a separate .cc file
// to avoid conflict when using the MKL based version.
#include <memory.h>

//! Copy data
template<>
inline void hwTMatrixN<double>::CopyData(void* dest, const void* src, int count) const
{
    if (IsReal())
        memcpy_s(dest, count * sizeof(double), src, count * sizeof(double));
    else
        memcpy_s(dest, count * sizeof(hwTComplex<double>), src, count * sizeof(hwTComplex<double>));
}
*/
// ****************************************************
//               Data Type, Ownership
// ****************************************************

//! Determine if the matrix contains real or complex data
//! (i.e. no non-zero imaginary components)
template<typename T1, typename T2>
bool hwTMatrixN<T1, T2>::IsRealData() const
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

// ****************************************************
//           hwTMatrix conversions functions
// ****************************************************
//! Convert hwTMatrixN to hwTMatrix
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::ConvertNDto2D(hwTMatrix<T1, T2>& target, bool copyData) const
{
    // This method repackages 1D or 2D data from an ND matrix into a 2D matrix
    int m = m_dim[0];
    int n = m_dim[1];

    for (int i = 2; i < m_dim.size(); ++i)
    {
        if (m_dim[i] > 2)
        {
            if (m == 1)
            {
                m = n;
                n = m_dim[i];
            }
            else if (n == 1)
            {
                n = m_dim[i];
            }
            else
            {
                throw hwMathException(HW_MATH_ERR_NOTALLOWED, 1);
            }
        }
    }

    if (copyData)
    {
        if (IsReal())
            target.Dimension(m, n, hwTMatrix<T1, T2>::REAL);
        else
            target.Dimension(m, n, hwTMatrix<T1, T2>::COMPLEX);

        if (m_size > 0)
        {
            if (m_real)
                CopyData(target.GetRealData(), m_real, m_size);
            else if (m_complex)
                CopyData(target.GetComplexData(), m_complex, m_size);
        }
    }
    else    // borrow data without taking ownership
    {
        if (IsReal())
        {
            target.Borrow(m, n, m_real_memory, m_real, hwTMatrix<T1, T2>::REAL);
        }
        else
        {
            target.Borrow(m, n, m_complex_memory, m_complex, hwTMatrix<T1, T2>::COMPLEX);
        }
    }
}

//! Convert hwTMatrixN to hwTMatrix
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::ConvertNDto2D(hwTMatrix<T1, T2>& target, bool copyData)
{
    // This method repackages 1D or 2D data from an ND matrix into a 2D matrix
    int m = m_dim[0];
    int n = m_dim[1];

    for (int i = 2; i < m_dim.size(); ++i)
    {
        if (m_dim[i] > 2)
        {
            if (m == 1)
            {
                m = n;
                n = m_dim[i];
            }
            else if (n == 1)
            {
                n = m_dim[i];
            }
            else
            {
                throw hwMathException(HW_MATH_ERR_NOTALLOWED, 1);
            }
        }
    }

    if (copyData)
    {
        if (IsReal())
            target.Dimension(m, n, hwTMatrix<T1, T2>::REAL);
        else
            target.Dimension(m, n, hwTMatrix<T1, T2>::COMPLEX);

        if (m_size > 0)
        {
            if (m_real)
                CopyData(target.GetRealData(), m_real, m_size);
            else if (m_complex)
                CopyData(target.GetComplexData(), m_complex, m_size);
        }
    }
    else    // borrow data, and transfer ownership (if applicable)
    {
        if (IsReal())
        {
            target.Borrow(m, n, m_real_memory, m_real, hwTMatrix<T1, T2>::REAL);
        }
        else
        {
            target.Borrow(m, n, m_complex_memory, m_complex, hwTMatrix<T1, T2>::COMPLEX);
        }

        if (m_bits.ownData)
        {
            OwnData(false);
            target.OwnData(true);
        }
        else
        {
            target.OwnData(false);
        }
    }
}

//! Convert hwTMatrix to hwTMatrixN
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Convert2DtoND(const hwTMatrix<T1, T2>& source, bool copyData)
{
    // This method repackages 1D or 2D data from a 2D matrix into a ND matrix
    if (!IsEmpty())
        MakeEmpty();

    m_dim[0] = source.M();
    m_dim[1] = source.N();

    if (copyData)
    {
        if (source.IsReal())
            Dimension(m_dim, REAL);
        else
            Dimension(m_dim, COMPLEX);

        if (m_size > 0)
        {
            if (m_real)
                CopyData(m_real, source.GetRealData(), m_size);
            else if (m_complex)
                CopyData(m_complex, source.GetComplexData(), m_size);
        }
    }
    else    // borrow data without taking ownership
    {
        OwnData(false);
        ComputeSize();

        if (source.IsReal())
        {
            m_real = const_cast<T1*> (source.GetRealData());
            m_bits.realData = 1;
        }
        else
        {
            m_complex = const_cast<T2*> (source.GetComplexData());
            m_bits.realData = 0;
        }
    }
}

//! Convert hwTMatrix to hwTMatrixN
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Convert2DtoND(hwTMatrix<T1, T2>& source, bool copyData)
{
    // This method repackages 1D or 2D data from a 2D matrix into a ND matrix
    if (!IsEmpty())
        MakeEmpty();

    m_dim[0] = source.M();
    m_dim[1] = source.N();

    if (copyData)
    {
        if (source.IsReal())
            Dimension(m_dim, REAL);
        else
            Dimension(m_dim, COMPLEX);

        if (m_size > 0)
        {
            if (m_real)
                CopyData(m_real, source.GetRealData(), m_size);
            else if (m_complex)
                CopyData(m_complex, source.GetComplexData(), m_size);
        }
    }
    else    // borrow data, and transfer ownership (if applicable)
    {
        ComputeSize();

        if (source.IsReal())
        {
            m_real_memory = source.GetRealMemory();
            m_real = source.GetRealData();
            m_bits.realData = 1;
        }
        else
        {
            m_complex_memory = source.GetComplexMemory();
            m_complex = source.GetComplexData();
            m_bits.realData = 0;
        }

        if (source.OwnData())
        {
            source.OwnData(false);
            m_bits.ownData = 1;
        }
        else
        {
            m_bits.ownData = 0;
        }
    }
}

//! Convert hwTMatrixN to hwTMatrix
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::ConvertNDto1D(hwTMatrix<T1, T2>& target) const
{
    // This method repackages ND data into a vector,
    // borrowing data without taking ownership
    if (IsReal())
    {
        target.Borrow(Size(), 1, m_real_memory, m_real, hwTMatrix<T1, T2>::REAL);
    }
    else
    {
        target.Borrow(Size(), 1, m_complex_memory, m_complex, hwTMatrix<T1, T2>::COMPLEX);
    }
}

//! Convert hwTMatrix to hwTMatrixN
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Convert1DtoND(hwTMatrix<T1, T2>& source, const std::vector<int>& dim)
{
    // This method repackages a vector into an ND matrix,
    // borrowring data, and transferring ownership (if applicable)
    if (!IsEmpty())
        MakeEmpty();

    m_dim = dim;
    ComputeSize();

    if (source.IsReal())
    {
        m_real_memory = source.GetRealMemory();
        m_real = source.GetRealData();
        m_bits.realData = 1;
    }
    else
    {
        m_complex_memory = source.GetComplexMemory();
        m_complex = source.GetComplexData();
        m_bits.realData = 0;
    }

    if (source.OwnData())
    {
        source.OwnData(false);
        m_bits.ownData = 1;
    }
    else
    {
        m_bits.ownData = 0;
    }
}

// ****************************************************
//                Set Matrix Dimensions
// ****************************************************
//! Dimension matrix, specifying real or complex
//! Applies to empty or populated matrix, does not preserve existing data
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Dimension(const std::vector<int>& dim, DataType dataType)
{
    if ((m_real && dataType == REAL) ||
        (m_complex && dataType == COMPLEX))
    {
        // nothing to do if dimensions have not changed
        int effDimSize = _min((int) m_dim.size(), (int) dim.size());
        bool sizeChange = false;

        for (int i = 0; i < effDimSize; ++i)
        {
            if (m_dim[i] != dim[i])
            {
                sizeChange = true;
                break;
            }
        }

        if (!sizeChange)
        {
            for (int i = effDimSize; i < m_dim.size(); ++i)
            {
                if (m_dim[i] != 1)
                {
                    sizeChange = true;
                    break;
                }
            }

            for (int i = effDimSize; i < dim.size(); ++i)
            {
                if (dim[i] != 1)
                {
                    sizeChange = true;
                    break;
                }
            }
        }

        if (!sizeChange)
            return;
    }

    if (m_real || m_complex)
    {
        MakeEmpty();
    }

    SetDimensions(dim);
    Allocate(dataType);

    if (!m_complex)
        m_bits.realData = 1;
    else
        m_bits.realData = 0;
}

//! Ignore high dimension singleton indices
template<typename T1, typename T2>
int hwTMatrixN<T1, T2>::RelevantNumberOfSlices(const std::vector<hwSliceArg>& sliceArg,
                                               bool matrixAssignment) const
{
    int relevantNum = (int) sliceArg.size();

    for (int i = relevantNum-1; i > m_dim.size()-1; --i)
    {
        if (sliceArg[i].IsScalar())
        {
            if (sliceArg[i].Scalar() == 0)
                --relevantNum;
            else
                break;
        }
        else if (sliceArg[i].IsColon())
        {
            if (matrixAssignment)   // a(slice) = matrix
                break;
            else
                --relevantNum;
        }
        else if (sliceArg[i].IsVector())
        {
            if (sliceArg[i].Vector().size() == 1)
            {
                if (sliceArg[i].Vector()[0] == 0)
                    --relevantNum;
                else
                    break;
            }
            else
            {
                break;
            }
        }
    }

    return relevantNum;
}

//! Resize the matrix, leaving the data type unchanged
//! Applies to populated matrix, can only resize an empty matrix if it remains empty
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Resize(const hwTMatrixN<T1, T2>& source, const std::vector<int>& dim, bool initZero)
{
    const std::vector<int>& sourceDim = source.Dimensions();
    int numDims = static_cast<int> (sourceDim.size());

    if (dim.size() != numDims)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    for (int i = 0; i < numDims; ++i)
    {
        if (dim[i] < 0)
            throw hwMathException(HW_MATH_ERR_ARRAYDIM, 1);
    }

    Dimension(dim, source.Type());

    if (IsEmpty())
        return;

    int numVecs = Size() / dim[0];
    int inc = 1;
    m_lhsMatrixIndex.resize(numDims);
    int maxZVdim = 0;
    bool zeroVector = false;

    if (IsReal())
    {
        T1* src = const_cast<T1*> (source.GetRealData());
        T1* dst = m_real;

        for (int i = 0; i < numVecs; ++i)
        {
            if (!zeroVector)
            {
                int srcIndex = source.Index(m_lhsMatrixIndex);

                if (sourceDim[0] < dim[0])
                {
                    // NOTE: If Resize is needed with non-trivial types then use a loop for
                    //       assignments and switch to dcopy in a partial template instantiation.
                    //       See CopyData for a similar case.
                    // dcopy_((int*) &sourceDim[0], src + srcIndex, &inc, dst, &inc);
                    memcpy_s(dst, sourceDim[0]*sizeof(T1), src + srcIndex, sourceDim[0]*sizeof(T1));

                    if (initZero)
                        memset(dst + sourceDim[0], 0, (dim[0] - sourceDim[0]) * sizeof(T1));
                }
                else
                {
                    // dcopy_((int*) &dim[0], src + srcIndex, &inc, dst, &inc);
                    memcpy_s(dst, dim[0]*sizeof(T1), src + srcIndex, dim[0]*sizeof(T1));
                }
            }
            else if (initZero)
            {
                memset(dst, 0, dim[0] * sizeof(T1));
            }

            // advance matrix pointers and indices
            dst += dim[0];

            for (int j = 1; j < numDims; ++j)
            {
                // increment index j if possible
                if (m_lhsMatrixIndex[j] < dim[j] - 1)
                {
                    ++m_lhsMatrixIndex[j];

                    if (m_lhsMatrixIndex[j] == sourceDim[j])
                    {
                        zeroVector = true;
                        maxZVdim   = _max(maxZVdim, j);
                    }

                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_lhsMatrixIndex[j] = 0;

                if (j == maxZVdim)
                {
                    zeroVector = false;
                    maxZVdim   = 0;
                }
            }
        }
    }
    else    // complex
    {
        T2* src = const_cast<T2*> (source.GetComplexData());
        T2* dst = m_complex;

        for (int i = 0; i < numVecs; ++i)
        {
            if (!zeroVector)
            {
                int srcIndex = source.Index(m_lhsMatrixIndex);

                if (sourceDim[0] < dim[0])
                {
                    // zcopy_((int*) &sourceDim[0], (complexD*) src + srcIndex, &inc, (complexD*) dst, &inc);
                    memcpy_s(dst, sourceDim[0]*sizeof(T2), src + srcIndex, sourceDim[0]*sizeof(T2));

                    if (initZero)
                        memset(dst + sourceDim[0], 0, (dim[0] - sourceDim[0]) * sizeof(T2));
                }
                else
                {
                    // zcopy_((int*) &dim[0], (complexD*) src + srcIndex, &inc, (complexD*) dst, &inc);
                    memcpy_s(dst, dim[0]*sizeof(T2), src + srcIndex, dim[0]*sizeof(T2));
                }
            }
            else if (initZero)
            {
                memset(dst, 0, dim[0] * sizeof(T2));
            }

            // advance matrix pointers and indices
            dst += dim[0];

            for (int j = 1; j < numDims; ++j)
            {
                // increment index j if possible
                if (m_lhsMatrixIndex[j] < dim[j] - 1)
                {
                    ++m_lhsMatrixIndex[j];

                    if (m_lhsMatrixIndex[j] == sourceDim[j])
                    {
                        zeroVector = true;
                        maxZVdim = _max(maxZVdim, j);
                    }

                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_lhsMatrixIndex[j] = 0;

                if (j == maxZVdim)
                {
                    zeroVector = false;
                    maxZVdim = 0;
                }
            }
        }
    }
}

//! Change the dimensions of a matrix while maintaining the same number of elements
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Reshape(const std::vector<int>& dim)
{
    SetDimensions(dim);

    // check new size
    int size = 1;
    int emptyDim = -1;
    
    for (int i = 0; i < m_dim.size(); ++i)
    {
        // a single element of dim can be set to -1 to allow the function
        // to set it as needed if an acceptable size exists
        if (m_dim[i] == -1)
        {
            if (emptyDim != -1)
                throw hwMathException(HW_MATH_ERR_MATRIXRESHAPE1, 1);

            emptyDim = i;
            continue;
        }
        else if (m_dim[i] < 0)
            throw hwMathException(HW_MATH_ERR_ARRAYDIM, 1);

        size *= m_dim[i];
    }

    if (emptyDim != -1)
    {
        if (size)
        {
            m_dim[emptyDim] = m_size / size;
            size *= m_dim[emptyDim];
        }
        else
        {
            m_dim[emptyDim] = 0;
        }
    }

    if (size != m_size)
        throw hwMathException(HW_MATH_ERR_MATRIXRESHAPE2, 1);
}

//! Replicate an ND matrix in each dimension
template<>
inline void hwTMatrixN<double>::Repmat(const hwTMatrixN<double>& source, const std::vector<int>& reps)
{
    // function needs to be optimized
    bool justCopy = true;

    for (int i = 0; i < reps.size(); ++i)
    {
        if (reps[i] != 1)
        {
            justCopy = false;
            break;
        }
    }

    if (justCopy)
    {
        (*this) = source;
    }

    std::vector<int> dimR;
    const hwTMatrixN<double>& base = source;
    hwTMatrixN<double>::DataType type = base.Type();

    const std::vector<int>& dimB = base.Dimensions();
    int min = _min(static_cast<int> (reps.size()), static_cast<int> (dimB.size()));

    // get output matrix dimensions and type
    for (int i = 0; i < min; ++i)
        dimR.push_back(dimB[i] * reps[i]);

    for (int i = min; i < dimB.size(); ++i)
        dimR.push_back(dimB[i]);

    for (int i = min; i < reps.size(); ++i)
        dimR.push_back(reps[i]);

    //// populate result
    hwTMatrixN<double>* result = this;
    result->Dimension(dimR, type);

    if (base.Size() == 1)
    {
        if (base.IsReal())
            result->SetElements(base(0));
        else
            result->SetElements(base.z(0));
    }
    else
    {
        const std::vector<int>& dimB = base.Dimensions();
        std::vector<int> base_index(dimB.size());
        int numElemsToCopy = dimB[0];
        int numBlocks = base.Size() / numElemsToCopy;

        if (result->IsReal() && !result->IsEmpty())
        {
            // copy base to result
            const double* src = base.GetRealData();
            double* dest = result->GetRealData();

            for (int j = 0; j < numBlocks; ++j)
            {
                if (numElemsToCopy == 1)
                {
                    (*dest) = (*src);
                }
                else
                {
                    memcpy_s(dest, numElemsToCopy * sizeof(double), src, numElemsToCopy * sizeof(double));
                }

                // advance base indices
                for (int k = 1; k < dimB.size(); ++k)
                {
                    // increment index k if possible
                    if (base_index[k] != dimB[k] - 1)
                    {
                        ++base_index[k];
                        break;
                    }

                    // index k is maxed out, so reset and continue to k+1
                    base_index[k] = 0;
                }

                src = base.GetRealData() + base.Index(base_index);
                dest = result->GetRealData() + result->Index(base_index);
            }

            // replicate the copied base in all dimensions
            numBlocks = base.Size();
            numElemsToCopy = 1;
            std::vector<int> result_index(dimR.size());

            for (int i = 0; i < reps.size(); ++i)
            {
                if (i < dimB.size())
                {
                    numBlocks /= dimB[i];
                    numElemsToCopy *= dimB[i];
                }

                double* src = result->GetRealData();

                for (int j = 0; j < numBlocks; ++j)
                {
                    double* dest = src;

                    if (numElemsToCopy == 1)
                    {
                        for (int k = 1; k < reps[i]; ++k)
                        {
                            ++dest;
                            (*dest) = (*src);
                        }
                    }
                    else
                    {
                        for (int k = 1; k < reps[i]; ++k)
                        {
                            dest += numElemsToCopy;
                            memcpy_s(dest, numElemsToCopy * sizeof(double), src, numElemsToCopy * sizeof(double));
                            src = dest;
                        }
                    }

                    // advance result indices
                    for (int k = i + 1; k < dimB.size(); ++k)
                    {
                        // increment index k if possible
                        if (result_index[k] < dimB[k] - 1)
                        {
                            ++result_index[k];
                            break;
                        }

                        // index k is maxed out, so reset and continue to k+1
                        result_index[k] = 0;
                    }

                    src = result->GetRealData() + result->Index(result_index);
                }

                numElemsToCopy *= reps[i];
            }
        }
        else if (!result->IsEmpty())  // complex
        {
            // copy base to result
            const hwTComplex<double>* src = base.GetComplexData();
            hwTComplex<double>* dest = result->GetComplexData();

            for (int j = 0; j < numBlocks; ++j)
            {
                if (numElemsToCopy == 1)
                {
                    (*dest) = (*src);
                }
                else
                {
                    memcpy_s(dest, numElemsToCopy * sizeof(hwTComplex<double>), src, numElemsToCopy * sizeof(hwTComplex<double>));
                }

                // advance base indices
                for (int k = 1; k < dimB.size(); ++k)
                {
                    // increment index k if possible
                    if (base_index[k] != dimB[k] - 1)
                    {
                        ++base_index[k];
                        break;
                    }

                    // index k is maxed out, so reset and continue to k+1
                    base_index[k] = 0;
                }

                src = base.GetComplexData() + base.Index(base_index);
                dest = result->GetComplexData() + result->Index(base_index);
            }

            // replicate the copied base in all dimensions
            numBlocks = base.Size();
            numElemsToCopy = 1;
            std::vector<int> result_index(dimR.size());

            for (int i = 0; i < reps.size(); ++i)
            {
                if (i < dimB.size())
                {
                    numBlocks /= dimB[i];
                    numElemsToCopy *= dimB[i];
                }

                hwTComplex<double>* src = result->GetComplexData();

                for (int j = 0; j < numBlocks; ++j)
                {
                    hwTComplex<double>* dest = src;

                    if (numElemsToCopy == 1)
                    {
                        for (int k = 1; k < reps[i]; ++k)
                        {
                            ++dest;
                            (*dest) = (*src);
                        }
                    }
                    else
                    {
                        for (int k = 1; k < reps[i]; ++k)
                        {
                            dest += numElemsToCopy;
                            memcpy_s(dest, numElemsToCopy * sizeof(hwTComplex<double>), src, numElemsToCopy * sizeof(hwTComplex<double>));
                            src = dest;
                        }
                    }

                    // advance result indices
                    for (int k = i + 1; k < dimB.size(); ++k)
                    {
                        // increment index k if possible
                        if (result_index[k] < dimB[k] - 1)
                        {
                            ++result_index[k];
                            break;
                        }

                        // index k is maxed out, so reset and continue to k+1
                        result_index[k] = 0;
                    }

                    src = result->GetComplexData() + result->Index(result_index);
                }

                numElemsToCopy *= reps[i];
            }
        }
    }
}

//! Reorder matrix dimensions, a generalized transpose
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Repmat(const hwTMatrixN<T1, T2>& source, const std::vector<int>& reps)
{
    throw hwMathException(HW_MATH_ERR_NOTIMPLEMENT);
}

//! Replicate a 2D matrix in each dimension
template<>
inline void hwTMatrixN<double>::Repmat(const hwTMatrix<double>& source, const std::vector<int>& reps)
{
    hwTMatrixN<double> base;
    base.Convert2DtoND(source, false);
    Repmat(base, reps);
}

//! Reorder matrix dimensions, a generalized transpose
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Repmat(const hwTMatrix<T1, T2>& source, const std::vector<int>& reps)
{
    throw hwMathException(HW_MATH_ERR_NOTIMPLEMENT);
}

//! Reorder matrix dimensions, a generalized transpose
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Permute(const hwTMatrixN<T1, T2>& source, const std::vector<int>& permuteVec)
{
    PermuteCheck(permuteVec);

    // dimension permuted matrix
    MakeEmpty();
    int numDims = static_cast<int> (permuteVec.size());
    int sourceDims = static_cast<int> (source.m_dim.size());

    if (sourceDims > numDims)
        throw hwMathException(HW_MATH_ERR_PERMVEC1);

    m_dim.resize(numDims);

    for (int i = 0; i < numDims; ++i)
    {
        if (permuteVec[i] < sourceDims)
            m_dim[i] = source.m_dim[permuteVec[i]];
        else
            m_dim[i] = 1;
    }

    Allocate(source.Type());
	m_bits.realData = source.m_bits.realData;

    // copy data
    m_rhsMatrixIndex.clear();
    m_rhsMatrixIndex.resize(numDims);

    if (m_real)
    {
        if (m_size)
            m_real[0] = source.m_real[0];

        for (int i = 1; i < m_size; ++i)
        {
            // advance rhs matrix indices
            for (int j = 0; j < numDims; ++j)
            {
                int k = permuteVec[j];

                // increment index k if possible
                if (k < sourceDims)
                {
                    if (m_rhsMatrixIndex[k] < source.m_dim[k] - 1)
                    {
                        ++m_rhsMatrixIndex[k];
                        break;
                    }

                    // index k is maxed out, so reset and continue to k+1
                    m_rhsMatrixIndex[k] = 0;
                }
            }

            m_real[i] = source(m_rhsMatrixIndex);
        }
    }
    else if (m_complex)
    {
        if (m_size)
            m_complex[0] = source.m_complex[0];

        for (int i = 1; i < m_size; ++i)
        {
            // advance rhs matrix indices
            for (int j = 0; j < numDims; ++j)
            {
                int k = permuteVec[j];

                // increment index k if possible
                if (k < sourceDims)
                {
                    if (m_rhsMatrixIndex[k] < source.m_dim[k] - 1)
                    {
                        ++m_rhsMatrixIndex[k];
                        break;
                    }

                    // index k is maxed out, so reset and continue to k+1
                    m_rhsMatrixIndex[k] = 0;
                }
            }

            m_complex[i] = source.z(m_rhsMatrixIndex);
        }
    }

    // discard any singleton dimensions that may have been created
    while (m_dim.size() > 2 && m_dim.back() == 1)
    {
        m_dim.pop_back();
    }
}

//! Verify permutation vector validity
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::PermuteCheck(const std::vector<int>& permuteVec)
{
    int numDims = static_cast<int> (permuteVec.size());
    bool* hasDim = new bool[numDims];

    for (int i = 0; i < numDims; ++i)
        hasDim[i] = false;

    for (int i = 0; i < numDims; ++i)
    {
        int dim = permuteVec[i];

        if (dim < 0)
            throw hwMathException(HW_MATH_ERR_PERMVEC2);

        if (dim >= numDims)
            throw hwMathException(HW_MATH_ERR_PERMVEC3);

        if (hasDim[dim])
            throw hwMathException(HW_MATH_ERR_PERMVEC4);

        hasDim[dim] = true;
    }

    delete [] hasDim;
}

// ****************************************************
//                  Matrix Properties
// ****************************************************

//! Determine if the matrix is a vector
template<typename T1, typename T2>
bool hwTMatrixN<T1, T2>::IsVector() const
{
    if (!m_real && !m_complex)
        return false;

    int numDim     = static_cast<int> (m_dim.size());
    int numNonOnes = 0;

    for (int i = 0; i < numDim; ++i)
    {
        if (m_dim[i] != 1)
        {
            if (++numNonOnes == 2)
                return false;
        }
    }

    return true;
}

//! Determine if the matrix is empty or a vector
template<typename T1, typename T2>
bool hwTMatrixN<T1, T2>::IsEmptyOrVector() const
{
    if (IsEmpty())
        return true;

    return IsVector();
}

//! Bounds check - upper and lower
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::BoundsCheckRHS(const std::vector<int>& indexVec) const
{
    int min = _min((int) m_dim.size(), (int) indexVec.size());

    for (int i = 0; i < min; ++i)
    {
        if (indexVec[i] < 0)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);

        if (indexVec[i] > m_dim[i]-1)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);
    }

    for (int i = min; i < indexVec.size(); ++i)
    {
        if (indexVec[i] != 0)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);
    }
}

//! Bounds check - upper and lower
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::BoundsCheckRHS(const std::vector<hwSliceArg>& sliceArg, int numSlices) const
{
    if (numSlices == 0)
        numSlices = RelevantNumberOfSlices(sliceArg, true);

    int common = _min((int) m_dim.size(), numSlices);

    for (int i = 0; i < common; ++i)
    {
        if (sliceArg[i].IsScalar())
        {
            if (sliceArg[i].Scalar() > m_dim[i] - 1)
                throw hwMathException(HW_MATH_ERR_SLICE_INDEX, 1);
        }
        else if (sliceArg[i].IsColon())
        {
        }
        else if (sliceArg[i].IsVector())
        {
            for (int j = 0; j < sliceArg[i].Vector().size(); ++j)
            {
                if (sliceArg[i].Vector()[j] < 0)
                    throw hwMathException(HW_MATH_ERR_SLICE_INDEX, 1);

                if (sliceArg[i].Vector()[j] > m_dim[i] - 1)
                    throw hwMathException(HW_MATH_ERR_SLICE_INDEX, 1);
            }
        }
    }

    for (int i = common; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar())
        {
            if (sliceArg[i].Scalar() != 0)
                throw hwMathException(HW_MATH_ERR_SLICE_INDEX, 1);
        }
        else if (sliceArg[i].IsColon())
        {
        }
        else if (sliceArg[i].IsVector())
        {
            for (int j = 0; j < sliceArg[i].Vector().size(); ++j)
            {
                if (sliceArg[i].Vector()[j] != 0)
                    throw hwMathException(HW_MATH_ERR_SLICE_INDEX, 1);
            }
        }
    }
}

//! Grow a matrix
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::GrowLHSMatrix(const std::vector<int>& newDim)
{
    // create matrix with expanded dimensions
    int numDims = (int) m_dim.size();

    if (newDim.size() < numDims)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    std::vector<int> oldDim = m_dim;
    m_dim = newDim;
    ComputeSize();

    if (m_real)
    {
        char* m_real_mem_old = m_real_memory;
        T1* m_real_old = m_real;
        hwTMatrixN<T1, T2> tempMatrix(oldDim, m_real_old, REAL);
        m_real_memory = nullptr;
        m_real = nullptr;
        Allocate(REAL);

        // copy zeros to the expanded matrix
        SetElements(0.0);

        // copy old data to the expanded matrix
        CopyMatrixLHS(tempMatrix);

        if (m_bits.ownData)
        {
            if (m_real_mem_old)
                delete[] m_real_mem_old;
            else if (m_real_old)
                delete[] m_real_old;
        }
        else
        {
            m_bits.ownData = 1;
        }
    }
    else if (m_complex)
    {
        char* m_complex_mem_old = m_complex_memory;
        T2* m_complex_old = m_complex;
        hwTMatrixN<T1, T2> tempMatrix(oldDim, m_complex_old, COMPLEX);
        m_complex_memory = nullptr;
        m_complex = nullptr;
        Allocate(COMPLEX);

        // copy zeros to the expanded matrix
        SetElements(0.0);

        // copy old data to the expanded matrix
        CopyMatrixLHS(tempMatrix);

        if (m_bits.ownData)
        {
            if (m_complex_mem_old)
                delete[] m_complex_mem_old;
            else if (m_complex_old)
                delete[] m_complex_old;
        }
        else
        {
            m_bits.ownData = 1;
        }
    }
    else
    {
        // cannot allocate an empty matrix with GrowLHSMatrix
        // verify that at least one indexVec element is 0
        std::vector<int>::iterator it;

        it = find(m_dim.begin(), m_dim.end(), 0);

        if (it != m_dim.end())
    	    throw hwMathException(HW_MATH_ERR_NOTALLOWED);
    }
}

// ****************************************************
//               Indexing Functions
// ****************************************************

//! Return the index vector corresponding to a single index
template<typename T1, typename T2>
std::vector<int> hwTMatrixN<T1, T2>::IndexVector(int index) const
{
    if (index < 0 || index >= m_size)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    std::vector<int> indexVec(m_dim.size());

    for (int i = 0; i < indexVec.size(); ++i)
    {
        int numItems = index / m_dim[i];            // number of complete items in the dimension
        indexVec[i] = index - numItems * m_dim[i];  // number of remaining items in the dimension
        index = numItems;
    }

    return indexVec;
}

//! Return the single index corresponding to an index vector
template<typename T1, typename T2>
int hwTMatrixN<T1, T2>::Index(const std::vector<int>& indexVec) const
{
    int numDim = static_cast<int> (indexVec.size());

    for (int i = numDim-1; i > m_dim.size()-1; --i)
    {
        if (indexVec[i] != 0)
            throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
        else
            --numDim;
    }

    int pos = indexVec[numDim-1];

    for (int i = numDim-2; i > -1; --i)
        pos = pos * m_dim[i] + indexVec[i];

    return pos;
}

//! Return the stride in a dimension of interest
template<typename T1, typename T2>
int hwTMatrixN<T1, T2>::Stride(int dim) const
{
    int numDims = static_cast<int> (m_dim.size());
    int stride = 1;

    if (dim < 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

    dim = _min(dim, numDims);

    for (int i = 1; i <= dim; ++i)
    {
        stride *= m_dim[i - 1];
    }

    return stride;
}

// ****************************************************
//         Access Functions for Real Elements
// ****************************************************

//! Return a reference to the real data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
T1& hwTMatrixN<T1, T2>::operator()(const std::vector<int>& indexVec)
{
    return m_real[Index(indexVec)];
}

//! Return a const reference to the real data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
const T1& hwTMatrixN<T1, T2>::operator()(const std::vector<int>& indexVec) const
{
    return m_real[Index(indexVec)];
}

//! Set every matrix element to a specified real value
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::SetElements(T1 real)
{
    int size = Size();

    if (m_real)
    {
        for (int i = 0; i < size; ++i)
            m_real[i] = real;
    }

    // the function does not reallocate a complex matrix as real
    // the user must manage that if desired
    if (m_complex)
    {
        for (int i = 0; i < size; ++i)
            m_complex[i] = real;
    }
}

// ****************************************************
//         Access Functions for Complex Elements
// ****************************************************

//! Return a reference to the complex data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
T2& hwTMatrixN<T1, T2>::z(const std::vector<int>& indexVec)
{
    return m_complex[Index(indexVec)];
}

//! Return a const reference to the complex data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
const T2& hwTMatrixN<T1, T2>::z(const std::vector<int>& indexVec) const
{
    return m_complex[Index(indexVec)];
}

//! Set every matrix element to a specified real value
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::SetElements(const T2& cmplx)
{
    if (cmplx.Imag() == (T1) 0)
    {
        SetElements(cmplx.Real());
        return;
    }

    int size = Size();

    if (m_real)
        MakeComplex();

    if (m_complex)
    {
        for (int i = 0; i < size; ++i)
            m_complex[i] = cmplx;
    }
}

// ****************************************************
//             Real / Complex Conversions
// ****************************************************

//! Convert the matrix to complex
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::MakeComplex()
{
    if (m_complex || !m_real)      // already complex, or empty
        return;

    // Note: an empty matrix is real by convention, but there is
    // no reason to trigger an error. The client will still have
    // to dimension the matrix somewhere.

    Allocate(COMPLEX);

    int size = Size();

    for (int i = 0; i < size; ++i)
    {
        m_complex[i].Real() = m_real[i];
        m_complex[i].Imag() = (T1) 0;
    }

    if (m_real_memory)
    {
        delete [] m_real_memory;
        m_real_memory = nullptr;
    }
    else if (m_real)
    {
        delete[] m_real;
    }

    m_real = nullptr;
    m_bits.realData = 0;
    return;
}

//! Pack real and imaginary components into a complex matrix
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::PackComplex(const hwTMatrixN<T1, T2>& real, const hwTMatrixN<T1, T2>* imag)
{
    if (!real.IsReal())
    {
        throw hwMathException(HW_MATH_ERR_ARRAYTYPE, 1);
    }

    if (imag)
    {
        if (real.m_dim != imag->m_dim)
            throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        if (!imag->IsReal())
            throw hwMathException(HW_MATH_ERR_ARRAYTYPE, 2);
    }

    Dimension(real.m_dim, COMPLEX);

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
}

//! Unpack real and imaginary components from a complex matrix
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::UnpackComplex(hwTMatrixN<T1, T2>* real, hwTMatrixN<T1, T2>* imag) const
{
    if (IsReal())
    {
        throw hwMathException(HW_MATH_ERR_ARRAYTYPE, 0);
    }

    if (real)
    {
        try
        {
            real->Dimension(m_dim, REAL);
        }
        catch (hwMathException& except)
        {
            except.Status().SetArg1(1);
            throw;
        }
    }

    if (imag)
    {
        try
        {
            imag->Dimension(m_dim, REAL);
        }
        catch (hwMathException& except)
        {
            except.Status().SetArg1(2);
            throw;
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
}

// ****************************************************
//                 Slice Operations
// ****************************************************

//! Read a matrix slice from the calling object, as if the calling object is being
//! sliced on the right hand side of an equals sign
//! lhs_matrix = matrix(slice args)
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::SliceRHS(const std::vector<hwSliceArg>& sliceArg,
                                  hwTMatrixN<T1, T2>& lhsMatrix) const
{
    if (sliceArg.size() == 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    // ignore high singleton dimensions
    int numSlices = RelevantNumberOfSlices(sliceArg, false);

    // perform implicit RHS reshape with reduced slices 
    if (numSlices < m_dim.size())
    {
        std::vector<int> newDim(numSlices);

        for (int i = 0; i < numSlices-1; ++i)
            newDim[i] = m_dim[i];

        newDim[numSlices-1] = -1;
        void* dataPtr;

        if (IsReal())
            dataPtr = (void*) m_real;
        else
            dataPtr = (void*) m_complex;

        hwTMatrixN<T1, T2> reshaped(m_dim, dataPtr, Type());

        try
        {
            if (numSlices == 1 && sliceArg[0].IsColon())    // handle a(:)
            {
                newDim.push_back(1);
                std::vector<hwSliceArg> newSliceArg;
                newSliceArg.push_back(hwSliceArg());
                newSliceArg.push_back(0);
                reshaped.Reshape(newDim);
                reshaped.SliceRHS(newSliceArg, lhsMatrix);
                return;
            }

            reshaped.Reshape(newDim);
            reshaped.SliceRHS(sliceArg, lhsMatrix);
        }
        catch (hwMathException& except)
        {
            except.Status().ResetArgs();
            throw;
        }

        return;
    }

    // check RHS bounds
    try
    {
        BoundsCheckRHS(sliceArg, numSlices);
    }
    catch (hwMathException& except)
    {
        except.Status().ResetArgs();
        throw;
    }

    // dimension LHS
    int lhsNumDims;

    if (numSlices > 2)
    {
        for (lhsNumDims = numSlices-1; lhsNumDims > 1; --lhsNumDims)
        {
            // find last non scalar slice
            if (!sliceArg[lhsNumDims].IsScalar())
            {
                ++lhsNumDims;
                break;
            }
        }
    }
    else
    {
        lhsNumDims = 2;
    }

    lhsNumDims = _max(lhsNumDims, 2);

    std::vector<int> lhsMatrixDim(lhsNumDims);

    for (int i = 0; i < lhsNumDims; ++i)
    {
        if (sliceArg[i].IsScalar())
            lhsMatrixDim[i] = 1;
        else if (sliceArg[i].IsColon())
            lhsMatrixDim[i] = m_dim[i];
        else if (sliceArg[i].IsVector())
            lhsMatrixDim[i] = (int) sliceArg[i].Vector().size();
    }

    try
    {
        lhsMatrix.Dimension(lhsMatrixDim, Type());
    }
    catch (hwMathException& except)
    {
        except.Status().ResetArgs();
        throw;
    }

    int size = lhsMatrix.Size();

    if (size == 0)
        return;

    // determine the best approach, depending on whether the slicing involves
    // contiguous or equally spaced memory.
    // 0. leading contiguous dimensions (singleton dimensions, multiple colons)
    // 1. at least one colon
    // 2. at least one continguous vector
    // 3. everything else
    int S_case;
    int keyDim[4]  {-1, -1, -1, -1};
    int keySize[3] { 1,  1,  1};

    // analyze S_case = 0
    for (int j = 0; j < numSlices; ++j)
    {
        // keyDim[0]  = last contiguous dimension
        // keySize[0] = size of the leading contiguous dimensions
        if (sliceArg[j].IsScalar())
        {
            if (m_dim[j] != 1)
            {
                // keyDim[0] = j - 1;
                break;
            }
        }
        else if (sliceArg[j].IsColon())
        {
            keyDim[0] = j;
            keySize[0] *= m_dim[j];
        }
        else if (sliceArg[j].IsVector())
        {
            // keyDim[0] = j - 1;
            break;
        }
    }

    // analyze S_case = 1,2
    for (int j = 0; j < numSlices; ++j)
    {
        // keyDim[j]  = dimension of largest contiguous vector or colon
        // keySize[j] = size of dimension keyDim[j]
        if (sliceArg[j].IsColon())
        {
            if (m_dim[j] > keySize[1])
            {
                keyDim[1] = j;
                keySize[1] = m_dim[j];
            }
        }
        else if (sliceArg[j].IsVector())
        {
            int vecLength = static_cast<int> (sliceArg[j].Vector().size());
            int indx = sliceArg[j].Vector()[0];
            bool contigVec = true;

            for (int i = 1; i < vecLength; ++i)
            {
                if (sliceArg[j].Vector()[i] != ++indx)
                {
                    contigVec = false;
                    break;
                }
            }

            if (contigVec && vecLength > keySize[2])
            {
                keyDim[2] = j;
                keySize[2] = vecLength;
            }
        }
    }

    if (keySize[0] >= keySize[1] && keySize[0] >= keySize[2])
        S_case = 0;
    else if (keySize[1] >= keySize[2])
        S_case = 1;
    else if (keySize[2] > 1)
        S_case = 2;
    else
        S_case = 3;

    // simulate nested loops to iterate over the lhsMatrix elements
    // in order of contiguous memory location
    m_lhsMatrixIndex.clear();
    m_lhsMatrixIndex.resize(numSlices);
    m_rhsMatrixIndex.resize(numSlices);

    for (int i = 0; i < numSlices; ++i)
    {
        // set the rhsMatrix indices to the first index in each slice
        if (sliceArg[i].IsScalar())
            m_rhsMatrixIndex[i] = sliceArg[i].Scalar();
        else if (sliceArg[i].IsColon())
            m_rhsMatrixIndex[i] = 0;
        else if (sliceArg[i].IsVector())
            m_rhsMatrixIndex[i] = sliceArg[i].Vector()[0];
    }

    int startIndx;

    switch (S_case)
    {
    case 0:
        startIndx = keyDim[0] + 1;
        break;
    case 1:
    case 2:
    case 3:
        startIndx = 0;
    }

    for (int i = 0; i < size; ++i)
    {
        switch (S_case)
        {
        case 0:
        case 1:
        case 2:
        {
            int startLHS = lhsMatrix.Index(m_lhsMatrixIndex);
            int startRHS = Index(m_rhsMatrixIndex);
            int strideLHS = (S_case == 0) ? 1 : lhsMatrix.Stride(keyDim[S_case]);
            int strideRHS = (S_case == 0) ? 1 : Stride(keyDim[S_case]);

            if (IsReal())
                CopyData(lhsMatrix.m_real + startLHS, strideLHS, m_real + startRHS, strideRHS, keySize[S_case]);
            else
                CopyData(lhsMatrix.m_complex + startLHS, strideLHS, m_complex + startRHS, strideRHS, keySize[S_case]);

            i += keySize[S_case] - 1;
            break;
        }
        case 3:
        {
            if (sliceArg[0].IsScalar())
            {
                if (lhsMatrix.IsReal())
                    lhsMatrix(i) = (*this)(m_rhsMatrixIndex);
                else
                    lhsMatrix.z(i) = this->z(m_rhsMatrixIndex);
            }
            else if (sliceArg[0].IsVector())
            {
                int vecLen = static_cast<int> (sliceArg[0].Vector().size());

                m_rhsMatrixIndex[0] = 0;
                int topOfDim = Index(m_rhsMatrixIndex);
                int k = i;

                for (int j = 0; j < vecLen; ++j)
                {
                    int pos = topOfDim + sliceArg[0].Vector()[j];

                    if (IsReal())
                        lhsMatrix(k++) = (*this)(pos);
                    else
                        lhsMatrix.z(k++) = this->z(pos);
                }

                m_rhsMatrixIndex[0] = sliceArg[0].Vector()[0];

                i += vecLen - 1;
            }

            break;
        }
        }

        // advance lhs and rhs matrix indices
        for (int j = startIndx; j < numSlices; ++j)
        {
            if (sliceArg[j].IsScalar() || j == keyDim[S_case])
            {
                continue;
            }
            else if (sliceArg[j].IsColon())
            {
                // increment index j if possible
                if (m_rhsMatrixIndex[j] < m_dim[j] - 1)
                {
                    ++m_lhsMatrixIndex[j];
                    ++m_rhsMatrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_lhsMatrixIndex[j] = 0;
                m_rhsMatrixIndex[j] = 0;
            }
            else // if (sliceArg[j].IsVector())
            {
                // increment index j if possible
                if (m_lhsMatrixIndex[j] < (int)sliceArg[j].Vector().size() - 1)
                {
                    ++m_lhsMatrixIndex[j];
                    m_rhsMatrixIndex[j] = sliceArg[j].Vector()[m_lhsMatrixIndex[j]];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_lhsMatrixIndex[j] = 0;
                m_rhsMatrixIndex[j] = sliceArg[j].Vector()[0];
            }
        }
    }
}

//! Write to a matrix slice of the calling object, as if the calling object is being
//! sliced on the left hand side of an equals sign
//! matrix(slice args) = rhs_matrix
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::SliceLHS(const std::vector<hwSliceArg>& sliceArg,
                                  const hwTMatrixN<T1, T2>& rhsMatrix)
{
    if (sliceArg.size() == 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    if (IsReal() && !rhsMatrix.IsReal())
    {
        MakeComplex();
    }
    else if (!IsReal() && rhsMatrix.IsReal())
    {
        hwTMatrixN<T1, T2> temp;
        temp.PackComplex(rhsMatrix, nullptr);
        return SliceLHS(sliceArg, temp);
    }

    // handle empty rhsMatrix
    if (rhsMatrix.m_dim.size() == 2 && rhsMatrix.m_dim[0] == 0 && rhsMatrix.m_dim[1] == 0)
    {
        try
        {
            DeleteSlice(sliceArg);
            return;
        }
        catch (hwMathException& except)
        {
            except.Status().ResetArgs();
            throw;
        }
    }

    // ignore singleton higher dimensions
    int numSlices = RelevantNumberOfSlices(sliceArg, true);

    // perform implicit LHS reshape with reduced slices 
    if (numSlices < m_dim.size())
    {
        std::vector<int> newDim(numSlices);

        for (int i = 0; i < numSlices-1; ++i)
            newDim[i] = m_dim[i];

        newDim[numSlices-1] = -1;
        void* dataPtr;

        if (IsReal())
            dataPtr = (void*) m_real;
        else
            dataPtr = (void*) m_complex;

        hwTMatrixN<T1, T2> reshaped(m_dim, dataPtr, Type());

        try
        {
            if (numSlices == 1)
            {
                if (sliceArg[0].IsColon())    // handle a(:)
                {
                    newDim.push_back(1);
                    std::vector<hwSliceArg> newSliceArg;
                    newSliceArg.push_back(hwSliceArg());
                    newSliceArg.push_back(0);
                    reshaped.Reshape(newDim);
                    reshaped.SliceLHS(newSliceArg, rhsMatrix);
                    return;
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
                }
            }

            reshaped.Reshape(newDim);
            reshaped.SliceLHS(sliceArg, rhsMatrix);
        }
        catch (hwMathException& except)
        {
            except.Status().ResetArgs();
            throw;
        }

        return;
    }

    // map RHS dimensions to LHS dimensions
    bool newMatrix = false;

    if (m_dim.size() == 2 && m_dim[0] == 0 && m_dim[1] == 0)
        newMatrix = true;

    std::vector<int> newDim(numSlices);
    std::vector<int> dimMap(numSlices);     // rhsDim = dimMap[lhsDim]
    int lhsDim = 0;
    int rhsDim = 0;
    int numColons = 0;

    for (int i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsColon())
            ++numColons;
    }

    for (rhsDim = 0; rhsDim < rhsMatrix.m_dim.size(); ++rhsDim)
    {
        if (lhsDim >= newDim.size())
            break;

        if (sliceArg[lhsDim].IsScalar())
        {
            if (newMatrix)
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;
            else if (lhsDim < m_dim.size())
                newDim[lhsDim] = _max(m_dim[lhsDim], sliceArg[lhsDim].Scalar() + 1);
            else
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;

            if (rhsMatrix.m_dim[rhsDim] == 1)
            {
                dimMap[lhsDim++] = rhsDim;
            }
            else
            {
                dimMap[lhsDim++] = -1;
                --rhsDim;
            }
        }
        else if (sliceArg[lhsDim].IsColon())
        {
            if (newMatrix)
            {
                if (numColons == rhsMatrix.m_dim.size() || rhsMatrix.m_dim[rhsDim] != 1)
                {
                    newDim[lhsDim] = rhsMatrix.m_dim[rhsDim];
                    dimMap[lhsDim++] = rhsDim;
                }
            }
            else if (rhsMatrix.m_dim[rhsDim] == m_dim[lhsDim])
            {
                newDim[lhsDim] = rhsMatrix.m_dim[rhsDim];
                dimMap[lhsDim++] = rhsDim;
            }
            else if (rhsMatrix.m_dim[rhsDim] == 1)
            {
                // effectively applies squeeze(rhsMatrix)
            }
            else if (m_dim[lhsDim] == 1)
            {
                // effectively applies squeeze(lhsMatrix)
                newDim[lhsDim] = 1;
                dimMap[lhsDim++] = -1;
                --rhsDim;
            }
            else
            {
                throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
            }
        }
        else // sliceArg[lhsDim].IsVector()
        {
            int maxVectorDim = -1;

            for (int i = 0; i < sliceArg[lhsDim].Vector().size(); ++i)
            {
                if (sliceArg[lhsDim].Vector()[i] < 0)
                    throw hwMathException(HW_MATH_ERR_SLICE_INDEX, 1);

                maxVectorDim = _max(sliceArg[lhsDim].Vector()[i] + 1, maxVectorDim);
            }

            if (newMatrix)
            {
                newDim[lhsDim] = maxVectorDim;
                dimMap[lhsDim++] = rhsDim;
            }
            else if (rhsMatrix.m_dim[rhsDim] == sliceArg[lhsDim].Vector().size() &&
                     lhsDim < m_dim.size())
            {
                newDim[lhsDim] = _max(m_dim[lhsDim], maxVectorDim);
                dimMap[lhsDim++] = rhsDim;
            }
            else if (rhsMatrix.m_dim[rhsDim] == 1)
            {
                // effectively applies squeeze(rhsMatrix)
            }
            else
            {
                throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
            }
        }
    }

    for ( ; rhsDim < rhsMatrix.m_dim.size(); ++rhsDim)
    {
        if (rhsMatrix.m_dim[rhsDim] != 1)
            throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
    }
        
    for ( ; lhsDim < numSlices; ++lhsDim)
    {
        if (sliceArg[lhsDim].IsScalar())
        {
            if (lhsDim < m_dim.size())
                newDim[lhsDim] = _max(m_dim[lhsDim], sliceArg[lhsDim].Scalar() + 1);
            else
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;

            dimMap[lhsDim] = -1;
        }
        else if (sliceArg[lhsDim].IsColon())
        {
            if (!newMatrix && lhsDim < m_dim.size() && m_dim[lhsDim] != 1)
                throw hwMathException(HW_MATH_ERR_SLICE_INDEX);

            newDim[lhsDim] = 1;
            dimMap[lhsDim] = -1;
        }
        else // sliceArg[lhsDim].IsVector()
        {
            throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
        }
    }

    // discard any singleton dimensions that may have been created
    // due to excess colon and "1" slices
    while (newDim.size() > 2 && newDim.back() == 1)
    {
        newDim.pop_back();
        dimMap.pop_back();
        numSlices--;
    }

    // dimension the matrix
    if (newMatrix)
    {
        m_dim = newDim;

        if (rhsMatrix.IsReal())
        {
            Allocate(REAL);
        }
        else
        {
            Allocate(COMPLEX);
            m_bits.realData = 0;
        }

        SetElements(0.0);
    }
    else
    {
        bool resize = false;

        if (newDim.size() > m_dim.size())
        {
            resize = true;
        }
        else
        {
            for (int i = 0; i < newDim.size(); ++i)
            {
                if (newDim[i] > m_dim[i])
                {
                    resize = true;
                    break;
                }
            }
        }

        if (resize)
        {
            GrowLHSMatrix(newDim);
        }
    }

    int size = rhsMatrix.Size();

    if (size == 0)
        return;

    if (m_size < size)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    // determine the best approach, depending on whether the slicing involves
    // contiguous or equally spaced memory.
    // 0. leading contiguous dimensions (singleton dimensions, multiple colons)
    // 1. at least one colon
    // 2. at least one continguous vector
    // 3. everything else
    int S_case;
    int keyDim[4]{ -1, -1, -1, -1 };
    int keySize[3]{ 1,  1,  1 };

    // analyze S_case = 0
    for (int j = 0; j < numSlices; ++j)
    {
        // keyDim[0]  = last contiguous dimension
        // keySize[0] = size of the leading contiguous dimensions
        if (sliceArg[j].IsScalar())
        {
            if (m_dim[j] != 1)
            {
                // keyDim[0] = j - 1;
                break;
            }
        }
        else if (sliceArg[j].IsColon())
        {
            keyDim[0] = j;
            keySize[0] *= m_dim[j];
        }
        else if (sliceArg[j].IsVector())
        {
            // keyDim[0] = j - 1;
            break;
        }
    }

    // analyze S_case = 1,2
    for (int j = 0; j < numSlices; ++j)
    {
        // keyDim[j]  = dimension of largest contiguous vector or colon
        // keySize[j] = size of dimension keyDim[j]
        if (sliceArg[j].IsColon())
        {
            if (m_dim[j] > keySize[1])
            {
                keyDim[1] = j;
                keySize[1] = m_dim[j];
            }
        }
        else if (sliceArg[j].IsVector())
        {
            int vecLength = static_cast<int> (sliceArg[j].Vector().size());
            int indx = sliceArg[j].Vector()[0];
            bool contigVec = true;

            for (int i = 1; i < vecLength; ++i)
            {
                if (sliceArg[j].Vector()[i] != ++indx)
                {
                    contigVec = false;
                    break;
                }
            }

            if (contigVec && vecLength > keySize[2])
            {
                keyDim[2] = j;
                keySize[2] = vecLength;
            }
        }
    }

    if (keySize[0] >= keySize[1] && keySize[0] >= keySize[2])
        S_case = 0;
    else if (keySize[1] >= keySize[2])
        S_case = 1;
    else if (keySize[2] > 1)
        S_case = 2;
    else
        S_case = 3;

    // simulate nested loops to iterate over the rhsMatrix elements
    // in order of contiguous memory location, copying blocks where possible
    m_rhsMatrixIndex.clear();
    m_rhsMatrixIndex.resize(numSlices);
    m_lhsMatrixIndex.clear();
    m_lhsMatrixIndex.resize(numSlices);

    for (int i = 0; i < numSlices; ++i)
    {
        // set the lhsMatrix indices to the first index in each slice
        if (sliceArg[i].IsScalar())
            m_lhsMatrixIndex[i] = sliceArg[i].Scalar();
        else if (sliceArg[i].IsColon())
            m_lhsMatrixIndex[i] = 0;
        else if (sliceArg[i].IsVector())
            m_lhsMatrixIndex[i] = sliceArg[i].Vector()[0];
    }

    int startIndx = (S_case == 0) ? keyDim[0] + 1 : 0;
    int strideLHS;
    int strideRHS;

    if (S_case != 3)
    {
        strideLHS = (S_case == 0) ? 1 : Stride(keyDim[S_case]);
        strideRHS = (S_case == 0) ? 1 : rhsMatrix.Stride(dimMap[keyDim[S_case]]);
    }

    for (int i = 0; i < size; ++i)
    {
        switch (S_case)
        {
        case 0:
        case 1:
        case 2:
        {
            int startLHS = Index(m_lhsMatrixIndex);
            int startRHS = rhsMatrix.Index(m_rhsMatrixIndex);

            if (IsReal())
                CopyData(m_real + startLHS, strideLHS, rhsMatrix.m_real + startRHS, strideRHS, keySize[S_case]);
            else
                CopyData(m_complex + startLHS, strideLHS, rhsMatrix.m_complex + startRHS, strideRHS, keySize[S_case]);

            i += keySize[S_case] - 1;
            break;
        }
        case 3:
        {
            if (sliceArg[0].IsScalar())
            {
                if (IsReal())
                    (*this)(m_lhsMatrixIndex) = rhsMatrix(i);
                else
                    this->z(m_lhsMatrixIndex) = rhsMatrix.z(i);
            }
            else if (sliceArg[0].IsVector())
            {
                int vecLen = static_cast<int> (sliceArg[0].Vector().size());

                m_lhsMatrixIndex[0] = 0;
                int topOfDim = Index(m_lhsMatrixIndex);
                int k = i;

                for (int j = 0; j < vecLen; ++j)
                {
                    int pos = topOfDim + sliceArg[0].Vector()[j];

                    if (IsReal())
                        (*this)(pos) = rhsMatrix(k++);
                    else
                        this->z(pos) = rhsMatrix.z(k++);
                }

                m_lhsMatrixIndex[0] = sliceArg[0].Vector()[0];

                i += vecLen - 1;
            }

            break;
        }
        }

        // advance rhs matrix indices
        for (int j = startIndx; j < numSlices; ++j)
        {
            if (sliceArg[j].IsScalar() || j == keyDim[S_case])
            {
                continue;
            }
            else if (sliceArg[j].IsColon())
            {
                // increment index j if possible
                if (m_lhsMatrixIndex[j] < (int) m_dim[j]-1)
                {
                    ++m_rhsMatrixIndex[dimMap[j]];
                    ++m_lhsMatrixIndex[j];

                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_rhsMatrixIndex[dimMap[j]] = 0;
                m_lhsMatrixIndex[j] = 0;
            }
            else if (sliceArg[j].IsVector())
            {
                // increment index j if possible
                if (m_rhsMatrixIndex[dimMap[j]] < (int) sliceArg[j].Vector().size()-1)
                {
                    ++m_rhsMatrixIndex[dimMap[j]];
                    m_lhsMatrixIndex[j] = sliceArg[j].Vector()[m_rhsMatrixIndex[dimMap[j]]];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_rhsMatrixIndex[dimMap[j]] = 0;
                m_lhsMatrixIndex[j] = sliceArg[j].Vector()[0];
            }
        }
    }
}

//! Write to a matrix slice of the calling object, as if the calling object is being
//! sliced on the left hand side of an equals sign
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::SliceLHS(const std::vector<hwSliceArg>& sliceArg, T1 real)
{
    if (sliceArg.size() == 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    if (!IsReal())
    {
        try
        {
            SliceLHS(sliceArg, hwTComplex<T1>(real, 0.0));
        }
        catch (hwMathException& except)
        {
            except.Status().ResetArgs();
            throw;
        }

        return;
    }

    // ignore high singleton dimensions
    int numSlices = RelevantNumberOfSlices(sliceArg, false);

    // manage single indexing assigment
    if (numSlices == 1 && sliceArg[0].IsScalar())
    {
        int index = sliceArg[0].Scalar();

        if (index > m_size - 1)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

        (*this)(index) = real;

        return;
    }

    // perform implicit LHS reshape with reduced slices 
    if (numSlices < m_dim.size())
    {
        std::vector<int> newDim(numSlices);

        for (int i = 0; i < numSlices - 1; ++i)
            newDim[i] = m_dim[i];

        newDim[numSlices - 1] = -1;
        void* dataPtr;

        if (IsReal())
            dataPtr = (void*)m_real;
        else
            dataPtr = (void*)m_complex;

        hwTMatrixN<T1, T2> reshaped(m_dim, dataPtr, Type());

        try
        {
            if (numSlices == 1 && sliceArg[0].IsColon())    // handle a(:)
            {
                newDim.push_back(1);
                std::vector<hwSliceArg> newSliceArg;
                newSliceArg.push_back(hwSliceArg());
                newSliceArg.push_back(0);
                reshaped.Reshape(newDim);
                reshaped.SliceLHS(newSliceArg, real);
                return;
            }

            reshaped.Reshape(newDim);
            reshaped.SliceLHS(sliceArg, real);

            if (reshaped.m_real != m_real)  // reallocation occurred
            {
                delete[] m_real;
                m_real = reshaped.m_real;
                reshaped.m_bits.ownData = 0;

                for (int i = 0; i < numSlices - 1; ++i)
                    m_dim[i] = reshaped.m_dim[i];

                ComputeSize();
            }
        }
        catch (hwMathException& except)
        {
            except.Status().ResetArgs();
            throw;
        }

        return;
    }

    // map RHS dimensions to LHS dimensions
    bool newMatrix = false;

    if (m_dim.size() == 2 && m_dim[0] == 0 && m_dim[1] == 0)
        newMatrix = true;

    std::vector<int> newDim(numSlices);
    int lhsDim = 0;

    for (lhsDim = 0; lhsDim < numSlices; ++lhsDim)
    {
        if (sliceArg[lhsDim].IsScalar())
        {
            if (newMatrix)
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;
            else if (lhsDim < m_dim.size())
                newDim[lhsDim] = _max(m_dim[lhsDim], sliceArg[lhsDim].Scalar() + 1);
            else
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;
        }
        else if (sliceArg[lhsDim].IsColon())
        {
            if (newMatrix)
            {
                newDim[lhsDim] = 1;
            }
            else
            {
                newDim[lhsDim] = m_dim[lhsDim];
            }
        }
        else // sliceArg[lhsDim].IsVector()
        {
            int maxVectorDim = -1;

            for (int i = 0; i < sliceArg[lhsDim].Vector().size(); ++i)
            {
                if (sliceArg[lhsDim].Vector()[i] < 0)
                    throw hwMathException(HW_MATH_ERR_SLICE_INDEX, 1);

                maxVectorDim = _max(sliceArg[lhsDim].Vector()[i] + 1, maxVectorDim);
            }

            if (newMatrix)
            {
                newDim[lhsDim] = maxVectorDim;
            }
            else if (lhsDim < m_dim.size())
            {
                newDim[lhsDim] = _max(m_dim[lhsDim], maxVectorDim);
            }
            else
            {
                throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
            }
        }
    }

    for (; lhsDim < numSlices; ++lhsDim)
    {
        if (sliceArg[lhsDim].IsScalar())
        {
            if (lhsDim < m_dim.size())
                newDim[lhsDim] = _max(m_dim[lhsDim], sliceArg[lhsDim].Scalar() + 1);
            else
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;
        }
        else if (sliceArg[lhsDim].IsColon())
        {
            newDim[lhsDim] = 1;
        }
        else // sliceArg[lhsDim].IsVector()
        {
            throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
        }
    }

    // discard any singleton dimensions that may have been created
    // due to excess colon and "1" slices
    while (newDim.size() > 2 && newDim.back() == 1)
    {
        newDim.pop_back();
        numSlices--;
    }

    // dimension the matrix
    if (newMatrix)
    {
        m_dim = newDim;
        Allocate(REAL);
        SetElements(0.0);
    }
    else
    {
        bool resize = false;

        if (newDim.size() > m_dim.size())
        {
            resize = true;
        }
        else
        {
            for (int i = 0; i < newDim.size(); ++i)
            {
                if (newDim[i] > m_dim[i])
                {
                    resize = true;
                    break;
                }
            }
        }

        if (resize)
        {
            GrowLHSMatrix(newDim);
        }
    }

    if (m_size == 0)
        return;

    int size = 1;

    m_rhsMatrixIndex.clear();
    m_rhsMatrixIndex.resize(numSlices);

    for (int i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar()) {}
        else if (sliceArg[i].IsColon())
            size *= m_dim[i];
        else if (sliceArg[i].IsVector())
            size *= (int)sliceArg[i].Vector().size();
    }

    if (size == 0)
        return;

    // determine the best approach, depending on whether the slicing involves
    // contiguous or equally spaced memory.
    // 0. leading contiguous dimensions (singleton dimensions, multiple colons)
    // 1. at least one colon
    // 2. at least one continguous vector
    // 3. everything else
    int S_case;
    int keyDim[4]{ -1, -1, -1, -1 };
    int keySize[3]{ 1,  1,  1 };

    // analyze S_case = 0
    for (int j = 0; j < numSlices; ++j)
    {
        // keyDim[0]  = last contiguous dimension
        // keySize[0] = size of the leading contiguous dimensions
        if (sliceArg[j].IsScalar())
        {
            if (m_dim[j] != 1)
            {
                // keyDim[0] = j - 1;
                break;
            }
        }
        else if (sliceArg[j].IsColon())
        {
            keyDim[0] = j;
            keySize[0] *= m_dim[j];
        }
        else if (sliceArg[j].IsVector())
        {
            // keyDim[0] = j - 1;
            break;
        }
    }

    // analyze S_case = 1,2
    for (int j = 0; j < numSlices; ++j)
    {
        // keyDim[j]  = dimension of largest contiguous vector or colon
        // keySize[j] = size of dimension keyDim[j]
        if (sliceArg[j].IsColon())
        {
            if (m_dim[j] > keySize[1])
            {
                keyDim[1] = j;
                keySize[1] = m_dim[j];
            }
        }
        else if (sliceArg[j].IsVector())
        {
            int vecLength = static_cast<int> (sliceArg[j].Vector().size());
            int indx = sliceArg[j].Vector()[0];
            bool contigVec = true;

            for (int i = 1; i < vecLength; ++i)
            {
                if (sliceArg[j].Vector()[i] != ++indx)
                {
                    contigVec = false;
                    break;
                }
            }

            if (contigVec && vecLength > keySize[2])
            {
                keyDim[2] = j;
                keySize[2] = vecLength;
            }
        }
    }

    if (keySize[0] >= keySize[1] && keySize[0] >= keySize[2])
        S_case = 0;
    else if (keySize[1] >= keySize[2])
        S_case = 1;
    else if (keySize[2] > 1)
        S_case = 2;
    else
        S_case = 3;

    // simulate nested loops to iterate over the rhsMatrix elements
    // in order of contiguous memory location, copying blocks where possible
    m_lhsMatrixIndex.resize(numSlices);

    for (int i = 0; i < numSlices; ++i)
    {
        // set the lhsMatrix indices to the first index in each slice
        if (sliceArg[i].IsScalar())
            m_lhsMatrixIndex[i] = sliceArg[i].Scalar();
        else if (sliceArg[i].IsColon())
            m_lhsMatrixIndex[i] = 0;
        else if (sliceArg[i].IsVector())
            m_lhsMatrixIndex[i] = sliceArg[i].Vector()[0];
    }

    int startIndx = (S_case == 0) ? keyDim[0] + 1 : 0;
    int strideLHS;

    if (S_case != 3)
    {
        strideLHS = (S_case == 0) ? 1 : Stride(keyDim[S_case]);
    }

    for (int i = 0; i < size; ++i)
    {
        switch (S_case)
        {
        case 0:
        case 1:
        case 2:
        {
            int startLHS = Index(m_lhsMatrixIndex);

            if (IsReal())
            {
                T1* dataPtr = m_real + startLHS;

                for (int k = 0; k < keySize[S_case]; ++k)
                {
                    *dataPtr = real;
                    dataPtr += strideLHS;
                }
            }
            else
            {
                T2* dataPtr = m_complex + startLHS;

                for (int k = 0; k < keySize[S_case]; ++k)
                {
                    *dataPtr = real;
                    dataPtr += strideLHS;
                }
            }

            i += keySize[S_case] - 1;

            break;
        }
        case 3:
        {
            if (sliceArg[0].IsScalar())
            {
                if (IsReal())
                    (*this)(m_lhsMatrixIndex) = real;
                else
                    this->z(m_lhsMatrixIndex) = real;
            }
            else if (sliceArg[0].IsVector())
            {
                int vecLen = static_cast<int> (sliceArg[0].Vector().size());

                m_lhsMatrixIndex[0] = 0;
                int topOfDim = Index(m_lhsMatrixIndex);

                for (int j = 0; j < vecLen; ++j)
                {
                    int pos = topOfDim + sliceArg[0].Vector()[j];

                    if (IsReal())
                        (*this)(pos) = real;
                    else
                        this->z(pos) = real;
                }

                m_lhsMatrixIndex[0] = sliceArg[0].Vector()[0];

                i += vecLen - 1;
                break;
            }
        }
        }

        // advance rhs matrix indices, as if the RHS is real*ones(slice_dims)
        // there is probably a better way to do this
        for (int j = startIndx; j < numSlices; ++j)
        {
            if (sliceArg[j].IsScalar() || j == keyDim[S_case])
            {
                continue;
            }
            else if (sliceArg[j].IsColon())
            {
                // increment index j if possible
                if (m_lhsMatrixIndex[j] < (int)m_dim[j] - 1)
                {
                    ++m_lhsMatrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_lhsMatrixIndex[j] = 0;
            }
            else if (sliceArg[j].IsVector())
            {
                // increment index j if possible
                int vecSize = (int)sliceArg[j].Vector().size();

                if (m_rhsMatrixIndex[j] < vecSize - 1)
                {
                    ++m_rhsMatrixIndex[j];
                    m_lhsMatrixIndex[j] = sliceArg[j].Vector()[m_rhsMatrixIndex[j]];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_rhsMatrixIndex[j] = 0;
                m_lhsMatrixIndex[j] = sliceArg[j].Vector()[0];
            }
        }
    }
}

//! Write to a matrix slice of the calling object, as if the calling object is being
//! sliced on the left hand side of an equals sign
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::SliceLHS(const std::vector<hwSliceArg>& sliceArg, const T2& cmplx)
{
    if (sliceArg.size() == 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    try
    {
        if (cmplx.Imag() == (T1) 0)
        {
            if (IsReal())
            {
                SliceLHS(sliceArg, cmplx.Real());
                return;
            }
        }
        else if (IsReal())
        {
            MakeComplex();
        }
    }
    catch (hwMathException& except)
    {
        except.Status().ResetArgs();
        throw;
    }

    // ignore high singleton dimensions
    int numSlices = RelevantNumberOfSlices(sliceArg, false);

    // manage single indexing assigment
    if (numSlices == 1 && sliceArg[0].IsScalar())
    {
        int index = sliceArg[0].Scalar();

        if (index > m_size - 1)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

        this->z(index) = cmplx;

        return;
    }

    // perform implicit LHS reshape with reduced slices 
    if (numSlices < m_dim.size())
    {
        std::vector<int> newDim(numSlices);

        for (int i = 0; i < numSlices - 1; ++i)
            newDim[i] = m_dim[i];

        newDim[numSlices - 1] = -1;
        void* dataPtr;

        if (IsReal())
            dataPtr = (void*)m_real;
        else
            dataPtr = (void*)m_complex;

        hwTMatrixN<T1, T2> reshaped(m_dim, dataPtr, Type());

        try
        {
            if (numSlices == 1 && sliceArg[0].IsColon())    // handle a(:)
            {
                newDim.push_back(1);
                std::vector<hwSliceArg> newSliceArg;
                newSliceArg.push_back(hwSliceArg());
                newSliceArg.push_back(0);
                reshaped.Reshape(newDim);
                reshaped.SliceLHS(newSliceArg, cmplx);

                if (reshaped.m_complex != m_complex)
                {
                    delete[] m_complex;
                    m_complex = reshaped.m_complex;
                    reshaped.m_bits.ownData = 0;
                }

                return;
            }

            reshaped.Reshape(newDim);
            reshaped.SliceLHS(sliceArg, cmplx);
        }
        catch (hwMathException& except)
        {
            except.Status().ResetArgs();
            throw;
        }

        return;
    }

    // map RHS dimensions to LHS dimensions
    bool newMatrix = false;

    if (m_dim.size() == 2 && m_dim[0] == 0 && m_dim[1] == 0)
        newMatrix = true;

    std::vector<int> newDim(numSlices);
    int lhsDim = 0;

    for (lhsDim = 0; lhsDim < numSlices; ++lhsDim)
    {
        if (sliceArg[lhsDim].IsScalar())
        {
            if (newMatrix)
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;
            else if (lhsDim < m_dim.size())
                newDim[lhsDim] = _max(m_dim[lhsDim], sliceArg[lhsDim].Scalar() + 1);
            else
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;
        }
        else if (sliceArg[lhsDim].IsColon())
        {
            if (newMatrix)
            {
                newDim[lhsDim] = 1;
            }
            else
            {
                newDim[lhsDim] = m_dim[lhsDim];
            }
        }
        else // sliceArg[lhsDim].IsVector()
        {
            int maxVectorDim = -1;

            for (int i = 0; i < sliceArg[lhsDim].Vector().size(); ++i)
            {
                if (sliceArg[lhsDim].Vector()[i] < 0)
                    throw hwMathException(HW_MATH_ERR_SLICE_INDEX, 1);

                maxVectorDim = _max(sliceArg[lhsDim].Vector()[i] + 1, maxVectorDim);
            }

            if (newMatrix)
            {
                newDim[lhsDim] = maxVectorDim;
            }
            else if (lhsDim < m_dim.size())
            {
                newDim[lhsDim] = _max(m_dim[lhsDim], maxVectorDim);
            }
            else
            {
                throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
            }
        }
    }

    for (; lhsDim < numSlices; ++lhsDim)
    {
        if (sliceArg[lhsDim].IsScalar())
        {
            if (lhsDim < m_dim.size())
                newDim[lhsDim] = _max(m_dim[lhsDim], sliceArg[lhsDim].Scalar() + 1);
            else
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;
        }
        else if (sliceArg[lhsDim].IsColon())
        {
            newDim[lhsDim] = 1;
        }
        else // sliceArg[lhsDim].IsVector()
        {
            throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
        }
    }

    // discard any singleton dimensions that may have been created
    // due to excess colon and "1" slices
    while (newDim.size() > 2 && newDim.back() == 1)
    {
        newDim.pop_back();
        numSlices--;
    }

    // dimension the matrix
    if (newMatrix)
    {
        m_dim = newDim;
        Allocate(COMPLEX);
        SetElements(0.0);
    }
    else
    {
        bool resize = false;

        if (newDim.size() > m_dim.size())
        {
            resize = true;
        }
        else
        {
            for (int i = 0; i < newDim.size(); ++i)
            {
                if (newDim[i] > m_dim[i])
                {
                    resize = true;
                    break;
                }
            }
        }

        if (resize)
        {
            GrowLHSMatrix(newDim);
        }
    }

    if (m_size == 0)
        return;

    int size = 1;

    m_rhsMatrixIndex.clear();
    m_rhsMatrixIndex.resize(numSlices);

    for (int i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar()) {}
        else if (sliceArg[i].IsColon())
            size *= m_dim[i];
        else if (sliceArg[i].IsVector())
            size *= (int)sliceArg[i].Vector().size();
    }

    if (size == 0)
        return;

    // determine the best approach, depending on whether the slicing involves
    // contiguous or equally spaced memory.
    // 0. leading contiguous dimensions (singleton dimensions, multiple colons)
    // 1. at least one colon
    // 2. at least one continguous vector
    // 3. everything else
    int S_case;
    int keyDim[4]{ -1, -1, -1, -1 };
    int keySize[3]{ 1,  1,  1 };

    // analyze S_case = 0
    for (int j = 0; j < numSlices; ++j)
    {
        // keyDim[0]  = last contiguous dimension
        // keySize[0] = size of the leading contiguous dimensions
        if (sliceArg[j].IsScalar())
        {
            if (m_dim[j] != 1)
            {
                // keyDim[0] = j - 1;
                break;
            }
        }
        else if (sliceArg[j].IsColon())
        {
            keyDim[0] = j;
            keySize[0] *= m_dim[j];
        }
        else if (sliceArg[j].IsVector())
        {
            // keyDim[0] = j - 1;
            break;
        }
    }

    // analyze S_case = 1,2
    for (int j = 0; j < numSlices; ++j)
    {
        // keyDim[j]  = dimension of largest contiguous vector or colon
        // keySize[j] = size of dimension keyDim[j]
        if (sliceArg[j].IsColon())
        {
            if (m_dim[j] > keySize[1])
            {
                keyDim[1] = j;
                keySize[1] = m_dim[j];
            }
        }
        else if (sliceArg[j].IsVector())
        {
            int vecLength = static_cast<int> (sliceArg[j].Vector().size());
            int indx = sliceArg[j].Vector()[0];
            bool contigVec = true;

            for (int i = 1; i < vecLength; ++i)
            {
                if (sliceArg[j].Vector()[i] != ++indx)
                {
                    contigVec = false;
                    break;
                }
            }

            if (contigVec && vecLength > keySize[2])
            {
                keyDim[2] = j;
                keySize[2] = vecLength;
            }
        }
    }

    if (keySize[0] >= keySize[1] && keySize[0] >= keySize[2])
        S_case = 0;
    else if (keySize[1] >= keySize[2])
        S_case = 1;
    else if (keySize[2] > 1)
        S_case = 2;
    else
        S_case = 3;

    // simulate nested loops to iterate over the rhsMatrix elements
    // in order of contiguous memory location, copying blocks where possible
    m_lhsMatrixIndex.resize(numSlices);

    for (int i = 0; i < numSlices; ++i)
    {
        // set the lhsMatrix indices to the first index in each slice
        if (sliceArg[i].IsScalar())
            m_lhsMatrixIndex[i] = sliceArg[i].Scalar();
        else if (sliceArg[i].IsColon())
            m_lhsMatrixIndex[i] = 0;
        else if (sliceArg[i].IsVector())
            m_lhsMatrixIndex[i] = sliceArg[i].Vector()[0];
    }

    int startIndx = (S_case == 0) ? keyDim[0] + 1 : 0;
    int strideLHS;

    if (S_case != 3)
    {
        strideLHS = (S_case == 0) ? 1 : Stride(keyDim[S_case]);
    }

    for (int i = 0; i < size; ++i)
    {
        switch (S_case)
        {
        case 0:
        case 1:
        case 2:
        {
            int startLHS = Index(m_lhsMatrixIndex);

            T2* dataPtr = m_complex + startLHS;

            for (int k = 0; k < keySize[S_case]; ++k)
            {
                *dataPtr = cmplx;
                dataPtr += strideLHS;
            }

            i += keySize[S_case] - 1;
            break;
        }
        case 3:
        {
            if (sliceArg[0].IsScalar())
            {
                this->z(m_lhsMatrixIndex) = cmplx;
            }
            else if (sliceArg[0].IsVector())
            {
                int vecLen = static_cast<int> (sliceArg[0].Vector().size());

                m_lhsMatrixIndex[0] = 0;
                int topOfDim = Index(m_lhsMatrixIndex);

                for (int j = 0; j < vecLen; ++j)
                {
                    int pos = topOfDim + sliceArg[0].Vector()[j];
                    this->z(pos) = cmplx;
                }

                m_lhsMatrixIndex[0] = sliceArg[0].Vector()[0];

                i += vecLen - 1;
            }

            break;
        }
        }

        // advance rhs matrix indices, as if the RHS is real*ones(slice_dims)
        // there is probably a better way to do this
        for (int j = startIndx; j < numSlices; ++j)
        {
            if (sliceArg[j].IsScalar() || j == keyDim[S_case])
            {
                continue;
            }
            else if (sliceArg[j].IsColon())
            {
                // increment index j if possible
                if (m_lhsMatrixIndex[j] < (int)m_dim[j] - 1)
                {
                    ++m_lhsMatrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_lhsMatrixIndex[j] = 0;
            }
            else if (sliceArg[j].IsVector())
            {
                // increment index j if possible
                int vecSize = (int)sliceArg[j].Vector().size();

                if (m_rhsMatrixIndex[j] < vecSize - 1)
                {
                    ++m_rhsMatrixIndex[j];
                    m_lhsMatrixIndex[j] = sliceArg[j].Vector()[m_rhsMatrixIndex[j]];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_rhsMatrixIndex[j] = 0;
                m_lhsMatrixIndex[j] = sliceArg[j].Vector()[0];
            }
        }
    }
}

// ****************************************************
//               Arithmetic Operators
// ****************************************************
//! Implement the == operator
template<typename T1, typename T2>
bool hwTMatrixN<T1, T2>::operator==(const hwTMatrixN<T1, T2>& A) const
{
    int effDimSize = _min((int) m_dim.size(), (int) A.m_dim.size());

    for (int i = 0; i < effDimSize; ++i)
    {
        if (m_dim[i] != A.m_dim[i])
            return false;
    }

    for (int i = effDimSize; i < m_dim.size(); ++i)
    {
        if (m_dim[i] != 1)
            return false;
    }

    for (int i = effDimSize; i < A.m_dim.size(); ++i)
    {
        if (A.m_dim[i] != 1)
            return false;
    }

    if (IsReal())
    {
        if (A.IsReal())
        {
            for (int i = 0; i < m_size; ++i)
            {
                if (m_real[i] != A.m_real[i])
                    return false;
            }
        }
        else
        {
            for (int i = 0; i < m_size; ++i)
            {
                if (A.m_complex[i] != m_real[i])
                    return false;
            }
        }
    }
    else if (m_complex)
    {
        if (A.IsReal())
        {
            for (int i = 0; i < m_size; ++i)
            {
                if (m_complex[i] != A.m_real[i])
                    return false;
            }
        }
        else
        {
            for (int i = 0; i < m_size; ++i)
            {
                if (m_complex[i] != A.m_complex[i])
                    return false;
            }
        }
    }

    return true;
}
    
//! Implement the != operator
template<typename T1, typename T2>
bool hwTMatrixN<T1, T2>::operator!=(const hwTMatrixN<T1, T2>& A) const
{
    return !((*this) == A);
}

// ****************************************************
//                 Private Utilities
// ****************************************************

//! Set Dimensions
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::SetDimensions(const std::vector<int>& dim)
{
    int size = (int) dim.size();

    if (!size)
    {
        MakeEmpty();
    }
    else if (size == 1)
    {
        m_dim.resize(2);
        m_dim[0] = dim[0];
        m_dim[1] = 1;
    }
    else
    {
        int numDims = 2;

        // ignore trailing dimensions greater than 2 of size 1
        for (int i = (int) size-1; i > 1; --i)
        {
            if (dim[i] != 1)
            {
                numDims = i+1;
                break;
            }
        }

        m_dim.resize(numDims);

        for (int i = 0; i < numDims; ++i)
            m_dim[i] = dim[i];
    }
}

//! Compute size
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::ComputeSize()
{
    if (m_dim.empty())
    {
        m_size = 0;
        return;
    }

    if (m_dim[0] < 0)
    {
        m_size = -1;
        return;
    }

    if (m_dim[0] == 0)
    {
        m_size = 0;
        return;
    }

    constexpr int maxSize = (std::numeric_limits<int>::max)();
    
    m_size = m_dim[0];

    for (int i = 1; i < m_dim.size(); ++i)
    {
        if (m_dim[i] < 0)
        {
            m_size = -1;
            return;
        }

        if (m_dim[i] == 0)
        {
            m_size = 0;
            return;
        }

        // detect overflow
        if (m_dim[i] >= maxSize / m_size)
        {
            m_size = -1;
            return;
        }

        m_size *= m_dim[i];
    }
}

//! Compute memory capacity for the matrix data
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::SetCapacity(DataType dataType)
{
    if (m_size > m_capacity)
    {
        m_capacity = _max(2 * m_capacity, m_size);
    }
    else if (m_size < 0)
    {
        m_capacity = -10;    // force allocaction failure
    }
}

//! Allocate memory for the matrix data
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Allocate(DataType dataType)
{
    ComputeSize();
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
            throw hwMathException(HW_MATH_ERR_ALLOCFAILED);
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
            throw hwMathException(HW_MATH_ERR_ALLOCFAILED);
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

//! Delete the matrix data
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Deallocate()
{
    if (m_bits.ownData)
    {
        if (m_real_memory)
            delete[] m_real_memory;
        else if (m_real)
            delete[] m_real;	    // ownership of external data was assumed

        if (m_complex_memory)
            delete[] m_complex_memory;
        else if (m_complex)
            delete[] m_complex;    // ownership of external data was assumed
    }
}

//! Set matrix to empty condition
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::MakeEmpty()
{
    Deallocate();
    m_capacity = 0;
    m_real_memory = nullptr;
    m_complex_memory = nullptr;
    m_real = nullptr;
    m_complex = nullptr;
    m_size = 0;
    m_bits.ownData = 1;
    m_bits.realData = 1;
    m_dim.resize(2);
    m_dim[0] = 0;
    m_dim[1] = 0;
}

//! Copy function
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Copy(const hwTMatrixN<T1, T2>& source)
{
    if (IsEmpty())
    {
        Dimension(source.m_dim, source.Type());
        m_bits.ownData = 1;

        if (m_real)
            CopyData(m_real, source.m_real, m_size);
        else if (m_complex)
            CopyData(m_complex, source.m_complex, m_size);
    }
    else if (m_dim == source.m_dim && m_bits.realData == source.m_bits.realData)
    {
        if (m_real)
            CopyData(m_real, source.m_real, m_size);
        else if (m_complex)
            CopyData(m_complex, source.m_complex, m_size);
    }
    else
    {
        MakeEmpty();
        Copy(source);
    }
}

template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::CopyData(void* dest, const void* src, int count) const
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

//! Copy data
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::CopyData(void* dest, int stride_dest,
                                  const void* src, int stride_src, int count) const
{
    // memcpy cannot be used here
    if (m_real)
    {
        T1* destTemp = reinterpret_cast<T1*>(dest);
        const T1* srcTemp = reinterpret_cast<T1*>(const_cast<void*>(src));

        for (int i = 0; i < count; ++i)
            destTemp[i * stride_dest] = srcTemp[i * stride_src];
    }
    else if (m_complex)
    {
        T2* destTemp = reinterpret_cast<T2*>(dest);
        const T2* srcTemp = reinterpret_cast<T2*>(const_cast<void*>(src));

        for (int i = 0; i < count; ++i)
            destTemp[i * stride_dest] = srcTemp[i * stride_src];
    }
}

//! Delete a matrix slice
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::DeleteSlice(const std::vector<hwSliceArg>& sliceArg)
{
    try
    {
        BoundsCheckRHS(sliceArg);   // make client do this
    }
    catch (hwMathException& except)
    {
        except.Status().ResetArgs();
        throw;
    }

    // verify rhs matrix dimensions
    // only one non-colon is allowed
    int nonColonArg = -1;

    for (int i = 0; i < m_dim.size(); ++i)
    {
        if (!sliceArg[i].IsColon())
        {
            if (nonColonArg == -1)
                nonColonArg = i;
            else
                throw hwMathException(HW_MATH_ERR_SLICE_NUMCOLON);
        }
    }

    if (nonColonArg == -1)
    {
        MakeEmpty();
        return;
    }

    // delete dimension nonColonArg
    std::vector<int> newDim(m_dim);

    if (sliceArg[nonColonArg].IsScalar())
        newDim[nonColonArg] -= 1;
    else if (sliceArg[nonColonArg].IsVector())
    {
        // need to remove any repeated elements in sliceArg[nonColonArg].Vector()
        std::vector<int> tempVec = GetUniqueVec(sliceArg[nonColonArg].Vector());
        newDim[nonColonArg] -= (int) tempVec.size();
    }

    // create new slice arguments
    std::vector<hwSliceArg> newSliceIndex(sliceArg);
    std::vector<int> newVec(newDim[nonColonArg]);
    std::vector<int> tempVec;

    if (sliceArg[nonColonArg].IsScalar())
    {
        for (int i = 0; i < newVec.size(); ++i)
        {
            if (i < sliceArg[nonColonArg].Scalar())
                newVec[i] = i;
            else
                newVec[i] = i+1;
        }
    }
    else if (sliceArg[nonColonArg].IsVector())
    {
        tempVec = GetUniqueVec(sliceArg[nonColonArg].Vector());
        int skip = 0;

        for (int i = 0; i < newVec.size(); ++i)
        {
            if (skip != tempVec.size())
            {
                if (i+skip < tempVec[skip])
                {
                    newVec[i] = i + skip;
                }
                else
                {
                    --i;
                    ++skip;
                }
            }
            else
                newVec[i] = i + skip;
        }
    }

    newSliceIndex.erase(newSliceIndex.begin()+nonColonArg);
    newSliceIndex.insert(newSliceIndex.begin()+nonColonArg, newVec);

    // get matrix with reduced dimensions
    hwTMatrixN<T1, T2> reducedMatrix(newDim, REAL);

    try
    {
        SliceRHS(newSliceIndex, reducedMatrix);
    }
    catch (hwMathException& except)
    {
        except.Status().ResetArgs();
        throw;
    }

    // transfer reduced matrix to *this
    Transfer(reducedMatrix);
}

//! Write a contiguous matrix to the calling object, as if the
//! calling object is on the left hand side of an equals sign
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::CopyMatrixLHS(const hwTMatrixN<T1, T2>& rhsMatrix)
{
    // determine the best approach, depending on whether the slicing involves
    // contiguous or equally spaced memory.
    // 0. leading contiguous dimensions
    // 1. largest equally spaced dimension
    int S_case;
    int keyDim[2]{ -1, -1 };
    int keySize[2]{ 1,  1 };

    for (int i = 0; i < rhsMatrix.m_dim.size(); ++i)
    {
        if (i == keyDim[0] + 1 && m_dim[i] == rhsMatrix.m_dim[i])
        {
            keyDim[0] = i;
            keySize[0] *= rhsMatrix.m_dim[i];
        }

        if (m_dim[i] >= keySize[1])
        {
            keyDim[1] = i;
            keySize[1] = rhsMatrix.m_dim[i];
        }
    }

    if (keySize[0] >= keySize[1])
        S_case = 0;
    else
        S_case = 1;

    // simulate nested loops to iterate over the rhsMatrix elements
    // in order of contiguous memory location, copying blocks where possible
    m_lhsMatrixIndex.clear();
    m_rhsMatrixIndex.clear();
    m_lhsMatrixIndex.resize(m_dim.size());
    m_rhsMatrixIndex.resize(rhsMatrix.m_dim.size());

    int size = rhsMatrix.Size();
    int startIndx = (S_case == 0) ? keyDim[0] + 1 : 0;
    int strideLHS = (S_case == 0) ? 1 : Stride(keyDim[S_case]);
    int strideRHS = (S_case == 0) ? 1 : rhsMatrix.Stride(keyDim[S_case]);

    for (int i = 0; i < size; ++i)
    {
        int startLHS = Index(m_lhsMatrixIndex);
        int startRHS = rhsMatrix.Index(m_rhsMatrixIndex);

        if (IsReal())
            CopyData(m_real + startLHS, strideLHS, rhsMatrix.m_real + startRHS, strideRHS, keySize[S_case]);
        else
            CopyData(m_complex + startLHS, strideLHS, rhsMatrix.m_complex + startRHS, strideRHS, keySize[S_case]);

        i += keySize[S_case] - 1;

        // advance rhs matrix indices
        for (int j = startIndx; j < rhsMatrix.m_dim.size(); ++j)
        {
            if (j == keyDim[S_case])
            {
                continue;
            }

            // increment index j if possible
            if (m_rhsMatrixIndex[j] < (int)rhsMatrix.m_dim[j] - 1)
            {
                ++m_lhsMatrixIndex[j];
                ++m_rhsMatrixIndex[j];
                break;
            }

            // index j is maxed out, so reset and continue to j+1
            m_lhsMatrixIndex[j] = 0;
            m_rhsMatrixIndex[j] = 0;
        }
    }
}

//! Transfer contents from another object
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Transfer(hwTMatrixN<T1, T2>& source)
{
    Deallocate();
    m_bits = source.m_bits;
    m_dim = source.m_dim;
    m_size = source.m_size;
    m_capacity = source.m_capacity;
    m_real_memory = source.m_real_memory;
    m_complex_memory = source.m_complex_memory;
    m_real = source.m_real;
    m_complex = source.m_complex;

    source.m_real_memory = nullptr;
    source.m_complex_memory = nullptr;
    source.m_real = nullptr;
    source.m_complex = nullptr;

    source.MakeEmpty();
}

