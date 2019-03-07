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
#include <hwMathException.h>
#include <hwSliceArg.h>
#include <tmpl/hwTComplex.h>

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
    : m_pos(-1), m_refCount(1)
{
    m_bits.ownData = 0;
    MakeEmpty();
}

//! Construct dimensioned matrix
template<typename T1, typename T2>
hwTMatrixN<T1, T2>::hwTMatrixN(const std::vector<int>& dim, DataType dataType)
    :  m_pos(-1), m_refCount(1), m_capacity(0)
{
    SetDimensions(dim);
    m_bits.ownData = 1;
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
    : m_refCount(1)
{
    m_bits = source.m_bits;
    m_dim = source.m_dim;
    m_size = source.m_size;
    m_pos = source.m_pos;

    if (m_bits.ownData)
    {
        m_real = nullptr;
        m_complex = nullptr;
        m_capacity = source.m_capacity;

        try
        {
            Copy(source);
        }
        catch (hwMathException& except)
        {
            except.Status().ResetArgs();
            throw;
        }
    }
    else
    {
        m_real = source.m_real;
        m_complex = source.m_complex;
        m_capacity = 0;
    }
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

    m_bits = rhs.m_bits;
    m_dim = rhs.m_dim;
    // m_size = rhs.m_size;
    m_real = nullptr;
    m_complex = nullptr;
    m_pos = rhs.m_pos;
    m_capacity = rhs.m_capacity;

    try
    {
        Copy(rhs);
    }
    catch (hwMathException& except)
    {
        except.Status().ResetArgs();
        throw;
    }

    return *this;
}

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
void hwTMatrixN<T1, T2>::ConvertNDto2D(hwTMatrix<T1, T2>& target) const
{
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

    if (IsReal())
        target.Dimension(m, n, hwTMatrix<T1, T2>::REAL);
    else
        target.Dimension(m, n, hwTMatrix<T1, T2>::COMPLEX);

    if (m_size > 0)
    {
        if (m_real)
            CopyData(target.GetRealData(), m_size, m_real, m_size);
        else if (m_complex)
            CopyData(target.GetComplexData(), m_size, m_complex, m_size);
    }
}

//! Convert hwTMatrix to hwTMatrixN
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Convert2DtoND(const hwTMatrix<T1, T2>& source)
{
    if (!IsEmpty())
        MakeEmpty();

    m_dim[0] = source.M();
    m_dim[1] = source.N();

    if (source.IsReal())
        Dimension(m_dim, REAL);
    else
        Dimension(m_dim, COMPLEX);

    if (m_size > 0)
    {
        if (m_real)
            CopyData(m_real, m_size, source.GetRealData(), m_size);
        else if (m_complex)
            CopyData(m_complex, m_size, source.GetComplexData(), m_size);
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
        int i;
        int effDimSize = _min((int) m_dim.size(), (int) dim.size());
        bool sizeChange = false;

        for (i = 0; i < effDimSize; ++i)
        {
            if (m_dim[i] != dim[i])
            {
                sizeChange = true;
                break;
            }
        }

        if (!sizeChange)
        {
            for (i = effDimSize; i < m_dim.size(); ++i)
            {
                if (m_dim[i] != 1)
                {
                    sizeChange = true;
                    break;
                }
            }

            for (i = effDimSize; i < dim.size(); ++i)
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

//! Change the dimensions of a matrix while maintaining the same number of elements
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Reshape(const std::vector<int>& dim)
{
    SetDimensions(dim);

    // check new size
    int size = 1;
    int emptyDim = -1;
    size_t maxSize = std::numeric_limits<int>::max();
    
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
                    if (m_rhsMatrixIndex[k] < (int) source.m_dim[k]-1)
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
                    if (m_rhsMatrixIndex[k] < (int) source.m_dim[k]-1)
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
    int i;
    int numDims = static_cast<int> (permuteVec.size());
    bool* hasDim = new bool[numDims];

    for (i = 0; i < numDims; ++i)
        hasDim[i] = false;

    for (i = 0; i < numDims; ++i)
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

    if (i != numDims)
        throw hwMathException(HW_MATH_ERR_PERMVEC3);
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

//! Bounds check - lower only, since LHS ops allow resizing
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::BoundsCheckLHS(const std::vector<int>& indexVec) const
{
    int i;
    int numDim = (int) m_dim.size();

    for (i = 0; i < numDim; ++i)
    {
        if (indexVec[i] < 0)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);
    }

    for (i = numDim; i < indexVec.size(); ++i)
    {
        if (indexVec[i] < 0)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);
    }
}

//! Bounds check - upper and lower
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::BoundsCheckRHS(const std::vector<int>& indexVec) const
{
    int i;
    int min = _min((int) m_dim.size(), (int) indexVec.size());

    for (i = 0; i < min; ++i)
    {
        if (indexVec[i] < 0)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);

        if (indexVec[i] > m_dim[i]-1)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);
    }

    for (i = min; i < indexVec.size(); ++i)
    {
        if (indexVec[i] != 0)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX, 1);
    }
}

//! Bounds check - lower only, since LHS ops allow resizing
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::BoundsCheckLHS(const std::vector<hwSliceArg>& sliceArg, int numSlices) const
{
    if (numSlices == 0)
        numSlices = RelevantNumberOfSlices(sliceArg, true);     // Evaluator support - deal with it later

    int i;
    int common = _min((int) m_dim.size(), numSlices);

    for (i = 0; i < common; ++i)
    {
        if (sliceArg[i].IsScalar())
        {
            if (sliceArg[i].Scalar() < 0)
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
            }
        }
    }

    for (i = common; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar())
        {
            if (sliceArg[i].Scalar() < 0)
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
            }
        }
    }
}

//! Bounds check - upper and lower
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::BoundsCheckRHS(const std::vector<hwSliceArg>& sliceArg, int numSlices) const
{
    if (numSlices == 0)
        numSlices = RelevantNumberOfSlices(sliceArg, true);

    int i;
    int common = _min((int) m_dim.size(), numSlices);

    for (i = 0; i < common; ++i)
    {
        if (sliceArg[i].IsScalar())
        {
            if (sliceArg[i].Scalar() < 0)
                throw hwMathException(HW_MATH_ERR_SLICE_INDEX, 1);

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

    for (i = common; i < numSlices; ++i)
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

//! Need to grow a matrix, for use with LHS indexing
template<typename T1, typename T2>
bool hwTMatrixN<T1, T2>::NeedToGrowLHS(const std::vector<int>& indexVec)
{
    int i;

    for (i = 0; i < m_dim.size(); ++i)
    {
        if (indexVec[i] > m_dim[i] - 1)
            return true;
    }

    for (i = (int) m_dim.size(); i < indexVec.size(); ++i)
    {
        if (indexVec[i] > 0)
            return true;
    }

    return false;
}

//! Need to grow a matrix, for use with LHS slicing
template<typename T1, typename T2>
bool hwTMatrixN<T1, T2>::NeedToGrowLHS(const std::vector<hwSliceArg>& sliceArg, int numSlices) const
{
    int i;
    int common = _min(numSlices, (int) m_dim.size());

    for (i = 0; i < common; ++i)
    {
        if (sliceArg[i].IsScalar())
        {
            if (sliceArg[i].Scalar() > m_dim[i] - 1)
                return true;
        }
        else if (sliceArg[i].IsColon())
        {
            if (m_dim[i] == 0)
                return true;
        }
        else if (sliceArg[i].IsVector())
        {
            for (int j = 0; j < sliceArg[i].Vector().size(); ++j)
            {
                if (sliceArg[i].Vector()[j] > m_dim[i] - 1)
                    return true;
            }
        }
    }

    for (i = common; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar())
        {
            if (sliceArg[i].Scalar() > -1)
                return true;
        }
        else if (sliceArg[i].IsColon())
        {
            return true;
        }
        else if (sliceArg[i].IsVector())
        {
            for (int j = 0; j < sliceArg[i].Vector().size(); ++j)
            {
                if (sliceArg[i].Vector()[j] > -1)
                    return true;
            }
        }
    }

    return false;
}

//! Grow a matrix
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::GrowLHSMatrix(const std::vector<int>& indexVec)
{
    // create matrix with expanded dimensions
    int numDims = (int) m_dim.size();

    if (indexVec.size() != numDims)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    std::vector<int> newDim(numDims);

    for (int i = 0; i < numDims; ++i)
    {
        if (indexVec[i] > m_dim[i] - 1)
            newDim[i] = indexVec[i] + 1;
        else
            newDim[i] = m_dim[i];
    }

    // this section should be in SetCapacity(), but getting it there
    // requires rethinking the allocation/copy functionality.
    // simply change dimensions if capacity is sufficient, and only
    // zero or singleton dimensions beyond the first change are allowed
    // and can also be changed (to ensure a simple contiguous copy).
    int firstChangedDim = -1;
    bool capacityIsOK = true;

    for (int i = 0; i < numDims; ++i)
    {
        if (newDim[i] != m_dim[i])  // find first changed dimension
        {
            firstChangedDim = i;
            break;
        }
    }

    if (firstChangedDim != -1)
    {
        for (int i = firstChangedDim+1; i < numDims; ++i)
        {
            if (m_dim[i] > 1) // zero or singleton dimensions only
            {
                capacityIsOK = false;
                break;
            }
        }
    }

    if (newDim == m_dim)
        return;

    std::vector<int> oldDim = m_dim;
    m_dim = newDim;
    ComputeSize();

    if (m_size <= m_capacity)
    {
        if (capacityIsOK && firstChangedDim > -1)
        {
            // set newly utilized elements to zero
            std::vector<hwSliceArg> sliceArgs;

            for (int i = 0; i < firstChangedDim; ++i)
                sliceArgs.push_back(hwSliceArg());

            for (int i = firstChangedDim; i < numDims; ++i)
            {
                if (m_dim[i] == 1)
                {
                    sliceArgs.push_back(0);
                }
                else if (m_dim[i] == oldDim[i] + 1)
                {
                    sliceArgs.push_back(oldDim[i]);
                }
                else
                {
                    std::vector<int> vec;

                    for (int j = oldDim[i]; j < m_dim[i]; ++j)
                        vec.push_back(j);

                    sliceArgs.push_back(vec);
                }
            }

            SliceLHS(sliceArgs, 0.0);
            return;
        }
    }
    else if (m_size < 0)
    {
        throw hwMathException(HW_MATH_ERR_ALLOCFAILED);
    }

    // reallocate because capacity was insufficient
    if (m_real)
    {
        T1* m_real_old = m_real;
        hwTMatrixN<T1, T2> tempMatrix(oldDim, m_real_old, REAL);
        m_real = nullptr;
        Allocate(REAL);

        // copy zeros to the expanded matrix
        SetElements(0.0);

        // copy old data to the expanded matrix
        CopyMatrixLHS(tempMatrix);

        if (m_bits.ownData)
            delete [] m_real_old;
        else
            m_bits.ownData = 1;
    }
    else if (m_complex)
    {
        T2* m_complex_old = m_complex;
        hwTMatrixN<T1, T2> tempMatrix(oldDim, m_complex_old, COMPLEX);
        m_complex = nullptr;
        Allocate(COMPLEX);

        // copy zeros to the expanded matrix
        SetElements(0.0);

        // copy old data to the expanded matrix
        CopyMatrixLHS(tempMatrix);

        if (m_bits.ownData)
            delete [] m_complex_old;
        else
            m_bits.ownData = 1;
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

//! Grow a matrix
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::GrowLHSMatrix(const std::vector<hwSliceArg>& sliceArg, int& numSlices,
                                       const hwTMatrixN<T1, T2>* rhsMatrix, bool rule2)
{
    // create matrix with expanded dimensions
    int i;
    int common = _min(numSlices, (int) m_dim.size());
    int rhsIndex = 0;
    std::vector<int> newDim(m_dim.size());

    for (i = 0; i < common; ++i)
    {
        if (sliceArg[i].IsScalar())
        {
            if (sliceArg[i].Scalar() > m_dim[i] - 1)
                newDim[i] = sliceArg[i].Scalar() + 1;
            else
                newDim[i] = m_dim[i];
        }
        else if (sliceArg[i].IsColon())
        {
            if (rhsMatrix)
            {
                // get next non-singleton RHS dimension
                while (rule2 && rhsIndex < rhsMatrix->m_dim.size() && rhsMatrix->m_dim[rhsIndex] == 1)
                    ++rhsIndex;

                if (rhsIndex < rhsMatrix->m_dim.size())
                {
                    if (m_dim[i] == 0)
                        newDim[i] = rhsMatrix->m_dim[rhsIndex];
                    else if (m_dim[i] == rhsMatrix->m_dim[rhsIndex])
                        newDim[i] = m_dim[i];
                    else
                        throw hwMathException(HW_MATH_ERR_SLICE_INDEX);

                    ++rhsIndex;
                }
                else if (m_dim[i] < 2)
                {
                    newDim[i] = 1;
                }
                else
                {
                    throw hwMathException(HW_MATH_ERR_SLICE_INDEX);
                }
            }
            else
            {
                if (m_dim[i] == 0)
                    newDim[i] = 1;
                else
                    newDim[i] = m_dim[i];
            }
        }
        else if (sliceArg[i].IsVector())
        {
            if (rhsMatrix)
            {
                // get next non-singleton RHS dimension
                while (rule2 && rhsIndex < rhsMatrix->m_dim.size() && rhsMatrix->m_dim[rhsIndex] == 1)
                    ++rhsIndex;
            }

            int j;
            int max = m_dim[i] - 1;

            for (j = 0; j < sliceArg[i].Vector().size(); ++j)
            {
                if (sliceArg[i].Vector()[j] > max)
                    max = sliceArg[i].Vector()[j];
            }

            if (max > m_dim[i] - 1)
                newDim[i] = max + 1;
            else
                newDim[i] = m_dim[i];

            if (rhsMatrix)
                ++rhsIndex;
        }
    }

    for (i = common; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar())
        {
            newDim.push_back(sliceArg[i].Scalar() + 1);
        }
        else if (sliceArg[i].IsColon())
        {
            if (rhsMatrix)
            {
                // get next non-singleton RHS dimension
                while (rule2 && rhsIndex < rhsMatrix->m_dim.size() && rhsMatrix->m_dim[rhsIndex] == 1)
                    ++rhsIndex;

                if (rhsIndex < rhsMatrix->m_dim.size())
                    newDim.push_back(rhsMatrix->m_dim[rhsIndex]);
                else
                    newDim.push_back(1);
            }
            else
            {
                newDim.push_back(1);
            }

            if (rhsMatrix)
                ++rhsIndex;
        }
        else if (sliceArg[i].IsVector())
        {
            if (rhsMatrix)
            {
                // get next non-singleton RHS dimension
                while (rule2 && rhsIndex < rhsMatrix->m_dim.size() && rhsMatrix->m_dim[rhsIndex] == 1)
                    ++rhsIndex;
            }

            int j;
            int max = -1;

            for (j = 0; j < sliceArg[i].Vector().size(); ++j)
            {
                if (sliceArg[i].Vector()[j] > max)
                    max = sliceArg[i].Vector()[j];
            }

            newDim.push_back(max + 1);

            if (rhsMatrix)
                ++rhsIndex;
        }
    }

    for (i = numSlices; i < m_dim.size(); ++i)
        newDim[i] = m_dim[i];

    // discard any singleton dimensions that may have been created
    // due to excess colon and "1" slices
    while (newDim.size() > 2 && newDim.back() == 1)
    {
        newDim.pop_back();
        numSlices--;
    }

    // this section should be in SetCapacity(), but getting it there
    // requires rethinking the allocation/copy functionality.
    // simply change dimensions if capacity is sufficient, and only
    // zero or singleton dimensions beyond the first change are allowed
    // and can also be changed (to ensure a simple contiguous copy).
    int numDims = (int) m_dim.size();
    int firstChangedDim = -1;
    bool capacityIsOK = true;

    for (int i = 0; i < numDims; ++i)
    {
        if (newDim[i] != m_dim[i])  // find first changed dimension
        {
            firstChangedDim = i;
            break;
        }
    }

    if (firstChangedDim != -1)
    {
        for (int i = firstChangedDim+1; i < numDims; ++i)
        {
            if (m_dim[i] > 1) // zero or singleton dimensions only
            {
                capacityIsOK = false;
                break;
            }
        }
    }

    if (newDim == m_dim)
        return;

    std::vector<int> oldDim = m_dim;
    m_dim = newDim;
    ComputeSize();

    if (m_size <= m_capacity)
    {
        if (capacityIsOK && firstChangedDim > -1)
        {
            // set newly utilized elements to zero
            std::vector<hwSliceArg> sliceArgs;

            for (int i = 0; i < firstChangedDim; ++i)
                sliceArgs.push_back(hwSliceArg());

            for (int i = firstChangedDim; i < numDims; ++i)
            {
                if (m_dim[i] == 1)
                {
                    sliceArgs.push_back(0);
                }
                else if (m_dim[i] == oldDim[i] + 1)
                {
                    sliceArgs.push_back(oldDim[i]);
                }
                else
                {
                    std::vector<int> vec(m_dim[i]-oldDim[i]);

                    for (int j = oldDim[i]; j < m_dim[i]; ++j)
                        vec.push_back(j);

                    sliceArgs.push_back(vec);
                }
            }

            SliceLHS(sliceArgs, 0.0);
            return;
        }
    }
    else if (m_size < 0)
    {
        throw hwMathException(HW_MATH_ERR_ALLOCFAILED);
    }

    // reallocate because capacity was insufficient
    if (m_bits.realData)
    {
        T1* m_real_old = m_real;
        hwTMatrixN<T1, T2> tempMatrix(oldDim, m_real_old, REAL);
        m_real = nullptr;
        Allocate(REAL);

        // copy zeros to the expanded matrix
        SetElements(0.0);

        // copy old data to the expanded matrix
        CopyMatrixLHS(tempMatrix);

        if (m_bits.ownData)
            delete [] m_real_old;
        else
            m_bits.ownData = 1;
    }
    else if (m_complex)
    {
        T2* m_complex_old = m_complex;
        hwTMatrixN<T1, T2> tempMatrix(oldDim, m_complex_old, COMPLEX);
        m_complex = nullptr;
        Allocate(COMPLEX);

        // copy zeros to the expanded matrix
        SetElements(0.0);

        // copy old data to the expanded matrix
        CopyMatrixLHS(tempMatrix);

        if (m_bits.ownData)
            delete [] m_complex_old;
        else
            m_bits.ownData = 1;
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

// ****************************************************
//         Access Functions for Real Elements
// ****************************************************

//! Return a reference to the real data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
T1& hwTMatrixN<T1, T2>::operator()(const std::vector<int>& indexVec)
{
    m_pos = Index(indexVec);

    return m_real[m_pos];
}

//! Return a const reference to the real data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
const T1& hwTMatrixN<T1, T2>::operator()(const std::vector<int>& indexVec) const
{
    m_pos = Index(indexVec);

    return m_real[m_pos];
}

//! Set every matrix element to a specified real value
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::SetElements(T1 real)
{
    int i;
    int size = Size();

    if (m_real)
    {
        for (i = 0; i < size; ++i)
            m_real[i] = real;
    }

    // the function does not reallocate a complex matrix as real
    // the user must manage that if desired
    if (m_complex)
    {
        for (i = 0; i < size; ++i)
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
    int numDim = _min((int) m_dim.size(), (int) indexVec.size());

    m_pos = indexVec[numDim-1];

    for (int i = numDim-2; i > -1; --i)
        m_pos = m_pos * m_dim[i] + indexVec[i];

    return m_complex[m_pos];
}

//! Return a const reference to the complex data element at the specified indices
//! Client is responsible for bound checking
template<typename T1, typename T2>
const T2& hwTMatrixN<T1, T2>::z(const std::vector<int>& indexVec) const
{
    int numDim = _min((int) m_dim.size(), (int) indexVec.size());

    m_pos = indexVec[numDim-1];

    for (int i = numDim-2; i > -1; --i)
        m_pos = m_pos * m_dim[i] + indexVec[i];

    return m_complex[m_pos];
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

    int i;
    int size = Size();

    for (i = 0; i < size; ++i)
    {
        m_complex[i].Real() = m_real[i];
        m_complex[i].Imag() = (T1) 0;
    }

    if (m_real)
    {
        delete [] m_real;
        m_real = nullptr;
    }

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

    int i, j;

    // ignore high singleton dimensions
    int numSlices = RelevantNumberOfSlices(sliceArg, false);

    if (numSlices > m_dim.size())
        throw hwMathException(HW_MATH_ERR_SLICE_INDEX);

    // perform implicit RHS reshape with reduced slices 
    if (numSlices < m_dim.size())
    {
        std::vector<int> newDim(numSlices);

        for (i = 0; i < numSlices-1; ++i)
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

    for (i = 0; i < lhsNumDims; ++i)
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

    if (m_size < size)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    if (size == 0)
        return;

    // determine slice at which discontiguity occurs
    int discontiguity = numSlices;

    for (j = 0; j < numSlices; ++j)
    {
        if (sliceArg[j].IsScalar())
        {
            if (m_dim[j] != 1)
            {
                discontiguity = j;
                break;
            }
        }
        else if (sliceArg[j].IsColon())
        {
        }
        else // if (sliceArg[j].IsVector())
        {
            discontiguity = j;
            break;
        }
    }

    // set the rhsMatrix indices to the first index in each slice
    m_rhsMatrixIndex.resize(numSlices);

    for (i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar())
            m_rhsMatrixIndex[i] = sliceArg[i].Scalar();
        else if (sliceArg[i].IsColon())
            m_rhsMatrixIndex[i] = 0;
        else if (sliceArg[i].IsVector())
            m_rhsMatrixIndex[i] = sliceArg[i].Vector()[0];
    }

    // simulate nested loops to iterate over the lhsMatrix elements
    // in order of contiguous memory location, copying blocks where possible
    m_lhsMatrixIndex.clear();
    m_lhsMatrixIndex.resize(numSlices);

    for (i = 0; i < size; ++i)
    {
        // copy data up to discontguity
        if (discontiguity == 0)
        {
            if (lhsMatrix.IsReal())
                lhsMatrix(i) = (*this)(m_rhsMatrixIndex);
            else
                lhsMatrix.z(i) = this->z(m_rhsMatrixIndex);
        }
        else
        {
            SetMemoryPosition(m_rhsMatrixIndex);    // update m_pos
            CopyBlockRHS(i, discontiguity, lhsMatrix);
            --i;
        }

        // advance lhs matrix indices
        for (j = discontiguity; j < numSlices; ++j)
        {
            if (sliceArg[j].IsScalar())
            {
                continue;
            }
            else if (sliceArg[j].IsColon())
            {
                // increment index j if possible
                if (m_rhsMatrixIndex[j] < (int) m_dim[j]-1)
                {
                    ++m_rhsMatrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                m_rhsMatrixIndex[j] = 0;
            }
            else // if (sliceArg[j].IsVector())
            {
                // increment index j if possible
                if (m_lhsMatrixIndex[j] < (int) sliceArg[j].Vector().size()-1)
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
            if (numSlices == 1 && sliceArg[0].IsColon())    // handle a(:)
            {
                newDim.push_back(1);
                std::vector<hwSliceArg> newSliceArg;
                newSliceArg.push_back(hwSliceArg());
                newSliceArg.push_back(0);
                reshaped.Reshape(newDim);
                reshaped.SliceLHS(newSliceArg, rhsMatrix);
                return;
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

    // check LHS bounds and data type, modify as needed
    try
    {
        BoundsCheckLHS(sliceArg, numSlices);
    }
    catch (hwMathException& except)
    {
        except.Status().ResetArgs();
        throw;
    }

    // count the LHS slicing colons pick the initial slicing rule
    // Rule 1: map the RHS dimensions to the LHS colons
    // Rule 2: map the non-singleton RHS dimensions to the LHS colons
    // Rule 3: attempt Rule2 if Rule 1 fails
    bool rule2;
    int numColons = 0;

    for (int i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsColon())
            ++numColons;
    }

    if (numColons >= rhsMatrix.m_dim.size())
    {
        rule2 = false;
    }
    else // if (numColons < rhsMatrix.m_dim.size())
    {
        rule2 = true;
    }

    // populate slices
    if (m_dim.size() == 2 && m_dim[0] == 0 && m_dim[1] == 0)
    {
        // dimension the matrix
        int rhsIndex = 0;
        std::vector<int> newDim(numSlices);

        for (int i = 0; i < numSlices; ++i)
        {
            if (sliceArg[i].IsScalar())
            {
                newDim[i] = sliceArg[i].Scalar() + 1;

                if (!rule2)
                    rule2 = true;
            }
            else if (sliceArg[i].IsColon())
            {
                // get next non-singleton RHS dimension
                while (rule2 && rhsIndex < rhsMatrix.m_dim.size() && rhsMatrix.m_dim[rhsIndex] == 1)
                    ++rhsIndex;

                if (rhsIndex < rhsMatrix.m_dim.size())
                    newDim[i] = rhsMatrix.m_dim[rhsIndex];
                else
                    newDim[i] = 1;

                ++rhsIndex;
            }
            else if (sliceArg[i].IsVector())
            {
                // get next non-singleton RHS dimension
                while (rule2 && rhsIndex < rhsMatrix.m_dim.size() && rhsMatrix.m_dim[rhsIndex] == 1)
                    ++rhsIndex;

                int max = -1;

                for (int j = 0; j < sliceArg[i].Vector().size(); ++j)
                {
                    if (sliceArg[i].Vector()[j] > max)
                        max = sliceArg[i].Vector()[j];
                }

                newDim[i] = max + 1;
                ++rhsIndex;

                if (!rule2)
                    rule2 = true;
            }
        }

        // discard any singleton dimensions that may have been created
        // due to excess colon and "1" slices
        while (newDim.size() > 2 && newDim.back() == 1)
        {
            newDim.pop_back();
            numSlices--;
        }

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
        if (NeedToGrowLHS(sliceArg, numSlices))
        {
            try
            {
                GrowLHSMatrix(sliceArg, numSlices, &rhsMatrix, rule2);
            }
            catch (hwMathException& except)
            {
                except.Status().ResetArgs();
                throw;
            }
        }

        if (!rhsMatrix.IsReal())
        {
            try
            {
                if (IsReal())
                    MakeComplex();
            }
            catch (hwMathException& except)
            {
                except.Status().ResetArgs();
                throw;
            }
        }

        // verify that slice op is well defined
        int rhsIndex = 0;

        for (int i = 0; i < numSlices; ++i)
        {
            if (sliceArg[i].IsScalar())
            {
                if (rhsIndex < rhsMatrix.m_dim.size())
                {
                    if (rhsMatrix.m_dim[rhsIndex] == 1)
                    {
                        ++rhsIndex;
                        continue;
                    }

                    if (!rule2)
                    {
                        rule2 = true;
                    }
                }
            }
            else if (sliceArg[i].IsColon())
            {
                if (rule2)
                {
                    // get next non-singleton RHS dimension
                    while (rhsIndex < rhsMatrix.m_dim.size() && rhsMatrix.m_dim[rhsIndex] == 1)
                        ++rhsIndex;
                }

                if (rhsIndex < rhsMatrix.m_dim.size())
                {
                    if (rhsMatrix.m_dim[rhsIndex] > m_dim[i])
                    {
                        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
                    }

                    ++rhsIndex;
                }
                else
                {
                    if (m_dim[i] != 1)
                        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
                }
            }
            else if (sliceArg[i].IsVector())
            {
                if (rhsMatrix.m_dim[rhsIndex] == sliceArg[i].Vector().size())
                {
                    if (rule2)
                    {
                        // get next non-singleton RHS dimension
                        while (rhsIndex < rhsMatrix.m_dim.size() && rhsMatrix.m_dim[rhsIndex] == 1)
                            ++rhsIndex;
                    }

                    if (rhsIndex < rhsMatrix.m_dim.size())
                    {
                        if (sliceArg[i].Vector().size() != rhsMatrix.m_dim[rhsIndex])
                            throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
                    }
                    else
                    {
                        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
                    }

                    ++rhsIndex;
                    continue;
                }

                if (!rule2)
                {
                    rule2 = true;
                }
            }
        }

        // ensure that any remaining RHS higher dimensions are singletons
        while (rhsIndex < rhsMatrix.m_dim.size() && rhsMatrix.m_dim[rhsIndex] == 1)
            ++rhsIndex;

        if (rhsIndex != rhsMatrix.m_dim.size())
            throw hwMathException(HW_MATH_ERR_INVALIDINPUT);
    }

    // determine slice at which discontiguity occurs
    int discontiguity = numSlices;

    for (int j = 0; j < numSlices; ++j)
    {
        if (sliceArg[j].IsScalar())
        {
            if (m_dim[j] != 1)
            {
                discontiguity = j;
                break;
            }
        }
        else if (sliceArg[j].IsColon())
        {
        }
        else // if (sliceArg[j].IsVector())
        {
            discontiguity = j;
            break;
        }
    }

    // set the lhsMatrix indices to the first index in each slice
    m_lhsMatrixIndex.resize(numSlices);

    for (int i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar())
            m_lhsMatrixIndex[i] = sliceArg[i].Scalar();
        else if (sliceArg[i].IsColon())
            m_lhsMatrixIndex[i] = 0;
        else if (sliceArg[i].IsVector())
            m_lhsMatrixIndex[i] = sliceArg[i].Vector()[0];
    }

    // simulate nested loops to iterate over the rhsMatrix elements
    // in order of contiguous memory location, copying blocks where possible
    int size = rhsMatrix.Size();

    m_rhsMatrixIndex.clear();
    m_rhsMatrixIndex.resize(numSlices);

    if (m_size < size)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    if (size == 0)
        return;

    if (!IsReal() && rhsMatrix.IsReal())
        discontiguity = 0;   // need dcopy to manage this

    for (int i = 0; i < size; ++i)
    {
        // copy data up to discontguity
        if (discontiguity == 0)
        {
            if (IsReal())
                (*this)(m_lhsMatrixIndex) = rhsMatrix(i);
            else if (rhsMatrix.IsReal())
                this->z(m_lhsMatrixIndex) = rhsMatrix(i);
            else
                this->z(m_lhsMatrixIndex) = rhsMatrix.z(i);
        }
        else
        {
            SetMemoryPosition(m_lhsMatrixIndex);    // update m_pos
            CopyBlockLHS(i, discontiguity, rhsMatrix);
            --i;
        }

        // advance rhs matrix indices
        for (int j = discontiguity; j < numSlices; ++j)
        {
            if (sliceArg[j].IsScalar())
            {
                continue;
            }
            else if (sliceArg[j].IsColon())
            {
                // increment index j if possible
                if (m_lhsMatrixIndex[j] < (int) m_dim[j]-1)
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
                if (m_rhsMatrixIndex[j] < (int) sliceArg[j].Vector().size()-1)
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
void hwTMatrixN<T1, T2>::SliceLHS(const std::vector<hwSliceArg>& sliceArg, T1 real)
{
    if (sliceArg.size() == 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    int i, j;

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

        if (index < 0)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

        if (index > m_size - 1)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

        (*this)(index) = real;

        return;
    }

    // perform implicit LHS reshape with reduced slices 
    if (numSlices < m_dim.size())
    {
        std::vector<int> newDim(numSlices);

        for (i = 0; i < numSlices-1; ++i)
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
                reshaped.SliceLHS(newSliceArg, real);
                return;
            }

            reshaped.Reshape(newDim);
            reshaped.SliceLHS(sliceArg, real);

            if (reshaped.m_real != m_real)  // reallocation occurred
            {
                delete [] m_real;
                m_real = reshaped.m_real;
                reshaped.m_bits.ownData = 0;

                for (i = 0; i < numSlices-1; ++i)
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

    // check LHS bounds, modify as needed
    try
    {
        BoundsCheckLHS(sliceArg, numSlices);
    }
    catch (hwMathException& except)
    {
        except.Status().ResetArgs();
        throw;
    }

    if (m_dim.size() == 2 && m_dim[0] == 0 && m_dim[1] == 0)
    {
        std::vector<int> newDim(numSlices);

        for (int i = 0; i < numSlices; ++i)
        {
            if (sliceArg[i].IsScalar())
            {
                newDim[i] = sliceArg[i].Scalar() + 1;
            }
            else if (sliceArg[i].IsColon())
            {
                newDim[i] = 1;
            }
            else if (sliceArg[i].IsVector())
            {
                int max = -1;

                for (int j = 0; j < sliceArg[i].Vector().size(); ++j)
                {
                    if (sliceArg[i].Vector()[j] > max)
                        max = sliceArg[i].Vector()[j];
                }

                newDim[i] = max + 1;
            }
        }

        m_dim = newDim;
        Allocate(REAL);
        SetElements(0.0);
    }
    else if (NeedToGrowLHS(sliceArg, numSlices))
    {
        try
        {
            GrowLHSMatrix(sliceArg, numSlices);
        }
        catch (hwMathException& except)
        {
            except.Status().ResetArgs();
            throw;
        }
    }

    // determine slice at which discontiguity occurs
    int discontiguity = numSlices;

    for (j = 0; j < numSlices; ++j)
    {
        if (sliceArg[j].IsScalar())
        {
            if (m_dim[j] != 1)
            {
                discontiguity = j;
                break;
            }
        }
        else if (sliceArg[j].IsColon())
        {
        }
        else // if (sliceArg[j].IsVector())
        {
            discontiguity = j;
            break;
        }
    }

    // set the lhsMatrix indices to the first index
    // in each slice, then assign the first rhs matrix element
    m_lhsMatrixIndex.resize(numSlices);

    for (i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar())
            m_lhsMatrixIndex[i] = sliceArg[i].Scalar();
        else if (sliceArg[i].IsColon())
            m_lhsMatrixIndex[i] = 0;
        else if (sliceArg[i].IsVector())
            m_lhsMatrixIndex[i] = sliceArg[i].Vector()[0];
    }

    // simulate nested loops to iterate over the rhsMatrix elements
    // in order of contiguous memory location, copying blocks where possible
    int size = 1;    

    m_rhsMatrixIndex.clear();
    m_rhsMatrixIndex.resize(numSlices);

    for (i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar()) {}
        else if (sliceArg[i].IsColon())
            size *= m_dim[i];
        else if (sliceArg[i].IsVector())
            size *= (int) sliceArg[i].Vector().size();
    }

    if (m_size < size)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    for (i = 0; i < size; ++i)
    {
        if (discontiguity == 0)
        {
            if (IsReal())
                (*this)(m_lhsMatrixIndex) = real;
            else
                this->z(m_lhsMatrixIndex) = real;
        }
        else
        {
            SetMemoryPosition(m_lhsMatrixIndex);    // update m_pos
            CopyBlockLHS(i, discontiguity, real);
            --i;
        }

        // advance rhs matrix indices, as if the RHS is real*ones(slice_dims)
        // there is probably a better way to do this
        for (j = discontiguity; j < numSlices; ++j)
        {
            if (sliceArg[j].IsScalar())
            {
                continue;
            }
            else if (sliceArg[j].IsColon())
            {
                // increment index j if possible
                if (m_lhsMatrixIndex[j] < (int) m_dim[j]-1)
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
                int vecSize = (int) sliceArg[j].Vector().size();

                if (m_lhsMatrixIndex[j] < (int) sliceArg[j].Vector()[vecSize-1])
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

    int i, j;

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

        if (index < 0)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

        if (index > m_size - 1)
            throw hwMathException(HW_MATH_ERR_INVALIDINDEX);

        this->z(index) = cmplx;
        return;
    }

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
                    delete [] m_complex;
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

    // check LHS bounds, modify as needed
    try
    {
        BoundsCheckLHS(sliceArg, numSlices);
    }
    catch (hwMathException& except)
    {
        except.Status().ResetArgs();
        throw;
    }

    if (m_dim.size() == 2 && m_dim[0] == 0 && m_dim[1] == 0)
    {
        std::vector<int> newDim(numSlices);

        for (int i = 0; i < numSlices; ++i)
        {
            if (sliceArg[i].IsScalar())
            {
                newDim[i] = sliceArg[i].Scalar() + 1;
            }
            else if (sliceArg[i].IsColon())
            {
                newDim[i] = 1;
            }
            else if (sliceArg[i].IsVector())
            {
                int max = -1;

                for (int j = 0; j < sliceArg[i].Vector().size(); ++j)
                {
                    if (sliceArg[i].Vector()[j] > max)
                        max = sliceArg[i].Vector()[j];
                }

                newDim[i] = max + 1;
            }
        }

        m_dim = newDim;
        Allocate(COMPLEX);
        m_bits.realData = 0;
        SetElements(0.0);
    }
    else if (NeedToGrowLHS(sliceArg, numSlices))
    {
        try
        {
            GrowLHSMatrix(sliceArg, numSlices);
        }
        catch (hwMathException& except)
        {
            except.Status().ResetArgs();
            throw;
        }
    }

    // determine slice at which discontiguity occurs
    int discontiguity = numSlices;

    for (j = 0; j < numSlices; ++j)
    {
        if (sliceArg[j].IsScalar())
        {
            if (m_dim[j] != 1)
            {
                discontiguity = j;
                break;
            }
        }
        else if (sliceArg[j].IsColon())
        {
        }
        else // if (sliceArg[j].IsVector())
        {
            discontiguity = j;
            break;
        }
    }

    // set the lhsMatrix indices to the first index
    // in each slice, then assign the first rhs matrix element
    m_lhsMatrixIndex.resize(numSlices);

    for (i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar())
            m_lhsMatrixIndex[i] = sliceArg[i].Scalar();
        else if (sliceArg[i].IsColon())
            m_lhsMatrixIndex[i] = 0;
        else if (sliceArg[i].IsVector())
            m_lhsMatrixIndex[i] = sliceArg[i].Vector()[0];
    }

    // simulate nested loops to iterate over the rhsMatrix elements
    // in order of contiguous memory location, copying blocks where possible
    int size = 1;    

    m_rhsMatrixIndex.clear();
    m_rhsMatrixIndex.resize(numSlices);

    for (i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsScalar()) {}
        else if (sliceArg[i].IsColon())
            size *= m_dim[i];
        else if (sliceArg[i].IsVector())
            size *= (int) sliceArg[i].Vector().size();
    }

    if (m_size < size)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

    for (i = 0; i < size; ++i)
    {
        if (discontiguity == 0)
        {
            this->z(m_lhsMatrixIndex) = cmplx;
        }
        else
        {
            SetMemoryPosition(m_lhsMatrixIndex);    // update m_pos
            CopyBlockLHS(i, discontiguity, cmplx);
            --i;
        }

        // advance rhs matrix indices, as if the RHS is cmplx*ones(slice_dims)
        // there is probably a better way to do this
        for (j = discontiguity; j < numSlices; ++j)
        {
            if (sliceArg[j].IsScalar())
            {
                continue;
            }
            else if (sliceArg[j].IsColon())
            {
                // increment index j if possible
                if (m_lhsMatrixIndex[j] < (int) m_dim[j]-1)
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
                int vecSize = (int) sliceArg[j].Vector().size();

                if (m_lhsMatrixIndex[j] < (int) sliceArg[j].Vector()[vecSize-1])
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
    int i;
    int effDimSize = _min((int) m_dim.size(), (int) A.m_dim.size());

    for (i = 0; i < effDimSize; ++i)
    {
        if (m_dim[i] != A.m_dim[i])
            return false;
    }

    for (i = effDimSize; i < m_dim.size(); ++i)
    {
        if (m_dim[i] != 1)
            return false;
    }

    for (i = effDimSize; i < A.m_dim.size(); ++i)
    {
        if (A.m_dim[i] != 1)
            return false;
    }

    if (IsReal())
    {
        if (A.IsReal())
        {
            for (i = 0; i < m_size; ++i)
            {
                if (m_real[i] != A.m_real[i])
                    return false;
            }
        }
        else
        {
            for (i = 0; i < m_size; ++i)
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
            for (i = 0; i < m_size; ++i)
            {
                if (m_complex[i] != A.m_real[i])
                    return false;
            }
        }
        else
        {
            for (i = 0; i < m_size; ++i)
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
        int i;
        int numDims = 2;

        // ignore trailing dimensions greater than 2 of size 1
        for (i = (int) size-1; i > 1; --i)
        {
            if (dim[i] != 1)
            {
                numDims = i+1;
                break;
            }
        }

        m_dim.resize(numDims);

        for (i = 0; i < numDims; ++i)
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

    size_t maxSize = std::numeric_limits<int>::max();
    
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
        if (m_dim[i] >= maxSize / (size_t) m_size)
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
        m_capacity = -1;    // force allocaction failure
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
        if (m_real)
            delete [] m_real;

        if (m_complex)
            delete [] m_complex;
    }
}

//! Set matrix to empty condition
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::MakeEmpty()
{
    Deallocate();
    m_capacity = 0;
    m_real = nullptr;
    m_complex = nullptr;
    m_size = 0;
    m_bits.ownData = 1;
    m_bits.realData = 1;
    m_pos = -1;
    m_dim.resize(2);
    m_dim[0] = 0;
    m_dim[1] = 0;
}

//! Set memory position corresponding to an index vector
template<typename T1, typename T2>
int hwTMatrixN<T1, T2>::SetMemoryPosition(const std::vector<int>& indexVec) const
{
    int numDim = _min((int) m_dim.size(), (int) indexVec.size());

    m_pos = indexVec[numDim-1];

    for (int i = numDim-2; i > -1; --i)
        m_pos = m_pos * m_dim[i] + indexVec[i];

    return m_pos;
}

//! Copy function
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::Copy(const hwTMatrixN<T1, T2>& source)
{
    bool sizeChange = (m_size != source.m_size);

    m_size = source.m_size;

    if (m_bits.ownData)
    {
        if (source.m_real)
        {
            if (!m_real)
            {
                Allocate(REAL);
            }
            else if (sizeChange)
            {
                delete [] m_real;
                m_real = nullptr;
                Allocate(REAL);
            }
        }
        else if (m_real)
        {
            delete [] m_real;
            m_real = nullptr;
        }

        if (source.m_complex)
        {
            if (!m_complex)
            {
                Allocate(COMPLEX);
            }
            else if (sizeChange)
            {
                delete [] m_complex;
                m_complex = nullptr;
                Allocate(COMPLEX);
            }
        }
        else if (m_complex)
        {
            delete [] m_complex;
            m_complex = nullptr;
        }
    }
    else if (sizeChange)
    {
        // copying data is not allowed when the target matrix does not
        // own its data and there is a change of dimension
        MakeEmpty();
        throw hwMathException(HW_MATH_ERR_NOTALLOWED);
    }

    if (m_size > 0)
    {
        if (m_real)
            CopyData(m_real, m_size, source.m_real, m_size);
        else if (m_complex)
            CopyData(m_complex, m_size, source.m_complex, m_size);
    }
}

//! Copy data
#ifdef _BLAS_LAPACK_h
    template<>
    inline void hwTMatrixN<double>::CopyData(void* dest, int arraySize, const void* src, int count) const
    {
        int inc = 1;

        if (IsReal())
            dcopy_((int*) &count, (double*) src, &inc, (double*) dest, &inc);
        else
            zcopy_((int*) &count, (complexD*) src, &inc, (complexD*) dest, &inc);
    }
#else
    #include <memory.h>

    template<>
    inline void hwTMatrixN<double>::CopyData(void* dest, int arraySize, const void* src, int count) const
    {
        if (IsReal())
            memcpy_s(dest, arraySize * sizeof(double), src, count * sizeof(double));
        else
            memcpy_s(dest, arraySize * sizeof(hwTComplex<double>), src, count * sizeof(hwTComplex<double>));
    }
#endif

template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::CopyData(void* dest, int arraySize, const void* src, int count) const
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
    int i;
    int nonColonArg = -1;

    for (i = 0; i < m_dim.size(); ++i)
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
        for (i = 0; i < newVec.size(); ++i)
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

        for (i = 0; i < newVec.size(); ++i)
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

//! Read a contiguous block from the calling object, as if the calling
//! object is being sliced on the right hand side of an equals sign
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::CopyBlockRHS(int& pos, int sliceArg, hwTMatrixN<T1, T2>& lhsMatrix) const
{
    int numVals = 1;
    int numArgs = _min(sliceArg, (int) lhsMatrix.m_dim.size());

    for (int i = 0; i < numArgs; ++i)
        numVals *= lhsMatrix.m_dim[i];

    // the lhsMatrix block begins at location pos
    // the lhsMatrix block begins at location m_pos, which must be set prior to the function call
    // switch memcpy_s to dcopy from BLAS
    if (numVals > 0)
    {
        if (m_real)
            CopyData(lhsMatrix.m_real+pos, numVals, m_real+m_pos, numVals);
        else if (m_complex)
            CopyData(lhsMatrix.m_complex+pos, numVals, m_complex+m_pos, numVals);
    }

    pos += numVals;
}

//! Write a contiguous block to the calling object, as if the calling
//! object is being sliced on the left hand side of an equals sign
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::CopyBlockLHS(int& pos, int sliceArg, const hwTMatrixN<T1, T2>& rhsMatrix)
{
    int numVals = 1;
    int numArgs = _min(sliceArg, (int) m_dim.size());

    for (int i = 0; i < numArgs; ++i)
        numVals *= m_dim[i];

    // the rhsMatrix block begins at location pos
    // the lhsMatrix block begins at location m_pos, which must be set prior to the function call
    if (numVals > 0)
    {
        if (m_real)
            CopyData(m_real+m_pos, numVals, rhsMatrix.m_real+pos, numVals);
        else if (m_complex)
            CopyData(m_complex+m_pos, numVals, rhsMatrix.m_complex+pos, numVals);
    }

    pos += numVals;
}

//! Write a contiguous block to the calling object, as if the calling
//! object is being sliced on the left hand side of an equals sign
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::CopyBlockLHS(int& pos, int sliceArg, T1 real)
{
    int numVals = 1;

    for (int i = 0; i < sliceArg; ++i)
       numVals *= m_dim[i];

    // the lhsMatrix block begins at location m_pos, which must be set prior to the function call
    for (int i = 0; i < numVals; ++i)
       (m_real+m_pos)[i] = real;

    pos += numVals;
}

//! Write a contiguous block to the calling object, as if the calling
//! object is being sliced on the left hand side of an equals sign
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::CopyBlockLHS(int& pos, int sliceArg, const T2& cmplx)
{
    int numVals = 1;

    for (int i = 0; i < sliceArg; ++i)
       numVals *= m_dim[i];

    // the lhsMatrix block begins at location m_pos, which must be set prior to the function call
    for (int i = 0; i < numVals; ++i)
       (m_complex+m_pos)[i] = cmplx;

    pos += numVals;
}

//! Write a contiguous matrix to the calling object, as if the
//! calling object is on the left hand side of an equals sign
template<typename T1, typename T2>
void hwTMatrixN<T1, T2>::CopyMatrixLHS(const hwTMatrixN<T1, T2>& rhsMatrix)
{
    // determine dimension at which discontiguity occurs
    int discontiguity = (int) rhsMatrix.m_dim.size();

    for (int i = 0; i < rhsMatrix.m_dim.size(); ++i)
    {
        if (m_dim[i] != rhsMatrix.m_dim[i])
        {
            discontiguity = i;
            break;
        }
    }

    // simulate nested loops to iterate over the rhsMatrix elements
    // in order of contiguous memory location, copying blocks where possible
    int size = rhsMatrix.Size();

    m_rhsMatrixIndex.clear();
    m_rhsMatrixIndex.resize(m_dim.size());

    for (int i = 0; i < size; ++i)
    {
        if (discontiguity == 0)
        {
            if (IsReal())
                (*this)(m_rhsMatrixIndex) = rhsMatrix(i);
            else if (rhsMatrix.IsReal())
                this->z(m_rhsMatrixIndex) = rhsMatrix(i);
            else
                this->z(m_rhsMatrixIndex) = rhsMatrix.z(i);
        }
        else
        {
            SetMemoryPosition(m_rhsMatrixIndex);    // update m_pos
            CopyBlockLHS(i, discontiguity, rhsMatrix);
            --i;
        }

        // advance rhs matrix indices
        for (int j = discontiguity; j < rhsMatrix.m_dim.size(); ++j)
        {
            // increment index j if possible
            if (m_rhsMatrixIndex[j] < (int) rhsMatrix.m_dim[j]-1)
            {
                ++m_rhsMatrixIndex[j];
                break;
            }

            // index j is maxed out, so reset and continue to j+1
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
    m_real = source.m_real;
    m_complex = source.m_complex;
    m_pos = source.m_pos;

    source.m_real = nullptr;
    source.m_complex = nullptr;

    source.MakeEmpty();
}

