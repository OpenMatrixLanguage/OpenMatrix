#include "hwMatrixN.h"
#include "Currency.h"

template<> inline
void hwTMatrixN<Currency, void*>::SliceLHS(const std::vector<hwSliceArg>& sliceArg, Currency real);

//! Write a contiguous matrix to the calling object, as if the
//! calling object is on the left hand side of an equals sign
template<> inline
void hwTMatrixN<Currency, void*>::CopyMatrixLHS(const hwTMatrixN<Currency, void*>& rhsMatrix)
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

            CopyData(m_real + startLHS, strideLHS, rhsMatrix.m_real + startRHS, strideRHS, keySize[S_case]);
            i += keySize[S_case] - 1;

            // advance rhs matrix indices
            for (int j = startIndx; j < rhsMatrix.m_dim.size(); ++j)
            {
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
}

//! Grow a matrix
template<> inline
void hwTMatrixN<Currency, void*>::GrowLHSMatrix(const std::vector<int>& newDim)
{
	// create matrix with expanded dimensions
	int numDims = (int)m_dim.size();

	if (newDim.size() < numDims)
		throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

	std::vector<int> oldDim = m_dim;
	m_dim = newDim;
	ComputeSize();

	if (m_real)
	{
		char* m_real_mem_old = m_real_memory;
        Currency* m_real_old = m_real;
		hwTMatrixN<Currency, void*> tempMatrix(oldDim, m_real_old, REAL);
		m_real_memory = nullptr;
		m_real = nullptr;
		Allocate(REAL);

		// copy zeros to the expanded matrix
		// SetElements(0.0);

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

//! Write to a matrix slice of the calling object, as if the calling object is being
//! sliced on the left hand side of an equals sign
template<> inline
void hwTMatrixN<Currency, void*>::SliceLHS(const std::vector<hwSliceArg>& sliceArg, Currency real)
{
    if (sliceArg.size() == 0)
        throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

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

        for (int i = 0; i < numSlices - 1; ++i)
            newDim[i] = m_dim[i];

        newDim[numSlices - 1] = -1;
        void* dataPtr = (void*)m_real;

        hwTMatrixN<Currency, void*> reshaped(m_dim, dataPtr, Type());

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
            if (sliceArg[lhsDim].Scalar() < 0)
                throw hwMathException(HW_MATH_ERR_SLICE_INDEX);

            if (newMatrix)
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;
            else if (lhsDim < m_dim.size())
                newDim[lhsDim] = _max(m_dim[lhsDim], static_cast<int> (sliceArg[lhsDim].Scalar()) + 1);
            else
                newDim[lhsDim] = static_cast<int> (sliceArg[lhsDim].Scalar()) + 1;
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
                newDim[lhsDim] = _max(m_dim[lhsDim], static_cast<int> (sliceArg[lhsDim].Scalar() + 1));
            else
                newDim[lhsDim] = static_cast<int> (sliceArg[lhsDim].Scalar() + 1);
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
        // SetElements(0.0);
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
            Currency* dataPtr = m_real + startLHS;

            for (int k = 0; k < keySize[S_case]; ++k)
            {
                *dataPtr = real;
                dataPtr += strideLHS;
            }

            break;
        }
        case 3:
        {
            if (sliceArg[0].IsScalar())
            {
                (*this)(m_lhsMatrixIndex) = real;
            }
            else if (sliceArg[0].IsVector())
            {
                int vecLen = static_cast<int> (sliceArg[0].Vector().size());

                m_lhsMatrixIndex[0] = 0;
                int topOfDim = Index(m_lhsMatrixIndex);

                for (int j = 0; j < vecLen; ++j)
                {
                    int pos = topOfDim + sliceArg[0].Vector()[j];
                    (*this)(pos) = real;
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
//! matrix(slice args) = rhs_matrix
template<> inline
void hwTMatrixN<Currency, void*>::SliceLHS(const std::vector<hwSliceArg>& sliceArg,
	const hwTMatrixN<Currency, void*>& rhsMatrix)
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

        for (int i = 0; i < numSlices - 1; ++i)
            newDim[i] = m_dim[i];

        newDim[numSlices - 1] = -1;
        void* dataPtr;

        if (IsReal())
            dataPtr = (void*)m_real;
        else
            dataPtr = (void*)m_complex;

        hwTMatrixN<Currency, void*> reshaped(m_dim, dataPtr, Type());

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

    // map RHS dimensions to LHS dimensions
    bool newMatrix = false;

    if (m_dim.size() == 2 && m_dim[0] == 0 && m_dim[1] == 0)
        newMatrix = true;

    std::vector<int> newDim(numSlices);
    std::vector<int> dimMap(numSlices);     // rhsDim = dimMap[lhsDim]
    int lhsDim = 0;
    int numColons = 0;

    for (int i = 0; i < numSlices; ++i)
    {
        if (sliceArg[i].IsColon())
            ++numColons;
    }

    for (int rhsDim = 0; rhsDim < rhsMatrix.m_dim.size(); ++rhsDim)
    {
        if (lhsDim >= sliceArg.size())
            throw hwMathException(HW_MATH_ERR_SLICE_INDEX);

        if (sliceArg[lhsDim].IsScalar())
        {
            if (sliceArg[lhsDim].Scalar() < 0)
                throw hwMathException(HW_MATH_ERR_SLICE_INDEX);

            if (newMatrix)
                newDim[lhsDim] = sliceArg[lhsDim].Scalar() + 1;
            else if (lhsDim < m_dim.size())
                newDim[lhsDim] = _max(m_dim[lhsDim], static_cast<int> (sliceArg[lhsDim].Scalar()) + 1);
            else
                newDim[lhsDim] = static_cast<int> (sliceArg[lhsDim].Scalar()) + 1;

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

    for (; lhsDim < numSlices; ++lhsDim)
    {
        if (sliceArg[lhsDim].IsScalar())
        {
            if (lhsDim < m_dim.size())
                newDim[lhsDim] = _max(m_dim[lhsDim], static_cast<int> (sliceArg[lhsDim].Scalar() + 1));
            else
                newDim[lhsDim] = static_cast<int> (sliceArg[lhsDim].Scalar() + 1);

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
        Allocate(REAL);
        // SetElements(0.0);
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

            CopyData(m_real + startLHS, strideLHS, rhsMatrix.m_real + startRHS, strideRHS, keySize[S_case]);

            i += keySize[S_case] - 1;
            break;
        }
        case 3:
        {
            if (sliceArg[0].IsScalar())
            {
                (*this)(m_lhsMatrixIndex) = rhsMatrix(i);
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
                    (*this)(pos) = rhsMatrix(k++);
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
                if (m_lhsMatrixIndex[j] < (int)m_dim[j] - 1)
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
                if (m_rhsMatrixIndex[dimMap[j]] < (int)sliceArg[j].Vector().size() - 1)
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
