#include "hwMatrixN.h"
#include "Currency.h"

template<> inline
void hwTMatrixN<Currency, void*>::SliceLHS(const std::vector<hwSliceArg>& sliceArg, Currency real);

// !Write a contiguous block to the calling object, as if the calling
//! object is being sliced on the left hand side of an equals sign
template<> inline
void hwTMatrixN<Currency, void*>::CopyBlockLHS(int& pos, int sliceArg, const hwTMatrixN<Currency, void*>& rhsMatrix)
{
	int numVals = 1;
	int numArgs = _min(sliceArg, (int)m_dim.size());

	for (int i = 0; i < numArgs; ++i)
		numVals *= m_dim[i];

	// the rhsMatrix block begins at location pos
	// the lhsMatrix block begins at location m_pos, which must be set prior to the function call
	if (numVals > 0)
	{
		if (m_real)
			CopyData(m_real + m_pos, numVals, rhsMatrix.m_real + pos, numVals);
	}

	pos += numVals;
}

//! Write a contiguous matrix to the calling object, as if the
//! calling object is on the left hand side of an equals sign
template<> inline
void hwTMatrixN<Currency, void*>::CopyMatrixLHS(const hwTMatrixN<Currency, void*>& rhsMatrix)
{
	// determine dimension at which discontiguity occurs
	int discontiguity = (int)rhsMatrix.m_dim.size();

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
			(*this)(m_rhsMatrixIndex) = rhsMatrix(i);
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
			if (m_rhsMatrixIndex[j] < (int)rhsMatrix.m_dim[j] - 1)
			{
				++m_rhsMatrixIndex[j];
				break;
			}

			// index j is maxed out, so reset and continue to j+1
			m_rhsMatrixIndex[j] = 0;
		}
	}
}

//! Grow a matrix
template<> inline
void hwTMatrixN<Currency, void*>::GrowLHSMatrix(const std::vector<hwSliceArg>& sliceArg, int& numSlices,
	const hwTMatrixN<Currency, void*>* rhsMatrix, bool rule2)
{
	// create matrix with expanded dimensions
	int i;
	int common = _min(numSlices, (int)m_dim.size());
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
	int numDims = (int)m_dim.size();
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
		for (int i = firstChangedDim + 1; i < numDims; ++i)
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
					std::vector<int> vec(m_dim[i] - oldDim[i]);

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
		Currency* m_real_old = m_real;
		hwTMatrixN<Currency, void*> tempMatrix(oldDim, m_real_old, REAL);
		m_real = nullptr;
		Allocate(REAL);

		// copy old data to the expanded matrix
		CopyMatrixLHS(tempMatrix);

		if (m_bits.ownData)
			delete[] m_real_old;
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

//! Write to a matrix slice of the calling object, as if the calling object is being
//! sliced on the left hand side of an equals sign
template<> inline
void hwTMatrixN<Currency, void*>::SliceLHS(const std::vector<hwSliceArg>& sliceArg, Currency real)
{
	if (sliceArg.size() == 0)
		throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

	int i, j;

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

		for (i = 0; i < numSlices - 1; ++i)
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

				for (i = 0; i < numSlices - 1; ++i)
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
		//SetElements(0.0);
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
			size *= (int)sliceArg[i].Vector().size();
	}

	if (m_size < size)
		throw hwMathException(HW_MATH_ERR_INVALIDINPUT);

	for (i = 0; i < size; ++i)
	{
		if (discontiguity == 0)
		{
			(*this)(m_lhsMatrixIndex) = real;
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

				if (m_lhsMatrixIndex[j] < (int)sliceArg[j].Vector()[vecSize - 1])
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

		Allocate(REAL);
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

	for (int i = 0; i < size; ++i)
	{
		// copy data up to discontguity
		if (discontiguity == 0)
		{
			(*this)(m_lhsMatrixIndex) = rhsMatrix(i);
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
				if (m_rhsMatrixIndex[j] < (int)sliceArg[j].Vector().size() - 1)
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




