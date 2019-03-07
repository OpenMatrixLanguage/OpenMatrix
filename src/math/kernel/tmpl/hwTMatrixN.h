/**
* @file hwTMatrixN.h
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

#ifndef _hwTMatrixN_h
#define _hwTMatrixN_h

#include <vector>
#include <tmpl/hwTMatrix.h>

//! Forward declarations
class hwSliceArg;

//! hwTMatrixN class definition
template < typename T1, typename T2 = hwTComplex<T1> > class hwTMatrixN
{
public:
    // ****************************************************
    //                       Enums
    // ****************************************************

    //! Define flags that specify whether a matrix has real or complex data
    enum DataType { REAL, COMPLEX };

    // ****************************************************
    //       Construction / Destruction / Assignment
    // ****************************************************

    //! Construct empty matrix
    hwTMatrixN();
    //! Construct dimensioned matrix
    hwTMatrixN(const std::vector<int>& dim, DataType dataType);
    //! Construct dimensioned matrix with external data
    hwTMatrixN(const std::vector<int>& dim, void* data, DataType dataType);
    //! Copy constructor
    hwTMatrixN(const hwTMatrixN<T1, T2>& source);

    //! Destructor
    ~hwTMatrixN();
    //! Implement the = operator
    hwTMatrixN<T1, T2>& operator=(const hwTMatrixN<T1, T2>& rhs);

    // ****************************************************
    //           hwTMatrix conversions functions
    // ****************************************************

    //! Convert hwTMatrixN to hwTMatrix 
    void ConvertNDto2D(hwTMatrix<T1, T2>& target) const;
    //! Convert hwTMatrix to hwTMatrixN 
    void Convert2DtoND(const hwTMatrix<T1, T2>& source);

    // ****************************************************
    //               Data Type, Ownership
    // ****************************************************

    //! Determine if the matrix contains real or complex data
    bool IsReal() const { return (m_complex ? false : true); }
    //! Determine if the matrix contains real or complex data
    //! (i.e. no non-zero imaginary components)
    bool IsRealData() const;
    //! Return the data type
    DataType Type() const { return (m_bits.realData ? REAL : COMPLEX); }
    //! Set the data ownership
    void OwnData(bool status) { m_bits.ownData = (int) status; }
    //! Determine the data ownership
    bool OwnData() const { return (m_bits.ownData ? true : false); }

    // ****************************************************
    //                Set Matrix Dimensions
    // ****************************************************

    //! Dimension matrix, specifying real or complex
    //! Applies to empty or populated matrix, does not preserve existing data
    void Dimension(const std::vector<int>& dim, DataType dataType);
    //! Change the dimensions of a matrix while maintaining the same number of elements
    void Reshape(const std::vector<int>& dim);
    //! Reorder matrix dimensions, a generalized transpose
    void Permute(const hwTMatrixN<T1, T2>& source, const std::vector<int>& permuteVec);
    //! Verify permutation vector validity
    static void PermuteCheck(const std::vector<int>& permuteVec);

    // ****************************************************
    //                  Matrix Properties
    // ****************************************************

    //! Return a reference to the dimension vector
    const std::vector<int>& Dimensions() const { return m_dim; }
    //! Return the number of elements
    int Size() const { return m_size; }
    //! Determine if the matrix is empty
    bool IsEmpty() const { return ((m_real || m_complex) ? false : true); }
    //! Determine if the matrix is a vector
    bool IsVector() const;
    //! Determine if the matrix is empty or a vector
    bool IsEmptyOrVector() const;
    //! Bounds check - lower only, since LHS ops allow resizing
    void BoundsCheckLHS(const std::vector<int>& indexVec) const;
    //! Bounds check - upper and lower
    void BoundsCheckRHS(const std::vector<int>& indexVec) const;
    //! Bounds check - lower only, since LHS ops allow resizing
    void BoundsCheckLHS(const std::vector<hwSliceArg>& sliceArg, int numSlices = 0) const;
    //! Bounds check - upper and lower
    void BoundsCheckRHS(const std::vector<hwSliceArg>& sliceArg, int numSlices = 0) const;
    //! Need to grow a matrix, for use with LHS indexing
    bool NeedToGrowLHS(const std::vector<int>& indexVec);
    //! Need to grow a matrix, for use with LHS slicing
    bool NeedToGrowLHS(const std::vector<hwSliceArg>& sliceArg, int numSlices) const;
    //! Grow a matrix
    void GrowLHSMatrix(const std::vector<int>& indexVec);
    //! Grow a matrix
    void GrowLHSMatrix(const std::vector<hwSliceArg>& sliceArg, int& numSlices,
                       const hwTMatrixN<T1, T2>* rhsMatrix = NULL, bool rule2 = false);

    // ****************************************************
    //               Indexing Functions
    // ****************************************************

    //! Return the index vector corresponding to a single index
    std::vector<int> IndexVector(int index) const;
    //! Return the single index corresponding to an index vector
    int Index(const std::vector<int>& indexVec) const;

    // ****************************************************
    //         Access Functions for Real Elements
    // ****************************************************

    //! Return the T1* data vector pointer
    const T1* GetRealData() const { return m_real; }
    //! Return the T1* data vector pointer
    T1* GetRealData() { return m_real; }
    //! Return a reference to the real data element at the specified single index
    T1& operator()(int index) { return m_real[index]; }
    //! Return a const reference to the real data element at the specified single index
    const T1& operator()(int index) const { return m_real[index]; }
    //! Return a reference to the real data element at the specified indices
    T1& operator()(const std::vector<int>& indexVec);
    //! Return a const reference to the real data element at the specified indices
    const T1& operator()(const std::vector<int>& indexVec) const;
    //! Set every matrix element to a specified real value
    void SetElements(T1 real);

    // ****************************************************
    //        Access Functions for Complex Elements
    // ****************************************************

    //! Return the T2* data vector pointer
    const T2* GetComplexData() const { return m_complex; };
    //! Return the T2* data vector pointer
    T2* GetComplexData() { return m_complex; }
    //! Return a reference to the complex data element at the specified single index
    T2& z(int index) { return m_complex[index]; }
    //! Return a const reference to the complex data element at the specified single index
    const T2& z(int index) const { return m_complex[index]; }
    //! Return a reference to the complex data element at the specified indices
    T2& z(const std::vector<int>& indexVec);
    //! Return a const reference to the complex data element at the specified indices
    const T2& z(const std::vector<int>& indexVec) const;
    //! Set every matrix element to a specified real value
    void SetElements(const T2& cmplx);

    // ****************************************************
    //             Real / Complex Conversions
    // ****************************************************

    //! Convert real matrix to complex
    void MakeComplex();
    //! Pack real and imaginary components into a complex matrix
    void PackComplex(const hwTMatrixN<T1, T2>& real, const hwTMatrixN<T1, T2>* imag = nullptr);
    //! Unpack real and imaginary components from a complex matrix
    void UnpackComplex(hwTMatrixN<T1, T2>* real, hwTMatrixN<T1, T2>* imag) const;

    // ****************************************************
    //                  Slice Operations
    // ****************************************************

    //! Read a matrix slice from the calling object, as if the calling object is being
    //! sliced on the the right hand side of an equals sign
    void SliceRHS(const std::vector<hwSliceArg>& sliceArg, hwTMatrixN<T1, T2>& lhsMatrix) const;
    //! Write to a matrix slice of the calling object, as if the calling object is being
    //! sliced on the the left hand side of an equals sign
    void SliceLHS(const std::vector<hwSliceArg>& sliceArg, const hwTMatrixN<T1, T2>& rhsMatrix);
    //! Write to a matrix slice of the calling object, as if the calling object is being
    //! sliced on the the left hand side of an equals sign
    void SliceLHS(const std::vector<hwSliceArg>& sliceArg, T1 real);
    //! Write to a matrix slice of the calling object, as if the calling object is being
    //! sliced on the the left hand side of an equals sign
    void SliceLHS(const std::vector<hwSliceArg>& sliceArg, const T2& cmplx);

    // ****************************************************
    //               Arithmetic Operators
    // ****************************************************
    
    //! Implement the == operator
    bool operator==(const hwTMatrixN<T1, T2>& A) const;
    //! Implement the != operator
    bool operator!=(const hwTMatrixN<T1, T2>& A) const;

    // ****************************************************
    //              Copy on Write functions
    // ****************************************************

    void IncrRefCount() { m_refCount++; }
    void DecrRefCount() { m_refCount--; }
    void ResetRefCount() { m_refCount = 1; }
    unsigned int GetRefCount() const { return m_refCount; }
    bool IsMatrixShared() const { return m_refCount != 1; }

private:
    // ****************************************************
    //                   Data Members
    // ****************************************************

    struct MatrixInfo
    {
        mutable unsigned short
        //! Bit indicating data type (0 = complex, 1 = real)
        realData : 1,
        //! Bit indicating data ownership (0 = external data, 1 = object owns data)
        ownData: 1,
        //! Bits for future use (0 = TBD, 1 = TBD)
        forRent: 14;
    } m_bits;

    //! Vector of dimensions
    std::vector<int> m_dim;
    //! Utilized size
    int m_size;
    //! Allocated memory capacity (utilized memory may be smaller)
    int m_capacity;
    //! Contiguous block of memory to store data of type T1
    T1* m_real;
    //! Contiguous block of memory to store data of type T2, which is hwTComplex<T> by default
    T2* m_complex;
    //! The cached memory index of the last accessed data element
    mutable int m_pos;
    //! Copy On Write reference counter
    int m_refCount;
    //! Utility vector for left hand side indexing
    mutable std::vector<int> m_lhsMatrixIndex;
    //! Utility vector for right hand side indexing
    mutable std::vector<int> m_rhsMatrixIndex;

    // ****************************************************
    //                 Private Utilities
    // ****************************************************

    //! Set Dimensions
    void SetDimensions(const std::vector<int>& dim);
    //! Compute size
    void ComputeSize();
    //! Compute memory capacity for the matrix data
    void SetCapacity(DataType dataType);
    //! Allocate memory for the matrix data
    void Allocate(DataType dataType);
    //! Delete the matrix data
    void Deallocate();
    //! Set matrix to empty condition
    void MakeEmpty();
    //! Set memory position corresponding to an index vector
    int SetMemoryPosition(const std::vector<int>& indexVec) const;
    //! Copy matrix data from a source
    void Copy(const hwTMatrixN<T1, T2>& source);
    //! Copy data
    void CopyData(void* dest, int arraySize, const void* src, int count) const;
    //! Ignore high dimension singleton indices
    int RelevantNumberOfSlices(const std::vector<hwSliceArg>& sliceArg,
                               bool matrixAssignment) const;
    //! Delete a matrix slice
    void DeleteSlice(const std::vector<hwSliceArg>& sliceArg);
    //! Read a contiguous block from the calling object, as if the calling
    //! object is being sliced on the the right hand side of an equals sign
    void CopyBlockRHS(int& pos, int sliceArg, hwTMatrixN<T1, T2>& lhsMatrix) const;
    //! Write a contiguous block to the calling object, as if the calling
    //! object is being sliced on the the left hand side of an equals sign
    void CopyBlockLHS(int& pos, int sliceArg, const hwTMatrixN<T1, T2>& rhsMatrix);
    //! Write a contiguous block to the calling object, as if the calling
    //! object is being sliced on the the left hand side of an equals sign
    void CopyBlockLHS(int& pos, int sliceArg, T1 real);
    //! Write a contiguous block to the calling object, as if the calling
    //! object is being sliced on the the left hand side of an equals sign
    void CopyBlockLHS(int& pos, int sliceArg, const T2& cmplx);
    //! Write a contiguous block to the calling object, as if the calling
    //! object is being sliced on the the left hand side of an equals sign
    void CopyMatrixLHS(const hwTMatrixN<T1, T2>& rhsMatrix);
    //! Transfer contents from another object
    void Transfer(hwTMatrixN<T1, T2>& source);
};

//! template implementation file
#include <tmpl/hwTMatrixN.cc>

#endif // _hwTMatrixN_h
