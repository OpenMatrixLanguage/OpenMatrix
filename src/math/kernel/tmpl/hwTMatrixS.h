/**
* @file hwTMatrixS.h
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

#ifndef _hwTMatrixS_h
#define _hwTMatrixS_h

#include <vector>
#include "tmpl/hwTMatrix.h"

//*******************************************************************
//                    MKL sparse prototypes
//*******************************************************************
#define MKL_INT int

//! Forward declarations
class hwSliceArg;

//! hwTMatrixS class definition
template < typename T1, typename T2 = hwTComplex<T1> > class hwTMatrixS
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
    hwTMatrixS();
    //! Construct sparse from dense
    explicit hwTMatrixS(const hwTMatrix<T1, T2>& A);
    //! Construct sparse from vector representation
    hwTMatrixS(const std::vector<int>&  ivec,
               const std::vector<int>&  jvec,
               const hwTMatrix<T1, T2>& V,
               int                      m = -1,
               int                      n = -1,
               const char*              option = nullptr);
    //! Construct real sparse from MKL representation
    hwTMatrixS(int            nRows,
               int            nCols,
               const MKL_INT* pBeginCount,
               const MKL_INT* pEndCount,
               const MKL_INT* pRowNum,
               const T1*      pValues);
    //! Construct complex sparse from MKL representation
    hwTMatrixS(int            nRows,
               int            nCols,
               const MKL_INT* pBeginCount,
               const MKL_INT* pEndCount,
               const MKL_INT* pRowNum,
               const T2*      pValues);
    //! Copy constructor
    hwTMatrixS(const hwTMatrixS<T1, T2>& source);
    //! Destructor
    ~hwTMatrixS();
    //! Implement the = operator
    hwTMatrixS<T1, T2>& operator=(const hwTMatrixS<T1, T2>& rhs);

    // ****************************************************
    //               Data Type, Ownership
    // ****************************************************

    //! Determine if the matrix contains real or complex data
    bool IsReal() const { return m_values.IsReal(); };
    //! Determine if the matrix contains real or complex data
    //! (i.e. no non-zero imaginary components)
    bool IsRealData() const { return m_values.IsRealData(); };
    //! Return the data type
    DataType Type() const { return m_values.Type(); }

    // ****************************************************
    //                Set Matrix Dimensions
    // ****************************************************

    //! Change the dimensions of a matrix while maintaining the same number of elements
    void Reshape(int m, int n);

    // ****************************************************
    //                  Matrix Properties
    // ****************************************************

    //! Return the number of rows
    int M() const { return m_nRows; }
    //! Return the number of columns
    int N() const { return m_nCols; }
    //! Return the number of elements
    int Size() const { return m_nCols * m_nRows; }
    //! Return the number of non-zero elements
    int NNZ() const { return m_values.Size(); }
    //! Determine if the matrix is empty
    bool IsEmpty() const { return !(m_nCols && m_nRows); }
    //! Determine if the matrix is 0x0
    bool Is0x0() const { return (m_nCols == 0 && m_nRows == 0); }
    //! Determine if the matrix is a vector
    bool IsVector() const { return (m_values.Size() && (m_nCols == 1 || m_nRows == 1) ? true : false); }

    // ****************************************************
    //                 Dense Functions
    // ****************************************************

    //! Create dense from sparse
    void Full(hwTMatrix<T1, T2>& A) const;
    //! Find the real data for a storage index
    void NZinfo(int index, int& row, int& col, T1& value) const;
    //! Find the complex data for a storage index
    void NZinfo(int index, int& row, int& col, T2& value) const;
    //! Find the non-zero data
    void NZinfo(int first, int last, std::vector<int>& row,
                std::vector<int>& col, hwTMatrix<T1, T2>& value) const;
    //! Find the single indicies of the non-zero data
    void NZinfo(int first, int last, std::vector<int>& index) const;

    // ****************************************************
    //         Access Functions for Real Elements
    // ****************************************************

    //! Return the T1* data vector pointer
    const T1* GetRealData() const { return m_values.GetRealData(); }
    //! Return the T1* data vector pointer
    T1* GetRealData() { return m_values.GetRealData(); }
    //! Return a reference to the real data element at the specified single index
    T1& operator()(int index);
    //! Return the real data element at the specified single index
    T1 operator()(int index) const;
    //! Return a reference to the real data element at the specified indices
    T1& operator()(int i, int j);
    //! Return the real data element at the specified indices
    T1 operator()(int i, int j) const;
    //! Set every matrix element to a specified real value
    void SetElements(T1 real) { m_values.SetElements(real); }
    //! Zero the element at the specified single index
    void ZeroElement(int index);
    //! Zero the element at the specified specified indices
    void ZeroElement(int i, int j);

    // ****************************************************
    //        Access Functions for Complex Elements
    // ****************************************************

    //! Return the T2* data vector pointer
    const T2* GetComplexData() const { return m_values.GetComplexData(); }
    //! Return the T2* data vector pointer
    T2* GetComplexData() { return m_values.GetComplexData(); }
    //! Return a reference to the complex data element at the specified single index
    T2& z(int index);
    //! Return the complex data element at the specified single index
    T2 z(int index) const;
    //! Return a reference to the complex data element at the specified indices
    T2& z(int i, int j);
    //! Return the complex data element at the specified indices
    T2 z(int i, int j) const;
    //! Set every matrix element to a specified real value
    void SetElements(const T2& cmplx) { m_values.SetElements(cmplx); }

    // ****************************************************
    //             Real / Complex Conversions
    // ****************************************************

    //! Convert real matrix to complex
    void MakeComplex() { m_values.MakeComplex(); }
    //! Pack real and imaginary components into a complex matrix
    void PackComplex(const hwTMatrixS<T1, T2>& real, const hwTMatrixS<T1, T2>* imag = nullptr);
    //! Unpack real and imaginary components from a complex matrix
    void UnpackComplex(hwTMatrixS<T1, T2>* real, hwTMatrixS<T1, T2>* imag) const;

    // ****************************************************
    //                  Slice Operations
    // ****************************************************

    //! Read a matrix slice from the calling object, as if the calling object is being
    //! sliced on the the right hand side of an equals sign
    void SliceRHS(const std::vector<hwSliceArg>& sliceArg, hwTMatrixS<T1, T2>& lhsMatrix) const;
    //! Write to a matrix slice of the calling object, as if the calling object is being
    //! sliced on the the left hand side of an equals sign
    void SliceLHS(const std::vector<hwSliceArg>& sliceArg, const hwTMatrix<T1, T2>& rhsMatrix);
    //! Write to a matrix slice of the calling object, as if the calling object is being
    //! sliced on the the left hand side of an equals sign
    void SliceLHS(const std::vector<hwSliceArg>& sliceArg, const hwTMatrixS<T1, T2>& rhsMatrix);
    //! Write to a matrix slice of the calling object, as if the calling object is being
    //! sliced on the the left hand side of an equals sign
    void SliceLHS(const std::vector<hwSliceArg>& sliceArg, T1 real);
    //! Write to a matrix slice of the calling object, as if the calling object is being
    //! sliced on the the left hand side of an equals sign
    void SliceLHS(const std::vector<hwSliceArg>& sliceArg, const T2& cmplx);

    // ****************************************************
    //                   Submatrix Operations
    // ****************************************************

    //! Write a submatrix source to the calling object, starting at the specified location
    void WriteSubmatrix(int startRow, int startCol, const hwTMatrixS<T1, T2>& source);
    //! Write a submatrix source to the calling object, starting at the specified location
    void WriteSubmatrix(int startRow, int startCol, const hwTMatrix<T1, T2>& source);
    //! Concatenate two matrices vertically
    void ConcatVertical(const hwTMatrixS<T1, T2>& A, const hwTMatrixS<T1, T2>& B);
    //! Concatenate two matrices horizontally
    void ConcatHorizontal(const hwTMatrixS<T1, T2>& A, const hwTMatrixS<T1, T2>& B);
    //! Expand and assign memory of a dimensioned sparse matrix, row-wise
    void FillEmptyDimensionedRows(int startRow, const hwTMatrixS<T1, T2>& B);
    //! Expand and assign memory of a dimensioned sparse matrix, column-wise
    void FillEmptyDimensionedColumns(int startCol, const hwTMatrixS<T1, T2>& B);

    // ****************************************************
    //           Misc Linear Algebra Operations
    // ****************************************************

    //! Set the diagonal to ones, with zeros elsewhere
    // void Identity();
    //! Generate a diagonal matrix or extract a diagonal from a matrix
    void Diag(const hwTMatrixS<T1, T2>& source, int k);
    //! Generate a diagonal matrix from a vector
    void Diag(const hwTMatrixS<T1, T2>& source, int m, int n);
    //! Transpose the matrix in place
    void Transpose();
    //! Transpose the matrix argument
    void Transpose(const hwTMatrixS<T1, T2>& source);
    //! Conjugate the matrix in place
    void Conjugate();
    //! Conjugate the matrix argument
    void Conjugate(const hwTMatrixS<T1, T2>& source);
    //! Transpose and conjugate the matrix in place
    void Hermitian();
    //! Transpose and conjugate the matrix argument
    void Hermitian(const hwTMatrixS<T1, T2>& source);
    //! Sum along a direction
    void Sum(const hwTMatrixS<T1, T2>& source, bool cols);
    //! Find maximum along a direction
    void Max(hwTMatrixS<T1, T2>& val, hwTMatrix<int, hwTComplex<int>>* rc, bool cols) const;
    //! Find minimum along a direction
    void Min(hwTMatrixS<T1, T2>& val, hwTMatrix<int, hwTComplex<int>>* rc, bool cols) const;

    // ****************************************************
    //               Arithmetic Operations
    // ****************************************************

    //! Add two sparse matrices, sum.Add(A,B)
    void Add(const hwTMatrixS<T1, T2>& A, const hwTMatrixS<T1, T2>& B);
    //! Subtract two sparse matrices, diff.Subtr(A,B)
    void Subtr(const hwTMatrixS<T1, T2>& A, const hwTMatrixS<T1, T2>& B);
    //! Multiply a matrix and a real number
    void Mult(const hwTMatrixS<T1, T2>& A, T1 real);
    //! Multiply a matrix and a complex number
    void Mult(const hwTMatrixS<T1, T2>& A, const T2& cmplx);
    //! Divide a matrix by a real number
    void Divide(const hwTMatrixS<T1, T2>& A, T1 real);
    //! Divide a matrix by a complex number
    void Divide(const hwTMatrixS<T1, T2>& A, const T2& cmplx);
    //! Negate a sparse matrix
    void Negate(const hwTMatrixS<T1, T2>& A);

    // ****************************************************
    //               Arithmetic Operators
    // ****************************************************
    
    //! Implement the == operator
    bool operator==(const hwTMatrixS<T1, T2>& A) const;
    //! Implement the != operator
    bool operator!=(const hwTMatrixS<T1, T2>& A) const;

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

    //! Copy On Write reference counter
    int m_refCount;

    //! The number of rows
    int m_nRows;
    //! The number of columns
    int m_nCols;

    //! The non-zero values
    hwTMatrix<T1, T2> m_values;
    //! The row numbers of the non-zero values
    std::vector<MKL_INT> m_rows;
    //! The cummulative non-zero value count at the Beginning of the columns
    std::vector<MKL_INT> m_pointerB;
    //! The cummulative non-zero value count at the End of the columns
    std::vector<MKL_INT> m_pointerE;

    // TODO: combine m_pointerB and m_pointerE into a single vector
public:
    const MKL_INT* rows() const { return m_rows.data(); }
    const MKL_INT* pointerB() const { return m_pointerB.data(); }
    const MKL_INT* pointerE() const { return m_pointerE.data(); }

    // ****************************************************
    //                 Private Utilities
    // ****************************************************
private:
    //! Set matrix to empty condition
    void MakeEmpty();
    //! Remove stored zeros
    void RemoveStoredZeros();
};

//! template implementation file
#include <tmpl/hwTMatrixS.cc>

#endif // _hwTMatrixS_h
