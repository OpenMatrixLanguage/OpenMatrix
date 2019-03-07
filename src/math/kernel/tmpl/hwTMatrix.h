/**
* @file hwTMatrix.h
* @date April 2012
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

#ifndef _hwTMatrix_h
#define _hwTMatrix_h

//! Forward declarations
class hwMathStatus;
template < typename T > class hwTComplex;

// -------------------------------------------------------------------
//! \class hwTMatrix hwTMatrix.h math/core/hwTMatrix.h
//! \brief A class to represent an MxN matrix of types T1 and T2
//!
//! Features of this matrix class are :
//! 1. Data is stored in contiguous memory to allow easy interface with third party libraries.
//! 2. Storage can be either by rows or by columns and can switch between the two at a cost.
//! 3. Complex numbers are stored as consecutive real and imaginary pairs in the contiguous memory.
//! 4. A real matrix becomes complex when an operation requires complex data.
//!
//! NOTE : If you require a complex matrix with separate real and imaginary arrays, please
//!        contact HMath development team with the use case for a strategy as to how to handle it.
// -------------------------------------------------------------------

//! hwTMatrix class definition
template < typename T1, typename T2 = hwTComplex<T1> > class hwTMatrix
{
public:
    // ****************************************************
    //                       Enums
    // ****************************************************

    //! Define flags that specify whether a matrix has real or complex data
    enum DataType { REAL, COMPLEX };
    //! Define flags to be used in ::operator() to specify whether to step
    //! through the elements by row or by column when using a single index
    enum IndexDir { BY_ROW, BY_COL };

    // ****************************************************
    //       Construction / Destruction / Assignment
    // ****************************************************

    //! Construct empty matrix
    hwTMatrix();
    //! Construct vector
    hwTMatrix(int size, DataType dataType);
    //! Construct matrix
    hwTMatrix(int m, int n, DataType dataType);
    //! Construct vector with external data
    hwTMatrix(int size, void* data, DataType dataType);
    //! Construct matrix with external data
    hwTMatrix(int m, int n, void* data, DataType dataType);
    //! Copy constructor
    hwTMatrix(const hwTMatrix<T1, T2>& source);
    //! Destructor
    ~hwTMatrix();
    //! Implement the = operator
    hwTMatrix<T1, T2>& operator=(const hwTMatrix<T1, T2>& rhs);

    // ****************************************************
    //             Data Type and Ownership
    // ****************************************************

    //! Determine if the matrix contains real or complex data memory
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

    //! Dimension matrix as a vector, specifying real or complex
    //! Applies to empty or populated matrix, does not preserve existing data
    hwMathStatus Dimension(int size, DataType dataType);
    //! Dimension matrix for double indexing, specifying real or complex
    //! Applies to empty or populated matrix, does not preserve existing data
    hwMathStatus Dimension(int m, int n, DataType dataType);
    //! Resize the matrix as a vector, leaving the data type unchanged
    //! Applies to populated matrix, can only resize an empty matrix if it remains empty
    hwMathStatus Resize(int size, bool initZero = false);
    //! Resize the matrix for double indexing, leaving the data type unchanged
    //! Applies to populated matrix, can only resize an empty matrix if it remains empty
    hwMathStatus Resize(int m, int n, bool initZero = false);
    //! Change the dimensions of a matrix while maintaining the same number of elements
    hwMathStatus Reshape(int m, int n);

    // ****************************************************
    //                  Matrix Properties
    // ****************************************************

    //! Return the number of rows
    int M() const { return m_nRows; }
    //! Return the number of columns
    int N() const { return m_nCols; }
    //! Return the number of elements
    int Size() const { return m_nCols * m_nRows; }
    //! Determine if the matrix is empty
    bool IsEmpty() const { return ((m_real || m_complex) ? false : true); }
	//! Determine if the matrix is 0x0
	bool Is0x0() const { return ((m_nCols == 0) && (m_nRows == 0)); }
    //! Determine if the matrix is a vector
    bool IsVector() const { return ((m_real||m_complex)&&(m_nCols==1||m_nRows==1) ? true : false); }
    //! Determine if the matrix is empty or a vector
    bool IsEmptyOrVector() const { return ((m_nCols<2 || m_nRows<2) ? true : false); }
    //! Determine if the matrix is square
    bool IsSquare() const { return (m_nCols == m_nRows ? true : false); }
    //! Determine if the matrix is symmetric
    bool IsSymmetric(T1 tol = (T1) 0) const;
    //! Determine if the matrix is Hermitian
    bool IsHermitian(T1 tol = (T1) 0) const;
    //! Determine if the matrix contains non-finite elements
    bool IsFinite() const;
    //! Determine if the matrix contains non-finite elements
    hwMathStatus IsFinite(hwTMatrix<bool>& R) const;
    //! Determine if the matrix contains non-finite elements
    hwMathStatus IsFinite(hwTMatrix<double>& R) const;

    // ****************************************************
    //         Access Functions for Real Elements
    // ****************************************************

    //! Return the T1* data vector pointer
    const T1* GetRealData() const { return m_real; }
    //! Return the T1* data vector pointer
    T1* GetRealData() { return m_real; }
    //! Return a reference to the real vector element at the specified single index
    T1& operator()(int index) { return m_real[index]; }
    //! Return a const reference to the real data element at the specified single index
    const T1& operator()(int index) const { return m_real[index]; }
    //! Return a reference to the real matrix element at the specified single index
    T1& operator()(int index, IndexDir dir);
    //! Return a const reference to the real data element at the specified single index
    const T1& operator()(int index, IndexDir dir) const;
    //! Return a reference to the real matrix element at the specified index pair
    T1& operator()(int i, int j);
    //! Return a const reference to the real data element at the specified index pair
    const T1& operator()(int i, int j) const;
    //! Set the element at the specified single index to a real value
    void SetElement(int index, T1 real, IndexDir dir = BY_COL);
    //! Set the element at the specified index pair to a real value
    void SetElement(int i, int j, T1 real);
    //! Set every matrix element to a specified real value
    void SetElements(T1 real);

    // ****************************************************
    //        Access Functions for Complex Elements
    // ****************************************************

    //! Return the T2* data vector pointer
    const T2* GetComplexData() const { return m_complex; };
    //! Return the T2* data vector pointer
    T2* GetComplexData() { return m_complex; };
    //! Return a reference to the complex vector element at the specified single index
    T2& z(int index) { return m_complex[index]; }
    //! Return a const reference to the complex data element at the specified single index
    const T2& z(int index) const { return m_complex[index]; }
    //! Return a reference to the complex matrix element at the specified single index
    T2& z(int index, IndexDir dir);
    //! Return a const reference to the complex data element at the specified single index
    const T2& z(int index, IndexDir dir) const;
    //! Return a reference to the complex matrix element at the specified index pair
    T2& z(int i, int j);
    //! Return a const reference to the complex data element at the specified index pair
    const T2& z(int i, int j) const;
    //! Set the element at the specified single index to a complex value
    void SetElement(int index, T2 cplx, IndexDir dir = BY_COL);
    //! Set the element at the specified index pair to a complex value
    void SetElement(int i, int j, T2 cplx);
    //! Set every matrix element to a specified complex value
    hwMathStatus SetElements(T2 cplx);

    // ****************************************************
    //             Real / Complex Conversions
    // ****************************************************

    //! Convert real matrix to complex
    hwMathStatus MakeComplex();
    //! Pack real and imaginary components into a complex matrix
    hwMathStatus PackComplex(const hwTMatrix<T1, T2>& real, const hwTMatrix<T1, T2>* imag = 0);
    //! Unpack real and imaginary components from a complex matrix
    hwMathStatus UnpackComplex(hwTMatrix<T1, T2>* real, hwTMatrix<T1, T2>* imag) const;

    // ****************************************************
    //                    Vector Operations
    // ****************************************************

    //! Insert a vector at the specified location of a source vector
    //! and store the result in the calling object
    hwMathStatus InsertElements(const hwTMatrix<T1, T2>& source, int startElem, const hwTMatrix<T1, T2>& elems);
    //! Insert a set zero elements at the specified location of a source vector
    //! and store the result in the calling object
    hwMathStatus InsertElements(const hwTMatrix<T1, T2>& source, int startElem, int numElems = 1);
    //! Delete a number of elements from a source vector at the specified location
    //! and store the result in the calling object
    hwMathStatus DeleteElements(const hwTMatrix<T1, T2>& source, int startElem, int numElems = 1);
    //! Delete a number of elements from the calling object at the specified location
    hwMathStatus DeleteElements(int startElem, int numElems = 1);

    // ****************************************************
    //                     Row Operations
    // ****************************************************

    //! Create a single vector from the rows of a matrix
    hwMathStatus ConcatRows(const hwTMatrix<T1, T2>& source);
    //! Insert a set of row vectors at the specified location of a source matrix
    //! and store the result in the calling object
    hwMathStatus InsertRows(const hwTMatrix<T1, T2>& source, int startRow, const hwTMatrix<T1, T2>& rows);
    //! Insert a set of zeroed rows at the specified location of a source matrix
    //! and store the result in the calling object
    hwMathStatus InsertRows(const hwTMatrix<T1, T2>& source, int startRow, int numRows = 1);
    //! Delete a number of rows from a source matrix at the specified location
    //! and store the result in the calling object
    hwMathStatus DeleteRows(const hwTMatrix<T1, T2>& source, int startRow, int numRows = 1);
    //! Delete a number of rows from the calling object at the specified location
    hwMathStatus DeleteRows(int startRow, int numRows = 1);
    //! Read the row at the specified location
    hwMathStatus ReadRow(int rowNum, hwTMatrix<T1, T2>& row) const;
    //! Write the row at the specified location
    hwMathStatus WriteRow(int rowNum, const hwTMatrix<T1, T2>& row);

    // ****************************************************
    //                    Column Operations
    // ****************************************************

    //! Create a single vector from the columns of a matrix
    hwMathStatus ConcatColumns(const hwTMatrix<T1, T2>& source);
    //! Insert a set of column vectors at the specified location of a source matrix
    //! and store the result in the calling object
    hwMathStatus InsertColumns(const hwTMatrix<T1, T2>& source, int startCol, const hwTMatrix<T1, T2>& cols);
    //! Insert a set of zeroed columns at the specified location of a source matrix
    //! and store the result in the calling object
    hwMathStatus InsertColumns(const hwTMatrix<T1, T2>& source, int startCol, int numCols = 1);
    //! Delete a number of columns from a source matrix at the specified location
    //! and store the result in the calling object
    hwMathStatus DeleteColumns(const hwTMatrix<T1, T2>& source, int startCol, int numCols = 1);
    //! Delete a number of columns from the calling object at the specified location
    hwMathStatus DeleteColumns(int startCol, int numCols = 1);
    //! Read the column at the specified location
    hwMathStatus ReadColumn(int colNum, hwTMatrix<T1, T2>& col) const;
    //! Write the column at the specified location
    hwMathStatus WriteColumn(int colNum, const hwTMatrix<T1, T2>& col);

    // ****************************************************
    //                   Submatrix Operations
    // ****************************************************

    //! Write a matrix to the calling object, starting at the specified location
    hwMathStatus WriteSubmatrix(int startRow, int startCol, const hwTMatrix<T1, T2>& source);
    //! Write a vector to the calling object, starting at the specified location
    hwMathStatus WriteSubmatrix(int startElem, const hwTMatrix<T1, T2>& source);

    // ****************************************************
    //           Misc Linear Algebra Operations
    // ****************************************************

    //! Set the diagonal to ones, with zeros elsewhere
    void Identity();
    //! Generate a diagonal matrix or extract a diagonal from a matrix
    hwMathStatus Diag(const hwTMatrix<T1, T2>& source, int k);
    //! Transpose the matrix in place
    hwMathStatus Transpose();
    //! Transpose the matrix argument
    hwMathStatus Transpose(const hwTMatrix<T1, T2>& source);
    //! Conjugate the matrix in place
    void Conjugate();
    //! Conjugate the matrix argument
    void Conjugate(const hwTMatrix<T1, T2>& source);
    //! Transpose and conjugate the matrix in place
    void Hermitian();
    //! Transpose and conjugate the matrix argument
    void Hermitian(const hwTMatrix<T1, T2>& source);
    //! Invert a matrix
    hwMathStatus Inverse(const hwTMatrix<double>& source);
    //! Matrix pseudo-inversion (Moore-Penrose)
    hwMathStatus Pinv(const hwTMatrix<double>& source);
    //! Real determinant
    hwMathStatus Determinant(double& det) const;
    //! Complex determinant
    hwMathStatus Determinant(hwTComplex<double>& det) const;
    //! Condition number
    hwMathStatus Cond(double& condNum) const;
    //! Reciprocal condition number estimate of a matrix
    hwMathStatus RCond(double& rCondNum) const;
    //! Rank
    hwMathStatus Rank(int& rank) const;
    //! Rank with a tolerance interval
    hwMathStatus Rank(double tol, int& rank) const;
    //! Matrix norm
    hwMathStatus Norm(double& norm, int p = 2) const;
    //! Matrix norm
    hwMathStatus Norm(double& norm, const char* type = "fro") const;
    //! Normalize a vector
    hwMathStatus Normalize();
    //! Matrix exponential (not a matrix of exponentials)
    hwMathStatus MatExp(const hwTMatrix<double>& power);
    //! Raise the object matrix to an integer power
    hwMathStatus Power(const hwTMatrix<double>& base, int power);
    //! Raise the object matrix (symmetric) to a possibly non integer power
    hwMathStatus Power(const hwTMatrix<double>& base, double power);
    //! Raise the object matrix (symmetric) to a complex power
    hwMathStatus Power(const hwTMatrix<double>& base, const hwTComplex<double>& power);

    // ****************************************************
    //               Arithmetic Operations
    // ****************************************************

    //! Add two matrices, sum.Add(A,B)
    hwMathStatus Add(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B);
    //! Add a matrix and a real number, sum.Add(A,real)
    hwMathStatus Add(const hwTMatrix<T1, T2>& A, T1 real);
    //! Add a matrix and a complex number, sum.Add(A,cmplx)
    hwMathStatus Add(const hwTMatrix<T1, T2>& A, const T2& cmplx);
    //! Add a matrix to the calling object
    hwMathStatus AddEquals(const hwTMatrix<T1, T2>& A);
    //! Add a real number to the calling object
    void AddEquals(T1 real);
    //! Add a complex number to the calling object
    hwMathStatus AddEquals(const T2& cmplx);
    //! Subtract one matrix from another
    hwMathStatus Subtr(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B);
    //! Subtract a real number from a matrix
    hwMathStatus Subtr(const hwTMatrix<T1, T2>& A, T1 real);
    //! Subtract a complex number from a matrix
    hwMathStatus Subtr(const hwTMatrix<T1, T2>& A, const T2& cmplx);
    //! Subtract a matrix from a real number
    hwMathStatus Subtr(T1 real, const hwTMatrix<T1, T2>& A);
    //! Subtract a matrix from a complex number
    hwMathStatus Subtr(const T2& cmplx, const hwTMatrix<T1, T2>& A);
    //! Negate a matrix
    hwMathStatus Negate(const hwTMatrix<T1, T2>& A);
    //! Subtract a matrix from the calling object
    hwMathStatus SubtrEquals(const hwTMatrix<T1, T2>& A);
    //! Subtract a real number from the calling object
    void SubtrEquals(T1 real);
    //! Subtract a complex number from the calling object
    hwMathStatus SubtrEquals(const T2& cmplx);
    //! Multiply two matrices
    hwMathStatus Mult(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B);
    //! Multiply a matrix and a real number
    hwMathStatus Mult(const hwTMatrix<T1, T2>& A, T1 real);
    //! Multiply a matrix and a complex number
    hwMathStatus Mult(const hwTMatrix<T1, T2>& A, const T2& cmplx);
    //! Multiply the calling object by a matrix
    hwMathStatus MultEquals(const hwTMatrix<T1, T2>& A);
    //! Multiply the calling object by a real number
    void MultEquals(T1 real);
    //! Multiply the calling object by a complex number
    hwMathStatus MultEquals(const T2& cmplx);
    //! Divide one matrix by another on the left side
    hwMathStatus DivideLeft(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B);
    //! Divide one matrix by another on the right side
    hwMathStatus DivideRight(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B);
    //! Divide a matrix by a real number
    hwMathStatus Divide(const hwTMatrix<T1, T2>& A, T1 real);
    //! Divide a matrix by a complex number
    hwMathStatus Divide(const hwTMatrix<T1, T2>& A, const T2& cmplx);
    //! Divide a real number by each element of a matrix
    hwMathStatus Divide(T1 real, const hwTMatrix<T1, T2>& A);
    //! Divide a complex number by each element of a matrix
    hwMathStatus Divide(const T2& cmplx, const hwTMatrix<T1, T2>& A);
    //! Divide the calling object by a matrix
    hwMathStatus DivideEquals(const hwTMatrix<T1, T2>& A);
    //! Divide the calling object by a real number
    void DivideEquals(T1 real);
    //! Divide the calling object by a complex number
    hwMathStatus DivideEquals(const T2& cmplx);
    //! Multiply two matrices element by element
    hwMathStatus MultByElems(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B);
    //! Divide two matrices element by element
    hwMathStatus DivideByElems(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B);
    //! Raise each matrix element to the specified real power
    hwMathStatus PowerByElems(const hwTMatrix<T1, T2>& Base, T1 power);
    //! Raise each matrix element to the specified complex power
    hwMathStatus PowerByElems(const hwTMatrix<T1, T2>& Base, T2 power);
    //! Raise each matrix element to the power in each element of another matrix
    hwMathStatus PowerByElems(const hwTMatrix<T1, T2>& Base, const hwTMatrix<T1, T2>& Power);
    //! Raise a real number to the power in each element of a matrix
    hwMathStatus PowerByElems(T1 base, const hwTMatrix<T1, T2>& Power);
    //! Raise a complex number to the power in each element of a matrix
    hwMathStatus PowerByElems(T2 base, const hwTMatrix<T1, T2>& Pow);

    // ****************************************************
    //         Decomposition and Solver Functions
    // ****************************************************

    //! Solve linear system AX=B, where X = *this
    //! uses LU decomposition when square, and QR otherwise
    hwMathStatus LSolve(const hwTMatrix<double>& A, const hwTMatrix<double>& B);
    //! Solve symmetric positive definite linear system AX=B with Cholesky decomposition,
    //! where X = *this
    hwMathStatus LSolveSPD(const hwTMatrix<double>& A, const hwTMatrix<double>& B);
    //! Solve tridigonal linear system AX=B with LU decomposition, where X = *this
    hwMathStatus LSolveT(const hwTMatrix<double>& D, const hwTMatrix<double>& DL,
                         const hwTMatrix<double>& DU, const hwTMatrix<double>& B);
    //! Solve band linear system AX=B with LU decomposition, where X = *this
    hwMathStatus LSolveB(const hwTMatrix<double>& A, int kl, int ku,
                         const hwTMatrix<double>& B);
    //! Solve SPD tridigonal linear system AX=B with Cholesky decomposition, where X = *this
    hwMathStatus LSolveSPDT(const hwTMatrix<double>& D, const hwTMatrix<double>& E,
                            const hwTMatrix<double>& B);
    //! Solve SPD band linear system AX=B with Cholesky decomposition, where X = *this
    hwMathStatus LSolveSPDB(const hwTMatrix<double>& A, int kd, const hwTMatrix<double>& B);
    //! Solve symmetric indefinite linear system AX=B, where X = *this
    hwMathStatus LSolveSI(const hwTMatrix<double>& A, const hwTMatrix<double>& B);
    //! Solve linear system AX=B with QR decomposition, where X = *this
    hwMathStatus QRSolve(const hwTMatrix<double>& A, const hwTMatrix<double>& B);
    //! Solve real linear system AX=B with SVD, where X = *this
    hwMathStatus SVDSolve(const hwTMatrix<double>& A, const hwTMatrix<double>& B);
    //! LU decomposition (PA=LU, where A = *this)
    hwMathStatus LU(hwTMatrix<int>& P, hwTMatrix<double>& L, hwTMatrix<double>& U) const;
    //! Cholesky decomposition (A=TT', where A = *this, and T is triangular)
    hwMathStatus Csky(hwTMatrix<double>& T, bool upper = false) const;
    //! Eigen decomposition (AV=DV, where A = *this)
    hwMathStatus Eigen(bool balance, hwTMatrix<double>* V, hwTMatrix<double>& D) const;
    //! Eigen decomposition - Symmetric/Hermitian  (AV=DV, where A = *this)
    hwMathStatus EigenSH(hwTMatrix<double>* V, hwTMatrix<double>& D) const;
    //! Eigen decomposition of a symmetric tridigonal matrix
    static hwMathStatus EigenST(const hwTMatrix<double>& D, const hwTMatrix<double>& E,
                                hwTMatrix<double>& W, hwTMatrix<double>& Z);
    //! Eigen decomposition of a symmetric band matrix
    hwMathStatus EigenSB(int kd, hwTMatrix<double>& W, hwTMatrix<double>& Z) const;
    //! Generalized Eigen decomposition
    hwMathStatus Eigen(const hwTMatrix<double>& A, const hwTMatrix<double>& B, hwTMatrix<double>& V, hwTMatrix<double>& D) const;
    //! Singular Value Decomposition (A=U*S*VT, where A = *this)
    hwMathStatus SVD(int flagSvd, hwTMatrix<double>* U, hwTMatrix<double>& S, hwTMatrix<double>* V) const;
    //! QR decomposition (A=QR, where A = *this)
    hwMathStatus QR(hwTMatrix<double>& Q, hwTMatrix<double>& R) const;
    // Schur decomposition
    hwMathStatus Schur(bool real, hwTMatrix<double>& U, hwTMatrix<double>& T) const;
    //! Balance matrix
    hwMathStatus Balance(bool noperm, hwTMatrix<double>& S, hwTMatrix<double>& P, hwTMatrix<double>& B) const;

    // ****************************************************
    //               Arithmetic Operators
    // ****************************************************

    //! Implement the == operator
    bool operator==(const hwTMatrix<T1, T2>& A) const;
    //! Determine if the matrix is equal to another
    bool IsEqual(const hwTMatrix<T1, T2>& A, T1 tol = (T1) 0) const;
    //! Implement the != operator
    bool operator!=(const hwTMatrix<T1, T2>& A) const;
    //! Implement the + operator with a matrix argument
    hwTMatrix<T1, T2> operator+(const hwTMatrix<T1, T2>& A) const;
    //! Implement the + operator with a real argument
    hwTMatrix<T1, T2> operator+(T1 real) const;
    //! Implement the + operator with a complex argument
    hwTMatrix<T1, T2> operator+(const T2& cmplx) const;
    //! Implement the += operator with a matrix argument
    hwTMatrix<T1, T2>& operator+=(const hwTMatrix<T1, T2>& A);
    //! Implement the += operator with a real argument
    hwTMatrix<T1, T2>& operator+=(T1 real);
    //! Implement the += operator with a complex argument
    hwTMatrix<T1, T2>& operator+=(const T2& cmplx);
    //! Implement the - operator with a matrix argument
    hwTMatrix<T1, T2> operator-(const hwTMatrix<T1, T2>& A) const;
    //! Implement the - operator with a real argument
    hwTMatrix<T1, T2> operator-(T1 real) const;
    //! Implement the - operator with a complex argument
    hwTMatrix<T1, T2> operator-(const T2& cmplx) const;
    //! Implement the - operator as urnary prefix
    hwTMatrix<T1, T2> operator-() const;
    //! Implement the -= operator with a matrix argument
    hwTMatrix<T1, T2>& operator-=(const hwTMatrix<T1, T2>& A);
    //! Implement the -= operator with a real argument
    hwTMatrix<T1, T2>& operator-=(T1 real);
    //! Implement the -= operator with a complex argument
    hwTMatrix<T1, T2>& operator-=(const T2& cmplx);
    //! Implement the * operator with a matrix argument
    hwTMatrix<T1, T2> operator*(const hwTMatrix<T1, T2>& A) const;
    //! Implement the * operator with a real argument
    hwTMatrix<T1, T2> operator*(T1 real) const;
    //! Implement the * operator with a complex argument
    hwTMatrix<T1, T2> operator*(const T2& cmplx) const;
    //! Implement the *= operator with a matrix argument
    hwTMatrix<T1, T2>& operator*=(const hwTMatrix<T1, T2>& A);
    //! Implement the *= operator with a real argument
    hwTMatrix<T1, T2>& operator*=(T1 real);
    //! Implement the *= operator with a complex argument
    hwTMatrix<T1, T2>& operator*=(const T2& cmplx);
    //! Implement the / operator with a matrix argument
    hwTMatrix<T1, T2> operator/(const hwTMatrix<T1, T2>& A) const;
    //! Implement the / operator with a real argument
    hwTMatrix<T1, T2> operator/(T1 real) const;
    //! Implement the / operator with a complex argument
    hwTMatrix<T1, T2> operator/(const T2& cmplx) const;
    //! Implement the /= operator with a matrix argument
    hwTMatrix<T1, T2>& operator/=(const hwTMatrix<T1, T2>& A);
    //! Implement the /= operator with a real argument
    hwTMatrix<T1, T2>& operator/=(T1 real);
    //! Implement the /= operator with a complex argument
    hwTMatrix<T1, T2>& operator/=(const T2& cmplx);

    // ****************************************************
    //                 Vector Operations
    // ****************************************************

    //! L2 vector norm squared
    hwMathStatus L2NormSq(T1& normSq) const;
    //! L2 vector norm
    hwMathStatus L2Norm(double& norm) const;
    //! Dot product of two real vectors
    static hwMathStatus Dot(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B, T1& dot);
    //! Dot product of two vectors, either real or complex
    static hwMathStatus Dot(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B, T2& dot);

    //! Kronecker product of two matrices
    hwMathStatus Kronecker(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B);

    //! Cross product of two vectors
    hwMathStatus Cross(const hwTMatrix<T1, T2>& A, const hwTMatrix<T1, T2>& B);
    //! Linear convolution of two vectors in the time domain
    hwMathStatus ConvLin(const hwTMatrix<T1, T2>& X, const hwTMatrix<T1, T2>& Y);
    //! Linear correlation of two vectors in the time domain
    hwMathStatus CorrLin(const hwTMatrix<T1, T2>& X, const hwTMatrix<T1, T2>& Y);
    //! Linear correlation of matrix columns in the time domain
    hwMathStatus CorrLin(const hwTMatrix<T1, T2>& X);
    //! 2D convolution of two matrices
    hwMathStatus Conv2D(const hwTMatrix<T1, T2>& X, const hwTMatrix<T1, T2>& Y);
    //! 2D convolution of a matrix with a column vector and a row vector
    hwMathStatus Conv2D(const hwTMatrix<T1, T2>& col, const hwTMatrix<T1, T2>& row, const hwTMatrix<T1, T2>& X);

    // ****************************************************
    //           Magnitude / Phase Operations
    // ****************************************************

    //! Compute the magnitudes of a matrix of values
    hwMathStatus Abs(const hwTMatrix<double>& A);
    //! Compute the squared magnitudes of a matrix of values
    hwMathStatus AbsSq(const hwTMatrix<double>& A);
    //! Compute the phases of a matrix of values
    hwMathStatus Phase(const hwTMatrix<double>& A);
    //! Compute the phases of a matrix of values
    hwMathStatus Hypot(const hwTMatrix<double>& A, const hwTMatrix<double>& B);
    //! Unwrap a vector of phases
    hwMathStatus UnwrapVec(const hwTMatrix<double>& phase, double tol);
    //! Unwrap a matrix of phases
    hwMathStatus UnwrapMat(const hwTMatrix<double>& phase, int dim, double tol);

    // ****************************************************
    //             Specialized Format Support
    // ****************************************************

    // The following functions support band and symmetric
    // band matrices.
    //! Dimension a band matrix
    hwMathStatus DimensionBandMatrix(int n, int kl, int ku, bool initZero = true);
    //! Dimension a symmetric band matrix
    hwMathStatus DimensionSymBandMatrix(int n, int kd, bool initZero = true);
    //! Set a band matrix element
    hwMathStatus SetBandMatrixElem(int i, int j, double value, int kl, int ku);
    //! Set a symmetric band matrix element
    hwMathStatus SetSymBandMatrixElem(int i, int j, double value, int kd);

    // ****************************************************
    //                   COW Support
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
        realData: 1,
        //! Bit indicating data ownership (0 = external data, 1 = object owns data)
        ownData: 1,
        //! Bits for future use (0 = TBD, 1 = TBD)
        unused: 14;
    } m_bits;

    //! The number of rows
    int m_nRows;
    //! The number of columns
    int m_nCols;
    //! Allocated memory capacity (utilized memory may be smaller)
    int m_capacity;
    //! Pointer to unaligned memory for new T1
    char* m_real_memory;
    //! Pointer to unaligned memory for new T2
    char* m_complex_memory;
    //! Pointer to aligned memory to store data of type T1
    T1* m_real;
    //! Pointer to aligned memory to store data of type T2,
    //! which is hwTComplex<T1> by default
    T2* m_complex;
    //! Copy On Write reference counter
    int m_refCount;

    // ****************************************************
    //                 Private Utilities
    // ****************************************************

    //! Compute memory capacity for the matrix data
    void SetCapacity(DataType dataType);
    //! Allocate memory for the matrix data
    void Allocate(DataType dataType);
    //! Delete the matrix data
    void Deallocate();
    //! Release real memory
    void FreeMemory(char*& pMemory, T1*& pReal);
    //! Release complex memory
    void FreeMemory(char*& pMemory, T2*& pComplex);
    //! Set matrix to empty condition
    void MakeEmpty();
    //! Copy matrix data from a source
    hwMathStatus Copy(const hwTMatrix<T1, T2>& source);
    //! Copy a real submatrix from another matrix to *this
    void CopyBlock(const T1* real, int m, int n, int row1, int row2,
                   int col1, int col2, int ii, int jj);
    //! Copy a complex submatrix from another matrix to *this
    void CopyBlock(const T2* cmplx, int m, int n, int row1, int row2,
                   int col1, int col2, int ii, int jj);
    //! Copy data
    void CopyData(void* dest, int arraySize, const void* src, int count);
    //! Set a submatrix of *this to zeros
    void ZeroBlock(int row1, int row2, int col1, int col2);

    // ****************************************************
    //         Decomposition and Solver Functions
    // ****************************************************

    //! Real LU decomposition (PA = LU)
    hwMathStatus RealLU(hwTMatrix<double>& L, hwTMatrix<double>& U, hwTMatrix<int>& P ) const;
    //! Complex LU decomposition (PA = LU)
    hwMathStatus ComplexLU(hwTMatrix<double>& L, hwTMatrix<double>& U, hwTMatrix<int>& P ) const;
    //! Real asymmetric Eigen decomposition with balance option
    hwMathStatus EigenDecompReal(bool balance, hwTMatrix<double>* V, hwTMatrix<double>& D) const;
    //! Real symmetric Eigen decomposition
    hwMathStatus EigenDecompRealSymmetric(hwTMatrix<double>& V, hwTMatrix<double>& D) const;
    //! Complex non-Hermitian Eigen decomposition with balance option
    hwMathStatus EigenDecompComplex(bool balance, hwTMatrix<double>* V, hwTMatrix<double>& D) const;
    //! Complex Hermitian Eigen decomposition
    hwMathStatus EigenDecompComplexHermitian(hwTMatrix<double>& V, hwTMatrix<double>& D) const;
    //! Balance real matrix
    hwMathStatus BalanceReal(bool noperm, hwTMatrix<double>& S, hwTMatrix<double>& P, hwTMatrix<double>& B) const;
    //! Balance complex matrix
    hwMathStatus BalanceComplex(bool noperm, hwTMatrix<double>& S, hwTMatrix<double>& P, hwTMatrix<double>& B) const;
    //! Generalized real asymmetric Eigen decomposition
    static hwMathStatus GeneralizedEigenDecompReal(const hwTMatrix<double>& A, const hwTMatrix<double>& B,
                        hwTMatrix<double>& V, hwTMatrix<double>& D);
    //! Generalized complex non-Hermitian Eigen decomposition
    static hwMathStatus GeneralizedEigenDecompComplex(const hwTMatrix<double>& A,
                        const hwTMatrix<double>& B, hwTMatrix<double>& V, hwTMatrix<double>& D);
    //! Generalized real symmetric Eigen decomposition
    static hwMathStatus GeneralizedEigenDecompRealSymmetric(const hwTMatrix<double>& A,
                        const hwTMatrix<double>& B, hwTMatrix<double>& V, hwTMatrix<double>& D);
    //! Generalized complex Hermitian Eigen decomposition
    static hwMathStatus GeneralizedEigenDecompComplexHermitian(const hwTMatrix<double>& A,
                        const hwTMatrix<double>& B, hwTMatrix<double>& V, hwTMatrix<double>& D);
    //! enum to define Singular Value decomposition type
    enum SVDtype
    {
        SVD_STD_FLAG,   // full SVD
        SVD_ECON_FLAG,  // economy size with S with reduced VT
        SVD_ECON_0_FLAG // economy size with S with reduced U
    };
    //! Real Singular Value Decomposition
    hwMathStatus RealSVD(const SVDtype& type, hwTMatrix<double>* U, hwTMatrix<double>& S, hwTMatrix<double>* VT) const;
    //! Complex Singular Value Decomposition
    hwMathStatus ComplexSVD(const SVDtype& type, hwTMatrix<double>* U, hwTMatrix<double>& S, hwTMatrix<double>* VT) const;
    //! Real matrix pseudo-inversion (Moore-Penrose)
    hwMathStatus RealPinv(const hwTMatrix<double>& source);
    //! Complex matrix pseudo-inversion (Moore-Penrose)
    hwMathStatus ComplexPinv(const hwTMatrix<double>& source);
    //! Real QR decomposition
    hwMathStatus RealQR(hwTMatrix<double>& Q, hwTMatrix<double>& R) const;
    //! Complex QR decomposition
    hwMathStatus ComplexQR(hwTMatrix<double>& Q, hwTMatrix<double>& R) const;
    // Real Schur decomposition
    hwMathStatus RealSchur(hwTMatrix<double>& U, hwTMatrix<double>& T) const;
    // Complex Schur decomposition
    hwMathStatus ComplexSchur(hwTMatrix<double>& U, hwTMatrix<double>& T) const;
    //! Reciprocal condition number estimate of a real matrix
    hwMathStatus RealRCond(double& rCondNum) const;
    //! Reciprocal condition number estimate of a complex matrix
    hwMathStatus ComplexRCond(double& rCondNum) const;
};

//! template implementation file
#include <tmpl/hwTMatrix.cc>

//! template utility function file
#include <utl/hwTMatrixUtil.cc>

#endif // _hwTMatrix_h
