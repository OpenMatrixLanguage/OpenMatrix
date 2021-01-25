/**
* @file BuiltInFuncsMKL.h
* @date September 2021
* Copyright (C) 2016-2021 Altair Engineering, Inc.
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

#ifndef __BUILTINFUNCSMKL__
#define __BUILTINFUNCSMKL__

// Begin defines/includes
#include "EvaluatorInt.h"

// End defines/includes

//------------------------------------------------------------------------------
//!
//! \brief Class for built-in functions implementing MKL based functions
//!
//------------------------------------------------------------------------------
class OMLDLL_DECLS BuiltInFuncsMKL
{
public:
    //!
    //! Destructor
    //!
    ~BuiltInFuncsMKL() {}
    //!
    //! Returns cosine of input in radians [cos command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Cos(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns sine of input in radians [sin command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Sin(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns tangent of input in radians [tan command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Tan(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns secant of input in radians [sec command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Sec(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns cosecant of input in radians [csc command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Csc(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns cotangent of input in radians [cot command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Cot(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns cosine of input in degrees [cosd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool CosD(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns sine of input in degrees [sind command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool SinD(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns tangent of input in degrees [tand command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool TanD(EvaluatorInterface          eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns secant of input in degrees [secd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool SecD(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns cosecant of input in degrees [cscd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool CscD(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns cotangent of input in degrees [cotd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool CotD(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cosine of input in radians [acos command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCos(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse sine of input in radians [asin command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aSin(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse tangent of input in radians [atan command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aTan(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse secant of input in radians [asec command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aSec(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cosecant of input in radians [acsc command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCsc(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cotangent of input in radians [acot command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCot(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cosine of input in degrees [acosd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCosD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse sine of input in degrees [asind command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aSinD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse tangent of input in degrees [atand command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aTanD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse secant of input in degrees [asecd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aSecD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cosecant of input in degrees [acscd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCscD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cotangent of input in degrees [acotd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCotD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse tangent of input in radians [atan2 command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aTan2(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse tangent of input in degrees [atan2d command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aTan2D(EvaluatorInterface           eval,
                       const std::vector<Currency>& inputs,
                       std::vector<Currency>&       outputs);
    //!
    //! Returns hyperbolic cosine of input [cosh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Cosh(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns hyperbolic sine of input [sinh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Sinh(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns hyperbolic tangent of input [tanh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Tanh(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse hyperbolic sine of input [asinh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aSinh(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns exponential of input [exp command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Exp(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns inverse hyperbolic cosine of input [acosh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCosh(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse hyperbolic tangent of input [atanh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aTanh(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);

    //!
    //! Returns natural logarithm of input [log command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Log(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Returns logarithm base 2 of input [log2 command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Log2(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);

    //!
    //! Returns logarithm base 10 of input [log10 command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Log10(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);

    //!
    //! Returns square root of input [sqrt command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Sqrt(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);

    //!
    //! Returns magnitude of input [abs command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Abs(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Returns angle of input [abs command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Arg(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Rounds input to the nearest integer [round command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Round(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);

    //!
    //! Rounds input toward inf [ceil command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Ceil(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);

    //!
    //! Rounds input input toward -inf [floor command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Floor(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);

    //!
    //! Rounds input toward zero [fix command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Fix(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Computes division remainder [rem command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Rem(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Computes modulo function [mod command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Mod(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Returns element-wise maximum of two inputs [max command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Max(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Returns element-wise minimum of two inputs [min command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Min(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Conjugate of a complex
    //! \param A   Full input
    //! \param hyp Full output
    //!
    static void Conj(const hwMatrix& A,
                     hwMatrix&       conj);

    //!
    //! Hypotenuse of paired matrix elements (hyp = sqrt(A^2 + B^2)
    //! \param A   Full input
    //! \param B   Full input
    //! \param hyp Full output
    //!
    static void Hypot(const hwMatrix& A,
                      const hwMatrix& B,
                      hwMatrix&       hyp);

    //!
    //! Add two full matrices (sum = A + B)
    //! \param A   Full input
    //! \param B   Full input
    //! \param sum Full output
    //!
    static void Add(const hwMatrix& A,
                    const hwMatrix& B,
                    hwMatrix&       sum);

    //!
    //! Subtract two full matrices (diff = A - B)
    //! \param A    Full input
    //! \param B    Full input
    //! \param diff Full output
    //!
    static void Subtr(const hwMatrix& A,
                      const hwMatrix& B,
                      hwMatrix& diff);

    //!
    //! Multiply two full matrices element-wise (prod = A .* B)
    //! \param A    Full input
    //! \param B    Full input
    //! \param prod Full output
    //!
    static void MultByElems(const hwMatrix& A,
                            const hwMatrix& B,
                            hwMatrix&       prod);

    //!
    //! Divide two full matrices element-wise (quot = A ./ B)
    //! \param A    Full input
    //! \param B    Full input
    //! \param quot Full output
    //!
    static void DivideByElems(const hwMatrix& A,
                              const hwMatrix& B,
                              hwMatrix&       quot);

    //!
    //! Power operation applied element-wise (B = A .^ P)
    //! \param A    base input
    //! \param P    power input
    //! \param B    output
    //!
    static void PowerByElems(const hwMatrix& A,
                             const hwMatrix& P,
                             hwMatrix&       B);

    //!
    //! Power operation applied element-wise (B = A .^ P)
    //! \param A    base input
    //! \param P    power input
    //! \param B    output
    //!
    static void PowerByElems(const hwMatrix& A,
                             double          P,
                             hwMatrix&       B);

    //!
    //! Power operation applied element-wise (B = A .^ P)
    //! \param A    base input
    //! \param P    power input
    //! \param B    output
    //!
    static void PowerByElems(const hwMatrix&  A,
                             const hwComplex& P,
                             hwMatrix&        B);

    //!
    //! Add two sparse matrices (sum = A + B)
    //! \param A   Sparse input
    //! \param B   Sparse input
    //! \param sum Sparse output
    //!
    static void SparseAdd(const hwMatrixS& A,
                          const hwMatrixS& B,
                          hwMatrixS&       sum);

    //!
    //! Multiply a sparse matrix by a full matrix (prod = A * B)
    //! \param A    Sparse input
    //! \param B    Full input
    //! \param prod Full output
    //!
    static void SparseMult(const hwMatrixS& A,
                           const hwMatrix&  B,
                           hwMatrix&        prod);

    //!
    //! Multiply a full matrix by a sparse matrix (prod = A * B)
    //! \param A    Full input
    //! \param B    Sparse input
    //! \param prod Full output
    //!
    static void SparseMult(const hwMatrix&  A,
                           const hwMatrixS& B,
                           hwMatrix&        prod);

    //!
    //! Multiply two sparse matrices (prod = A * B)
    //! \param A    Sparse input
    //! \param B    Sparse input
    //! \param prod Sparse output
    //!
    // Multiply a sparse matrix by a sparse matrix (prod = A * B)
    static void SparseMult(const hwMatrixS& A,
                           const hwMatrixS& B,
                           hwMatrixS&       prod);

    //!
    //! Divide a full matrix by a sparse matrix on the left side (Q = A \ B)
    //! \param A Sparse input
    //! \param B Full input
    //! \param Q Full output
    //!
    static void SparseDivideLeft(const hwMatrixS& A,
                                 const hwMatrix&  B,
                                 hwMatrix&        Q);

    //!
    //! Divide a full matrix by a sparse matrix on the right side (Q = A / B)
    //! \param A Full input
    //! \param B Sparse input
    //! \param Q Full output
    //!
    static void SparseDivideRight(const hwMatrix&  A,
                                  const hwMatrixS& B,
                                   hwMatrix&       Q);

    //!
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    // Set pivot threshold for MKL sparse matrix division
    static bool oml_mkl_sdpt(EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs);

private:
    //!
    //! Constructor
    //!
    BuiltInFuncsMKL() {}
};

#endif // __BUILTINFUNCSMKL__
