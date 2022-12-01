/**
* @file MKL_Tutilities.cxx
* @date March 2022
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

// Begin defines/includes
#include "../hwComplex.h"
#include "../hwMatrix.h"
#include "../hwMatrixS.h"
#include "../hwMathException.h"
#include "mkl_vml_functions.h"
#include "mkl_spblas.h"
#include "mkl_pardiso.h"

//------------------------------------------------------------------------------
// Returns cosine of input in radians
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::Cos(const hwMatrix& X,
                                        hwMatrix&       C)
{
    // return ElemFunc1<cos, hwComplex::cos, vdCos, vzCos>(eval, inputs, outputs);
    hwMathStatus status = C.Dimension(X.M(), X.N(), X.Type());

    if (!status.IsOk())
        throw hwMathException(status.GetMsgCode());

    int size = X.Size();

    if (X.IsReal())
    {
        const double* input = X.GetRealData();
        double* output = C.GetRealData();

        vdCos(size, input, output);
    }
    else
    {
        hwComplex* input = const_cast<hwComplex*> (X.GetComplexData());
        hwComplex* output = C.GetComplexData();
        MKL_Complex16* inputMKL = reinterpret_cast<MKL_Complex16*> (input);
        MKL_Complex16* outputMKL = reinterpret_cast<MKL_Complex16*> (output);

        vzCos(size, inputMKL, outputMKL);
    }
}
//------------------------------------------------------------------------------
// Returns cosine of input in radians
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::Cos(const hwMatrix& X)
{
    hwMatrix C;

    Cos(X, C);
    return C;
}
//------------------------------------------------------------------------------
// Returns sine of input in radians
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::Sin(const hwMatrix& X,
                                        hwMatrix&       S)
{
    // return ElemFunc1<sin, hwComplex::sin, vdSin, vzSin>(eval, inputs, outputs);
    hwMathStatus status = S.Dimension(X.M(), X.N(), X.Type());

    if (!status.IsOk())
        throw hwMathException(status.GetMsgCode());

    int size = X.Size();

    if (X.IsReal())
    {
        const double* input = X.GetRealData();
        double* output = S.GetRealData();

        vdSin(size, input, output);
    }
    else
    {
        hwComplex* input = const_cast<hwComplex*> (X.GetComplexData());
        hwComplex* output = S.GetComplexData();
        MKL_Complex16* inputMKL = reinterpret_cast<MKL_Complex16*> (input);
        MKL_Complex16* outputMKL = reinterpret_cast<MKL_Complex16*> (output);

        vzSin(size, inputMKL, outputMKL);
    }
}
//------------------------------------------------------------------------------
// Returns cosine of input in radians
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::Sin(const hwMatrix& X)
{
    hwMatrix S;

    Sin(X, S);
    return S;
}
//------------------------------------------------------------------------------
// Returns exponential of input
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::Exp(const hwMatrix& X,
                                        hwMatrix&       E)
{
    // return ElemFunc1<exp, hwComplex::exp, vdExp, vzExp>(eval, inputs, outputs);
    hwMathStatus status = E.Dimension(X.M(), X.N(), X.Type());

    if (!status.IsOk())
        throw hwMathException(status.GetMsgCode());

    int size = X.Size();

    if (X.IsReal())
    {
        const double* input = X.GetRealData();
        double* output = E.GetRealData();

        vdExp(size, input, output);
    }
    else
    {
        hwComplex* input = const_cast<hwComplex*> (X.GetComplexData());
        hwComplex* output = E.GetComplexData();
        MKL_Complex16* inputMKL = reinterpret_cast<MKL_Complex16*> (input);
        MKL_Complex16* outputMKL = reinterpret_cast<MKL_Complex16*> (output);

        vzExp(size, inputMKL, outputMKL);
    }
}
//------------------------------------------------------------------------------
// Returns exponential of input
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::Exp(const hwMatrix& X)
{
    hwMatrix E;

    Exp(X, E);
    return E;
}
//------------------------------------------------------------------------------
// Add two full matrices, sum = A + B
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::Add(const hwMatrix& A,
                                        const hwMatrix& B,
                                        hwMatrix&       sum)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (B.M() != m || B.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.IsReal() && B.IsReal())
    {
        const double* pA = A.GetRealData();
        const double* pB = B.GetRealData();
        hwMathStatus status = sum.Dimension(m, n, hwMatrix::REAL);
        double* pC = sum.GetRealData();
        vdAdd(sizeA, pA, pB, pC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pB = const_cast<hwComplex*> (B.GetComplexData());
        hwMathStatus status = sum.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pC = sum.GetComplexData();
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
        MKL_Complex16* pMKL_C = reinterpret_cast<MKL_Complex16*> (pC);
        vzAdd(sizeA, pMKL_A, pMKL_B, pMKL_C);
    }
    else if (!A.IsReal())
    {
        hwMatrix BC;
        BC.PackComplex(B);
        return Add(A, BC, sum);
    }
    else    // !B.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        return Add(AC, B, sum);
    }
}
//------------------------------------------------------------------------------
// Add two full matrices, sum = A + B
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::Add(const hwMatrix& A,
                                            const hwMatrix& B)
{
    hwMatrix sum;

    Add(A, B, sum);
    return sum;
}
//------------------------------------------------------------------------------
// Subtract two full matrices, diff = A - B
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::Subtr(const hwMatrix& A,
                                          const hwMatrix& B,
                                          hwMatrix&       diff)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (B.M() != m || B.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.IsReal() && B.IsReal())
    {
        const double* pA = A.GetRealData();
        const double* pB = B.GetRealData();
        hwMathStatus status = diff.Dimension(m, n, hwMatrix::REAL);
        double* pC = diff.GetRealData();
        vdSub(sizeA, pA, pB, pC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pB = const_cast<hwComplex*> (B.GetComplexData());
        hwMathStatus status = diff.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pC = diff.GetComplexData();
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
        MKL_Complex16* pMKL_C = reinterpret_cast<MKL_Complex16*> (pC);
        vzSub(sizeA, pMKL_A, pMKL_B, pMKL_C);
    }
    else if (!A.IsReal())
    {
        hwMatrix BC;
        BC.PackComplex(B);
        Subtr(A, BC, diff);
    }
    else    // !B.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        Subtr(AC, B, diff);
    }
}
//------------------------------------------------------------------------------
// Subtract two full matrices, diff = A - B
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::Subtr(const hwMatrix& A,
                                              const hwMatrix& B)
{
    hwMatrix diff;

    Subtr(A, B, diff);
    return diff;
}
//------------------------------------------------------------------------------
// Multiply two full matrices element-wise, prod = A .* B
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::MultByElems(const hwMatrix& A,
                                                const hwMatrix& B,
                                                hwMatrix&       prod)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (B.M() != m || B.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.IsReal() && B.IsReal())
    {
        const double* pA = A.GetRealData();
        const double* pB = B.GetRealData();
        hwMathStatus status = prod.Dimension(m, n, hwMatrix::REAL);
        double* pC = prod.GetRealData();
        vdMul(sizeA, pA, pB, pC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pB = const_cast<hwComplex*> (B.GetComplexData());
        hwMathStatus status = prod.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pC = prod.GetComplexData();
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
        MKL_Complex16* pMKL_C = reinterpret_cast<MKL_Complex16*> (pC);
        vzMul(sizeA, pMKL_A, pMKL_B, pMKL_C);
    }
    else if (!A.IsReal())
    {
        hwMatrix BC;
        BC.PackComplex(B);
        MultByElems(A, BC, prod);
    }
    else    // !B.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        MultByElems(AC, B, prod);
    }
}
//------------------------------------------------------------------------------
// Multiply two full matrices element-wise, prod = A .* B
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::MultByElems(const hwMatrix& A,
                                                    const hwMatrix& B)
{
    hwMatrix prod;

    MultByElems(A, B, prod);
    return prod;
}
//------------------------------------------------------------------------------
// Power operation applied element-wise (B = A .^ P)
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::PowerByElems(const hwMatrix& A,
                                                 const hwMatrix& P,
                                                 hwMatrix&       B)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (P.M() != m || P.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.IsReal() && P.IsReal())
    {
        int count = sizeA;
        const double* pA = A.GetRealData();
        const double* pP = P.GetRealData();
        bool A_gt_zero = true;

        while (count--)
        {
            if (*pA <= 0.0)
            {
                A_gt_zero = false;
                break;
            }

            ++pA;
        }

        pA = A.GetRealData();

        if (A_gt_zero)
        {
            hwMathStatus status = B.Dimension(m, n, hwMatrix::REAL);
            double* pB = B.GetRealData();
            vdPowr(sizeA, pA, pP, pB);
        }
        else
        {
            int count = sizeA;
            bool P_is_int = true;

            while (count--)
            {
                if (!IsInteger(*pP).IsOk())
                {
                    P_is_int = false;
                    break;
                }

                ++pP;
            }

            pP = P.GetRealData();

            if (P_is_int)
            {
                pA = A.GetRealData();
                hwMathStatus status = B.Dimension(m, n, hwMatrix::REAL);
                double* pB = B.GetRealData();
                vdPow(sizeA, pA, pP, pB);
            }
            else
            {
                hwMatrix AC;
                AC.PackComplex(A);
                hwMatrix PC;
                PC.PackComplex(P);
                PowerByElems(AC, PC, B);
            }
        }
    }
    else if (!A.IsReal() && !P.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pP = const_cast<hwComplex*> (P.GetComplexData());
        hwMathStatus status = B.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pB = B.GetComplexData();
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        MKL_Complex16* pMKL_P = reinterpret_cast<MKL_Complex16*> (pP);
        MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
        vzPow(sizeA, pMKL_A, pMKL_P, pMKL_B);
    }
    else if (!A.IsReal())
    {
        hwMatrix PC;
        PC.PackComplex(P);
        PowerByElems(A, PC, B);
    }
    else    // !P.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        PowerByElems(AC, P, B);
    }
}
//------------------------------------------------------------------------------
// Power operation applied element-wise (B = A .^ P)
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::PowerByElems(const hwMatrix&  A,
                                                 const hwComplex& P,
                                                 hwMatrix&        B)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (A.IsReal())
    {
        hwMatrix AC;
        AC.PackComplex(A);
        PowerByElems(AC, P, B);
    }
    else
    {
        if (P == 2.0)
        {
            MultByElems(A, A, B);
        }
        else
        {
            hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
            hwMathStatus status = B.Dimension(m, n, hwMatrix::COMPLEX);
            hwComplex* pB = B.GetComplexData();
            MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
            const MKL_Complex16* pMKL_P = reinterpret_cast<const MKL_Complex16*> (&P);
            MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
            vzPowx(sizeA, pMKL_A, *pMKL_P, pMKL_B);
        }
    }
}
//------------------------------------------------------------------------------
// Power operation applied element-wise (B = A .^ P)
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::PowerByElems(const hwMatrix& A,
                                                 double          P,
                                                 hwMatrix&       B)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (A.IsReal())
    {
        const double* pA = A.GetRealData();

        if (IsInteger(P).IsOk())
        {
            hwMathStatus status = B.Dimension(m, n, hwMatrix::REAL);
            double* pB = B.GetRealData();

            if (P == 2.0)
                vdSqr(sizeA, pA, pB);
            else
                vdPowx(sizeA, pA, P, pB);
        }
        else
        {
            int count = sizeA;
            bool A_ge_zero = true;

            while (count--)
            {
                if (*pA < 0.0)
                {
                    A_ge_zero = false;
                    break;
                }

                ++pA;
            }

            pA = A.GetRealData();

            if (A_ge_zero)
            {
                hwMathStatus status = B.Dimension(m, n, hwMatrix::REAL);
                double* pB = B.GetRealData();
                vdPowx(sizeA, pA, P, pB);
            }
            else
            {
                hwMatrix AC;
                AC.PackComplex(A);
                hwComplex PC(P, 0.0);
                PowerByElems(AC, PC, B);
            }
        }
    }
    else    // !A.IsReal()
    {
        hwComplex PC(P, 0.0);
        PowerByElems(A, PC, B);
    }
}
//------------------------------------------------------------------------------
// Power operation applied element-wise (B = A .^ P)
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::PowerByElems(const hwMatrix& A,
                                                     double          P)
{
    hwMatrix B;

    PowerByElems(A, P, B);
    return B;
}
//------------------------------------------------------------------------------
// Power operation applied element-wise (B = A .^ P)
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::PowerByElems(const hwMatrix& A,
                                                     const hwComplex& P)
{
    hwMatrix B;

    PowerByElems(A, P, B);
    return B;
}
//------------------------------------------------------------------------------
// Power operation applied element-wise (B = A .^ P)
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::PowerByElems(const hwMatrix& A,
                                                     const hwMatrix& P)
{
    hwMatrix B;

    PowerByElems(A, P, B);
    return B;
}
//------------------------------------------------------------------------------
// Divide two full matrices element-wise, prod = A ./ B
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::DivideByElems(const hwMatrix& A,
                                                  const hwMatrix& B,
                                                  hwMatrix&       quot)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (B.M() != m || B.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.IsReal() && B.IsReal())
    {
        const double* pA = A.GetRealData();
        const double* pB = B.GetRealData();
        hwMathStatus status = quot.Dimension(m, n, hwMatrix::REAL);
        double* pC = quot.GetRealData();
        vdDiv(sizeA, pA, pB, pC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pB = const_cast<hwComplex*> (B.GetComplexData());
        hwMathStatus status = quot.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pC = quot.GetComplexData();
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        MKL_Complex16* pMKL_B = reinterpret_cast<MKL_Complex16*> (pB);
        MKL_Complex16* pMKL_C = reinterpret_cast<MKL_Complex16*> (pC);
        vzDiv(sizeA, pMKL_A, pMKL_B, pMKL_C);
    }
    else if (!A.IsReal())
    {
        hwMatrix BC;
        BC.PackComplex(B);
        DivideByElems(A, BC, quot);
    }
    else    // !B.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        DivideByElems(AC, B, quot);
    }
}
//------------------------------------------------------------------------------
// Divide two full matrices element-wise, prod = A ./ B
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::DivideByElems(const hwMatrix& A,
                                        const hwMatrix& B)
{
    hwMatrix quot;
    
    DivideByElems(A, B, quot);
    return quot;
}
//------------------------------------------------------------------------------
// Rounds input toward -inf
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::Floor(const hwMatrix& A,
                                          hwMatrix&       F)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int size = A.Size();

    hwMathStatus status = F.Dimension(m, n, A.Type());

    if (A.IsReal())
    {
        const double* input = A.GetRealData();
        double* output = F.GetRealData();

        vdFloor(size, input, output);
    }
    else
    {
        hwComplex* input = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* output = F.GetComplexData();
        double* dbl_in = reinterpret_cast<double*> (input);
        double* dbl_out = reinterpret_cast<double*> (output);

        vdFloor(2 * size, dbl_in, dbl_out);
    }
}
//------------------------------------------------------------------------------
// Rounds input toward -inf
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::Floor(const hwMatrix& A)
{
    hwMatrix F;
    Floor(A, F);
    return F;
}
//------------------------------------------------------------------------------
// Computes modulo function
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::Mod(const hwMatrix& A,
                                        const hwMatrix& B,
                                        hwMatrix&       mod)
{
    if (!A.IsReal())
        throw hwMathException(HW_MATH_ERR_COMPLEX, 1);

    if (!B.IsReal())
        throw hwMathException(HW_MATH_ERR_COMPLEX, 2);

    if (A.M() != B.M() || A.N() != B.N())
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    hwMathStatus status = mod.Dimension(A.M(), A.N(), hwMatrix::REAL);
    // vdRemainder(size, A, B, mod);
    // vdRemainder does not follow our convention
    hwMatrix temp = Floor(DivideByElems(A, B));
    mod = Subtr(A, MultByElems(temp, B));
}
//------------------------------------------------------------------------------
// Computes modulo function
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::Mod(const hwMatrix& A,
                                            const hwMatrix& B)
{
    hwMatrix M;
    Mod(A, B, M);
    return M;
}
//------------------------------------------------------------------------------
// Computes modulo function
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::Mod(const hwMatrix& A,
                                        double          b,
                                        hwMatrix&       mod)
{
    hwMatrix temp(A.M(), A.N(), hwMatrix::REAL);
    temp.SetElements(b);

    Mod(A, temp, mod);
}
//------------------------------------------------------------------------------
// Computes modulo function
//------------------------------------------------------------------------------
template<>
inline hwMatrix MKL_Tutilities<double>::Mod(const hwMatrix& A,
                                            double          b)
{
    hwMatrix M;
    Mod(A, b, M);
    return M;
}
//------------------------------------------------------------------------------
// Conjugate of a complex matrix
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::Conj(const hwMatrix& A,
                                         hwMatrix& conj)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (!A.IsReal())
    {
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        MKL_Complex16* pMKL_A = reinterpret_cast<MKL_Complex16*> (pA);
        hwMathStatus status = conj.Dimension(m, n, hwMatrix::COMPLEX);
        hwComplex* pC = conj.GetComplexData();
        MKL_Complex16* pMKL_C = reinterpret_cast<MKL_Complex16*> (pC);
        vzConj(sizeA, pMKL_A, pMKL_C);
    }
    else
    {
        conj = A;
    }
}
//------------------------------------------------------------------------------
// Hypotenuse of paired matrix elements, hyp = sqrt(A^2 + B^2
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::Hypot(const hwMatrix& A,
                                          const hwMatrix& B,
                                          hwMatrix& hyp)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    int sizeA = A.Size();

    if (A.IsReal() && B.IsReal())
    {
        const double* pA = A.GetRealData();
        const double* pB = B.GetRealData();
        hwMathStatus status = hyp.Dimension(m, n, hwMatrix::REAL);
        double* pC = hyp.GetRealData();
        vdHypot(sizeA, pA, pB, pC);
    }
}
//------------------------------------------------------------------------------
// Add two sparse matrices, sum = A + B
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::SparseAdd(const hwMatrixS& A,
                                              const hwMatrixS& B,
                                              hwMatrixS&       sum)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    hwMathStatus status;

    if (B.M() != m || B.N() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.NNZ() == 0)
    {
        sum = B;
        return;
    }
    else if (B.NNZ() == 0)
    {
        sum = A;
        return;
    }

    MKL_INT* prA = const_cast<MKL_INT*> (A.rows());
    MKL_INT* pbA = const_cast<MKL_INT*> (A.pointerB());
    MKL_INT* peA = const_cast<MKL_INT*> (A.pointerE());
    MKL_INT* prB = const_cast<MKL_INT*> (B.rows());
    MKL_INT* pbB = const_cast<MKL_INT*> (B.pointerB());
    MKL_INT* peB = const_cast<MKL_INT*> (B.pointerE());

    MKL_INT* prC = nullptr;
    MKL_INT* pbC = nullptr;
    MKL_INT* peC = nullptr;

    // MKL matrix and description
    sparse_matrix_t A_MKL = NULL;
    sparse_matrix_t B_MKL = NULL;
    sparse_matrix_t C_MKL = NULL;
    sparse_status_t mkl_status;

    if (A.IsReal() && B.IsReal())
    {
        double* pVA = const_cast<double*> (A.GetRealData());
        double* pVB = const_cast<double*> (B.GetRealData());
        double* pVC = nullptr;

        mkl_status = mkl_sparse_d_create_csr(&A_MKL, SPARSE_INDEX_BASE_ZERO,
            n, m, pbA, peA, prA, pVA);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_d_create_csr(&B_MKL, SPARSE_INDEX_BASE_ZERO,
            n, m, pbB, peB, prB, pVB);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_d_add(SPARSE_OPERATION_NON_TRANSPOSE, A_MKL, 1.0, B_MKL, &C_MKL);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

        mkl_status = mkl_sparse_d_export_csr(C_MKL, &indexing, &n, &m, &pbC, &peC, &prC, &pVC);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sum = hwMatrixS(m, n, pbC, peC, prC, pVC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pVA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pVB = const_cast<hwComplex*> (B.GetComplexData());
        hwComplex* pVC = nullptr;
        struct _MKL_Complex16 one;

        one.real = 1.0;
        one.imag = 0.0;

        mkl_status = mkl_sparse_z_create_csr(&A_MKL, SPARSE_INDEX_BASE_ZERO,
                                             n, m, pbA, peA, prA,
                                             reinterpret_cast<MKL_Complex16*> (pVA));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_z_create_csr(&B_MKL, SPARSE_INDEX_BASE_ZERO,
                                             n, m, pbB, peB, prB,
                                             reinterpret_cast<MKL_Complex16*> (pVB));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_z_add(SPARSE_OPERATION_NON_TRANSPOSE, A_MKL, one, B_MKL, &C_MKL);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

        mkl_status = mkl_sparse_z_export_csr(C_MKL, &indexing, &n, &m, &pbC, &peC, &prC,
            reinterpret_cast<MKL_Complex16**> (&pVC));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sum = hwMatrixS(m, n, pbC, peC, prC, pVC);
    }
    else if (!A.IsReal())
    {
        hwMatrixS BC;
        BC.PackComplex(B);
        return SparseAdd(A, BC, sum);
    }
    else    // !B.IsReal()
    {
        hwMatrixS AC;
        AC.PackComplex(A);
        return SparseAdd(AC, B, sum);
    }

    mkl_sparse_destroy(A_MKL);
    mkl_sparse_destroy(B_MKL);
    mkl_sparse_destroy(C_MKL);
}
//------------------------------------------------------------------------------
// Multiply a sparse matrix by a full matrix, prod = A * B
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::SparseMult(const hwMatrixS& A,
                                               const hwMatrix&  B,
                                               hwMatrix&        prod)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    MKL_INT k = B.N();
    hwMathStatus status;

    if (B.M() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.NNZ() == 0 || B.IsEmpty())
    {
        status = prod.Dimension(m, k, hwMatrix::REAL);
        prod.SetElements(0.0);
        return;
    }

    MKL_INT* pr = const_cast<MKL_INT*> (A.rows());
    MKL_INT* pb = const_cast<MKL_INT*> (A.pointerB());
    MKL_INT* pe = const_cast<MKL_INT*> (A.pointerE());

    // MKL matrix and description
    sparse_matrix_t A_MKL = NULL;
    struct matrix_descr DSC;
    DSC.type = SPARSE_MATRIX_TYPE_GENERAL;
    sparse_status_t mkl_status;

    if (A.IsReal() && B.IsReal())
    {
        status = prod.Dimension(m, k, hwMatrix::REAL);

        // perform MKL operations
        double* pV = const_cast<double*> (A.GetRealData());

        if (k == 1)
        {
            // populate MKL matrix
            mkl_status = mkl_sparse_d_create_csc(&A_MKL, SPARSE_INDEX_BASE_ZERO,
                                                 m, n, pb, pe, pr, pV);

            if (mkl_status != SPARSE_STATUS_SUCCESS)
                throw hwMathException(HW_MATH_ERR_UNKNOWN);

            mkl_status = mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE,
                                         1.0, A_MKL, DSC, B.GetRealData(),
                                         0.0, prod.GetRealData());
        }
        else
        {
            // SPARSE_INDEX_BASE_ZERO is not supported in MKL for this
            // function, so convert to SPARSE_INDEX_BASE_ONE
            int nnz = A.NNZ();
            MKL_INT* pb1 = new MKL_INT[n];
            MKL_INT* pe1 = new MKL_INT[n];
            MKL_INT* pr1 = new MKL_INT[nnz];

            for (int i = 0; i < n; ++i)
            {
                pb1[i] = pb[i] + 1;
                pe1[i] = pe[i] + 1;
            }

            for (int i = 0; i < nnz; ++i)
                pr1[i] = pr[i] + 1;

            // populate MKL matrix
            mkl_status = mkl_sparse_d_create_csc(&A_MKL, SPARSE_INDEX_BASE_ONE,
                                                 m, n, pb1, pe1, pr1, pV);

            if (mkl_status != SPARSE_STATUS_SUCCESS)
                throw hwMathException(HW_MATH_ERR_UNKNOWN);

            mkl_status = mkl_sparse_d_mm(SPARSE_OPERATION_NON_TRANSPOSE,
                                         1.0, A_MKL, DSC, SPARSE_LAYOUT_COLUMN_MAJOR,
                                         B.GetRealData(), k, n, 0.0,
                                         prod.GetRealData(), m);

            delete[] pb1;
            delete[] pe1;
            delete[] pr1;
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        status = prod.Dimension(m, k, hwMatrix::COMPLEX);

        // perform MKL operations
        hwComplex* pV = const_cast<hwComplex*> (A.GetComplexData());
        hwTComplex<double>* pB = const_cast<hwComplex*> (B.GetComplexData());
        MKL_Complex16* pP = reinterpret_cast<MKL_Complex16*> (prod.GetComplexData());
        struct _MKL_Complex16 one;
        struct _MKL_Complex16 zero;

        one.real = 1.0;
        one.imag = 0.0;
        zero.real = 0.0;
        zero.imag = 0.0;

        if (k == 1)
        {
            // populate MKL matrix
            mkl_status = mkl_sparse_z_create_csc(&A_MKL, SPARSE_INDEX_BASE_ZERO,
                                                 m, n, pb, pe, pr,
                                                 reinterpret_cast<MKL_Complex16*> (pV));

            if (mkl_status != SPARSE_STATUS_SUCCESS)
                throw hwMathException(HW_MATH_ERR_UNKNOWN);

            mkl_status = mkl_sparse_z_mv(SPARSE_OPERATION_NON_TRANSPOSE,
                                         one, A_MKL, DSC,
                                         reinterpret_cast<MKL_Complex16*> (pB), zero, pP);
        }
        else
        {
            // SPARSE_INDEX_BASE_ZERO is not supported in MKL for this
            // function, so convert to SPARSE_INDEX_BASE_ONE
            int nnz = A.NNZ();
            MKL_INT* pb1 = new MKL_INT[n];
            MKL_INT* pe1 = new MKL_INT[n];
            MKL_INT* pr1 = new MKL_INT[nnz];

            for (int i = 0; i < n; ++i)
            {
                pb1[i] = pb[i] + 1;
                pe1[i] = pe[i] + 1;
            }

            for (int i = 0; i < nnz; ++i)
                pr1[i] = pr[i] + 1;

            // populate MKL matrix
            mkl_status = mkl_sparse_z_create_csc(&A_MKL, SPARSE_INDEX_BASE_ONE,
                                                 m, n, pb1, pe1, pr1,
                                                 reinterpret_cast<MKL_Complex16*> (pV));

            if (mkl_status != SPARSE_STATUS_SUCCESS)
                throw hwMathException(HW_MATH_ERR_UNKNOWN);

            mkl_status = mkl_sparse_z_mm(SPARSE_OPERATION_NON_TRANSPOSE,
                                         one, A_MKL, DSC, SPARSE_LAYOUT_COLUMN_MAJOR,
                                         reinterpret_cast<MKL_Complex16*> (pB), k, n,
                                         zero, pP, m);

            delete[] pb1;
            delete[] pe1;
            delete[] pr1;
        }
    }
    else if (!A.IsReal())
    {
        hwMatrix C;
        C.PackComplex(B);
        return SparseMult(A, C, prod);
    }
    else    // !B.IsReal()
    {
        hwMatrixS C;
        C.PackComplex(A);
        return SparseMult(C, B, prod);
    }

    if (mkl_status != SPARSE_STATUS_SUCCESS)
        throw hwMathException(HW_MATH_ERR_UNKNOWN);

    mkl_sparse_destroy(A_MKL);
}
//------------------------------------------------------------------------------
// Multiply a full matrix by a sparse matrix, prod = A * B
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::SparseMult(const hwMatrix&  A,
                                               const hwMatrixS& B,
                                               hwMatrix&        prod)
{
    // prod = A * B will be computed using trans(prod) = trans(B) * trans(A)
    // while interpreting prod and A as having row major storage

    // A: row major matrix
    MKL_INT rows_A = A.N();
    MKL_INT cols_A = A.M();

    // B: column major sparse matrix
    MKL_INT rows_B = B.M();
    MKL_INT cols_B = B.N();

    if (rows_A != rows_B)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    // prod: row major matrix
    MKL_INT rows_P = cols_B;
    MKL_INT cols_P = cols_A;

    if (B.NNZ() == 0 || A.IsEmpty())
    {
        hwMathStatus status = prod.Dimension(cols_P, rows_P, hwMatrix::REAL);
        prod.SetElements(0.0);
        return;
    }

    // compute prod
    MKL_INT* pr = const_cast<MKL_INT*> (B.rows());
    MKL_INT* pb = const_cast<MKL_INT*> (B.pointerB());
    MKL_INT* pe = const_cast<MKL_INT*> (B.pointerE());
    sparse_matrix_t B_MKL = NULL;
    struct matrix_descr DSC;
    DSC.type = SPARSE_MATRIX_TYPE_GENERAL;
    sparse_status_t mkl_status;

    MKL_INT columns = cols_P;
    MKL_INT ldx = cols_A;
    MKL_INT ldp = cols_P;

    if (A.IsReal() && B.IsReal())
    {
        hwMathStatus status = prod.Dimension(cols_P, rows_P, hwMatrix::REAL);

        // populate MKL matrix and compute
        double* pB = const_cast<double*> (B.GetRealData());

        mkl_status = mkl_sparse_d_create_csc(&B_MKL, SPARSE_INDEX_BASE_ZERO,
                                             rows_B, cols_B, pb, pe, pr, pB);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        double alpha = 1.0;
        double beta = 0.0;

        if (cols_A == 1)
        {
            mkl_status = mkl_sparse_d_mv(SPARSE_OPERATION_TRANSPOSE,
                                         alpha, B_MKL, DSC, A.GetRealData(),
                                         beta, prod.GetRealData());
        }
        else
        {
            mkl_status = mkl_sparse_d_mm(SPARSE_OPERATION_TRANSPOSE,
                                         alpha, B_MKL, DSC, SPARSE_LAYOUT_ROW_MAJOR,
                                         A.GetRealData(), columns, ldx, beta,
                                         prod.GetRealData(), ldp);
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwMathStatus status = prod.Dimension(cols_P, rows_P, hwMatrix::COMPLEX);

        // populate MKL matrix and compute
        hwComplex* pA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pB = const_cast<hwComplex*> (B.GetComplexData());
        MKL_Complex16* pP = reinterpret_cast<MKL_Complex16*> (prod.GetComplexData());

        // populate MKL matrix
        mkl_status = mkl_sparse_z_create_csc(&B_MKL, SPARSE_INDEX_BASE_ZERO,
                                             rows_B, cols_B, pb, pe, pr,
                                             reinterpret_cast<MKL_Complex16*> (pB));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        struct _MKL_Complex16 alpha;
        struct _MKL_Complex16 beta;

        alpha.real = 1.0;
        alpha.imag = 0.0;
        beta.real = 0.0;
        beta.imag = 0.0;

        if (cols_A == 1)
        {
            mkl_status = mkl_sparse_z_mv(SPARSE_OPERATION_TRANSPOSE,
                                         alpha, B_MKL, DSC, reinterpret_cast<MKL_Complex16*> (pA), beta, pP);
        }
        else
        {
            mkl_status = mkl_sparse_z_mm(SPARSE_OPERATION_TRANSPOSE,
                                         alpha, B_MKL, DSC, SPARSE_LAYOUT_ROW_MAJOR,
                                         reinterpret_cast<MKL_Complex16*> (pA), columns,
                                         ldx, beta, pP, ldp);
        }
    }
    else if (!A.IsReal())
    {
        hwMatrixS BC;
        BC.PackComplex(B);
        return SparseMult(A, BC, prod);
    }
    else    // !B.IsReal()
    {
        hwMatrix AC;
        AC.PackComplex(A);
        return SparseMult(AC, B, prod);
    }

    if (mkl_status != SPARSE_STATUS_SUCCESS)
        throw hwMathException(HW_MATH_ERR_UNKNOWN);

    mkl_sparse_destroy(B_MKL);
}
//------------------------------------------------------------------------------
// Multiply two sparse matrices, prod = A * B
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::SparseMult(const hwMatrixS& A,
                                               const hwMatrixS& B,
                                               hwMatrixS&       prod)
{
    // get dimensions info
    MKL_INT m = A.M();
    MKL_INT n = A.N();
    MKL_INT k = B.N();
    hwMathStatus status;

    if (B.M() != n)
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE);

    if (A.NNZ() == 0 || B.NNZ() == 0)
    {
        std::vector<int> ivec;
        std::vector<int> jvec;
        hwMatrix vals;
        prod = hwMatrixS(ivec, jvec, vals, m, k);
        return;
    }

    MKL_INT* prA = const_cast<MKL_INT*> (A.rows());
    MKL_INT* pbA = const_cast<MKL_INT*> (A.pointerB());
    MKL_INT* peA = const_cast<MKL_INT*> (A.pointerE());
    MKL_INT* prB = const_cast<MKL_INT*> (B.rows());
    MKL_INT* pbB = const_cast<MKL_INT*> (B.pointerB());
    MKL_INT* peB = const_cast<MKL_INT*> (B.pointerE());

    MKL_INT* prC = nullptr;
    MKL_INT* pbC = nullptr;
    MKL_INT* peC = nullptr;

    // MKL matrix and description
    sparse_matrix_t A_MKL = NULL;
    sparse_matrix_t B_MKL = NULL;
    sparse_matrix_t C_MKL = NULL;
    sparse_status_t mkl_status;

    if (A.IsReal() && B.IsReal())
    {
        double* pVA = const_cast<double*> (A.GetRealData());
        double* pVB = const_cast<double*> (B.GetRealData());
        double* pVC = nullptr;

        mkl_status = mkl_sparse_d_create_csc(&A_MKL, SPARSE_INDEX_BASE_ZERO,
                                             m, n, pbA, peA, prA, pVA);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_d_create_csc(&B_MKL, SPARSE_INDEX_BASE_ZERO,
                                             n, k, pbB, peB, prB, pVB);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, A_MKL, B_MKL, &C_MKL);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

        mkl_status = mkl_sparse_d_export_csc(C_MKL, &indexing, &n, &k, &pbC, &peC, &prC, &pVC);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        prod = hwMatrixS(n, k, pbC, peC, prC, pVC);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        hwComplex* pVA = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* pVB = const_cast<hwComplex*> (B.GetComplexData());
        hwComplex* pVC = nullptr;

        mkl_status = mkl_sparse_z_create_csc(&A_MKL, SPARSE_INDEX_BASE_ZERO,
                                             m, n, pbA, peA, prA,
                                             reinterpret_cast<MKL_Complex16*> (pVA));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_z_create_csc(&B_MKL, SPARSE_INDEX_BASE_ZERO,
                                             n, k, pbB, peB, prB, reinterpret_cast<MKL_Complex16*> (pVB));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        mkl_status = mkl_sparse_spmm(SPARSE_OPERATION_NON_TRANSPOSE, A_MKL, B_MKL, &C_MKL);

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        sparse_index_base_t indexing = SPARSE_INDEX_BASE_ZERO;

        mkl_status = mkl_sparse_z_export_csc(C_MKL, &indexing, &n, &k, &pbC, &peC, &prC,
                                             reinterpret_cast<MKL_Complex16**> (&pVC));

        if (mkl_status != SPARSE_STATUS_SUCCESS)
            throw hwMathException(HW_MATH_ERR_UNKNOWN);

        prod = hwMatrixS(n, k, pbC, peC, prC, pVC);
    }
    else if (!A.IsReal())
    {
        hwMatrixS BC;
        BC.PackComplex(B);
        return SparseMult(A, BC, prod);
    }
    else    // !B.IsReal()
    {
        hwMatrixS AC;
        AC.PackComplex(A);
        return SparseMult(AC, B, prod);
    }

    mkl_sparse_destroy(A_MKL);
    mkl_sparse_destroy(B_MKL);
    mkl_sparse_destroy(C_MKL);
}

//------------------------------------------------------------------------------
// Set MKL iparam9
//------------------------------------------------------------------------------
static int iparm9 = 13;

template<>
inline void MKL_Tutilities<double>::SetIparam9(int _iparm9)
{
    iparm9 = _iparm9;
}
//------------------------------------------------------------------------------
// Divide a full matrix by a sparse matrix on the left side, Q = A \ B
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::SparseDivideLeft(const hwMatrixS& A,
                                                     const hwMatrix&  B,
                                                     hwMatrix&        Q)
{
    int nRows = A.M();
    int nCols = A.N();
    int n = B.N();

    if (nRows != nCols)
        throw hwMathException(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (nRows != B.M())
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (B.IsEmpty())
    {
        hwMathStatus status = Q.Dimension(nCols, n, hwMatrix::REAL);
        Q.SetElements(0.0);
        return;
    }

    if (A.NNZ() == 0)
    {
        throw hwMathException(HW_MATH_ERR_SINGMATRIX);
    }

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void* pt[64];
    // Pardiso control parameters.
    MKL_INT iparm[64];

    // Setup Pardiso control parameters
    for (MKL_INT i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;        // No solver default
    iparm[1] = 2;        // Fill-in reordering from METIS
    iparm[3] = 0;        // No iterative-direct algorithm
    iparm[4] = 0;        // No user fill-in reducing permutation
    iparm[5] = 0;        // Write solution into x
    iparm[6] = 0;        // Not in use
    iparm[7] = 2;        // Max numbers of iterative refinement steps
    iparm[8] = 0;        // Not in use
    iparm[9] = iparm9;   // Perturb the pivot elements with 1E-13
    iparm[10] = 1;        // Use nonsymmetric permutation and scaling MPS
    iparm[11] = 2;        // Transpose solve for systems in CSC format
    iparm[12] = 1;        // Maximum weighted matching algorithm is switched-on (default for non-symmetric)
    iparm[13] = 0;        // Output: Number of perturbed pivots
    iparm[14] = 0;        // Not in use
    iparm[15] = 0;        // Not in use
    iparm[16] = 0;        // Not in use
    iparm[17] = 0;        // Output: Number of nonzeros in the factor LU (default: -1)
    iparm[18] = 0;        // Output: Mflops for LU factorization (default: -1)
    iparm[19] = 0;        // Output: Numbers of CG Iterations
    iparm[34] = 1;        // Zero based indexing

    MKL_INT maxfct = 1;   // Maximum number of numerical factorizations.
    MKL_INT mnum = 1;   // Which factorization to use.
    MKL_INT msglvl = 0;   // Do not print statistical information in file
    MKL_INT error = 0;   // Initialize error flag
    MKL_INT phase;
    double ddum;          // Auxiliary double dummy
    MKL_INT idum;         // Auxiliary integer dummy

    // Initialize the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver.
    for (MKL_INT i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }

    // define the non-zero structure of the matrix
    int* pr = const_cast<int*> (A.rows());
    int* pe = const_cast<int*> (A.pointerE());

    std::vector<int> colCount(nCols + 1);  // TODO: rework, merging pointerB and pointerE

    for (int ii = 0; ii < nCols; ++ii)
        colCount[ii + 1] = pe[ii];

    MKL_INT* ia = colCount.data();
    MKL_INT* ja = pr;

    // choose matrix case
    if (A.IsReal() && B.IsReal())
    {
        MKL_INT mtype = 11;   // Real unsymmetric matrix

        hwMathStatus status = Q.Dimension(nCols, n, hwMatrix::REAL);

        if (!status.IsOk())
            throw hwMathException(status.GetMsgCode());

        double* a = const_cast<double*> (A.GetRealData());
        double* b = const_cast<double*> (B.GetRealData());
        double* x = Q.GetRealData();

        // Reordering and Symbolic Factorization. This step also allocates
        // all memory that is necessary for the factorization
        phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &n, iparm, &msglvl,
                &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during symbolic factorization: %d", error);
        }

        // Numerical factorization
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &n, iparm, &msglvl,
                &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during numerical factorization: %d", error);
        }

        // Solution phase
        phase = 33;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
            &nRows, a, ia, ja, &idum, &n, iparm, &msglvl, b, x, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during solution: %d", error);
        }

        // Termination and release of memory
        phase = -1;

        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, &ddum, ia, ja, &idum, &n,
                iparm, &msglvl, &ddum, &ddum, &error);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        MKL_INT mtype = 13;   // Complex unsymmetric matrix

        hwMathStatus status = Q.Dimension(nCols, n, hwMatrix::COMPLEX);

        if (!status.IsOk())
            throw hwMathException(status.GetMsgCode());

        hwComplex* a = const_cast<hwComplex*> (A.GetComplexData());
        hwComplex* b = const_cast<hwComplex*> (B.GetComplexData());
        hwComplex* x = Q.GetComplexData();

        // Reordering and Symbolic Factorization. This step also allocates
        // all memory that is necessary for the factorization
        phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &n, iparm, &msglvl,
                &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during symbolic factorization: %d", error);
        }

        // Numerical factorization
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &n, iparm, &msglvl,
                &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during numerical factorization: %d", error);
        }

        // Solution phase
        phase = 33;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &n, iparm, &msglvl, b, x, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &n,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during solution: %d", error);
        }

        // Termination and release of memory
        phase = -1;

        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, &ddum, ia, ja, &idum, &n,
                iparm, &msglvl, &ddum, &ddum, &error);
    }
    else if (!A.IsReal())
    {
        hwMatrix C;
        C.PackComplex(B);
        return SparseDivideLeft(A, C, Q);
    }
    else    // !B.IsReal()
    {
        hwMatrixS C;
        C.PackComplex(A);
        return SparseDivideLeft(C, B, Q);
    }
}
//------------------------------------------------------------------------------
// Divide a full matrix by a sparse matrix on the right side, Q = A / B
//------------------------------------------------------------------------------
template<>
inline void MKL_Tutilities<double>::SparseDivideRight(const hwMatrix&   A,
                                                      const hwMatrixS& B,
                                                      hwMatrix&        Q)
{
    int nRows = B.M();
    int nCols = B.N();
    int m = A.M();

    if (nRows != nCols)
        throw hwMathException(HW_MATH_ERR_MTXNOTSQUARE, 1);

    if (nCols != A.N())
        throw hwMathException(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    if (A.IsEmpty())
    {
        hwMathStatus status = Q.Dimension(m, nRows, hwMatrix::REAL);
        Q.SetElements(0.0);
        return;
    }

    if (B.NNZ() == 0)
    {
        throw hwMathException(HW_MATH_ERR_SINGMATRIX);
    }

    // Internal solver memory pointer pt,
    // 32-bit: int pt[64]; 64-bit: long int pt[64]
    // or void *pt[64] should be OK on both architectures
    void* pt[64];
    // Pardiso control parameters.
    MKL_INT iparm[64];

    // Setup Pardiso control parameters
    for (MKL_INT i = 0; i < 64; i++)
    {
        iparm[i] = 0;
    }
    iparm[0] = 1;        // No solver default
    iparm[1] = 2;        // Fill-in reordering from METIS
    iparm[3] = 0;        // No iterative-direct algorithm
    iparm[4] = 0;        // No user fill-in reducing permutation
    iparm[5] = 0;        // Write solution into x
    iparm[6] = 0;        // Not in use
    iparm[7] = 2;        // Max numbers of iterative refinement steps
    iparm[8] = 0;        // Not in use
    iparm[9] = iparm9;   // Perturb the pivot elements with 1E-13
    iparm[10] = 1;       // Use nonsymmetric permutation and scaling MPS
    iparm[11] = 0;       // Non-transpose solve for systems in CSC format
    iparm[12] = 1;       // Maximum weighted matching algorithm is switched-on (default for non-symmetric)
    iparm[13] = 0;       // Output: Number of perturbed pivots
    iparm[14] = 0;       // Not in use
    iparm[15] = 0;       // Not in use
    iparm[16] = 0;       // Not in use
    iparm[17] = 0;       // Output: Number of nonzeros in the factor LU (default: -1)
    iparm[18] = 0;       // Output: Mflops for LU factorization (default: -1)
    iparm[19] = 0;       // Output: Numbers of CG Iterations
    iparm[34] = 1;       // Zero based indexing

    MKL_INT maxfct = 1;  // Maximum number of numerical factorizations.
    MKL_INT mnum = 1;    // Which factorization to use.
    MKL_INT msglvl = 0;  // Do not print statistical information in file
    MKL_INT error = 0;   // Initialize error flag
    MKL_INT phase;
    double ddum;         // Auxiliary double dummy
    MKL_INT idum;        // Auxiliary integer dummy

    // Initialize the internal solver memory pointer. This is only
    // necessary for the FIRST call of the PARDISO solver.
    for (MKL_INT i = 0; i < 64; i++)
    {
        pt[i] = 0;
    }

    // define the non-zero structure of the matrix
    int* pr = const_cast<int*> (B.rows());
    int* pe = const_cast<int*> (B.pointerE());

    std::vector<int> colCount(nCols + 1);  // TODO: rework, merging pointerB and pointerE

    for (int ii = 0; ii < nCols; ++ii)
        colCount[ii + 1] = pe[ii];

    MKL_INT* ia = colCount.data();
    MKL_INT* ja = pr;

    hwMatrix AT;

    hwMathStatus status = AT.Transpose(A);

    if (!status.IsOk())
        throw hwMathException(status.GetMsgCode());

    // choose matrix case
    if (A.IsReal() && B.IsReal())
    {
        MKL_INT mtype = 11;   // Real unsymmetric matrix

        status = Q.Dimension(nRows, m, hwMatrix::REAL);

        if (!status.IsOk())
            throw hwMathException(status.GetMsgCode());

        double* a = const_cast<double*> (B.GetRealData());
        double* b = const_cast<double*> (AT.GetRealData());
        double* x = Q.GetRealData();

        // Reordering and Symbolic Factorization. This step also allocates
        // all memory that is necessary for the factorization
        phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl,
                &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during symbolic factorization: %d", error);
        }

        // Numerical factorization
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl,
                &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during numerical factorization: %d", error);
        }

        // Solution phase
        phase = 33;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl, b, x, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during solution: %d", error);
        }

        status = Q.Transpose();

        // Termination and release of memory
        phase = -1;

        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, &ddum, ia, ja, &idum, &m,
                iparm, &msglvl, &ddum, &ddum, &error);
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        MKL_INT mtype = 13;   // Complex unsymmetric matrix

        status = Q.Dimension(nRows, m, hwMatrix::COMPLEX);

        if (!status.IsOk())
            throw hwMathException(status.GetMsgCode());

        hwComplex* a = const_cast<hwComplex*> (B.GetComplexData());
        hwComplex* b = const_cast<hwComplex*> (AT.GetComplexData());
        hwComplex* x = Q.GetComplexData();

        // Reordering and Symbolic Factorization. This step also allocates
        // all memory that is necessary for the factorization
        phase = 11;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl,
                &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during symbolic factorization: %d", error);
        }

        // Numerical factorization
        phase = 22;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl,
                &ddum, &ddum, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during numerical factorization: %d", error);
        }

        // Solution phase
        phase = 33;
        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, a, ia, ja, &idum, &m, iparm, &msglvl, b, x, &error);

        if (error != 0)
        {
            phase = -1;
            PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                    &nRows, &ddum, ia, ja, &idum, &m,
                    iparm, &msglvl, &ddum, &ddum, &error);

            throw hwMathException(HW_MATH_ERR_UNKNOWN);
            // printf("\nERROR during solution: %d", error);
        }

        status = Q.Transpose();

        // Termination and release of memory
        phase = -1;

        PARDISO(pt, &maxfct, &mnum, &mtype, &phase,
                &nRows, &ddum, ia, ja, &idum, &m,
                iparm, &msglvl, &ddum, &ddum, &error);
    }
    else if (!A.IsReal())
    {
        hwMatrixS C;
        C.PackComplex(B);
        return SparseDivideRight(A, C, Q);
    }
    else    // !B.IsReal()
    {
        hwMatrix C;
        C.PackComplex(A);
        return SparseDivideRight(C, B, Q);
    }
}
