/**
* @file MatrixNUtils.h
* @date January 2018
* Copyright (C) 2018 Altair Engineering, Inc.  
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

#ifndef __MATRIX_N_UTILS_H__
#define __MATRIX_N_UTILS_H__

// Begin defines/includes
#include "Hml2Dll.h"

#include "BuiltInFuncsUtils.h"
#include "Currency.h"
#include "EvaluatorInt.h"

typedef bool (*OML_func1)(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);
typedef Currency (ExprTreeEvaluator::*OML_func2)(const Currency& op);
typedef Currency (ExprTreeEvaluator::*OML_func3)(const Currency& op1, const Currency& op2);
typedef Currency (ExprTreeEvaluator::*OML_func4)(const Currency& op1, const Currency& op2, int op);
typedef bool (*OML_func5)(const Currency& op1, const Currency& op2);
typedef bool (*OML_func6)(const Currency& op1, const Currency& op2, const Currency& tol);

// Apply a function to each element of an ND matrix. The syntax is oml_func(ND).
// Examples: cos(ND), assert(ND), rat(ND, tol)
HML2DLL_DECLS
bool oml_MatrixNUtil1(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                      const OML_func1 oml_func);

// Apply a function to each element pair from the arguments, with at least one ND matrix. The syntax is oml_func(ND1, ND2).
// Example: atan2(ND1, ND2)
bool oml_MatrixNUtil2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                      const OML_func1 oml_func);

// Apply a function to each vector of an ND matrix in the specified dimension. The
// syntax is oml_func(ND, dim). This utility applies to functions for which the
// output for each vector input is a scalar.
// Examples: sum(ND, dim), max(ND, [], dim)
HML2DLL_DECLS
bool oml_MatrixNUtil3(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                      const OML_func1 oml_func, int dimArg = 0);

// Apply a function to each vector of an ND matrix in the specified dimension. The
// syntax is oml_func(ND, dim). This utility applies to functions for which the
// output for each vector input is a vector.
// Examples: cumsum(ND, dim), sort(ND, dim, mode)
HML2DLL_DECLS
bool oml_MatrixNUtil4(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                      const OML_func1 oml_func, int dimArg = 0, int ndArg = 1);

// Apply a function to an argument list that contains vectors stored as ND matrices. The
// syntax is oml_func(..., ND, ...). Output vectors will be in the same dimension as the first
// argument if it is a vector.
// Examples: not yet in use
HML2DLL_DECLS
bool oml_MatrixNVecs(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                     const OML_func1 oml_func);

// ND support for dot and cross functions
bool oml_MatrixN_VecProd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                         const OML_func1 oml_func, int vecLength);

// ND support for sum function
bool oml_MatrixN_sum(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                     int dimArg = 0);

// ND support for cumsum function
bool oml_MatrixN_cumsum(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                        int dimArg = 0);

// ND support for diff function
bool oml_MatrixN_diff(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                      const OML_func1 oml_func);

// ND support for circshift function
bool oml_MatrixN_circshift(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                           const OML_func1 oml_func);

// Apply an operator function to each element of an ND matrix.
// Examples: -ND, not(ND)
Currency oml_MatrixNUtil6(const Currency& op, OML_func2 oml_func);

// Apply an arithmetic operator function to each element pair from the arguments, with at least one ND matrix.
// Examples: ND1 + ND2, times(ND1, ND2)
Currency oml_MatrixNUtil7(const Currency& op1, const Currency& op2, OML_func3 oml_func);

// Apply a comparison operator function to each element pair from the arguments, with at least one ND matrix.
// Examples: ND1 < ND2, eq(ND, ND2)
Currency oml_MatrixNUtil8(const Currency& op1, const Currency& op2, int op, OML_func4 oml_func);

// Apply a bool function to a pair of arguments, with at least one ND matrix.
// Examples: assert(ND1, ND2), isequal(ND1, ND2)
bool oml_MatrixNUtil9(const Currency& op1, const Currency& op2, OML_func5 oml_func);

// Apply a bool function to a pair of arguments, with at least one ND matrix and a tolerance.
// Example: isequal(ND1, ND2, tol)
bool oml_MatrixNUtil10(const Currency& op1, const Currency& op2, const Currency& tol, OML_func6 oml_func);

// Apply a function to each vector of a 2D matrix in the specified dimension. The
// syntax is oml_func(2D, dim). This utility applies to functions for which the
// output for each vector is a scalar.
// Examples: any(2D, dim), all(2D, dim)
HML2DLL_DECLS
Currency oml_MatrixUtil(EvaluatorInterface& eval, const hwMatrix* mtx, int dim, double (*func)(EvaluatorInterface&, const hwMatrix*));

HML2DLL_DECLS
void checkMathStatus(EvaluatorInterface& eval, hwMathStatus stat);  // also declared in BuiltInFuncs.h, review later

// 2D support for sort function
template <typename T1, typename T2>
inline Currency oml_Matrix_sort(EvaluatorInterface& eval, const hwTMatrix<T1, T2>* mtx, int dim,
                    hwTMatrix<T1, T2>* (*func)(EvaluatorInterface&, const hwTMatrix<T1, T2>*, std::pair<int&, hwMatrix*>*),
                    std::pair<int&, hwMatrix*>* index_data)
{
    // todo: evolve function to be more parallel to oml_MatrixNUtil4, at least the signature
    int i;
    hwTMatrix<T1, T2>* result;
    if (dim == 1)
    {
        hwTMatrix<T1, T2>* col = new hwTMatrix<T1, T2>;
        Currency colcur(col);
        for (i = 0; i < mtx->N(); ++i)
        {
            hwMathStatus stat = mtx->ReadColumn(i, *col);
            if (!stat.IsOk())
            {
                if (i)
                    delete result;
                BuiltInFuncsUtils::CheckMathStatus(eval, stat);
            }

            hwTMatrix<T1, T2>* temp = (*func)(eval, col, index_data);
            if (!i)
                result = new hwTMatrix<T1, T2>(temp->Size(), mtx->N(), hwTMatrix<T1, T2>::REAL);

            writeCol(eval, result, temp, i);
            delete temp;
        }
        return i ? result : new hwTMatrix<T1, T2>(0, 0, hwTMatrix<T1, T2>::REAL);
    }
    else // dim == 2
    {
        hwTMatrix<T1, T2>* row = new hwTMatrix<T1, T2>;
        Currency rowcur(row);
        for (i = 0; i < mtx->M(); ++i)
        {
            hwMathStatus stat = mtx->ReadRow(i, *row);
            if (!stat.IsOk())
            {
                if (i)
                    delete result;
                BuiltInFuncsUtils::CheckMathStatus(eval, stat);
            }

            hwTMatrix<T1, T2>* temp = (*func)(eval, row, index_data);
            if (!i)
                result = new hwTMatrix<T1, T2>(mtx->M(), temp->Size(), hwTMatrix<T1, T2>::REAL);

            writeRow(eval, result, temp, i);
            delete temp;
        }
        return i ? result : new hwTMatrix<T1, T2>(0, 0, hwTMatrix<T1, T2>::REAL);
    }
}

#endif // __MATRIX_N_UTILS_H__
