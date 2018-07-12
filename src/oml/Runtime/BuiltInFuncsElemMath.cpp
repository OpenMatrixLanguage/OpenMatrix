/**
* @file BuiltInFuncsElemMath.cpp
* @date July 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language (“OpenMatrix”) software.
* Open Source License Information:
* OpenMatrix is free software. You can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* OpenMatrix is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
* You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
* 
* Commercial License Information: 
* For a copy of the commercial license terms and conditions, contact the Altair Legal Department at Legal@altair.com and in the subject line, use the following wording: Request for Commercial License Terms for OpenMatrix.
* Altair’s dual-license business model allows companies, individuals, and organizations to create proprietary derivative works of OpenMatrix and distribute them - whether embedded or bundled with other software - under a commercial license agreement.
* Use of Altair’s trademarks and logos is subject to Altair's trademark licensing policies.  To request a copy, email Legal@altair.com and in the subject line, enter: Request copy of trademark and logo usage policy.
*/

// Begin defines/includes
#include "BuiltInFuncsElemMath.h"

#include <cassert>
#include <climits>
#include <iomanip>
#include <cmath>

#include "BuiltInFuncsUtils.h"
#include "OML_Error.h"
#include "MatrixNUtils.h"

#include "hwMatrix.h"
#include "hwMatrixN.h"
#include "math/kernel/GeneralFuncs.h"

typedef hwTMatrix< double,  hwTComplex<double> > hwMatrix;
typedef hwTMatrixN<double, hwTComplex<double> > hwMatrixN;

// End defines/includes

//------------------------------------------------------------------------------
// True if given complex numbers are comparable
//------------------------------------------------------------------------------
static bool complexLessThanTol(const hwComplex &cplx1, const hwComplex &cplx2, double tol)
{
    double comp = cplx1.Mag() - cplx2.Mag();
    if (comp < -tol)
        return true;
    else if (comp > tol)
        return false;

    comp = cplx1.Arg() - cplx2.Arg();
    
    return comp < -tol;
}
//------------------------------------------------------------------------------
// True if first complex number is less than second
//------------------------------------------------------------------------------
static bool complexLessThan(const hwComplex &cplx1, const hwComplex &cplx2)
{
    return complexLessThanTol(cplx1, cplx2, 0.0);
}
//------------------------------------------------------------------------------
// Helper function for Unique, gets a deque of uniue doubles
//------------------------------------------------------------------------------
std::deque<double> BuiltInFuncsElemMath::UniqueHelperRealMtxN(const hwMatrixN* mtx)
{
    if (!mtx ) return std::deque<double>(); // Nothing to sort
    assert(mtx->IsReal());

    bool hasNegInf = false;
    bool hasPosInf = false;
    bool hasNan    = false;
    int  numElem   = mtx->Size();

    std::deque<double> vals;
    for (int i = 0; i < numElem; ++i)
    {
        double val = (*mtx)(i);
        // Check for valid field width values
        if (IsNegInf_T (val))  
            hasNegInf = true;
        else if (IsInf_T(val))
            hasPosInf = true;
        else if (IsNaN_T(val))
            hasNan = true;
        else
            vals.push_back(val);
    }

    if (!vals.empty() && vals.size() > 1)
    {
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    }

    if (hasNegInf)
        vals.insert(vals.begin(), std::numeric_limits<double>::infinity() * -1);

    if (hasPosInf)
        vals.push_back(std::numeric_limits<double>::infinity());

    if (hasNan)
        vals.push_back(std::numeric_limits<double>::quiet_NaN());

    return vals;
}
//------------------------------------------------------------------------------
// Helper function for Unique, gets a deque of uniue complex numbers
//------------------------------------------------------------------------------
std::deque<hwComplex> BuiltInFuncsElemMath::UniqueHelperComplexMtxN(const hwMatrixN* mtx)
{
    if (!mtx ) return std::deque<hwComplex>(); // Nothing to sort
    assert(!mtx->IsReal());

    bool hasNegInf = false;
    bool hasPosInf = false;
    bool hasNan    = false;
    int  numElem   = mtx->Size();

    std::deque<hwComplex> vals;
    for (int i = 0; i < numElem; ++i)
    {
        hwComplex val = mtx->z(i);
        double    rval = val.Real();
        double    ival = val.Imag();

        // Check for valid field width values
        if (IsNegInf_T (rval) || IsNegInf_T (ival))
            hasNegInf = true;
        else if (IsInf_T(rval) || IsInf_T (ival))
            hasPosInf = true;
        else if (IsNaN_T(rval) || IsNaN_T(ival))
            hasNan = true;
        else
            vals.push_back(val);
    }

    if (!vals.empty() && vals.size() > 1)
    {
        std::sort(vals.begin(), vals.end(), &complexLessThan);
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    }

    if (hasNegInf)
        vals.insert(vals.begin(), 
                    hwComplex(std::numeric_limits<double>::infinity()*-1, 0));

    if (hasPosInf)
        vals.push_back(hwComplex(std::numeric_limits<double>::infinity(), 0));

    if (hasNan)
        vals.push_back(hwComplex(std::numeric_limits<double>::quiet_NaN(), 0));
    return vals;
}
//------------------------------------------------------------------------------
// Helper function for unique command for ND matrices
//------------------------------------------------------------------------------
void BuiltInFuncsElemMath::UniqueHelperFuncMtxN(EvaluatorInterface&    eval,
                                                const Currency&        in,
                                                bool                   forward,
                                                bool                   outputIdx,
                                                bool                   inputIdx,
                                                std::vector<Currency>& outputs)
{
    const hwMatrixN* x = in.MatrixN();
    if (!x || x->Size() <= 0)
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        if (outputIdx)
            outputs.push_back(EvaluatorInterface::allocateMatrix());
        if (inputIdx)
            outputs.push_back(EvaluatorInterface::allocateMatrix());
        return;
    }
    int mtxSize = x->Size();

    if (mtxSize == 1)
    {
        outputs.push_back(in);
        if (outputIdx)
            outputs.push_back(1);
        if (inputIdx)
            outputs.push_back(1);
        return;
    }
    BuiltInFuncsElemMath funcs;
    bool singlerow = funcs.IsSingleRowND(x);

    if (x->IsReal())
    {
        std::deque<double> y (funcs.UniqueHelperRealMtxN(x));

        outputs.push_back(BuiltInFuncsUtils::ContainerToMatrix(y, singlerow));
        if (outputIdx)
            outputs.push_back(funcs.GetMatrixNIndices(y, x, forward));

        if (inputIdx)
            outputs.push_back(funcs.GetValueIndicesND(y, x));
    }
    else
    {
        std::deque<hwComplex> y (funcs.UniqueHelperComplexMtxN(x));
        outputs.push_back(BuiltInFuncsUtils::ContainerToMatrix(y, singlerow));
        if (outputIdx)
            outputs.push_back(funcs.GetMatrixNIndices(x, y, forward));

        if (inputIdx)
            outputs.push_back(funcs.GetValueIndicesND(x, y));
    }
}
//------------------------------------------------------------------------------
// Gets indices of occurences of unsorted values in a real matrix
//------------------------------------------------------------------------------
Currency BuiltInFuncsElemMath::GetMatrixNIndices(const std::deque<double>& vals, 
                                                 const hwMatrixN*          in, 
                                                 bool                      forward)
{
    assert(in);
    if (!in || in->Size() <= 0) 
        return EvaluatorInterface::allocateMatrix(0, 0, hwMatrix::REAL);

    std::vector<int> dims = in->Dimensions();

    int mtxsize = in->Size();
    int valsize = static_cast<int>(vals.size());
    int rows    = valsize;
    int cols    = 1;
    if (!dims.empty() && dims[0] == 1)
    {
        rows = 1;
        cols = valsize;
    }

    hwMatrix* indices = EvaluatorInterface::allocateMatrix(rows, cols, hwMatrix::REAL);

    for (int i = 0; i < valsize; ++i)
    {
        double tofind   = vals[i];
        bool   isNegInf = false; // Handle -Inf/Inf/Nan
        bool   isPosInf = false;
        bool   isNan    = false;
            
        if (IsNegInf_T (tofind))
            isNegInf = true;
        else if (IsInf_T (tofind))
            isPosInf = true;
        else if (IsNaN_T(tofind))
            isNan = true;

        int index = forward ? 0 : mtxsize - 1;
        while (1)
        {
            if ((forward  && index > mtxsize - 1) ||
                (!forward && index < 0))
                break;

            double val = (*in)(index);
            if ((isNegInf && IsNegInf_T (val)) ||
                (isPosInf && IsInf_T(val))     ||
                (isNan    && IsNaN_T(val))     ||
                (AreEqual(val, tofind)))
            {
                (*indices)(i) = index + 1;
                break;
            }
            if (forward)
                index++;
            else
                index--;
        }
    }
    return Currency(indices);
}
//------------------------------------------------------------------------------
// Gets indices of occurences of sorted values in a real matrix
//------------------------------------------------------------------------------
Currency BuiltInFuncsElemMath::GetValueIndicesND(const std::deque<double>& vals, 
                                                 const hwMatrixN*          in)
{
    if (!in || vals.empty()) 
        return EvaluatorInterface::allocateMatrix(0, 0, hwMatrix::REAL);
        
    size_t valsize = vals.size();

    int mtxsize = in->Size();
    int rows    = mtxsize;
    int cols    = 1;
    if (IsSingleRowND(in))
    {
        rows = 1;
        cols = mtxsize;
    }

    hwMatrix* indices = EvaluatorInterface::allocateMatrix(rows, cols, hwMatrix::REAL);

    for (int i = 0; i < mtxsize; ++i)
    {
        double ref      = (*in)(i);
        bool   isNegInf = false; // Handle -Inf/Inf/Nan
        bool   isPosInf = false;
        bool   isNan    = false;
            
        if (IsNegInf_T (ref))
            isNegInf = true;
        else if (IsInf_T (ref))
            isPosInf = true;
        else if (IsNaN_T(ref))
            isNan = true;

        for (size_t j = 0; j < valsize; ++j)
        {
            double val = vals[j];
            if ((isNegInf && IsNegInf_T (val)) ||
                (isPosInf && IsInf_T(val))     ||
                (isNan    && IsNaN_T(val))     ||
                 AreEqual(val, ref))
            {
                (*indices)(i) = static_cast<int>(j) + 1;
                break;
            }
        }
    }
    return Currency(indices);
}
//------------------------------------------------------------------------------
// Gets indices of occurences of values in a complex matrix => x = y(j)
//------------------------------------------------------------------------------
Currency BuiltInFuncsElemMath::GetValueIndicesND(const hwMatrixN*             x,
                                                 const std::deque<hwComplex>& y)
{
    if (!x || x->Size() <= 0 || y.empty()) return EvaluatorInterface::allocateMatrix();

    int matsize = x->Size();
    int rows    = matsize;
    int cols    = 1;
    if (IsSingleRowND(x))
    {
        rows = 1;
        cols = matsize;
    }

    hwMatrix* indices = EvaluatorInterface::allocateMatrix(rows, cols, hwMatrix::REAL);

    size_t valsize = y.size();
    for (int i = 0; i < matsize; ++i)
    {
        hwComplex ref  = x->z(i);
        double    rref = ref.Real();
        double    iref = ref.Imag();

        bool   isNegInf = false; // Handle -Inf/Inf/Nan
        bool   isPosInf = false;
        bool   isNan    = false;
            
        if (IsNegInf_T (rref) || IsNegInf_T (iref))
            isNegInf = true;
        else if (IsInf_T (rref) || IsInf_T(iref))
            isPosInf = true;
        else if (IsNaN_T(rref) || IsNaN_T(iref))
            isNan = true;

        for (size_t j = 0; j < valsize; ++j)
        {
            hwComplex val  = y[j];
            double    rval = val.Real();
            double    ival = val.Imag();

            if ((isNegInf && (IsNegInf_T (rval) || IsNegInf_T(ival))) ||
                (isPosInf && (IsInf_T    (rval) || IsInf_T(ival)))    ||
                (isNan   &&  (IsNaN_T    (rval) || IsNaN_T(ival)))    ||
                val == ref)
                (*indices)(i) = static_cast<int>(j) + 1;
        }
    }

    return Currency(indices);
}
//------------------------------------------------------------------------------
// True if given ND matrix has a single row 
//------------------------------------------------------------------------------
bool BuiltInFuncsElemMath::IsSingleRowND(const hwMatrixN* mtx) const
{
    if (!mtx) return false;

    std::vector<int> dims = mtx->Dimensions();
    return (!dims.empty() && dims[0] == 1);
}
//------------------------------------------------------------------------------
// Gets indices of occurences of matrix elements(x) in values => y = x(i)
//------------------------------------------------------------------------------
Currency BuiltInFuncsElemMath::GetMatrixNIndices(const hwMatrixN*             x,
                                                 const std::deque<hwComplex>& y,
                                                 bool                         forward)
{
    if (!x || x->Size() <= 0 || y.empty()) return EvaluatorInterface::allocateMatrix();

    int mtxsize = x->Size();
    int valsize = static_cast<int>(y.size());
    int rows    = valsize;
    int cols    = 1;
    if (IsSingleRowND(x))
    {
        rows = 1;
        cols = valsize;
    }
    hwMatrix* indices = EvaluatorInterface::allocateMatrix(rows, cols, hwMatrix::REAL);

    for (int i = 0; i < valsize; ++i)
    {
        hwComplex tofind = y[i];
        double    rref   = tofind.Real();
        double    iref   = tofind.Imag();

        bool   isNegInf = false; // Handle -Inf/Inf/Nan
        bool   isPosInf = false;
        bool   isNan    = false;
            
        if (IsNegInf_T (rref) || IsNegInf_T (iref))
            isNegInf = true;
        else if (IsInf_T (rref) || IsInf_T(iref))
            isPosInf = true;
        else if (IsNaN_T(rref) || IsNaN_T(iref))
            isNan = true;

        int index = forward ? 0 : mtxsize - 1;
        while (1)
        {
            if ((forward  && index >= mtxsize - 1) ||
                (!forward && index < 0))
                break;

            hwComplex val  = x->z(index);
            double    rval = val.Real();
            double    ival = val.Imag();

            if ((isNegInf && (IsNegInf_T (rval) || IsNegInf_T(ival))) ||
                (isPosInf && (IsInf_T    (rval) || IsInf_T(ival)))    ||
                (isNan   &&  (IsNaN_T    (rval) || IsNaN_T(ival)))    ||
                val == tofind)
            {
                (*indices)(i) = index + 1;
                break;
            }

            if (forward)
                index++;
            else
                index--;
        }
    }
    return Currency(indices);
}
//------------------------------------------------------------------------------
// Helper function for unique command for 2D matrices
//------------------------------------------------------------------------------
void BuiltInFuncsElemMath::UniqueHelperFuncMtx(EvaluatorInterface&    eval,
                                               const Currency&        x,
                                               bool                   forward,
                                               bool                   outputIdx,
                                               bool                   inputIdx,
                                               std::vector<Currency>& outputs)
{
    const hwMatrix* mtx = x.Matrix();
    if (!mtx || mtx->IsEmpty())  // Handle empty matrices
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        if (outputIdx)
            outputs.push_back(EvaluatorInterface::allocateMatrix());
        if (inputIdx)
            outputs.push_back(EvaluatorInterface::allocateMatrix());
        return;
    }
    else if (mtx->Size() == 1)
    {
        outputs.push_back(x);
        if (outputIdx)
            outputs.push_back(1);
        if (inputIdx)
            outputs.push_back(1);
        return;
    }

    BuiltInFuncsElemMath funcs;
    if (mtx->IsReal())
    {
        std::deque<double> y (funcs.UniqueHelperRealMtx(mtx));

        Currency out = BuiltInFuncsUtils::ContainerToMatrix(y, mtx->M() == 1);
        if (x.IsString())
            out.SetMask(Currency::MASK_STRING);
        outputs.push_back(out);

        if (outputIdx)
            outputs.push_back(funcs.GetMatrixIndices(mtx, y, forward));

        if (inputIdx)
            outputs.push_back(funcs.GetValueIndices(mtx, y));
        return;
    }

    std::deque<hwComplex> y (funcs.UniqueHelperComplexMtx(mtx));

    outputs.push_back(BuiltInFuncsUtils::ContainerToMatrix(y, mtx->M() == 1));
    if (outputIdx)
        outputs.push_back(funcs.GetMatrixIndices(mtx, y, forward));

    if (inputIdx)
        outputs.push_back(funcs.GetValueIndices(mtx, y));

}
//------------------------------------------------------------------------------
// Helper function for Unique with complex Matrix
//------------------------------------------------------------------------------
std::deque<hwComplex> BuiltInFuncsElemMath::UniqueHelperComplexMtx(const hwMatrix* mtx)
{
    assert (mtx);
    assert(!mtx->IsReal());

    bool hasNegInf = false;
    bool hasPosInf = false;
    bool hasNan    = false;
    int  numElem   = mtx->Size();

    std::deque<hwComplex> y;
    for (int i = 0; i < numElem; ++i)
    {
        hwComplex val = mtx->z(i);
        double    rval = val.Real();
        double    ival = val.Imag();

        if (IsNegInf_T (rval) || IsNegInf_T (ival))  // Handle -Inf/Inf/Nan
            hasNegInf = true;
        else if (IsInf_T(rval) || IsInf_T (ival))
            hasPosInf = true;
        else if (IsNaN_T(rval) || IsNaN_T(ival))
            hasNan = true;
        else
            y.push_back(val);
    }

    if (!y.empty() && y.size() > 1)
    {
        std::sort(y.begin(), y.end(), &complexLessThan);
        y.erase(std::unique(y.begin(), y.end()), y.end());
    }

    if (hasNegInf)
        y.insert(y.begin(), 
                 hwComplex(std::numeric_limits<double>::infinity()*-1, 0));

    if (hasPosInf)
        y.push_back(hwComplex(std::numeric_limits<double>::infinity(), 0));

    if (hasNan)
        y.push_back(hwComplex(std::numeric_limits<double>::quiet_NaN(), 0));

    return y;
}
//------------------------------------------------------------------------------
// Gets indices of occurences of matrix elements(x) in values => y = x(i)
//------------------------------------------------------------------------------
Currency BuiltInFuncsElemMath::GetMatrixIndices(const hwMatrix*           x, 
                                                const std::deque<double>& y, 
                                                bool                      forward)
{
    assert(x);
    if (!x || y.empty()) return EvaluatorInterface::allocateMatrix();

    int       mtxsize    = x->Size();
    int       valsize    = static_cast<int>(y.size());
    hwMatrix* indices = (x->M() == 1) ?
        EvaluatorInterface::allocateMatrix(1, valsize, hwMatrix::REAL) :
        EvaluatorInterface::allocateMatrix(valsize, 1, hwMatrix::REAL) ;

    for (int i = 0; i < valsize; ++i)
    {
        double tofind   = y[i];
        bool   isNegInf = false; // Handle -Inf/Inf/Nan
        bool   isPosInf = false;
        bool   isNan    = false;
            
        if (IsNegInf_T (tofind))
            isNegInf = true;
        else if (IsInf_T (tofind))
            isPosInf = true;
        else if (IsNaN_T(tofind))
            isNan = true;

        int index = forward ? 0 : mtxsize - 1;
        while (1)
        {
            if ((forward  && index > mtxsize - 1) ||
                (!forward && index < 0))
                break;

            double val = (*x)(index);
            if ((isNegInf && IsNegInf_T (val)) ||
                (isPosInf && IsInf_T(val))     ||
                (isNan    && IsNaN_T(val))     ||
                (AreEqual(val, tofind)))
            {
                (*indices)(i) = index + 1;
                break;
            }
            if (forward)
                index++;
            else
                index--;
        }
    }
    return Currency(indices);
}
//------------------------------------------------------------------------------
// Gets indices of occurences of matrix elements(x) in values => y = x(i)
//------------------------------------------------------------------------------
Currency BuiltInFuncsElemMath::GetMatrixIndices(const hwMatrix*              x, 
                                                const std::deque<hwComplex>& y,
                                                bool                         forward)
{
    assert(x);
    if (!x || y.empty()) return EvaluatorInterface::allocateMatrix();

    int       mtxsize = x->Size();
    int       valsize = static_cast<int>(y.size());
    hwMatrix* indices = (x->M() == 1) ?
        EvaluatorInterface::allocateMatrix(1, valsize, hwMatrix::REAL) :
        EvaluatorInterface::allocateMatrix(valsize, 1, hwMatrix::REAL) ;

    for (int i = 0; i < valsize; ++i)
    {
        hwComplex tofind = y[i];
        double    rref   = tofind.Real();
        double    iref   = tofind.Imag();

        bool   isNegInf = false; // Handle -Inf/Inf/Nan
        bool   isPosInf = false;
        bool   isNan    = false;
            
        if (IsNegInf_T (rref) || IsNegInf_T (iref))
            isNegInf = true;
        else if (IsInf_T (rref) || IsInf_T(iref))
            isPosInf = true;
        else if (IsNaN_T(rref) || IsNaN_T(iref))
            isNan = true;

        int index = forward ? 0 : mtxsize - 1;
        while (1)
        {
            if ((forward  && index >= mtxsize) ||
                (!forward && index < 0))
                break;

            hwComplex val  = x->z(index);
            double    rval = val.Real();
            double    ival = val.Imag();

            if ((isNegInf && (IsNegInf_T (rval) || IsNegInf_T(ival))) ||
                (isPosInf && (IsInf_T    (rval) || IsInf_T(ival)))    ||
                (isNan   &&  (IsNaN_T    (rval) || IsNaN_T(ival)))    ||
                val == tofind)
            {
                (*indices)(i) = index + 1;
                break;
            }

            if (forward)
                index++;
            else
                index--;
        }
    }
    return Currency(indices);
}
//------------------------------------------------------------------------------
// Gets indices of occurences of values in a complex matrix => x = y(j)
//------------------------------------------------------------------------------
Currency BuiltInFuncsElemMath::GetValueIndices(const hwMatrix*              x,
                                               const std::deque<hwComplex>& y)
{
    if (!x || y.empty()) return EvaluatorInterface::allocateMatrix();

    int matsize = x->Size();
    hwMatrix *indices = (x->M() == 1) ?
        EvaluatorInterface::allocateMatrix(1, matsize, hwMatrix::REAL):
        EvaluatorInterface::allocateMatrix(matsize, 1, hwMatrix::REAL);

    size_t valsize = y.size();
    for (int i = 0; i < matsize; ++i)
    {
        hwComplex ref  = x->z(i);
        double    rref = ref.Real();
        double    iref = ref.Imag();

        bool   isNegInf = false; // Handle -Inf/Inf/Nan
        bool   isPosInf = false;
        bool   isNan    = false;
            
        if (IsNegInf_T (rref) || IsNegInf_T (iref))
            isNegInf = true;
        else if (IsInf_T (rref) || IsInf_T(iref))
            isPosInf = true;
        else if (IsNaN_T(rref) || IsNaN_T(iref))
            isNan = true;

        for (size_t j = 0; j < valsize; ++j)
        {
            hwComplex val  = y[j];
            double    rval = val.Real();
            double    ival = val.Imag();

            if ((isNegInf && (IsNegInf_T (rval) || IsNegInf_T(ival))) ||
                (isPosInf && (IsInf_T    (rval) || IsInf_T(ival)))    ||
                (isNan   &&  (IsNaN_T    (rval) || IsNaN_T(ival)))    ||
                val == ref)
                (*indices)(i) = static_cast<int>(j) + 1;
        }
    }

    return Currency(indices);
}
//------------------------------------------------------------------------------
// Gets indices of occurences of values in a complex matrix => x = y(j)
//------------------------------------------------------------------------------
Currency BuiltInFuncsElemMath::GetValueIndices(const hwMatrix*           x,
                                               const std::deque<double>& y)
{
    if (!x || y.empty()) return EvaluatorInterface::allocateMatrix();

    int matsize = x->Size();
    hwMatrix *indices = (x->M() == 1) ?
        EvaluatorInterface::allocateMatrix(1, matsize, hwMatrix::REAL):
        EvaluatorInterface::allocateMatrix(matsize, 1, hwMatrix::REAL);

    size_t valsize = y.size();
    for (int i = 0; i < matsize; ++i)
    {
        double ref      = (*x)(i);
        bool   isNegInf = false; // Handle -Inf/Inf/Nan
        bool   isPosInf = false;
        bool   isNan    = false;
            
        if (IsNegInf_T (ref))
            isNegInf = true;
        else if (IsInf_T (ref))
            isPosInf = true;
        else if (IsNaN_T(ref))
            isNan = true;

        for (size_t j = 0; j < valsize; ++j)
        {
            double val = y[j];
            if ((isNegInf && IsNegInf_T (val)) ||
                (isPosInf && IsInf_T(val))     ||
                (isNan    && IsNaN_T(val))     ||
                 AreEqual(val, ref))
            {
                (*indices)(i) = static_cast<int>(j) + 1;
                break;
            }
        }
    }
    return Currency(indices);
}
//------------------------------------------------------------------------------
// Helper function for Unique, gets a deque of unique doubles
//------------------------------------------------------------------------------
std::deque<double> BuiltInFuncsElemMath::UniqueHelperRealMtx(const hwMatrix* mtx)
{
    assert (mtx);
    assert(mtx->IsReal());

    bool hasNegInf = false;
    bool hasPosInf = false;
    bool hasNan    = false;
    int  numElem   = mtx->Size();

    std::deque<double> vals;
    for (int i = 0; i < numElem; ++i)
    {
        double val = (*mtx)(i);

        if (IsNegInf_T (val))  // Handle -Inf/Inf/Nan
            hasNegInf = true;
        else if (IsInf_T(val))
            hasPosInf = true;
        else if (IsNaN_T(val))
            hasNan = true;
        else
            vals.push_back(val);
    }

    if (!vals.empty() && vals.size() > 1)
    {
        std::sort(vals.begin(), vals.end());
        vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
    }

    if (hasNegInf)
        vals.insert(vals.begin(), std::numeric_limits<double>::infinity() * -1);

    if (hasPosInf)
        vals.push_back(std::numeric_limits<double>::infinity());

    if (hasNan)
        vals.push_back(std::numeric_limits<double>::quiet_NaN());

    return vals;
}
//------------------------------------------------------------------------------
// Returns true after flipping matrix
//------------------------------------------------------------------------------
bool BuiltInFuncsElemMath::Flip(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    size_t nargin = (inputs.empty()) ? 0 : inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur1 = inputs[0];

    if (cur1.IsMatrix() || cur1.IsString())
    {
        BuiltInFuncsElemMath funcs;
        const hwMatrix* inmtx = cur1.Matrix();
        if (nargin == 1)
        {
            int dim = (inmtx->IsVector()) ? 0 : 1;
            hwMatrix* outmtx = funcs.FlipHelper(inmtx, dim);
            Currency out (outmtx);
            out.SetMask(cur1.GetMask());
            outputs.push_back(out);
            return true;
        }
        Currency cur2 = inputs[1];
        if (!cur2.IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DATA);

        int dim = static_cast<int>(cur2.Scalar());

        if (dim != 1 && dim != 2)
            throw OML_Error(OML_ERR_UNSUPPORTDIM, 2);

        hwMatrix* outmtx = funcs.FlipHelper(inmtx, dim);
        Currency out (outmtx);
        out.SetMask(cur1.GetMask());

        outputs.push_back(out);
    }
    else if (cur1.IsNDMatrix())
    {
        if (nargin == 1)
        {
            return oml_MatrixNUtil4(eval, inputs, outputs, BuiltInFuncsElemMath::Flip);
        }
        else
        {
            return oml_MatrixNUtil4(eval, inputs, outputs, BuiltInFuncsElemMath::Flip, 2);
        }
    }
    else if (cur1.IsScalar() || cur1.IsComplex())
    {
        outputs.push_back(cur1);
    }
    else
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_TYPE);
    }

    return true;
}
//------------------------------------------------------------------------------
// Helper function for flip methods which returns flipped matrix 
//------------------------------------------------------------------------------
hwMatrix* BuiltInFuncsElemMath::FlipHelper(const hwMatrix* in, int dim)
{
    assert(in);
    if (!in)
        return EvaluatorInterface::allocateMatrix();

    int       m   = in->M();
    int       n   = in->N();
    hwMatrix* out = EvaluatorInterface::allocateMatrix(m, n, in->Type());
    assert(out);

    bool isreal = in->IsReal();

    if (in->IsVector() && dim == 0)
    {
        int matrixSize = in->Size();
        for (int i = 0; i < matrixSize; ++i)
        {
            if (isreal)
                (*out)(i) = (*in)(matrixSize - 1 - i);
            else
                out->z(i) = in->z(matrixSize - 1 - i);
        }
        return out;
    }
    switch (dim)
    {
        case 2:
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < m; ++i)
                {
                    if (isreal)
                        (*out)(i, j) =  (*in)(i, n - 1 - j);
                    else
                        out->z(i, j) =  in->z(i, n - 1 - j);

                }
            }
            break;
            
        case 1:
        default:
            for (int j = 0; j < n; ++j)
            {
                for (int i = 0; i < m; ++i)
                {
                    if (isreal)
                        (*out)(i, j) = (*in)(m - 1 - i, j);
                    else
                        out->z(i, j) = in->z(m - 1 - i, j);
                }
            }
            break;
    }
    return out;
}
//------------------------------------------------------------------------------
// Returns true after flipping matrix from left to right
//------------------------------------------------------------------------------
bool BuiltInFuncsElemMath::Fliplr(EvaluatorInterface           eval,
                                  const std::vector<Currency>& inputs,
                                  std::vector<Currency>&       outputs)
{
    if (inputs.empty())
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> newinputs(inputs);
    newinputs.push_back(2);

    return Flip(eval, newinputs, outputs);
}
//------------------------------------------------------------------------------
// Returns true after flipping matrix from up to down
//------------------------------------------------------------------------------
bool BuiltInFuncsElemMath::Flipud(EvaluatorInterface           eval,
                                  const std::vector<Currency>& inputs,
                                  std::vector<Currency>&       outputs)
{
    if (inputs.empty())
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> newinputs(inputs);
    newinputs.push_back(1);

    return Flip(eval, newinputs, outputs);
}
