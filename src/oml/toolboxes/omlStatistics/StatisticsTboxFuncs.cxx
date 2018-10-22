/**
* @file StatisticsTboxFuncs.cxx
* @date January 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#include "StatisticsTboxFuncs.h"

#include <cassert>
#include <memory>  // For std::unique_ptr

#include "BuiltInFuncsUtils.h"
#include "StructData.h"
#include "MatrixNUtils.h"

#include "DistributionFuncs.h"
#include "StatisticsFuncs.h"
#include "StatUtilFuncs.h"
#include "StatisticsTests.h"

#define STATAN "StatisticalAnalysis"
#define TBOXVERSION 2019.0

static hwMersenneTwisterState* twister = nullptr;

// Helper functions

// Template function used in functions in BuiltInFuncs.cpp
template <hwMathStatus (*func)(const hwMatrix&, double&)>
double callOnVector(EvaluatorInterface& eval, const hwMatrix* vec)
{
    double val;
    hwMathStatus mstat = (*func)(*vec, val);
    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
    return val;
}

// Creates an instance of hwMersenneTwisterState
void CreateTwister();
// Helper method to get dimensions
void GetDims(EvaluatorInterface&          eval, 
             const std::vector<Currency>& currencies, 
             size_t                       offset, 
             int*                         m, 
             int*                         n);
// Helper method to get alpha and dimension options
bool ReadOptionsAlphaDim(EvaluatorInterface&          eval, 
                         const std::vector<Currency>& inputs, 
                         int                          i, 
                         double*                      alpha, 
                         int*                         dim);
// Helper method to determine output type for random number function dimension arguments
bool RNG_areDimArgsND(EvaluatorInterface&          eval,
                      const std::vector<Currency>& inputs,
                      int                          firstDimArg);
// Helper method to determine 2D output dimensions of random number function
void RNG_numRowsAndCols(const EvaluatorInterface&    eval,
                        const std::vector<Currency>& inputs,
                        int                          firstDimArg,
                        int&                         m,
                        int&                         n);
// Helper method to determine ND output dimensions of random number function
void RNG_dimensionVecAndSize(const EvaluatorInterface&    eval,
                             const std::vector<Currency>& inputs,
                             int                          firstDimArg,
                             std::vector<int>&            dims,
                             int&                         numVals);

//------------------------------------------------------------------------------
// Entry point which registers Statistics functions with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("pdf",  &OmlPdf,    FunctionMetaData(-2, -1, STATAN));
    eval.RegisterBuiltInFunction("cdf",  &OmlCdf,    FunctionMetaData(-2, -1, STATAN));
    eval.RegisterBuiltInFunction("icdf", &OmlInvcdf, FunctionMetaData(-2, -1, STATAN));
    eval.RegisterBuiltInFunction("erf",  &OmlErf,    FunctionMetaData(1, 1, STATAN));

    eval.RegisterBuiltInFunction("rand",   &OmlRand,   FunctionMetaData(-1, 1, STATAN));
    eval.RegisterBuiltInFunction("randn",  &OmlRandn,  FunctionMetaData(-1, 1, STATAN));
    eval.RegisterBuiltInFunction("random", &OmlRandom, FunctionMetaData(-2, -1, STATAN));

    eval.RegisterBuiltInFunction("betapdf", &OmlBetapdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("betacdf", &OmlBetacdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("betainv", &OmlBetainv, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("betarnd", &OmlBetarnd, FunctionMetaData(-3, 1, STATAN));
    eval.RegisterBuiltInFunction("betafit", &OmlBetafit, FunctionMetaData(1, 4, STATAN));

    eval.RegisterBuiltInFunction("chi2pdf", &OmlChi2pdf, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("chi2cdf", &OmlChi2cdf, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("chi2inv", &OmlChi2inv, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("chi2rnd", &OmlChi2rnd, FunctionMetaData(-2, 1, STATAN));

    eval.RegisterBuiltInFunction("exppdf", &OmlExppdf, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("expcdf", &OmlExpcdf, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("expinv", &OmlExpinv, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("exprnd", &OmlExprnd, FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("expfit", &OmlExpfit, FunctionMetaData(1, 2, STATAN));

    eval.RegisterBuiltInFunction("fpdf", &OmlFpdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("fcdf", &OmlFcdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("finv", &OmlFinv, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("frnd", &OmlFrnd, FunctionMetaData(-3, 1, STATAN));

    eval.RegisterBuiltInFunction("gampdf", &OmlGampdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("gamcdf", &OmlGamcdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("gaminv", &OmlGaminv, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("gamrnd", &OmlGamrnd, FunctionMetaData(-3, 1, STATAN));
    eval.RegisterBuiltInFunction("gamfit", &OmlGamfit, FunctionMetaData(1, 4, STATAN));

    eval.RegisterBuiltInFunction("lognpdf", &OmlLognpdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("logncdf", &OmlLogncdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("logninv", &OmlLogninv, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("lognrnd", &OmlLognrnd, FunctionMetaData(-3, 1, STATAN));
    eval.RegisterBuiltInFunction("lognfit", &OmlLognfit, FunctionMetaData(1, 4, STATAN));

    eval.RegisterBuiltInFunction("normpdf", &OmlNormpdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("normcdf", &OmlNormcdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("norminv", &OmlNorminv, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("normrnd", &OmlNormrnd, FunctionMetaData(-3, 1, STATAN));
    eval.RegisterBuiltInFunction("normfit", &OmlNormfit, FunctionMetaData(1,4, STATAN));
    
    eval.RegisterBuiltInFunction("poisspdf", &OmlPoisspdf, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("poisscdf", &OmlPoisscdf, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("poissinv", &OmlPoissinv, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("poissfit", &OmlPoissfit, FunctionMetaData(1, 2, STATAN));
    eval.RegisterBuiltInFunction("poissrnd", &OmlPoissrnd, FunctionMetaData(-2, 2, STATAN));

    eval.RegisterBuiltInFunction("tpdf", &OmlTpdf, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("tcdf", &OmlTcdf, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("tinv", &OmlTinv, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("trnd", &OmlTrnd, FunctionMetaData(-2, 1, STATAN));

    eval.RegisterBuiltInFunction("unifpdf", &OmlUnifpdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("unifcdf", &OmlUnifcdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("unifinv", &OmlUnifinv, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("unifrnd", &OmlUnifrnd, FunctionMetaData(-3, 1, STATAN));
    eval.RegisterBuiltInFunction("unifit",  &OmlUnifit,  FunctionMetaData(1,4, STATAN));

    eval.RegisterBuiltInFunction("wblpdf", &OmlWeibullpdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("wblcdf", &OmlWeibullcdf, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("wblinv", &OmlWeibullinv, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("wblrnd", &OmlWeibullrnd, FunctionMetaData(-3, 1, STATAN));
    eval.RegisterBuiltInFunction("wblfit", &OmlWeibullfit, FunctionMetaData(1, 4, STATAN));

    eval.RegisterBuiltInFunction("ttest",    &OmlTtest,    FunctionMetaData(-3, 3, STATAN));
    eval.RegisterBuiltInFunction("ttest2",   &OmlTtest2,   FunctionMetaData(-3, 3, STATAN));
    eval.RegisterBuiltInFunction("vartest",  &OmlChi2test, FunctionMetaData(-3, 3, STATAN));
    eval.RegisterBuiltInFunction("vartest2", &OmlFtest,    FunctionMetaData(-3, 3, STATAN));
    eval.RegisterBuiltInFunction("ztest",    &OmlZtest,    FunctionMetaData(-4, 3, STATAN));

    eval.RegisterBuiltInFunction("rms",      &OmlRms,      FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("skewness", &OmlSkewness, FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("var",      &OmlVariance, FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("std",      &OmlStd,      FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("median",   &OmlMedian,   FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("meandev",  &OmlMeandev,  FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("mean",     &OmlMean,     FunctionMetaData(3, 1, STATAN));
    eval.RegisterBuiltInFunction("cov",      &OmlCov,      FunctionMetaData(1, 1, STATAN));
    eval.RegisterBuiltInFunction("corr",     &OmlCorr,     FunctionMetaData(1, 1, STATAN));
    eval.RegisterBuiltInFunction("detrend",  &OmlDetrend,  FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("polyfit",  &OmlPolyfit,  FunctionMetaData(4, 4, STATAN));
   
    eval.RegisterBuiltInFunction("regress",  &OmlMultiregress, FunctionMetaData(-3, 5, STATAN));
    eval.RegisterBuiltInFunction("bbdesign", &OmlBBdoe,        FunctionMetaData(1, 1, STATAN));
    eval.RegisterBuiltInFunction("fullfact", &OmlFulldoe,      FunctionMetaData(1, 1, STATAN));

    return 1;
}

//------------------------------------------------------------------------------
// Computess uniform distribution probability density function values 
//------------------------------------------------------------------------------
bool OmlUnifpdf(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();

    if (nargin == 1)
    {
        newinputs.push_back(0.0);   // Lower bound
        newinputs.push_back(1.0);   // Upper bound
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
   
    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMathStatus mstat = UnifPDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = UnifPDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double    b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = UnifPDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = UnifPDF(x, realval(a, i), realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = UnifPDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = UnifPDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = UnifPDF(realval(x, i), a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = UnifPDF(realval(x, i), realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = UnifPDF(realval(x, i), realval(a, i), 
                                                 realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess uniform distribution probability density function values 
//------------------------------------------------------------------------------
bool OmlUnifcdf(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(0.0);   // Lower bound
        newinputs.push_back(1.0);   // Upper bound
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
   
    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMathStatus mstat = UnifCDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = UnifCDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double    b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = UnifCDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = UnifCDF(x, realval(a, i), realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = UnifCDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = UnifCDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = UnifCDF(realval(x, i), a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = UnifCDF(realval(x, i), realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = UnifCDF(realval(x, i), realval(a, i), 
                                                 realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess uniform distribution inverse cumulative distribution function values
//------------------------------------------------------------------------------
bool OmlUnifinv(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(0.0);   // Lower bound
        newinputs.push_back(1.0);   // Upper bound
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (newinputs.size() != 3)
        throw OML_Error(OML_ERR_NUMARGIN);
   
    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMathStatus mstat = UnifInvCDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = UnifInvCDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double    b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = UnifInvCDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = UnifInvCDF(x, realval(a, i), realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = UnifInvCDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = UnifInvCDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = UnifInvCDF(realval(x, i), a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = UnifInvCDF(realval(x, i), realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = UnifInvCDF(realval(x, i), realval(a, i), 
                                                    realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generates random data from a uniform distribution [unifrnd]
//------------------------------------------------------------------------------
bool OmlUnifrnd(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m = -1;
    int n = -1;

    GetDims(eval, inputs, 2, &m, &n);
    CreateTwister();

    // Get (or generate) the value or matrix
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    bool isScalarA = cur1.IsScalar();
    bool isScalarB = cur2.IsScalar();

    if (isScalarA && isScalarB)
    {
        double a = cur1.Scalar();
        double b = cur2.Scalar();

        if (m == 1 && n == 1)  // Only a single random number is needed
        {
            double result;
            hwMathStatus mstat = UnifRnd(a, b, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat  = UnifRnd(a, b, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
    }
    else
    {
        if (!isScalarA && (!cur1.IsMatrix() || !cur1.Matrix()->IsReal()))
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 1, OML_VAR_PARAMETER);

        hwMatrix* A = (isScalarA) ?
            EvaluatorInterface::allocateMatrix(m, n, cur1.Scalar()) :
            EvaluatorInterface::allocateMatrix(cur1.Matrix());

        if (!isScalarB && (!cur2.IsMatrix() || !cur2.Matrix()->IsReal()))
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 2, OML_VAR_PARAMETER);

        hwMatrix* B = (isScalarB) ?
            EvaluatorInterface::allocateMatrix(m, n, cur2.Scalar()) :
            EvaluatorInterface::allocateMatrix(cur2.Matrix());

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
        hwMathStatus mstat  = UnifRnd(*A, *B, twister, nullptr, *result);
        delete A;
        delete B;
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Fit a uniform distribution to a data sample [unifit]
//------------------------------------------------------------------------------
bool OmlUnifit(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur = inputs[0];

    if (!cur.IsMatrix() && !cur.IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = cur.ConvertToMatrix();
    assert(mtx);
    
    double       ahat    = 0;
    double       bhat    = 0;
    hwMatrix*    aCI     = EvaluatorInterface::allocateMatrix();
    hwMatrix*    bCI     = EvaluatorInterface::allocateMatrix();
    hwMathStatus mstatus = UnifFit(*mtx, ahat, bhat, aCI, bCI);
                           
    BuiltInFuncsUtils::CheckMathStatus(eval, mstatus);

    outputs.push_back(ahat);
    outputs.push_back(bhat);
    outputs.push_back(aCI);
    outputs.push_back(bCI);

    return true;
}
//------------------------------------------------------------------------------
// Computess normal distribution probability density function values [normpdf]
//------------------------------------------------------------------------------
bool OmlNormpdf(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(0.0);
        newinputs.push_back(1.0);
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
   
    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMathStatus mstat = NormPDF(x, mu, sigma, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma      = cur3.Matrix();
                int             sigmasize  = sigma->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         sigma->M(), sigma->N(), hwMatrix::REAL);

                for (int i = 0; i < sigmasize; ++i)
                {
                    hwMathStatus mstat = NormPDF(x, mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
            }
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu     = cur2.Matrix();
            int             musize = mu->Size();

            if (cur3.IsScalar())
            {
                double    sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = NormPDF(x, realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(mu, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = NormPDF(x, realval(mu, i), realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    sigma      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = NormPDF((*x)(i), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = NormPDF(x->z(i).Real(), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = NormPDF(realval(x, i), mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu = cur2.Matrix();

            if (!sameSize(x, mu))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = NormPDF(realval(x, i), realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = NormPDF(realval(x, i), realval(mu, i), 
                                                 realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess normal distribution cumulative distribution values [normcdf]
//------------------------------------------------------------------------------
bool OmlNormcdf(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(0.0);
        newinputs.push_back(1.0);
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
   
    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMathStatus mstat = NormCDF(x, mu, sigma, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma      = cur3.Matrix();
                int             sigmasize  = sigma->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         sigma->M(), sigma->N(), hwMatrix::REAL);

                for (int i = 0; i < sigmasize; ++i)
                {
                    hwMathStatus mstat = NormCDF(x, mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu     = cur2.Matrix();
            int             musize = mu->Size();

            if (cur3.IsScalar())
            {
                double    sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = NormCDF(x, realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(mu, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = NormCDF(x, realval(mu, i), realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    sigma      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = NormCDF((*x)(i), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = NormCDF(x->z(i).Real(), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = NormCDF(realval(x, i), mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu = cur2.Matrix();

            if (!sameSize(x, mu))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = NormCDF(realval(x, i), realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = NormCDF(realval(x, i), realval(mu, i), 
                                                 realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess normal distribution inverse cumulative distribution values [norminv]
//------------------------------------------------------------------------------
bool OmlNorminv(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,  
                std::vector<Currency>&       outputs)
{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(0.0);
        newinputs.push_back(1.0);
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
   
    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMathStatus mstat = NormInvCDF(x, mu, sigma, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma      = cur3.Matrix();
                int             sigmasize  = sigma->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         sigma->M(), sigma->N(), hwMatrix::REAL);

                for (int i = 0; i < sigmasize; ++i)
                {
                    hwMathStatus mstat = NormInvCDF(x, mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu     = cur2.Matrix();
            int             musize = mu->Size();

            if (cur3.IsScalar())
            {
                double    sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = NormInvCDF(x, realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(mu, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = NormInvCDF(x, realval(mu, i), realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    sigma      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = NormInvCDF((*x)(i), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = NormInvCDF(x->z(i).Real(), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = NormInvCDF(realval(x, i), mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu = cur2.Matrix();

            if (!sameSize(x, mu))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = NormInvCDF(realval(x, i), realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = NormInvCDF(realval(x, i), realval(mu, i), 
                                                    realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generates random data from a normal distribution [normrnd]
//------------------------------------------------------------------------------
bool OmlNormrnd(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    int firstDimArg = 3;
    bool NDout = RNG_areDimArgsND(eval, inputs, firstDimArg);

    if (nargin == 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (nargin > 2)
    {
        if (!inputs[0].IsScalar())
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, static_cast<int>(nargin), OML_VAR_DIMS);
        }

        if (!inputs[1].IsScalar())
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 2, static_cast<int>(nargin), OML_VAR_DIMS);
        }
    }

    if (nargin)
    {
        if (inputs[0].IsScalar())
        {
            if (inputs[1].IsNDMatrix())
            {
                // convert matrix to 2D
                const hwMatrixN* matrix = inputs[1].MatrixN();
                const std::vector<int>& dims = matrix->Dimensions();
                std::vector<hwSliceArg> sliceArgs;
                sliceArgs.push_back(hwSliceArg());
                hwMatrixN slice;
                matrix->SliceRHS(sliceArgs, slice);

                hwMatrix* slice2D = new hwMatrix;
                slice.ConvertNDto2D(*slice2D);

                // call computation function
                std::vector<Currency> inputs2;
                std::vector<Currency> outputs2;

                inputs2.push_back(inputs[0]);
                inputs2.push_back(slice2D);
                bool retv = OmlNormrnd(eval, inputs2, outputs2);

                // convert output to ND
                if (outputs2.size() && outputs2[0].IsMatrix())
                {
                    // convert from vector to ND
                    hwMatrixN* outMatrix = new hwMatrixN;
                    outMatrix->Convert2DtoND(*outputs2[0].Matrix());
                    outMatrix->Reshape(dims);

                    Currency out(outMatrix);
                    outputs.push_back(out);
                }

                return retv;
            }

            if (!inputs[1].IsScalar() && !inputs[1].IsMatrix())
            {
                throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
            }
        }
        else if (inputs[0].IsMatrix())
        {
            if (inputs[1].IsNDMatrix())
            {
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);
            }

            if (!inputs[1].IsScalar() && !inputs[1].IsMatrix())
            {
                throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
            }
        }
        else if (inputs[0].IsNDMatrix())
        {
            if (inputs[1].IsScalar())
            {
                // convert matrix to 2D
                const hwMatrixN* matrix = inputs[0].MatrixN();
                const std::vector<int>& dims = matrix->Dimensions();
                std::vector<hwSliceArg> sliceArgs;
                sliceArgs.push_back(hwSliceArg());
                hwMatrixN slice;
                matrix->SliceRHS(sliceArgs, slice);

                hwMatrix* slice2D = new hwMatrix;
                slice.ConvertNDto2D(*slice2D);

                // call computation function
                std::vector<Currency> inputs2;
                std::vector<Currency> outputs2;

                inputs2.push_back(slice2D);
                inputs2.push_back(inputs[1]);
                bool retv = OmlNormrnd(eval, inputs2, outputs2);

                // convert output to ND
                if (outputs2.size() && outputs2[0].IsMatrix())
                {
                    // convert from vector to ND
                    hwMatrixN* outMatrix = new hwMatrixN;
                    outMatrix->Convert2DtoND(*outputs2[0].Matrix());
                    outMatrix->Reshape(dims);

                    Currency out(outMatrix);
                    outputs.push_back(out);
                }

                return retv;
            }

            if (inputs[1].IsNDMatrix())
            {
                const hwMatrixN* matrix1 = inputs[0].MatrixN();
                const hwMatrixN* matrix2 = inputs[1].MatrixN();
                const std::vector<int>& dims = matrix1->Dimensions();

                if (matrix2->Dimensions() != dims)
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

                // convert matrices to 2D
                std::vector<hwSliceArg> sliceArgs;
                sliceArgs.push_back(hwSliceArg());
                hwMatrixN slice1;
                hwMatrixN slice2;
                matrix1->SliceRHS(sliceArgs, slice1);
                matrix2->SliceRHS(sliceArgs, slice2);

                hwMatrix* slice2D_1 = new hwMatrix;
                hwMatrix* slice2D_2 = new hwMatrix;
                slice1.ConvertNDto2D(*slice2D_1);
                slice2.ConvertNDto2D(*slice2D_2);

                // call computation function
                std::vector<Currency> inputs2;
                std::vector<Currency> outputs2;

                inputs2.push_back(slice2D_1);
                inputs2.push_back(slice2D_2);

                bool retv = OmlNormrnd(eval, inputs2, outputs2);

                // convert output to ND
                if (outputs2.size() && outputs2[0].IsMatrix())
                {
                    // convert from vector to ND
                    hwMatrixN* outMatrix = new hwMatrixN;
                    outMatrix->Convert2DtoND(*outputs2[0].Matrix());
                    outMatrix->Reshape(dims);
                    outputs.push_back(outMatrix);
                }

                return retv;
            }

            if (inputs[1].IsMatrix())
            {
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);
            }

            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 1, OML_VAR_TYPE);
        }
    }

    CreateTwister();

    // Generate the value or matrix
    if (!NDout)     // 2D case
    {
        if (nargin > 4)
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, static_cast<int>(nargin), OML_VAR_DIMS);

        if (!nargin)
        {
            double mu    = 0.0;
            double sigma = 1.0;
            double result;
            hwMathStatus mstat = NormRnd(mu, sigma, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            int m;
            int n;

            RNG_numRowsAndCols(eval, inputs, 3, m, n);

            Currency cur1 = inputs[0];
            Currency cur2 = inputs[1];
            bool isScalarM = cur1.IsScalar();
            bool isScalarS = cur2.IsScalar();

            if (isScalarM && isScalarS)
            {
                double mu    = cur1.Scalar();
                double sigma = cur2.Scalar();

                if (m == 1 && n == 1)  // Only a single random number is needed
                {
                    double result;
                    hwMathStatus mstat = NormRnd(mu, sigma, twister, nullptr, result);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    outputs.push_back(result);
                }
                else
                {
                    hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
                    hwMathStatus mstat = NormRnd(mu, sigma, twister, nullptr, *result);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    outputs.push_back(result);
                }
            }
            else
            {
                if (!isScalarM && (!cur1.IsMatrix() || !cur1.Matrix()->IsReal()))
                    throw OML_Error(OML_ERR_SCALAR_REALMTX, 1, OML_VAR_PARAMETER);

                hwMatrix* Mu = (isScalarM) ?
                    EvaluatorInterface::allocateMatrix(m, n, cur1.Scalar()) :
                EvaluatorInterface::allocateMatrix(cur1.Matrix());

                if (!isScalarS && (!cur2.IsMatrix() || !cur2.Matrix()->IsReal()))
                    throw OML_Error(OML_ERR_SCALAR_REALMTX, 2, OML_VAR_PARAMETER);

                hwMatrix* Sigma = (isScalarS) ?
                    EvaluatorInterface::allocateMatrix(m, n, cur2.Scalar()) :
                EvaluatorInterface::allocateMatrix(cur2.Matrix());

                hwMatrix*    result = EvaluatorInterface::allocateMatrix();
                hwMathStatus mstat  = NormRnd(*Mu, *Sigma, twister, nullptr, *result);
                delete Mu;
                delete Sigma;
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(result);
            }
        }
    }
    else    // ND case
    {
        std::vector<int> dims;
        int numVals = 1;

        RNG_dimensionVecAndSize(eval, inputs, 3, dims, numVals);

        Currency cur1 = inputs[0];
        Currency cur2 = inputs[1];
        bool isScalarM = cur1.IsScalar();
        bool isScalarS = cur2.IsScalar();

        if (isScalarM && isScalarS)
        {
            double mu    = cur1.Scalar();
            double sigma = cur2.Scalar();

            hwMatrix temp(numVals, 1, hwMatrix::REAL);
            hwMathStatus status = NormRnd(mu, sigma, twister, nullptr, temp);
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
            hwMatrixN* result = EvaluatorInterface::allocateMatrixN();
            result->Convert2DtoND(temp);
            result->Reshape(dims);
            outputs.push_back(result);
        }
        else
        {
            throw OML_Error(OML_ERR_ARRAYSIZE); // should never happen
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Fits a normal distribution to the given data sample [normfit]
//------------------------------------------------------------------------------
bool OmlNormfit(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur = inputs[0];

    if (!cur.IsMatrix() && !cur.IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = cur.ConvertToMatrix();
    assert(mtx);
    
    double       muhat    = 0;
    double       sigmahat = 0;
    hwMatrix*    muCI     = EvaluatorInterface::allocateMatrix();
    hwMatrix*    sigmaCI  = EvaluatorInterface::allocateMatrix();
    hwMathStatus mstatus  = NormFit(*mtx, muhat, sigmahat, muCI, sigmaCI);
                           
    BuiltInFuncsUtils::CheckMathStatus(eval, mstatus);

    outputs.push_back(muhat);
    outputs.push_back(sigmahat);
    outputs.push_back(muCI);
    outputs.push_back(sigmaCI);

    return true;
}
//------------------------------------------------------------------------------
// Computess beta distribution probability distribution function values [betapdf]
//------------------------------------------------------------------------------
bool OmlBetapdf(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 3) 
        throw OML_Error(OML_ERR_NUMARGIN);

    double   r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    Currency cur3 = inputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double       b     = cur3.Scalar();
                hwMathStatus mstat = BetaPDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = BetaPDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();                
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = BetaPDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = BetaPDF(x, realval(a, i), 
                                                 realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(x->M(), 
                                   x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = BetaPDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = BetaPDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(x->M(), 
                                   x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = BetaPDF(realval(x, i), a, realval(b, i), r); 
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(x->M(), 
                                   x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = BetaPDF(realval(x, i), 
                                                 realval(a, i), b, r); 

                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = BetaPDF(realval(x, i), realval(a, i),
                                         realval(b, i), r); 
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess beta distribution cumulative distribution values [betacdf]
//------------------------------------------------------------------------------
bool OmlBetacdf(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 3) 
        throw OML_Error(OML_ERR_NUMARGIN);
   
    double   r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    Currency cur3 = inputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMathStatus mstat = BetaCDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = BetaCDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double    b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = BetaCDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = BetaCDF(x, realval(a, i), realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = BetaCDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = BetaCDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = BetaCDF(realval(x, i), a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = BetaCDF(realval(x, i), realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = BetaCDF(realval(x, i), realval(a, i), 
                                                 realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess beta distribution inverse cumulative distribution values [betainv]
//------------------------------------------------------------------------------
bool OmlBetainv(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 3) 
        throw OML_Error(OML_ERR_NUMARGIN);

    double   r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    Currency cur3 = inputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double       b     = cur3.Scalar();
                hwMathStatus mstat = BetaInvCDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {   
                    hwMathStatus mstat = BetaInvCDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = BetaInvCDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = BetaInvCDF(x, realval(a, i), 
                                         realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = BetaInvCDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = BetaInvCDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = BetaInvCDF(realval(x, i), a, 
                                         realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = BetaInvCDF(realval(x, i), 
                                                    realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = BetaInvCDF(realval(x, i), realval(a, i), 
                                                    realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generates random data from a beta distribution [betarnd]
//------------------------------------------------------------------------------
bool OmlBetarnd(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m = -1;
    int n = -1;

    GetDims(eval, inputs, 2, &m, &n);
    CreateTwister();

    // Get (or generate) the value or matrix
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    bool isScalarA = cur1.IsScalar();
    bool isScalarB = cur2.IsScalar();

    if (isScalarA && isScalarB)
    {
        double a = cur1.Scalar();
        double b = cur2.Scalar();

        if (m == 1 && n == 1)  // Only a single random number is needed
        {
            double result;
            hwMathStatus mstat = BetaRnd(a, b, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat  = BetaRnd(a, b, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
    }
    else
    {
        if (!isScalarA && (!cur1.IsMatrix() || !cur1.Matrix()->IsReal()))
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 1, OML_VAR_PARAMETER);

        hwMatrix* A = (isScalarA) ?
            EvaluatorInterface::allocateMatrix(m, n, cur1.Scalar()) :
            EvaluatorInterface::allocateMatrix(cur1.Matrix());

        if (!isScalarB && (!cur2.IsMatrix() || !cur2.Matrix()->IsReal()))
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 2, OML_VAR_PARAMETER);

        hwMatrix* B = (isScalarB) ?
            EvaluatorInterface::allocateMatrix(m, n, cur2.Scalar()) :
            EvaluatorInterface::allocateMatrix(cur2.Matrix());

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
        hwMathStatus mstat  = BetaRnd(*A, *B, twister, nullptr, *result);
        delete A;
        delete B;
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Fits a beta distribution to the given data sample [betafit]
//------------------------------------------------------------------------------
bool OmlBetafit(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur = inputs[0];

    if (!cur.IsMatrix() && !cur.IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = cur.ConvertToMatrix();
    assert(mtx);
    
    double       ahat    = 0;
    double       bhat    = 0;
    hwMatrix*    aCI     = EvaluatorInterface::allocateMatrix();
    hwMatrix*    bCI     = EvaluatorInterface::allocateMatrix();
    hwMathStatus mstatus = BetaFit(*mtx, ahat, bhat, aCI, bCI);
                           
    BuiltInFuncsUtils::CheckMathStatus(eval, mstatus);

    outputs.push_back(ahat);
    outputs.push_back(bhat);
    outputs.push_back(aCI);
    outputs.push_back(bCI);

    return true;
}
//------------------------------------------------------------------------------
// Computes gamma distribution probability density values [gampdf]
//------------------------------------------------------------------------------
bool OmlGampdf(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)

{
    if (inputs.size() != 3) 
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    Currency cur3 = inputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMathStatus mstat = GammaPDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = GammaPDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = GammaPDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = GammaPDF(x, realval(a, i), 
                        realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = GammaPDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = GammaPDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = inputs[2].Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = GammaPDF(realval(x, i), a, 
                                                  realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = GammaPDF(realval(x, i), 
                                                  realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = GammaPDF(realval(x, i), 
                                         realval(a, i), realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computes gamma distribution cumulative distribution values [gamcdf]
//------------------------------------------------------------------------------
bool OmlGamcdf(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 3) 
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    Currency cur3 = inputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMathStatus mstat = GammaCDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = GammaCDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = GammaCDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = GammaCDF(x, realval(a, i), 
                        realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = GammaCDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = GammaCDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = GammaCDF(realval(x, i), a, 
                                                  realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = GammaCDF(realval(x, i), 
                                                  realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = GammaCDF(realval(x, i), 
                                                  realval(a, i), realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computes gamma distribution inverse cumulative distribution values [gaminv]
//------------------------------------------------------------------------------
bool OmlGaminv(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 3) 
        throw OML_Error(OML_ERR_NUMARGIN);

    double   r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    Currency cur3 = inputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double       b     = cur3.Scalar();
                hwMathStatus mstat = GammaInvCDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = GammaInvCDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = GammaInvCDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = GammaInvCDF(x, realval(a, i), 
                        realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = GammaInvCDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = GammaInvCDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = inputs[2].Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = GammaInvCDF(realval(x, i), a, 
                                                     realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = GammaInvCDF(realval(x, i), 
                                                     realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = GammaInvCDF(realval(x, i), 
                                         realval(a, i), realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generate random data from a gamma distribution [gamrnd]
//------------------------------------------------------------------------------
bool OmlGamrnd(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m = -1;
    int n = -1;

    GetDims(eval, inputs, 2, &m, &n);
    CreateTwister();

    // Get (or generate) the value or matrix
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    bool isScalarA = cur1.IsScalar();
    bool isScalarB = cur2.IsScalar();

    if (isScalarA && isScalarB)
    {
        double a = cur1.Scalar();
        double b = cur2.Scalar();

        if (m == 1 && n == 1)  // Only a single random number is needed
        {
            double result;
            hwMathStatus mstat = GammaRnd(a, b, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat  = GammaRnd(a, b, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
    }
    else
    {
        if (!isScalarA && (!cur1.IsMatrix() || !cur1.Matrix()->IsReal()))
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 1, OML_VAR_PARAMETER);

        hwMatrix* A = (isScalarA) ?
            EvaluatorInterface::allocateMatrix(m, n, cur1.Scalar()) :
            EvaluatorInterface::allocateMatrix(cur1.Matrix());

        if (!isScalarB && (!cur2.IsMatrix() || !cur2.Matrix()->IsReal()))
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 2, OML_VAR_PARAMETER);

        hwMatrix* B = (isScalarB) ?
            EvaluatorInterface::allocateMatrix(m, n, cur2.Scalar()) :
            EvaluatorInterface::allocateMatrix(cur2.Matrix());

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
        hwMathStatus mstat  = GammaRnd(*A, *B, twister, nullptr, *result);
        delete A;
        delete B;
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Fits a gamma distribution to the given data sample [gamfit]
//------------------------------------------------------------------------------
bool OmlGamfit(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur = inputs[0];

    if (!cur.IsMatrix() && !cur.IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = cur.ConvertToMatrix();
    assert(mtx);
    
    double       ahat    = 0;
    double       bhat    = 0;
    hwMatrix*    aCI     = EvaluatorInterface::allocateMatrix();
    hwMatrix*    bCI     = EvaluatorInterface::allocateMatrix();
    hwMathStatus mstatus = GammaFit(*mtx, ahat, bhat, aCI, bCI);
                           
    BuiltInFuncsUtils::CheckMathStatus(eval, mstatus);

    outputs.push_back(ahat);
    outputs.push_back(bhat);
    outputs.push_back(aCI);
    outputs.push_back(bCI);

    return true;
}
//------------------------------------------------------------------------------
// Computess exponential distribution probability density values[exppdf]
//------------------------------------------------------------------------------
bool OmlExppdf(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double       a     = cur2.Scalar();
            hwMathStatus mstat = ExpPDF(x, a, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a      = cur2.Matrix();
            int             asize  = a->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < asize; ++i)
            {
                hwMathStatus mstat = ExpPDF(x, realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double    a      = cur2.Scalar();
            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = ExpPDF((*x)(i), a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = ExpPDF(x->z(i).Real(), a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                hwMathStatus mstat = ExpPDF(realval(x, i), realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess exponential distribution cumulative distribution values [expcdf]
//------------------------------------------------------------------------------
bool OmlExpcdf(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double       a     = cur2.Scalar();
            hwMathStatus mstat = ExpCDF(x, a, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a      = cur2.Matrix();
            int             asize  = a->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < asize; ++i)
            {
                hwMathStatus mstat = ExpCDF(x, realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double    a      = cur2.Scalar();
            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = ExpCDF((*x)(i), a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = ExpCDF(x->z(i).Real(), a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                hwMathStatus mstat = ExpCDF(realval(x, i), realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess exponential distribution inverse cumulative distribution values [expinv]
//------------------------------------------------------------------------------
bool OmlExpinv(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double       a     = cur2.Scalar();
            hwMathStatus mstat = ExpInvCDF(x, a, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a      = cur2.Matrix();
            int             asize  = a->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < asize; ++i)
            {
                hwMathStatus mstat = ExpInvCDF(x, realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double    a      = cur2.Scalar();
            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = ExpInvCDF((*x)(i), a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = ExpInvCDF(x->z(i).Real(), a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                hwMathStatus mstat = ExpInvCDF(realval(x, i), realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generate random data from an exponential distribution [exprnd]
//------------------------------------------------------------------------------
bool OmlExprnd(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m = -1;
    int n = -1;

    GetDims(eval, inputs, 1, &m, &n);
    CreateTwister();

    // Get (or generate) the value or matrix
    Currency cur1 = inputs[0];
    bool isScalarA = cur1.IsScalar();

    if (isScalarA)
    {
        double a = cur1.Scalar();

        if (m == 1 && n == 1)  // Only a single random number is needed
        {
            double result;
            hwMathStatus mstat = ExpRnd(a, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat  = ExpRnd(a, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
    }
    else
    {
        if (!cur1.IsMatrix() || !cur1.Matrix()->IsReal())
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 1, OML_VAR_PARAMETER);

        const hwMatrix* A      = cur1.Matrix();
        hwMatrix*       result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
        hwMathStatus    mstat  = ExpRnd(*A, twister, nullptr, *result);
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Fits an exponential distribution to the given data sample [expfit]
//------------------------------------------------------------------------------
bool OmlExpfit(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur = inputs[0];

    if (!cur.IsMatrix() && !cur.IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = cur.ConvertToMatrix();
    assert(mtx);
    
    double       ahat    = 0;
    hwMatrix*    aCI     = EvaluatorInterface::allocateMatrix();
    hwMathStatus mstatus = ExpFit(*mtx, ahat, aCI);
                           
    BuiltInFuncsUtils::CheckMathStatus(eval, mstatus);

    outputs.push_back(ahat);
    outputs.push_back(aCI);

    return true;
}
//------------------------------------------------------------------------------
// Computess chi-squared distribution probability density values [chi2pdf]
//------------------------------------------------------------------------------
bool OmlChi2pdf(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());
            hwMathStatus mstat = Chi2PDF(x, n, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n      = cur2.Matrix();
            int             nsize  = n->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < nsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = Chi2PDF(x, static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = Chi2PDF((*x)(i), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = Chi2PDF(x->z(i).Real(), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n = cur2.Matrix();

            if (!sameSize(x, n))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = Chi2PDF(realval(x, i), 
                                     static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess chi-squared distribution cumulative distribution values [chi2cdf]
//------------------------------------------------------------------------------
bool OmlChi2cdf(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());
            hwMathStatus mstat = Chi2CDF(x, n, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n      = cur2.Matrix();
            int             nsize  = n->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < nsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = Chi2CDF(x, static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = Chi2CDF((*x)(i), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = Chi2CDF(x->z(i).Real(), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n = cur2.Matrix();

            if (!sameSize(x, n))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = Chi2CDF(realval(x, i), 
                                     static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess chi-squared distribution inverse cumulative distribution values [chi2inv]
//------------------------------------------------------------------------------
bool OmlChi2inv(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());
            hwMathStatus mstat = Chi2InvCDF(x, n, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n      = cur2.Matrix();
            int             nsize  = n->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < nsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = Chi2InvCDF(x, 
                                     static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = Chi2InvCDF((*x)(i), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = Chi2InvCDF(x->z(i).Real(), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n = cur2.Matrix();

            if (!sameSize(x, n))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = Chi2InvCDF(realval(x, i), 
                                     static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generates random data from a chi-squared distribution [chi2rnd]
//------------------------------------------------------------------------------
bool OmlChi2rnd(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m = -1;
    int n = -1;

    GetDims(eval, inputs, 1, &m, &n);
    CreateTwister();

    // Get (or generate) the value or matrix
    Currency cur1 = inputs[0];
    bool isIntegerN = cur1.IsInteger();

    if (isIntegerN)
    {
        int ndof = static_cast<int>(cur1.Scalar());

        if (m == 1 && n == 1)  // Only a single random number is needed
        {
            double result;
            hwMathStatus mstat = Chi2Rnd(ndof, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat  = Chi2Rnd(ndof, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
    }
    else
    {
        if (cur1.IsMatrix() && cur1.Matrix()->IsReal())
        {
            const hwMatrix* temp = cur1.Matrix();
            int size = temp->Size();
            hwMatrixI N(m, n, hwMatrixI::REAL);

            for (int i = 0; i < size; ++i)
            {
                if (isint((*temp)(i)))
                    N(i) = static_cast<int>((*temp)(i));
                else
                    throw OML_Error(OML_ERR_INTEGER_INTMTX, 1, OML_VAR_PARAMETER);
            }

            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat  = Chi2Rnd(N, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_INTEGER_INTMTX, 1, OML_VAR_PARAMETER);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computess Student t distribution probability density values [tpdf]
//------------------------------------------------------------------------------
bool OmlTpdf(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());
            hwMathStatus mstat = T_PDF(x, n, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n      = cur2.Matrix();
            int             nsize  = n->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < nsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = T_PDF(x, static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = T_PDF((*x)(i), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = T_PDF(x->z(i).Real(), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n = cur2.Matrix();

            if (!sameSize(x, n))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = T_PDF(realval(x, i), 
                                     static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess Student t distribution cumulative distribution function values [tcdf]
//------------------------------------------------------------------------------
bool OmlTcdf(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());
            hwMathStatus mstat = T_CDF(x, n, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n      = cur2.Matrix();
            int             nsize  = n->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < nsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = T_CDF(x, static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = T_CDF((*x)(i), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = T_CDF(x->z(i).Real(), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n = cur2.Matrix();

            if (!sameSize(x, n))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = T_CDF(realval(x, i), 
                                     static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess Student t distribution inverse cumulative distribution values [tinv]
//------------------------------------------------------------------------------
bool OmlTinv(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());
            hwMathStatus mstat = T_InvCDF(x, n, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n      = cur2.Matrix();
            int             nsize  = n->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < nsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = T_InvCDF(x, static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

            int n = static_cast<int>(cur2.Scalar());

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = T_InvCDF((*x)(i), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = T_InvCDF(x->z(i).Real(), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* n = cur2.Matrix();

            if (!sameSize(x, n))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               n->M(), n->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                if (!isint(realval(n, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

                hwMathStatus mstat = T_InvCDF(realval(x, i), 
                                     static_cast<int>(realval(n, i)), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generates random data from a Student t distribution [trnd]
//------------------------------------------------------------------------------
bool OmlTrnd(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)

{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m = -1;
    int n = -1;

    GetDims(eval, inputs, 1, &m, &n);
    CreateTwister();

    // Get (or generate) the value or matrix
    Currency cur1 = inputs[0];
    bool isIntegerN = cur1.IsInteger();

    if (isIntegerN)
    {
        int ndof = static_cast<int>(cur1.Scalar());

        if (m == 1 && n == 1)  // Only a single random number is needed
        {
            double result;
            hwMathStatus mstat = T_Rnd(ndof, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat  = T_Rnd(ndof, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
    }
    else
    {
        if (cur1.IsMatrix() && cur1.Matrix()->IsReal())
        {
            const hwMatrix* temp = cur1.Matrix();
            int size = temp->Size();
            hwMatrixI N(m, n, hwMatrixI::REAL);

            for (int i = 0; i < size; ++i)
            {
                if (isint((*temp)(i)))
                    N(i) = static_cast<int>((*temp)(i));
                else
                    throw OML_Error(OML_ERR_INTEGER_INTMTX, 1, OML_VAR_PARAMETER);
            }

            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat  = T_Rnd(N, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_INTEGER_INTMTX, 1, OML_VAR_PARAMETER);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes F distribution probability density function values [fpdf]
//------------------------------------------------------------------------------
bool OmlFpdf(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs,
             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 3) 
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    Currency cur3 = inputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

            int m = static_cast<int>(cur2.Scalar());

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int n = static_cast<int>(cur3.Scalar());
                hwMathStatus mstat = F_PDF(x, m, n, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n      = cur3.Matrix();
                int             nsize  = n->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         n->M(), n->N(), hwMatrix::REAL);

                for (int i = 0; i < nsize; ++i)
                {
                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_PDF(x, m, 
                        static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* m     = cur2.Matrix();
            int             msize = m->Size();

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int n = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   m->M(), m->N(), hwMatrix::REAL);

                for (int i = 0; i < msize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_PDF(x, 
                                         static_cast<int>(realval(m, i)), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n = cur3.Matrix();

                if (!sameSize(m, n))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(m->M(), m->N(), hwMatrix::REAL);

                for (int i = 0; i < msize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_PDF(x, 
                        static_cast<int>(realval(m, i)), 
                        static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

            int m = static_cast<int>(cur2.Scalar());

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int       n      = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = F_PDF((*x)(i), m, n, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = F_PDF(x->z(i).Real(), m, n, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n = inputs[2].Matrix();

                if (!sameSize(x, n))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_PDF(realval(x, i), m, 
                                         static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* m = cur2.Matrix();

            if (!sameSize(x, m))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int       n      = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_PDF(realval(x, i), 
                                         static_cast<int>(realval(m, i)), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n = cur3.Matrix();

                if (!sameSize(x, n))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_PDF(realval(x, i), 
                        static_cast<int>(realval(m, i)),
                        static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computes F distribution cumulative distribution values [fcdf]
//------------------------------------------------------------------------------
bool OmlFcdf(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 3) 
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    Currency cur3 = inputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

            int m = static_cast<int>(cur2.Scalar());

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int n = static_cast<int>(cur3.Scalar());
                hwMathStatus mstat = F_CDF(x, m, n, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n      = cur3.Matrix();
                int             nsize  = n->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         n->M(), n->N(), hwMatrix::REAL);

                for (int i = 0; i < nsize; ++i)
                {
                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_CDF(x, m, 
                        static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* m     = cur2.Matrix();
            int             msize = m->Size();

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int n = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   m->M(), m->N(), hwMatrix::REAL);

                for (int i = 0; i < msize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_CDF(x, 
                                         static_cast<int>(realval(m, i)), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n = cur3.Matrix();

                if (!sameSize(m, n))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(m->M(), m->N(), hwMatrix::REAL);

                for (int i = 0; i < msize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_CDF(x, 
                        static_cast<int>(realval(m, i)), 
                        static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

            int m = static_cast<int>(cur2.Scalar());

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int       n      = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = F_CDF((*x)(i), m, n, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = F_CDF(x->z(i).Real(), m, n, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n = inputs[2].Matrix();

                if (!sameSize(x, n))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_CDF(realval(x, i), m, 
                                         static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* m = cur2.Matrix();

            if (!sameSize(x, m))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int       n      = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_CDF(realval(x, i), 
                                         static_cast<int>(realval(m, i)), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n = cur3.Matrix();

                if (!sameSize(x, n))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_CDF(realval(x, i), 
                        static_cast<int>(realval(m, i)),
                        static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computes F distribution inverse cumulative distribution values [finv]
//------------------------------------------------------------------------------
bool OmlFinv(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 3) 
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    Currency cur3 = inputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

            int m = static_cast<int>(cur2.Scalar());

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int n = static_cast<int>(cur3.Scalar());
                hwMathStatus mstat = F_InvCDF(x, m, n, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n      = cur3.Matrix();
                int             nsize  = n->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         n->M(), n->N(), hwMatrix::REAL);

                for (int i = 0; i < nsize; ++i)
                {
                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_InvCDF(x, m, 
                        static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* m     = cur2.Matrix();
            int             msize = m->Size();

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int n = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   m->M(), m->N(), hwMatrix::REAL);

                for (int i = 0; i < msize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_InvCDF(x, 
                                         static_cast<int>(realval(m, i)), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n = cur3.Matrix();

                if (!sameSize(m, n))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(m->M(), m->N(), hwMatrix::REAL);

                for (int i = 0; i < msize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_InvCDF(x, 
                        static_cast<int>(realval(m, i)), 
                        static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            if (!cur2.IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

            int m = static_cast<int>(cur2.Scalar());

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int       n      = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = F_InvCDF((*x)(i), m, n, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = F_InvCDF(x->z(i).Real(), m, n, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n = inputs[2].Matrix();

                if (!sameSize(x, n))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_InvCDF(realval(x, i), m, 
                                         static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* m = cur2.Matrix();

            if (!sameSize(x, m))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int       n      = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_InvCDF(realval(x, i), 
                                         static_cast<int>(realval(m, i)), n, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* n = cur3.Matrix();

                if (!sameSize(x, n))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint(realval(m, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);

                    if (!isint(realval(n, i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                    hwMathStatus mstat = F_InvCDF(realval(x, i), 
                        static_cast<int>(realval(m, i)),
                        static_cast<int>(realval(n, i)), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generate random data from an F distribution [frnd]
//------------------------------------------------------------------------------
bool OmlFrnd(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m = -1;
    int n = -1;

    GetDims(eval, inputs, 2, &m, &n);
    CreateTwister();

    // Get (or generate) the value or matrix
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    bool isIntegerM = cur1.IsInteger();
    bool isIntegerN = cur2.IsInteger();

    if (isIntegerM && isIntegerN)
    {
        int mdof = static_cast<int>(cur1.Scalar());
        int ndof = static_cast<int>(cur2.Scalar());

        if (m == 1 && n == 1)  // Only a single random number is needed
        {
            double result;
            hwMathStatus mstat = F_Rnd(mdof, ndof, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat  = F_Rnd(mdof, ndof, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
    }
    else
    {
        hwMatrixI* M = nullptr;

        if (cur1.IsInteger())
        {
            int mdof = static_cast<int>(cur1.Scalar());
            M = EvaluatorInterface::allocateMatrix(m, n, mdof);
        }
        else if (cur1.IsMatrix() && cur1.Matrix()->IsReal())
        {
            const hwMatrix* temp = cur1.Matrix();
            int size = temp->Size();
            M = EvaluatorInterface::allocateMatrix(m, n, hwMatrixI::REAL);

            for (int i = 0; i < size; ++i)
            {
                if (isint((*temp)(i)))
                    (*M)(i) = static_cast<int>((*temp)(i));
                else
                    throw OML_Error(OML_ERR_INTEGER_INTMTX, 1, OML_VAR_PARAMETER);
            }
        }
        else
            throw OML_Error(OML_ERR_INTEGER_INTMTX, 1, OML_VAR_PARAMETER);


        hwMatrixI* N = nullptr;

        if (cur2.IsInteger())
        {
            int ndof = static_cast<int>(cur2.Scalar());
            N = EvaluatorInterface::allocateMatrix(m, n, ndof);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsReal())
        {
            const hwMatrix* temp = cur2.Matrix();
            int size = temp->Size();
            N = EvaluatorInterface::allocateMatrix(m, n, hwMatrixI::REAL);

            for (int i = 0; i < size; ++i)
            {
                if (isint((*temp)(i)))
                    (*N)(i) = static_cast<int>((*temp)(i));
                else
                    throw OML_Error(OML_ERR_INTEGER_INTMTX, 2, OML_VAR_PARAMETER);
            }
        }
        else
            throw OML_Error(OML_ERR_INTEGER_INTMTX, 2, OML_VAR_PARAMETER);

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
        hwMathStatus mstat  = F_Rnd(*M, *N, twister, nullptr, *result);
        delete M;
        delete N;
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes lognormal distribution probability density values [lognpdf]
//------------------------------------------------------------------------------
bool OmlLognpdf(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(0.0);
        newinputs.push_back(1.0);
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMathStatus mstat = LogNormPDF(x, mu, sigma, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma      = cur3.Matrix();
                int             sigmasize  = sigma->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         sigma->M(), sigma->N(), hwMatrix::REAL);

                for (int i = 0; i < sigmasize; ++i)
                {
                    hwMathStatus mstat = LogNormPDF(x, mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu     = cur2.Matrix();
            int             musize = mu->Size();

            if (cur3.IsScalar())
            {
                double    sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = LogNormPDF(x, realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(mu, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = LogNormPDF(x, realval(mu, i), realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    sigma      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = LogNormPDF((*x)(i), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = LogNormPDF(x->z(i).Real(), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = LogNormPDF(realval(x, i), mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu = cur2.Matrix();

            if (!sameSize(x, mu))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = LogNormPDF(realval(x, i), realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = LogNormPDF(realval(x, i), realval(mu, i), 
                                                 realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computess lognormal distribution cumulative distribution values [logncdf]
//------------------------------------------------------------------------------
bool OmlLogncdf(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)

{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(0.0);
        newinputs.push_back(1.0);
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
   
    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMathStatus mstat = LogNormCDF(x, mu, sigma, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma      = cur3.Matrix();
                int             sigmasize  = sigma->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         sigma->M(), sigma->N(), hwMatrix::REAL);

                for (int i = 0; i < sigmasize; ++i)
                {
                    hwMathStatus mstat = LogNormCDF(x, mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu     = cur2.Matrix();
            int             musize = mu->Size();

            if (cur3.IsScalar())
            {
                double    sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = LogNormCDF(x, realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(mu, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = LogNormCDF(x, realval(mu, i), realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    sigma      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = LogNormCDF((*x)(i), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = LogNormCDF(x->z(i).Real(), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = LogNormCDF(realval(x, i), mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu = cur2.Matrix();

            if (!sameSize(x, mu))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = LogNormCDF(realval(x, i), realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = LogNormCDF(realval(x, i), realval(mu, i), 
                                                 realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computes lognormal distribution inverse cummulative distribution values [logninv]
//------------------------------------------------------------------------------
bool OmlLogninv(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)

{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(0.0);
        newinputs.push_back(1.0);
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
   
    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMathStatus mstat = LogNormInvCDF(x, mu, sigma, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma      = cur3.Matrix();
                int             sigmasize  = sigma->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         sigma->M(), sigma->N(), hwMatrix::REAL);

                for (int i = 0; i < sigmasize; ++i)
                {
                    hwMathStatus mstat = LogNormInvCDF(x, mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu     = cur2.Matrix();
            int             musize = mu->Size();

            if (cur3.IsScalar())
            {
                double    sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = LogNormInvCDF(x, realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(mu, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), hwMatrix::REAL);

                for (int i = 0; i < musize; ++i)
                {
                    hwMathStatus mstat = LogNormInvCDF(x, realval(mu, i), realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double mu = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    sigma      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = LogNormInvCDF((*x)(i), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = LogNormInvCDF(x->z(i).Real(), mu, sigma, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = LogNormInvCDF(realval(x, i), mu, realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* mu = cur2.Matrix();

            if (!sameSize(x, mu))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = LogNormInvCDF(realval(x, i), realval(mu, i), sigma, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* sigma = cur3.Matrix();

                if (!sameSize(x, sigma))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = LogNormInvCDF(realval(x, i), realval(mu, i), 
                                                 realval(sigma, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generate random data from a lognormal distribution [lognrnd]
//------------------------------------------------------------------------------
bool OmlLognrnd(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m = -1;
    int n = -1;

    GetDims(eval, inputs, 2, &m, &n);
    CreateTwister();

    // Get (or generate) the value or matrix
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    bool isScalarM = cur1.IsScalar();
    bool isScalarS = cur2.IsScalar();

    if (isScalarM && isScalarS)
    {
        double mu    = cur1.Scalar();
        double sigma = cur2.Scalar();

        if (m == 1 && n == 1)  // Only a single random number is needed
        {
            double result;
            hwMathStatus mstat = LogNormRnd(mu, sigma, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat = LogNormRnd(mu, sigma, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
    }
    else
    {
        if (!isScalarM && (!cur1.IsMatrix() || !cur1.Matrix()->IsReal()))
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 1, OML_VAR_PARAMETER);

        hwMatrix* Mu = (isScalarM) ?
            EvaluatorInterface::allocateMatrix(m, n, cur1.Scalar()) :
            EvaluatorInterface::allocateMatrix(cur1.Matrix());

        if (!isScalarS && (!cur2.IsMatrix() || !cur2.Matrix()->IsReal()))
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 2, OML_VAR_PARAMETER);

        hwMatrix* Sigma = (isScalarS) ?
            EvaluatorInterface::allocateMatrix(m, n, cur2.Scalar()) :
            EvaluatorInterface::allocateMatrix(cur2.Matrix());

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
        hwMathStatus mstat  = LogNormRnd(*Mu, *Sigma, twister, nullptr, *result);
        delete Mu;
        delete Sigma;
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Fits a lognormal distribution to the given data sample [lognfit]
//------------------------------------------------------------------------------
bool OmlLognfit(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur = inputs[0];

    if (!cur.IsMatrix() && !cur.IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = cur.ConvertToMatrix();
    assert(mtx);
    
    double       muhat    = 0;
    double       sigmahat = 0;
    hwMatrix*    muCI     = EvaluatorInterface::allocateMatrix();
    hwMatrix*    sigmaCI  = EvaluatorInterface::allocateMatrix();
    hwMathStatus mstatus  = LogNormFit(*mtx, muhat, sigmahat, muCI, sigmaCI);
                           
    BuiltInFuncsUtils::CheckMathStatus(eval, mstatus);

    outputs.push_back(muhat);
    outputs.push_back(sigmahat);
    outputs.push_back(muCI);
    outputs.push_back(sigmaCI);

    return true;
}
//------------------------------------------------------------------------------
// Computes Weibull distribution probability density values [wblpdf]
//------------------------------------------------------------------------------
bool OmlWeibullpdf(EvaluatorInterface           eval,
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs)
{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(1.0);
        newinputs.push_back(1.0);
    }
    else if (nargin == 2)
    {
        newinputs.push_back(1.0);
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
   
    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMathStatus mstat = WeibullPDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = WeibullPDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double    b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = WeibullPDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = WeibullPDF(x, realval(a, i), realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = WeibullPDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = WeibullPDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = WeibullPDF(realval(x, i), a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = WeibullPDF(realval(x, i), realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = WeibullPDF(realval(x, i), realval(a, i), 
                                                 realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computes Weibull distribution cumulative distribution values [wblcdf]
//------------------------------------------------------------------------------
bool OmlWeibullcdf(EvaluatorInterface           eval,
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs)

{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(1.0);
        newinputs.push_back(1.0);
    }
    else if (nargin == 2)
    {
        newinputs.push_back(1.0);
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMathStatus mstat = WeibullCDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = WeibullCDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double    b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = WeibullCDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = WeibullCDF(x, realval(a, i), realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = WeibullCDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = WeibullCDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = WeibullCDF(realval(x, i), a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = WeibullCDF(realval(x, i), realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = WeibullCDF(realval(x, i), realval(a, i), 
                                                 realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computes Weibull distribution inverse cumulative distribution values [wblinv]
//------------------------------------------------------------------------------
bool OmlWeibullinv(EvaluatorInterface           eval,
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs)
{
    std::vector<Currency> newinputs(inputs);
    size_t nargin = inputs.size();    

    if (nargin == 1)
    {
        newinputs.push_back(1.0);
        newinputs.push_back(1.0);
    }
    else if (nargin == 2)
    {
        newinputs.push_back(1.0);
    }
    else if (newinputs.size() != 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    double   r;
    Currency cur1 = newinputs[0];
    Currency cur2 = newinputs[1];
    Currency cur3 = newinputs[2];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMathStatus mstat = WeibullInvCDF(x, a, b, r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                outputs.push_back(r);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b      = cur3.Matrix();
                int             bsize  = b->Size();
                hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                         b->M(), b->N(), hwMatrix::REAL);

                for (int i = 0; i < bsize; ++i)
                {
                    hwMathStatus mstat = WeibullInvCDF(x, a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a     = cur2.Matrix();
            int             asize = a->Size();

            if (cur3.IsScalar())
            {
                double    b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = WeibullInvCDF(x, realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(a, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), hwMatrix::REAL);

                for (int i = 0; i < asize; ++i)
                {
                    hwMathStatus mstat = WeibullInvCDF(x, realval(a, i), realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                if (x->IsReal())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = WeibullInvCDF((*x)(i), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else if (x->IsRealData())
                {
                    for (int i = 0; i < xsize; ++i)
                    {
                        hwMathStatus mstat = WeibullInvCDF(x->z(i).Real(), a, b, r);
                        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                        (*result)(i) = r;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = WeibullInvCDF(realval(x, i), a, realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = WeibullInvCDF(realval(x, i), realval(a, i), b, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else if (cur3.IsMatrix() && cur3.Matrix()->IsRealData())
            {
                const hwMatrix* b = cur3.Matrix();

                if (!sameSize(x, b))
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3, OML_VAR_DIMS);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), hwMatrix::REAL);

                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = WeibullInvCDF(realval(x, i), realval(a, i), 
                                                 realval(b, i), r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }

                outputs.push_back(result);
            }
            else
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_PARAMETER);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generate random data from a Weibull distribution [wblrnd]
//------------------------------------------------------------------------------
bool OmlWeibullrnd(EvaluatorInterface           eval,
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m = -1;
    int n = -1;

    GetDims(eval, inputs, 2, &m, &n);
    CreateTwister();

    // Get (or generate) the value or matrix
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];
    bool isScalarA = cur1.IsScalar();
    bool isScalarB = cur2.IsScalar();

    if (isScalarA && isScalarB)
    {
        double a = cur1.Scalar();
        double b = cur2.Scalar();

        if (m == 1 && n == 1)  // Only a single random number is needed
        {
            double result;
            hwMathStatus mstat = WeibullRnd(a, b, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus mstat  = WeibullRnd(a, b, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
    }
    else
    {
        if (!isScalarA && (!cur1.IsMatrix() || !cur1.Matrix()->IsReal()))
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 1, OML_VAR_PARAMETER);

        hwMatrix* A = (isScalarA) ?
            EvaluatorInterface::allocateMatrix(m, n, cur1.Scalar()) :
            EvaluatorInterface::allocateMatrix(cur1.Matrix());

        if (!isScalarB && (!cur2.IsMatrix() || !cur2.Matrix()->IsReal()))
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 2, OML_VAR_PARAMETER);

        hwMatrix* B = (isScalarB) ?
            EvaluatorInterface::allocateMatrix(m, n, cur2.Scalar()) :
            EvaluatorInterface::allocateMatrix(cur2.Matrix());

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
        hwMathStatus mstat  = WeibullRnd(*A, *B, twister, nullptr, *result);
        delete A;
        delete B;
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Fits a Weibull distribution to the given data sample [wblfit]
//------------------------------------------------------------------------------
bool OmlWeibullfit(EvaluatorInterface           eval,
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur = inputs[0];

    if (!cur.IsMatrix() && !cur.IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = cur.ConvertToMatrix();
    assert(mtx);
    
    double       ahat    = 0;
    double       bhat    = 0;
    hwMatrix*    aCI     = EvaluatorInterface::allocateMatrix();
    hwMatrix*    bCI     = EvaluatorInterface::allocateMatrix();
    hwMathStatus mstatus = WeibullFit(*mtx, ahat, bhat, aCI, bCI);
                           
    BuiltInFuncsUtils::CheckMathStatus(eval, mstatus);

    outputs.push_back(ahat);
    outputs.push_back(bhat);
    outputs.push_back(aCI);
    outputs.push_back(bCI);

    return true;
}
//------------------------------------------------------------------------------
// Computes Poisson distribution probability density values [poisspdf]
//------------------------------------------------------------------------------
bool OmlPoisspdf(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)

{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        if (!cur1.IsInteger())
            throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_VALUE);

        int x = static_cast<int>(cur1.Scalar());

        if (cur2.IsScalar())
        {
            double       a     = cur2.Scalar();
            hwMathStatus mstat = PoissonPDF(x, a, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a      = cur2.Matrix();
            int             asize  = a->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < asize; ++i)
            {
                hwMathStatus mstat = PoissonPDF(x, realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double    a      = cur2.Scalar();
            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint((*x)(i)))
                        throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DATA);

                    hwMathStatus mstat = PoissonPDF(static_cast<int>((*x)(i)), 
                        a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    if (!isint(x->z(i).Real()))
                        throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DATA);

                    hwMathStatus mstat = PoissonPDF(static_cast<int>(x->z(i).Real()), 
                        a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                if (!isint(realval(x, i)))
                    throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DATA);

                hwMathStatus mstat = PoissonPDF(static_cast<int>(realval(x, i)),
                                     realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computes Poisson distribution cumulative distribution values [poisscdf]
//------------------------------------------------------------------------------
bool OmlPoisscdf(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)

{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    double r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double       a     = cur2.Scalar();
            hwMathStatus mstat = PoissonCDF(x, a, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a      = cur2.Matrix();
            int             asize  = a->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < asize; ++i)
            {
                hwMathStatus mstat = PoissonCDF(x, realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double    a      = cur2.Scalar();
            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = PoissonCDF((*x)(i), a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = PoissonCDF(x->z(i).Real(), a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                hwMathStatus mstat = PoissonCDF(realval(x, i), realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Computes Poisson distribution inverse cumulative distribution values [poissinv]
//------------------------------------------------------------------------------
bool OmlPoissinv(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    int r;
    Currency cur1 = inputs[0];
    Currency cur2 = inputs[1];

    if (cur1.IsScalar())
    {
        double x = cur1.Scalar();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();
            hwMathStatus mstat = PoissonInvCDF(x, a, r);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(r);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a      = cur2.Matrix();
            int             asize  = a->Size();
            hwMatrix*       result = EvaluatorInterface::allocateMatrix(
                                     a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < asize; ++i)
            {
                hwMathStatus mstat = PoissonInvCDF(x, realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (cur1.IsMatrix() && cur1.Matrix()->IsRealData())
    {
        const hwMatrix* x     = cur1.Matrix();
        int             xsize = x->Size();

        if (cur2.IsScalar())
        {
            double a = cur2.Scalar();
            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = PoissonInvCDF((*x)(i), a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else if (x->IsRealData())
            {
                for (int i = 0; i < xsize; ++i)
                {
                    hwMathStatus mstat = PoissonInvCDF(x->z(i).Real(), a, r);
                    BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                    (*result)(i) = r;
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
            }

            outputs.push_back(result);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsRealData())
        {
            const hwMatrix* a = cur2.Matrix();

            if (!sameSize(x, a))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(a->M(), a->N(), hwMatrix::REAL);

            for (int i = 0; i < xsize; ++i)
            {
                hwMathStatus mstat = PoissonInvCDF(realval(x, i), realval(a, i), r);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
                (*result)(i) = r;
            }

            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_PARAMETER);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Generate random data from a Poisson distribution [poissrnd]
//------------------------------------------------------------------------------
bool OmlPoissrnd(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m = -1;
    int n = -1;

    GetDims(eval, inputs, 1, &m, &n);
    CreateTwister();

    // Get (or generate) the value or matrix
    Currency cur1 = inputs[0];
    bool isScalarA = cur1.IsScalar();

    if (isScalarA)
    {
        double a = cur1.Scalar();

        if (m == 1 && n == 1)  // Only a single random number is needed
        {
            double result;
            hwMathStatus mstat = PoissonRnd(a, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(result);
        }
        else
        {
            hwMatrixI    result(m, n, hwMatrixI::REAL);
            hwMathStatus mstat  = PoissonRnd(a, twister, nullptr, result);

            hwMatrix* dresult = EvaluatorInterface::allocateMatrix(result.M(), 
                result.N(), hwMatrix::REAL);

            int dsize = result.Size();

            for (int i = 0; i < dsize; ++i)
                (*dresult)(i) = static_cast<double>(result(i));

            BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            outputs.push_back(dresult);
        }
    }
    else
    {
        if (!cur1.IsMatrix() || !cur1.Matrix()->IsReal())
            throw OML_Error(OML_ERR_SCALAR_REALMTX, 1, OML_VAR_PARAMETER);

        const hwMatrix* A      = cur1.Matrix();
        hwMatrixI       result(m, n, hwMatrixI::REAL);
        hwMathStatus    mstat  = PoissonRnd(*A, twister, nullptr, result);
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);

        hwMatrix* dresult = EvaluatorInterface::allocateMatrix(result.M(), 
            result.N(), hwMatrix::REAL);

        int dsize = result.Size();

        for (int i = 0; i < dsize; ++i)
            (*dresult)(i) = static_cast<double>(result(i));

        outputs.push_back(dresult);
    }

    return true;
}
//------------------------------------------------------------------------------
// Fits a Poisson distribution to the given data sample [poissfit]
//------------------------------------------------------------------------------
bool OmlPoissfit(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency cur = inputs[0];

    if (!cur.IsMatrix() && !cur.IsScalar())
        throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = cur.ConvertToMatrix();
    assert(mtx);
    
    double       ahat    = 0;
    hwMatrix*    aCI     = EvaluatorInterface::allocateMatrix();
    hwMathStatus mstatus = PoissonFit(*mtx, ahat, aCI);
                           
    BuiltInFuncsUtils::CheckMathStatus(eval, mstatus);

    outputs.push_back(ahat);
    outputs.push_back(aCI);

    return true;
}
//------------------------------------------------------------------------------
// Generates random samples from a distribution [random]
//------------------------------------------------------------------------------
bool OmlRandom(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.size() < 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsString())
        throw OML_Error(OML_ERR_STRING, 1);

    std::string type = readOption(eval, inputs[0]);
    std::vector<Currency> newinputs(inputs.cbegin() + 1, inputs.cend());

    if (type == "beta")
        return OmlBetarnd(eval, newinputs, outputs);

    else if (type == "chi2" || type == "chisquare")
        return OmlChi2rnd(eval, newinputs, outputs);

    else if (type == "exp" || type == "exponential")
        return OmlExprnd(eval, newinputs, outputs);

    else if (type == "f")
        return OmlFrnd(eval, newinputs, outputs);

    else if (type == "gam" || type == "gamma")
        return OmlGamrnd(eval, newinputs, outputs);

    else if (type == "logn" || type == "lognormal")
        return OmlLognrnd(eval, newinputs, outputs);

    else if (type == "norm" || type == "normal")
        return OmlNormrnd(eval, newinputs, outputs);

    else if (type == "poiss" || type == "poisson")
        return OmlPoissrnd(eval, newinputs, outputs);

    else if (type == "t")
        return OmlTrnd(eval, newinputs, outputs);

    else if (type == "unif" || type == "uniform")
        return OmlUnifrnd(eval, newinputs, outputs);

    else if (type == "wbl" || type == "weibull")
        return OmlWeibullrnd(eval, newinputs, outputs);

    throw OML_Error(HW_ERROR_INVALIDOPTION(type));
    return false;
}
//------------------------------------------------------------------------------
// Computess inverse cumulative distribution function values [icdf]
//------------------------------------------------------------------------------
bool OmlInvcdf(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (!inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsString())
        throw OML_Error(OML_ERR_STRING, 1);

    std::string type = readOption(eval, inputs[0]);
    std::vector<Currency> newinputs(inputs.cbegin() + 1, inputs.cend());

    if (type == "beta")
        return OmlBetainv(eval, newinputs, outputs);

    else if (type == "chi2" || type == "chisquare")
        return OmlChi2inv(eval, newinputs, outputs);

    else if (type == "exp" || type == "exponential")
        return OmlExpinv(eval, newinputs, outputs);

    else if (type == "f")
        return OmlFinv(eval, newinputs, outputs);

    else if (type == "gam" || type == "gamma")
        return OmlGaminv(eval, newinputs, outputs);

    else if (type == "logn" || type == "lognormal")
        return OmlLogninv(eval, newinputs, outputs);

    else if (type == "norm" || type == "normal")
        return OmlNorminv(eval, newinputs, outputs);

    else if (type == "poiss" || type == "poisson")
        return OmlPoissinv(eval, newinputs, outputs);

    else if (type == "t")
        return OmlTinv(eval, newinputs, outputs);

    else if (type == "unif" || type == "uniform")
        return OmlUnifinv(eval, newinputs, outputs);

    else if (type == "wbl" || type == "weibull")
        return OmlWeibullinv(eval, newinputs, outputs);

    throw OML_Error(HW_ERROR_INVALIDOPTION(type));
    return false;
}
//------------------------------------------------------------------------------
// Computess cumulative distribution function values [cdf]
//------------------------------------------------------------------------------
bool OmlCdf(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs, 
            std::vector<Currency>&       outputs)
{
    if (!inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsString())
        throw OML_Error(OML_ERR_STRING, 1);

    std::string type = readOption(eval, inputs[0]);
    std::vector<Currency> newinputs(inputs.cbegin() + 1, inputs.cend());

    if (type == "beta")
        return OmlBetacdf(eval, newinputs, outputs);

    else if (type == "chi2" || type == "chisquare")
        return OmlChi2cdf(eval, newinputs, outputs);
    
    else if (type == "exp" || type == "exponential")
        return OmlExpcdf(eval, newinputs, outputs);
    
    else if (type == "f")
        return OmlFcdf(eval, newinputs, outputs);
   
    else if (type == "gam" || type == "gamma")
        return OmlGamcdf(eval, newinputs, outputs);
    
    else if (type == "logn" || type == "lognormal")
        return OmlLogncdf(eval, newinputs, outputs);
    
    else if (type == "norm" || type == "normal")
        return OmlNormcdf(eval, newinputs, outputs);
    
    else if (type == "poiss" || type == "poisson")
        return OmlPoisscdf(eval, newinputs, outputs);
    
    else if (type == "t")
        return OmlTcdf(eval, newinputs, outputs);
    
    else if (type == "unif" || type == "uniform")
        return OmlUnifcdf(eval, newinputs, outputs);
    
    else if (type == "wbl" || type == "weibull")
        return OmlWeibullcdf(eval, newinputs, outputs);

    throw OML_Error(HW_ERROR_INVALIDOPTION(type));
    return false;
}
//------------------------------------------------------------------------------
// Computess probability density function values [pdf]
//------------------------------------------------------------------------------
bool OmlPdf(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs, 
            std::vector<Currency>&       outputs)
{
    if (!inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsString())
        throw OML_Error(OML_ERR_STRING, 1);

    std::string type = readOption(eval, inputs[0]);
    std::vector<Currency> newinputs(inputs.cbegin() + 1, inputs.cend());

    if (type == "beta")
        return OmlBetapdf(eval, newinputs, outputs);

    else if (type == "chi2" || type == "chisquare")
        return OmlChi2pdf(eval, newinputs, outputs);
    
    else if (type == "exp" || type == "exponential")
        return OmlExppdf(eval, newinputs, outputs);
    
    else if (type == "f")
        return OmlFpdf(eval, newinputs, outputs);
    
    else if (type == "gam" || type == "gamma")
        return OmlGampdf(eval, newinputs, outputs);
    
    else if (type == "logn" || type == "lognormal")
        return OmlLognpdf(eval, newinputs, outputs);
    
    else if (type == "norm" || type == "normal")
        return OmlNormpdf(eval, newinputs, outputs);
    
    else if (type == "poiss" || type == "poisson")
        return OmlPoisspdf(eval, newinputs, outputs);

    else if (type == "t")
        return OmlTpdf(eval, newinputs, outputs);

    else if (type == "unif" || type == "uniform")
        return OmlUnifpdf(eval, newinputs, outputs);

    else if (type == "wbl" || type == "weibull")
        return OmlWeibullpdf(eval, newinputs, outputs);

    throw OML_Error(HW_ERROR_INVALIDOPTION(type));
    return false;
}
//------------------------------------------------------------------------------
// Generates uniform random values on the interval (0,1) [rand]
//------------------------------------------------------------------------------
bool OmlRand(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    int firstDimArg;

    // manage seed/state, or determine output object type: 2D or ND
    if (nargin == 0)
    {
        firstDimArg = -1;   // no dimension args
    }
    else if (nargin == 1)
    {
        const Currency& input1 = inputs[0];

        if (input1.IsScalar() || input1.IsMatrix())
        {
            firstDimArg = 1;
        }
        else if (input1.IsString())
        {
            hwMatrix* state = nullptr;

            if (input1.StringVal() == "state")
            {
                state = new hwMatrix(625, 1, hwMatrix::REAL);
                assert(state);
                if (!state)
                {
                    hwMathStatus status(HW_MATH_ERR_ALLOCFAILED);
                    throw OML_Error(status);
                }

                CreateTwister();

                for (int i = 0; i < 624; ++i)
                    (*state)(i) = static_cast<double>(twister->mt[i]);

                (*state)(624) = static_cast<double>(twister->mti);

                outputs.push_back(state);
                return true;
            }
            else
            {
                throw OML_Error(OML_ERR_OPTIONVAL, 1, OML_VAR_STRING);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_INTVECSTR, 1, OML_VAR_DIMS);
        }
    }
    else if (nargin == 2)
    {
        const Currency& input1 = inputs[0];
        const Currency& input2 = inputs[1];

        if (input1.IsInteger())
        {
            firstDimArg = 1;
        }
        else if (input1.IsString())
        {
            if (input1.StringVal() == "seed")
            {
                if (input2.IsInteger())
                {
                    double seed = input2.Scalar();

                    if (!islonglong(seed) || seed < 0)
                        throw OML_Error(OML_ERR_NATURALNUM, 2);

                    if (!twister)
                    {
                        twister = new hwMersenneTwisterState(static_cast<unsigned long>(seed));
                        twister->m_initialized = true;
                    }
                    else
                        twister->Initialize(static_cast<unsigned long>(seed));
                }
                else
                {
                    throw OML_Error(OML_ERR_NATURALNUM, 2);
                }
            }
            else if (input1.StringVal() == "state")
            {
                if (input2.IsMatrix())
                {
                    const hwMatrix* state = input2.Matrix();

                    if (state->Size() != 625)
                        throw OML_Error(OML_ERR_OPTIONVAL, 2);

                    if (!twister)
                    {
                        twister = new hwMersenneTwisterState();
                        twister->m_initialized = true;
                    }

                    for (int i = 0; i < 624; ++i)
                        twister->mt[i] = (unsigned long) (*state)(i);

                    twister->mti = static_cast<int>((*state)(624));
                }
                else
                {
                    throw OML_Error(OML_ERR_MATRIX, 1);
                }
            }
            else
            {
                throw OML_Error(OML_ERR_OPTIONVAL, 1, OML_VAR_STRING);
            }

            return true;
        }
        else
        {
            throw OML_Error(OML_ERR_INTVECSTR, 2, OML_VAR_PARAMETER);
        }
    }
    else
    {
        firstDimArg = 1;
    }
    
    bool NDout = RNG_areDimArgsND(eval, inputs, firstDimArg);

    CreateTwister();

    if (!NDout)     // 2D case
    {
        int m;
        int n;
    
        RNG_numRowsAndCols(eval, inputs, firstDimArg, m, n);

        if (m == 1 && n == 1)
        {
            double result;
            hwMathStatus status = UnifRnd(0.0, 1.0, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus status = UnifRnd(0.0, 1.0, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
            outputs.push_back(result);
        }
    }
    else    // ND case
    {
        std::vector<int> dims;
        int numVals = 1;

        RNG_dimensionVecAndSize(eval, inputs, firstDimArg, dims, numVals);

        hwMatrix temp(numVals, 1, hwMatrix::REAL);
        hwMathStatus status = UnifRnd(0.0, 1.0, twister, nullptr, temp);
        BuiltInFuncsUtils::CheckMathStatus(eval, status);
        hwMatrixN* result = EvaluatorInterface::allocateMatrixN();
        result->Convert2DtoND(temp);
        result->Reshape(dims);
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Generates standard normal random values with mean 0 and variance 1 [randn]
//------------------------------------------------------------------------------
bool OmlRandn(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    int firstDimArg;

    // manage seed/state, or determine output object type: 2D or ND
    if (nargin == 0)
    {
        firstDimArg = -1;   // no dimension args
    }
    else if (nargin == 1)
    {
        const Currency& input1 = inputs[0];

        if (input1.IsScalar() || input1.IsMatrix())
        {
            firstDimArg = 1;
        }
        else if (input1.IsString())
        {
            hwMatrix* state = nullptr;

            if (input1.StringVal() == "state")
            {
                state = new hwMatrix(625, 1, hwMatrix::REAL);
                assert(state);
                if (!state)
                {
                    hwMathStatus status(HW_MATH_ERR_ALLOCFAILED);
                    throw OML_Error(status);
                }

                CreateTwister();

                for (int i = 0; i < 624; ++i)
                    (*state)(i) = static_cast<double>(twister->mt[i]);

                (*state)(624) = static_cast<double>(twister->mti);

                outputs.push_back(state);
                return true;
            }
            else
            {
                throw OML_Error(OML_ERR_OPTIONVAL, 1, OML_VAR_STRING);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_INTVECSTR, 1, OML_VAR_DIMS);
        }
    }
    else if (nargin == 2)
    {
        const Currency& input1 = inputs[0];
        const Currency& input2 = inputs[1];

        if (input1.IsInteger())
        {
            firstDimArg = 1;
        }
        else if (input1.IsString())
        {
            if (input1.StringVal() == "seed")
            {
                if (input2.IsInteger())
                {
                    double seed = input2.Scalar();

                    if (!islonglong(seed) || seed < 0)
                        throw OML_Error(OML_ERR_NATURALNUM, 2);

                    if (!twister)
                    {
                        twister = new hwMersenneTwisterState(static_cast<unsigned long>(seed));
                        twister->m_initialized = true;
                    }
                    else
                        twister->Initialize(static_cast<unsigned long>(seed));
                }
                else
                {
                    throw OML_Error(OML_ERR_NATURALNUM, 2);
                }
            }
            else if (input1.StringVal() == "state")
            {
                if (input2.IsMatrix())
                {
                    const hwMatrix* state = input2.Matrix();

                    if (state->Size() != 625)
                        throw OML_Error(OML_ERR_OPTIONVAL, 2);

                    if (!twister)
                    {
                        twister = new hwMersenneTwisterState();
                        twister->m_initialized = true;
                    }

                    for (int i = 0; i < 624; ++i)
                        twister->mt[i] = (unsigned long) (*state)(i);

                    twister->mti = static_cast<int>((*state)(624));
                }
                else
                {
                    throw OML_Error(OML_ERR_MATRIX, 1);
                }
            }
            else
            {
                throw OML_Error(OML_ERR_OPTIONVAL, 1, OML_VAR_STRING);
            }

            return true;
        }
        else
        {
            throw OML_Error(OML_ERR_INTVECSTR, 2, OML_VAR_PARAMETER);
        }
    }
    else
    {
        firstDimArg = 1;
    }
    
    bool NDout = RNG_areDimArgsND(eval, inputs, firstDimArg);

    CreateTwister();

    if (!NDout)     // 2D case
    {
        int m;
        int n;
    
        RNG_numRowsAndCols(eval, inputs, firstDimArg, m, n);

        if (m == 1 && n == 1)
        {
            double result;
            hwMathStatus status = NormRnd(0.0, 1.0, twister, nullptr, result);
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
            outputs.push_back(result);
        }
        else
        {
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            hwMathStatus status = NormRnd(0.0, 1.0, twister, nullptr, *result);
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
            outputs.push_back(result);
        }
    }
    else    // ND case
    {
        std::vector<int> dims;
        int numVals = 1;

        RNG_dimensionVecAndSize(eval, inputs, firstDimArg, dims, numVals);

        hwMatrix temp(numVals, 1, hwMatrix::REAL);
        hwMathStatus status = NormRnd(0.0, 1.0, twister, nullptr, temp);
        BuiltInFuncsUtils::CheckMathStatus(eval, status);
        hwMatrixN* result = EvaluatorInterface::allocateMatrixN();
        result->Convert2DtoND(temp);
        result->Reshape(dims);
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Gets the error function of given scalar or real matrix [erf]
//------------------------------------------------------------------------------
bool OmlErf(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs, 
            std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency input = inputs[0];

    if (input.IsScalar())
    {
        outputs.push_back(ErrorFunc(input.Scalar()));
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* m = input.Matrix();
        hwMatrix* out = EvaluatorInterface::allocateMatrix(m->M(), m->N(), hwMatrix::REAL);

        if (m->IsReal())
        {
            for (int i = 0; i < m->Size(); i++)
                (*out)(i) = ErrorFunc((*m)(i));
        }
        else
        {
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);
        }

        outputs.push_back(out);
    }
    else
        throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Hypothesis test for the mean of a sample with unknown standard deviation 
//------------------------------------------------------------------------------
bool OmlTtest(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs)
{
    int    dim    = 0;
    double alpha  = 0.05;
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 6 || (nargin % 2))
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_DATA);

    if (!inputs[1].IsScalar())
        throw OML_Error(OML_ERR_SCALAR, 2, OML_VAR_VALUE);

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    double mu = inputs[1].Scalar();
    bool dimset = ReadOptionsAlphaDim(eval, inputs, 2, &alpha, &dim);

    if (!dimset)
        dim = (data->M() == 1) ? 2 : 1;

    if (dim > 2)
        OML_Error(OML_ERR_UNSUPPORTDIM, 3);

    double prob;
    bool   reject;
    int    max = (dim == 1) ? data->N() : data->M();

    hwMatrix* result = EvaluatorInterface::allocateMatrix(1, max, hwMatrix::REAL);
    hwMatrix* p      = EvaluatorInterface::allocateMatrix(1, max, hwMatrix::REAL);
    hwMatrix* CI     = EvaluatorInterface::allocateMatrix(2, max, hwMatrix::REAL);

    hwMatrix* piece    = EvaluatorInterface::allocateMatrix();
    hwMatrix* interval = EvaluatorInterface::allocateMatrix();

    for (int i = 0; i < max; ++i)
    {
        hwMathStatus mstat;
        
        if (dim == 1)
            mstat = data->ReadColumn(i, *piece);
        else
            mstat = data->ReadRow(i, *piece);

        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);

        mstat = TTest(*piece, mu, reject, prob, *interval, alpha);
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);

        (*result)(i) = static_cast<double>(reject);
        (*p)(i) = prob;
        (*CI)(0, i) = (*interval)(0);
        (*CI)(1, i) = (*interval)(1);
    }

    if (dim == 2)
    {
        result->Transpose();
        p->Transpose();
        CI->Transpose();
    }

    Currency rescur(result);
    rescur.SetMask(Currency::MASK_LOGICAL);

    outputs.push_back(rescur);
    outputs.push_back(p);
    outputs.push_back(CI);
    
    return true;
}
//------------------------------------------------------------------------------
// Executes an f-test for non-equal variances [vartest command]
//------------------------------------------------------------------------------
bool OmlChi2test(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 6 || (nargin % 2))
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_DATA);

    if (!inputs[1].IsScalar())
        throw OML_Error(OML_ERR_SCALAR, 2, OML_VAR_VALUE);

    int    dim    = 0;
    double alpha  = 0.05;
    bool   dimset = ReadOptionsAlphaDim(eval, inputs, 2, &alpha, &dim);

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    double          mu   = inputs[1].Scalar();

    if (!dimset)
        dim = data->M() == 1 ? 2 : 1;

    if (dim > 2)
        OML_Error(OML_ERR_UNSUPPORTDIM, 3);

    double prob;
    bool   reject;

    int max = (dim == 1) ? data->N() : data->M();

    hwMatrix* result = EvaluatorInterface::allocateMatrix(1, max, hwMatrix::REAL);
    hwMatrix* p      = EvaluatorInterface::allocateMatrix(1, max, hwMatrix::REAL);
    hwMatrix* CI     = EvaluatorInterface::allocateMatrix(2, max, hwMatrix::REAL);

    hwMatrix* piece    = EvaluatorInterface::allocateMatrix();
    hwMatrix* interval = EvaluatorInterface::allocateMatrix();

    for (int i = 0; i < max; ++i)
    {
        hwMathStatus mstat;

        if (dim == 1)
            mstat = data->ReadColumn(i, *piece);
        else
            mstat = data->ReadRow(i, *piece);

        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);

        mstat = ChiSqTest(*piece, mu, reject, prob, *interval, alpha);
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);

        (*result)(i) = static_cast<double>(reject);
        (*p)(i)      = prob;
        (*CI)(0, i)  = (*interval)(0);
        (*CI)(1, i)  = (*interval)(1);
    }

    if (dim == 2)
    {
        result->Transpose();
        p->Transpose();
        CI->Transpose();
    }

    Currency rescur(result);
    rescur.SetMask(Currency::MASK_LOGICAL);

    outputs.push_back(rescur);
    outputs.push_back(p);
    outputs.push_back(CI);
    
    return true;
}
//------------------------------------------------------------------------------
// Hypothesis test for the variances of two samples [vartest2 command]
//------------------------------------------------------------------------------
bool OmlFtest(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs)
{
    int    dim    = 1;
    double alpha  = 0.05;
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 6 || (nargin % 2))
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 2, OML_VAR_DATA);

    const hwMatrix* data1 = inputs[0].ConvertToMatrix();
    const hwMatrix* data2 = inputs[1].ConvertToMatrix();
    bool dimset = ReadOptionsAlphaDim(eval, inputs, 2, &alpha, &dim);
    double prob;
    bool   reject;
    bool   firstvec  = false;
    bool   secondvec = false;

    if (dim > 2)
        throw OML_Error(OML_ERR_UNSUPPORTDIM, 3);

    if (!dimset)
    {
        firstvec = data1->IsVector();
        secondvec = data2->IsVector();
        if (!(firstvec && secondvec))
            dim = (data1->N() == data2->N()) ? 1 : 2;
    }

    int max;
    int size1;
    int size2;

    if (dim == 1)
    {
        size1 = firstvec  ? 1 : data1->N();
        size2 = secondvec ? 1 : data2->N();
    }
    else
    {
        size1 = firstvec  ? 1 : data1->M();
        size2 = secondvec ? 1 : data2->M();
    }

    if (dimset)
    {
        if (size1 == 1)
            firstvec = true;

        if (size2 == 1)
            secondvec = true;
    }

    if (firstvec)
    {
        max = size2;
    }
    else
    {
        if (secondvec || size1 == size2)
            max = size1;
        else
            throw OML_Error(HW_ERROR_INPMATSAMESIZEALONGSPECDIMEN);
    }

    hwMatrix* result = EvaluatorInterface::allocateMatrix(1, max, hwMatrix::REAL);
    hwMatrix* p      = EvaluatorInterface::allocateMatrix(1, max, hwMatrix::REAL);
    hwMatrix* CI     = EvaluatorInterface::allocateMatrix(2, max, hwMatrix::REAL);

    hwMatrix* piece1 = firstvec ?  EvaluatorInterface::allocateMatrix(data1) : 
                                   EvaluatorInterface::allocateMatrix();
    hwMatrix* piece2 = secondvec ? EvaluatorInterface::allocateMatrix(data2) : 
                                   EvaluatorInterface::allocateMatrix();
    hwMatrix* interval = EvaluatorInterface::allocateMatrix();

    for (int i = 0; i < max; ++i)
    {
        hwMathStatus mstat;

        if (dim == 1)
        {
            if (!firstvec)
            {
                mstat = data1->ReadColumn(i, *piece1);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            }
            if (!secondvec)
            {
                mstat = data2->ReadColumn(i, *piece2);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            }
        }
        else
        {
            if (!firstvec)
            {
                mstat = data1->ReadRow(i, *piece1);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            }
            if (!secondvec)
            {
                mstat = data2->ReadRow(i, *piece2);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            }
        }

        mstat = FTest(*piece1, *piece2, reject, prob, *interval, alpha);
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);

        (*result)(i) = static_cast<double>(reject);
        (*p)(i) = prob;
        (*CI)(0, i) = (*interval)(0);
        (*CI)(1, i) = (*interval)(1);
    }
    
    if (dim == 2 || (!dimset && data1->M() == 1))
    {
        result->Transpose();
        p->Transpose();
        CI->Transpose();
    }

    Currency rescur(result);
    rescur.SetMask(Currency::MASK_LOGICAL);

    outputs.push_back(rescur);
    outputs.push_back(p);
    outputs.push_back(CI);
    return true;
}
//------------------------------------------------------------------------------
// Hypothesis test for the means of two samples with unknown and equal standard
// deviations [ttest2 command]
//------------------------------------------------------------------------------
bool OmlTtest2(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    int    dim    = 1;
    double alpha  = 0.05;
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 6 || (nargin % 2))
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 2, OML_VAR_DATA);

    const hwMatrix* data1 = inputs[0].ConvertToMatrix();
    const hwMatrix* data2 = inputs[1].ConvertToMatrix();
    bool dimset = ReadOptionsAlphaDim(eval, inputs, 2, &alpha, &dim);
    double prob;
    bool   reject;
    bool   firstvec  = false;
    bool   secondvec = false;

    if (dim > 2)
        throw OML_Error(OML_ERR_UNSUPPORTDIM, 3);

    if (!dimset)
    {
        firstvec = data1->IsVector();
        secondvec = data2->IsVector();
        if (!(firstvec && secondvec))
            dim = (data1->N() == data2->N()) ? 1 : 2;
    }

    int max;
    int size1;
    int size2;

    if (dim == 1)
    {
        size1 = firstvec  ? 1 : data1->N();
        size2 = secondvec ? 1 : data2->N();
    }
    else
    {
        size1 = firstvec  ? 1 : data1->M();
        size2 = secondvec ? 1 : data2->M();
    }

    if (dimset)
    {
        if (size1 == 1)
            firstvec = true;

        if (size2 == 1)
            secondvec = true;
    }

    if (firstvec)
    {
        max = size2;
    }
    else
    {
        if (secondvec || size1 == size2)
            max = size1;
        else
            throw OML_Error(HW_ERROR_INPMATSAMESIZEALONGSPECDIMEN);
    }

    hwMatrix* result = EvaluatorInterface::allocateMatrix(1, max, hwMatrix::REAL);
    hwMatrix* p      = EvaluatorInterface::allocateMatrix(1, max, hwMatrix::REAL);
    hwMatrix* CI     = EvaluatorInterface::allocateMatrix(2, max, hwMatrix::REAL);

    hwMatrix* piece1 = firstvec  ? EvaluatorInterface::allocateMatrix(data1) : 
                                   EvaluatorInterface::allocateMatrix();
    hwMatrix* piece2 = secondvec ? EvaluatorInterface::allocateMatrix(data2) :
                                   EvaluatorInterface::allocateMatrix();
    hwMatrix* interval = EvaluatorInterface::allocateMatrix();

    for (int i = 0; i < max; ++i)
    {
        hwMathStatus mstat;
        if (dim == 1)
        {
            if (!firstvec)
            {
                mstat = data1->ReadColumn(i, *piece1);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            }
            if (!secondvec)
            {
                mstat = data2->ReadColumn(i, *piece2);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            }
        }
        else
        {
            if (!firstvec)
            {
                mstat = data1->ReadRow(i, *piece1);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            }
            if (!secondvec)
            {
                mstat = data2->ReadRow(i, *piece2);
                BuiltInFuncsUtils::CheckMathStatus(eval, mstat);
            }
        }

        mstat = TTest2(*piece1, *piece2, reject, prob, *interval, alpha);
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);

        (*result)(i) = static_cast<double>(reject);
        (*p)(i) = prob;
        (*CI)(0, i) = (*interval)(0);
        (*CI)(1, i) = (*interval)(1);
    }
    
    if (dim == 2 || (!dimset && data1->M() == 1))
    {
        result->Transpose();
        p->Transpose();
        CI->Transpose();
    }

    Currency rescur(result);
    rescur.SetMask(Currency::MASK_LOGICAL);

    outputs.push_back(rescur);
    outputs.push_back(p);
    outputs.push_back(CI);
    return true;
}
//------------------------------------------------------------------------------
// Hypothesis test for the mean of a sample with known standard deviation [ztest] 
//------------------------------------------------------------------------------
bool OmlZtest(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs)
{
    int    dim   = 0;
    double alpha = 0.05;
    size_t nargin = inputs.size();

    if (nargin < 3 || nargin > 7 || !(nargin % 2))
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_DATA);

    if (!inputs[1].IsScalar())
        throw OML_Error(OML_ERR_SCALAR, 2, OML_VAR_VALUE);

    if (!inputs[2].IsScalar())
        throw OML_Error(OML_ERR_SCALAR, 3, OML_VAR_VALUE);

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    double mu = inputs[1].Scalar();
    double sigma = inputs[2].Scalar();

    bool dimset = ReadOptionsAlphaDim(eval, inputs, 3, &alpha, &dim);

    if (!dimset)
        dim = data->M() == 1 ? 2 : 1;

    if (dim > 2)
        OML_Error(OML_ERR_UNSUPPORTDIM, 3);

    bool reject;
    double prob;

    int max = (dim == 1)? data->N() : data->M();

    hwMatrix* result = EvaluatorInterface::allocateMatrix(1, max, hwMatrix::REAL);
    hwMatrix* p      = EvaluatorInterface::allocateMatrix(1, max, hwMatrix::REAL);
    hwMatrix* CI     = EvaluatorInterface::allocateMatrix(2, max, hwMatrix::REAL);

    hwMatrix* piece    = EvaluatorInterface::allocateMatrix();
    hwMatrix* interval = EvaluatorInterface::allocateMatrix();

    for (int i = 0; i < max; ++i)
    {
        hwMathStatus mstat = (dim == 1) ? data->ReadColumn(i, *piece) :
                                          data->ReadRow(i, *piece);
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);

        mstat = ZTest(*piece, mu, sigma, reject, prob, *interval, alpha);
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);

        (*result)(i) = static_cast<double>(reject);
        (*p)(i)      = prob;
        (*CI)(0, i)  = (*interval)(0);
        (*CI)(1, i)  = (*interval)(1);
    }

    if (dim == 2)
    {
        result->Transpose();
        p->Transpose();
        CI->Transpose();
    }

    Currency rescur(result);
    rescur.SetMask(Currency::MASK_LOGICAL);

    outputs.push_back(rescur);
    outputs.push_back(p);
    outputs.push_back(CI);
    
    return true;
}
//------------------------------------------------------------------------------
// Helper method to get alpha and dimension options
//------------------------------------------------------------------------------
bool ReadOptionsAlphaDim(EvaluatorInterface&          eval, 
                         const std::vector<Currency>& inputs, 
                         int                          i, 
                         double*                      alpha, 
                         int*                         dim)
{
    bool alphaset = false;
    bool dimset   = false;
    size_t nargin = inputs.size();

    while (i < nargin)
    {
        std::string opt = readOption(eval, inputs[i++]);

        if (!inputs[i].IsScalar())
            throw OML_Error(OML_ERR_REAL, i+1, OML_VAR_VALUE);

        double val = inputs[i++].Scalar();

        if (opt == "alpha")
        {
            if (alphaset)
                throw OML_Error(HW_ERROR_ALREADYSETALPHA);
            alphaset = true;
            *alpha = val;
        }
        else if (opt == "dim")
        {
            if (dimset)
                throw OML_Error(HW_ERROR_ALREADYSETDIM);
            dimset = true;
            if (isint(val) && val > 0)
                *dim = static_cast<int>(val);
            else
                throw OML_Error(OML_ERR_POSINTEGER, i-1, OML_VAR_DIM);
        }
        else
            throw OML_Error(HW_ERROR_INVALIDOPTION(opt));
    }
    return dimset;
}
//------------------------------------------------------------------------------
// Helper method to get dimensions
//------------------------------------------------------------------------------
void GetDims(EvaluatorInterface&          eval, 
             const std::vector<Currency>& currencies, 
             size_t                       offset, 
             int*                         m, 
             int*                         n)
{
    *m = *n = -1;
    
    size_t endOffset = currencies.size();

    for (size_t i = offset; i < endOffset; ++i)
    {
        if (currencies[i].IsString())
            endOffset = i;
    }

    if (endOffset > offset)
    {
        const std::vector<Currency> dims(currencies.begin() + offset, currencies.begin() + endOffset);

        try
        {
            getDimensionsFromInput(dims, m, n);
        }
        catch (OML_Error& err)
        {
            if (err.Arg1() != -1)
                err.Arg1(err.Arg1() + static_cast<int>(offset));

            if (err.Arg2() != -1)
                err.Arg2(err.Arg2() + static_cast<int>(offset));

            throw err;
        }
    }

    for (size_t i = 0; i < offset; ++i)
    {
        if (currencies[i].IsMatrix())
        {
            const hwMatrix* mtx = currencies[i].Matrix();
            if (mtx->Size() != 1)
            {
                if (*m == -1)
                {
                    *m = mtx->M();
                    *n = mtx->N();
                }
                else if (!(mtx->M() == *m && mtx->N() == *n))
                {
                    throw OML_Error(HW_ERROR_ALLINPMATCHSPECSIZE);
                }
            }
        }
    }
    
    if (*m == -1)
        *m = *n = 1;
}
//------------------------------------------------------------------------------
// Computess root mean square values [rms]
//------------------------------------------------------------------------------
bool OmlRms(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs, 
            std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsNDMatrix())
    {
        if (nargin == 1)
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlRms);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlRms, 2);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    int dim;

    if (nargin > 1)
    {
        if (!inputs[1].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DIM);

        dim = (int) inputs[1].Scalar();
    }
    else
    {
        dim = data->M() == 1 ? 2 : 1;
    }

    outputs.push_back(oml_MatrixUtil(eval, data, dim, &callOnVector<&RMS>));
    return true;
}
//------------------------------------------------------------------------------
// Computess skewness values [skewness]
//------------------------------------------------------------------------------
bool OmlSkewness(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsNDMatrix())
    {
        if (nargin < 2)
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlSkewness);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlSkewness, 3);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    int dim;
    bool biasCorrected;

    if (nargin > 1)
    {
        if (inputs[1].IsScalar())
        {
            if (inputs[1].Scalar() == 0.0)
                biasCorrected = true;
            else if (inputs[1].Scalar() == 1.0)
                biasCorrected = false;
            else
                throw OML_Error(OML_ERR_FLAG_01, 2);
        }
        else if (inputs[1].IsEmpty())
        {
            biasCorrected = false;  // default is opposite of var and std
        }
        else
        {
            throw OML_Error(OML_ERR_FLAG_01, 2);
        }
    }
    else
    {
        biasCorrected = false;  // default is opposite of var and std
    }

    if (nargin > 2)
    {
        if (!inputs[2].IsPositiveInteger())
            throw OML_Error(OML_ERR_UNSUPPORTDIM, 3);

        dim = static_cast<int>(inputs[2].Scalar());

        if (dim !=1 && dim != 2)
            throw OML_Error(OML_ERR_UNSUPPORTDIM, 3);
    }
    else
        dim = data->M() == 1 ? 2 : 1;

    if (data->IsEmpty())
    {
        outputs.push_back(!biasCorrected ? 0.0 : std::numeric_limits<double>::quiet_NaN());
        return true;
    }

    Currency result = oml_MatrixUtil(eval, data, dim, &callOnVector<&Skewness>);

    if (biasCorrected)
    {
        int n;

        if (dim == 1)
            n = data->M();
        else if (dim == 2)
            n = data->N();
        else
            n = 1;

        if (n < 2)
            result.GetWritableMatrix()->SetElements(std::numeric_limits<double>::quiet_NaN());
        else
        {
            double coef = sqrt(static_cast<double>(n * (n - 1))) / (n - 2);
            result.GetWritableMatrix()->MultEquals(coef);
        }
    }

    outputs.push_back(result);
    return true;
}
//------------------------------------------------------------------------------
// Computess variance values [var]
//------------------------------------------------------------------------------
bool OmlVariance(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsNDMatrix())
    {
        if (nargin == 1)
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlVariance);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlVariance, 3);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    bool sampleStat;    // false = population statistic
    int dim;

    if (!data->IsReal())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    if (nargin > 1)
    {
        if (inputs[1].IsScalar())
        {
            if (inputs[1].Scalar() == 0)
                sampleStat = true;
            else if (inputs[1].Scalar() == 1)
                sampleStat = false;
            else
                throw OML_Error(OML_ERR_FLAG_01, 2);
        }
        else if (inputs[1].IsEmpty())
        {
            sampleStat = true;
        }
        else
        {
            throw OML_Error(OML_ERR_FLAG_01, 2);
        }
    }
    else
    {
        sampleStat = true;
    }

    if (nargin > 2)
    {
        if (!inputs[2].IsPositiveInteger())
            throw OML_Error(OML_ERR_UNSUPPORTDIM, 3);

        dim = (int) inputs[2].Scalar();
    }
    else
    {
        if (data->IsVector())
            dim = -1;
        else
            dim = 1;
    }

    hwMathStatus status;

    if (dim == -1)
    {
        double var;
        status = Variance(*data, var, sampleStat);
        outputs.push_back(var);
    }
    else if (dim == 1)
    {
        const double* colPtr;
        hwMatrix* var;

        if (data->M())
        {
            var = new hwMatrix(1, data->N(), hwMatrix::REAL);

            for (int i = 0; i < data->N(); ++i)
            {
                colPtr = &((*data)(0, i));
                hwMatrix col(data->M(), (void*) colPtr, hwMatrix::REAL);
                status = Variance(col, (*var)(0, i), sampleStat);
            }
        }
        else
        {
            if (data->N())
                var = new hwMatrix(1, data->N(), hwMatrix::REAL);
            else
                var = new hwMatrix(1, 1, hwMatrix::REAL);

            var->SetElements(std::numeric_limits<double>::quiet_NaN());
        }

        outputs.push_back(var);
    }
    else if (dim == 2)
    {
        hwMatrix row;
        hwMatrix* var;

        if (data->N())
        {
            var = new hwMatrix(data->M(), 1, hwMatrix::REAL);

            for (int i = 0; i < data->M(); ++i)
            {
                status = data->ReadRow(i, row);
                status = Variance(row, (*var)(i, 0), sampleStat);
            }
        }
        else
        {
            if (data->M())
                var = new hwMatrix(data->M(), 1, hwMatrix::REAL);
            else
                var = new hwMatrix(1, 1, hwMatrix::REAL);

            var->SetElements(std::numeric_limits<double>::quiet_NaN());
        }

        outputs.push_back(var);
    }
    else
    {
        hwMatrix* var = new hwMatrix(data->M(), data->N(), hwMatrix::REAL);
        var->SetElements(0.0);

        outputs.push_back(var);
    }

    return true;
}
//------------------------------------------------------------------------------
// Gets the standard deviation of the given input [std]
//------------------------------------------------------------------------------
bool OmlStd(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs, 
            std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsNDMatrix())
    {
        if (nargin < 2)
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlStd);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlStd, 3);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    bool sampleStat;    // false = population statistic
    int dim;

    if (!data->IsReal())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    if (nargin > 1)
    {
        if (inputs[1].IsScalar())
        {
            if (inputs[1].Scalar() == 0)
                sampleStat = true;
            else if (inputs[1].Scalar() == 1)
                sampleStat = false;
            else
                throw OML_Error(OML_ERR_FLAG_01, 2);
        }
        else if (inputs[1].IsEmpty())
        {
            sampleStat = true;
        }
        else
        {
            throw OML_Error(OML_ERR_FLAG_01, 2);
        }
    }
    else
    {
        sampleStat = true;
    }

    if (nargin > 2)
    {
        if (!inputs[2].IsPositiveInteger())
            throw OML_Error(OML_ERR_UNSUPPORTDIM, 3);

        dim = (int) inputs[2].Scalar();
    }
    else
    {
        if (data->IsVector())
            dim = -1;
        else
            dim = 1;
    }

    hwMathStatus status;

    if (dim == -1)
    {
        double std;
        status = StdDev(*data, std, sampleStat);
        outputs.push_back(std);
    }
    else if (dim == 1)
    {
        const double* colPtr;
        hwMatrix* std;

        if (data->M())
        {
            std = new hwMatrix(1, data->N(), hwMatrix::REAL);

            for (int i = 0; i < data->N(); ++i)
            {
                colPtr = &((*data)(0, i));
                hwMatrix col(data->M(), (void*) colPtr, hwMatrix::REAL);
                status = StdDev(col, (*std)(0, i), sampleStat);
            }
        }
        else
        {
            if (data->N())
                std = new hwMatrix(1, data->N(), hwMatrix::REAL);
            else
                std = new hwMatrix(1, 1, hwMatrix::REAL);

            std->SetElements(std::numeric_limits<double>::quiet_NaN());
        }

        outputs.push_back(std);
    }
    else if (dim == 2)
    {
        hwMatrix row;
        hwMatrix* std;

        if (data->N())
        {
            std = new hwMatrix(data->M(), 1, hwMatrix::REAL);

            for (int i = 0; i < data->M(); ++i)
            {
                status = data->ReadRow(i, row);
                status = StdDev(row, (*std)(i, 0), sampleStat);
            }
        }
        else
        {
            if (data->M())
                std = new hwMatrix(data->M(), 1, hwMatrix::REAL);
            else
                std = new hwMatrix(1, 1, hwMatrix::REAL);

            std->SetElements(std::numeric_limits<double>::quiet_NaN());
        }

        outputs.push_back(std);
    }
    else
    {
        hwMatrix* std = new hwMatrix(data->M(), data->N(), hwMatrix::REAL);
        std->SetElements(0.0);

        outputs.push_back(std);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computess median values [median]
//------------------------------------------------------------------------------
bool OmlMedian(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsNDMatrix())
    {
        if (nargin == 1)
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlMedian);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlMedian, 2);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* data = inputs[0].ConvertToMatrix();

    if (data->IsReal() && !data->IsRealData())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    int dim;

    if (nargin > 1)
    {
        if (!inputs[1].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DIM);

        dim = static_cast<int> (inputs[1].Scalar());
    }
    else
    {
        dim = data->M() == 1 ? 2 : 1;
    }

    outputs.push_back(oml_MatrixUtil(eval, data, dim, &callOnVector<&Median>));
    return true;
}
//------------------------------------------------------------------------------
// Computess mean absolute deviation values [meandev]
//------------------------------------------------------------------------------
bool OmlMeandev(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsNDMatrix())
    {
        if (nargin == 1)
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlMeandev);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlMeandev, 2);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* data = inputs[0].ConvertToMatrix();

    if (data->IsReal() && !data->IsRealData())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    int dim;

    if (nargin > 1)
    {
        if (!inputs[1].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DIM);

        dim = static_cast<int> (inputs[1].Scalar());
    }
    else
    {
        dim = data->M() == 1 ? 2 : 1;
    }

    outputs.push_back(oml_MatrixUtil(eval, data, dim, &callOnVector<&AvgDev>));
    return true;
}
//------------------------------------------------------------------------------
// Gets the mean of the given input [mean]
//------------------------------------------------------------------------------
bool OmlMean(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsNDMatrix())
    {
        if (nargin == 1)
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlMean);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlMean, 2);
    }

    enum
    {  
        Arithmetic,
        Geometric
    } mean_type = Arithmetic;

    if (nargin > 1 && inputs.back().IsString())
    {
        --nargin;
        std::string opt = inputs.back().StringVal();
        if (opt == "a")
            mean_type = Arithmetic;
        else if (opt == "g")
            mean_type = Geometric;
        else
            throw OML_Error(HW_ERROR_INVALIDOPTION(opt));
    }
    else if (nargin > 2)
        throw OML_Error(OML_ERR_STRING, static_cast<int>(nargin));

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1);

    const hwMatrix *data = inputs[0].ConvertToMatrix();
    int dim;

    if (nargin == 1)
    {
        dim = data->M() == 1 ? 2 : 1;
    }
    else // nargin > 1
    {
        if (!inputs[1].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DIM);

        dim = static_cast<int>(inputs[1].Scalar());
    }

    if (mean_type == Arithmetic)
    {
        if (inputs.back().IsString())
            oml_sum(eval, std::vector<Currency>(inputs.cbegin(), inputs.cend() - 1), outputs);
        else
            oml_sum(eval, inputs, outputs);

        if (inputs[0].IsEmpty())
        {
            outputs[0] = std::numeric_limits<double>::quiet_NaN();
        }
        else if (outputs[0].IsMatrix())
        {
            double divby = 1;

            if (dim == 1)
                divby = data->M();
            else if (dim == 2)
                divby = data->N();
            
            outputs[0].GetWritableMatrix()->DivideEquals(divby);
        }
    }
    else if (mean_type == Geometric)
    {
        const hwMatrix* data = inputs[0].ConvertToMatrix();
        outputs.push_back(oml_MatrixUtil(eval, data, dim, &callOnVector<&GeoMean>));
    }

    return true;
}
//------------------------------------------------------------------------------
// Computess covariances [cov]
//------------------------------------------------------------------------------
bool OmlCov(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs, 
            std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* m1 = inputs[0].ConvertToMatrix();
    const hwMatrix* m2 = NULL;
    bool sampleStat;

    if (nargin == 1)
    {
        sampleStat = true;
    }
    else if (nargin == 2)
    {
        if (inputs[1].IsScalar())
        {
            if (inputs[1].Scalar() == 0)
                sampleStat = true;
            else if (inputs[1].Scalar() == 1)
                sampleStat = false;
            else
                throw OML_Error(OML_ERR_FLAG_01, 2);
        }
        else if (inputs[1].IsMatrix())
        {
            sampleStat = true;
            m2 = inputs[1].ConvertToMatrix();
        }
        else
        {
            throw OML_Error(OML_ERR_MATRIX, 2); // need better choice
        }
    }
    else // nargin == 3
    {
        if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);

        m2 = inputs[1].ConvertToMatrix();

        if (!inputs[2].IsScalar())
            throw OML_Error(OML_ERR_FLAG_01, 3);

        if (inputs[2].Scalar() == 0)
            sampleStat = true;
        else if (inputs[2].Scalar() == 1)
            sampleStat = false;
        else
            throw OML_Error(OML_ERR_FLAG_01, 3);
    }

    hwMatrix* cov = EvaluatorInterface::allocateMatrix();

    if (!m2)
    {
        if (m1->IsVector())
        {
            double cov;
            BuiltInFuncsUtils::CheckMathStatus(eval, Variance(*m1, cov, sampleStat));
            outputs.push_back(cov);
        }
        else
        {
            hwMatrix* cov = EvaluatorInterface::allocateMatrix();
            BuiltInFuncsUtils::CheckMathStatus(eval, Cov(*m1, *cov, sampleStat));
            outputs.push_back(cov);
        }
    }
    else
    {
        if (m1->IsVector() && m2->IsVector())
        {
            double cov;
            BuiltInFuncsUtils::CheckMathStatus(eval, Covariance(*m1, *m2, cov, sampleStat));
            outputs.push_back(cov);
        }
        else
        {
            std::unique_ptr<hwMatrix> cov(EvaluatorInterface::allocateMatrix());
            BuiltInFuncsUtils::CheckMathStatus(eval, Cov(*m1, *m2, *cov, sampleStat));
            outputs.push_back(cov.release());
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Computess correlation coefficients [corr]
//------------------------------------------------------------------------------
bool OmlCorr(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* m1         = inputs[0].ConvertToMatrix();
    bool            sampleStat = true;
    std::unique_ptr<hwMatrix> corr(EvaluatorInterface::allocateMatrix());

    if (nargin == 1)
    {
        if (m1->IsVector())
        {
            outputs.push_back(1.0);
        }
        else
        {
            std::unique_ptr<hwMatrix> corr(EvaluatorInterface::allocateMatrix());
            BuiltInFuncsUtils::CheckMathStatus(eval, Corr(*m1, *corr));
            outputs.push_back(corr.release());
        }
    }
    else
    {
        if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);

        const hwMatrix* m2 = inputs[1].ConvertToMatrix();

        if (m1->IsVector() && m2->IsVector())
        {
            double cov;
            double var1;
            double var2;
            BuiltInFuncsUtils::CheckMathStatus(eval, Covariance(*m1, *m2, cov));
            BuiltInFuncsUtils::CheckMathStatus(eval, Variance(*m1, var1));
            BuiltInFuncsUtils::CheckMathStatus(eval, Variance(*m2, var2));
            outputs.push_back(cov / sqrt(var1 * var2));
        }
        else
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, Corr(*m1, *m2, *corr));
            outputs.push_back(corr.release());
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Performs multiple linear regression using the model y = X * beta + e [regress]
//------------------------------------------------------------------------------
bool OmlMultiregress(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs, 
                     std::vector<Currency>&       outputs)
{
    int nargin = eval.GetNarginValue();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    int nargout = eval.GetNargoutValue();
    if (nargout > 5)
        throw OML_Error(OML_ERR_NUMARGOUT);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 2, OML_VAR_DATA);

    const hwMatrix* y = inputs[0].ConvertToMatrix();
    const hwMatrix* X = inputs[1].ConvertToMatrix();
    double alpha;
    hwMathStatus status;

    if (nargin == 3)
    {
        if (!inputs[2].IsScalar())
            throw OML_Error(OML_ERR_SCALAR, 3, OML_VAR_VALUE);

        alpha = inputs[2].Scalar();
    }
    else
    {
        alpha = 0.05;
    }

    std::unique_ptr<hwMatrix> b(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> bci(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> r(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> rci(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> stats(EvaluatorInterface::allocateMatrix());

    if (nargout < 2)
        status = MultiRegress(*y, *X, *b, alpha);
    else if (nargout == 2)
        status = MultiRegress(*y, *X, *b, alpha, bci.get());
    else if (nargout == 3)
        status = MultiRegress(*y, *X, *b, alpha, bci.get(), r.get());
    else if (nargout == 4)
        status = MultiRegress(*y, *X, *b, alpha, bci.get(), r.get(), rci.get());
    else if (nargout == 5)
        status = MultiRegress(*y, *X, *b, alpha, bci.get(), r.get(), rci.get(), stats.get());

    if (!status.IsOk())
    {
        if (status.GetArg1() == 4)
            status.SetArg1(3);

        if (status.IsWarning())
            BuiltInFuncsUtils::SetWarning(eval, status.GetMessageString());
        else
            throw OML_Error(status);
    }

    outputs.push_back(b.release());

    if (nargout > 1)
        outputs.push_back(bci.release());
    if (nargout > 2)
        outputs.push_back(r.release());
    if (nargout > 3)
        outputs.push_back(rci.release());
    if (nargout == 5)
        outputs.push_back(stats.release());
    
    return true;
}
//------------------------------------------------------------------------------
// Generates a design matrix for Box-Behnken with a given number of factors [bbdesign]
//------------------------------------------------------------------------------
bool OmlBBdoe(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsPositiveInteger())
        throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_VALUE);

    int n = static_cast<int>(inputs[0].Scalar());

    hwMatrixI mtxi;
    BuiltInFuncsUtils::CheckMathStatus(eval, BBDoe(n, mtxi));
    
    hwMatrix* mtxd = EvaluatorInterface::allocateMatrix(mtxi.M(), mtxi.N(), 
                                                        hwMatrix::REAL);
    
    for (int i = 0; i < mtxi.Size(); ++i)
        (*mtxd)(i) = (double) mtxi(i);

    outputs.push_back(mtxd);
    return true;
}
//------------------------------------------------------------------------------
// Generates full factorial design matrices [fullfact]
//------------------------------------------------------------------------------
bool OmlFulldoe(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = inputs[0].ConvertToMatrix();

    if (mtx->IsReal())
    {
        if (mtx->IsEmptyOrVector())
        {
            if (isint(mtx))
            {
                hwMatrixI mtxi(mtx->M(), mtx->N(), hwMatrixI::REAL);

                for (int i = 0; i < mtxi.Size(); ++i)
                    mtxi(i) = (int) (*mtx)(i);

                hwMatrixI outi;
                BuiltInFuncsUtils::CheckMathStatus(eval, FullDoe(mtxi, outi));
                hwMatrix* out = EvaluatorInterface::allocateMatrix(outi.M(), outi.N(), hwMatrix::REAL);

                for (int i = 0; i < outi.Size(); ++i)
                    (*out)(i) = (double) outi(i);

                outputs.push_back(out);
                return true;
            }
            else
                throw OML_Error(OML_ERR_INTEGER, 1);
        }
        else
            throw OML_Error(OML_ERR_VECTOR, 1);
    }
    else
        throw OML_Error(OML_ERR_REALVECTOR, 1);
}
//------------------------------------------------------------------------------
// Removes the mean or best fit line from a data vector [detrend]
//------------------------------------------------------------------------------
bool OmlDetrend(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_DATA);

    const hwMatrix* mtx = inputs[0].ConvertToMatrix();
    std::string     opt = "linear";

    if (nargin > 1)
    {
        // opt = readOption(eval, inputs[1]);
        if (!inputs[1].IsString())
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_PARAMETER);

        opt = inputs[1].StringVal();

        if (opt != "constant" && opt != "linear")
            throw OML_Error(OML_ERR_OPTIONVAL, 2, OML_VAR_PARAMETER);
    }

    std::unique_ptr<hwMatrix> ymod(EvaluatorInterface::allocateMatrix());
    BuiltInFuncsUtils::CheckMathStatus(eval, Detrend(*mtx, opt == "linear" ? nullptr : opt.c_str(), *ymod));

    outputs.push_back(ymod.release());
    return true;
}
//------------------------------------------------------------------------------
// Fits a polynomial to a set of paired data [polyfit]
//------------------------------------------------------------------------------
bool OmlPolyfit(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin != 3)
        throw OML_Error(OML_ERR_NUMARGIN);
    
    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_VECTOR, 2, OML_VAR_DATA);

    if (!inputs[2].IsPositiveInteger())
        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_ORDER);

    const hwMatrix* x       = inputs[0].ConvertToMatrix();
    const hwMatrix* y       = inputs[1].ConvertToMatrix();
    int             order   = static_cast<int>(inputs[2].Scalar());
    size_t          nargout = eval.GetNargoutValue();

    std::unique_ptr<hwMatrix> coef(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> yest      = NULL;
    std::unique_ptr<hwMatrix> transform = NULL;

    if (nargout < 2)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, PolyCurveFit(*x, *y, order, *coef, NULL, NULL));
    }
    else if (nargout == 2)
    {
        yest.reset(EvaluatorInterface::allocateMatrix());

        BuiltInFuncsUtils::CheckMathStatus(eval, PolyCurveFit(*x, *y, order, *coef, NULL, yest.get()));
    }
    else if (nargout == 3)
    {
        double mu;
        double sigma;
        yest.reset(EvaluatorInterface::allocateMatrix());
        transform.reset(EvaluatorInterface::allocateMatrix(1, 2, hwMatrix::REAL));

        BuiltInFuncsUtils::CheckMathStatus(eval, PolyCurveFit(*x, *y, order, *coef, NULL, yest.get(), true, &mu, &sigma));

        (*transform)(0) = mu;
        (*transform)(1) = sigma;
    }
    else
    {
        throw OML_Error(OML_ERR_NUMARGOUT);
    }

    // reverse coefficients for decreasing order
    int    index;
    double temp;

    for (int i = 0; i < coef->Size()/2; i++)
    {
        index = coef->Size()-1-i;
        temp = (*coef)(i);
        (*coef)(i) = (*coef)(index);
        (*coef)(index) = temp;
    }

    // switch to row vectors
    BuiltInFuncsUtils::CheckMathStatus(eval, coef->Transpose());  

    // return outputs
    outputs.push_back(coef.release());

    if (nargout >= 2)
    {
        Currency outStruct;
        outStruct.MakeStruct();

        // push normr into structure
        double norm;
        hwMatrix residuals;
        residuals = (*y) - (*yest);
        residuals.L2Norm(norm);
        outStruct.Struct()->SetValue(0, -1, "normr", norm);

        // push yest into structure
        outStruct.Struct()->SetValue(0, -1, "yf", yest.release());

        outputs.push_back(outStruct);
    }

    if (nargout == 3)
        outputs.push_back(transform.release());

    return true;
}
//------------------------------------------------------------------------------
// Creates an instance of hwMersenneTwisterState
//------------------------------------------------------------------------------
void CreateTwister()
{
    if (!twister)
        twister = new hwMersenneTwisterState(0);
}
//------------------------------------------------------------------------------
// Determines output type for random number function dimension arguments
//------------------------------------------------------------------------------
bool RNG_areDimArgsND(EvaluatorInterface&          eval,
                      const std::vector<Currency>& inputs,
                      int                          firstDimArg)
{
    size_t nargin = inputs.size();

    if (firstDimArg < 1 || firstDimArg > nargin)
    {
        return false;
    }
    else if (nargin == firstDimArg)
    {
        if (inputs[firstDimArg-1].IsMatrix())
        {
            const hwMatrix* matrix = inputs[firstDimArg-1].Matrix();
            if (!matrix->IsReal())
                throw OML_Error(OML_ERR_NNINTVECTOR, firstDimArg, OML_VAR_DIMS);

            if (!matrix->IsVector())
                throw OML_Error(OML_ERR_NNINTVECTOR, firstDimArg, OML_VAR_DIMS);

            if (matrix->Size() > 2)
                return true;
            else
                return false;
        }
        else if (!inputs[firstDimArg-1].IsScalar())
        {
            throw OML_Error(OML_ERR_POSINTEGER);
        }
    }

    // check for all scalars
    for (int i = firstDimArg; i < nargin; ++i)
    {
        if (!inputs[i].IsScalar())
            throw OML_Error(OML_ERR_ARRAYSIZE);
    }

    if (nargin > firstDimArg + 1)
        return true;
    else
        return false;
}
//------------------------------------------------------------------------------
// Determines 2D output dimensions of random number function
//------------------------------------------------------------------------------
void RNG_numRowsAndCols(const EvaluatorInterface&    eval,
                        const std::vector<Currency>& inputs,
                        int                          firstDimArg,
                        int&                         m,
                        int&                         n)
{
    size_t nargin = inputs.size();

    if (firstDimArg < 1 || firstDimArg > nargin)
    {
        m = 1;
        n = 1;
    }
    else if (nargin == firstDimArg)
    {
        if (inputs[firstDimArg-1].IsInteger())
        {
            m = static_cast<int>(inputs[firstDimArg-1].Scalar());
            n = m;

            if (m < 0)
                throw OML_Error(OML_ERR_NATURALNUM, firstDimArg, OML_VAR_DIM);
        }
        else if (inputs[firstDimArg-1].IsMatrix())
        {
            const hwMatrix* dims = inputs[firstDimArg-1].Matrix();

            m = static_cast<int>(realval(dims, 0));
            n = static_cast<int>(realval(dims, 1));

            if (m < 0)
                throw OML_Error(OML_ERR_NATURALNUM, firstDimArg, OML_VAR_DIMS);

            if (n < 0)
                throw OML_Error(OML_ERR_NATURALNUM, firstDimArg, OML_VAR_DIMS);
        }
        else
        {
            throw OML_Error(OML_ERR_NATURALNUM, firstDimArg, OML_VAR_DIMS);
        }
    }
    else // nargin == firstDimArg+1
    {
        m = static_cast<int>(inputs[firstDimArg-1].Scalar());
        n = static_cast<int>(inputs[firstDimArg].Scalar());

        if (m < 0)
            throw OML_Error(OML_ERR_NATURALNUM, firstDimArg, OML_VAR_DIM);

        if (n < 0)
            throw OML_Error(OML_ERR_NATURALNUM, firstDimArg+1, OML_VAR_DIM);
    }
}
//------------------------------------------------------------------------------
// Determines ND output dimensions of random number function
//------------------------------------------------------------------------------
void RNG_dimensionVecAndSize(const EvaluatorInterface&    eval,
                             const std::vector<Currency>& inputs,
                             int                          firstDimArg,
                             std::vector<int>&            dims,
                             int&                         numVals)
{
    size_t nargin = inputs.size();

    if (nargin == firstDimArg)
    {
        const hwMatrix* indx = inputs[firstDimArg-1].Matrix();

        for (int i = 0; i < indx->Size(); ++i)
        {
            double dimD = (*indx)(i);

            if (dimD < 0)
                throw OML_Error(OML_ERR_NATURALNUM, firstDimArg, OML_VAR_DIMS);

            int dim = static_cast<int>(dimD);
            dims.push_back(dim);
            numVals *= dim;
        }
    }
    else
    {
        for (int i = firstDimArg-1; i < nargin; ++i)
        {
            const Currency& cur = inputs[i];

            if (!cur.IsInteger())
                throw OML_Error(OML_ERR_NATURALNUM, i+1, OML_VAR_DIM);

            int dim = static_cast<int>(cur.Scalar());

            if (dim < 0)
                throw OML_Error(OML_ERR_NATURALNUM, i+1, OML_VAR_DIM);

            dims.push_back(dim);
            numVals *= dim;
        }
    }
}

#if 0 // Commented code
template <hwMathStatus (*func)(const hwMatrix&, hwMersenneTwisterState*, unsigned long*, hwMatrix&)>
inline bool oml_rndfunc_m(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{

    size_t nargin = inputs.size();
    
    if (nargin < 1 || nargin > 7)
        throw OML_Error(OML_ERR_NUMARGIN);

    int m;
    int n;
    // to avoid checking dimensions of the matrix input
    std::vector<Currency> temp(inputs.begin() + 1, inputs.end());
    GetDims(eval, temp, 0, &m, &n);
    CreateTwister();

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_SCALARVECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = inputs[0].ConvertToMatrix();
    std::unique_ptr<hwMatrix> outmtx(EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL));

    checkMathStatus(eval, (*func)(*mtx, twister, nullptr, *outmtx));
    outputs.push_back(outmtx.release());
    return true;
}
#define rndfunc_m(func)          oml_rndfunc_m<func>
#define realfunc_dm_switch(func) oml_realfunc_dm_switch<func>
bool oml_realfunc_dm_switch(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
    hwMathStatus (*func)(const double&, const hwMatrix&, double&))
{
    if (inputs.size() != 2) 
        throw OML_Error(OML_ERR_NUMARGIN);

    const hwMatrix* mtx = inputs[0].ConvertToMatrix();

    if (inputs[1].IsScalar())
    {
        double result;
        checkMathStatus(eval, (*func)(inputs[1].Scalar(), *mtx, result));
        outputs.push_back(result);
    }
    else if (inputs[1].IsMatrix() && inputs[1].Matrix()->IsRealData())
    {
        const hwMatrix* data = inputs[1].Matrix();
        std::unique_ptr<hwMatrix> result(EvaluatorInterface::allocateMatrix(data->M(), data->N(), hwMatrix::REAL));
        double temp;
        for (int i = 0; i < data->Size(); ++i)
        {
            checkMathStatus(eval, (*func)(realval(data, i), *mtx, temp));
            (*result)(i) = temp;
        }
        outputs.push_back(result.release());
    }
    else
        throw OML_Error(OML_ERR_REAL, 2);

    return true;
}
void getRndArg(const Currency& input, hwMatrix*& arg, std::unique_ptr<hwMatrix>& holder, int m, int n)
{
    if (input.IsScalar())
    {
        arg = EvaluatorInterface::allocateMatrix(m, n, input.Scalar());
        holder.reset(arg);
    }
    else if (input.IsMatrix() && input.Matrix()->IsReal())
    {
        arg = EvaluatorInterface::allocateMatrix(input.Matrix());
    }
    else
        throw OML_Error(OML_ERR_REAL);
}
void getRndArg(const Currency& input, hwMatrixI*& arg, std::unique_ptr<hwMatrixI>& holder, int m, int n)
{
    if (input.IsScalar())
    {
        double value = input.Scalar();

        if (!IsInteger(value).IsOk())
            throw OML_Error(OML_ERR_INTEGER, -1, OML_VAR_VALUE);

        arg = EvaluatorInterface::allocateMatrix(m, n, (int) value);
    }
    else if (input.IsMatrix() && input.Matrix()->IsReal())
    {
        const hwMatrix* temp = input.Matrix();
        arg = EvaluatorInterface::allocateMatrix(m, n, hwMatrixI::REAL);
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                if (isint((*temp)(i,j)))
                    (*arg)(i,j) = (int) (*temp)(i,j);
                else
                    throw OML_Error(OML_ERR_INTEGER);
            }
        }
    }
    else
        throw OML_Error(OML_ERR_INTEGER);

    holder.reset(arg);
}

#endif
//------------------------------------------------------------------------------
// Returns toolbox version
//------------------------------------------------------------------------------
double GetToolboxVersion(EvaluatorInterface eval)
{
    return TBOXVERSION;
}
