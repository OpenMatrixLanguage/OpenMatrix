/**
* @file StatisticsTboxFuncs.cxx
* @date January 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#include "StatisticsTboxFuncs.h"

#include <cassert>
#include <memory>  // For std::unique_ptr

#include "BuiltInFuncs.h"
#include "BuiltInFuncsUtils.h"
#include "BuiltInFuncsData.h"
#include "StructData.h"
#include "MatrixNUtils.h"

#include "hwMatrix.h"
#include "hwMatrixN.h"
#include "OML_Error.h"
#include "hwMersenneTwister.h"
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
    eval.RegisterBuiltInFunction("betainc", &OmlBetacdf, FunctionMetaData(3, 1, STATAN));

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
    eval.RegisterBuiltInFunction("gammainc", &OmlGammaInc, FunctionMetaData(-3, 1, STATAN));

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

    eval.RegisterBuiltInFunction("rms",      &OmlRms,      FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("skewness", &OmlSkewness, FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("kurtosis", &OmlKurtosis, FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("var",      &OmlVariance, FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("std",      &OmlStd,      FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("median",   &OmlMedian,   FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("quantile", &OmlQuantile, FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("histc",    &OmlHistC,    FunctionMetaData(-3, -2, STATAN));
    eval.RegisterBuiltInFunction("meandev",  &OmlMeandev,  FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("mad",      &OmlMAD,      FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("mean",     &OmlMean,     FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("mode",     &OmlMode,     FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("movmean",  &OmlMovMean,  FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("cov",      &OmlCov,      FunctionMetaData(1, 1, STATAN));
    eval.RegisterBuiltInFunction("corr",     &OmlCorr,     FunctionMetaData(1, 1, STATAN));
    eval.RegisterBuiltInFunction("detrend",  &OmlDetrend,  FunctionMetaData(2, 1, STATAN));
    eval.RegisterBuiltInFunction("polyfit",  &OmlPolyfit,  FunctionMetaData(4, 4, STATAN));
    eval.RegisterBuiltInFunction("nchoosek", &OmlNchooseK, FunctionMetaData(2, 1, STATAN));

    eval.RegisterBuiltInFunction("regress",  &OmlMultiregress, FunctionMetaData(-3, 5, STATAN));
    eval.RegisterBuiltInFunction("randperm", &OmlRandperm,     FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("bbdesign", &OmlBBdoe,        FunctionMetaData(1, 1, STATAN));
    eval.RegisterBuiltInFunction("fullfact", &OmlFulldoe,      FunctionMetaData(1, 1, STATAN));

    eval.RegisterBuiltInFunction("nanmax",    &OmlNanMax,      FunctionMetaData(-2, -2, STATAN));
    eval.RegisterBuiltInFunction("nanmin",    &OmlNanMin,      FunctionMetaData(-2, -2, STATAN));
    eval.RegisterBuiltInFunction("nansum",    &OmlNanSum,      FunctionMetaData(-2, -1, STATAN));
    eval.RegisterBuiltInFunction("nanmean",   &OmlNanMean,     FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("nanstd",    &OmlNanStd,      FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("nanvar",    &OmlNanVar,      FunctionMetaData(-2, 1, STATAN));
    eval.RegisterBuiltInFunction("nanmedian", &OmlNanMedian,   FunctionMetaData(-2, 1, STATAN));

    return 1;
}

//------------------------------------------------------------------------------
// Computes uniform distribution probability density function values 
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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
// Computes uniform distribution probability density function values 
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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
// Computes uniform distribution inverse cumulative distribution function values
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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
// Computes normal distribution probability density function values [normpdf]
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
                                         sigma->M(), sigma->N(), true);

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
                                   mu->M(), mu->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
// Computes normal distribution cumulative distribution values [normcdf]
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
                                         sigma->M(), sigma->N(), true);

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
                                   mu->M(), mu->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
// Computes normal distribution inverse cumulative distribution values [norminv]
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
                                         sigma->M(), sigma->N(), true);

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
                                   mu->M(), mu->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
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
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, static_cast<int>(nargin));

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
                    hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
// Computes beta distribution probability distribution function values [betapdf]
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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), true);

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
                                   x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(x->M(), 
                                   x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(x->M(), 
                                   x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
// Computes beta distribution cumulative distribution values [betacdf]
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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
// Computes beta distribution inverse cumulative distribution values [betainv]
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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double    b      = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
// Computes incomplete gamma function values [gammainc]
//------------------------------------------------------------------------------
bool OmlGammaInc(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin != 2 && nargin != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> inputs2;
    inputs2.push_back(inputs[0]);
    inputs2.push_back(inputs[1]);
    inputs2.push_back(1.0);

    if (!OmlGamcdf(eval, inputs2, outputs))
        return false;

    if (nargin == 3)
    {
        if (!inputs[2].IsString())
            throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);

        std::string str = inputs[2].StringVal();

        if (str == "upper")
        {
            // switch tail
            if (outputs[0].IsScalar())
            {
                double value = outputs[0].Scalar();
                outputs[0] = 1.0 - value;
            }
            else if (outputs[0].IsMatrix())
            {
                hwMatrix* matrix = outputs[0].GetWritableMatrix();
                int size = matrix->Size();

                for (int i = 0; i < size; ++i)
                {
                    (*matrix)(i) = 1.0 - (*matrix)(i);
                }
            }
            else // ND matrix
            {
                hwMatrixN* matrix = outputs[0].GetWritableMatrixN();
                int size = matrix->Size();

                for (int i = 0; i < size; ++i)
                {
                    (*matrix)(i) = 1.0 - (*matrix)(i);
                }
            }
        }
        else if (str != "lower")
        {
            throw OML_Error(OML_ERR_OPTIONVAL, 3, OML_VAR_STRING);
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes exponential distribution probability density values[exppdf]
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
                                     a->M(), a->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               a->M(), a->N(), true);

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
// Computes exponential distribution cumulative distribution values [expcdf]
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
                                     a->M(), a->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               a->M(), a->N(), true);

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
// Computes exponential distribution inverse cumulative distribution values [expinv]
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
                                     a->M(), a->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               a->M(), a->N(), true);

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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
        hwMatrix*       result = EvaluatorInterface::allocateMatrix(m, n, true);
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
// Computes chi-squared distribution probability density values [chi2pdf]
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
                                     n->M(), n->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(n->M(), n->N(), true);

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
// Computes chi-squared distribution cumulative distribution values [chi2cdf]
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
                                     n->M(), n->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               n->M(), n->N(), true);

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
// Computes chi-squared distribution inverse cumulative distribution values [chi2inv]
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
                                     n->M(), n->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               n->M(), n->N(), true);

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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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

            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
// Computes Student t distribution probability density values [tpdf]
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
                                     n->M(), n->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               n->M(), n->N(), true);

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
// Computes Student t distribution cumulative distribution function values [tcdf]
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
                                     n->M(), n->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               n->M(), n->N(), true);

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
// Computes Student t distribution inverse cumulative distribution values [tinv]
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
                                     n->M(), n->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               n->M(), n->N(), true);

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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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

            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
                                         n->M(), n->N(), true);

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
                                   m->M(), m->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(m->M(), m->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int       n      = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                                         n->M(), n->N(), true);

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
                                   m->M(), m->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(m->M(), m->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int       n      = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                                         n->M(), n->N(), true);

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
                                   m->M(), m->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(m->M(), m->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                if (!cur3.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

                int       n      = static_cast<int>(cur3.Scalar());
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
            M = EvaluatorInterface::allocateMatrixI(m, n, mdof);
        }
        else if (cur1.IsMatrix() && cur1.Matrix()->IsReal())
        {
            const hwMatrix* temp = cur1.Matrix();
            int size = temp->Size();
            M = EvaluatorInterface::allocateMatrixI(m, n, true);

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
            N = EvaluatorInterface::allocateMatrixI(m, n, ndof);
        }
        else if (cur2.IsMatrix() && cur2.Matrix()->IsReal())
        {
            const hwMatrix* temp = cur2.Matrix();
            int size = temp->Size();
            N = EvaluatorInterface::allocateMatrixI(m, n, true);

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

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
                                         sigma->M(), sigma->N(), true);

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
                                   mu->M(), mu->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
// Computes lognormal distribution cumulative distribution values [logncdf]
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
                                         sigma->M(), sigma->N(), true);

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
                                   mu->M(), mu->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                                         sigma->M(), sigma->N(), true);

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
                                   mu->M(), mu->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   mu->M(), mu->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double sigma = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                                         b->M(), b->N(), true);

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
                                   a->M(), a->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 2, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   a->M(), a->N(), true);

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
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (cur3.IsScalar())
            {
                double b = cur3.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);

                hwMatrix* result = EvaluatorInterface::allocateMatrix(
                                   x->M(), x->N(), true);

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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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

        hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
                                     a->M(), a->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(
                               a->M(), a->N(), true);

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
                                     a->M(), a->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(a->M(), a->N(), true);

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
                                     a->M(), a->N(), true);

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
                               x->M(), x->N(), true);

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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(a->M(), a->N(), true);

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
                result.N(), true);

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
            result.N(), true);

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
// Computes inverse cumulative distribution function values [icdf]
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
// Computes cumulative distribution function values [cdf]
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
// Computes probability density function values [pdf]
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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
            hwMatrix*    result = EvaluatorInterface::allocateMatrix(m, n, true);
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
// Computes error function values [erf]
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
        hwMatrix* out = EvaluatorInterface::allocateMatrix(m->M(), m->N(), true);

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

    hwMatrix* result = EvaluatorInterface::allocateMatrix(1, max, true);
    hwMatrix* p      = EvaluatorInterface::allocateMatrix(1, max, true);
    hwMatrix* CI     = EvaluatorInterface::allocateMatrix(2, max, true);

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

        if (!mstat.IsOk())
        {
            if (mstat.GetArg1() == 6)
            {
                if (nargin > 3 && inputs[2].StringVal() == "alpha")
                    mstat.SetArg1(4);
            }
        }

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

    hwMatrix* result = EvaluatorInterface::allocateMatrix(1, max, true);
    hwMatrix* p      = EvaluatorInterface::allocateMatrix(1, max, true);
    hwMatrix* CI     = EvaluatorInterface::allocateMatrix(2, max, true);

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

        if (!mstat.IsOk())
        {
            if (mstat.GetArg1() == 6)
            {
                if (nargin > 3 && inputs[2].StringVal() == "alpha")
                    mstat.SetArg1(4);
            }
        }

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

    hwMatrix* result = EvaluatorInterface::allocateMatrix(1, max, true);
    hwMatrix* p      = EvaluatorInterface::allocateMatrix(1, max, true);
    hwMatrix* CI     = EvaluatorInterface::allocateMatrix(2, max, true);

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

        if (!mstat.IsOk())
        {
            if (mstat.GetArg1() == 6)
            {
                if (nargin > 3 && inputs[2].StringVal() == "alpha")
                    mstat.SetArg1(4);
            }
        }

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

    hwMatrix* result = EvaluatorInterface::allocateMatrix(1, max, true);
    hwMatrix* p      = EvaluatorInterface::allocateMatrix(1, max, true);
    hwMatrix* CI     = EvaluatorInterface::allocateMatrix(2, max, true);

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

        if (!mstat.IsOk())
        {
            if (mstat.GetArg1() == 6)
            {
                if (nargin > 3 && inputs[2].StringVal() == "alpha")
                    mstat.SetArg1(4);
            }
        }

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

    hwMatrix* result = EvaluatorInterface::allocateMatrix(1, max, true);
    hwMatrix* p      = EvaluatorInterface::allocateMatrix(1, max, true);
    hwMatrix* CI     = EvaluatorInterface::allocateMatrix(2, max, true);

    hwMatrix* piece    = EvaluatorInterface::allocateMatrix();
    hwMatrix* interval = EvaluatorInterface::allocateMatrix();

    for (int i = 0; i < max; ++i)
    {
        hwMathStatus mstat = (dim == 1) ? data->ReadColumn(i, *piece) :
                                          data->ReadRow(i, *piece);
        BuiltInFuncsUtils::CheckMathStatus(eval, mstat);

        mstat = ZTest(*piece, mu, sigma, reject, prob, *interval, alpha);

        if (!mstat.IsOk())
        {
            if (mstat.GetArg1() == 7)
            {
                if (nargin > 4 && inputs[3].StringVal() == "alpha")
                    mstat.SetArg1(5);
            }
        }

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
// Computes root mean square values [rms]
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
// Computes skewness values [skewness]
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
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlSkewness, 23);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    int dim;
    bool correctBias;

    if (!data->IsReal())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    if (nargin > 1)
    {
        if (inputs[1].IsScalar())
        {
            if (inputs[1].Scalar() == 0.0)
                correctBias = true;
            else if (inputs[1].Scalar() == 1.0)
                correctBias = false;
            else
                throw OML_Error(OML_ERR_FLAG_01, 2);
        }
        else if (inputs[1].IsEmpty())
        {
            correctBias = false;  // default is opposite of var and std
        }
        else
        {
            throw OML_Error(OML_ERR_FLAG_01, 2);
        }
    }
    else
    {
        correctBias = false;  // default is opposite of var and std
    }

    if (nargin > 2)
    {
        if (!inputs[2].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 3);

        dim = (int)inputs[2].Scalar();
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
        double skew;
        status = Skewness(*data, skew, correctBias);
        outputs.push_back(skew);
    }
    else if (dim == 1)
    {
        const double* colPtr;
        hwMatrix* skew;

        if (data->M())
        {
            skew = new hwMatrix(1, data->N(), hwMatrix::REAL);

            for (int i = 0; i < data->N(); ++i)
            {
                colPtr = &((*data)(0, i));
                hwMatrix col(data->M(), (void*)colPtr, hwMatrix::REAL);
                status = Skewness(col, (*skew)(0, i), correctBias);
            }
        }
        else
        {
            if (data->N())
                skew = new hwMatrix(1, data->N(), hwMatrix::REAL);
            else
                skew = new hwMatrix(1, 1, hwMatrix::REAL);

            skew->SetElements(std::numeric_limits<double>::quiet_NaN());
        }

        outputs.push_back(skew);
    }
    else if (dim == 2)
    {
        hwMatrix row;
        hwMatrix* skew;

        if (data->N())
        {
            skew = new hwMatrix(data->M(), 1, hwMatrix::REAL);

            for (int i = 0; i < data->M(); ++i)
            {
                status = data->ReadRow(i, row);
                status = Skewness(row, (*skew)(i, 0), correctBias);
            }
        }
        else
        {
            if (data->M())
                skew = new hwMatrix(data->M(), 1, hwMatrix::REAL);
            else
                skew = new hwMatrix(1, 1, hwMatrix::REAL);

            skew->SetElements(std::numeric_limits<double>::quiet_NaN());
        }

        outputs.push_back(skew);
    }
    else
    {
        hwMatrix* skew = new hwMatrix(data->M(), data->N(), hwMatrix::REAL);
        skew->SetElements(std::numeric_limits<double>::quiet_NaN());

        outputs.push_back(skew);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes kurtosis values [kurtosis]
//------------------------------------------------------------------------------
bool OmlKurtosis(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsNDMatrix())
    {
        if (nargin < 2)
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlKurtosis);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlKurtosis, 23);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    int dim;
    bool correctBias;

    if (!data->IsReal())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    if (nargin > 1)
    {
        if (inputs[1].IsScalar())
        {
            if (inputs[1].Scalar() == 0.0)
                correctBias = true;
            else if (inputs[1].Scalar() == 1.0)
                correctBias = false;
            else
                throw OML_Error(OML_ERR_FLAG_01, 2);
        }
        else if (inputs[1].IsEmpty())
        {
            correctBias = false;  // default is opposite of var and std
        }
        else
        {
            throw OML_Error(OML_ERR_FLAG_01, 2);
        }
    }
    else
    {
        correctBias = false;  // default is opposite of var and std
    }

    if (nargin > 2)
    {
        if (!inputs[2].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 3);

        dim = (int)inputs[2].Scalar();
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
        double kurt;
        status = Kurtosis(*data, kurt, correctBias);
        outputs.push_back(kurt);
    }
    else if (dim == 1)
    {
        const double* colPtr;
        hwMatrix* kurt;

        if (data->M())
        {
            kurt = new hwMatrix(1, data->N(), hwMatrix::REAL);

            for (int i = 0; i < data->N(); ++i)
            {
                colPtr = &((*data)(0, i));
                hwMatrix col(data->M(), (void*)colPtr, hwMatrix::REAL);
                status = Kurtosis(col, (*kurt)(0, i), correctBias);
            }
        }
        else
        {
            if (data->N())
                kurt = new hwMatrix(1, data->N(), hwMatrix::REAL);
            else
                kurt = new hwMatrix(1, 1, hwMatrix::REAL);

            kurt->SetElements(std::numeric_limits<double>::quiet_NaN());
        }

        outputs.push_back(kurt);
    }
    else if (dim == 2)
    {
        hwMatrix row;
        hwMatrix* kurt;

        if (data->N())
        {
            kurt = new hwMatrix(data->M(), 1, hwMatrix::REAL);

            for (int i = 0; i < data->M(); ++i)
            {
                status = data->ReadRow(i, row);
                status = Kurtosis(row, (*kurt)(i, 0), correctBias);
            }
        }
        else
        {
            if (data->M())
                kurt = new hwMatrix(data->M(), 1, hwMatrix::REAL);
            else
                kurt = new hwMatrix(1, 1, hwMatrix::REAL);

            kurt->SetElements(std::numeric_limits<double>::quiet_NaN());
        }

        outputs.push_back(kurt);
    }
    else
    {
        hwMatrix* kurt = new hwMatrix(data->M(), data->N(), hwMatrix::REAL);
        kurt->SetElements(std::numeric_limits<double>::quiet_NaN());

        outputs.push_back(kurt);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes variance values [var]
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
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlVariance, 23);
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
            throw OML_Error(OML_ERR_POSINTEGER, 3);

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
// Computes standard deviation values [std]
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
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlStd, 23);
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
            throw OML_Error(OML_ERR_POSINTEGER, 3);

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
// Computes median values [median]
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
// Computes quantile values [quantile]
//------------------------------------------------------------------------------
bool OmlQuantile(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    const hwMatrix* P = nullptr;

    if (nargin > 1)
    {
        if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
            throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_TYPE);

        P = inputs[1].ConvertToMatrix();
    }

    if (nargin == 1 || P->Is0x0())
    {
        // default quantile request
        hwMatrix* P = EvaluatorInterface::allocateMatrix(1, 5, true);
        (*P)(0) = 0.0;
        (*P)(1) = 0.25;
        (*P)(2) = 0.50;
        (*P)(3) = 0.75;
        (*P)(4) = 1.0;

        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);
        inputs2.push_back(P);

        if (nargin > 1)
        {
            for (std::vector<Currency>::const_iterator it = inputs.begin() + 2; it != inputs.end(); it++)
                inputs2.push_back(*it);
        }

        return OmlQuantile(eval, inputs2, outputs);
    }

    int dim = -1;

    if (nargin > 2)
    {
        if (inputs[2].IsPositiveInteger())
        {
            dim = static_cast<int> (inputs[2].Scalar());
        }
        else if (inputs[2].IsMatrix())
        {
            if (!inputs[2].Matrix()->Is0x0())
                throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DATA);
        }
        else
        {
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DATA);
        }
    }

    if (dim == -1)
    {
        if (inputs[0].IsMatrix() || inputs[0].IsScalar())
        {
            const hwMatrix* mtx = inputs[0].ConvertToMatrix();

            if (!mtx->IsReal())
                throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_TYPE);

            if (mtx->M() == 1)
                dim = 2;
            else
                dim = 1;
        }
        else if (inputs[0].IsNDMatrix())
        {
            const hwMatrixN* mtx = inputs[0].MatrixN();
            const std::vector<int>& dims = mtx->Dimensions();

            if (!mtx->IsReal())
                throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_TYPE);

            for (int i = 0; i < dims.size(); ++i)
            {
                if (dims[i] != 1)
                {
                    dim = i + 1;
                    break;
                }
            }
        }
        else
        {
            throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
        }
    }

    // ensure that inputs[0] is sorted
    static bool unsorted = true;

    if (dim != 1)
    {
        // permute matrix
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);
        oml_ndims(eval, inputs2, outputs);
        int numDim = static_cast<int>(outputs[0].Scalar());

        if (dim < 1 || dim > numDim)
        {
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);
        }

        outputs.clear();
        hwMatrix* permvec = new hwMatrix(numDim, 1, hwMatrix::REAL);

        for (int i = 1; i < numDim; ++i)
            (*permvec)(i) = i + 1;

        (*permvec)(0) = dim;
        (*permvec)(dim - 1) = 1;
        Currency permCur(permvec);
        inputs2.push_back(permCur);
        oml_permute(eval, inputs2, outputs);

        // compute permuted quantiles
        inputs2.clear();
        inputs2.push_back(outputs[0]);  // permuted data
        inputs2.push_back(inputs[1]);   // percentiles
        inputs2.push_back(1);           // dimension 1

        if (nargin == 4)
            inputs2.push_back(inputs[3]);   // method

        outputs.clear();

        OmlQuantile(eval, inputs2, outputs);

        if (unsorted && inputs2[0].Matrix()->IsVector())
        {
            // unsorted has been reset, so preparing to exit
            // check special vector case
            hwMatrix* Q = outputs[0].GetWritableMatrix();

            if (P->M() != Q->M())
                Q->Transpose();
        }
        else
        {
            // ipermute quantiles
            inputs2.clear();
            inputs2.push_back(outputs[0]);
            inputs2.push_back(permCur);
            outputs.clear();
            oml_ipermute(eval, inputs2, outputs);
        }

        return true;
    }

    if (inputs[0].IsNDMatrix())
    {
        const hwMatrixN* ndmat = inputs[0].MatrixN();
        std::vector<int> dims = ndmat->Dimensions();

        // reshape to 2D matrix
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);
        inputs2.push_back(dims[0]);
        inputs2.push_back(Currency());
        oml_reshape(eval, inputs2, outputs);

        // compute quantiles
        inputs2.clear();
        inputs2.push_back(outputs[0]);

        for (std::vector<Currency>::const_iterator it = inputs.begin() + 1; it != inputs.end(); it++)
            inputs2.push_back(*it);

        outputs.clear();
        OmlQuantile(eval, inputs2, outputs);

        if (unsorted && inputs2[0].Matrix()->IsVector())
        {
            // unsorted has been reset, so preparing to exit
            // check special vector case
            hwMatrix* Q = outputs[0].GetWritableMatrix();

            if (P->M() != Q->M())
                Q->Transpose();
        }
        else
        {
            // reshape quantiles
            const hwMatrix* P = inputs[1].ConvertToMatrix();
            hwMatrix* Q = outputs[0].GetWritableMatrix();

            inputs2.clear();
            inputs2.push_back(outputs[0]);  // Q
            dims[0] = P->Size();

            for (int i = 0; i < dims.size(); ++i)
                inputs2.push_back(dims[i]);

            outputs.clear();
            oml_reshape(eval, inputs2, outputs);
        }

        return true;
    }

    if (unsorted)
    {
        // sort x
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);           // push_back(x)
        inputs2.push_back(dim);
        inputs2.push_back("ascend");
        outputs.clear();
        oml_sort(eval, inputs2, outputs);   // sort(x)
        inputs2.clear();
        inputs2.push_back(outputs[0]);
        outputs.clear();

        for (std::vector<Currency>::const_iterator it = inputs.begin()+1; it != inputs.end(); it++)
            inputs2.push_back(*it);

        try
        {
            unsorted = false;
            return OmlQuantile(eval, inputs2, outputs);
        }
        catch (OML_Error&)
        {
            unsorted = true;    // reset
            throw;
        }
        catch (hwMathException&)
        {
            unsorted = true;    // reset
            throw;
        }

        outputs.clear();
    }

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
    {
        unsorted = true;    // reset
        throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DATA);
    }

    const hwMatrix* data = inputs[0].ConvertToMatrix();
    int method = 5;

    if (nargin > 3)
    {
        if (!inputs[3].IsPositiveInteger())
        {
            unsorted = true;    // reset
            throw OML_Error(OML_ERR_POSINTEGER, 4, OML_VAR_DIM);
        }

        method = static_cast<int> (inputs[3].Scalar());
    }

    hwMatrix* Q = EvaluatorInterface::allocateMatrix();
    hwMathStatus status = Quantile(*data, *P, method, *Q);

    unsorted = true;    // reset
    BuiltInFuncsUtils::CheckMathStatus(eval, status);
    outputs.push_back(Q);

    return true;
}

void SortedLookup(const double* edges, int numEdges, const double* values, int numValues,
                  double* count, double* idx, int stride)
{
    int idx_e = 0;
    int idx_v = 0;
    int mem_v = 0;
    int mem_c = 0;

    while (idx_v < numValues && idx_e < numEdges - 1)
    {
        if (values[mem_v] < edges[0])
        {
            if (idx)
            {
                idx[mem_v] = 0;
            }

            ++idx_v;
            mem_v += stride;
            continue;
        }

        break;
    }

    while (idx_v < numValues && idx_e < numEdges - 1)
    {
        if (values[mem_v] < edges[idx_e + 1])
        {
            if (idx)
            {
                idx[mem_v] = idx_e + 1;
            }

            ++idx_v;
            mem_v += stride;
            ++count[mem_c];
            continue;
        }

        ++idx_e;
        mem_c += stride;
    }

    while (idx_v < numValues && idx_e == numEdges - 1)
    {
        if (values[mem_v] == edges[idx_e])
        {
            if (idx)
            {
                idx[mem_v] = numEdges;
            }

            ++idx_v;
            mem_v += stride;
            ++count[mem_c];
            continue;
        }

        break;
    }

    while (idx && idx_v < numValues)
    {
        idx[mem_v] = 0;
        ++idx_v;
        mem_v += stride;
    }
}

//------------------------------------------------------------------------------
// Compute histogram counts [histc]
//------------------------------------------------------------------------------
bool OmlHistC(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    int nargout = eval.GetNargoutValue();

    if (nargout > 2)
        throw OML_Error(OML_ERR_NUMARGOUT);

    // get edges vector
    if (!inputs[1].IsMatrix())
    {
        throw OML_Error(OML_ERR_MATRIX, 2);
    }

    const hwMatrix* edges = inputs[1].Matrix();

    if (!edges->IsReal())
    {
        throw OML_Error(OML_ERR_REAL, 2);
    }

    if (!edges->IsVector())
    {
        throw OML_Error(OML_ERR_VECTOR, 2);
    }

    int numEdges = edges->Size();

    // get dimension input
    int dim = -1;

    if (nargin == 3)
    {
        if (inputs[2].IsPositiveInteger())
            dim = static_cast<int> (inputs[2].Scalar()) - 1;
    }

    // sort data
    std::vector<Currency> inputs2;
    std::vector<Currency> outputs2;
    inputs2.push_back(inputs[0]);

    if (dim != -1)
    {
        inputs2.push_back(dim + 1);
    }

    inputs2.push_back("ascend");
    oml_sort(eval, inputs2, outputs2);   // sort(x)

    if (outputs2[0].IsMatrix())
    {
        const hwMatrix* data = outputs2[0].Matrix();

        if (!data->IsReal())
        {
            throw OML_Error(OML_ERR_REAL, 1);
        }

        if (dim == -1)
        {
            // first non-singleton
            if (data->M() != 1 || data->N() == 1)
                dim = 0;
            else
                dim = 1;
        }
        else if (dim != 0 && dim != 1)
        {
            throw OML_Error(OML_ERR_INVALID_RANGE, 3, OML_VAR_DIM);
        }

        hwMatrix* count = EvaluatorInterface::allocateMatrix();
        hwMatrix* indx = nullptr;
        double* indxVec = nullptr;
        int m = data->M();
        int n = data->N();

        if (dim == 0)
            count->Dimension(numEdges, n, hwMatrix::REAL);
        else
            count->Dimension(m, numEdges, hwMatrix::REAL);

        count->SetElements(0.0);

        if (nargout == 2)
        {
            indx = new hwMatrix(m, n, hwMatrix::REAL);
        }

        std::vector<int> dims = { m, n };
        int numVecs = (dim == 0) ? n : m;
        int stride = (dim == 0) ? 1 : m;
        int row = 0;
        int col = 0;

        for (int i = 0; i < numVecs; ++i)
        {
            // set the matrix indices to the first index in each slice
            int start = col * m + row;
            const double* real = data->GetRealData() + start;

            if (nargout == 2)
            {
                indxVec = indx->GetRealData() + start;
            }

            start = col * count->M() + row;
            double* count_vec = count->GetRealData() + start;

            // perform op
            SortedLookup(edges->GetRealData(), numEdges, real, dims[dim],
                         count_vec, indxVec, stride);

            // advance slice indices
            if (dim == 0)
                ++col;
            else
                ++row;
        }

        outputs.push_back(count);

        if (nargout == 2)
        {
            // unsort the indices
            hwMatrix* sortIndx = outputs2[1].GetWritableMatrix();
            hwMatrix* indx2 = new hwMatrix(m, n, hwMatrix::REAL);

            if (dim == 0)
            {
                for (int j = 0; j < n; ++j)
                {
                    for (int i = 0; i < m; ++i)
                    {
                        (*indx2)(static_cast<int>((*sortIndx)(i, j)) - 1, j) = (*indx)(i, j);
                    }
                }
            }
            else
            {
                for (int j = 0; j < n; ++j)
                {
                    for (int i = 0; i < m; ++i)
                    {
                        (*indx2)(i, static_cast<int>((*sortIndx)(i, j)) - 1) = (*indx)(i, j);
                    }
                }
            }

            delete indx;
            outputs.push_back(indx2);
        }
    }
    else if (outputs2[0].IsNDMatrix())
    {
        const hwMatrixN* data = outputs2[0].MatrixN();

        if (!data->IsReal())
        {
            throw OML_Error(OML_ERR_REAL, 1);
        }

        const std::vector<int>& dims = data->Dimensions();
        int numDim = static_cast<int> (dims.size());

        if (dim == -1)
        {
            // first non-singleton
            for (int i = 0; i < numDim; ++i)
            {
                if (dims[i] != 1)
                {
                    dim = i;
                    break;
                }
            }
        }
        else if (dim < 0 || dim > numDim - 1)
        {
            throw OML_Error(OML_ERR_INVALID_RANGE, 3, OML_VAR_DIM);
        }

        hwMatrixN* count = EvaluatorInterface::allocateMatrixN();
        hwMatrixN* indx = nullptr;
        double* indxVec = nullptr;
        std::vector<int> countDims = dims;

        countDims[dim] = numEdges;
        count->Dimension(countDims, hwMatrixN::REAL);
        count->SetElements(0.0);

        if (nargout == 2)
        {
            indx = new hwMatrixN(dims, hwMatrixN::REAL);
        }

        int numVecs = data->Size() / dims[dim];
        int stride = data->Stride(dim);
        std::vector<int> matrixIndex(numDim);

        for (int i = 0; i < numVecs; ++i)
        {
            // set the matrix indices to the first index in each slice
            int start = data->Index(matrixIndex);
            const double* real = data->GetRealData() + start;

            if (nargout == 2)
            {
                indxVec = indx->GetRealData() + start;
            }

            start = count->Index(matrixIndex);
            double* count_vec = count->GetRealData() + start;

            // perform op
            SortedLookup(edges->GetRealData(), numEdges, real, dims[dim],
                         count_vec, indxVec, stride);

            // advance slice indices
            for (int j = 0; j < numDim; ++j)
            {
                if (j == dim)
                    continue;

                // increment index j if possible
                if (matrixIndex[j] < static_cast<int> (dims[j]) - 1)
                {
                    ++matrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                matrixIndex[j] = 0;
            }
        }

        outputs.push_back(count);

        if (nargout == 2)
        {
            // unsort the indices
            hwMatrixN* sortIndx = outputs2[1].GetWritableMatrixN();
            hwMatrixN* indx2 = new hwMatrixN(dims, hwMatrixN::REAL);
            int size = indx->Size();

            for (int i = 0; i < size; ++i)
            {
                double rhs = (*indx)(matrixIndex);
                int save_idx = matrixIndex[dim];
                int sort_idx = static_cast<int>((*sortIndx)(matrixIndex)) - 1;

                matrixIndex[dim] = sort_idx;
                (*indx2)(matrixIndex) = rhs;
                matrixIndex[dim] = save_idx;

                // advance slice indices
                for (int j = 0; j < numDim; ++j)
                {
                    // increment index j if possible
                    if (matrixIndex[j] < static_cast<int> (dims[j]) - 1)
                    {
                        ++matrixIndex[j];
                        break;
                    }

                    // index j is maxed out, so reset and continue to j+1
                    matrixIndex[j] = 0;
                }
            }

            delete indx;
            outputs.push_back(indx2);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes mean absolute deviation values [meandev]
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
// Computes mean or median absolute deviation values [mad]
//------------------------------------------------------------------------------
bool OmlMAD(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs,
            std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    int opt = 0;

    if (nargin > 1)
    {
        if (!inputs[1].IsInteger())
            throw OML_Error(OML_ERR_FLAG_01, 2, OML_VAR_OPTION);

        opt = static_cast<int> (inputs[1].Scalar());

        if (opt != 0 && opt != 1)
            throw OML_Error(OML_ERR_FLAG_01, 2, OML_VAR_OPTION);
    }

    if (opt == 0)
    {
        std::vector<Currency> inputs2;
        inputs2.push_back(inputs[0]);

        if (nargin > 2)
            inputs2.push_back(inputs[2]);

        return OmlMeandev(eval, inputs2, outputs);
    }

    if (inputs[0].IsNDMatrix())
    {
        if (nargin == 1)
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlMAD);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, OmlMAD, 23);
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* data = inputs[0].ConvertToMatrix();

    if (data->IsReal() && !data->IsRealData())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    int dim;

    if (nargin > 2)
    {
        if (!inputs[2].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);

        dim = static_cast<int> (inputs[2].Scalar());
    }
    else
    {
        dim = data->M() == 1 ? 2 : 1;
    }

    outputs.push_back(oml_MatrixUtil(eval, data, dim, &callOnVector<&MedianDev>));
    return true;
}
//------------------------------------------------------------------------------
// Computes mean values [mean]
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
        if (inputs[1].IsPositiveInteger())
            dim = static_cast<int>(inputs[1].Scalar());
        else if (inputs[1].IsMatrix() && inputs[1].Matrix()->Is0x0())
            dim = data->M() == 1 ? 2 : 1;
        else
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DIM);
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
        else if (outputs[0].IsScalar())
        {
            double mean = outputs[0].Scalar();
            outputs.clear();
            outputs.push_back(mean / data->Size());
        }
        else if (outputs[0].IsComplex())
        {
            hwComplex mean = outputs[0].Complex();
            outputs.clear();
            outputs.push_back(mean / data->Size());
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
// Mode helper
//------------------------------------------------------------------------------
void ModeHelper(double* data, int n, int stride, double* modeData)
{
    int idx = 0;
    modeData[0] = 1;

    for (int i = stride; i < n * stride; i += stride)
    {
        if (data[i] == data[i - stride])
        {
            ++modeData[idx];
        }
        else
        {
            idx = i;
            modeData[idx] = 1;
        }
    }
}
//------------------------------------------------------------------------------
// Computes mode values [mode]
//------------------------------------------------------------------------------
bool OmlMode(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs,
             std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    int nargout = eval.GetNargoutValue();

    if (nargin != 1 && nargin != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargout > 3)
        throw OML_Error(OML_ERR_NUMARGOUT);

    // get dimension on which to operate
    int dim;

    if (nargin == 2)
    {
        if (!inputs[1].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

        dim = static_cast<int> (inputs[1].Scalar());
    }
    else if (inputs[0].IsMatrix() || inputs[0].IsScalar())
    {
        const hwMatrix* mtx = inputs[0].ConvertToMatrix();

        if (!mtx->IsReal())
            throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_TYPE);

        if (mtx->M() == 1)
            dim = 2;
        else
            dim = 1;
    }
    else if (inputs[0].IsNDMatrix())
    {
        const hwMatrixN* mtx = inputs[0].MatrixN();
        const std::vector<int>& dims = mtx->Dimensions();

        if (!mtx->IsReal())
            throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_TYPE);

        for (int i = 0; i < dims.size(); ++i)
        {
            if (dims[i] != 1)
            {
                dim = i;
                break;
            }
        }
    }
    else
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
    }

    // sort
    std::vector<Currency> outputs2;
    oml_sort(eval, inputs, outputs2);

    Currency& cur = outputs2[0];

    // compute based on matrix type
    int stride;
    int length;
    int vecDelta;
    int numVecs;

    std::vector<Currency> inputsMax;
    std::vector<Currency> outputsMax;

    if (cur.IsMatrix() || cur.IsScalar())
    {
        // collect mode data
        hwMatrix* mtx = cur.GetWritableMatrix();
        double* data = mtx->GetRealData();
        int m = mtx->M();
        int n = mtx->N();
        hwMatrix* modeData = EvaluatorInterface::allocateMatrix(m, n, true);
        modeData->SetElements(0.0);
        double* mData = modeData->GetRealData();

        if (dim == 1)
        {
            numVecs = n;
            vecDelta = m;
            length = m;
            stride = 1;
        }
        else if (dim == 2)
        {
            numVecs = m;
            vecDelta = 1;
            length = n;
            stride = m;
        }
        else
        {
            outputs.push_back(inputs[0]);

            hwMatrix* onesM = EvaluatorInterface::allocateMatrix(m, n, true);
            onesM->SetElements(1.0);
            outputs.push_back(onesM);

            std::vector<Currency> outputs2;
            BuiltInFuncsData::Num2Cell(eval, inputs, outputs2);
            outputs.push_back(outputs2[0]);

            return true;
        }

        for (int i = 0; i < numVecs; ++i)
        {
            ModeHelper(data, length, stride, mData);
            data += vecDelta;
            mData += vecDelta;
        }

        // compute mode frequency
        inputsMax.push_back(modeData);
        inputsMax.push_back(EvaluatorInterface::allocateMatrix());
        inputsMax.push_back(dim);
        outputsMax = eval.DoMultiReturnFunctionCall(oml_max, inputsMax, 3, 1, true);

        Currency& maxC = outputsMax[0];

        // populate mode and cell outputs
        hwMatrix* mode = nullptr;
        HML_CELLARRAY* cArray = nullptr;

        if (dim == 1)
        {
            mode = EvaluatorInterface::allocateMatrix(1, n, true);
            cArray = EvaluatorInterface::allocateCellArray(1, n);
        }
        else    // dim == 2
        {
            mode = EvaluatorInterface::allocateMatrix(m, 1, true);
            cArray = EvaluatorInterface::allocateCellArray(m, 1);
        }

        data = mtx->GetRealData();
        mData = modeData->GetRealData();

        for (int i = 0; i < numVecs; ++i)
        {
            hwMatrix* modeVals = EvaluatorInterface::allocateMatrix(1, 1, true);
            (*cArray)(i) = Currency(modeVals);
            int numModes = 0;
            double maxVal;

            if (maxC.IsMatrix())
            {
                maxVal = (*maxC.GetWritableMatrix())(i);
            }
            else    // scalar
            {
                maxVal = maxC.Scalar();
            }

            for (int k = 0; k < length * stride; k += stride)
            {
                if (mData[k] == maxVal)
                {
                    ++numModes;
                    hwMathStatus status = modeVals->Resize(numModes, 1);
                    (*modeVals)(numModes - 1) = data[k];
                }
            }

            data += vecDelta;
            mData += vecDelta;
            (*mode)(i) = (*modeVals)(0);
        }

        outputs.push_back(mode);
        outputs.push_back(outputsMax[0]);
        outputs.push_back(cArray);
    }
    else // if (inputs[0].IsNDMatrix())
    {
        // collect mode data
        hwMatrixN* mtx = cur.GetWritableMatrixN();
        const std::vector<int>& dims = mtx->Dimensions();
        int numDim = static_cast<int> (dims.size());
        hwMatrixN* modeData = new hwMatrixN(dims, hwMatrixN::REAL);
        modeData->SetElements(0.0);

        if (dim > dims.size())
        {
            outputs.push_back(inputs[0]);

            hwMatrixN* onesM = EvaluatorInterface::allocateMatrixN(dims, true);
            onesM->SetElements(1.0);
            outputs.push_back(onesM);

            std::vector<Currency> outputs2;
            BuiltInFuncsData::Num2Cell(eval, inputs, outputs2);
            outputs.push_back(outputs2[0]);

            return true;
        }

        --dim;  // switch to zero based
        length = dims[dim];
        numVecs = mtx->Size() / length;
        stride = mtx->Stride(dim);

        std::vector<int> rhsMatrixIndex(numDim);

        for (int i = 0; i < numVecs; ++i)
        {
            // set the rhsMatrix indices to the first index in each slice
            int start = mtx->Index(rhsMatrixIndex);
            double* data = mtx->GetRealData() + start;
            double* mData = modeData->GetRealData() + start;

            ModeHelper(data, length, stride, mData);

            // advance slice indices
            for (int j = 0; j < numDim; ++j)
            {
                if (j == dim)
                    continue;

                // increment index j if possible
                if (rhsMatrixIndex[j] < static_cast<int> (dims[j]) - 1)
                {
                    ++rhsMatrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                rhsMatrixIndex[j] = 0;
            }
        }

        // compute mode frequency
        inputsMax.push_back(modeData);
        inputsMax.push_back(EvaluatorInterface::allocateMatrix());
        inputsMax.push_back(dim + 1);
        outputsMax = eval.DoMultiReturnFunctionCall(oml_max, inputsMax, 3, 1, true);
        Currency& maxC = outputsMax[0];

        // populate mode and cell outputs
        std::vector<int> modeDims = dims;
        modeDims[dim] = 1;

        hwMatrixN* mode = new hwMatrixN(modeDims, hwMatrixN::REAL);
        HML_ND_CELLARRAY* cArray = EvaluatorInterface::allocateNDCellArray(modeDims);

        rhsMatrixIndex.clear();
        rhsMatrixIndex.resize(numDim);

        for (int i = 0; i < numVecs; ++i)
        {
            // set the rhsMatrix indices to the first index in each slice
            int start = mtx->Index(rhsMatrixIndex);
            double* data = mtx->GetRealData() + start;
            double* mData = modeData->GetRealData() + start;

            hwMatrix* modeVals = EvaluatorInterface::allocateMatrix(1, 1, true);
            (*cArray)(i) = Currency(modeVals);
            int numModes = 0;
            double maxVal;

            if (maxC.IsNDMatrix())
            {
                maxVal = (*maxC.GetWritableMatrixN())(i);
            }
            else if (maxC.IsMatrix())
            {
                maxVal = (*maxC.GetWritableMatrix())(i);
            }
            else    // scalar
            {
                maxVal = maxC.Scalar();
            }

            for (int k = 0; k < length * stride; k += stride)
            {
                if (mData[k] == maxVal)
                {
                    ++numModes;
                    hwMathStatus status = modeVals->Resize(numModes, 1);
                    (*modeVals)(numModes - 1) = data[k];
                }
            }

            (*mode)(i) = (*modeVals)(0);

            // advance slice indices
            for (int j = 0; j < numDim; ++j)
            {
                if (j == dim)
                    continue;

                // increment index j if possible
                if (rhsMatrixIndex[j] < static_cast<int> (dims[j]) - 1)
                {
                    ++rhsMatrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                rhsMatrixIndex[j] = 0;
            }
        }

        outputs.push_back(mode);
        outputs.push_back(outputsMax[0]);
        outputs.push_back(cArray);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes moving mean values [movmean]
//------------------------------------------------------------------------------
bool OmlMovMean(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin == 4 || nargin > 5)
        throw OML_Error(OML_ERR_NUMARGIN);

    int dim = -1;
    int na = -1;
    int nb = -1;

    if (nargin > 1)
    {
        if (inputs[1].IsPositiveInteger())
        {
            int wlen = static_cast<int> (inputs[1].Scalar());

            if (wlen % 2 == 0)
            {
                nb = wlen / 2;
                na = nb - 1;
            }
            else
            {
                na = nb = (wlen - 1) / 2;
            }
        }
        else if (inputs[1].IsVector())
        {
            const hwMatrix* data = inputs[1].Matrix();

            if (data->Size() != 2)
            {
            }

            nb = static_cast<int> ((*data)(0));
            na = static_cast<int> ((*data)(1));
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARVECTOR, 2, OML_VAR_DIM);
        }
    }

    if (nargin > 2)
    {
        if (inputs[2].IsPositiveInteger())
            dim = static_cast<int>(inputs[2].Scalar());
        else if (!inputs[2].IsMatrix() || !inputs[2].Matrix()->Is0x0())
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);
    }

    if (dim == -1)
    {
        if (inputs[0].IsMatrix() || inputs[0].IsScalar() || inputs[0].IsComplex())
        {
            const hwMatrix* mtx = inputs[0].ConvertToMatrix();

            if (mtx->M() == 1)
                dim = 2;
            else
                dim = 1;
        }
        else if (inputs[0].IsNDMatrix())
        {
            const hwMatrixN* mtx = inputs[0].MatrixN();
            const std::vector<int>& dims = mtx->Dimensions();

            for (int i = 0; i < dims.size(); ++i)
            {
                if (dims[i] != 1)
                {
                    dim = i + 1;
                    break;
                }
            }
        }
    }

    std::string endproperty = "shrink";
    hwComplex userVal(std::numeric_limits<double>::quiet_NaN(), 0.0);

    if (nargin == 5)
    {
        if (!inputs[3].IsString())
            throw OML_Error(OML_ERR_STRING, 4, OML_VAR_TYPE);

        if (inputs[3].StringVal() != "Endpoints")
            throw OML_Error(OML_ERR_OPTION, 4);

        if (inputs[4].IsString())
        {
            endproperty = inputs[4].StringVal();
        }
        else if (inputs[4].IsScalar())
        {
            endproperty = "userval";
            userVal = inputs[4].Scalar();
        }
        else if (inputs[4].IsComplex())
        {
            endproperty = "userval";
            userVal = inputs[4].Complex();
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARSTRING, 5, OML_VAR_TYPE);
        }
    }

    if (inputs[0].IsMatrix() || inputs[0].IsScalar() || inputs[0].IsComplex())
    {
        const hwMatrix* matrix = inputs[0].ConvertToMatrix();
        hwMatrixN matrixN;
        matrixN.Convert2DtoND(*matrix, false);

        hwMatrixN xBarN;
        hwMathStatus status = MovMean(matrixN, nb, na, dim-1, endproperty, userVal, xBarN);

        if (!status.IsOk())
        {
            int arg1 = status.GetArg1();

            if (arg1 > 2 && arg1 < 5)
                status.SetArg1(arg1 - 1);

            BuiltInFuncsUtils::CheckMathStatus(eval, status);
        }

        hwMatrix* xBar = EvaluatorInterface::allocateMatrix();
        xBarN.ConvertNDto2D(*xBar, false);
        outputs.push_back(xBar);
    }
    else if (inputs[0].IsNDMatrix())
    {
        const hwMatrixN* matrix = inputs[0].MatrixN();
        std::unique_ptr<hwMatrixN> xBar(EvaluatorInterface::allocateMatrixN());

        hwMathStatus status = MovMean(*matrix, nb, na, dim-1, endproperty, userVal, *xBar);

        if (!status.IsOk())
        {
            int arg1 = status.GetArg1();

            if (arg1 > 2 && arg1 < 5)
                status.SetArg1(arg1 - 1);

            BuiltInFuncsUtils::CheckMathStatus(eval, status);
        }

        outputs.push_back(xBar.release());
    }
    else
    {
        throw OML_Error(OML_ERR_MATRIX, 1);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes covariances [cov]
//------------------------------------------------------------------------------
bool OmlCov(EvaluatorInterface           eval,
            const std::vector<Currency>& inputs, 
            std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_DATA);

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
            throw OML_Error(OML_ERR_SCALARMATRIX, 2);
        }
    }
    else // nargin == 3
    {
        if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
            throw OML_Error(OML_ERR_REALMATRIX, 2, OML_VAR_DATA);

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
// Computes correlation coefficients [corr]
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
// Utility reassigns NaN values to an input value
//------------------------------------------------------------------------------
void NanAssignmentUtil(hwMatrix& data, double value, hwMatrixI* count,
                       int& dim)
{
    if (dim == -1)
    {
        // first non-singleton
        if (data.M() != 1 || data.N() == 1)
            dim = 0;
        else
            dim = 1;
    }
    else if (dim != 0 && dim != 1)
    {
        return;
    }

    if (count)
    {
        if (dim == 0)
            count->Dimension(1, data.N(), hwMatrixI::REAL);
        else
            count->Dimension(data.M(), 1, hwMatrixI::REAL);
    }

    int m = data.M();
    int n = data.N();
    std::vector<int> dims = { m, n };
    int numVecs = (dim == 0) ? n : m;
    int stride = (dim == 0) ? 1 : m;
    int row = 0;
    int col = 0;

    for (int i = 0; i < numVecs; ++i)
    {
        // set the matrix indices to the first index in each slice
        int start = col * m + row;
        int total = dims[dim];  // total non-NaN elements

        if (data.IsReal())
        {
            double* real = data.GetRealData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                // replace and count NaN values
                if (IsNaN_T(*real))
                {
                    *real = value;
                    --total;
                }

                real += stride;
            }
        }
        else
        {
            hwComplex* cplx = data.GetComplexData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                // replace and count NaN values
                if (IsNaN_T(*cplx))
                {
                    *cplx = value;
                    --total;
                }

                cplx += stride;
            }
        }

        if (count)
            (*count)(row, col) = total;

        // advance slice indices
        if (dim == 0)
            ++col;
        else
            ++row;
    }
}
//------------------------------------------------------------------------------
// Utility reassigns NaN values to an input value
//------------------------------------------------------------------------------
void NanAssignmentUtil(hwMatrixN& data, double value, hwMatrixNI* count,
                       int& dim)
{
    const std::vector<int>& dims = data.Dimensions();
    int numDim = static_cast<int> (dims.size());

    if (dim == -1)
    {
        // first non-singleton
        for (int i = 0; i < numDim; ++i)
        {
            if (dims[i] != 1)
            {
                dim = i;
                break;
            }
        }
    }
    else if (dim < 0 || dim > numDim - 1)
    {
        return;
    }

    if (count)
    {
        std::vector<int> countDims = dims;
        countDims[dim] = 1;
        count->Dimension(countDims, hwMatrixNI::REAL);
    }

    int numVecs = data.Size() / dims[dim];
    int stride = data.Stride(dim);
    std::vector<int> matrixIndex(numDim);

    for (int i = 0; i < numVecs; ++i)
    {
        // set the matrix indices to the first index in each slice
        int start = data.Index(matrixIndex);
        int total = dims[dim];  // total non-NaN elements

        if (data.IsReal())
        {
            double* real = data.GetRealData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                // replace and count NaN values
                if (IsNaN_T(*real))
                {
                    *real = value;
                    --total;
                }

                real += stride;
            }
        }
        else
        {
            hwComplex* cplx = data.GetComplexData() + start;

            for (int j = 0; j < dims[dim]; ++j)
            {
                // replace and count NaN values
                if (IsNaN_T(*cplx))
                {
                    *cplx = value;
                    --total;
                }

                cplx += stride;
            }
        }

        if (count)
            (*count)(matrixIndex) = total;

        // advance slice indices
        for (int j = 0; j < numDim; ++j)
        {
            if (j == dim)
                continue;

            // increment index j if possible
            if (matrixIndex[j] < static_cast<int> (dims[j]) - 1)
            {
                ++matrixIndex[j];
                break;
            }

            // index j is maxed out, so reset and continue to j+1
            matrixIndex[j] = 0;
        }
    }
}
//------------------------------------------------------------------------------
// Computes min values, excluding NaN [nanmin]
//------------------------------------------------------------------------------
bool OmlNanMin(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs)
{
    return oml_min(eval, inputs, outputs);
}
//------------------------------------------------------------------------------
// Computes max values, excluding NaN [nanmax]
//------------------------------------------------------------------------------
bool OmlNanMax(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs)
{
    return oml_max(eval, inputs, outputs);
}
//------------------------------------------------------------------------------
// Computes sum values, excluding NaN [nansum]
//------------------------------------------------------------------------------
bool OmlNanSum(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    int dim = -1;

    if (nargin == 2)
    {
        if (inputs[1].IsPositiveInteger())
            dim = static_cast<int> (inputs[1].Scalar()) - 1;
    }

    std::vector<Currency> inputs2;
    Currency copy = inputs[0].GetWritableCurrency();

    if (inputs[0].IsMatrix())
    {
        hwMatrix* m = copy.GetWritableMatrix();
        NanAssignmentUtil(*m, 0.0, nullptr, dim);
        inputs2.push_back(copy);

        if (nargin == 2)
        {
            inputs2.push_back(inputs[1]);
        }

        oml_sum(eval, inputs2, outputs);
    }
    else if (inputs[0].IsNDMatrix())
    {
        hwMatrixN* m = copy.GetWritableMatrixN();
        NanAssignmentUtil(*m, 0.0, nullptr, dim);
        inputs2.push_back(copy);

        if (nargin == 2)
        {
            inputs2.push_back(inputs[1]);
        }

        oml_sum(eval, inputs2, outputs);
    }
    else if (copy.IsScalar())
    {
        if (IsNaN_T(copy.Scalar()))
        {
            outputs.push_back(0.0);
        }
        else
        {
            outputs.push_back(copy);
        }
    }
    else if (copy.IsComplex())
    {
        if (IsNaN_T(copy.Complex()))
        {
            outputs.push_back(0.0);
        }
        else
        {
            outputs.push_back(copy);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, 1, OML_VAR_DATA);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes mean values, excluding NaN [nanmean]
//------------------------------------------------------------------------------
bool OmlNanMean(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    int dim = -1;

    if (nargin == 2)
    {
        if (inputs[1].IsPositiveInteger())
            dim = static_cast<int> (inputs[1].Scalar()) - 1;
    }

    std::vector<Currency> inputs2;
    Currency copy = inputs[0].GetWritableCurrency();

    if (copy.IsMatrix())
    {
        hwMatrix* m = copy.GetWritableMatrix();

        if (!m->IsReal())
        {
            throw OML_Error(OML_ERR_REAL, 1);
        }

        hwMatrixI count;
        NanAssignmentUtil(*m, 0.0, &count, dim);
        inputs2.push_back(copy);

        if (nargin == 2)
        {
            inputs2.push_back(inputs[1]);
        }

        oml_sum(eval, inputs2, outputs);

        if (outputs[0].IsMatrix())
        {
            hwMatrix* mean = outputs[0].GetWritableMatrix();
            int size = mean->Size();

            for (int i = 0; i < size; ++i)
            {
                if (count(i))
                    (*mean)(i) /= static_cast<int> (count(i));
                else
                    (*mean)(i) = std::numeric_limits<double>::quiet_NaN();
            }
        }
        else // outputs[0].IsScalar()
        {
            double mean = outputs[0].Scalar();
            
            if (count(0))
                mean /= count(0);
            else
                mean = std::numeric_limits<double>::quiet_NaN();
            
            outputs.clear();
            outputs.push_back(mean);
        }
    }
    else if (copy.IsNDMatrix())
    {
        hwMatrixN* m = copy.GetWritableMatrixN();

        if (!m->IsReal())
        {
            throw OML_Error(OML_ERR_REAL, 1);
        }

        hwMatrixNI count;
        NanAssignmentUtil(*m, 0.0, &count, dim);
        inputs2.push_back(copy);

        if (nargin == 2)
        {
            inputs2.push_back(inputs[1]);
        }

        oml_sum(eval, inputs2, outputs);

        if (outputs[0].IsNDMatrix())
        {
            hwMatrixN* mean = outputs[0].GetWritableMatrixN();
            int size = mean->Size();

            for (int i = 0; i < size; ++i)
            {
                if (count(i))
                    (*mean)(i) /= static_cast<int> (count(i));
                else
                    (*mean)(i) = std::numeric_limits<double>::quiet_NaN();
            }
        }
        else if (outputs[0].IsMatrix())
        {
            // this happens if dim = 3
            hwMatrix* mean = outputs[0].GetWritableMatrix();
            int size = mean->Size();

            for (int i = 0; i < size; ++i)
            {
                if (count(i))
                    (*mean)(i) /= static_cast<int> (count(i));
                else
                    (*mean)(i) = std::numeric_limits<double>::quiet_NaN();
            }
        }
        else // outputs[0].IsScalar()
        {
            double mean = outputs[0].Scalar() / count(0);

            if (count(0))
                mean /= count(0);
            else
                mean = std::numeric_limits<double>::quiet_NaN();

            outputs.clear();
            outputs.push_back(mean);
        }
    }
    else if (copy.IsScalar())
    {
        if (IsNaN_T(copy.Scalar()))
        {
            outputs.push_back(0.0);
        }
        else
        {
            outputs.push_back(inputs[0]);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, 1, OML_VAR_DATA);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes standard deviation values, excluding NaN [nanstd]
//------------------------------------------------------------------------------
bool OmlNanStd(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs)
{
    OmlNanVar(eval, inputs, outputs);

    if (outputs[0].IsScalar())
    {
        double std = sqrt(outputs[0].Scalar());
        outputs.clear();
        outputs.push_back(std);
    }
    else if (outputs[0].IsMatrix())
    {
        hwMatrix* var = outputs[0].GetWritableMatrix();
        int size = var->Size();

        for (int i = 0; i < size; ++i)
            (*var)(i) = sqrt((*var)(i));
    }
    else // outputs[0].IsNDMatrix()
    {
        hwMatrixN* var = outputs[0].GetWritableMatrixN();
        int size = var->Size();

        for (int i = 0; i < size; ++i)
            (*var)(i) = sqrt((*var)(i));
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes variance values, excluding NaN [nanvar]
//------------------------------------------------------------------------------
bool OmlNanVar(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    bool sampleStat;    // false = population statistic

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

    int dim = -1;

    if (nargin == 3)
    {
        if (inputs[2].IsPositiveInteger())
            dim = static_cast<int> (inputs[2].Scalar()) - 1;
    }

    Currency copy = inputs[0].GetWritableCurrency();

    if (copy.IsMatrix())
    {
        hwMatrix* data = copy.GetWritableMatrix();

        if (!data->IsReal())
        {
            throw OML_Error(OML_ERR_REAL, 1);
        }

        hwMatrixI count;
        NanAssignmentUtil(*data, 0.0, &count, dim);

        if (dim != 0 && dim != 1)
        {
            outputs.push_back(copy);
            return true;
        }

        int m = data->M();
        int n = data->N();
        hwMatrix* var = EvaluatorInterface::allocateMatrix();

        if (dim == 0)
            var->Dimension(1, n, hwMatrix::REAL);
        else
            var->Dimension(m, 1, hwMatrix::REAL);

        std::vector<int> dims = { m, n };
        int numVecs = (dim == 0) ? n : m;
        int stride = (dim == 0) ? 1 : m;
        int row = 0;
        int col = 0;

        for (int i = 0; i < numVecs; ++i)
        {
            // set the matrix indices to the first index in each slice
            int start = col * m + row;

            // shift data and compensate for replaced NaN values
            double* real = data->GetRealData() + start;
            int m = count(i);
            int n = dims[dim];
            double data_zero;

            for (int j = 0; j < n; ++j)
            {
                if (!IsNaN_T(*(real + j)))
                {
                    data_zero = *(real + j);
                    break;
                }
            }

            double sum   = 0.0;
            double sumSq = 0.0;
            double value;

            for (int j = 0; j < n; ++j)
            {
                value  = *real - data_zero;     // shift to reduce overflow risk
                sum   += value;
                sumSq += value * value;
                real  += stride;
            }

            sum   += (n - m) * data_zero;
            sumSq -= (n - m) * data_zero * data_zero;

            if (sampleStat)
            {
                if (m != 1)
                {
                    (*var)(i) = ((double)m * sumSq - sum * sum) / ((double)m * (double)(m - 1));
                }
                else // (m == 1)
                {
                    (*var)(i) = 0.0;
                }
            }
            else // population variance
            {
                (*var)(i) = ((double)m * sumSq - sum * sum) / ((double)(m * m));
            }

            // advance slice indices
            if (dim == 0)
                ++col;
            else
                ++row;
        }

        outputs.clear();
        outputs.push_back(var);
    }
    else if (copy.IsNDMatrix())
    {
        hwMatrixN* data = copy.GetWritableMatrixN();

        if (!data->IsReal())
        {
            throw OML_Error(OML_ERR_REAL, 1);
        }

        hwMatrixNI count;
        NanAssignmentUtil(*data, 0.0, &count, dim);

        const std::vector<int>& dims = data->Dimensions();
        int numDim = static_cast<int> (dims.size());

        if (dim > numDim - 1)
        {
            outputs.push_back(copy);
            return true;
        }

        const std::vector<int>& varDims = count.Dimensions();
        hwMatrixN* var = new hwMatrixN(varDims, hwMatrixN::REAL);

        int numVecs = data->Size() / dims[dim];
        int stride = data->Stride(dim);
        std::vector<int> matrixIndex(numDim);

        for (int i = 0; i < numVecs; ++i)
        {
            // set the matrix indices to the first index in each slice
            int start = data->Index(matrixIndex);

            // shift data and compensate for replaced NaN values
            double* real = data->GetRealData() + start;
            int m = count(i);
            int n = dims[dim];
            double data_zero;

            for (int j = 0; j < n; ++j)
            {
                if (!IsNaN_T(*(real + j)))
                {
                    data_zero = *(real + j);
                    break;
                }
            }

            double sum   = 0.0;
            double sumSq = 0.0;
            double value;

            for (int j = 0; j < n; ++j)
            {
                value  = *real - data_zero;        // shift to reduce overflow risk
                sum   += value;
                sumSq += value * value;
                real  += stride;
            }

            sum   += (n - m) * data_zero;
            sumSq -= (n - m) * data_zero * data_zero;

            if (sampleStat)
            {
                if (m != 1)
                {
                    (*var)(i) = ((double)m * sumSq - sum * sum) / ((double)m * (double)(m - 1));
                }
                else // (m == 1)
                {
                    (*var)(i) = 0.0;
                }
            }
            else // population variance
            {
                (*var)(i) = ((double)m * sumSq - sum * sum) / ((double)(m * m));
            }

            // advance slice indices
            for (int j = 0; j < numDim; ++j)
            {
                if (j == dim)
                    continue;

                // increment index j if possible
                if (matrixIndex[j] < static_cast<int> (dims[j]) - 1)
                {
                    ++matrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                matrixIndex[j] = 0;
            }
        }

        outputs.clear();
        outputs.push_back(var);
    }
    else if (copy.IsScalar())
    {
        if (IsNaN_T(copy.Scalar()))
        {
            outputs.push_back(0.0);
        }
        else
        {
            outputs.push_back(copy);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, 1, OML_VAR_DATA);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes median values, excluding NaN [nanmedian]
//------------------------------------------------------------------------------
bool OmlNanMedian(EvaluatorInterface           eval,
                  const std::vector<Currency>& inputs,
                  std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    int dim = -1;

    if (nargin == 2)
    {
        if (inputs[1].IsPositiveInteger())
            dim = static_cast<int> (inputs[1].Scalar()) - 1;
    }

    std::vector<Currency> inputs2;
    Currency copy = inputs[0].GetWritableCurrency();

    if (copy.IsMatrix())
    {
        hwMatrix* data = copy.GetWritableMatrix();

        if (!data->IsReal())
        {
            throw OML_Error(OML_ERR_REAL, 1);
        }

        hwMatrixI count;
        NanAssignmentUtil(*data, std::numeric_limits<double>::infinity(), &count, dim);

        if (dim != 0 && dim != 1)
        {
            outputs.push_back(inputs[0]);
            return true;
        }

        inputs2.push_back(copy);

        if (nargin == 2)
        {
            inputs2.push_back(inputs[1]);
        }

        oml_sort(eval, inputs2, outputs);

        data = outputs[0].GetWritableMatrix();
        int m = data->M();
        int n = data->N();

        hwMatrix* median = EvaluatorInterface::allocateMatrix();

        if (dim == 0)
            median->Dimension(1, n, hwMatrix::REAL);
        else
            median->Dimension(m, 1, hwMatrix::REAL);

        std::vector<int> dims = { m, n };
        int numVecs = (dim == 0) ? n : m;
        int stride = (dim == 0) ? 1 : m;
        int row = 0;
        int col = 0;

        for (int i = 0; i < numVecs; ++i)
        {
            // set the matrix indices to the first index in each slice
            int start = col * m + row;

            // compute median value
            int cnt = count(i);
            double* real = data->GetRealData() + start;

            if (cnt == 0)
            {
                (*median)(i) = std::numeric_limits<double>::quiet_NaN();
            }
            else if (cnt % 2)
            {
                (*median)(i) = *(real + cnt / 2 * stride);
            }
            else
            {
                (*median)(i)  = *(real + (cnt / 2 - 1) * stride);
                (*median)(i) += *(real + (cnt / 2) * stride);
                (*median)(i) /= 2.0;
            }

            // advance slice indices
            if (dim == 0)
                ++col;
            else
                ++row;
        }

        outputs.clear();
        outputs.push_back(median);
    }
    else if (copy.IsNDMatrix())
    {
        hwMatrixN* data = copy.GetWritableMatrixN();

        if (!data->IsReal())
        {
            throw OML_Error(OML_ERR_REAL, 1);
        }

        hwMatrixNI count;
        NanAssignmentUtil(*data, std::numeric_limits<double>::infinity(), &count, dim);

        const std::vector<int>& dims = data->Dimensions();
        int numDim = static_cast<int> (dims.size());

        if (dim > numDim - 1)
        {
            outputs.push_back(inputs[0]);
            return true;
        }

        inputs2.push_back(copy);

        if (nargin == 2)
        {
            inputs2.push_back(inputs[1]);
        }

        oml_sort(eval, inputs2, outputs);

        data = outputs[0].GetWritableMatrixN();

        std::vector<int> medianDims = count.Dimensions();
        hwMatrixN* median = new hwMatrixN(medianDims, hwMatrixN::REAL);

        int numVecs = data->Size() / dims[dim];
        int stride = data->Stride(dim);
        std::vector<int> matrixIndex(numDim);

        for (int i = 0; i < numVecs; ++i)
        {
            // set the matrix indices to the first index in each slice
            int start = data->Index(matrixIndex);

            // compute median values
            int cnt = count(i);
            double* real = data->GetRealData() + start;

            if (cnt == 0)
            {
                (*median)(i) = std::numeric_limits<double>::quiet_NaN();
            }
            else if (cnt % 2)
            {
                (*median)(i) = *(real + cnt / 2 * stride);
            }
            else
            {
                (*median)(i) = *(real + (cnt / 2 - 1) * stride);
                (*median)(i) += *(real + (cnt / 2) * stride);
                (*median)(i) /= 2.0;
            }

            // advance slice indices
            for (int j = 0; j < numDim; ++j)
            {
                if (j == dim)
                    continue;

                // increment index j if possible
                if (matrixIndex[j] < static_cast<int> (dims[j]) - 1)
                {
                    ++matrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                matrixIndex[j] = 0;
            }
        }

        outputs.clear();
        outputs.push_back(median);
    }
    else if (copy.IsScalar())
    {
        if (IsNaN_T(copy.Scalar()))
        {
            outputs.push_back(0.0);
        }
        else
        {
            outputs.push_back(copy);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, 1, OML_VAR_DATA);
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
// Generates a random permutation vector [randperm]
//------------------------------------------------------------------------------
bool OmlRandperm(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs)
{
    int nargin = eval.GetNarginValue();

    if (nargin != 1 && nargin != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsPositiveInteger() && !inputs[0].IsPositiveInteger64() &&
        (!inputs[0].IsScalar() || inputs[0].Scalar() != 0.0) && !inputs[0].IsEmpty())
        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DATA);

    hwMathStatus status;
    CreateTwister();

    if (inputs[0].IsPositiveInteger() || (inputs[0].IsScalar() && inputs[0].Scalar() == 0.0))
    {
        int max = static_cast<int> (inputs[0].Scalar());
        int numPts = max;

        if (nargin == 2)
        {
            if (!inputs[1].IsInteger())
                throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DATA);

            numPts = static_cast<int> (inputs[1].Scalar());
        }

        hwMatrixI permVecI;
        status = RandPerm(max, numPts, twister, permVecI);

        if (status.IsOk())
        {
            hwMatrix* permVec = EvaluatorInterface::allocateMatrix(1, numPts, true);

            for (int i = 0; i < numPts; ++i)
            {
                (*permVec)(i) = static_cast<double> (permVecI(i));
            }

            outputs.push_back(permVec);
        }
    }
    else if (inputs[0].IsEmpty())
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix(1, 0, true));
    }
    else
    {
        if (!inputs[1].IsPositiveInteger())
        {
            if (inputs[1].IsPositiveInteger64())
                throw OML_Error(hwMathStatus(HW_MATH_ERR_ALLOCFAILED));  // alloc failure
            else
                throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DATA);
        }

        int64_t max = static_cast<int64_t> (inputs[0].Scalar());
        int numPts = static_cast<int> (inputs[1].Scalar());

        hwMatrixI64 permVecI;
        status = RandPerm(max, numPts, twister, permVecI);

        if (status.IsOk())
        {
            hwMatrix* permVec = EvaluatorInterface::allocateMatrix(1, numPts, true);

            for (int i = 0; i < numPts; ++i)
            {
                (*permVec)(i) = static_cast<double> (permVecI(i));
            }

            outputs.push_back(permVec);
        }
    }

    if (!status.IsOk())
    {
        if (status.GetArg1() == 1)
        {
            status.SetArg1(status.GetArg2());
            status.SetArg2(-1);
        }

        if (status.GetArg1() == 2)
            status.SetArg1(1);
        else if (status.GetArg1() == 3)
            status.SetArg1(2);
        else
            status.ResetArgs();
    }

    BuiltInFuncsUtils::CheckMathStatus(eval, status);

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
                                                        true);
    
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
                hwMatrix* out = EvaluatorInterface::allocateMatrix(outi.M(), outi.N(), true);

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
        transform.reset(EvaluatorInterface::allocateMatrix(1, 2, true));

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
// Compute binomial combinations  [nchoosek]
//------------------------------------------------------------------------------
bool OmlNchooseK(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[1].IsInteger())
        throw OML_Error(OML_ERR_NNINTVECTOR, 1, OML_VAR_TYPE);

    int k = static_cast<int> (inputs[1].Scalar());

    if (inputs[0].IsPositiveInteger())
    {
        int n = static_cast<int> (inputs[0].Scalar());
        double nCk;

        BuiltInFuncsUtils::CheckMathStatus(eval, NchooseK(n, k, nCk));
        outputs.push_back(nCk);
    }
    else if (inputs[0].IsMatrix())
    {
        const hwMatrix* N = inputs[0].Matrix();
        std::unique_ptr<hwMatrix> combos(EvaluatorInterface::allocateMatrix());

        BuiltInFuncsUtils::CheckMathStatus(eval, NchooseK(*N, k, *combos));
        outputs.push_back(combos.release());
    }
    else
    {
        throw OML_Error(OML_ERR_POSINTEGER_VEC, 1, OML_VAR_TYPE);
    }

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
        if (!inputs[i-1].IsScalar())
            throw OML_Error(OML_ERR_NATURALNUM, i);
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

    if (!nargin)
    {
        m = 1;
        n = 1;
    }
    else if (firstDimArg < 1 || firstDimArg > nargin)
    {
        if (inputs[0].IsScalar())
        {
            if (inputs[1].IsScalar())
            {
                m = 1;
                n = 1;
            }
            else if (inputs[1].IsMatrix())
            {
                m = inputs[1].Matrix()->M();
                n = inputs[1].Matrix()->N();
            }
        }
        else if (inputs[1].IsScalar())
        {
            if (inputs[0].IsMatrix())
            {
                m = inputs[0].Matrix()->M();
                n = inputs[0].Matrix()->N();
            }
        }
        else
        {
            m = inputs[0].Matrix()->M();
            n = inputs[0].Matrix()->N();

            if (inputs[1].Matrix()->M() != m)
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

            if (inputs[1].Matrix()->N() != n)
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
        }
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

            if (!isint(realval(dims, 0)))
                throw OML_Error(OML_ERR_NATURALNUM, firstDimArg, OML_VAR_DIMS);

            if (!isint(realval(dims, 1)))
                throw OML_Error(OML_ERR_NATURALNUM, firstDimArg, OML_VAR_DIMS);

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
        if (!inputs[firstDimArg - 1].IsScalar())
        {
            throw OML_Error(OML_ERR_NATURALNUM, firstDimArg, OML_VAR_DIMS);
        }

        if (!inputs[firstDimArg].IsScalar())
        {
            throw OML_Error(OML_ERR_NATURALNUM, firstDimArg+1, OML_VAR_DIMS);
        }

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
        arg = EvaluatorInterface::allocateMatrix(m, n, true);
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
