/**
* @file OptimizationFuncs.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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

#include "hwOptimizationFuncs.h"

//------------------------------------------------------------------------------
// Total sum of squares of squares
//------------------------------------------------------------------------------
static
hwMathStatus TotalSumOfSquares(const hwMatrix& data,
                               double&         sst)
{
    // This is a copy from statistics/_private/StatUtilFuncs.cxx
    // \todo: Eliminate when restructuring
    if (data.IsEmpty())
    {
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);
    }
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int n = data.Size();

    double sum   = 0.0;
    double sumSq = 0.0;
    double value;

    for (int i = 0; i < n; i++)
    {
        value  = data(i);
        sum   += value;
        sumSq += value * value;
    }

    sst = (n > 1) ? 
          sumSq - sum * sum / (double) n :
          0.0;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Nonlinear curve fit which works on surfaces and curves
//------------------------------------------------------------------------------
hwMathStatus NLCurveFit(const LSqFitFunc pRespFunc,
                        const LSqFitFunc pJacFunc,
                        hwMatrix&        P,
                        const hwMatrix&  X,
                        const hwMatrix&  y,
                        const hwMatrix*  lowerBound,
                        const hwMatrix*  upperBound,
                        int&             maxIter,
                        int&             maxFuncEval,
                        hwMatrix*        stats,
                        hwMatrix*        yEst,
                        double           tolf,
                        double           tolx,
                        hwMatrix*        objHist,
                        hwMatrix*        designHist,
                        const hwMatrix*  userData)
{
    hwGenericFuncFitter fittedFunc(pRespFunc, pJacFunc, P, X, y,
                                   lowerBound, upperBound, maxIter,
                                   maxFuncEval, tolf, tolx, userData);

    fittedFunc.RequestDesignHist(designHist);
    fittedFunc.RequestObjectiveHist(objHist);

    hwMathStatus status = fittedFunc.Compute();

    if (!status.IsOk())
    {
        switch (status.GetArg1())
        {
            case 10: status.SetArg1(12); break;
            case 11: status.SetArg1(13); break;
            default: break;
        }

        if (status == HW_MATH_ERR_UNDERDETSYS_E)
        {
            status.SetMsgCode(HW_MATH_ERR_UNDERDETSYS_P);
            return status;
        }

        if (!status.IsWarning() && !status.IsInfoMsg())
        {
            return status;
        }
    }

    fittedFunc.GetParams(P);
    maxIter = fittedFunc.Iterations();
    maxFuncEval -= fittedFunc.FunctionEvals();

    // compute statistics and fitted points
    int n = y.Size();

    if (stats)
    {
        double sse = 0.0;
        double sst;
        double value;

        hwMathStatus status2 = stats->Dimension(2, hwMatrix::REAL);
        if (!status2.IsOk())
        {
            if (status2.GetArg1() == 0)
            {
                status2.SetArg1(10);
            }
            else
            {
                status2.ResetArgs();
            }
            return status2;
        }

        if (yEst)
        {
            status2 = yEst->Dimension(y.M(), y.N(), hwMatrix::REAL);
            if (!status2.IsOk())
            {
                if (status2.GetArg1() == 0)
                {
                    status2.SetArg1(11);
                }
                else
                {
                    status2.ResetArgs();
                }
                return status2;
            }

            status2 = pRespFunc(P, X, userData, *yEst);
            if (!status2.IsOk())
            {
                return status2;
            }

            for (int i = 0; i < n; i++)
            {
                value = y(i) - (*yEst)(i);
                sse  += value * value;
            }
        }
        else
        {
            hwMatrix temp_y(y.M(), y.N(), hwMatrix::REAL);

            status2 = pRespFunc(P, X, userData, temp_y);
            if (!status2.IsOk())
            {
                return status2;
            }

            for (int i = 0; i < n; i++)
            {
                value = y(i) - temp_y(i);
                sse  += value * value;
            }
        }

        status2 = TotalSumOfSquares(y, sst);
        if (!status2.IsOk())
        {
            return status2;
        }

        (*stats)(0) = sse / n;              // mse

        if (IsZero(sst, 1.0e-12))
        {
            (*stats)(1) = 1.0;
        }
        else
        {
            value = 1.0 - sse / sst;
            (*stats)(1) = (value > 0.0) ?
                          sqrt(value) : // rho
                          0.0;
        }
    }
    else if (yEst)
    {
        hwMathStatus status2 = yEst->Dimension(y.M(), y.N(), hwMatrix::REAL);
        if (!status2.IsOk())
        {
            if (status2.GetArg1() == 0)
            {
                status2.SetArg1(11);
            }
            else
            {
                status2.ResetArgs();
            }
            return status2;
        }

        status2 = pRespFunc(P, X, userData, *yEst);
        if (!status2.IsOk())
        {
            return status2;
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Solve a nonlinear system of equations
//------------------------------------------------------------------------------
hwMathStatus NLSolve(const LSqFitFunc pRespFunc,
                     const LSqFitFunc pJacFunc,
                     hwMatrix&        P,
                     double&          minVal,
                     int&             maxIter,
                     int&             maxFuncEval,
                     int              numEqns_,
                     double           tolf,
                     double           tolx,
                     hwMatrix*        objHist,
                     hwMatrix*        designHist,
                     const hwMatrix*  userData)
{
    int numEqns = (numEqns_ == -1) ?
                  P.Size() :    // solve system
                  numEqns_ ;    // minimize over determined system

    hwGenericFuncFitter fittedFunc(pRespFunc, pJacFunc, P, numEqns, maxIter,
                                   maxFuncEval, tolf, tolx, userData);

    fittedFunc.RequestDesignHist(designHist);
    fittedFunc.RequestObjectiveHist(objHist);

    hwMathStatus status = fittedFunc.Compute();
    if (!status.IsOk())
    {
        switch (status.GetArg1())
        {
            case 4: status.SetArg1(7); break;
            case 7: status.SetArg1(8); break;
            case 8: status.SetArg1(9); break;
            default: break;
        }

        if (status.GetArg2() == 4)
        {
            status.SetArg2(7);
        }

        if (!status.IsWarning() && !status.IsInfoMsg())
        {
            return status;
        }
    }

    fittedFunc.GetParams(P);
    minVal = fittedFunc.ObjFuncVal();
    maxIter = fittedFunc.Iterations();
    maxFuncEval -= fittedFunc.FunctionEvals();

    if (!status.IsWarning())
    {
        if (tolf)
        {
            if (minVal > 10.0 * tolf)
            {
                status(HW_MATH_WARN_NOSOLUTION);
            }
        }
        else
        {
            if (minVal > 10.0 * tolx)
            {
                status(HW_MATH_WARN_NOSOLUTION);
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Utility to update Brent history
//------------------------------------------------------------------------------
static void BrentHistory(hwMatrix* designHist,
                         hwMatrix* objHist,
                         double    a,
                         double    b,
                         double    fa,
                         double    fb,
                         int       iter)
{
    hwMathStatus status;

    if (objHist)
    {
        if (objHist->IsEmpty())
            status = objHist->Dimension(iter, 2, hwMatrix::REAL);
        else
            status = objHist->Resize(iter, 2);

        (*objHist)(iter-1, 0) = fa;
        (*objHist)(iter-1, 1) = fb;
    }

    if (designHist)
    {
        if (designHist->IsEmpty())
            status = designHist->Dimension(iter, 2, hwMatrix::REAL);
        else
            status = designHist->Resize(iter, 2);

        (*designHist)(iter-1, 0) = a;
        (*designHist)(iter-1, 1) = b;
    }
}
//------------------------------------------------------------------------------
// Find the minimum of a function without derivatives using Brent's method
//------------------------------------------------------------------------------
hwMathStatus Brent(const BrentFunc f,
                   double&         a,
                   double&         b,
                   double&         fa,
                   double&         fb,
                   double&         min,
                   double&         fmin,
                   int&            maxIter,
                   int&            maxFuncEval, 
                   double          tol,
                   hwMatrix*       objHist,
                   hwMatrix*       designHist)
{
  // An approximation x to the point where f attains a minimum on the interval
  // (a,b)  is determined.
  // The method used is a combination of golden section searchand successive
  // parabolic interpolation. Convergence is never much slower than that for a  
  // fibonacci search. If f has a continuous second derivative which is positive 
  // at the minimum (which is not at a or b), then convergence is superlinear, 
  // and usually of the order of about  1.324. The function f is never evaluated 
  // at two points closer together than eps*abs(fmin) + (tol/3), where eps is  
  // approximately the square root of the relative machine precision. If f is a 
  // unimodal function and the computed values of f are always unimodal when
  // separated by at least eps*abs(x) + (tol/3), then fmin approximates the 
  // abcissa of the global minimum of f on the interval a, b with an error less 
  // than 3*eps*abs(fmin) + tol. If f is not unimodal, then fmin may approximate 
  // a local, but perhaps non-global, minimum to the same accuracy. This function
  // is a slightly modified version of the algol 60 procedure localmin given in 
  // Richard Brent, algorithms for minimization without derivatives, 
  // Prentice-Hall, Inc. (1973).

    // c is the squared inverse of the golden ratio
    double c = 0.5 * (3.0 - sqrt(5.0));

    // eps is approximately the square root of the relative machine precision
    double eps = sqrt(MACHEP);

    // Initialization
    double v = a + c * (b - a);
    double w = v;
    double x = v;
    double e = 0.0;
    double fx;

    hwMathStatus status = f(x, fx);

    if (objHist || designHist)
    {
        status = f(a, fa);
        status = f(b, fb);
    }

    double fv       = fx;
    double fw       = fx;
    double d;
    double xm;
    double p;
    double q;
    double r;
    double tol1;
    double tol2;
    double u;
    double fu;
    bool   parabolicStep;
    int numFunEvals = 1;
    int numIters    = 0;

    // main loop starts here
    while (numIters < maxIter && numFunEvals < maxFuncEval)
    {
        if (objHist || designHist)
            BrentHistory(designHist, objHist, a, b, fa, fb, numIters+1);

        xm   = (a + b) / 2.0;
        tol1 = eps * fabs(x) + tol / 3.0;
        tol2 = 2.0 * tol1;

        // check stopping criterion
        if (fabs(x - xm) + 0.5 * (b - a) <= tol2)
        {
            break;
        }

        if (fabs(e) > tol1)
        {
            // fit parabola
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);

            if (q > 0.0)
            {
                p = -p;
            }

            q = fabs(q);
            r = e;
            e = d;

            if ( (fabs(p) >= fabs(0.5 * q * r)) ||
                 (p <= q * (a - x))             ||
                 (p >= q * (b - x)))
            {
                parabolicStep = false;
            }
            else
            {
                parabolicStep = true;
            }
        }
        else
        {
            parabolicStep = false;
        }

        if (parabolicStep)
        {
            // parabolic interpolation step
            d = p / q;
            u = x + d;

            // f must not be evaluated too close to ax or bx
            if (u - a < tol2)
            {
                d = sign(tol1, xm - x);
            }

            if (b - u < tol2)
            {
                d = sign(tol1, xm - x);
            }
        }
        else
        {
            // golden-section step
            e = (x >= xm) ? a - x : b - x;
            d = c * e;
        }

        // f must not be evaluated too close to x
        if (fabs(d) < tol1)
        {
            u = x + sign(tol1, d);
        }
        else
        {
            u = x + d;
        }

        status = f(u, fu);
        numFunEvals++;

        // update a, b, v, w, and x
        if (fu < fx)
        {
            if (u < x)
            {
                b  = x;
                fb = fx;
            }
            else
            {
                a  = x;
                fa = fx;
            }

            v  = w;
            w  = x;
            x  = u;
            fv = fw;
            fw = fx;
            fx = fu;
        }
        else
        {
            if (u < x)
            {
                a  = u;
                fa = fu;
            }
            else
            {
                b  = u;
                fb = fu;
            }

            if (fu <= fw || w == x)
            {
                v  = w;
                w  = u;
                fv = fw;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v  = u;
                fv = fu;
            }
        }

        numIters++;
    }

    min  = x;
    fmin = fx;

    if (numIters == maxIter)
    {
        status(HW_MATH_WARN_MAXITERATE);
    }
    else if (numFunEvals >= maxFuncEval)
    {
        status(HW_MATH_WARN_MAXFUNCEVAL);
    }

    maxIter     = numIters;
    maxFuncEval = numFunEvals;

    return status;
}
//------------------------------------------------------------------------------
// Find zero of a function without derivatives using the TOMS 748 algorithm
//------------------------------------------------------------------------------
hwMathStatus F_zero(const TOMS748Func f, 
                    double&           a, 
                    double&           b, 
                    double&           fa, 
                    double&           fb,
                    double&           root, 
                    double&           froot, 
                    int&              maxIter,
                    int&              maxFuncEval,
                    double            tol,
                    hwMatrix*         objHist,
                    hwMatrix*         designHist)
{
    return TOMS748(f, a, b, fa, fb, root, froot, maxIter, maxFuncEval, tol,
                   objHist, designHist);
}
//------------------------------------------------------------------------------
// Utility to update NelderMead history
//------------------------------------------------------------------------------
static void NelderMeadHistory(hwMatrix*       designHist,
                              hwMatrix*       objHist,
                              const hwMatrix& x,
                              double          fx,
                              int             iter)
{
    int size = x.Size();
    hwMathStatus status;

    if (objHist)
    {
        if (objHist->IsEmpty())
            status = objHist->Dimension(iter, 1, hwMatrix::REAL);
        else
            status = objHist->Resize(iter, 1);

        (*objHist)(iter-1, 0) = fx;
    }

    if (designHist)
    {
        if (designHist->IsEmpty())
            status = designHist->Dimension(iter, size, hwMatrix::REAL);
        else
            status = designHist->Resize(iter, size);

        for (int i = 0; i < size; ++i)
            (*designHist)(iter-1, i) = x(i);
    }
}
//------------------------------------------------------------------------------
// Finds minimum of a function without derivatives using Nelder-Mead simplex method
//------------------------------------------------------------------------------
hwMathStatus NelderMead(const NelderMeadFunc f, 
                        hwMatrix&            optPoint, 
                        double&              minVal,
                        int&                 maxIter,
                        int&                 maxFuncEval,
                        double               tolx,
                        hwMatrix*            objHist,
                        hwMatrix*            designHist)
{
    hwMathStatus status;
    bool needTranspose = false;

    if (!optPoint.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!optPoint.IsVector())
    {
        return status(HW_MATH_ERR_VECTOR, 2);
    }

    if (optPoint.N() > 1)
    {
        optPoint.Transpose();
        needTranspose = true;
    }

    int        dim = optPoint.Size();
    hwMatrix   fx(dim + 1, hwMatrix::REAL);
    hwMatrix** vertex = new hwMatrix*[dim+1];   // the simplex

    // initialize the simplex
    for (int i = 0; i < dim+1; ++i)
    {
        vertex[i] = new hwMatrix(optPoint);
    }

    for (int i = 0; i < dim; ++i)
    {
        (*vertex[i])(i) += 1.0;
    }

    // evaulate the function at each vertex
    int numFunEvals = dim + 1;

    for (int i = 0; i < dim + 1; i++)
    {
        status = f(*vertex[i], fx(i));

        if (!status.IsOk())
        {
            numFunEvals = 0;
            break;
        }
    }

    // search
    double* min  = nullptr;
    double* max  = nullptr;
    double* max2 = nullptr;
    double next;
    double norm1;
    double norm2;
    double scale;

    hwMatrix midpoint(dim, hwMatrix::REAL);
    hwMatrix line(dim,     hwMatrix::REAL);
    hwMatrix nextv(dim,    hwMatrix::REAL);

    hwMatrix* minv        = nullptr;
    hwMatrix* maxv        = nullptr;
    int       numIters    = 0;

    while (numIters < maxIter && numFunEvals < maxFuncEval)
    {
        // find simplex extremes
        if (fx(0) > fx(1))
        {
            max = &fx(0);
            min = &fx(1);
            max2 = min;
            maxv = vertex[0];
            minv = vertex[1];
        }
        else
        {
            max = &fx(1);
            min = &fx(0);
            max2 = min;
            maxv = vertex[1];
            minv = vertex[0];
        }

        for (int i = 2; i < dim+1; ++i)
        {
            if (fx(i) <= *min)
            {
                min = &fx(i);
                minv = vertex[i];
            }
            else if (fx(i) > *max)
            {
                max2 = max;
                max = &fx(i);
                maxv = vertex[i];
            }
            else if (fx(i) > *max2)
            {
                max2 = &fx(i);
            }
        }

        // check convergence
        status = (*maxv-*minv).L2Norm(norm1);
        status = maxv->L2Norm(norm2);

        if (objHist || designHist)
            NelderMeadHistory(designHist, objHist, *minv, *min, numIters+1);

        if (norm1 < norm2 * tolx + 0.1 * tolx)
        {
            break;
        }

        // set simplex direction
        midpoint.SetElements(0.0);

        for (int i = 0; i < dim + 1; i++)
        {
            if (vertex[i] != maxv)
            {
                midpoint += *vertex[i];
            }
        }

        midpoint /= dim;
        line = (*maxv) - midpoint;

        // update simplex
        scale = -1.0;   // reflect
        nextv = midpoint + scale * line;
        ++numFunEvals;
        status = f(nextv, next);

        if (!status.IsOk())
        {
            break;
        }

        if (next < *max)
        {
            *max    = next;
            (*maxv) = nextv; 
        }

        if (*max < *min)
        {
            scale = -2.0;   // extend
            nextv = midpoint + scale * line;
            ++numFunEvals;

            status = f(nextv, next);
            if (!status.IsOk())
            {
                break;
            }

            if (next < *max)
            {
                *max = next;
                (*maxv) = nextv; 
            }
        }
        else if (*max >= *max2)
        {
            scale = 0.5;    // contract
            nextv = midpoint + scale * line;
            ++numFunEvals;

            status = f(nextv, next);
            if (!status.IsOk())
            {
                break;
            }

            if (next < *max)
            {
                *max = next;
                (*maxv) = nextv; 
            }
            else
            {
                // shrink
                for (int i = 0; i < dim+1; ++i)
                {
                    if (vertex[i] != minv)
                    {
                        *vertex[i] = scale * (*vertex[i] + *minv);
                        ++numFunEvals;

                        status = f(*vertex[i], fx(i));
                        if (!status.IsOk())
                        {
                            break;
                        }
                    }
                }
            }
        }

        ++numIters;
    }

    if (status.IsOk())
    {
        if (numIters >= maxIter)
        {
            status.ResetArgs();
            status(HW_MATH_WARN_MAXITERATE);
        }
        else if (numFunEvals >= maxFuncEval)
        {
            status.ResetArgs();
            status(HW_MATH_WARN_MAXFUNCEVAL);
        }

        minVal   = *min;
        optPoint = *minv;
    }

    maxIter     = numIters;
    maxFuncEval = numFunEvals;

    if (needTranspose)
        optPoint.Transpose();

    for (int i = 0; i < dim+1; ++i)
    {
        delete vertex[i];
    }

    delete [] vertex;
    return status;
}
//------------------------------------------------------------------------------
// Find the unconstrained minimum of a multivariate function
//------------------------------------------------------------------------------
hwMathStatus FMinUncon(const UnConMinObjFunc  pRespFunc, 
                       const UnConMinGradFunc pGradFunc,
                       hwMatrix&              P, 
                       double&                minVal, 
                       int&                   maxIter,
                       int&                   maxFuncEval,
                       double                 tolf,
                       double                 tolx, 
                       hwMatrix*              objHist, 
                       hwMatrix*              designHist,
                       const hwMatrix*        userData)
{
    hwUnconstrMinimizer minimizer(pRespFunc, pGradFunc, P, maxIter, maxFuncEval,
                                  tolf, tolx, userData);

    minimizer.RequestDesignHist(designHist);
    minimizer.RequestObjectiveHist(objHist);

    hwMathStatus status = minimizer.Compute();

    if (!status.IsOk())
    {
        switch (status.GetArg1())
        {
            case 4: status.SetArg1(5); break;
            case 5: status.SetArg1(6); break;
            case 6: status.SetArg1(7); break;
            case 7: status.SetArg1(8); break;
            default: break;
        }

        if (!status.IsWarning() && !status.IsInfoMsg())
        {
            return status;
        }
    }

    minimizer.GetParams(P);
    minVal = minimizer.ObjFuncVal();
    maxIter = minimizer.Iterations();
    maxFuncEval -= minimizer.FunctionEvals();

    return status;
}
