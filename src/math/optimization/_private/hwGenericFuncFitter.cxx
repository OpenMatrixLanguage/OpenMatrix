/**
* @file hwGenericFuncFitter.cxx
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

#include "hwGenericFuncFitter.h"

//------------------------------------------------------------------------------
// Constructor - Find parameters when fitting equations of the form y=f(x)
//------------------------------------------------------------------------------
hwGenericFuncFitter::hwGenericFuncFitter(const LSqFitFunc pObjFunc,
    const LSqFitFunc pJacFunc,
    const hwMatrix& P,
    const hwMatrix& X_,
    const hwMatrix& y_,
    const hwMatrix* lowerBound,
    const hwMatrix* upperBound,
    int             maxIter,
    int             maxFuncEval,
    double          tolf,
    double          tolx,
    const hwMatrix* userData)
    : hwGaussNewtLSqFit(P, y_.Size(), maxIter, maxFuncEval, tolf, tolx),
      m_lowerBound(lowerBound),
      m_upperBound(upperBound),
      X(nullptr)
{
    if (!m_status.IsOk())
    {
        switch (m_status.GetArg1())
        {
            case 1: m_status.SetArg1(3); break;
            case 2: m_status.SetArg1(5); break;
            case 3: m_status.SetArg1(8); break;
            case 4: m_status.SetArg1(9); break;
            case 5: m_status.SetArg1(10); break;
            case 6: m_status.SetArg1(11); break;
            default: break;
        }

        if (m_status.GetArg2() == 2)
        {
            m_status.SetArg2(5);
        }
        return;
    }

    if (!pObjFunc)
    {
        m_status(HW_MATH_ERR_NULLPOINTER, 1);
        return;
    }

    int numVars = P.Size();

    if (m_lowerBound && m_upperBound)
    {
        if (!m_lowerBound->IsReal())
        {
            m_status(HW_MATH_ERR_COMPLEX, 6);
            return;
        }

        if (!m_lowerBound->IsVector())
        {
            m_status(HW_MATH_ERR_VECTOR, 6);
            return;
        }

        if (!m_upperBound->IsReal())
        {
            m_status(HW_MATH_ERR_COMPLEX, 7);
            return;
        }

        if (!m_upperBound->IsVector())
        {
            m_status(HW_MATH_ERR_VECTOR, 7);
            return;
        }

        if (m_lowerBound->Size() != numVars)
        {
            m_status(HW_MATH_ERR_ARRAYSIZE, 2, 6);
            return;
        }

        if (m_upperBound->Size() != numVars)
        {
            m_status(HW_MATH_ERR_ARRAYSIZE, 2, 7);
            return;
        }

        for (int i = 0; i < numVars; i++)
        {
            if ((*m_lowerBound)(i) > P(i))
                m_status(HW_MATH_ERR_INVALIDINTERVAL, 2, 6);
        }

        for (int i = 0; i < numVars; i++)
        {
            if (P(i) > (*m_upperBound)(i))
                m_status(HW_MATH_ERR_INVALIDINTERVAL, 2, 7);
        }
    }
    else if (m_lowerBound || m_upperBound)
    {
        m_status(HW_MATH_ERR_INVALIDINTERVAL, 6, 7);   // needs to be better
        return;
    }

    if (X_.IsEmpty())
    {
        m_status(HW_MATH_ERR_EMPTYMATRIX, 4);
        return;
    }

    if (!X_.IsReal())
    {
        m_status(HW_MATH_ERR_COMPLEX, 4);
        return;
    }

    if (!y_.IsReal())
    {
        m_status(HW_MATH_ERR_COMPLEX, 5);
        return;
    }

    if (!y_.IsVector())
    {
        m_status(HW_MATH_ERR_VECTOR, 5);
        return;
    }

    if (X_.M() != y_.Size() && X_.N() != y_.Size())
    {
        m_status(HW_MATH_ERR_ARRAYSIZE, 4, 5);
        return;
    }

    m_pObjFunc = pObjFunc;
    m_pJacFunc = pJacFunc;
    X = &X_;
    y = &y_;
    m_userData = userData;

    if ((m_lowerBound && !m_lowerBound->IsEmpty()) || (m_upperBound && !m_upperBound->IsEmpty()))
    {
        m_status = D.Dimension(m_numParams, m_numParams, hwMatrix::REAL);
        m_status = DinvSqrt.Dimension(m_numParams, m_numParams, hwMatrix::REAL);
        m_status = M.Dimension(m_numParams, m_numParams, hwMatrix::REAL);

        D.Identity();
        DinvSqrt.Identity();
    }
}
//------------------------------------------------------------------------------
// Constructor - Find parameters when solving f(P)=0
// -----------------------------------------------------------------------------
hwGenericFuncFitter::hwGenericFuncFitter(const LSqFitFunc pObjFunc,
                                         const LSqFitFunc pJacFunc,
                                         const hwMatrix&  P,
                                         int              numEqns,
                                         int              maxIter,
                                         int              maxFuncEval,
                                         double           tolf,
                                         double           tolx,
                                         const hwMatrix*  userData)
    : hwGaussNewtLSqFit(P, numEqns, maxIter, maxFuncEval, tolf, tolx),
      m_lowerBound(nullptr),
      m_upperBound(nullptr),
      X(nullptr)
{
    if (!m_status.IsOk())
    {
        switch (m_status.GetArg1())
        {
            case 1: m_status.SetArg1(3); break;
            case 2: m_status.SetArg1(4); break;
            case 3: m_status.SetArg1(5); break;
            case 4: m_status.SetArg1(6); break;
            case 5: m_status.SetArg1(7); break;
            case 6: m_status.SetArg1(8); break;
            default: break;
        }

        if (m_status.GetArg2() == 2)
        {
            m_status.SetArg2(4);
        }
        return;
    }

    if (!pObjFunc)
    {
        m_status(HW_MATH_ERR_NULLPOINTER, 1);
        return;
    }

    if (numEqns < 1)
    {
        m_status(HW_MATH_ERR_NONPOSITIVE, 4);
        return;
    }

    m_pObjFunc = pObjFunc;
    m_pJacFunc = pJacFunc;
    X = new hwMatrix;           // dummy, not used with this constructor
    y = nullptr;                // not used with this constructor
    m_userData = userData;
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwGenericFuncFitter::~hwGenericFuncFitter()
{
    if (X)
    {
        if (X->M() == 0)    // dummy matrix from second constructor
            delete X;
    }
}
//------------------------------------------------------------------------------
// Compute parameters that minimize the objective function
//------------------------------------------------------------------------------
hwMathStatus hwGenericFuncFitter::Compute()
{
    if (!m_lowerBound || !m_upperBound)
    {
        m_status = hwTrustRegionMinimizer::Compute();
        return m_status;
    }

    if (!m_status.IsOk())
    {
        return m_status;
    }

    // get initial estimates for parameters and objective function
    hwMatrix Pcand(m_numParams, 1, hwMatrix::REAL);   // candidate parameter estimates
    GetInitialEstimate(Pcand);

    double f; // Objective function values
    EvalObjectiveFunc(f);

    if (!m_status.IsOk())
    {
        return m_status;
    }

    UpdateHistory(Pcand, f);

    // allocate
    hwMatrix Pstep(m_numParams, 1, hwMatrix::REAL);   // incremental parameter updates
    hwMatrix Pnew(m_numParams, 1, hwMatrix::REAL);    // candidate parameter estimates
    hwMatrix Gstep(m_numParams, 1, hwMatrix::REAL);   // unconstrained gradient step
    hwMatrix Nstep(m_numParams, 1, hwMatrix::REAL);   // unconstrained Newton step

    bool converged = false;
    bool newEstimate;
    bool trBoundary;

    double a;
    double b;
    double c;
    double f_new;                 // objective function values
    double m_new;                 // objective function values
    double gradLength;
    double GstepLength = -1.0;
    double NstepLength;
    double PstepLength;
    double PcandLength;
    double trRadius = 1.0;
    double trRadiusMin = 1.0e-14;
    double trRadiusMax = 10000.0;
    double trRatio = 0.0;
    double numer;
    double denom;
    double zero1;
    double zero2;

    // minimze f
    int i;
    for (i = 0; i < m_maxIter; i++)
    {
        if (i == 0) // step 1: manage the gradient
        {
            EvalGradient();
            if (!m_status.IsOk())
            {
                return m_status;
            }
            newEstimate = true;

            if (m_numFuncEvals <= 0)
            {
                m_status(HW_MATH_WARN_MAXFUNCEVAL);
                break;
            }
        }
        else if (trRatio > 0.25)
        {
            EvalGradient();
            if (m_status.IsOk())
            {
                newEstimate = true;
            }
            else
            {
                // reject new step and update trust region radius
                // Pstep retains previous value
                SetParams(Pcand);

                if (trRadius <= trRadiusMin)
                {
                    m_status(HW_MATH_INFO_SMALLTRUST);
                    break;
                }

                trRadius = _max(0.25 * PstepLength, trRadiusMin);
                newEstimate = false;

                m_status = (D * G).L2Norm(gradLength);
                if (!m_status.IsOk())
                {
                    return m_status;
                }
            }

            if (m_numFuncEvals <= 0)
            {
                m_status(HW_MATH_WARN_MAXFUNCEVAL);
                break;
            }
        }
        else
        {
            // reject new step and update trust region radius
            // Pstep, G retain previous values
            SetParams(Pcand);

            if (trRadius <= trRadiusMin)
            {
                m_status(HW_MATH_INFO_SMALLTRUST);
                break;
            }
            trRadius = _max(0.25 * PstepLength, trRadiusMin);
            newEstimate = false;
        }

        // step 2: manage the Hessian
        if (newEstimate)
        {
            ScaleMatrix();

            // check for gradient convergence
            m_status = (D * G).L2Norm(gradLength);
            if (!m_status.IsOk())
            {
                return m_status;
            }

            if (i == 0)
            {
                if (gradLength == 0.0)
                {
                    converged = true;
                    m_status(HW_MATH_INFO_TOLXCONV);
                    break;
                }
            }

            EvalHessian();
            EvalNewtonStep(Nstep);

            if (!m_status.IsOk())
            {
                if (i == 0)
                {
                    return m_status(HW_MATH_ERR_NOLOCALMIN);
                }

                // reject new step and update trust region radius
                // Pstep retains previous value
                SetParams(Pcand);

                if (trRadius <= trRadiusMin)
                {
                    return m_status(HW_MATH_INFO_SMALLTRUST);
                }

                trRadius = _max(0.25 * PstepLength, trRadiusMin);
                newEstimate = false;

                m_status = (D * G).L2Norm(gradLength);

                if (!m_status.IsOk())
                {
                    return m_status;
                }

                EvalNewtonStep(Nstep);
            }

            m_status = (Nstep).L2Norm(NstepLength);

            if (!m_status.IsOk())
            {
                return m_status;
            }
        }

        if (newEstimate)
        {
            if (i > 0)
            {
                // accept new step and update trust region radius
                UpdateHistory(Pnew, f_new);

                if (fabs(f_new) > 1.0e-12)
                {
                    if (fabs(f - f_new) < m_tolf * fabs(f))
                    {
                        f = f_new;
                        SetParams(Pnew);
                        converged = true;
                        m_status(HW_MATH_INFO_TOLFCONV_R);
                        break;
                    }
                }
                else if (fabs(f - f_new) < m_tolf)
                {
                    f = f_new;
                    SetParams(Pnew);
                    converged = true;
                    m_status(HW_MATH_INFO_TOLFCONV);
                    break;
                }

                f = f_new;
                Pcand = Pnew;

                if (trBoundary == true)
                {
                    if (trRatio > 0.75)
                    {
                        trRadius = _min(2.0 * trRadius, trRadiusMax);
                    }
                }
            }
        }

        // step 3: dogleg analysis
        if (NstepLength <= trRadius)
        {
            Pstep = Nstep;
            trBoundary = false;
        }
        else
        {
            // compute unconstrained steepest descent (gradient) step
            if (newEstimate || GstepLength == -1.0)
            {
                EvalSteepDescentStep(Gstep);

                if (!m_status.IsOk())
                {
                    if (m_status == HW_MATH_ERR_DIVIDEZERO)
                    {
                        return m_status(HW_MATH_ERR_NOTCONVERGE);
                    }
                    return m_status;
                }

                m_status = Gstep.L2Norm(GstepLength);
                if (!m_status.IsOk())
                {
                    return m_status;
                }
            }

            if (GstepLength >= trRadius)
            {
                // use Cauchy point
                Pstep = (D * G) * (trRadius / gradLength);
                trBoundary = true;
            }
            else
            {
                // use dogleg
                m_status = (Nstep - Gstep).L2NormSq(a);
                if (!m_status.IsOk())
                {
                    return m_status;
                }

                m_status = hwMatrix::Dot(Gstep, Nstep - Gstep, b);
                if (!m_status.IsOk())
                {
                    return m_status;
                }

                b *= 2.0;

                m_status = Gstep.L2NormSq(c);
                if (!m_status.IsOk())
                {
                    return m_status;
                }

                c -= trRadius * trRadius;
                quadraticRoots(a, b, c, zero1, zero2);

                if (zero1 >= 0.0 && zero1 <= 1.0)
                {
                    Pstep = (Gstep + (Nstep - Gstep) * zero1);
                }
                else
                {
                    Pstep = (Gstep + (Nstep - Gstep) * zero2);
                }
            }
        }

        EnforceBounds(Pcand, Pstep);

        // step 4: check for step convergence
        m_status = Pstep.L2Norm(PstepLength);
        if (!m_status.IsOk())
        {
            return m_status;
        }

        m_status = Pcand.L2Norm(PcandLength);
        if (!m_status.IsOk())
        {
            return m_status;
        }

        if (PcandLength > 1.0e-12)
        {
            if (PstepLength < m_tolx * PcandLength)
            {
                converged = true;
                m_status(HW_MATH_INFO_TOLXCONV_R);
                Pnew = Pcand - Pstep;
                SetParams(Pnew);
                break;
            }
        }
        else if (PstepLength < m_tolx)
        {
            converged = true;
            m_status(HW_MATH_INFO_TOLXCONV);
            Pnew = Pcand - Pstep;
            SetParams(Pnew);
            break;
        }

        // step 5: prepare next iteration
        EvalQuadModelFunc(f, Pstep, m_new);  // estimate objective function update
        if (!m_status.IsOk())
        {
            return m_status;
        }

        if (m_new < f)                              // update must be improvement 
        {
            // compute new candidate parameters
            Pnew = Pcand - Pstep;
            SetParams(Pnew);

            // update trust region ratio
            EvalObjectiveFunc(f_new);

            if (m_status.IsOk())
            {
                numer = f - f_new;
                denom = f - m_new;

                if (numer > 0.0)
                {
                    if (numer < 1.0e+8 * denom)     // blow up check
                    {
                        if (numer > 1.0e-8 || denom > 1.0e-8)
                        {
                            trRatio = numer / denom;
                        }
                        else
                        {
                            trRatio = 1.0;          // close to 0/0
                        }
                    }
                    else if (numer < 1.0e-8)
                    {
                        trRatio = 1.0;
                    }
                    else
                    {
                        trRatio = 0.0;
                    }
                }
                else
                {
                    trRatio = 0.0;
                }
            }
            else
            {
                trRatio = 0.0;
            }
        }
        else if (m_new == f)
        {
            break;
        }
        else
        {
            trRatio = 0.0;
        }
    }

    m_numIter = i + 1;
    m_objFuncVal = f;

    if (!converged)
    {
        if (i == m_maxIter)
        {
            SetParams(Pcand);   // return to last history entry
            m_status(HW_MATH_WARN_MAXITERATE);
        }
        else if (m_status == HW_MATH_WARN_MAXFUNCEVAL)
        {
            SetParams(Pcand);   // return to last history entry
        }
        else
        {
            m_status(HW_MATH_WARN_NOTCONVERGE);
        }
    }

    return m_status;
}
//------------------------------------------------------------------------------
// Evaluate fitted function
//------------------------------------------------------------------------------
void hwGenericFuncFitter::EvalFittedFunc(hwMatrix& y_est)
{
    if (y_est.IsEmpty())
    {
        m_status(HW_MATH_ERR_EMPTYMATRIX, 1);
        return;
    }

    if (!y_est.IsReal())
    {
        m_status(HW_MATH_ERR_COMPLEX, 1);
        return;
    }

    if (!y_est.IsVector())
    {
        m_status(HW_MATH_ERR_VECTOR, 1);
        return;
    }

    if (m_numFuncEvals == 0)
    {
        return;
    }

    m_status = m_pObjFunc(P, *X, m_userData, y_est);
    --m_numFuncEvals;
}
//------------------------------------------------------------------------------
// Evaluate residuals
//------------------------------------------------------------------------------
void hwGenericFuncFitter::EvalResiduals(hwMatrix& residual)
{
    EvalFittedFunc(residual);

    if (!m_status.IsOk())
    {
        return;
    }

    if (y)
    {
        for (int i = 0; i < m_numEqns; i++)     // m_numEqns = m_numPnts
        {
            residual(i) -= (*y)(i);
        }
    }
}
//------------------------------------------------------------------------------
// Evaluate Jacobian matrix
//------------------------------------------------------------------------------
void hwGenericFuncFitter::EvalJacobian()
{
    if (m_pJacFunc)
    {
        // use analytical derivatives
        m_status = m_pJacFunc(P, *X, m_userData, J);
    }
    else
    {
        // use numerical derivatives
        double delta;
        double param;
        double temp;
        hwMatrix y_est(m_numEqns, hwMatrix::REAL);

        for (int j = 0; j < m_numParams; j++)
        {
            // find delta
            param = GetParam(j);
            delta = _max(m_numDeriv_eps * fabs(param), 1.0e-10); // 2.0 * MACHEP);
            temp  = param + delta;
            delta = temp - param;

            // compute central difference
            SetParam(j, param + delta);
            EvalFittedFunc(y_est);

            if (!m_status.IsOk())
            {
                if (m_status.GetArg1() != 111)
                    m_status.ResetArgs();

                return;
            }

            for (int i = 0; i < m_numEqns; i++)
                J(i, j) = y_est(i);

            SetParam(j, param - delta);
            EvalFittedFunc(y_est);

            if (!m_status.IsOk())
            {
                if (m_status.GetArg1() != 111)
                    m_status.ResetArgs();

                return;
            }

            for (int i = 0; i < m_numEqns; i++)
            {
                J(i, j) -= y_est(i);
                J(i, j) /= (2.0 * delta);
            }

            SetParam(j, param);         // restore original value
        }
    }
}

/*
Sources for bounded lsqcurvefit concepts:
1. J.B. Francisco, N. Krejic, M.Martínez, "An interior-point method for solving box-constrained underdetermined nonlinear systems"
   Journal of Computational and Applied Mathematics
   Volume 177, Issue 1, 1 May 2005, Pages 67-88

2. S. Bellavia, M. Macconi, S. Pieraccini, "Constrained Dogleg Methods for nonlinear systems with simple bounds"
   Computational Optimization and Applications, Vol. 53, pp. 771-794, 2012

3. A. Klug, "Affine-Scaling Methods for Nonlinear Minimization Problems and Nonlinear Systems of Equations with Bound Constraints"
   PhD Dissertation, Bavarian Julius Maximilians University, Wurzburg, 2006

4. S. Bellavia, M. Macconi, S. Pieraccini, "On affine scaling inexact dogleg methods for bound-constrained nonlinear systems"
   Optimization Methods and Software. ISSN 1055-6788, 30:2(2015), pp. 276-300.
*/

//------------------------------------------------------------------------------
// Evaluate approximate Hessian matrix
//------------------------------------------------------------------------------
void hwGenericFuncFitter::EvalHessian()
{
    hwGaussNewtLSqFit::EvalHessian();

    if (!m_lowerBound && !m_upperBound)
        return;

    hwMatrix Jd(m_numParams, m_numParams, hwMatrix::REAL);
    hwMatrix Gd(m_numParams, m_numParams, hwMatrix::REAL);

    Jd.SetElements(0.0);
    Gd.SetElements(0.0);

    for (int j = 0; j < m_numParams; j++)
    {
        if (G(j) < 0.0)
        {
            if (!m_upperBound || (*m_upperBound)(j) == std::numeric_limits<double>::infinity())
            {
                Jd(j, j) = 0.0;
            }
            else
            {
                Jd(j, j) = -1.0;
            }
        }
        else if (G(j) > 0.0)
        {
            if (!m_lowerBound || (*m_lowerBound)(j) == -std::numeric_limits<double>::infinity())
            {
                Jd(j, j) = 0.0;
            }
            else
            {
                Jd(j, j) = 1.0;
            }
        }
        else
        {
            Jd(j, j) = 0.0;
        }

        Gd(j, j) = G(j);
    }

    hwMatrix C;

    C = DinvSqrt * Gd * Jd * DinvSqrt;

    if (m_numEqns == m_numParams)
    {
        // unlikely case
        hwMatrix JT;
        JT.Transpose(J);
        B = JT * J;
    }

    M = B + C;
}
//------------------------------------------------------------------------------
// Evaluate steepest descent step
//------------------------------------------------------------------------------
void hwGenericFuncFitter::EvalSteepDescentStep(hwMatrix& Gstep)
{
    if (!m_lowerBound && !m_upperBound)
    {
        hwGaussNewtLSqFit::EvalSteepDescentStep(Gstep);
        return;
    }

    double numer;
    m_status = (D * G).L2NormSq(numer);      // = G * G;

    if (!m_status.IsOk())
    {
        return;
    }

    double denom;

    hwMatrix GT;
    m_status = GT.Transpose(G);
    denom = (GT * M * G)(0);  // = G * B * G;

    if (denom < 1.0e-12 * numer)
    {
        m_status(HW_MATH_ERR_DIVIDEZERO);
        return;
    }

    Gstep = G * (numer / denom);
}
//------------------------------------------------------------------------------
// Evaluate steepest descent step
//------------------------------------------------------------------------------
void hwGenericFuncFitter::EvalNewtonStep(hwMatrix& Nstep)
{
    if (!m_lowerBound && !m_upperBound)
    {
        hwGaussNewtLSqFit::EvalNewtonStep(Nstep);
        return;
    }

    // compute unconstrained Newton step
    m_status = Nstep.LSolveSPD(M, G);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
    }
}
//------------------------------------------------------------------------------
// Scale matrix
//------------------------------------------------------------------------------
void hwGenericFuncFitter::ScaleMatrix()
{
    for (int j = 0; j < m_numParams; j++)
    {
        double v;

        if (G(j) < 0.0)
        {
            if (!m_upperBound || (*m_upperBound)(j) == std::numeric_limits<double>::infinity())
            {
                v = 1.0;
            }
            else
            {
                v = (*m_upperBound)(j) - P(j);
            }
        }
        else if (G(j) > 0.0)
        {
            if (!m_lowerBound || (*m_lowerBound)(j) == -std::numeric_limits<double>::infinity())
            {
                v = 1.0;
            }
            else
            {
                v = P(j) - (*m_lowerBound)(j);
            }
        }
        else
        {
            if ((!m_upperBound || (*m_upperBound)(j) == std::numeric_limits<double>::infinity()) &&
                (!m_lowerBound || (*m_lowerBound)(j) == -std::numeric_limits<double>::infinity()))
            {
                v = 1.0;
            }
            else
            {
                v = _min((*m_upperBound)(j) - P(j), P(j) - (*m_lowerBound)(j));
            }
        }

        if (v != 0.0)
        {
            D(j, j) = v;
            DinvSqrt(j, j) = 1.0 / sqrt(v);
        }
        else
        {
            D(j, j) = 1.0;
            DinvSqrt(j, j) = 1.0;
        }
    }
}
//------------------------------------------------------------------------------
// Constrain step to simple bounds limits
//------------------------------------------------------------------------------
void hwGenericFuncFitter::EnforceBounds(const hwMatrix& Pcand, hwMatrix& Pstep)
{
    if (!m_lowerBound && !m_upperBound)
        return;

    bool truncate = false;

    // project onto the bounding box and then truncate to an interior point
    for (int i = 0; i < m_numParams; i++)
    {
        if (m_lowerBound && (Pcand(i) - Pstep(i) < (*m_lowerBound)(i)))
        {
            Pstep(i) = Pcand(i) - (*m_lowerBound)(i);
            truncate = true;
        }
        else if (m_upperBound && (Pcand(i) - Pstep(i) > (*m_upperBound)(i)))
        {
            Pstep(i) = Pcand(i) - (*m_upperBound)(i);
            truncate = true;
        }
    }

    if (truncate)
    {
        double norm;

        m_status = Pstep.L2Norm(norm);
        double alpha = _max(0.995, (1.0 - norm));
        Pstep *= alpha;
    }
}
