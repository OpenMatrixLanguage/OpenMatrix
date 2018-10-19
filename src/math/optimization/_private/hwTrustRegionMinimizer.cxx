/**
* @file hwTrustRegionMinimizer.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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
#include <hwTrustRegionMinimizer.h>

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwTrustRegionMinimizer::hwTrustRegionMinimizer(const hwMatrix& P_, 
                                               int             maxIter,
                                               int             maxFuncEval, 
                                               double          tolf, 
                                               double          tolx)
    : Obj_history    (nullptr),
      DV_history     (nullptr),
      m_maxIter      (maxIter),
      m_numFuncEvals (maxFuncEval),
      m_tolf         (tolf),
      m_tolx         (tolx),
      m_numIter      (0),
      m_numHistPnts  (0),
      m_objFuncVal   (0.0)
{
    if (P_.IsEmpty())
    {
        m_status(HW_MATH_ERR_EMPTYMATRIX, 1);
        return;
    }

    if (!P_.IsReal())
    {
        m_status(HW_MATH_ERR_COMPLEX, 1);
        return;
    }

    if (!P_.IsVector())
    {
        m_status(HW_MATH_ERR_VECTOR, 1);
        return;
    }

    if (maxIter < 0)
    {
        m_status(HW_MATH_ERR_NONPOSINT, 2);
        return;
    }

    if (maxFuncEval < 0)
    {
        m_status(HW_MATH_ERR_NONPOSINT, 3);
        return;
    }

    if (tolf <= 0.0)
    {
        m_status(HW_MATH_ERR_NONPOSITIVE, 4);
        return;
    }

    if (tolx <= 0.0)
    {
        m_status(HW_MATH_ERR_NONPOSITIVE, 5);
        return;
    }

    m_numDeriv_eps = pow(MachPrecision(1.0), 1.0/3.0);
    m_numParams    = P_.Size();

    m_status = P.Dimension(m_numParams, 1, hwMatrix::REAL);	// parameters
    if (!m_status.IsOk())
    {
        return;
    }

    m_status = G.Dimension(m_numParams, 1, hwMatrix::REAL);	// gradient
    if (!m_status.IsOk())
    {
        return;
    }

    for (int i = 0; i < m_numParams; i++)   // allows P_ to be row or column
    {
        P(i) = P_(i);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwTrustRegionMinimizer::~hwTrustRegionMinimizer()
{
}
//------------------------------------------------------------------------------
// Compute parameters that minimize the objective function
//------------------------------------------------------------------------------
hwMathStatus hwTrustRegionMinimizer::Compute()
{
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
    bool GstepSet = false;

    double a;
    double b;
    double c;
    double f_new;                 // objective function values
    double m_new;                 // objective function values
    double gradLength;
    double GstepLength;
    double NstepLength;
    double PstepLength;
    double PcandLength;
    double trRadius    = 1.0;
    double trRadiusMin = 1.0e-14;
    double trRadiusMax = 10000.0;
    double trRatio     = 0.0;
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

                trRadius    = _max(0.25 * PstepLength, trRadiusMin);
                newEstimate = false;
                ResetGradient();

                m_status = G.L2Norm(gradLength);
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
            trRadius    = _max(0.25 * PstepLength, trRadiusMin);
            newEstimate = false;
        }

        // step 2: manage the Hessian
        if (newEstimate)
        {
            // first check for gradient convergence
            m_status = G.L2Norm(gradLength);
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
                    // return m_status(HW_MATH_ERR_NOLOCALMIN);
                }

                trRadius    = _max(0.25 * PstepLength, trRadiusMin);
                newEstimate = false;
                ResetGradient();

                m_status = G.L2Norm(gradLength);

                if (!m_status.IsOk())
                {
                    return m_status;
                }
                EvalNewtonStep(Nstep);
            }

            m_status = Nstep.L2Norm(NstepLength);

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

                if (fabs(f_new) > _max(1.0e-6 * m_tolf, 1.0e-20))
                {
                    if (fabs(f - f_new) < m_tolf * fabs(f))
                    {
                        f = f_new;
                        SetParams(Pnew);

                        if (m_tolf < 1.0e-20 || gradLength * m_tolx < 100.0 * m_tolf)   // sanity check on gradient
                        {
                            converged = true;
                            m_status(HW_MATH_INFO_TOLFCONV_R);
                            break;
                        }
                    }
                }
                else
                {
                    if (fabs(f - f_new) < m_tolf)
                    {
                        f = f_new;
                        SetParams(Pnew);

                        if (m_tolf < 1.0e-20 || gradLength * m_tolx < 100.0 * m_tolf)   // sanity check on gradient
                        {
                            converged = true;
                            m_status(HW_MATH_INFO_TOLFCONV);
                            break;
                        }
                    }
                }

                f     = f_new;
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
            if (newEstimate || !GstepSet)
            {
                EvalSteepDescentStep(Gstep);
                GstepSet = true;

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
                Pstep = G * (trRadius / gradLength);
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

            trBoundary = true;
        }

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

        if (PcandLength > _max(1.0e-6 * m_tolx, 1.0e-20))
        {
            if (PstepLength < m_tolx * PcandLength)
            {
                if (m_tolf < 1.0e-20 || gradLength * PstepLength < 100.0 * m_tolf)   // sanity check on gradient
                {
                    converged = true;
                    m_status(HW_MATH_INFO_TOLXCONV_R);
                    Pnew = Pcand - Pstep;
                    SetParams(Pnew);
                    break;
                }
                else if (PstepLength <= trRadiusMin)
                {
                    converged = true;
                    m_status(HW_MATH_INFO_TOLXCONV);
                    Pnew = Pcand - Pstep;
                    SetParams(Pnew);
                    break;
                }
            }
        }
        else
        {
            if (PstepLength < m_tolx)
            {
                if (m_tolf < 1.0e-20 || gradLength * PstepLength < 100.0 * m_tolf)   // sanity check on gradient
                {
                    converged = true;
                    m_status(HW_MATH_INFO_TOLXCONV);
                    Pnew = Pcand - Pstep;
                    SetParams(Pnew);
                    break;
                }
                else if (PstepLength <= trRadiusMin)
                {
                    converged = true;
                    m_status(HW_MATH_INFO_TOLXCONV);
                    Pnew = Pcand - Pstep;
                    SetParams(Pnew);
                    break;
                }
            }
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

    m_numIter    = i + 1;
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
// Get parameter at the specified index
//------------------------------------------------------------------------------
double hwTrustRegionMinimizer::GetParam(int index) const
{
    // override in derived class if parameters are not stored as a vector
    return P(index);
}
//------------------------------------------------------------------------------
// Get parameters as a vector
//------------------------------------------------------------------------------
void hwTrustRegionMinimizer::GetParams(hwMatrix& Pcand)
{
    // override in derived class if parameters are not stored as a vector
    // Pcand = P;
    if (!Pcand.IsVector())
    {
        m_status(HW_MATH_ERR_VECTOR);
        return;
    }

    if (Pcand.Size() != m_numParams)
    {
        m_status(HW_MATH_ERR_ARRAYSIZE);
        return;
    }

    for (int i = 0; i < m_numParams; i++)   // allows Pcand to be row or column
    {
        Pcand(i) = P(i);
    }
}
//------------------------------------------------------------------------------
// Set parameter at the specified index
//------------------------------------------------------------------------------
void hwTrustRegionMinimizer::SetParam(int index,
                                      double value)
{
    // override in derived class if parameters are not stored as a vector
    P(index) = value;
}
//------------------------------------------------------------------------------
// Set parameters as a vector
//------------------------------------------------------------------------------
void hwTrustRegionMinimizer::SetParams(const hwMatrix& Pcand)
{
    // override in derived class if parameters are not stored as a vector
    P = Pcand;
}
//------------------------------------------------------------------------------
// Get initial estimate of the parameter vector
//------------------------------------------------------------------------------
void hwTrustRegionMinimizer::GetInitialEstimate(hwMatrix& Pcand)
{
    // override in derived class if parameters are not stored as a vector
    Pcand = P;
}
//------------------------------------------------------------------------------
// Update iteration history
//------------------------------------------------------------------------------
void hwTrustRegionMinimizer::UpdateHistory(const hwMatrix& Pcand,
                                           double          f)
{
    // make this a pure virtual function if residual history is desired
    if (DV_history || Obj_history)
    {
        if (DV_history)
        {
            if (DV_history->IsEmpty())
            {
                DV_history->Dimension(m_numParams, m_numHistPnts+1, hwMatrix::REAL);
            }
            else
            {
                DV_history->Resize(m_numParams, m_numHistPnts+1, true);
            }
            DV_history->WriteColumn(m_numHistPnts, Pcand);
        }

        if (Obj_history)
        {
            if (Obj_history->IsEmpty())
            {
                Obj_history->Dimension(1, m_numHistPnts+1, hwMatrix::REAL);
            }
            else
            {
                Obj_history->Resize(1, m_numHistPnts+1, true);
            }
            (*Obj_history)(m_numHistPnts) = f;
        }

        m_numHistPnts++;
    }
}
