/**
* @file hwUnconstrMinimizer.cxx
* @date June, 2007
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

#include <hwUnconstrMinimizer.h>

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwUnconstrMinimizer::hwUnconstrMinimizer(const UnConMinObjFunc  pObjFunc,
                                         const UnConMinGradFunc pGradFunc,
                                         const hwMatrix&        P, 
                                         int                    maxIter,
                                         int                    maxFuncEval, 
                                         double                 tolf, 
                                         double                 tolx,
                                         const hwMatrix*        userData)
    : hwTrustRegionMinimizer(P, maxIter, maxFuncEval, tolf, tolx)
{
    if (!m_status.IsOk())
    {
        switch (m_status.GetArg1())
        {
            case 1:  m_status.SetArg1(3); break;
            case 2:  m_status.SetArg1(4); break;
            case 3:  m_status.SetArg1(5); break;
            case 4:  m_status.SetArg1(6); break;
            case 5:  m_status.SetArg1(7); break;
            default: break;
        }
        return;
    }

    if (!pObjFunc)
    {
        m_status(HW_MATH_ERR_NULLPOINTER, 1);
        return;
    }

    initStep    = true;
    m_pObjFunc  = pObjFunc;
    m_pGradFunc = pGradFunc;
    m_userData  = userData;
    Gk.Dimension(m_numParams, 1, hwMatrix::REAL);	         // gradient at previous step
    Gb.Dimension(m_numParams, 1, hwMatrix::REAL);		     // backup of Gk
    Pk.Dimension(m_numParams, 1, hwMatrix::REAL);	         // parameter estimates at previous step
    B.Dimension(m_numParams,  m_numParams, hwMatrix::REAL);	 // approximate Hessian
    Bk.Dimension(m_numParams, m_numParams, hwMatrix::REAL);  // Hessian approx at previous step
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwUnconstrMinimizer::~hwUnconstrMinimizer()
{
}
//------------------------------------------------------------------------------
// Evaluate objective function
//------------------------------------------------------------------------------
void hwUnconstrMinimizer::EvalObjectiveFunc(double& result)
{
    if (m_numFuncEvals == 0)
        return;

    m_status = m_pObjFunc(P, m_userData, result);
    --m_numFuncEvals;
}
//------------------------------------------------------------------------------
// Evaluate gradient vector
//------------------------------------------------------------------------------
void hwUnconstrMinimizer::EvalGradient()
{
    if (m_pGradFunc)
    {
        // use analytical derivatives
        m_status = m_pGradFunc(P, m_userData, G);
    }
    else
    {
        // use numerical derivatives
        double delta;
        double param;
        double result;
        double temp;

        for (int i = 0; i < m_numParams; i++)
        {
            // find delta
            param = GetParam(i);
            delta = _max(m_numDeriv_eps * fabs(param), 1.0e-6);
            temp = param + delta;
            delta = temp - param;

            // compute central difference
            SetParam(i, param + delta);
            EvalObjectiveFunc(result);

            if (!m_status.IsOk())
            {
                return;
            }
            G(i) = result;

            SetParam(i, param - delta);
            EvalObjectiveFunc(result);

            if (!m_status.IsOk())
            {
                return;
            }
            G(i) -= result;
            G(i) /= (2.0 * delta);

            SetParam(i, param);         // restore original value
        }
    }
}
//------------------------------------------------------------------------------
// Evaluate steepest descent step
//------------------------------------------------------------------------------
void hwUnconstrMinimizer::EvalSteepDescentStep(hwMatrix& Gstep)
{
    double numer;
    m_status = hwMatrix::Dot(G, G, numer);
    if (!m_status.IsOk())
    {
        return;
    }

    double denom;
    m_status = hwMatrix::Dot(G, B * G, denom);
    if (!m_status.IsOk())
    {
        return;
    }

    if (denom < 1.0e-12 * numer)
    {
        m_status(HW_MATH_ERR_DIVIDEZERO);
        return;
    }

    Gstep = G * (numer / denom);
}
//------------------------------------------------------------------------------
// Evaluate approximate Hessian matrix
//------------------------------------------------------------------------------
void hwUnconstrMinimizer::EvalHessian()
{
    if (initStep == true)
    {
        // initialize BFGS approximation
        B.Identity();
        Gk = G;
    }
    else
    {
        // damped BFGS
        hwMatrix yk = G - Gk;
        hwMatrix ykT;
        m_status = ykT.Transpose(yk);
        if (!m_status.IsOk())
        {
            return;
        }

        hwMatrix sk = P - Pk;
        hwMatrix skT;
        m_status = skT.Transpose(sk);
        if (!m_status.IsOk())
        {
            return;
        }

        hwMatrix Hk;
        m_status = Hk.Inverse(Bk);
        if (!m_status.IsOk())
        {
            return;
        }

        // enforce Wolfe condition
        double skTyk = (skT * yk)(0);
        if (skTyk < 1.0e-15)
        {
            m_status(HW_MATH_ERR_INVALIDINPUT);
            return;
        }

        double ykTHkyk;
                      
        hwMatrix skTBk;
        hwMatrix BkskskTB;
        hwMatrix wk;
        hwMatrix wkT;

        // update yk with damping
        // "Improved Damped Quasi-Newton Methods for Unconstrained
        // Optimization", Al_Baali and Grandinetti, August 2015
        hwMatrix Bksk = Bk * sk;
        double skTBksk = (skT * Bksk)(0);
        double bk      = skTBksk / skTyk;
        double bkr     = 1 / bk;
        hwMatrix Hkyk = Hk * yk;
        hwMatrix ykykT;
        ykTHkyk = (ykT * Hkyk)(0);
        double hk = ykTHkyk / skTyk;
        double thetak = (hk < 0.95) ? (1.0 / (1.0 - bk)) : 0.0;
        double ak     = (bk * hk - 1.0) * _max(fabs(thetak), 1.0);
        double sigma2 = _max(1.0 - 1.0 / ak, 0.5);
        double sigma3 = exp(1.0);
        double phi = 1.0;

        if (bkr < 1.0 - sigma2)
        {
            phi = sigma2 / (1.0 - bkr);
        }
        else if (bkr > 1.0 + sigma3)
        {
            phi = sigma3 / (bkr - 1.0);
        }
        else if (ak > sigma3)
        {
            phi = sqrt(sigma3 / ak);
        }

        yk = yk * phi + Bksk * (1.0 - phi);

        m_status = ykT.Transpose(yk);
        if (!m_status.IsOk())
        {
            return;
        }

        skTyk = (skT * yk)(0);

        // approximate Hessian
        m_status = skTBk.Transpose(Bksk);
        if (!m_status.IsOk())
        {
            return;
        }

        wk = yk / skTyk - Bksk / skTBksk;

        m_status = wkT.Transpose(wk);
        if (!m_status.IsOk())
        {
            return;
        }

        B = Bk - (Bksk * skTBk) / skTBksk + (yk * ykT) / skTyk
            + thetak * skTBksk * (wk * wkT);
    }

    // save for next interation
    Pk = P;
    Gb = Gk;
    Gk = G;
    Bk = B;
}
//------------------------------------------------------------------------------
// Evaluate Newton step
//------------------------------------------------------------------------------
void hwUnconstrMinimizer::EvalNewtonStep(hwMatrix& Nstep)
{
    if (initStep == true)
    {
        Nstep    = G;
        initStep = false;
        return;
    }

    m_status = Nstep.LSolveSPD(B, G);
}
//------------------------------------------------------------------------------
// Evaluate trust region quadratic model
//!
// -------------------------------------------------------------------
void hwUnconstrMinimizer::EvalQuadModelFunc(double          f,
                                            const hwMatrix& Pstep,
                                            double&         m_new)
{
    // quadratic model of objective function
    double linearTerm;
    m_status = hwMatrix::Dot(G, Pstep, linearTerm);
    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return;
    }

    double quadTerm;
    m_status = hwMatrix::Dot(Pstep, B * Pstep, quadTerm);   // Pstep * B * Pstep;
    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return;
    }

    m_new = f - linearTerm + 0.5 * quadTerm;    // minus due to sign of Pstep in Compute()
}
