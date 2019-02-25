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

#include <hwGenericFuncFitter.h>

//------------------------------------------------------------------------------
// Constructor - Find parameters when fitting equations of the form y=f(x)
//------------------------------------------------------------------------------
hwGenericFuncFitter::hwGenericFuncFitter(const LSqFitFunc pObjFunc,
                                         const LSqFitFunc pJacFunc,
                                         const hwMatrix&  P,
                                         const hwMatrix&  X_,
                                         const hwMatrix&  y_,
                                         int              maxIter,
                                         int              maxFuncEval,
                                         double           tolf,
                                         double           tolx,
                                         const hwMatrix*  userData)
    : hwGaussNewtLSqFit(P, y_.Size(), maxIter, maxFuncEval, tolf, tolx),
      X (nullptr)
{
    if (!m_status.IsOk())
    {
        switch (m_status.GetArg1())
        {
            case 1: m_status.SetArg1(3); break;
            case 2: m_status.SetArg1(5); break;
            case 3: m_status.SetArg1(6); break;
            case 4: m_status.SetArg1(7); break;
            case 5: m_status.SetArg1(8); break;
            case 6: m_status.SetArg1(9); break;
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
      X (nullptr)
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
