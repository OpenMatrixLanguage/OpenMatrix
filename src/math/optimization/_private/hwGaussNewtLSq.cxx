/**
* @file hwGaussNewtLSqFit.cxx
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
#include <hwGaussNewtLSq.h>
#include <hwMatrix.h>

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwGaussNewtLSqFit::hwGaussNewtLSqFit(const hwMatrix& P, 
                                     int             numEqns, 
                                     int             maxIter,
                                     int             maxFuncEval, 
                                     double          tolf, 
                                     double          tolx)
    : hwTrustRegionMinimizer(P, maxIter, maxFuncEval, tolf, tolx)
{
    // When finding best fit params for y=f(x), numPoints = numEqns
    if (!m_status.IsOk())
    {
        switch (m_status.GetArg1())
        {
            case 2:  m_status.SetArg1(3); break;
            case 3:  m_status.SetArg1(4); break;
            case 4:  m_status.SetArg1(5); break;
            case 5:  m_status.SetArg1(6); break;
            default: break;
        }
        return;
    }

    if (numEqns < P.Size())
    {
        m_status(HW_MATH_ERR_UNDERDETSYS_E, 1, 2);
        return;
    }

    m_numEqns = numEqns;

    F.Dimension(m_numEqns,    1,           hwMatrix::REAL);
    J.Dimension(m_numEqns,    m_numParams, hwMatrix::REAL);
    JT.Dimension(m_numParams, m_numEqns,   hwMatrix::REAL);

    if (m_numEqns != m_numParams) // approximate Hessian
    {
        B.Dimension(m_numParams, m_numParams, hwMatrix::REAL);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwGaussNewtLSqFit::~hwGaussNewtLSqFit()
{
}
//------------------------------------------------------------------------------
// Evaluate objective function
//------------------------------------------------------------------------------
void hwGaussNewtLSqFit::EvalObjectiveFunc(double& result)
{
    // Non-linear least squares literature commonly includes a 0.5 factor on
    // 'result' that is omitted here. See the effect on EvalQuadModelFunc() below.
    EvalResiduals(F);

    if (!m_status.IsOk())
    {
        return;
    }

    m_status = F.L2NormSq(result);
}
//------------------------------------------------------------------------------
// Evaluate gradient vector
//------------------------------------------------------------------------------
void hwGaussNewtLSqFit::EvalGradient()
{
    EvalJacobian();

    if (!m_status.IsOk())
    {
        return;
    }

    JT.Transpose(J);
    G = JT * F;
}
//------------------------------------------------------------------------------
// Evaluate steepest descent step
//------------------------------------------------------------------------------
void hwGaussNewtLSqFit::EvalSteepDescentStep(hwMatrix& Gstep)
{
    double numer;
    m_status = G.L2NormSq(numer);      // = G * G;

    if (!m_status.IsOk())
    {
        return;
    }

    double denom;
    m_status = (J*G).L2NormSq(denom);  // = G * B * G;

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
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
void hwGaussNewtLSqFit::EvalHessian()
{
    if (m_numEqns == m_numParams)
    {
        // algorithm simplifies to Newton-Raphson, so Hessian is not computed
    }
    else
    {
        B = JT * J; // compute approximate Hessian for Gauss-Newton method
    }

#if 0 // Commented code
    //if (!m_status.IsOk())
    //    m_status.ResetArgs();
#endif
}
//------------------------------------------------------------------------------
// Evaluate Newton step
//------------------------------------------------------------------------------
void hwGaussNewtLSqFit::EvalNewtonStep(hwMatrix& Nstep)
{
    // compute unconstrained Newton step
    m_status = (m_numEqns == m_numParams) ?
               Nstep.LSolve(J, F):
               Nstep.LSolveSPD(B, G);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
    }
}
//------------------------------------------------------------------------------
// Evaluate trust region quadratic model
//------------------------------------------------------------------------------
void hwGaussNewtLSqFit::EvalQuadModelFunc(double          f, 
                                          const hwMatrix& Pstep, 
                                          double&         m_new)
{
    // quadratic model of objective function
    hwMatrix temp = J * Pstep;

    double linearTerm;
    m_status = hwMatrix::Dot(G, Pstep, linearTerm);
    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return;
    }

    double quadTerm;
    m_status = temp.L2NormSq(quadTerm);     // = Pstep * B * Pstep;
    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return;
    }

    // If EvalObjectiveFunc() were to include a 0.5 factor (as is common in
    // the literature) this would be written as
    // m_new = f - linearTerm + 0.5 * quadTerm;
    m_new = f - 2.0 * linearTerm + quadTerm;    // minus due to sign of Pstep in Compute()
}

