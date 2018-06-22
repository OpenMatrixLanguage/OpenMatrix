/**
* @file  hwElliptic_Proto.cxx
* @date April 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
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
#include "hwElliptic_Proto.h"

#include <math.h>

#include "EllipticFuncs.h"
#include "GeneralFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwElliptic_Proto::hwElliptic_Proto(int    order, 
                                   double passEdgeDb,
                                   double stopEdgeDb)
	: hwLowPass_Proto(order)
{
    if (!m_status.IsOk())
    {
        return;
    }
    if (passEdgeDb <= 0.0)
    {
        m_status(HW_MATH_ERR_DB_SIGN, 2);
        return;
    }

    if (stopEdgeDb <= 0.0)
    {
        m_status(HW_MATH_ERR_DB_SIGN, 3);
        return;
    }

    if (stopEdgeDb < passEdgeDb + 0.001)
    {
        m_status(HW_MATH_ERR_FILTERRIPPLE, 2, 3);
        return;
    }

    double dbfac = log(10.0) / 10.0;
    m_epsilon = sqrt(exp(passEdgeDb * dbfac) - 1.0);

    if (IsZero(m_epsilon, 1.0e-12))
    {
        m_status(HW_MATH_ERR_FILTERSPEC_E, 2);
        return;
    }

    double delta = sqrt(exp(stopEdgeDb * dbfac) - 1.0);

    if (IsZero(delta, 1.0e-12))
    {
        m_status(HW_MATH_ERR_FILTERSPEC_E, 3);
        return;
    }

    double m1 = m_epsilon / delta;

    m1 *= m1;

    double m1p = 1.0 - m1;
    double Kk1 = ellpk(m1p);
    double Kpk1 = ellpk(m1);                    // Kp denotes K prime
    double q = exp(-PI * Kpk1 / (m_order * Kk1));
    m_k = cay(q);                               // m_k is the elliptic modulus

    LambdaPln();
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwElliptic_Proto::~hwElliptic_Proto()
{
}
//------------------------------------------------------------------------------
// Compute location of the real pole
//------------------------------------------------------------------------------
void hwElliptic_Proto::GetSPlaneInfo(double& poleReal) const
{
    poleReal = -m_sn / m_cn;
}
//------------------------------------------------------------------------------
// Compute real component and squared magnitude of the ith pole
//------------------------------------------------------------------------------
void hwElliptic_Proto::GetSPlaneInfo(int     i, 
                                     double& poleReal, 
                                     double& poleMagSq,
                                     double& zeroMagSq) const
{
    // complex conjugate pole pair
    double a = m_order - 1 - i - i;
    double b = a * m_Kk / m_order;	

    double sn;
    double cn;
    double dn;
    double phi;
    ellpj(b, m_m, sn, cn, dn, phi);

    // pole real component
    double r = m_k * sn * m_sn;
    b = m_cn * m_cn + r * r;
    a = -cn * dn * m_sn * m_cn / b;
    poleReal = a;

    // pole mag^2
    b = sn * m_dn / b;
    poleMagSq = ((a * a) + (b * b));

    // zero mag^2
    b = 1.0 / (m_k * sn);
    zeroMagSq = b * b;
}
//------------------------------------------------------------------------------
// Compute low pass ripple factor at DC
//------------------------------------------------------------------------------
double hwElliptic_Proto::GetRippleFactor() const
{
    if (m_order%2 == 0)
    {
        return 1.0 / sqrt(1.0 + m_epsilon * m_epsilon);
    }
    return 1.0;
}
//------------------------------------------------------------------------------
// Compute Jacobian elliptic information in the lambda plane
//------------------------------------------------------------------------------
void hwElliptic_Proto::LambdaPln()
{
    // This function is adapted from public domain code
    // by Stephen Moshier available at http://www.moshier.net/

    m_m  = m_k * m_k;
    m_Kk = ellpk(1.0 - m_m);

    double Kpk = ellpk(m_m);
    double q   = exp(-PI * m_order * Kpk / m_Kk);  // the nome of k1
    double m1  = cay(q);                           // modulus for the nome
    m1 *= m1;                                      // m1 is now properly the elliptic parameter

    double m1p = 1.0 - m1;                // m1p is the complementaty elliptic parameter
    double Kk1 = ellpk(m1p);
    double phi = atan(1.0 / m_epsilon);   // phi is the Jacobi amplitude
    double u = ellik(phi, m1p);
    u *= m_Kk / (m_order * Kk1);	      // Jacobian elliptic function argument (or, u *= Kpk / Kpk1)

    ellpj(u, 1.0 - m_m, m_sn, m_cn, m_dn, phi); // dummy use of phi
}

