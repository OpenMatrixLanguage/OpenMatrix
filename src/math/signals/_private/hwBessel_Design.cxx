/**
* @file hwBessel_Design.cxx
* @date May 2012
* Copyright (C) 2012-2018 Altair Engineering, Inc.  
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
#include "hwBessel_Design.h"

#include "hwMatrix.h"
#include "PolynomFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwBessel_Design::hwBessel_Design()
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwBessel_Design::~hwBessel_Design()
{
}
//------------------------------------------------------------------------------
// Design a low pass prototype filter with normalized -3db cutoff frequencies
//------------------------------------------------------------------------------
void hwBessel_Design::DesignPrototype(double  epsilon, 
                                      double  delta,
                                      double  omegaS2, 
                                      int&    order,
                                      double& omegaC1,
                                      double& omegaC2)
{
    // compute order and range of possible normalized omegaC values
    // omegaP1 = default pass band corner = 1.0
    // omegaP2 = alternate pass band corner > omegaP1
    // omegaS1 = default stop band corner < omegaS2, not computed
    // omegaS2 = alternate stop band corner, the input upper bound
    // omegaSP = omegaS2 / omegaP2 = omegaS1 / omegaP1, adjusted for integer filter size
    double omegaC;
    double omegaP;
    double omegaS;

    // find the order
    int k = 1;
    for (k = 1; k < 15; ++k)
    {
        SolveProtoCutoff(epsilon, delta, k, omegaC, omegaP, omegaS);
        if (!m_status.IsOk())
        {
            return;
        }

        if (omegaS/omegaP <= omegaS2)
        {
            order = k;
            break;
        }
    }

    if (k == 15)
    {
        m_status(HW_MATH_ERR_FILTERSPEC_E);
		return;
    }

    omegaC1 = omegaC / omegaP;
    omegaC2 = omegaC / omegaS * omegaS2;
}
//------------------------------------------------------------------------------
// Determine the prototype pass and stop corner frequencies for a given order
//------------------------------------------------------------------------------
void hwBessel_Design::SolveProtoCutoff(double  epsilon, 
                                       double  delta,
                                       int     order, 
                                       double& omegaC,
                                       double& omegaP, 
                                       double& omegaS)
{
    // generate H(s) denominator coefficients in descending order
    double factor;
    hwMatrix Hs_denom(order+1, hwMatrix::REAL);

    Hs_denom(0) = 1.0;

    for (int k = 0; k < order; ++k)
    {
        factor = (order+k+1)*(order-k)/(2.0*(k+1));	// factor = a(k+1) / a(k)
        Hs_denom(k+1) = factor * Hs_denom(k);
    }

    // generate H(s*) denominator coefficients
    hwMatrix Hs_conj(order+1, hwMatrix::REAL);

    for (int k = order; k > -1; k-=2)
    {
        Hs_conj(k) = Hs_denom(k);
    }

    for (int k = order-1; k > -1; k-=2)
    {
        Hs_conj(k) = -Hs_denom(k);
    }
    // generate |H(s)|^2 denominator coefficients
    hwMatrix Hs2_denom;

    m_status = Hs2_denom.ConvLin(Hs_denom, Hs_conj);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return;
    }

    // convert |H(s)|^2 denominator coefficients to |H(omega)|^2 coefficients
    for (int k = Hs2_denom.Size()-3; k > -1; k-=4)
    {
        Hs2_denom(k) = -Hs2_denom(k);
    }

    // convert |H(omega)|^2 denominator coefficients to |H(omega^2)|^2 coefficients
    hwMatrix Hs22_denom(order+1, hwMatrix::REAL);

    for (int k = 0; k <= order; ++k)
    {
        Hs22_denom(k) = Hs2_denom(2*k);
    }

    // solve |H(omega^2)|^2 = 1/2 to find angular pass frequency
    hwMatrix Roots;

    Hs22_denom(order) = -Hs22_denom(order);

    m_status = PolyRoots(Hs22_denom, Roots);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return;
    }

    for (int k = 0; k <= order; ++k)
    {
        // find the positive real root
        if (Roots.IsReal())
        {
            if (Roots(k) > 1.0e-8)
            {
                omegaC = sqrt(Roots(k));
                break;
            }
        }
        else
        {
            if (Roots.z(k).IsReal(1.0e-8))
            {
                if (Roots.z(k).Real() > 0.0)
                {
                    omegaC = sqrt(Roots.z(k).Real());
                    break;
                }
            }
        }
    }

    // solve |H(omega^2)|^2 = 1/(1+epsilon^2) to find angular pass cutoff frequency
    double epsilon2 = epsilon * epsilon;

    Hs22_denom(order) *= epsilon2;

    m_status = PolyRoots(Hs22_denom, Roots);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return;
    }

    for (int k = 0; k <= order; ++k)
    {
        // find the positive real root
        if (Roots.IsReal())
        {
            if (Roots(k) > 1.0e-8)
            {
                omegaP = sqrt(Roots(k));
                break;
            }
        }
        else
        {
            if (Roots.z(k).IsReal(1.0e-8))
            {
                if (Roots.z(k).Real() > 0.0)
                {
                    omegaP = sqrt(Roots.z(k).Real());
                    break;
                }
            }
        }
    }

    // solve |H(omega^2)|^2 = 1/(1+delta^2) to find angular stop frequency
    Hs22_denom(order) *= (delta * delta) / epsilon2;

    if (Hs22_denom(order) < -1.0e+12)
    {
        m_status(HW_MATH_ERR_FILTERSPEC_E);
		return;
    }

    m_status = PolyRoots(Hs22_denom, Roots);

    if (!m_status.IsOk())
    {
        m_status.ResetArgs();
        return;
    }

    for (int k = 0; k <= order; ++k)
    {
        // find the positive real root
        if (Roots.IsReal())
        {
            if (Roots(k) > 1.0e-8)
            {
                omegaS = sqrt(Roots(k));
                break;
            }
        }
        else
        {
            if (Roots.z(k).IsReal(1.0e-8))
            {
                if (Roots.z(k).Real() > 0.0)
                {
                    omegaS = sqrt(Roots.z(k).Real());
                    break;
                }
            }
        }
    }
}

