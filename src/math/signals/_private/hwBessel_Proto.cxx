/**
* @file hwBessel_Proto.cxx
* @date May 2012
* Copyright (C) 2012-2018 Altair Engineering, Inc.  
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
#include "hwBessel_Proto.h"

#include "PolynomFuncs.h"

#ifndef OS_WIN
  #include <math.h>
#endif

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwBessel_Proto::hwBessel_Proto(int         order,
                               const char* type)
	: hwLowPass_Proto(order)
{
    // generate H(s) denominator coefficients in descending order
    double factor;
    hwMatrix Hs_denom(m_order+1, hwMatrix::REAL);

    Hs_denom(0) = 1.0;

    for (int k = 0; k < m_order; ++k)
    {
        factor = (m_order+k+1)*(m_order-k)/(2.0*(k+1));	// factor = a(k+1) / a(k)
        Hs_denom(k+1) = factor * Hs_denom(k);
    }

    // compute scale factor    
    if (!strcmp(type, "besself"))
    {
        // normalize for the Butterworth asymptotic form
        m_scale = pow(Hs_denom(m_order), -1.0/m_order);
    }
    else if (!strcmp(type, "besself3"))
    {
        // normalize to set the cutoff frequency attenuation at -3dB
        // generate H(s*) denominator coefficients
        hwMatrix Hs_conj(m_order+1, hwMatrix::REAL);

        for (int k = m_order; k > -1; k-=2)
        {
            Hs_conj(k) = Hs_denom(k);
        }

        for (int k = m_order-1; k > -1; k-=2)
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
        hwMatrix Hs22_denom(m_order+1, hwMatrix::REAL);

        for (int k = 0; k <= m_order; ++k)
        {
            Hs22_denom(k) = Hs2_denom(2*k);
        }

        if (Hs22_denom(m_order) > 1.0e+12)
        {
            m_status(HW_MATH_WARN_FILTERSPEC_W);
        }

        // solve |H(omega^2)|^2 = 1/2 to find angular cutoff frequency
        Hs22_denom(m_order) = -Hs22_denom(m_order);

        hwMathStatus status = PolyRoots(Hs22_denom, m_roots);
        if (!status.IsOk())
        {
            m_status = status;
            m_status.ResetArgs();
            return;
        }

        for (int k = 0; k <= m_order; ++k)
        {
            // find the positive real root
            if (m_roots.IsReal())
            {
                if (m_roots(k) > 1.0e-8)
                {
                    m_scale = 1.0 / sqrt(m_roots(k));
                    break;
                }
            }
            else
            {
                if (m_roots.z(k).IsReal(1.0e-8))
                {
                    if (m_roots.z(k).Real() > 0.0)
                    {
                        m_scale = 1.0 / sqrt(m_roots.z(k).Real());
                        break;
                    }
                }
            }
        }
    }
    else
    {
        m_status(HW_MATH_ERR_FILTERTYPE, 2);
        return;
    }

    // get poles of H(s) denominator
    hwMathStatus status = PolyRoots(Hs_denom, m_roots);
    if (!status.IsOk())
    {
        m_status = status;
        m_status.ResetArgs();
        return;
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwBessel_Proto::~hwBessel_Proto()
{
}
//------------------------------------------------------------------------------
// Compute location of the real pole
//------------------------------------------------------------------------------
void hwBessel_Proto::GetSPlaneInfo(double& poleReal) const
{
    // real pole
    for (int k = 0; k < m_order; ++k)
    {
        // find the real pole in the left half plane
        if (m_roots.IsReal())
        {
            if (m_roots(k) < 1.0e-8)
            {
                poleReal = m_roots(k) * m_scale;
                break;
            }
        }
        else
        {
            if (m_roots.z(k).IsReal(1.0e-8))
            {
                if (m_roots.z(k).Real() < 0.0)
                {
                    poleReal = m_roots.z(k).Real() * m_scale;
                    break;
                }
            }
        }
    }
}
//------------------------------------------------------------------------------
// Compute real component and squared magnitude of the ith pole
//------------------------------------------------------------------------------
void hwBessel_Proto::GetSPlaneInfo(int     i,
                                   double& poleReal,
                                   double& poleMagSq) const
{
    // complex conjugate pole pair
    int count = 0;
    double real;
    double imag;

    for (int k = 0; k < m_order; ++k)
    {
        // find the ith left half plane pole with a positive imaginary component
        imag = m_roots.z(k).Imag();

        if (imag > 1.0e-8)
        {
            real = m_roots.z(k).Real();

            if (real < 0.0)
            {
                if (count == i)
                {
                    real *= m_scale;
                    imag *= m_scale;
                    poleReal = real;
                    poleMagSq = real*real + imag*imag;
                    break;
                }

                ++count;
            }
        }
    }
}
