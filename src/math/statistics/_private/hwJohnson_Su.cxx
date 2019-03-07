/**
* @file hwJohnson_Su.cxx
* @date May 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
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
#include "hwJohnson_Su.h"

#include <math.h>

#include "GeneralFuncs.h"
#include "hwNormal.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwJohnson_Su::hwJohnson_Su(double                  beta1, 
                           double                  beta2,
                           hwMersenneTwisterState* pMTState)
    : m_pNormal (nullptr)
{
    double lowerCaseOmega;
    double upperCaseOmega;

    ComputeOmegaValues(beta1, beta2, lowerCaseOmega, upperCaseOmega);

    m_delta = 1.0 / sqrt(log(lowerCaseOmega));
    m_gamma = upperCaseOmega * m_delta;

    m_mean = -sqrt(lowerCaseOmega) * sinh(upperCaseOmega);
    m_variance = 0.5*(lowerCaseOmega-1.0)*(lowerCaseOmega*cosh(2.0*upperCaseOmega)+1.0);

    m_pNormal = new hwNormal(pMTState);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwJohnson_Su::~hwJohnson_Su()
{
    delete m_pNormal;
}
//------------------------------------------------------------------------------
// Compute the omega parameters of the distribution
//------------------------------------------------------------------------------
void hwJohnson_Su::ComputeOmegaValues(double  beta1,
                                      double  beta2,
                                      double& lowerCaseOmega,
                                      double& upperCaseOmega)
{
    double A0;
    double A1;
    double A2;
    double A;
    double B;
    double C;
    double beta1_approx;
    double m;
    double omega2;
    double omega3; 
    double omega4;
    double omega5;
    double value1;
    double value2;
    double value4; 
    double value5;
    double zero1;
    double zero2;

    double omega = sqrt(sqrt(2.0 * beta2 - 2.8 * beta1 - 2.0) - 1.0); // = exp(1.0 / (delta * delta));
    double value3 = 2.0 * (beta2 - 3.0);

    while (true)
    {
        omega2 = omega * omega;
        omega3 = omega2 * omega;
        omega4 = omega2 * omega2;
        omega5 = omega3 * omega2;

        A0 = omega5 + 3.0 * omega4 + 6.0 * omega3 + 10.0 * omega2 + 9.0 * omega + 3;
        A1 = 8.0 * (omega4 + 3.0 * omega3 + 6.0 * omega2 + 7.0 * omega + 3.0);
        A2 = 8.0 * (omega3 + 3.0 * omega2 + 6.0 * omega + 6.0);

        value1 = omega - 1.0;
        value2 = omega + 1.0;
        A = A2 * value1 - 4.0 * value3;
        B = A1 * value1 - 4.0 * value3 * value2;
        C = A0 * value1 - value3 * value2 * value2;

        quadraticRoots(A, B, C, zero1, zero2);

        m = (zero1 > zero2) ? zero1 : zero2;

        value4 = 4.0 * (omega + 2.0) * m + 3.0 * value2 * value2;
        value5 = 2.0 * m + value2;
        value5 *= value5 * value5;
        beta1_approx = 0.5 * m * value1 * value4 * value4 / value5;

        value1 = beta2 - 0.5 * (omega4 + 2.0 * omega2 + 3.0);
        value2 = beta1 / beta1_approx * value1;
        A = 0.5;
        B = 1.0;
        C = value2 - beta2 + 1.5;

        quadraticRoots(A, B, C, zero1, zero2);

        omega = (zero1 > zero2) ? sqrt(zero1) : sqrt(zero2);

        if (fabs(beta1 - beta1_approx) < 1.0e-6)
        {
            break;
        }
    }

    lowerCaseOmega = omega;
    upperCaseOmega = log(sqrt(m/omega) + sqrt(m/omega + 1.0));
}
//------------------------------------------------------------------------------
// Returns probability density function
//------------------------------------------------------------------------------
double hwJohnson_Su::Pdf(double x)
{
    double value1 = sqrt(x * x + 1.0);
    double value2 = m_gamma + m_delta * log(x + value1);
    double value3 = m_pNormal->Pdf(value2);

    return m_delta / value1 * value3;
}
//------------------------------------------------------------------------------
// Returns cumulative density function
//------------------------------------------------------------------------------
double hwJohnson_Su::Cdf(double x)
{
    double value = m_gamma + m_delta * log(x + sqrt(x * x + 1.0));

    return m_pNormal->Cdf(value);
}
//------------------------------------------------------------------------------
// Returns inverse cumulative density function
//------------------------------------------------------------------------------
double hwJohnson_Su::CdfInv(double prob)
{
    double z = m_pNormal->CdfInv(prob);

    return sinh((z - m_gamma) / m_delta);
}
//------------------------------------------------------------------------------
// Compute the mean of the distribution
//------------------------------------------------------------------------------
double hwJohnson_Su::Mean()
{
    return m_mean;
}
//------------------------------------------------------------------------------
// Compute the variance of the distribution
//------------------------------------------------------------------------------
double hwJohnson_Su::Variance()
{
    return m_variance;
}
//------------------------------------------------------------------------------
// Compute the mode of the distribution
//------------------------------------------------------------------------------
double hwJohnson_Su::Mode()
{
    double x = sinh(-m_gamma / m_delta);
    double f;
    double fp;
    double value1;
    double value2;

    for (int i = 0; i < 20; i++)
    {
        value1 = x * x + 1.0;
        value2 = sqrt(value1);

        f = x / value2 + m_delta * (m_gamma + m_delta * log(x + value2));

        if (fabs(f) < 1.0e-10)
        {
            break;          // function convergence
        }

        fp = (1.0 + m_delta * m_delta * value1) / (value1 * value2);

        if (fabs(f) < 1.0e-10 * fabs(fp))
        {
            break;          // step convergence
        }
        x -= f / fp;
    }

    return x;
}
//------------------------------------------------------------------------------
// Compute a random deviate from the distribution
//------------------------------------------------------------------------------
double hwJohnson_Su::GetDeviate()
{
    double z = m_pNormal->GetDeviate();

    return sinh((z - m_gamma) / m_delta);
}
