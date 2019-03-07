/**
* @file hwPoisson.cxx
* @date June 2009
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
#include "hwPoisson.h"

#include "hwMatrix.h"
#include "hwUniform.h"
#include "SpecialFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwPoisson::hwPoisson(double lambda, hwMersenneTwisterState* pMTState)
    : m_pUniform (nullptr)
    , m_lambda   (lambda)
{
    if (lambda <= 0.0)
    {
        return;
    }

    if (pMTState)
    {
        m_pUniform = new hwUniform(pMTState);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwPoisson::~hwPoisson()
{
}
//------------------------------------------------------------------------------
// Compute the log of the factorial of the argument
//------------------------------------------------------------------------------
static double LogFactorial(double k)
{
    double logFact;

    if (k < 10.0)
    {
        double factorial;

        switch ((int) k)
        {
            case 0:  factorial = 1.0;      break;
            case 1:  factorial = 1.0;      break;
            case 2:  factorial = 2.0;      break;
            case 3:  factorial = 6.0;      break;
            case 4:  factorial = 24.0;     break;
            case 5:  factorial = 120.0;    break;
            case 6:  factorial = 720.0;    break;
            case 7:  factorial = 5040.0;   break;
            case 8:  factorial = 40320.0;  break;
            case 9:  factorial = 362880.0; break;
            default: break;
        }
        logFact = log(factorial);
    }
    else
    {
        logFact = GammaLog((double) k + 1.0);
    }
    return logFact;
}
//------------------------------------------------------------------------------
// Probability density function
//------------------------------------------------------------------------------
hwMathStatus hwPoisson::Pdf(double k, double& density)
{
    if (k < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 1);
    }
    double gammaLog = LogFactorial(k);

    density = exp(k * log(m_lambda) - m_lambda - gammaLog);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Returns cumulative density function
//------------------------------------------------------------------------------
double hwPoisson::Cdf(double k)
{
    if (k <= 0.0)
    {
        return 0.0;
    }

    double pdf;
    double cdf;

    Pdf(0, cdf);

    for (int i = 1; i <= static_cast<int>(k); i++)
    {
        Pdf(i, pdf);
        cdf += pdf;
    }

    return cdf;
}
//------------------------------------------------------------------------------
// Returns cumulative density function values from 0 to k
//------------------------------------------------------------------------------
hwMathStatus hwPoisson::Cdf(int k, hwMatrix& cdf)
{
    if (k < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }

    // compute vector of CDF values from 0 to k
    hwMathStatus status = cdf.Dimension(k + 1, hwMatrix::REAL);
    if (!status.IsOk())
    {
        status.SetArg1(2);
        return status;
    }

    Pdf(0, cdf(0));

    double pdf;
    for (int i = 1; i <= k; i++)
    {
        Pdf(i, pdf);
        cdf(i) = cdf(i-1) + pdf;
    }

    return status;
}
//------------------------------------------------------------------------------
// Returns inverse cumulative density function
//------------------------------------------------------------------------------
hwMathStatus hwPoisson::CdfInv(double prob, int& k)
{
    if (prob < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }
    if (prob > 0.999999)
    {
        return hwMathStatus(HW_MATH_ERR_BADRANGE, 1);
    }

    double pdf;
    double cdf;

    k = 0;
    Pdf(0, cdf);

    while (cdf < prob)
    {
        k++;
        Pdf(k, pdf);
        cdf += pdf;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Returns inverse cumulative density function values from 0 to k
//------------------------------------------------------------------------------
hwMathStatus hwPoisson::CdfInv(double prob, const hwMatrix& cdf, int& k)
{
    // compute inverse of CDF using vector of CDF values from 0 to k
    if (prob < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }
    if (prob > 1.0)
    {
        return hwMathStatus(HW_MATH_ERR_BADRANGE, 1);
    }

    int i    = 0;
    int size = cdf.Size();

    for (i = 0; i < size; i++)
    {
        if (prob > cdf(i))
        {
            k = i;
            break;
        }
    }

    if (i == size)
    {
        return hwMathStatus(HW_MATH_ERR_BADRANGE, 1);
    }
    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Compute a random deviate from the distribution
//------------------------------------------------------------------------------
int hwPoisson::GetDeviate()
{
    double g;       // deviate

    if (m_lambda < 10.0)
    {
        // accumulate exponential waiting times
        double p = 1.0;   // product
        double thresh = exp(-m_lambda);

        g = -1.0;

        do
        {
            g += 1.0;
            p *= m_pUniform->GetDeviate();
        }
        while (p > thresh);
    }
    else
    {
        // transformed rejection method with decomposition (Wolfgang Hormann)
        double c = m_lambda + 0.445;
        double b = 0.931 + 2.53 * sqrt(m_lambda);
        double a = -0.059 + 0.02483 * b;
        double alpha_r = 1.1239 + 1.1328 / (b - 3.4);
        double ur = 0.43;
        double vr = 0.9277 - 3.6224 / (b - 2.0);
        double u;
        double v;
        double lhs;
        double rhs;
        double value;

        do
        {
            v = m_pUniform->GetDeviate();

            if (v <= 2.0 * ur * vr)
            {
                u = v / vr - ur;
                value = 0.5 - fabs(u);
                g = (2.0 * a / value + b) * u + c;
                break;
            }

            if (v < vr)
            {
                u = v / vr - (ur + 0.5);

                if (u < 0.0)
                {
                    u = -0.5 - u;
                }
                else
                {
                    u = 0.5 - u;
                }

                v = vr * m_pUniform->GetDeviate();
            }
            else
            {
                u = m_pUniform->GetDeviate() - 0.5;
            }
            value = 0.5 - fabs(u);
            lhs = log(v * alpha_r / (a / (value*value) + b));
            g = (2.0 * a / value + b) * u + c;

            if (g > 0.0)
            {
                rhs = -m_lambda + floor(g) * log(m_lambda) - LogFactorial(floor(g));
            }
            else    // force loop to continue
            {
                lhs = 1.0;
                rhs = 0.0;
            }
        }
        while (lhs > rhs);
    }

    return static_cast<int> (floor(g));
}
