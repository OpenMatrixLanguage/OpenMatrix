/**
* @file hwWeibull.cxx
* @date October 2011
* Copyright (C) 2011-2018 Altair Engineering, Inc.  
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
#include "hwWeibull.h"

#include <math.h>

#include "Globals.h"
#include "hwUniform.h"
#include "SpecialFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwWeibull::hwWeibull(double                  alpha,
                     double                  beta,
                     hwMersenneTwisterState* pMTState)
    : m_pUniform (nullptr)
    , m_alpha    (alpha)
    , m_beta     (beta)
{
    if (alpha <= 0.0)
    {
        return;
    }

    if (beta <= 0.0)
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
hwWeibull::~hwWeibull()
{
    delete m_pUniform;
}
//------------------------------------------------------------------------------
// Probability density function
//------------------------------------------------------------------------------
double hwWeibull::Pdf(double x)
{
    if (x < 0.0)
    {
        return 0.0;
    }

	return m_beta/pow(m_alpha, m_beta) * pow(x, m_beta-1) * exp(-pow(x/m_alpha, m_beta));
}
//------------------------------------------------------------------------------
// Returns cumulative density function
//------------------------------------------------------------------------------
double hwWeibull::Cdf(double x)
{
    if (x < 0.0)
    {
        return 0.0;
    }
	return 1.0 - exp(-pow(x/m_alpha, m_beta));
}

//------------------------------------------------------------------------------
// Returns inverse cumulative density function
//------------------------------------------------------------------------------
double hwWeibull::CdfInv(double prob)
{
    if (prob >= 1.0)
    {
        return MAXNUM;  // should change to an error
    }
    return m_alpha * pow(-log(1.0-prob), 1.0/m_beta);
}
//------------------------------------------------------------------------------
// Compute the mean of the distribution
//------------------------------------------------------------------------------
double hwWeibull::Mean()
{
	return m_alpha * GammaFunc(1.0 + 1.0 / m_beta);
}
//------------------------------------------------------------------------------
// Compute the variance of the distribution
//------------------------------------------------------------------------------
double hwWeibull::Variance()
{
    double temp1 = GammaFunc(1.0 + 1.0 / m_beta);
    double temp2 = GammaFunc(1.0 + 2.0 / m_beta);

    return m_alpha * m_alpha * (temp2 - temp1 * temp1);
}
//------------------------------------------------------------------------------
// Compute the median of the distribution
//------------------------------------------------------------------------------
double hwWeibull::Median()
{
    return m_alpha * pow(log(2.0), 1.0/m_beta);
}
//------------------------------------------------------------------------------
// Compute the mode of the distribution
//------------------------------------------------------------------------------
double hwWeibull::Mode()
{
    return m_alpha * pow(log((m_beta-1.0)/m_beta), 1.0/m_beta);
}
//------------------------------------------------------------------------------
// Compute a random deviate from the distribution
//------------------------------------------------------------------------------
double hwWeibull::GetDeviate()
{
    double u = m_pUniform->GetDeviate();

	return m_alpha * pow(-log(u), 1.0/m_beta);
}
