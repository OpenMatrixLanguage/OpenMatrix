/**
* @file hwLogNormal.cxx
* @date June 2009
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
#include "hwLogNormal.h"

#ifndef OS_WIN
    #include <math.h>
#endif

#include "Globals.h"
#include "hwNormal.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwLogNormal::hwLogNormal(double                  mu,
                         double                  sigma,
                         hwMersenneTwisterState* pMTState)
    : m_pNormDist (nullptr)
    , m_mu        (mu)
    , m_sigma     (sigma)
{
    if (sigma < 0.0)
    {
        return;
    }
    m_pNormDist = new hwNormal(pMTState);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwLogNormal::~hwLogNormal()
{
    delete m_pNormDist;
}
//------------------------------------------------------------------------------
// Returns probability density function
//------------------------------------------------------------------------------
double hwLogNormal::Pdf(double x)
{
    if (x < 0.0)
    {
        return 0.0;
    }

    double value = (log(x) - m_mu) / m_sigma;

    value *= 0.5 * value;

    return exp(-value) / (sqrt(2.0*PI) * m_sigma * x);
}
//------------------------------------------------------------------------------
// Returns cumulative density function
//------------------------------------------------------------------------------
double hwLogNormal::Cdf(double x)
{
    if (x < 0.0)
    {
        return 0.0;
    }

    double value = (log(x) - m_mu) / m_sigma;

    return m_pNormDist->Cdf(value);
}
//------------------------------------------------------------------------------
// Returns inverse cumulative density function
//------------------------------------------------------------------------------
double hwLogNormal::CdfInv(double prob)
{
    if (prob <= 0.0)
    {
        return 0.0;
    }
    if (prob >= 1.0)
    {
        return 99999.0 * Mean();
    }

    double x = m_pNormDist->CdfInv(prob);

    return exp(x * m_sigma + m_mu);
}
//------------------------------------------------------------------------------
// Compute the mean of the distribution
//------------------------------------------------------------------------------
double hwLogNormal::Mean()
{
    return exp(m_mu + 0.5 * m_sigma * m_sigma);
}
//------------------------------------------------------------------------------
// Compute the variance of the distribution
//------------------------------------------------------------------------------
double hwLogNormal::Variance()
{
    return (exp(m_sigma*m_sigma) - 1.0) * exp(2.0 * m_mu + m_sigma * m_sigma);
}
//------------------------------------------------------------------------------
// Compute the median of the distribution
//------------------------------------------------------------------------------
double hwLogNormal::Median()
{
    return exp(m_mu);
}
//------------------------------------------------------------------------------
// Compute the mode of the distribution
//------------------------------------------------------------------------------
double hwLogNormal::Mode()
{
    return exp(m_mu - m_sigma * m_sigma);
}
//------------------------------------------------------------------------------
// Compute a random deviate from the distribution
//------------------------------------------------------------------------------
double hwLogNormal::GetDeviate()
{
    double x = m_pNormDist->GetDeviate();

    return exp(x * m_sigma + m_mu);
}
