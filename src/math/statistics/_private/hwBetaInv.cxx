/**
* @file hwBetaInv.cxx
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
#include "hwBetaInv.h"

#include <math.h>

#include "hwBeta.h"
#include "SpecialFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwBetaInv::hwBetaInv(double alpha, double beta, hwMersenneTwisterState* pMTState)
    : m_pBeta (nullptr)
    , m_alpha (alpha)
    , m_beta  (beta)
{
    if (alpha <= 0.0 || beta <= 0.0)
    {
        return;
    }
    m_pBeta = new hwBeta(m_alpha, m_beta, pMTState);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwBetaInv::~hwBetaInv()
{
    delete m_pBeta;
}
//------------------------------------------------------------------------------
// Returns probability density function
//------------------------------------------------------------------------------
double hwBetaInv::Pdf(double x)
{
    double nearZero = 1.0e-12;
    if (x < nearZero)
    {
        return 0.0;
    }

    double beta    = BetaFunc(m_alpha, m_beta);
    double density = 0.0;
    if (fabs(beta) > nearZero)
    {
        density = pow(x, m_alpha - 1.0) * pow(x + 1.0, -(m_alpha + m_beta)) / beta;
    }
    return density;
}
//------------------------------------------------------------------------------
// Returns cumulative density function
//------------------------------------------------------------------------------
double hwBetaInv::Cdf(double x)
{
    return m_pBeta->Cdf(x / (x + 1.0));
}
//------------------------------------------------------------------------------
// Returns inverse cumulative density function
//------------------------------------------------------------------------------
double hwBetaInv::CdfInv(double prob)
{
    double nearZero = 1.0e-12;
    double value    = m_pBeta->CdfInv(prob);

    if (1.0 - value < nearZero)
    {
        return 0.0;
    }

    return value / (1.0 - value);
}
//------------------------------------------------------------------------------
// Returns distribution mean
//------------------------------------------------------------------------------
double hwBetaInv::Mean()
{
    return m_alpha / (m_beta - 1.0);
}
//------------------------------------------------------------------------------
// Returns distribution variance
//------------------------------------------------------------------------------
double hwBetaInv::Variance()
{
    double value1 = m_alpha * (m_alpha + m_beta - 1.0);
    double value2 = (m_beta - 1.0) * (m_beta - 1.0) * (m_beta - 2.0);

    return value1 / value2;
}
//------------------------------------------------------------------------------
// Returns distribution mode
//------------------------------------------------------------------------------
double hwBetaInv::Mode()
{
    return (m_alpha - 1.0) / (m_beta + 1.0);
}
//------------------------------------------------------------------------------
// Returns random value
//------------------------------------------------------------------------------
double hwBetaInv::GetDeviate()
{
    double nearZero = 1.0e-12;
    double x        = m_pBeta->GetDeviate();

    if (x > 1.0 - nearZero)
    {
        return 0.0;
    }
    return x / (1.0 - x);
}
