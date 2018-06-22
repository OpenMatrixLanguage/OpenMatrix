/**
* @file hwPearson_VII.cxx
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
#include "hwPearson_VII.h"

#include <math.h>

#include "hwBetaInv.h"
#include "hwUniform.h"
#include "SpecialFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwPearson_VII::hwPearson_VII(double shape, hwMersenneTwisterState* pMTState)
    : m_pBetaInv (nullptr)
    , m_pUniform (nullptr)
    , m_shape    (shape)
{
    if (shape <= 0.5)
    {
        return;
    }

    m_ndof = 2.0 * m_shape - 1.0;

    m_pBetaInv = new hwBetaInv(0.5, m_shape - 0.5, pMTState);

    if (pMTState)
    {
        m_pUniform = new hwUniform(pMTState);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwPearson_VII::~hwPearson_VII()
{
    delete m_pBetaInv;
    delete m_pUniform;
}
//------------------------------------------------------------------------------
// Probability density function
//------------------------------------------------------------------------------
double hwPearson_VII::Pdf(double x)
{
    double value1 = BetaFunc(0.5, m_shape - 0.5) * sqrt(m_ndof);
    double value2 = pow(1.0 + x * x / m_ndof, -m_shape);

    return value2 / value1;
}
//------------------------------------------------------------------------------
// Cummulative density function
//------------------------------------------------------------------------------
double hwPearson_VII::Cdf(double x)
{
    if (x < 0.0)
    {
        return 0.5 * (1.0 - m_pBetaInv->Cdf(x * x / m_ndof));
    }

    return 0.5 * (1.0 + m_pBetaInv->Cdf(x * x / m_ndof));
}
//------------------------------------------------------------------------------
// Inverse cummulative density function
//------------------------------------------------------------------------------
double hwPearson_VII::CdfInv(double prob)
{
    if (prob < 0.5)
    {
        return -sqrt(m_pBetaInv->CdfInv(1.0 - 2.0 * prob) * m_ndof);
    }
    return sqrt(m_pBetaInv->CdfInv(2.0 * prob - 1.0) * m_ndof);
}
//------------------------------------------------------------------------------
// Compute the variance of the distribution
//------------------------------------------------------------------------------
double hwPearson_VII::Variance()
{
    return m_ndof / (m_ndof - 2.0);
}
//------------------------------------------------------------------------------
// Compute a random deviate from the distribution
//------------------------------------------------------------------------------
double hwPearson_VII::GetDeviate()
{
    double u = m_pUniform->GetDeviate();
    double b = sqrt(m_ndof * m_pBetaInv->GetDeviate());

    if (u > 0.5)
    {
        return b;
    }

    return -b;
}
