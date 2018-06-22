/**
* @file hwGammaInv.cxx
* @date May 2009
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
#include "hwGammaInv.h"

#include <math.h>

#include "Globals.h"
#include "hwGamma.h"
#include "SpecialFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwGammaInv::hwGammaInv(double alpha, hwMersenneTwisterState* pMTState)
    : m_pGamma (nullptr)
    , m_alpha  (alpha)
{
    if (alpha <= 0.0)
    {
        return;
    }

    m_pGamma = new hwGamma(m_alpha, pMTState);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwGammaInv::~hwGammaInv()
{
    delete m_pGamma;
}
//------------------------------------------------------------------------------
// Returns probability density function
//------------------------------------------------------------------------------
double hwGammaInv::Pdf(double x)
{
    if (x < 0.0)
    {
        return 0.0;
    }

    double gamma;
    if (m_alpha < MAXGAM)
    {
        gamma = GammaFunc(m_alpha);
    }
    else
    {
        gamma = GammaLog(m_alpha);
        gamma = (gamma < MINLOG) ? 0.0 : exp(gamma);
    }

    double nearZero = 1.0e-12;
    if (fabs(gamma) < nearZero)
    {
        return 0.0;
    }
    return pow(x, -(m_alpha + 1.0)) * exp(-1.0 / x) / gamma;
}
//------------------------------------------------------------------------------
// Returns cumulative density function
//------------------------------------------------------------------------------
double hwGammaInv::Cdf(double x)
{
    double nearZero = 1.0e-12;

    if (x < nearZero)
    {
        return 0.0;
    }
    return 1.0 - m_pGamma->Cdf(1.0 / x);
}
//------------------------------------------------------------------------------
// Returns inverse cumulative density function
//------------------------------------------------------------------------------
double hwGammaInv::CdfInv(double prob)
{
    double nearZero = 1.0e-12;
    double value    = m_pGamma->CdfInv(1.0 - prob);

    if (value < nearZero)
    {
        return 0.0;
    }

    return 1.0 / value;
}
//------------------------------------------------------------------------------
// Returns distribution mean
//------------------------------------------------------------------------------
double hwGammaInv::Mean()
{
    return 1.0 / (m_alpha - 1.0);
}
//------------------------------------------------------------------------------
// Returns distribution variance
//------------------------------------------------------------------------------
double hwGammaInv::Variance()
{
    return 1.0 / ((m_alpha - 1.0) * (m_alpha - 1.0) * (m_alpha - 2.0));
}
//------------------------------------------------------------------------------
// Returns distribution mode
//------------------------------------------------------------------------------
double hwGammaInv::Mode()
{
    return 1.0 / (m_alpha + 1.0);
}
//------------------------------------------------------------------------------
// Returns random value
//------------------------------------------------------------------------------
double hwGammaInv::GetDeviate()
{
    double nearZero = 1.0e-12;
    double g        = m_pGamma->GetDeviate();

    if (g < nearZero)
    {
        return 0.0;
    }
    return 1.0 / g;
}
