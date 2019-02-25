/**
* @file hw_F.cxx
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
#include "hw_F.h"

#include "hwBetaInv.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hw_F::hw_F(int m, int n, hwMersenneTwisterState* pMTState)
    : m_pBetaInv (nullptr)
{
	if (m < 1 || n < 1)
	{
		return;
	}
	
	m_scale = (double) n / (double) m;
	
	m_pBetaInv = new hwBetaInv(0.5 * m, 0.5 * n, pMTState);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hw_F::~hw_F()
{
    delete m_pBetaInv;
}
//------------------------------------------------------------------------------
// Returns probability density function
//------------------------------------------------------------------------------
double hw_F::Pdf(double x)
{
    return m_pBetaInv->Pdf(x / m_scale) / m_scale;
}
//------------------------------------------------------------------------------
// Returns cumulative density function
//------------------------------------------------------------------------------
double hw_F::Cdf(double x)
{
    return m_pBetaInv->Cdf(x / m_scale);
}
//------------------------------------------------------------------------------
// Returns inverse cumulative density function
//------------------------------------------------------------------------------
double hw_F::CdfInv(double prob)
{
    return m_scale * m_pBetaInv->CdfInv(prob);
}
//------------------------------------------------------------------------------
// Compute the mean of the distribution
//------------------------------------------------------------------------------
double hw_F::Mean()
{
    return m_scale * m_pBetaInv->Mean();
}
//------------------------------------------------------------------------------
// Returns distribution variance
//------------------------------------------------------------------------------
double hw_F::Variance()
{
    return m_scale * m_scale * m_pBetaInv->Variance();
}
//------------------------------------------------------------------------------
// Returns distribution mode
//------------------------------------------------------------------------------
double hw_F::Mode()
{
    return m_scale * m_pBetaInv->Mode();
}
//------------------------------------------------------------------------------
// Returns random value
//------------------------------------------------------------------------------
double hw_F::GetDeviate()
{
    return m_scale * m_pBetaInv->GetDeviate();
}
