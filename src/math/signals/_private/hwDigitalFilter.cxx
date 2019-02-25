/**
* @file  hwDigitalFilter.cxx
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
#include "hwDigitalFilter.h"

#include "hwMatrix.h"
#include "GeneralFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwDigitalFilter::hwDigitalFilter()
{
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwDigitalFilter::hwDigitalFilter(const hwMatrix& pNumerCoef, 
                                 const hwMatrix* pDenomCoef)
    : hwFilter(pNumerCoef, pDenomCoef)
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwDigitalFilter::~hwDigitalFilter()
{
}
//------------------------------------------------------------------------------
// Compute the filter response at angular frequency omega
//------------------------------------------------------------------------------
void hwDigitalFilter::Response(double omega, hwComplex& resp)
{
    if (!m_pNumerCoef)
    {
        resp.Set(0.0, 0.0); 
        m_status(HW_MATH_ERR_EMPTYMATRIX, 0);
        return;
    }

    // evaluate the transfer function where z = e^(j*omega)
    hwComplex numer (0, 0);
    hwComplex value;
    double*   numerCoef = m_pNumerCoef->GetRealData();

    for (int i = 0; i < m_pNumerCoef->Size(); ++i)
    {
        value.Set(cos(-omega * i), sin(-omega * i));
        numer += value * numerCoef[i];
    }

    if (m_pDenomCoef)
    {
        hwComplex denom (0, 0);
        double*   denomCoef = m_pDenomCoef->GetRealData();

        for (int i = 0; i < m_pDenomCoef->Size(); ++i)
        {
            value.Set(cos(-omega * i), sin(-omega * i));
            denom += value * denomCoef[i];
        }

        resp = numer / denom;
    }
    else
    {
        resp = numer;
    }

    if (omega < -PI || omega > PI)
    {
        m_status(HW_MATH_WARN_NYQUIST, 1);
    }
}
