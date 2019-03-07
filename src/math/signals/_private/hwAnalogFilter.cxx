/**
* @file hwAnalogFilter.cxx
* @date April 2009
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
#include "hwAnalogFilter.h"

#include "hwMatrix.h"
#include "GeneralFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwAnalogFilter::hwAnalogFilter()
{
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwAnalogFilter::hwAnalogFilter(const hwMatrix& pNumerCoef, 
                               const hwMatrix* pDenomCoef)
    : hwFilter(pNumerCoef, pDenomCoef)
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwAnalogFilter::~hwAnalogFilter()
{
}
//------------------------------------------------------------------------------
// Evaluate a polynomial in j*omega
//------------------------------------------------------------------------------
static void EvalPoly_s(const hwMatrix* poly,
                       double          omega,
                       hwComplex&      polyVal)
{
    int       k;
    int       size   = poly->Size();
    double    omega2 = omega * omega;
    double    omega4 = omega2 * omega2;
    double    sum[4];
    double*   pCoef;
    hwMatrix* P = (hwMatrix*) poly;     // overreide const
    double* pStart = P->GetRealData();

    // separate the sum into four components based on powers of j
    // each component is a polynomial of omega^4
    // temporarily ignore remaining factors of omega
    for (int i = 0; i < 4; ++i)
    {
        if (i < size)
        {
            pCoef = pStart + i;
            sum[i] = *pCoef;
            k = size - i - 4;

            while (k > 0)
            {
                pCoef += 4;
                sum[i] = sum[i] * omega4 + *pCoef;
                k -= 4;
            }
        }
        else
            break;
    }

    // account for remaining factors of omega and combine
    // components to form result
    if (size == 1)
    {
        polyVal.Set(sum[0], 0.0);
    }
    else if (size == 2)
    {
        polyVal.Set(sum[1], sum[0] * omega);
    }
    else if (size == 3)
    {
        polyVal.Set(sum[2] - sum[0] * omega2, sum[1] * omega);
    }
    else if (size % 4 == 0)
    {
        polyVal.Set(sum[3] - sum[1] * omega2, sum[2] * omega - sum[0] * omega2 * omega);
    }
    else if (size % 4 == 1)
    {
        polyVal.Set(sum[0] - sum[2] * omega2, sum[3] * omega - sum[1] * omega2 * omega);
    }
    else if (size % 4 == 2)
    {
        polyVal.Set(sum[1] - sum[3] * omega2, sum[0] * omega - sum[2] * omega2 * omega);
    }
    else 
    {
        polyVal.Set(sum[2] - sum[0] * omega2, sum[1] * omega - sum[3] * omega2 * omega);
    }
}
//------------------------------------------------------------------------------
// Compute the filter response at angular frequency omega
//------------------------------------------------------------------------------
void hwAnalogFilter::Response(double omega, hwComplex& resp)
{
    if (!m_pNumerCoef)
    {
        resp.Set(0.0, 0.0); 
        m_status(HW_MATH_ERR_EMPTYMATRIX, 0);
        return;
    }

    // evaluate the transfer function where s = j*omega
    hwComplex numer;
    EvalPoly_s(m_pNumerCoef, omega, numer);

    if (m_pDenomCoef)
    {
        hwComplex denom;

        EvalPoly_s(m_pDenomCoef, omega, denom);
        resp = numer / denom;
    }
    else
    {
        resp = numer;
    }
}
