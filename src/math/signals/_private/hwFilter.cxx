/**
* @file  hwFilter.cxx
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
#include "hwFilter.h"

#include "hwMatrix.h"
#include "GeneralFuncs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwFilter::hwFilter()
    : m_pNumerCoef(nullptr)
    , m_pDenomCoef(nullptr)
{
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwFilter::hwFilter(const hwMatrix& pNumerCoef, const hwMatrix* pDenomCoef)
    : m_pNumerCoef(nullptr)
    , m_pDenomCoef(nullptr)
{
    if (!pNumerCoef.IsReal())
    {
        m_status(HW_MATH_ERR_COMPLEX, 1);
        return;
    }

    if (!pNumerCoef.IsVector())
    {
        if (pNumerCoef.IsEmpty())
        {
            m_status(HW_MATH_ERR_EMPTYMATRIX, 1);
        }
        else
        {
            m_status(HW_MATH_ERR_VECTOR, 1);
        }
        return;
    }

    if (pDenomCoef)
    {
        if (!pDenomCoef->IsReal())
        {
            m_status(HW_MATH_ERR_COMPLEX, 2);
            return;
        }

        if (!pDenomCoef->IsVector())
        {
            m_status(HW_MATH_ERR_VECTOR, 2);
            return;
        }

        m_pDenomCoef = new hwMatrix(*pDenomCoef);
    }
    m_pNumerCoef = new hwMatrix(pNumerCoef);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwFilter::~hwFilter()
{
    if (m_pNumerCoef)
    {
        delete m_pNumerCoef;
    }
    if (m_pDenomCoef)
    {
        delete m_pDenomCoef;
    }
}
//------------------------------------------------------------------------------
// Set length of transfer function polynomials
//------------------------------------------------------------------------------
void hwFilter::SetSize(int size, bool IIR)
{
    if (m_pNumerCoef)
    {
        m_status = m_pNumerCoef->Resize(size);
        if (!m_status.IsOk())
        {
            return;
        }
    }
    else
    {
        m_pNumerCoef = new hwMatrix(size, hwMatrix::REAL);
    }

    if (IIR)
    {
        if (m_pDenomCoef)
        {
            m_status = m_pDenomCoef->Resize(size);
        }
        else
        {
            m_pDenomCoef = new hwMatrix(size, hwMatrix::REAL);
        }
    }
    else
    {
        if (m_pDenomCoef)
        {
            delete m_pDenomCoef;
            m_pDenomCoef = NULL;
        }
    }
}
//------------------------------------------------------------------------------
// Compute the response magnitude
//------------------------------------------------------------------------------
void hwFilter::RespMag(double omega, double& mag)
{
    hwComplex response;
    Response(omega, response);

    if (!m_status.IsOk())
    {
        if (m_status != HW_MATH_WARN_NYQUIST)
        {
            return;
        }
    }

    mag = response.Mag();
}
//------------------------------------------------------------------------------
// Compute the response phase
//------------------------------------------------------------------------------
void hwFilter::RespPhase(double omega, double& phase)
{
    hwComplex response;
    
    Response(omega, response);

    if (!m_status.IsOk())
    {
        if (m_status != HW_MATH_WARN_NYQUIST)
        {
            return;
        }
    }

    if (!IsZero(response.Mag(), 1.0e-12))
    {
        phase = response.Arg();
    }
    else
    {
        phase = 0.0;
    }
}
