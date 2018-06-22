/**
* @file  hwDigitalFilterGen_IIR.cxx
* @date April 2009
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
#include "hwDigitalFilterGen_IIR.h"

#include "hwMatrix.h"
#include "hwMathStatus.h"
#include "hwDigitalFilter.h"
#include "hwFilterSpecs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwDigitalFilterGen_IIR::hwDigitalFilterGen_IIR(const hwFilterSpecs& filterSpecs,
                                               hwDigitalFilter&     digitalFilter)
    : m_pFilterSpecs (nullptr)
{
    if (filterSpecs.Order() < 1 || filterSpecs.Order() > 14)
    {
        m_status(HW_MATH_ERR_FILTERORDERIRR, 1);
        return;
    }

    m_status = filterSpecs.Status();

    if (!m_status.IsOk())   // checked after order for API reasons
    {
        return;
    }

    if (filterSpecs.BandType() == LowPass || filterSpecs.BandType() == HighPass)
    {
        digitalFilter.SetSize(filterSpecs.Order() + 1);
    }
    else // BandPass || BandStop
    {
        digitalFilter.SetSize(2 * filterSpecs.Order() + 1);
    }

    if (!digitalFilter.Status().IsOk())
    {
        m_status = digitalFilter.Status();
        return;
    }

    m_pDigitalFilter = &digitalFilter;
    m_pFilterSpecs   = &filterSpecs;
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwDigitalFilterGen_IIR::~hwDigitalFilterGen_IIR()
{
}
//------------------------------------------------------------------------------
// Scales transfer function coeffients, sets gain factor and normalizes
//------------------------------------------------------------------------------
void hwDigitalFilterGen_IIR::ScaleCoefs()
{
    hwMatrix* numer = m_pDigitalFilter->GetNumerCoefs();
    hwMatrix* denom = m_pDigitalFilter->GetDenomCoefs();

    double rippleFactor = m_pDigitalFilter->GetRippleFactor();

    double scale2 = 1.0 / (*denom)(0);
    double scale1 = rippleFactor * scale2;

    for (int i = 0; i < numer->Size(); ++i)
    {
        (*numer)(i) *= scale1;
    }

    (*denom)(0) = 1.0;
    for (int i = 1; i < denom->Size(); ++i)
    {
        (*denom)(i) *= scale2;
    }
}
