/**
* @file hwAnalogFilterGen_IIR.cxx
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
#include "hwAnalogFilterGen_IIR.h"

#include "hwMatrix.h"
#include "hwAnalogFilter.h"
#include "hwFilterSpecs.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwAnalogFilterGen_IIR::hwAnalogFilterGen_IIR(const hwFilterSpecs& filterSpecs,
                                             hwAnalogFilter&      analogFilter)
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
        analogFilter.SetSize(filterSpecs.Order() + 1);
    }
    else // BandPass || BandStop
    {
        analogFilter.SetSize(2 * filterSpecs.Order() + 1);
    }

    if (!analogFilter.Status().IsOk())
    {
        m_status = analogFilter.Status();
        return;
    }

    m_pAnalogFilter = &analogFilter;
    m_pFilterSpecs = &filterSpecs;
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwAnalogFilterGen_IIR::~hwAnalogFilterGen_IIR()
{
}
//------------------------------------------------------------------------------
// Scale transfer function coeffients. Set gain factor and normalize
//------------------------------------------------------------------------------
void hwAnalogFilterGen_IIR::ScaleCoefs()
{
    hwMatrix* numer = m_pAnalogFilter->GetNumerCoefs();
    hwMatrix* denom = m_pAnalogFilter->GetDenomCoefs();

    double rippleFactor = m_pAnalogFilter->GetRippleFactor();
    double scale2       = 1.0 / (*denom)(0);
    double scale1       = rippleFactor * scale2;

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
