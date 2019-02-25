/**
* @file hwButterworth_z.cxx
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

#include "hwButterworth_z.h"

#include "hwFilterSpecs.h"
#include "hwDigitalFilterGen_AP.h"
#include "hwButterworth_Proto.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwButterworth_z::hwButterworth_z(int    order, 
                                 double lowCutoffFreq,
                                 double highCutoffFreq)
    : m_pButterProto (NULL)
{
	hwFilterSpecs         filterSpecs(order, lowCutoffFreq, highCutoffFreq, true);
    hwDigitalFilterGen_AP filterGen(filterSpecs, *this);

    m_status = filterGen.Status();
    if (!m_status.IsOk())
    {
        return;
    }

    m_pButterProto = new hwButterworth_Proto(order);
    if (!m_pButterProto->Status().IsOk())
    {
        m_status = m_pButterProto->Status();
        return;
    }

    switch (filterSpecs.BandType())
    {
        case LowPass:  filterGen.CreateLowPassFilter();  break;
        case HighPass: filterGen.CreateHighPassFilter(); break;
        case BandPass: filterGen.CreateBandPassFilter(); break;
        case BandStop: filterGen.CreateBandStopFilter(); break;
        default: break;
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwButterworth_z::~hwButterworth_z()
{
    if (m_pButterProto)
    {
        delete m_pButterProto;
    }
}
//------------------------------------------------------------------------------
// Compute location of the real pole
//------------------------------------------------------------------------------
void hwButterworth_z::GetSPlaneInfo(double omegaC, double& poleReal) const
{
    m_pButterProto->GetSPlaneInfo(poleReal);
    poleReal *= omegaC;
}
//------------------------------------------------------------------------------
// Compute real component and squared magnitude of the ith pole
//------------------------------------------------------------------------------
void hwButterworth_z::GetSPlaneInfo(int     i, 
                                    double  omegaC,
                                    double  omegaCSq, 
                                    double& poleReal,
                                    double& poleMagSq) const
{
    m_pButterProto->GetSPlaneInfo(i, poleReal, poleMagSq);
    poleReal *= omegaC;
    poleMagSq *= omegaCSq;
}
