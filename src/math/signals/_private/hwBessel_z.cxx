/**
* @file hwBessel_z.cxx
* @date June 2012
* Copyright (C) 2012-2018 Altair Engineering, Inc.  
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
#include "hwBessel_z.h"

#include "hwFilterSpecs.h"
#include "hwDigitalFilterGen_AP.h"
#include "hwBessel_Proto.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwBessel_z::hwBessel_z(int         order, 
                       double      lowCutoffFreq,
                       double      highCutoffFreq, 
                       const char* type)
{
	hwFilterSpecs filterSpecs(order, lowCutoffFreq, highCutoffFreq, true);
    m_pBesselProto = NULL;

    hwDigitalFilterGen_AP filterGen(filterSpecs, *this);

    m_status = filterGen.Status();
    if (!m_status.IsOk())
    {
        return;
    }

    m_pBesselProto = new hwBessel_Proto(order, type);
    m_status = m_pBesselProto->Status();

    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 2)
        {
            m_status.SetArg1(4);
        }
        if (!m_status.IsWarning())
        {
            return;
        }
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
hwBessel_z::~hwBessel_z()
{
    if (m_pBesselProto)
    {
        delete m_pBesselProto;
    }
}
//------------------------------------------------------------------------------
// Compute location of the real pole
//------------------------------------------------------------------------------
void hwBessel_z::GetSPlaneInfo(double omegaC, double& poleReal) const
{
    m_pBesselProto->GetSPlaneInfo(poleReal);
    poleReal *= omegaC;
}
//------------------------------------------------------------------------------
// Compute real component and squared magnitude of the ith pole
//------------------------------------------------------------------------------
void hwBessel_z::GetSPlaneInfo(int     i, 
                               double  omegaC,
                               double  omegaCSq,
                               double& poleReal, 
                               double& poleMagSq) const
{
    m_pBesselProto->GetSPlaneInfo(i, poleReal, poleMagSq);
    poleReal *= omegaC;
    poleMagSq *= omegaCSq;
}
