/**
* @file  hwElliptic_z.cxx
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

#include "hwElliptic_z.h"

#include "hwFilterSpecs.h"
#include "hwDigitalFilterGen_ZP.h"
#include "hwElliptic_Proto.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwElliptic_z::hwElliptic_z(int    order, 
                           double lowCutoffFreq,
                           double highCutoffFreq, 
                           double passEdgeDb,
                           double stopEdgeDb)
    : m_pElliptic_Proto (nullptr)
{
	hwFilterSpecs filterSpecs(order, lowCutoffFreq, highCutoffFreq, true);

    hwDigitalFilterGen_ZP filterGen(filterSpecs, *this);

    m_status = filterGen.Status();
    if (!m_status.IsOk())
    {
        return;
    }

    m_pElliptic_Proto = new hwElliptic_Proto(order, passEdgeDb, stopEdgeDb);
    m_status          = m_pElliptic_Proto->Status();

    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 2)
        {
            m_status.SetArg1(4);
        }
        else if (m_status.GetArg1() == 3)
        {
            m_status.SetArg1(5);
        }
        if (m_status.GetArg2() == 3)
        {
            m_status.SetArg2(5);
        }
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
hwElliptic_z::~hwElliptic_z()
{
}
//------------------------------------------------------------------------------
// Compute location of the real pole
//------------------------------------------------------------------------------
void hwElliptic_z::GetSPlaneInfo(double omegaC, double& poleReal) const
{
    m_pElliptic_Proto->GetSPlaneInfo(poleReal);
    poleReal *= omegaC;
}
//------------------------------------------------------------------------------
// Compute real component, squared magnitude of the ith pole and the squared
// magnitude of the ith zero
//------------------------------------------------------------------------------
void hwElliptic_z::GetSPlaneInfo(int     i, 
                                 double  omegaC,
                                 double  omegaCSq, 
                                 double& poleReal,
                                 double& poleMagSq, 
                                 double& zeroMagSq) const
{
    m_pElliptic_Proto->GetSPlaneInfo(i, poleReal, poleMagSq, zeroMagSq);
    poleReal *= omegaC;
    poleMagSq *= omegaCSq;
    zeroMagSq *= omegaCSq;
}
//------------------------------------------------------------------------------
// Compute low pass ripple factor at DC
//------------------------------------------------------------------------------
double hwElliptic_z::GetRippleFactor() const
{
    return m_pElliptic_Proto->GetRippleFactor();
}
