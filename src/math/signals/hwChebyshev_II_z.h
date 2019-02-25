/**
* @file  hwChebyshev_II_z.h
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
#ifndef _Signals_hwChebyshev_II_z_h
#define _Signals_hwChebyshev_II_z_h

#include "hwDigitalFilter.h"

class hwChebyshev_II_Proto; // forward declaration

//------------------------------------------------------------------------------
//!
//! \class hwChebyshev_II_z
//! \brief Chebyshev II digital filter class
//!
//------------------------------------------------------------------------------
class hwChebyshev_II_z : public hwDigitalFilter
{
public:
    //!
    //! Constructor
    //! \param order
    //! \param lowCutoffFreq
    //! \param highCutoffFreq
    //! \param stopEdgeDb
    //!
    hwChebyshev_II_z(int    order, 
                     double lowCutoffFreq,
                     double highCutoffFreq, 
                     double stopEdgeDb);
    //!
    //! Destructor
    //!
    virtual ~hwChebyshev_II_z();

    //!
    //! Compute location of the real pole
    //! \param omegaC
    //! \param poleReal real pole
    //!
    void GetSPlaneInfo(double  omegaC, 
                       double& poleReal) const;
    //!
    //! Compute real component/squared magnitude of the ith pole, and the
    //! squared magnitude of the ith zero
    //! \param i
    //! \param omegaC
    //! \param omegaCSq
    //! \param poleReal  Real component of ith pole
    //! \param poleMagSq Squared magnitude of the ith pole
    //! \param zeroMagSq Squared magnitude of the ith zero
    //!
    void GetSPlaneInfo(int i, 
                       double  omegaC, 
                       double  omegaCSq,
                       double& poleReal, 
                       double& poleMagSq, 
                       double& zeroMagSq) const;

protected:
    hwChebyshev_II_Proto* m_pCheby_II_Proto; //!< pointer to prototype filter
};

#endif // _Signals_hwChebyshev_II_z_h
