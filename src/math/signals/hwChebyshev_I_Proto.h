/**
* @file hwChebyshev_I_Proto.h
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
#ifndef _Signals_hwChebyshev_I_Proto_h
#define _Signals_hwChebyshev_I_Proto_h

#include "hwLowPass_Proto.h"

//------------------------------------------------------------------------------
//!
//! \class hwChebyshev_I_Proto
//! \brief Chebyshev I prototype low pass filter class. Generates s plane data 
//! for a low pass Chebyshev I prototype filter that is transformed into digital
//! and analog filters having the various passband types
//!
//------------------------------------------------------------------------------
class hwChebyshev_I_Proto : public hwLowPass_Proto
{
public:
    //!
    //! Constructor
    //! \param order
    //! \param passEdgeDb
    //!
    hwChebyshev_I_Proto(int order, double passEdgeDb);
    //!
    //! Destructor
    //!
    virtual ~hwChebyshev_I_Proto();

    //!
    //! Compute location of the real pole
    //! \param poleReal Real pole
    //!
    void GetSPlaneInfo(double& poleReal) const;
    //!
    //! Compute real component and squared magnitude of the ith pole
    //! \param i
    //! \param poleReal  Real component of ith pole
    //! \param poleMagSq Squared magnitude of the ith pole
    //!
    void GetSPlaneInfo(int     i, 
                       double& poleReal, 
                       double& poleMagSq) const;
    //!
    //! Compute ripple factor at DC
    //!
    double GetRippleFactor() const;

protected:
    double m_epsilon;  //!< constant in ripple calculation
    double m_sinh;     //!< constant in pole/zero calculation
    double m_cosh;     //!< constant in pole/zero calculation
};

#endif // _Signals_hwChebyshev_I_Proto_h
