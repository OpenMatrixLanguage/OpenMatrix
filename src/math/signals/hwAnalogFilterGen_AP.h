/**
* @file hwAnalogFilterGen_AP.h
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
#ifndef _Signals_hwAnalogFilter_AP_h
#define _Signals_hwAnalogFilter_AP_h

#include "hwAnalogFilterGen_IIR.h"

//------------------------------------------------------------------------------
//!
//! \brief Analog filters with only poles class. All Poles in the S plane, no zeros
//!
//------------------------------------------------------------------------------
class hwAnalogFilterGen_AP : public hwAnalogFilterGen_IIR
{
public:
    //!
    //! Constructor
    //! \param filterSpecs
    //! \param analogFilter
    //!
    hwAnalogFilterGen_AP(const hwFilterSpecs& filterSpecs, 
                         hwAnalogFilter&      analogFilter);
    //! 
    //! Destructor
    //!
    virtual ~hwAnalogFilterGen_AP();

    //!
    //! Compute transfer function polynomial coefficients for a low pass filter
    //!
    void CreateLowPassFilter();
    //!
    //! Compute transfer function polynomial coefficients for a high pass filter
    //!
    void CreateHighPassFilter();
    //!
    //! Compute transfer function polynomial coefficients for a band pass filter
    //!
    void CreateBandPassFilter();
    //!
    //! Compute transfer function polynomial coefficients for a band stop filter
    //!
    void CreateBandStopFilter();
};

#endif // _Signals_hwAnalogFilter_AP_h
