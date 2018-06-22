/**
* @file  hwFilterSpecs.h
* @date April 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language (�OpenMatrix�) software.
* Open Source License Information:
* OpenMatrix is free software. You can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* OpenMatrix is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
* You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
* 
* Commercial License Information: 
* For a copy of the commercial license terms and conditions, contact the Altair Legal Department at Legal@altair.com and in the subject line, use the following wording: Request for Commercial License Terms for OpenMatrix.
* Altair�s dual-license business model allows companies, individuals, and organizations to create proprietary derivative works of OpenMatrix and distribute them - whether embedded or bundled with other software - under a commercial license agreement.
* Use of Altair�s trademarks and logos is subject to Altair's trademark licensing policies.  To request a copy, email Legal@altair.com and in the subject line, enter: Request copy of trademark and logo usage policy.
*/
#ifndef _hwFilterSpecs_h
#define _hwFilterSpecs_h

#include "hwMathStatus.h"

//!
//! \enum hwFilterType
//!
enum hwFilterType 
{
    Unknown, 
    LowPass,
    HighPass, 
    BandPass, 
    BandStop 
};

//------------------------------------------------------------------------------
//!
//! \class hwFilterSpecs
//! \brief Filter base class specifying the type
//!
//------------------------------------------------------------------------------
class hwFilterSpecs
{
public:
    //! 
    //! Constructor
    //! \param order 
    //! \param lowCutoffFreq  Low cut off frequency
    //! \param highCutoffFreq High cut off frequency
    //! \param digital        True if digital
    //!
    hwFilterSpecs(int    order, 
                  double lowCutoffFreq,
                  double highCutoffFreq, 
                  bool   digital);
    //!
    //! Destructor
    //!
    virtual ~hwFilterSpecs();

    //!
    //! Gets the status
    //!
    const hwMathStatus& Status() const { return m_status; }
    //!
    //! Gets the order
    //!
    int Order() const { return m_order; }
    //! 
    //! Gets the lower corner frequency
    //!
    double GetLowerCorner() const { return m_lowerCorner; }
    //! 
    //! Gets the upper corner frequency
    //!
    double GetUpperCorner() const { return m_upperCorner; }
    //!
    //! Gets the band type
    hwFilterType BandType() const { return m_filterType; }

protected:
    int          m_order;         //!< Order
    double       m_lowerCorner;   //!< corner frequency - not used for LowPass
    double       m_upperCorner;   //!< corner frequency - not used for HighPass
    hwFilterType m_filterType;    //!< Band type
    hwMathStatus m_status;        //!< Status
};

#endif // _hwFilterSpecs_h
