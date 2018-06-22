/**
* @file  hwFilterManager.h
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
#ifndef _Signals_hwFilterManager_h
#define _Signals_hwFilterManager_h

#include "hwDigitalFilter.h"

//------------------------------------------------------------------------------
//!
//! \class hwFilterManager
//! \brief Digital filter management class
//!
//------------------------------------------------------------------------------
class hwFilterManager
{
public:
    //!
    //! Constructor
    //!
    hwFilterManager();
    //!
    //! Destructor
    //!
    virtual ~hwFilterManager();

    //!
    //! Returns status after creating filter with transfer function polynomials
    //! \param numerCoef
    //! \param denomCoef
    //!
    hwMathStatus CreateFilter(const hwMatrix& numerCoef, 
                              const hwMatrix* denomCoef = NULL);
    //!
    //! Returns status after applying filter to signal data
    //! \param inSignal  Input signal
    //! \param outSignal Output signal
    //! \param initCond  Initial conditions
    //!
    hwMathStatus ApplyFilter(const hwMatrix& inSignal, 
                             hwMatrix&       outSignal,
                             const hwMatrix* initCond = NULL);
    //!
    //! Returns status after applying filter to signal data in reverse
    //! \param inSignal  Input signal
    //! \param outSignal Output signal
    //! \param initCond  Initial conditions
    //!
    hwMathStatus ApplyReverseFilter(const hwMatrix& inSignal, 
                                    hwMatrix&       outSignal,
                                    const hwMatrix* initCond = NULL);
    //!
    //! Returns status after applying zero phase filter to signal data
    //! \param inSignal  Input signal
    //! \param outSignal Output signal
    //!
    hwMathStatus ApplyZeroPhaseFilter(const hwMatrix& inSignal,
                                      hwMatrix&       outSignal);
    //!
    //! Returns status and ups sample by a given factor
    //! \param inSignal  Input signal
    //! \param outSignal Output signal
    //! \param k         Up sample factor
    //! \param phase     Offset number of samples
    //!
    hwMathStatus UpSample(const hwMatrix& inSignal, 
                          hwMatrix&       outSignal,
                          int             k,
                          int             phase = 0);
    //!
    //! Returns status and downs sample by a given factor
    //! \param inSignal  Input signal
    //! \param outSignal Output signal
    //! \param k         Up sample factor
    //! \param phase     Offset number of samples
    //!
    hwMathStatus DownSample(const hwMatrix& inSignal, 
                            hwMatrix&       outSignal,
                            int             k,
                            int             phase = 0);

protected:
    hwDigitalFilter* m_pFilter; //!< Filter
};

#endif // _Signals_hwFilterManager_h
