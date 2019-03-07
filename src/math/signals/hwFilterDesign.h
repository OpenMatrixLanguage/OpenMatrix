/**
* @file  hwFilterDesign.h
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
#ifndef _Signals_hwFilterDesign_h
#define _Signals_hwFilterDesign_h

#include "hwMathStatus.h"

// forward declarations
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \class hwFilterDesign
//! \brief Filter design class
//!
//------------------------------------------------------------------------------
class hwFilterDesign
{
public:
    //!
    //! Constructor
    //!
    hwFilterDesign();
    //!
    //! Destructor
    //!
    virtual ~hwFilterDesign();

    //!
    //! Design an analog filter (low or high pass) and return the status
    //! \param passBandFreq
    //! \param stopBandFreq
    //! \param passEdgeDb
    //! \param stopEdgeDb
    //! \param order
    //! \param freqC
    //!
    hwMathStatus Analog(double    passBandFreq, 
                        double    stopBandFreq,
                        double    passEdgeDb, 
                        double    stopEdgeDb,
                        int&      order, 
                        hwMatrix& freqC);
    //!
    //! Design an analog filter (low or high pass) and return the status
    //! \param passBandFreq
    //! \param stopBandFreq
    //! \param passEdgeDb
    //! \param stopEdgeDb
    //! \param order
    //! \param freqC1
    //! \param freqC2
    //! \param cheby2
    //!
    hwMathStatus Analog(const hwMatrix& passBandFreq, 
                        const hwMatrix& stopBandFreq,
                        double          passEdgeDb, 
                        double          stopEdgeDb, 
                        int&            order,
                        hwMatrix&       freqC1, 
                        hwMatrix&       freqC2, 
                        bool            cheby2 = false);
    //!
    //! Design an digital filter (low or high pass) and return the status
    //! \param passBandFreq
    //! \param stopBandFreq
    //! \param passEdgeDb
    //! \param stopEdgeDb
    //! \param order
    //! \param freqC         Permissible range for the corner freqency
    //!
    hwMathStatus Digital(double    passBandFreq, 
                         double    stopBandFreq,
                         double    passEdgeDb, 
                         double    stopEdgeDb,
                         int&      order, 
                         hwMatrix& freqC);
    //!
    //! Design an digital filter (low or high pass) and return the status
    //! \param passBandFreq
    //! \param stopBandFreq
    //! \param passEdgeDb
    //! \param stopEdgeDb
    //! \param order
    //! \param freqC1       Permissible range for the first corner freqency
    //! \param freqC2       Permissible range for the second corner freqency
    //! \param cheby2
    //!
    hwMathStatus Digital(const hwMatrix& passBandFreq,
                         const hwMatrix& stopBandFreq,
                         double          passEdgeDb, 
                         double          stopEdgeDb, 
                         int&            order,
                         hwMatrix&       freqC1, 
                         hwMatrix&       freqC2, 
                         bool            cheby2 = false);

protected:
    hwMathStatus m_status;  //!< Status

    //!
    //! \enum Type
    //!
    enum Type 
    {
        Unknown,
        LowPass,
        HighPass,
        BandPass,
        BandStop
    };

    //!
    //! Design a prototype filter
    //! \param epsilon
    //! \param delta
    //! \param omegaSP
    //! \param order
    //! \param omegaC1
    //! \param omegaC2
    //!
    virtual void DesignPrototype(double  epsilon, 
                                 double  delta,
                                 double  omegaSP,
                                 int&    order, 
                                 double& omegaC1, 
                                 double& omegaC2) = 0;
};

#endif // _Signals_hwFilterDesign_h
