/**
* @file hwWindowFunc.h
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
#ifndef _Signals_WindowFunc_h
#define _Signals_WindowFunc_h

#include "hwMatrix.h"

// forward declarations
class hwMathStatus;

//------------------------------------------------------------------------------
//!
//! \class hwWindowFunc
//! \brief Window function base class
//!
//------------------------------------------------------------------------------
class hwWindowFunc
{
public:
    //!
    //! Constructor
    //! \param periodic True if periodic window type
    //!
    hwWindowFunc(bool periodic);
    //!
    //! Constructor
    //! \param numPoints Number of points
    //! \param periodic True if periodic window type 
    //!
    hwWindowFunc(int  numPoints, 
                 bool periodic);
    //!
    //! Constructor
    //! \param weight 
    //!
    hwWindowFunc(const hwMatrix& weight);
    //!
    //! Destructor
    //!
    virtual ~hwWindowFunc();

    //!
    //! Copy constructor
    //! \param src Source
    //!
    hwWindowFunc(const hwWindowFunc& src);
    //!
    //! Copy constructor
    //! \param src Source
    //!
    hwWindowFunc& operator=(const hwWindowFunc& src);

    //!
    //! Returns status after setting window length
    //! \param numPoints Number of points
    //!
    hwMathStatus SetSize(int numPoints);
    //!
    //! Window the signal data
    //! \param data_to_window
    //! \param start
    //! \param windowed_data
    //! \param bias
    hwMathStatus ApplyWindow(const hwMatrix& data_to_window, 
                             int             start,
                             hwMatrix&       windowed_data, 
                             bool            bias = false);
    //!
    //! Returns status after computing weights
    //! \param weight
    //!
    virtual hwMathStatus ComputeWeights(hwMatrix& weight) { return hwMathStatus(); }

protected:
    //!
    //! \enum Type
    //!
    enum Type 
    {  
        Periodic, 
        Symmetric 
    };
    Type     m_windowType;  //!< Window type
    hwMatrix m_weight;      //!< Window weights

private:
    //!
    //! Helper method for Copy
    //! \param src Source
    //!
    hwMathStatus Copy(const hwWindowFunc& copy);
    //!
    //! Computes bias correction constants
    //! \param signal
    //! \param start
    //! \param k1     Windowed data mean
    //! \param k2     Power correction factor
    //!
    void Bias(const hwMatrix& signal, 
              int             start, 
              double&         k1, 
              double&         k2);
};

#endif // _Signals_WindowFunc_h
