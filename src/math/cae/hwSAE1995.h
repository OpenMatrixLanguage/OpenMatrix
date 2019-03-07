/**
* @file hwSAE1995.h
* @date October, 2009
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
#ifndef _CAE_SAE1995_h
#define _CAE_SAE1995_h

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \class hwSAE1995
//! \brief SAE J211-1 (1995) filter class. See Appendix C.
//!
//------------------------------------------------------------------------------
class hwSAE1995
{
public:
    //!
    //! Constructor
    //!
    hwSAE1995();
    //!
    //! Destructor
    //!
    ~hwSAE1995();

    //!
    //! Returns status after applying filter to the signal data
    //! \param response
    //! \param sampFreq  Sampling frequency
    //! \param cfc
    //! \param stdpad    0(no padding), <0(pad with zeros), >0(pad using mirrors)
    //! \param direction
    //! \param output
    hwMathStatus Filter(const hwMatrix& response, 
                        double          sampFreq,
                        double          cfc, 
                        int             stdpad,
                        int             direction, 
                        hwMatrix&       output);

private:
    //!
    //! Returns status after padding source vector
    //! \param source     Source vector
    //! \param padmode    Pad type
    //! \param padcnt  
    //! \param source_pad Padded source vector
    hwMathStatus PadVec(const hwMatrix& source, 
                        int             padmode,
                        int             padcnt, 
                        hwMatrix&       source_pad);
    //!
    //! Performs the filtering
    //! \param response
    //! \param sampFreq  Sampling frequency
    //! \param cfc
    //! \param direction Direction of output
    //! \param output    Output
    void Filter1D(const hwMatrix& response, 
                  double          sampFreq, 
                  double          cfc,
                  int             direction, 
                  hwMatrix&       output);
};

#endif // _CAE_SAE1995_h
