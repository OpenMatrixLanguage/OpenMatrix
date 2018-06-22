/**
* @file  hwISO6487.h
* @date June 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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
#ifndef _CAE_ISO6487_h
#define _CAE_ISO6487_h

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \class hwISO6487
//! \brief ISO6487 filter class. See Appendix A.
//!
//------------------------------------------------------------------------------
class hwISO6487
{
public:
    //! 
    //! Constructor
    //!
    hwISO6487();
    //! 
    //! Destructor
    //!
    ~hwISO6487();

    //!
    //! Filters the signal and returns the status
    //! \param response
    //! \param sampFreq Sampling frequency
    //! \param cfc      Channel frequency class
    //! \param output
    //!
    hwMathStatus Filter(const hwMatrix& response, 
                        double          sampFreq, 
                        double          cfc,
                        hwMatrix&       output);

private:
    //!
    //! Returns the status after padding the signal
    //! \param response Given input
    //! \param padcnt   Pad count
    //! \param pad_resp Padded output
    //!
    hwMathStatus FillIt(const hwMatrix& response, 
                        int             padcnt, 
                        hwMatrix&       pad_resp);
    //!
    //! Returns the status after filtering the signal in place
    //! \param response Given signal
    //! \param delta  
    //! \param cfc      Channel frequency class
    //! \param dir      Direction
    //!
    void Filter(hwMatrix& response, double delta, double cfc, int dir);
};

#endif // _CAE_ISO6487_h
