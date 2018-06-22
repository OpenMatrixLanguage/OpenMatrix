/**
* @file hwPearson_VII.h
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
#ifndef _Stats_hwPearsson_VII_h
#define _Stats_hwPearsson_VII_h

#include "hwDistribution.h"

class hwMersenneTwisterState;   // forward declaration
class hwBetaInv;                // forward declaration
class hwUniform;                // forward declaration

//------------------------------------------------------------------------------
//!
//! \class hwPearson_VII
//! \brief Pearson type VII distribution class
//!
//------------------------------------------------------------------------------
class hwPearson_VII : public hwDistribution
{
public:
    //!
    //! Constructor
    //! \param shape
    //! \param pMTState
    //!
    hwPearson_VII(double                  shape,
                  hwMersenneTwisterState* pMTState = nullptr);
    //!
    //! Destructor
    //!
    virtual ~hwPearson_VII();

    //!
    //! Returns probability density function
    //! \param x Domain value
    //!
    double Pdf(double x);
    //!
    //! Returns cumulative density function
    //! \param x Domain value
    //!
    double Cdf(double x);
    //!
    //! Returns inverse cumulative density function
    //! \param prob Probability value
    //!
    double CdfInv(double prob);
    //!
    //! Returns distribution mean
    //!
    double Mean() { return 0.0; }
    //!
    //! Returns distribution variance
    //!
    double Variance();
    //!
    //! Returns distribution mode
    //!
    double Mode() { return 0.0; }
    //!
    //! Returns distribution median
    //!
    double Median() { return 0.0; }
    //!
    //! Returns random value
    //!
    double GetDeviate();

private:
    double     m_shape;    //!< shape parameter, m in most literature
    double     m_ndof;     //!< number of degrees of freedom (student t)
    hwBetaInv* m_pBetaInv; //!< pointer to inverse Beta distribution
    hwUniform* m_pUniform; //!< pointer to uniform distribution
};

#endif // _Stats_hwPearsson_VII_h
