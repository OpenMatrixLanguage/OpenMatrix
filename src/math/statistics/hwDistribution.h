/**
* @file hwDistribution.h
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
#ifndef _Stats_hwDistribution_h
#define _Stats_hwDistribution_h

#include "StatisticsExports.h"

//------------------------------------------------------------------------------
//!
//! \class hwDistribution
//! \brief Continuous distribution base class
//!
//------------------------------------------------------------------------------
class STATISTICS_DECLS hwDistribution
{
public:
    //!
    //! Constructor
    //!
    hwDistribution();
    //!
    //! Destructor
    //!
    virtual ~hwDistribution();

    //!
    //! Returns probability density function
    //! \param x Domain value
    //!
    virtual double Pdf(double x) = 0;
    //!
    //! Returns cumulative density function
    //! \param x Domain value
    //!
    virtual double Cdf(double x) = 0;
    //!
    //! Returns inverse cumulative density function
    //! \param prob Probability value
    //!
    virtual double CdfInv(double prob) = 0;
    //!
    //! Returns distribution mean
    //!
    virtual double Mean() = 0;
    //!
    //! Returns distribution variance
    //!
    virtual double Variance() = 0;
    //!
    //! Returns distribution median
    //!
    virtual double Median();
    //!
    //! Returns distribution mode
    //!
    virtual double Mode() = 0;
    //!
    //! Returns random value
    //!
    virtual double GetDeviate() = 0;
};

#endif // _Stats_hwDistribution_h
