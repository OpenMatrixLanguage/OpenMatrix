/**
* @file StatisticsTests.h
* @date June 2012
* Copyright (C) 2012-2018 Altair Engineering, Inc.  
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
#ifndef _Stats_StatisticsTests_h
#define _Stats_StatisticsTests_h

#include "StatisticsExports.h"

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \brief Statisitical tests
//!
//------------------------------------------------------------------------------

//!
//! z test
//! \param data
//! \param mu
//! \param sigma
//! \param reject
//! \param p
//! \param CI
//! \param alpha
//!
STATISTICS_DECLS hwMathStatus ZTest(const hwMatrix& data,
                                    double          mu,
                                    double          sigma,
                                    bool&           reject,
                                    double&         p,
                                    hwMatrix&       CI,
                                    double          alpha = 0.05);
//!
//! One sample t Test
//! \param data
//! \param mu
//! \param reject
//! \param p
//! \param CI
//! \param alpha
//!
STATISTICS_DECLS hwMathStatus TTest(const hwMatrix& data,
                                    double          mu,
                                    bool&           reject,
                                    double&         p,
                                    hwMatrix&       CI,
                                    double          alpha = 0.05);
//!
//! Two sample t Test
//! \param data1
//! \param data2
//! \param reject
//! \param p
//! \param CI
//! \param alpha
//!
STATISTICS_DECLS hwMathStatus TTest2(const hwMatrix& data1,
                                     const hwMatrix& data2,
                                     bool&           reject,
                                     double&         p,
                                     hwMatrix&       CI,
                                     double          alpha = 0.05);
//!
//! Chi-squared Test
//! \param data
//! \param var
//! \param reject
//! \param p
//! \param CI
//! \param alpha
//!
STATISTICS_DECLS hwMathStatus ChiSqTest(const hwMatrix& data,
                                        double          var,
                                        bool&           reject,
                                        double&         p,
                                        hwMatrix&       CI,
                                        double          alpha = 0.05);
//!
//! F Test
//! \param data1
//! \param data2
//! \param reject
//! \param p
//! \param CI
//! \param alpha
//!
STATISTICS_DECLS hwMathStatus FTest(const hwMatrix& data1,
                                    const hwMatrix& data2,
                                    bool&           reject,
                                    double&         p,
                                    hwMatrix&       CI,
                                    double          alpha = 0.05);

#endif // _Stats_StatisticsTests_h
