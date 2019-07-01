/**
* @file hwAdaptiveQuadrature.h
* @date December 2015
* Copyright (C) 2015-2019 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language ("OpenMatrix") software.
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
#ifndef _Calculus_AdaptiveQuadrature_h
#define _Calculus_AdaptiveQuadrature_h

#include "hwGaussQuadrature.h"

//!
//! \typedef hwMathStatus
//!
typedef hwMathStatus (*QuadFunc1)(const hwMatrix& x, hwMatrix& y);

//------------------------------------------------------------------------------
//! 
//! \class hwAdaptiveQuadrature
//! \brief Gaussian quadrature class 
//!
//------------------------------------------------------------------------------
class hwAdaptiveQuadrature
{
public:
    //!
    //! Constructor
    //! \param numPnts Number of points
    //!
    hwAdaptiveQuadrature(int numPnts);
    //!
    //! Destructor
    //!
    virtual ~hwAdaptiveQuadrature();
    //!
    //! Returns status after integrating from a to b and gets area. The method
    //! manages proper and improper integrals. Transformations are used for
    //! improper integrals
    //! \param pFunc  Input function  
    //! \param a      Input a
    //! \param b      Input b
    //! \param area   Output 
    //! \param count
    //! \param reltol Optional relative tolerance
    //! \param abs    Optional absolute tolerance
    //!
    hwMathStatus Compute(const QuadFunc1 pFunc, 
                         double          a, 
                         double          b,
                         double&         area, 
                         int&            count,
                         double          reltol = 1.0e-6, 
                         double          abstol = 1.0e-6);
private:
    int n;                      // the order of the first approximation
    hwGaussLegendre m_kernel1;  // the order n integration object 
    hwGaussLegendre m_kernel2;  // the order n+1 integration object

    //!
    //! Returns status and gets area after integrating from a to b. Used with
    //! improper integrals with either a=-Inf or b=Inf. a and b MUST have the 
    //! same sign. The integrand is transformed by u=1/log(abs(x)+1).
    //! \param pFunc  Input function  
    //! \param a      Input a
    //! \param b      Input b
    //! \param area   Output 
    //! \param count
    //! \param reltol Relative tolerance
    //! \param abs    Absolute tolerance
    //!
    hwMathStatus ComputeRLog(const QuadFunc1 pFunc, 
                             double          a, 
                             double          b,
                             double&         area, 
                             int&            count, 
                             double          reltol, 
                             double          abstol);
    //!
    //! Returns status and gets area after integrating from a to b. This method
    //! is used with improper integrals with an integrand that has a vertical
    //! asymptote at a. The integrand is transformed by u=sqrt(x-a).
    //! \param pFunc  Input function  
    //! \param a      Input a
    //! \param b      Input b
    //! \param area   Output 
    //! \param count
    //! \param reltol Relative tolerance
    //! \param abs    Absolute tolerance
    //!
    hwMathStatus ComputeSqrt1(const QuadFunc1 pFunc, 
                              double          a, 
                              double          b,
                              double&         area, 
                              int&            count, 
                              double          reltol, 
                              double          abstol);
    //!
    //! Returns status and gets area after integrating from a to b. This method 
    //! is used with improper integrals with an integrand that has a vertical
    //! asymptote at b. The integrand is transformed by u=sqrt(b-x).
    //! \param pFunc  Input function  
    //! \param a      Input a
    //! \param b      Input b
    //! \param area   Output 
    //! \param count
    //! \param reltol Relative tolerance
    //! \param abs    Absolute tolerance
    //!
    hwMathStatus ComputeSqrt2(const QuadFunc1 pFunc, 
                              double          a, 
                              double          b,
                              double&         area, 
                              int&            count, 
                              double          reltol, 
                              double          abstol);
};

#endif // _Calculus_AdaptiveQuadrature_h
