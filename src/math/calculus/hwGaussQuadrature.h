/**
* @file hwGaussQuadrature.h
* @date June 2007
* Copyright (C) 2007-2019 Altair Engineering, Inc.  
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
#ifndef _Calculus_GaussQuadrature_h
#define _Calculus_GaussQuadrature_h

#include "hwMatrix.h"

//!
//! \typedef hwMathStatus
//!
typedef hwMathStatus (*QuadFunc1)(const hwMatrix& x, hwMatrix& y);

//------------------------------------------------------------------------------
//! 
//! \class hwGaussQuadrature
//! \brief Gaussian quadrature class 
//!
//------------------------------------------------------------------------------
class hwGaussQuadrature
{
public:
    //!
    //! Constructor
    //!
    hwGaussQuadrature();
    //!
    //! Destructor
    //!
    virtual ~hwGaussQuadrature();

    //!
    //! Returns the status
    //!
    const hwMathStatus& GetStatus() const { return m_status; }
    //!
    //! Sets points
    //! \param n_ Number of points
    //!
    virtual void SetPnts(int n_) = 0;
    //!
    //! Returns status and gets area after integrating pFunc(z) from a to b
    //! \param pFunc Input function
    //! \param a     Input a
    //! \param b     Input b
    //! \param area  Output
    //!
    hwMathStatus Compute(const QuadFunc1 pFunc,
                         double          a,
                         double          b, 
                         double&         area);
    //!
    //! Returns status and gets area after integrating from a to b. This method
    //! is used with improper integrals with either a=-Inf or b=Inf. a and b 
    //! MUST have the same sign. The integrand is transformed by 
    //! u=1/log(abs(x)+1)
    //! \param pFunc Input function
    //! \param a     Input a
    //! \param b     Input b
    //! \param area  Output
    //!
    hwMathStatus ComputeRLog(const QuadFunc1 pFunc,
                             double          a,
                             double          b, 
                             double&         area);
    //!
    //! Returns status and gets area after integrating from a to b. This method
    //! is used with improper integrals with an integrand that has a vertical
    //! asymptote at a. The integrand is transformed by u=sqrt(x-a).  It is 
    //! required that a <= b
    //! \param pFunc Input function
    //! \param a     Input a
    //! \param b     Input b
    //! \param area  Output
    //!
    hwMathStatus ComputeSqrt1(const QuadFunc1 pFunc,
                              double          a,
                              double          b, 
                              double&         area);
    //!
    //! Returns status and gets area after integrating from a to b. This method
    //! is used with improper integrals with an integrand that has a vertical
    //! asymptote at b. The integrand is transformed by u=sqrt(b-x). It is
    //! required that a <= b
    //! \param pFunc Input function
    //! \param a     Input a
    //! \param b     Input b
    //! \param area  Output
    //!
    hwMathStatus ComputeSqrt2(const QuadFunc1 pFunc,
                              double          a,
                              double          b, 
                              double&         area);
protected:
    hwMathStatus m_status;           //!< status 
    int n;                           //!< number of points
    hwMatrix X;                      //!< locations
    hwMatrix W;                      //!< weights

    virtual void GetLocations() = 0; //!< Populates locations
    virtual void GetWeights()   = 0; //!< Populates weights
};
//------------------------------------------------------------------------------
//! 
//! \class hwGaussLegendre
//! \brief Provides exact results for polynomials of order up to 2n-1
//!
//------------------------------------------------------------------------------
class hwGaussLegendre : public hwGaussQuadrature
{
public:
    //!
    //! Constructor
    //! \param n_ Number of points
    //!
    hwGaussLegendre(int n_);
    //!
    //! Destructor
    //!
    virtual ~hwGaussLegendre();

    //!
    //! Sets points
    //! \param n Number of points
    //!
    void SetPnts(int n);

protected:
    void GetLocations(); //!< Populates locations
    void GetWeights();   //!< Populates weights
};
//------------------------------------------------------------------------------
//! 
//! \class hwGaussLobatto
//! \brief Provides exact results for polynomials of order up to 2n-3
//!
//------------------------------------------------------------------------------
class hwGaussLobatto : public hwGaussQuadrature
{
public:
    //!
    //! Constructor
    //! \param n_ Number of points
    //!
    hwGaussLobatto(int n_);
    //!
    //! Destructor
    //!
    virtual ~hwGaussLobatto();

    //!
    //! Sets points
    //! \param n Number of points
    //!
    void SetPnts(int n_);

protected:
    void GetLocations(); //!< Populates locations
    void GetWeights();   //!< Populates weights
};

#endif // _Calculus_GaussQuadrature_h
