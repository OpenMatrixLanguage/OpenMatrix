/**
* @file hwPearson_IV.h
* @date October 2014
* Copyright (C) 2014-2018 Altair Engineering, Inc.  
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
#ifndef _Stats_hwPearsson_IV_h
#define _Stats_hwPearsson_IV_h

#include "hwDistribution.h"

// forward declarations
class hwMersenneTwisterState;
class hwUniform;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \class hwPearson_IV
//! \brief Pearson type IV distribution class
//!
//------------------------------------------------------------------------------
class hwPearson_IV : public hwDistribution
{
public:
    //!
    //! Constructor
    //! \param m
    //! \param nu
    //! \param pMTState
    //!
    hwPearson_IV(double                  m, 
                 double                  nu,
                 hwMersenneTwisterState* pMTState = nullptr);
    //!
    //! Destructor
    //!
    virtual ~hwPearson_IV();

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
    double Mean() { return m_mean; }
    //!
    //! Returns distribution variance
    //!
    double Variance() { return m_variance; }
    //!
    //! Returns distribution mode
    //!
    double Mode();
    //!
    //! Returns random value
    //!
    double GetDeviate();

private:
    double     m_m;        //!< shape parameter
    double     m_nu;       //!< shape parameter
    double     m_mean;     //!< mean
    double     m_variance; //!< variance
    double     m_k;        //!< pdf normalization constant
    hwUniform* m_pUniform; //!< pointer to uniform distribution

    //!
    //! Compute abs(gamma(x+iy)/gamma(x))^2
    //! \param x
    //! \param y
    //!
    double gammar2(double x, 
                   double y);
    //! 
    //! Compute pdf normalization constant
    //!
    double type4norm();
    //!
    //! Gauss-Legender quadrature
    //! \param a
    //! \param b
    //! \param X
    //! \param W
    //!
    double GaussQuad(double          a, 
                     double          b, 
                     const hwMatrix& X, 
                     const hwMatrix& W);
};

#endif // _Stats_hwPearsson_IV_h
