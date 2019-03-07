/**
* @file hwPoisson.h
* @date June 2009
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
#ifndef _Stats_hwPoisson_h
#define _Stats_hwPoisson_h

// forward declarations
class hwMathStatus;
class hwMersenneTwisterState;
class hwUniform;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \class hwPoisson
//! \brief Poisson distribution class
//!
//------------------------------------------------------------------------------
class hwPoisson
{
public:
    //!
    //! Constructor
    //! \param lamda
    //! \param pMTState
    //!
    hwPoisson(double                  lambda,
              hwMersenneTwisterState* pMTState = nullptr);
    //!
    //! Destructor
    //!
    virtual ~hwPoisson();

    //!
    //! Returns probability density function
    //! \param k       Domain value
    //! \param density Probability density function value
    //!
    hwMathStatus Pdf(double  k,
                     double& density);
    //!
    //! Returns cumulative density function
    //! \param k Domain value
    //!
    double Cdf(double k);
    //!
    //! Returns cumulative density function values from 0 to k
    //! \param k   Domain value
    //! \param cdf Cumulative density function values from 0 to k
    //!
    hwMathStatus Cdf(int       k,
                     hwMatrix& cdf);
    //!
    //! Returns inverse cumulative density function
    //! \param prob Probability value
    //! \param k    Inverse cumulative density function
    //!
    hwMathStatus CdfInv(double prob,
                        int&   k);
    //!
    //! Returns inverse cumulative density function values from 0 to k
    //! \param prob Probability value
    //! \param cdf  Cumulative density function values
    //! \param k    Inverse cumulative density function values from 0 to k
    //!
    hwMathStatus CdfInv(double          prob,
                        const hwMatrix& cdf,
                        int&            k);
    //!
    //! Returns distribution mean
    //!
    double Mean() { return m_lambda; }
    //!
    //! Returns distribution variance
    //!
    double Variance() { return m_lambda; }
    //!
    //! Returns random value
    //!
    int GetDeviate();

private:
    double     m_lambda;   //!< shape parameter
    hwUniform* m_pUniform; //!< pointer to uniform distribution
};

#endif // _Stats_hwPoisson_h
