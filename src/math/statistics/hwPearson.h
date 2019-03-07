/**
* @file hwPearson.h
* @date May 2009
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
#ifndef _Stats_hwPearson_h
#define _Stats_hwPearson_h

#include "hwDistribution.h"

// forward declarations
class hwMathStatus;
class hwMersenneTwisterState;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//------------------------------------------------------------------------------
//!
//! \class hwPearson_IV
//! \brief Pearson family distribution class
//!
//------------------------------------------------------------------------------
class hwPearson : public hwDistribution
{
public:
    //!
    //! Constructor
    //! \param moment
    //! \param pMTState
    //!
    hwPearson(const hwMatrix&         moment,
              hwMersenneTwisterState* pMTState = nullptr);
    //!
    //! Constructor
    //! \param mean
    //! \param var
    //! \param skew
    //! \param kurt
    //! \param pMTState
    //!
    hwPearson(double                  mean,
              double                  var,
              double                  skew,
              double                  kurt,
              hwMersenneTwisterState* pMTState = nullptr);
    //!
    //! Destructor
    //!
    virtual ~hwPearson();

    //!
    //! Returns distribution type
    //!
    int Type() { return m_type; }
    //!
    //! Returns distribution parameters
    //! \param param Distribution parameters
    //!
    hwMathStatus GetParams(hwMatrix& param);
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
    //! Returns distribution median
    //!
    double Median();
    //!
    //! Returns distribution mode
    //!
    double Mode();
    //!
    //! Returns random value
    //!
    double GetDeviate();

private:
    int                     m_type;             //!< Pearson family type
    hwDistribution*         m_pStdDistribution; //!< pointer to standard distribution
    double                  m_mean;             //!< distribution mean
    double                  m_variance;         //!< distribution variance
    double                  m_skewness;         //!< distribution skewness
    double                  m_kurtosis;         //!< distribution kurtosis
    double                  m_shape1;           //!< distribution shape parameter
    double                  m_shape2;           //!< distribution shape parameter
    double                  m_scale;            //!< distribution scale parameter
    double                  m_offset;           //!< distribution location parameter
    bool                    m_bTailSwitch;      //!< switch for base distributions that only skew one way
    hwMersenneTwisterState* m_pMTState;         //!< pointer to MersenneTwister

    //!
    //! Determines distribution type
    //!
    void FindType();
    //!
    //! Constructs standard distribution
    //!
    void ConstructDistribution();
};

#endif // _Stats_hwPearson_h
