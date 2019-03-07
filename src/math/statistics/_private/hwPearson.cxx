/**
* @file hwPearson.cxx
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
#include "hwPearson.h"

#include "hwMatrix.h"
#include "hwNormal.h"
#include "hwBeta.h"
#include "hwPearson_III.h"
#include "hwPearson_IV.h"
#include "hwPearson_V.h"
#include "hwBetaInv.h"
#include "hwPearson_VII.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwPearson::hwPearson(const hwMatrix& moment, hwMersenneTwisterState* pMTState)
    : m_pMTState         (nullptr)
    , m_pStdDistribution (nullptr)
    , m_type             (-1)
    , m_bTailSwitch      (false)
{
    if (!moment.IsReal())
    {
        return;
    }

    if (!moment.IsVector())
    {
        return;
    }

    if (moment.Size() != 4)
    {
        return;
    }

    m_mean     = moment(0);
    m_variance = moment(1);
    m_skewness = moment(2);
    m_kurtosis = moment(3);
    m_pMTState = pMTState;

    FindType();
    if (m_type != -1)
    {
        ConstructDistribution();
    }
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwPearson::hwPearson(double                  mean,
                     double                  var,
                     double                  skew,
                     double                  kurt,
                     hwMersenneTwisterState* pMTState)
    : m_pMTState         (pMTState)
    , m_mean             (mean)
    , m_variance         (var)
    , m_skewness         (skew)
    , m_kurtosis         (kurt)
    , m_pStdDistribution (nullptr)
    , m_type             (-1)
    , m_bTailSwitch      (false)
{
    FindType();
    if (m_type != -1)
    {
        ConstructDistribution();
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwPearson::~hwPearson()
{
    delete m_pStdDistribution;
}
//------------------------------------------------------------------------------
// Compute distribution parameters from the moments
//------------------------------------------------------------------------------
hwMathStatus hwPearson::GetParams(hwMatrix& param)
{
    hwMathStatus status = param.Dimension(4, hwMatrix::REAL);
    if (!status.IsOk())
    {
        status.SetArg1(1);
        return status;
    }

    param(0) = m_shape1;
    param(1) = m_shape2;
    param(2) = m_scale;
    param(3) = m_offset;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Find the Pearson family type
//------------------------------------------------------------------------------
void hwPearson::FindType()
{
    double nearZero = 1.0e-8;
    double beta1    = m_skewness * m_skewness;
    double beta2    = m_kurtosis + 3.0;
    double denom    = 4.0 * (4.0 * beta2 - 3.0 * beta1) * (2.0 * beta2 - 3.0 * beta1 - 6.0);

    if (beta2 - beta1 - 1.0 < 0.0)
    {
        m_type = -1;
        return;
    }

    if (fabs(denom) < nearZero)
    {
        m_type = (beta1 < nearZero) ? 0 :     // Normal
                                      3 ;     // modified Gamma
    }
    else
    {
        double numer = beta1 * (beta2 + 3.0) * (beta2 + 3.0);
        double kappa = numer / denom;

        if (fabs(kappa) < nearZero)
        {
            if (fabs(beta2 - 3.0) < nearZero)
            {
                m_type = 0;						// Normal
            }
            else if (beta2 < 3.0)
            {
                m_type = 2;						// Symmetric Beta
            }
            else  //(beta2 > 3.0)
            {
                m_type = 7;						// (generalized Student T)
            }
        }
        else if (kappa < 0.0)
        {
            m_type = 1;							// Beta
        }
        else if (fabs(kappa - 1.0) < nearZero)
        {
            m_type = 5;							// modified Inverse Gamma
        }
        else if (kappa < 1.0)
        {
            m_type = 4;							// Johnson Su
        }
        else // kappa > 1.0
        {
            m_type = 6;							// Inverse Beta (generalized F)
        }
    }
}
//------------------------------------------------------------------------------
// Construct distribution for the Pearson family type
//------------------------------------------------------------------------------
void hwPearson::ConstructDistribution()
{
    // compute parameters and construct standard distribution
    double beta1 = m_skewness * m_skewness;
    double beta2 = m_kurtosis + 3.0;
    double r;
    double value1;
    double value2;

    m_bTailSwitch = false;

    switch(m_type)
    {
    case 0:		// Normal
        m_pStdDistribution = new hwNormal(m_pMTState);
        break;

    case 1:		// Beta
        r = 6.0*(beta2 - beta1 - 1.0) / (6.0 + 3.0*beta1 - 2.0*beta2);	// a + b
        value1 = sqrt(beta1 * (r + 2.0) * (r + 2.0) + 16.0 * (r + 1.0));
        value2 = m_skewness * (r + 2.0) / value1;
        m_shape1 = 0.5 * r * (1.0 - value2);				// a
        m_shape2 = 0.5 * r * (1.0 + value2);				// b
        m_pStdDistribution = new hwBeta(m_shape1, m_shape2, m_pMTState);
        break;

    case 2:		// Symmetric Beta
        r = 6.0 * (beta2 - 1.0) / (6.0 - 2.0 * beta2);	// a + b
        value1 = 4.0 * sqrt(r + 1.0);
        m_shape1 = 0.5 * r;								// a
        m_pStdDistribution = new hwBeta(m_shape1, m_shape1, m_pMTState);
        break;

    case 3:		// modified Gamma
        m_shape1 = 4.0 / beta1 - 1.0;						// a
        m_pStdDistribution = new hwPearson_III(m_shape1, m_pMTState);
        if (m_skewness < 0.0)
        {
            m_bTailSwitch = true;
        }
        break;

    case 4:		// Pearson type IV
        r = 6.0 * (beta2-beta1-1.0) / (2.0*beta2-3.0*beta1-6.0); // 2(m-1)
        m_shape1 = 0.5 * r + 1.0;    // m
        m_shape2 = -(r*(r-2)*sqrt(beta1)) / sqrt(16.0*(r-1.0)-beta1*(r-2.0)*(r-2.0));   // nu

        m_pStdDistribution = new hwPearson_IV(m_shape1, m_shape2, m_pMTState);
        if (m_skewness < 0.0)
        {
            m_bTailSwitch = true;
        }
        break;

    case 5:		// modified Inverse Gamma
        m_shape1 = 4.0 + (8.0 + 4.0 * sqrt(4.0 + beta1)) / beta1;
        m_pStdDistribution = new hwPearson_V(m_shape1, m_pMTState);
        if (m_skewness < 0.0)
        {
            m_bTailSwitch = true;
        }
        break;

    case 6:		// Inverse Beta (generalized F)
        r = 6.0 * (beta2 - beta1 - 1.0) / (6.0 + 3.0 * beta1 - 2.0 * beta2);
        value1 = sqrt(beta1 * (r + 2.0) * (r + 2.0) + 16.0 * (r + 1.0));
        value2 = fabs(m_skewness) * (r + 2.0) / value1;
        m_shape1 = 0.5 * r * (1.0 + value2);				// a
        m_shape2 = 1.0 - r;								// b
        m_pStdDistribution = new hwBetaInv(m_shape1, m_shape2, m_pMTState);
        if (m_skewness < 0.0)
        {
            m_bTailSwitch = true;
        }
        break;

    case 7:		// (generalized Student T)
        m_shape1 = 2.5 + 3.0 / (beta2 - 3.0);
        m_pStdDistribution = new hwPearson_VII(m_shape1, m_pMTState);
        break;

    default: m_pStdDistribution = nullptr; break;
    }

    m_scale = sqrt(m_variance / m_pStdDistribution->Variance());

    if (m_bTailSwitch)
    {
        m_scale = -m_scale;
    }

    m_offset = m_mean - m_pStdDistribution->Mean() * m_scale;
}
//------------------------------------------------------------------------------
// Probability density function
//------------------------------------------------------------------------------
double hwPearson::Pdf(double x)
{
    double xx = x;

    xx -= m_offset;
    xx /= m_scale;

    return m_pStdDistribution->Pdf(xx) / fabs(m_scale);
}
//------------------------------------------------------------------------------
// Cummulative density function
//------------------------------------------------------------------------------
double hwPearson::Cdf(double x)
{
    double xx = x;

    xx -= m_offset;
    xx /= m_scale;

    if (m_bTailSwitch)
    {
        return 1.0 - m_pStdDistribution->Cdf(xx);
    }

    return m_pStdDistribution->Cdf(xx);
}
//------------------------------------------------------------------------------
// Inverse cummulative density function
//------------------------------------------------------------------------------
double hwPearson::CdfInv(double prob)
{
    double x = (m_bTailSwitch) ? m_pStdDistribution->CdfInv(1.0 - prob) :
                                 m_pStdDistribution->CdfInv(prob);

    x *= m_scale;
    x += m_offset;

    return x;
}
//------------------------------------------------------------------------------
// Compute the median of the distribution
//------------------------------------------------------------------------------
double hwPearson::Median()
{
    double median = m_pStdDistribution->Median();

    median *= m_scale;
    median += m_offset;

    return median;
}
//------------------------------------------------------------------------------
// Compute the mode of the distribution
//------------------------------------------------------------------------------
double hwPearson::Mode()
{
    double mode = m_pStdDistribution->Mode();

    mode *= m_scale;
    mode += m_offset;

    return mode;
}
//------------------------------------------------------------------------------
// Compute a random deviate from the distribution
//------------------------------------------------------------------------------
double hwPearson::GetDeviate()
{
    double deviate = m_pStdDistribution->GetDeviate();

    deviate *= m_scale;
    deviate += m_offset;

    return deviate;
}
