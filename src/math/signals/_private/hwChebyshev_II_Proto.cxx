/**
* @file  hwChebyshev_II_Proto.h
* @date April 2009
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
#include "hwChebyshev_II_Proto.h"

#include <math.h>

#include "GeneralFuncs.h"


//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwChebyshev_II_Proto::hwChebyshev_II_Proto(int    order,
                                           double stopEdgeDb)
	: hwLowPass_Proto(order)
{
    if (!m_status.IsOk())
    {
        return;
    }

    if (stopEdgeDb <= 0.0)
    {
        m_status(HW_MATH_ERR_DB_SIGN, 2);
        return;
    }

    m_epsilon = 1.0 / sqrt(exp(stopEdgeDb * log(10.0) / 10.0) - 1.0);

    if (IsZero(m_epsilon, 1.0e-12))
    {
        m_status(HW_MATH_ERR_FILTERSPEC_E);
        return;
    }

    double value = asinh(1.0 / m_epsilon) / (double) m_order;

    m_sinh = sinh(value);
    m_cosh = cosh(value);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwChebyshev_II_Proto::~hwChebyshev_II_Proto()
{
}
//------------------------------------------------------------------------------
// Compute location of the real pole
//------------------------------------------------------------------------------
void hwChebyshev_II_Proto::GetSPlaneInfo(double& poleReal) const
{
    // real pole
    poleReal = -1.0 / m_sinh;
}
//------------------------------------------------------------------------------
// Compute real component/squared magnitude of the ith pole, and the squared 
// magnitude of the ith zero
//------------------------------------------------------------------------------
void hwChebyshev_II_Proto::GetSPlaneInfo(int     i, 
                                         double& poleReal,
                                         double& poleMagSq, 
                                         double& zeroMagSq) const
{
    // complex conjugate pole pair
    double angle = PI * (2*i + 1) / (double) m_order;
    double value = sin(0.5 * angle);

    poleReal = -m_sinh * value / (m_cosh * m_cosh - value * value);
    poleMagSq = 1.0 / (m_cosh * m_cosh - value * value);
    zeroMagSq = 2.0 / (1.0 + cos(angle));   // use half angle identity for cos^2()
}
