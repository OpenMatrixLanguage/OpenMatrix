/**
* @file hwAdaptiveQuadrature.cxx
* @date December 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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
#include "hwAdaptiveQuadrature.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwAdaptiveQuadrature::hwAdaptiveQuadrature(int numPnts)
    : m_kernel1(numPnts)
    , m_kernel2(numPnts+1)
    , n(numPnts)
{
}
//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
hwAdaptiveQuadrature::~hwAdaptiveQuadrature()
{
}
//------------------------------------------------------------------------------
// Returns status and gets area after integrating from a to b
//------------------------------------------------------------------------------
hwMathStatus hwAdaptiveQuadrature::Compute(const UnivarFunc pFunc, 
                                           double           a, 
                                           double           b,
                                           double&          area, 
                                           int&             count,
                                           double           reltol, 
                                           double           abstol)
{
    if (pFunc == NULL)
    {
        return hwMathStatus(HW_MATH_ERR_NULLPOINTER, 1);
    }
    if (reltol < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 6);
    }
    if (abstol < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 7);
    }
    if (a == b)
    {
        area = 0.0;
        return hwMathStatus();
    }
    if (fabs(a - b) < 1.0e-12)
    {
        return hwMathStatus(HW_MATH_ERR_QUADSTEPSIZE);
    }
    if (count > 10000)
    {
        return hwMathStatus(HW_MATH_ERR_FUNCTIONCOUNT);
    }

    hwMathStatus status;
    double       area1;
    double       area2;

    // manage improper integral
    // integration limits of either -Inf or Inf are managed
    // with a reciprocal log transform via ComputeRLog.
    // integration limits with an asymptote are managed
    // with a square root transformation via ComputeSqrt1,2.
    if (a == std::numeric_limits<double>::infinity())
    {
        if (b == std::numeric_limits<double>::infinity())
        {
            area = 0.0;
        }
        else
        {
            status = Compute(pFunc, b, a, area, count, reltol, abstol);

            if (!status.IsOk())
                return status;

            area = -area;
        }

        return status;
    }
    else if (-a == std::numeric_limits<double>::infinity())
    {
        if (b < 0.0)
        {
            status = ComputeRLog(pFunc, a, b, area, count, reltol, abstol);
        }
        else    // b >= 0.0
        {
            status = ComputeRLog(pFunc, a, -1.0, area2, count, reltol/2.0, abstol/2.0);

            if (!status.IsOk())
            {
                return status;
            }
            status = Compute(pFunc, -1.0, b, area1, count, reltol/2.0, abstol/2.0);

            if (!status.IsOk())
            {
                return status;
            }
            area = area1 + area2;
        }

        return status;
    }
    else if (b == std::numeric_limits<double>::infinity())
    {
        if (a <= 0.0)
        {
            // handle Inf by splitting interval and transforming second piece
            status = Compute(pFunc, a, 1.0, area1, count, reltol/2.0, abstol/2.0);

            if (!status.IsOk())
            {
                return status;
            }
            status = ComputeRLog(pFunc, 1.0, b, area2, count, reltol/2.0, abstol/2.0);

            if (!status.IsOk())
            {
                return status;
            }
            area = area1 + area2;
        }
        else    // a > 0.0
        {
            status = ComputeRLog(pFunc, a, b, area, count, reltol, abstol);
        }

        return status;
    }
    else if (-b == std::numeric_limits<double>::infinity())
    {
        if (-a == std::numeric_limits<double>::infinity())
        {
            area = 0.0;
        }
        else
        {
            status = Compute(pFunc, b, a, area, count, reltol, abstol);

            if (!status.IsOk())
            {
                return status;
            }
            area = -area;
        }

        return status;
    }
    else
    {
        double y1, y2;

        status = pFunc(a, y1);

        if (!status.IsOk())
        {
            return status;
        }
        status = pFunc(b, y2);

        if (!status.IsOk())
        {
            return status;
        }
        if (fabs(y1) == std::numeric_limits<double>::infinity())
        {
            if (a <= b)
            {
                status = ComputeSqrt1(pFunc, a, b, area, count, reltol, abstol);
            }
            else    // a > b
            {
                status = ComputeSqrt2(pFunc, b, a, area, count, reltol, abstol);
                area = -area;
            }

            if (!status.IsOk())
            {
                return status;
            }
            return status;
        }

        if (fabs(y2) == std::numeric_limits<double>::infinity())
        {
            if (a <= b)
            {
                status = ComputeSqrt2(pFunc, a, b, area, count, reltol, abstol);
            }
            else    // a > b
            {
                status = ComputeSqrt1(pFunc, b, a, area, count, reltol, abstol);
                area = -area;
            }

            if (!status.IsOk())
            {
                return status;
            }
            return status;
        }
    }

    // manage proper integral
    count += 2*n+1;

    status = m_kernel1.Compute(pFunc, a, b, area1);

    if (!status.IsOk())
    {
        return status;
    }
    status = m_kernel2.Compute(pFunc, a, b, area);

    if (!status.IsOk())
    {
        return status;
    }
    if (fabs(area - area1) > reltol * fabs(area) && fabs(area - area1) > abstol)
    {
        double mid = (a+b) / 2.0;
        status = Compute(pFunc, a, mid, area1, count, reltol/2.0, abstol/2.0);

        if (!status.IsOk())
        {
            return status;
        }
        status = Compute(pFunc, mid, b, area2, count, reltol/2.0, abstol/2.0);

        if (!status.IsOk())
        {
            return status;
        }
        area = area1 + area2;
    }

    return status;
}
#if 0
/*
//------------------------------------------------------------------------------
// Returns status and gets area after integrating from a to b. This method is
// is used with improper integrals with either a=-Inf or b=Inf. The integrand
// is transformed by u=1/x. The function has not been fully tested. It is to be
// retained should use cases be found that cause other methods to fail.
//------------------------------------------------------------------------------
hwMathStatus hwAdaptiveQuadrature::ComputeR(const UnivarFunc pFunc,
                                            double           a,
                                            double           b,
                                            double&          area,
                                            int&             count,
                                            double           reltol,
                                            double           abstol)
{
    double area1, area2;
    hwMathStatus status;

    if (fabs(b - a) < 1.0e-12 * fabs(a*b))
        return status(HW_MATH_ERR_QUADSTEPSIZE);

    if (count > 10000)
        return status(HW_MATH_ERR_FUNCTIONCOUNT);

    count += 2*n+1;

    status = m_kernel1.ComputeR(pFunc, a, b, area1);

    if (!status.IsOk())
        return status;

    status = m_kernel2.ComputeR(pFunc, a, b, area);

    if (!status.IsOk())
        return status;

    if (fabs(area - area1) > reltol * fabs(area) && fabs(area - area1) > abstol)
    {
        double mid = 2.0 / (1.0/a + 1.0/b);     // pre-transformed midpoint
        status = ComputeR(pFunc, a, mid, area1, count, reltol/2.0, abstol/2.0);

        if (!status.IsOk())
            return status;

        status = ComputeR(pFunc, mid, b, area2, count, reltol/2.0, abstol/2.0);

        if (!status.IsOk())
            return status;

        area = area1 + area2;
    }

    return status;
}
*/
#endif
//------------------------------------------------------------------------------
// Returns status and gets area after integrating from a to b
//------------------------------------------------------------------------------
hwMathStatus hwAdaptiveQuadrature::ComputeRLog(const UnivarFunc pFunc,
                                               double           a,
                                               double           b,
                                               double&          area,
                                               int&             count,
                                               double           reltol,
                                               double           abstol)
{
    bool rightSize;

    if (a > 0.0 && b > 0.0)
    {
        rightSize = true;

        if (fabs(1.0/log(b+1.0) - 1.0/log(a+1.0)) < 1.0e-12)
        {
            //return hwMathStatus(HW_MATH_ERR_QUADSTEPSIZE);
            area = 0.0;
            return hwMathStatus();
        }
    }
    else if (a < 0.0 && b < 0.0)
    {
        rightSize = false;

        if (fabs(1.0/log(-b+1.0) - 1.0/log(-a+1.0)) < 1.0e-12)
        {
            //return hwMathStatus(HW_MATH_ERR_QUADSTEPSIZE);
            area = 0.0;
            return hwMathStatus();
        }
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2, 3);
    }

    if (count > 10000)
    {
        return hwMathStatus(HW_MATH_ERR_FUNCTIONCOUNT);
    }
    count += 2*n+1;

    double       area1;
    hwMathStatus status = m_kernel1.ComputeRLog(pFunc, a, b, area1);
    if (!status.IsOk())
    {
        return status;
    }

    status = m_kernel2.ComputeRLog(pFunc, a, b, area);
    if (!status.IsOk())
    {
        return status;
    }

    if (fabs(area - area1) > reltol * fabs(area) && fabs(area - area1) > abstol)
    {
        double mid;
        if (rightSize)
        {
            mid = exp(2.0 / (1.0/log(a+1.0) + 1.0/log(b+1.0)))-1.0;     // pre-transformed midpoint
        }
        else
        {
            mid = 1.0-exp(2.0 / (1.0/log(-a+1.0) + 1.0/log(-b+1.0)));   // pre-transformed midpoint
        }

        status = ComputeRLog(pFunc, a, mid, area1, count, reltol/2.0, abstol/2.0);

        if (!status.IsOk())
        {
            return status;
        }

        double area2;
        status = ComputeRLog(pFunc, mid, b, area2, count, reltol/2.0, abstol/2.0);
        if (!status.IsOk())
        {
            return status;
        }
        area = area1 + area2;
    }

    return status;
}
//------------------------------------------------------------------------------
// Returns status and gets area after integrating from a to b
//------------------------------------------------------------------------------
hwMathStatus hwAdaptiveQuadrature::ComputeSqrt1(const UnivarFunc pFunc, 
                                                double           a,
                                                double           b,
                                                double&          area, 
                                                int&             count,
                                                double           reltol, 
                                                double           abstol)
{
    if (fabs(sqrt(b-a)) < 1.0e-12)
    {
        // return hwMathStatus(HW_MATH_ERR_QUADSTEPSIZE);
        area = 0.0;
        return hwMathStatus();
    }

    if (count > 10000)
    {
        return hwMathStatus(HW_MATH_ERR_FUNCTIONCOUNT);
    }
    count += 2*n+1;

    double       area1;
    hwMathStatus status = m_kernel1.ComputeSqrt1(pFunc, a, b, area1);
    if (!status.IsOk())
    {
        return status;
    }

    status = m_kernel2.ComputeSqrt1(pFunc, a, b, area);
    if (!status.IsOk())
    {
        return status;
    }

    if (fabs(area - area1) > reltol * fabs(area) && fabs(area - area1) > abstol)
    {
        double mid = (3.0*a + b) / 4.0;     // pre-transformed midpoint
        status = ComputeSqrt1(pFunc, a, mid, area1, count, reltol/2.0, abstol/2.0);

        if (!status.IsOk())
        {
            return status;
        }

        double area2;
        status = ComputeSqrt1(pFunc, mid, b, area2, count, reltol/2.0, abstol/2.0);
        if (!status.IsOk())
        {
            return status;
        }
        area = area1 + area2;
    }

    return status;
}
//------------------------------------------------------------------------------
// Returns status and gets area after integrating from a to b
//------------------------------------------------------------------------------
hwMathStatus hwAdaptiveQuadrature::ComputeSqrt2(const UnivarFunc pFunc, 
                                                double           a, 
                                                double           b,
                                                double&          area, 
                                                int&             count, 
                                                double           reltol, 
                                                double           abstol)
{
    if (fabs(sqrt(b-a)) < 1.0e-12)
    {
        // return hwMathStatus(HW_MATH_ERR_QUADSTEPSIZE);
        area = 0.0;
        return hwMathStatus();
    }

    if (count > 10000)
    {
        return hwMathStatus(HW_MATH_ERR_FUNCTIONCOUNT);
    }
    count += 2*n+1;

    double area1;
    hwMathStatus status = m_kernel1.ComputeSqrt2(pFunc, a, b, area1);
    if (!status.IsOk())
    {
        return status;
    }

    status = m_kernel2.ComputeSqrt2(pFunc, a, b, area);
    if (!status.IsOk())
    {
        return status;
    }

    if (fabs(area - area1) > reltol * fabs(area) && fabs(area - area1) > abstol)
    {
        double mid = (a + 3.0*b) / 4.0;     // pre-transformed midpoint
        status = ComputeSqrt2(pFunc, a, mid, area1, count, reltol/2.0, abstol/2.0);
        if (!status.IsOk())
        {
            return status;
        }

        double area2;
        status = ComputeSqrt2(pFunc, mid, b, area2, count, reltol/2.0, abstol/2.0);
        if (!status.IsOk())
        {
            return status;
        }

        area = area1 + area2;
    }

    return status;
}
