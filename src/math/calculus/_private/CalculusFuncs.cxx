/**
* @file CalculusFuncs.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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
#include "CalculusFuncs.h"
#include "hwAdaptiveQuadrature.h"

//------------------------------------------------------------------------------
// Computes forward differences and return status
//------------------------------------------------------------------------------
hwMathStatus Derivative(const hwMatrix& x,
                        const hwMatrix& y,
                        hwMatrix&       derivative)
{
    if (!x.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!x.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int size = x.Size();
    if (size < 2)
    {
        return hwMathStatus(HW_MATH_ERR_TOOFEWPOINTS, 1);
    }

    if (!y.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!y.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    if (y.Size() != size)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    hwMathStatus status = derivative.Dimension(x.M(), x.N(), hwMatrix::REAL);
    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(3);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    double delta_x;
    if (size > 1)
    {
        delta_x = x(1) - x(0);
        derivative(0) = (y(1) - y(0)) / delta_x;
    }

    for (int k = 1; k < size-1; k++)
    {
        delta_x = x(k+1) - x(k-1);
        derivative(k) = (y(k+1) - y(k-1)) / delta_x;
    }

    if (size > 2)
    {
        delta_x = x(size-1) - x(size-2);
        derivative(size-1) = (y(size-1) - y(size-2)) / delta_x;
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes integral using trapezoidal method and returns status
//------------------------------------------------------------------------------
hwMathStatus TrapZ(const hwMatrix& x,
                   const hwMatrix& y,
                   double&         integral)
{    
    if (!x.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!x.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (!y.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!y.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }

    int size = x.Size();
    if (y.Size() != size)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    integral = 0.0;
    for (int k = 1; k < size; k++)
        integral += 0.5 * (x(k) - x(k-1)) * (y(k) + y(k-1));

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes cumulative integral using trapezoidal method and returns status
//------------------------------------------------------------------------------
hwMathStatus CumTrapZ(const hwMatrix& x,
                      const hwMatrix& y,
                      hwMatrix&       integral)
{
    if (!x.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!x.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (!y.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!y.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }

    int size = x.Size();
    if (y.Size() != size)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    hwMathStatus status = integral.Dimension(x.M(), x.N(), hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(3);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    if (size)
    {
        integral(0) = 0.0;
    }

    for (int k = 1; k < size; k++)
    {
        integral(k) = integral(k-1) + 0.5 * (x(k) - x(k-1)) * (y(k) + y(k-1));
    }
    return status;
}
//------------------------------------------------------------------------------
// Computes integral using Gauss-Legendre quadrature and returns status
//------------------------------------------------------------------------------
hwMathStatus GaussLegendre(const UnivarFunc pFunc, 
                           double           a,
                           double           b, 
                           int              numPnts, 
                           double&          area)
{
    hwGaussLegendre quadObj(numPnts);
    hwMathStatus status = quadObj.GetStatus();

    if (!status.IsOk())
    {
        if (status.GetArg1() == 1)
        {
            status.SetArg1(4);
        }
        return status;
    }

    status = quadObj.Compute(pFunc, a, b, area);
    return status;
}
//------------------------------------------------------------------------------
// Computes integral using Gauss-Lobatto quadrature and returns status
//------------------------------------------------------------------------------
hwMathStatus GaussLobatto(const UnivarFunc pFunc, 
                          double           a, 
                          double           b,
                          int              numPnts, 
                          double&          area)
{
    hwGaussLobatto quadObj(numPnts);
    hwMathStatus status = quadObj.GetStatus();

    if (!status.IsOk())
    {
        if (status.GetArg1() == 1)
        {
            status.SetArg1(4);
        }
        return status;
    }

    status = quadObj.Compute(pFunc, a, b, area);
    return status;
}
//------------------------------------------------------------------------------
// Computes integral using adaptive quadrature and returns status
//------------------------------------------------------------------------------
hwMathStatus Quad(const UnivarFunc pFunc, 
                  double           a, 
                  double           b, 
                  double&          area,
                  int&             count, 
                  double           reltol, 
                  double           abstol)
{
    hwAdaptiveQuadrature quadObj(4);
    count = 0;

    hwMathStatus stat = quadObj.Compute(pFunc, a, b, area, count, reltol, abstol);
    return stat;
}
//------------------------------------------------------------------------------
// Helper for adaptive Simpson's rule
//------------------------------------------------------------------------------
static hwMathStatus AdaptiveSimpsonsRule(const UnivarFunc pFunc, 
                                         double           a, 
                                         double           b, 
                                         double           c,
                                         double           fa, 
                                         double           fb, 
                                         double           fc, 
                                         double&          area,
                                         int&             count, 
                                         double           tol)
{    
    double h = b - a;
    if (h == 0.0)
    {
        area = 0.0;
        return hwMathStatus();
    }
    if (h < 1.0e-12)
    {
        return hwMathStatus(HW_MATH_ERR_QUADSTEPSIZE);
    }

    double d = (a + c) / 2.0;
    double fd;
    hwMathStatus status = pFunc(d, fd);
    if (!status.IsOk())
    {
        return status;
    }

    double e = (c + b) / 2.0;
    double fe;
    status = pFunc(e, fe);
    if (!status.IsOk())
    {
        return status;
    }

    count       += 2;
    double area1 = (h/12.0) * (fa + 4.0*fd + fc);
    double area2 = (h/12.0) * (fc + 4.0*fe + fb);
    double sum   = area1 + area2;

    if (fabs(sum - area) <= 15.0 * tol)
    {
        area = sum + (sum - area) / 15.0;

        if (count > 10000)
        {
            return status(HW_MATH_ERR_FUNCTIONCOUNT);
        }
        return status;
    }

    status = AdaptiveSimpsonsRule(pFunc, a, c, d, fa, fc, fd, area1, count, tol/2.0);

    if (!status.IsOk())
    {
        return status;
    }

    status = AdaptiveSimpsonsRule(pFunc, c, b, e, fc, fb, fe, area2, count, tol/2.0);
    area   = area1 + area2;

    return status;
}
//------------------------------------------------------------------------------
// Computes integral using adaptive Simpson's rule and returns status
//------------------------------------------------------------------------------
hwMathStatus QuadV(const UnivarFunc pFunc, 
                   double           a, 
                   double           b, 
                   double&          area,
                   int&             count, 
                   double           tol)
{
    hwMathStatus status;

    if (b < a)
    {
        status = QuadV(pFunc, b, a, area, count, tol);
        area   = -area;
        return status;
    }

    double fa;
    status = pFunc(a, fa);

    if (fabs(fa) == std::numeric_limits<double>::infinity())
    {
        status = pFunc(a + MACHEP*(b-a), fa);
    }
    if (!status.IsOk())
    {
        return status;
    }

    double fb;
    status = pFunc(b, fb);

    if (fabs(fb) == std::numeric_limits<double>::infinity())
    {
        status = pFunc(b - MACHEP*(b-a), fb);
    }
    if (!status.IsOk())
    {
        return status;
    }

    double c = (a + b) / 2.0;
    double fc;
    status = pFunc(c, fc);
    if (!status.IsOk())
    {
        return status;
    }

    double h = b - a;
    count    = 3;
    area     = (h/6.0) * (fa + 4.0 * fc + fb);

    status = AdaptiveSimpsonsRule(pFunc, a, b, c, fa, fb, fc, area, count, tol);
    return status;
}
