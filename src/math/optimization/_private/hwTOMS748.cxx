/**
* @file hwTOMS748.cxx
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
#include "hwTOMS748.h"

#include "GeneralFuncs.h"

//------------------------------------------------------------------------------
// Determines the termination criterion
// tol = 2*(2*EPS*|B| + 10D-{NEPS}),  IF NEPS IS NOT 1000
// tol = 2*(2*EPS*|B|),               IF NEPS = 1000
//------------------------------------------------------------------------------
void tole(double  b,
          int     neps,
          double  eps,
          double& tol)
{
    if (neps == 1000)
    {
        tol = 0.0;
    }
    else
    {
        tol = 1.0;
        for (int i = 0; i < neps; ++i)
        {
            tol /= 10.0;
        }
    }

    tol += 2.0 * fabs(b) * eps;
    tol *= 2.0;
}
//------------------------------------------------------------------------------
// Gets a shrinked enclosing interval and updates termination criterion
//------------------------------------------------------------------------------
void Bracket(const TOMS748Func f, 
             double&           a, 
             double&           b, 
             double            c, 
             double&           fa, 
             double&           fb,
             double&           tol, 
             int               neps,
             double            eps, 
             double&           d,
             double&           fd)
{
    // Given the current enclosing interval [a, b] and a number c in [a, b], if
    // f(c) = 0, sets output a = c. Otherwise the new enclosing interval is
    // determined: [a. b] = [a, c] or [a, b] = [c, b]. The termination 
    // criterion corresponds to the new enclosing interval

    // Adjust c if (b - a) is very small or if c is very close to a or b
    tol *= 0.7;

    if (b - a <= 2.0 * tol)
    {
        c = a + 0.5 * (b - a);
    }
    else if (c <= a + tol)
    {
        c = a + tol;
    }
    else if (c >= b - tol)
    {
        c = b - tol;
    }

    double fc;
    hwMathStatus status = f(c, fc);

    if (fc == 0.0)
    {
        a  = c;
        fa = 0.0;
        d  = 0.0;
        fd = 0.0;
    }

    // If f(c) is not zero, get the new enclosing interval
    if (fa * fc < 0.0)
    {
        d  = b;
        fd = fb;
        b  = c;
        fb = fc;
    }
    else
    {
        d  = a;
        fd = fa;
        a  = c;
        fa = fc;
    }

    // Update termination criterion to new enclosing interval
    if (fabs(fb) <= fabs(fa))
    {
        tole(b, neps, eps, tol);
    }
    else
    {
        tole(a, neps, eps, tol);
    }
}
//------------------------------------------------------------------------------
// Indicates sign of the input, 1 if positive, -1 if negative and 0 otherwise
//------------------------------------------------------------------------------
int ISIGN(double x)
{
    if (x > 0.0)
    {
        return 1;
    }
    else if (x == 0.0)
    {
        return 0;
    }

    return -1;
}
//------------------------------------------------------------------------------
// Uses Newton steps to approximate 0 in (A,B) of the quadratic polynomial
// interpolating f(x) at a, b and d. Safeguard is used to avoid overflow
//------------------------------------------------------------------------------
void NewtonQuad(double  a, 
                double  b, 
                double  d, 
                double  fa, 
                double  fb, 
                double  fd,
                double& c, 
                int     k)
{
    // Initialization
    // Find the coefficeints of the quadratic polynomial
    int ierror = 0;

    double a0 = fa;
    double a1 = (fb - fa)  / (b - a);
    double a2 = ((fd - fb) / (d - b) - a1) / (d - a);

    double pc;
    double pdc;

    // Safeguard to avoid overflow
    do
    {
        if (a2 == 0.0 || ierror == 1)
        {
            c = a - a0 / a1;
            return;
        }

        // Determine the starting point of the Newton steps
        c = (ISIGN(a2) * ISIGN(fa) > 0) ? a : b;

        // Start the safeguarded Newton steps
        for (int i = 0; i < k; ++i)
        {
            if (ierror == 0)
            {
                pc  = a0 + (a1 + a2 * (c - b)) * (c - a);
                pdc = a1 + a2 * (2.0 * c - (a + b));

                if (pdc == 0.0)
                {
                    ierror = 1;
                }
                else
                {
                    c -= pc / pdc;
                }
            }
        }
    }
    while (ierror == 1);
}
//------------------------------------------------------------------------------
// Uses cubic inverse interpolation of f(x) at a, b, d and e to get approximate
// root of f(x). This procedure is a slight modification of Aitken-Neville
// algorithm for interpolation described by Stoer and Burlirsch in the
// Introduction to Numerical Analysis, Springer-Verlag, New York (1980)
// d and e lie outside the interval [a,b]
//------------------------------------------------------------------------------
void Pzero(double  a, 
           double  b, 
           double  d, 
           double  e, 
           double  fa, 
           double  fb, 
           double  fd,
           double  fe, 
           double& c)
{        
    double q11 = (d - e) * fd / (fe - fd);
    double q21 = (b - d) * fb / (fd - fb);
    double q31 = (a - b) * fa / (fb - fa);
    double d21 = (b - d) * fd / (fd - fb);
    double d31 = (a - b) * fb / (fb - fa);
    double q22 = (d21 - q11) * fb / (fe - fb);
    double q32 = (d31 - q21) * fa / (fd - fa);
    double d32 = (d31 - q21) * fd / (fd - fa);
    double q33 = (d32 - q22) * fa / (fe - fa);

    // calculate the output
    c = q31 + q32 + q33;
    c += a;
}
//------------------------------------------------------------------------------
// Finds zero of a function without derivates using TOMS 748 algorithm
//------------------------------------------------------------------------------
hwMathStatus TOMS748(const TOMS748Func f, 
                     double&           a, 
                     double&           b, 
                     double&           fa, 
                     double&           fb,
                     double&           root, 
                     double&           froot, 
                     double            eps, 
                     int&              numIters, 
                     int&              numFunEvals)
{
    // Subroutine rroot(nprob, neps, eps, a, b, root)
    // Finds either an exact solution or an approximate solution of the equation
    // f(x) = 0 in the interval [a, b]. At the beginning of each iteration. the
    // current enclosing interval is recorded as [a0, b0]. The first step is 
    // simply a secant step. Starting with the second iteration, three steps are
    // taken in each iteration. First two steps are either quadratic 
    // interpolation or cubic inverse interpolation. The third step is a double
    // size secant step. If the diameter of the enclosing interval obtained 
    // after those three steps is larger than 0.5 * (b0, a0), an additional
    // bisection step will be taken

    // Initialization. Number of iterations is set to 0. Calls function to get
    // initial function values f(a) and f(b). Dummy values used for e and fe

    int maxIters    = numIters;
    int maxFunEvals = numFunEvals;
    int neps        = 1000;

    double c;
    double d;
    double u;
    double fd;
    double fu;
    double a0;
    double b0;
    double tol;
    double prof;

    hwMathStatus status = f(a, fa);
    status = f(b, fb);

    numFunEvals = 2;
    numIters    = 0;

    double mu = 0.5;
    double e  = 1.0e5;
    double fe = 1.0e5;
    
    while (numIters < maxIters)
    {
        // Iteration starts. The enclosing interval before executing the 
        // iteration is recorded as [a0, b0]
        a0 = a;
        b0 = b;

        numIters++; // Update iteration number

        // Calculates termination criterion and stops if criterion is satisfied
        if (fabs(fb) < fabs(fa))
        {
            tole(b, neps, eps, tol);
        }
        else
        {
            tole(a, neps, eps, tol);
        }

        if (b - a <= tol)
        {
            break;
        }

        if (numIters == 1) // Take secant step for first iteration
        {
            c = a - fa / (fb - fa) * (b - a);

            // Calls Bracket to get a shrinked enclosing enclosing interval and
            // update the termination criterion. Stops if criterion is met or 
            // exact solution is obtained
            Bracket(f, a, b, c, fa, fb, tol, neps, eps, d, fd);
            ++numFunEvals;

            if (fa == 0.0 || b - a <= tol || numFunEvals == maxFunEvals)
                break;

            continue;
        }

        // Starting with the second iteration, in the first two steps, either
        // quadratic interpolation is used by calling NewtonQuad or the cubic
        // inverse interpolation is used by calling the method Pzero. If prof
        // is not equal to 0, then the four function values fa, fb, fd and fe
        // are distinct and hence pzero will be called
        prof = (fa-fb) * (fa-fd) * (fa-fe) * (fb-fd) * (fb-fe) * (fd-fe);

        if (numIters == 2 || prof == 0.0)
        {
            NewtonQuad(a, b, d, fa, fb, fd, c, 2);
        }
        else
        {
            Pzero(a, b, d, e, fa, fb, fd, fe, c);

            if ((c-a) * (c-b) >= 0.0)
                NewtonQuad(a, b, d, fa, fb, fd, c, 2);
        }

        e  = d;
        fe = fd;

        // Call Bracket to get a shrinked enclosing interval as well as to 
        // update the termination criterion. Execution is stopped if the 
        // criterion is met or the exact solution is met
        Bracket(f, a, b, c, fa, fb, tol, neps, eps, d, fd);
        ++numFunEvals;
        if (fa == 0.0 || (b-a) <= tol || numFunEvals == maxFunEvals)
        {
            break;
        }

        prof = (fa-fb) * (fa-fd) * (fa-fe) * (fb-fd) * (fb-fe) * (fd-fe);
        if (prof == 0.0)
        {
            NewtonQuad(a, b, d, fa, fb, fd, c, 3);
        }
        else
        {
            Pzero(a, b, d, e, fa, fb, fd, fe, c);

            if ((c-a) * (c-b) > 0.0)
            {
                NewtonQuad(a, b, d, fa, fb, fd, c, 3);
            }
        }

        // Bracket gets a shrinked enclosing interval and updates the termination
        // criterion. If the criterion is met or solution is obtained, execution
        // is quit
        Bracket(f, a, b, c, fa, fb, tol, neps, eps, d, fd);
        ++numFunEvals;
        if (fa == 0.0 || (b-a) <= tol || numFunEvals == maxFunEvals)
        {
            break;
        }
        e  = d;
        fe = fd;

        // Double size secant step
        if (fabs(fa) < fabs(fb))
        {
            u  = a;
            fu = fa;
        }
        else
        {
            u  = b;
            fu = fb;
        }

        c = u - 2.0 * fu / (fb-fa) * (b-a);

        if (fabs(c-u) > 0.5 * (b-a))
        {
            c = a + 0.5 * (b-a);
        }

        // Bracket gets a shrinked enclosing interval and updates the termination
        // criterion.
        Bracket(f, a, b, c, fa, fb, tol, neps, eps, d, fd);
        ++numFunEvals;
        if (fa == 0.0 || (b-a) <= tol || numFunEvals == maxFunEvals)
        {
            break;
        }

        // Performs additional bisection step is needed
        if (b-a < mu * (b0 - a0))
        {
            continue;
        }
        e  = d;
        fe = fd;

        // Bracket gets a shrinked enclosing interval and updates the termination
        // criterion.
        Bracket(f, a, b, a + 0.5*(b-a), fa, fb, tol, neps, eps, d, fd);
        ++numFunEvals;
        if (fa == 0.0 || b-a <= tol || numFunEvals == maxFunEvals)
        {
            break;
        }
    }

    root  = a;
    froot = fa;

    if (numIters == maxIters)
    {
        status(HW_MATH_WARN_MAXITERATE);
    }
    else if (numFunEvals == maxFunEvals)
    {
        status(HW_MATH_WARN_MAXFUNCEVAL);
    }
    return status;
}
