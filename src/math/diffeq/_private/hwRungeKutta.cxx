/**
* @file hwRungeKutta.cxx
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
#include "hwRungeKutta.h"

#include <math.h>
#include <stdlib.h>

#include <GeneralFuncs.h>

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwRungeKutta::hwRungeKutta(RKF45Fn_client  sysfunc, 
                           const hwMatrix& y,
                           double          relerr_, 
                           const hwMatrix* abserr_,
                           const hwMatrix* userData)
	: hwDiffEqSolver(y)
    , work (nullptr)
{
    if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 1)
        {
            m_status.SetArg1(2);
        }
        else
        {
            m_status.ResetArgs();
        }
        return;
    }

    if (!sysfunc)
    {
        m_status(HW_MATH_ERR_NULLPOINTER, 1);
        return;
    }

    m_pFunc = sysfunc;

    if (relerr_ <= 0.0)
    {
        m_status(HW_MATH_ERR_NONPOSITIVE, 3);
        return;
    }

    relerr = relerr_;

    int numEqns = m_y.Size();

    hwMathStatus status = abserr.Dimension(numEqns, hwMatrix::REAL);

    if (!status.IsOk())
    {
        m_status(HW_MATH_ERR_ALLOCFAILED);
        return;
    }

    if (abserr_)
    {
        if (abserr_->IsEmpty())
        {
            for (int i = 0; i < numEqns; i++)
            {
                abserr(i) = 1.0e-6;
            }
        }
        else
        {
            if (!abserr_->IsReal())
            {
                m_status(HW_MATH_ERR_COMPLEX, 4);
                return;
            }

            if (!abserr_->IsVector())
            {
                m_status(HW_MATH_ERR_VECTOR, 4);
                return;
            }

            if (abserr_->Size() == 1)
            {
                if ((*abserr_)(0) < 0.0)
                {
                    m_status(HW_MATH_ERR_NEGATIVE, 4);
                    return;
                }

                // abserr = (*abserr_)(0);
                for (int i = 0; i < numEqns; i++)
                {
                    abserr(i) = (*abserr_)(0);
                }
            }
            else if (abserr_->Size() == numEqns)
            {
                abserr = (*abserr_);

                for (int i = 0; i < numEqns; i++)
                {
                    if (abserr(i) < 0.0)
                    {
                        m_status(HW_MATH_ERR_NEGATIVE, 4);
                        return;
                    }
                }
            }
            else
            {
                m_status(HW_MATH_ERR_ARRAYSIZE, 2, 4);
                return;
            }
        }
    }
    else
    {
        for (int i = 0; i < numEqns; i++)
        {
            abserr(i) = 1.0e-6;
        }
    }

    work = new double[6*numEqns];   // single allocation - carried over from FORTRAN version
    yp = work;
    f1 = yp + numEqns;
    f2 = f1 + numEqns;
    f3 = f2 + numEqns;
    f4 = f3 + numEqns;
    f5 = f4 + numEqns;

    eps = 2.22044604925031308085e-16;   // = 2.0 * MACHEP
    u26 = 26.0*eps;

    m_userData = userData;
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwRungeKutta::~hwRungeKutta()
{
    if (work)
    {
        delete [] work;
    }
}
//------------------------------------------------------------------------------
// RungeKF45      (Runge-Kutta-Fehlberg 45 ODE solver)
// Ported from FORTRAN to C++ by Brian Shier 10/2006
// Arguments changed to double precision by jmb 10/19/95
//
// fehlberg fourth-fifth order runge-kutta method
// written by h.a.watts and l.f.shampine, Sandia Lboratories, Albuquerque, NM
//
// Rkf45 is primarily designed to solve non-stiff and mildly stiff differential 
// equations when derivative evaluations are inexpensive. Rkf45 should generally 
// not be used when the user is demanding high accuracy.
// Abstract:
// Rkf45 integrates a system of neqn first order ordinary differential equations 
// of the form
// dy(i)/dt = f(t,y(1),y(2),...,y(neqn))
// where the y(i) are given at t. Typically the subroutine is used to integrate 
// from t to tout but it can be used as a one-step integrator to advance the 
// solution a single step in the direction of tout. On return the parameters in 
// the call list are set for continuing the integration. The user has only to 
// call rkf45 again (and perhaps define a new value for tout).
// Rkf45 calls Fehl which computes an approximate solution over one step.
//
// Rkf45 uses the runge-kutta-fehlberg (4,5) method described in the reference 
// E. Fehlberg , low-order classical runge-kutta formulas with stepsize control, 
// nasa tr r-315 
//
// The performance of rkf45 is illustrated in the reference
// L. F. Shampine, H. A. Watts, S. Davenport, solving non-stiff ordinary 
// differential equations - the state of the art, Sandia Lboratories report 
// sand75-0182, to appear in SIAM review.
//
// The parameters represent
// pFunc -- function to evaluate derivatives yp(i)=dy(i)/dt
// neqn  -- number of equations to be integrated (now m_nEqn)
// y(*)  -- solution vector at t
// t     -- independent variable
// tout  -- output point at which solution is desired
// iflag -- indicator for status of integration 
// relerr, abserr -- relative and absolute error tolerances for local error test
// at each step the code requires that abs(local error) .le. relerr*abs(y) + 
// abserr for each component of the local error and solution vectors
//        
// First call to rkf45
// The user must provide storage in his calling program for y(neqn), declare 
// pFunc in an external statement, supply subroutine pFunc(t,y,yp) and 
// initialize the following parameters.
//
// neqn -- number of equations to be integrated.  (neqn .ge. 1)
// y(*) -- vector of initial conditions
// t    -- starting point of integration, must be a variable
// tout -- output point at which solution is desired. t=tout is allowed on the 
// first call only, in which case Rkf45 returns with iflag=2 if continuation is 
// possible.
// relerr, abserr -- relative and absolute local error tolerances which must be 
// non-negative. relerr must be a variable while abserr may be a constant. The 
// code should normally not be used with relative error control smaller than 
// about 1.e-8. To avoid limiting precision difficulties the code requires 
// relerr to be larger than an internally computed relative error parameter 
// which is machine dependent. In particular, pure absolute error is not 
// permitted. If a smaller than allowable value of relerr is attempted, rkf45 
// increases relerr appropriately and returns control to the user before 
// continuing the integration.
// iflag -- +1,-1  indicator to initialize the code for each new problem. Normal 
// input is +1. The user should set iflag=-1 only when one-step integrator 
// control is essential. In this case, Rkf45 attempts to advance the solution a 
// single step in the direction of tout each time it is called. Since this mode
// of operation results in extra computing overhead, it should be avoided unless 
// needed.
//
// Output from rkf45
// y(*)      -- solution at t
// t         -- last point reached in integration
// iflag = 2 -- integration reached tout. indicates successful return and is the 
//              Normal mode for continuing integration.
//       =-2 -- a single successful step in the direction of tout has been taken.
//              Normal mode for continuing integration one step at a time. 
//       = 3 -- integration was not completed because relative error tolerance 
//              was too small. relerr has been increased appropriately for 
//              continuing.
//       = 4 -- integration was not completed because more than 3000 derivative
//              evaluations were needed. This is approximately 500 steps.
//       = 5 -- Integration was not completed because solution vanished making a 
//              pure relative error test impossible. Must use non-zero abserr to 
//              continue. Using the one-step integration mode for one step is a 
//              good way to proceed.
//       = 6 -- integration was not completed because requested accuracy could
//              not be achieved using the smallest allowable stepsize. User must 
//              increase the error tolerance before continued integration can be 
//              attempted.
//       = 7 -- it is likely that rkf45 is inefficient for solving this problem. 
//              Too much output is restricting the natural stepsize choice. Use 
//              the one-step integrator mode.
//       = 8 -- invalid input parameters this indicator occurs if any of the 
//              following is satisfied -
//              neqn .le. 0
//              t=tout  and  iflag .ne. +1 or -1
//              relerr or abserr .lt. 0.
//              iflag .eq. 0  or  .lt. -2  or  .gt. 8
//       = 9 -- the user function failed because of invalid output vector length
//
// Subsequent calls to rkf45
// Rkf45 returns with all information needed to continue the integration. If the 
// integration reached tout, the user need only define a new tout and call rkf45 
// again. in the one-step integrator mode (iflag=-2) the user must keep in mind 
// that each step taken is in the direction of the current tout. upon reaching 
// tout (indicated by changing iflag to 2),the user must then define a new tout 
// and reset iflag to -2 to continue in the one-step integrator mode.
//
// If the integration was not completed but the user still wants to continue 
// (iflag=3,4 cases), he just calls rkf45 again. With iflag=3 the relerr 
// parameter has been adjusted appropriately for continuing the integration. In 
// the case of iflag=4 the function counter will be reset to 0 and another 3000 
// function evaluations are allowed.
//
// However,in the case iflag=5, the user must first alter the error criterion to 
// use a positive value of abserr before integration can proceed. if he does 
// not,execution is terminated.
//
// Also, in the case iflag=6, it is necessary for the user to reset iflag to 2 
// (or -2 when the one-step integration mode is being used) as well as 
// increasing either abserr,relerr or both before the integration can be 
// continued. If this is not done, execution will be terminated. the occurrence 
// of iflag=6 indicates a trouble spot (solution is changing rapidly,singularity 
// may be present) and it often is inadvisable to continue.

// If iflag=7 is encountered, the user should use the one-stepintegration mode 
// with the stepsize determined by the code or consider switching to the adams 
// codes de/step,intrp. if the user insists upon continuing the integration with 
// rkf45, he must reset iflag to 2 before calling rkf45 again. otherwise, 
// execution will be terminated.

// If iflag=8 is obtained, integration can not be continued unless the invalid 
// input parameters are corrected.

// It should be noted that the arrays work,iwork contain information required 
// for subsequent integration. accordingly, work and iwork should not be altered.

// Member variable notes
// rkf45 integrates a system of first order ordinary differential equations. The 
// arrays yp, f1, f2, f3, f4,and f5 and the variables h, savre, savae, nfe, kop,
// init, jflag,and kflag are used internally by the code and should not be 
// altered. Items of possible interest are
// yp  - derivative of solution vector at t
// h   - an appropriate stepsize to be used for the next step
// nfe - counter on the number of derivative function evaluations 

//------------------------------------------------------------------------------
// Performs an integration step. This function manages the step
//------------------------------------------------------------------------------
void hwRungeKutta::TakeStep(double& t, double tout, int& iflag)
{
    if (!m_pFunc)
    {
        m_status(HW_MATH_ERR_NULLPOINTER);
        return;
    }

    // check input parameters
    int mflag    = abs(iflag);
    int  numEqns = m_y.Size();

    // if (numEqns < 1 || relerr < 0.0 || abserr < 0.0 || mflag == 0 || mflag > 8)
    if (numEqns < 1 || mflag == 0 || mflag > 8)
    {
        m_status(HW_MATH_ERR_INVALIDINPUT);
        return;
    }


#if 0 // Commented code
    // the expense is controlled by restricting the number of function evaluations to be
    // approximately maxnfe.
    // maxnfe = 3000;     // corresponds to about 500 steps.
    // maxnfe = 15000;    // corresponds to about 2500 steps. it's not 1974 anymore. Now completely disabled.
#endif

    // check continuation possibilities
    if (t == tout && kflag != 3)
    {
        // invalid input
        m_status(HW_MATH_ERR_INVALIDINPUT);
        iflag = 8;
        return;
    }

    if (mflag != 2)
    {
        if (iflag == 3)
        {
            iflag = jflag; // reset flag value from previous call
            if (kflag == 3)
            {
                mflag = abs(iflag);
            }
        }
        else if (iflag == 4)
        {
            nfe = 0;       // reset function evaluation counter
            if (mflag != 2)
            {
                iflag = jflag;
                if (kflag == 3)
                {
                    mflag = abs(iflag);
                }
            }
        }
        else if (iflag == 5) // && abserr > 0.0)
        {
            iflag = jflag; // reset flag value from previous call
            if (kflag == 3)
            {
                mflag = abs(iflag);
            }
        }
    }
    else
    {
        // iflag = +2 or -2
        if (kflag == 3)
        {
            iflag = jflag; // reset flag value from previous call
            if (kflag == 3)
            {
                mflag = abs(iflag);
            }
        }
        else if (init == 0)
        {
            iflag = jflag;  // reset flag value from previous call
            if (kflag == 3)
            {
                mflag = abs(iflag);
            }
        }
        else if (kflag == 4)
        {
            nfe = 0; // reset function evaluation counter
            if (mflag != 2)
            {
                iflag = jflag; // reset flag value from previous call 
                if (kflag == 3)
                {
                    mflag = abs(iflag);
                }
            }
        }
        else if (kflag == 5) // && abserr == 0.0)
        {
            // integration cannot be continued
            m_status(HW_MATH_ERR_NOTCONVERGE);
            return;
        }
        else if (kflag == 6) // && relerr <= savre && abserr <= savae)
        {
            // integration cannot be continued
            m_status(HW_MATH_ERR_NOTCONVERGE);
            return;
        }
    }

    // save input iflag and set continuation flag value for subsequent input checking
    jflag = iflag;
    kflag = 0;

    // save relerr and abserr for checking input on subsequent calls 
    // savre = relerr;
    // savae = abserr;

    // remin is the minimum acceptable value of relerr. Attempts to obtain 
    // higher accuracy with this subroutine are usually very expensive and 
    // often unsuccessful.
    double remin = 1.0e-12;

    // restrict relative error tolerance to be at least as large as 2*eps+remin
    // to avoid limiting precision difficulties arising from impossible accuracy requests 
    double rer = 2.0 * eps + remin;

    if (relerr < rer)
    {
        // relative error tolerance too small
        relerr = rer * (1.0 + eps);
        iflag = 3;
        kflag = 3;
        return;
    }

    double dt = tout - t;

    // Advance y from t to tout if possible
    bool hfaild;
    bool output;
    int  k;
    double a;
    double ee;
    double eeoet;
    double esttol;
    double et;
    double hmin;
    double s;
    double scale;
    double tol;
    double toln;
    double ypk;
    double tmpreal;
    hwMatrix ae;

    // initialization --
    //      set initialization completion indicator,init
    //      set indicator for too many output points, kop
    //      evaluate initial derivatives
    //      set counter for function evaluations,nfe
    //      estimate starting stepsize

    if (mflag == 1)
    {
        init = 0;
        kop = 0;
        a = t;

        m_status = m_pFunc(a, m_y.GetRealData(), yp, m_userData);

        if (!m_status.IsOk())
        {
            iflag = 9;
            return;
        }

        nfe = 1;

        if (t == tout)
        {
            iflag = 2;
            return;
        }
    }

    if (init == 0 || t == tout)
    {
        init = 1;
        h = fabs(dt);
        toln = 0.0;

        for (k = 0; k < numEqns; k++)
        {
            tol = relerr * fabs(m_y(k)) + abserr(k);

            if (tol <= 0.0)
                continue;

            toln = tol;
            ypk = fabs(yp[k]);

            if (ypk * pow(h, 5) > tol)
            {
                h = pow(tol / ypk, 0.2);
            }
        }

        if (toln <= 0.0)
        {
            h = 0.0;
        }
        h = _max(h, u26 * _max(fabs(t), fabs(dt)));

        jflag = sign(2, iflag);
    }

    // set stepsize for integration in the direction from t to tout
    h = sign(h, dt);

    // test to see if rkf45 is being severely impacted by too many output points 
    // if (fabs(h) >= 2.0 * fabs(dt))
    if (fabs(h) > 10.0 * fabs(dt))   // it's not 1974 anymore, so relax condition
    {
       kop++;
    }

    if (kop == 100)
    {
        // unnecessary frequency of output 
        kop = 0;
        iflag = 7;
        return;
    }

    if (fabs(dt) <= u26 * fabs(t))
    {
        // if too close to output point, extrapolate and return 
        for (k = 0; k < numEqns; k++)
        {
            m_y(k) += dt * yp[k];
        }

        a = tout;
        m_pFunc(a, m_y.GetRealData(), yp, m_userData);
        nfe++;
        t = tout;
        iflag = 2;
        return;
    }

    // initialize output point indicator 
    output = false;

    // to avoid premature underflow in the error tolerance function, scale the error tolerances
    scale = 2.0 / relerr;
    ae = abserr * scale;

    // step by step integration
    do
    {
        hfaild = false;

        // set smallest allowable stepsize 
        hmin = u26 * fabs(t);

        // adjust stepsize if necessary to hit the output point. Look ahead two steps to avoid
        // drastic changes in the stepsize and thus lessen the impact of output points on the code.
        dt = tout - t;

        if (fabs(dt) < 2.0 * fabs(h))
        {
            if (fabs(dt) <= fabs(h))
            {
                // the next successful step will complete the integration to the output point
                output = true;
                h = dt;
            }
            else
            {
                h = 0.5 * dt;
            }
        }

        // core integrator for taking a single step
        /*
        the tolerances have been scaled to avoid premature underflow in computing the
        error tolerance function et. to avoid problems with zero crossings, relative
        error is measured using the average of the magnitudes of the solution at the
        beginning and end of a step. The error estimate formula has been grouped to
        control loss of significance. To distinguish the various arguments, h is not
        permitted to become smaller than 26 units of roundoff in t. Practical limits
        on the change in the stepsize are enforced to smooth the stepsize selection
        process and to avoid excessive chattering on problems having discontinuities.
        To prevent unnecessary failures, the code uses 9/10 the stepsize it estimates
        will succeed. After a step failure, the stepsize is not allowed to increase for 
        the next attempted step. this makes the code more efficient on problems having
        discontinuities and more effective in general since local extrapolation is
        being used and extra caution seems warranted.
        */

        do
        {
#if 0 // Commented code
            // test number of derivative function evaluations. If okay, try to advance
            // the integration from t to t+h
            /*
            if (nfe > maxnfe)
            {
                // too much work 
                iflag = 4;
                kflag = 4;
                m_status(HW_MATH_ERR_RK45WORKLOAD);
                return;s
            }
            */
#endif
            // advance an approximate solution over one step of length h 
            Fehld(t);
            nfe += 5;

            // compute and test allowable tolerances versus local error estimates and
            // remove scaling of tolerances. note that relative error is measured with
            // respect to the average of the magnitudes of the solution at the beginning
            // and end of the step.
            eeoet = 0.0;

            for (k = 0; k < numEqns; k++)
            {
                et = fabs(m_y(k)) + fabs(f1[k]) + ae(k);

                if (et > 0.0)
                {
                    ee = fabs((-2090.0 * yp[k] + (21970.0 * f3[k] - 15048.0 * f4[k])) + (22528.0 * f2[k] - 27360.0 * f5[k]));
                }
                else
                {
                    // inappropriate error tolerance 
                    iflag=5;
                    return;
                }

                eeoet = _max(eeoet, ee / et);
            }

            esttol = fabs(h) * eeoet * scale / 752400.0;

            if (esttol <= 1.0)
            {
                break;
            }
            // unsuccessful step
            // reduce the stepsize, try again
            // the decrease is limited to a factor of 1/10 
            hfaild = true;
            output = false;
            s = 0.1;

            if (esttol < 59049.0)
                s = 0.9 / pow(esttol, 0.2);

            h *= s;
        }
        while (fabs(h) > hmin);

        // requested error unattainable at smallest allowable stepsize 
        if (esttol > 1.0)
        {
            iflag = 6;
            kflag = 6;
            return;
        }

        // successful step 
        // store solution at t+h
        // and evaluate derivatives there 
        t += h;

        for (k = 0; k < numEqns; k++)
        {
            m_y(k) = f1[k];
        }

        a = t;
        m_pFunc(a, m_y.GetRealData(), yp, m_userData);
        nfe++;

        // choose next stepsize
        // the increase is limited to a factor of 5
        // if step failure has just occurred, next
        // stepsize is not allowed to increase
        s = 5.0;

        if (esttol > 1.889568e-4)
        {
            s = 0.9 / pow(esttol, 0.2);
        }

        tmpreal = 1.0;

        if (hfaild)
        {
            s = _min(s, tmpreal);
        }

        h = sign(_max(s * fabs(h), hmin), h);

        // end of core integrator, should we take another step?
        if (output)
        {
            break;
        }
    }
    while (iflag > 0);

    // integration successfully completed
    if (output)
    {
        // interval mode 
        t = tout;
        iflag = 2;
    }
    else
    {
        // one-step mode 
        iflag = -2;
    }
}
//------------------------------------------------------------------------------
// Fehld - fehlberg fourth-fifth order runge-kutta method
// fehl integrates a system of neqn first order ordinary differential equations 
// of the form
// dy(i)/dt=f(t,y(1),---,y(neqn))
// where the initial values y(i) and the initial derivatives yp(i) are specified 
// at the starting point t. fehl advances the solution over the fixed step h and 
// returns the fifth order (sixth order accurate locally) solution approximation 
// at t+h in array s(i). f1,---,f5 are arrays of dimension neqn which are needed 
// for internal storage. The formulas have been grouped to control loss of 
// significance. Fehl should be called with an h not smaller than 13 units of 
// roundoff in t so that the various independent arguments can be distinguished.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Performs an integration step. This is the core integrator
//------------------------------------------------------------------------------
void hwRungeKutta::Fehld(double t)
{
    int     numEqns = m_y.Size();
    double* s = f1;

    double  ch = h / 4.0;

    for (int k = 0; k < numEqns; k++)
    {
        f5[k] = m_y(k) + ch * yp[k];
    }

    m_pFunc(t + ch, f5, f1, m_userData);
    ch = 3.0 * h / 32.0;

    for (int k = 0; k < numEqns; k++)
    {
        f5[k] = m_y(k) + ch * (yp[k] + 3.0 * f1[k]);
    }

    m_pFunc(t + 3.0 * h / 8.0, f5, f2, m_userData);
    ch = h / 2197.0;

    for (int k = 0; k < numEqns; k++)
    {
        f5[k] = m_y(k) + ch * (1932.0 * yp[k] + (7296.0 * f2[k] - 7200.0 * f1[k]));
    }

    m_pFunc(t + 12.0 * h / 13.0, f5, f3, m_userData); 
    ch = h / 4104.0;

    for (int k = 0; k < numEqns; k++)
    {
        f5[k] = m_y(k) + ch * ((8341.0 * yp[k] - 845.0 * f3[k]) + (29440.0 * f2[k] - 32832.0 * f1[k]));
    }

    m_pFunc(t + h, f5, f4, m_userData);
    ch = h / 20520.0;

    for (int k = 0; k < numEqns; k++)
    {
        f1[k] = m_y(k) + ch * ((-6080.0 * yp[k] + (9295.0 * f3[k] - 5643.0 * f4[k])) + (41040.0 * f1[k] - 28352.0 * f2[k]));
    }
    m_pFunc(t + h / 2.0, f1, f5, m_userData);

    // compute approximate solution at t+h 
    ch = h / 7618050.0;

    for (int k = 0; k < numEqns; k++)
    {
        s[k] = m_y(k) + ch * ((902880.0 * yp[k] + (3855735.0 * f3[k] - 1371249.0 * f4[k])) + (3953664.0 * f2[k]+ 277020.0 * f5[k]));
    }
}
//------------------------------------------------------------------------------
// Returns true if successful
//------------------------------------------------------------------------------
bool hwRungeKutta::Success(int flag)
{
    return (flag == 2);
}
//------------------------------------------------------------------------------
// Returns true if execution can continue
//------------------------------------------------------------------------------
bool hwRungeKutta::Continue(int flag)
{
    if (flag != -2)
    {
        if (flag == 7)
        {
            m_status(HW_MATH_ERR_RK45STEPSIZE);
        }
        return false;
    }

    return true;
}
