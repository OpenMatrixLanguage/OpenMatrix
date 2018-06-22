/**
* @file hwRungeKutta.h
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
#ifndef _DiffEq_RungeKutta_h
#define _DiffEq_RungeKutta_h

#include "hwDiffEqSolver.h"

//!
//! \typedef hwMathStatus
//!
typedef hwMathStatus (*RKF45Fn_client)(double          t, 
                                       const double*   y, 
                                       double*         yp,
                                       const hwMatrix* userData);

//------------------------------------------------------------------------------
//!
//! \class hwRungeKutta
//! \brief Runge-Kutta-Fehlberg 45 ODE solver class
//!
//------------------------------------------------------------------------------
class hwRungeKutta : public hwDiffEqSolver
{
public:

    //!
    //! Constructor
    //! \param sysfunc
    //! \param y
    //! \param relerr_  Optional argument for relative error
    //! \param abserr_  Optional argument for absolute error
    //! \param userData Optional argument for user data
    //!
    hwRungeKutta(RKF45Fn_client  sysfunc, 
                 const hwMatrix& y,
                 double          relerr_  = 1.0e-3, 
                 const hwMatrix* abserr_  = nullptr,
                 const hwMatrix* userData = nullptr);
    //!
    //! Destructor
    //!
    virtual ~hwRungeKutta();

    //!
    //! Performs an integration step. This function manages the step
    //! \param t
    //! \param tout
    //! \param iflag
    //!
    void TakeStep(double& t, 
                  double  tout, 
                  int&    iflag);
    //!
    //! Sets the status flag
    //! \param flag A value of 2 indicates success
    //!
    bool Success(int flag);
    //!
    //! Assess continuation flag
    //! \param flag
    //!
    bool Continue(int flag);
    //!
    //! Returns the relative error
    //!
    double RelativeError() { return relerr; }

private:
    int             nfe;         //!< number of function evaluations
    int             kop;         //!< operation counter for excessive outputs
    int             init;        //!< initialization flag
    int             jflag;       //!< internal flag
    int             kflag;       //!< internal flag
    hwMatrix        abserr;      //!< absolute error
    double          relerr;      //!< relative error
    double          h;           //!< step size
    double          eps;         //!< precision variable
    double          u26;         //!< 26 * eps
    double*         work;        //!< work array
    double*         yp;          //!< solution derivative vector
    double*         f1;          //!< 1st function output vector
    double*         f2;          //!< 2nd function output vector
    double*         f3;          //!< 3rd function output vector
    double*         f4;          //!< 4th function output vector
    double*         f5;          //!< 5th function output vector
    const hwMatrix* m_userData;  //!< user data matrix
    RKF45Fn_client    m_pFunc;     //!< function pointer for ode system

protected:
    //!
    //! Performs an integration step. This is the core integrator
    //! \param t
    //! 
    void Fehld(double t);
};

#endif // _DiffEq_RungeKutta_h
