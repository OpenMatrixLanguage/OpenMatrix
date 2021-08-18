/**
* @file hwIdaWrap.h
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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
#ifndef _DiffEq_IDA_h
#define _DiffEq_IDA_h

#include "hwDiffEqSolver.h"

#include "ida/ida.h"
#include "sundials/sundials_dense.h"   // access dense SUNMatrix
#include "sunlinsol/sunlinsol_dense.h" // access dense SUNLinearSolver

//!
//! \typedef IDAResFn_client
//!
typedef int (*IDAResFn_client)(double  t, 
                               double* y, 
                               double* yp,
                               double* r, 
                               void*   res_data);
//!
//! \typedef IDARootFn_client
//!
typedef int (*IDARootFn_client)(double  t, 
                                double* y, 
                                double* yp,
                                double* gout, 
                                void*   g_data);
//!
//! \typedef IDADenseJacFn_client
//!
typedef int (*IDADenseJacFn_client)(long int N, 
                                    double   t, 
                                    double   cj,
                                    double*  y, 
                                    double*  yp, 
                                    double*  r, 
                                    void*    jac_data, 
                                    double** Jac);

#include "DiffEqExports.h"
DIFFEQ_DECLS void WriteEventDataIDA(double* isterminal, double* direction, int nrtfn);

//------------------------------------------------------------------------------
//! 
//! \class hwIdaWrap (IDA DAE solver)
//! \brief IDA wrapper class
//!
//------------------------------------------------------------------------------
class hwIdaWrap : public hwDiffEqSolver
{
public:
    //! Constructor
    //! \param sysfunc
    //! \param rootfunc
    //! \param nrtfn
    //! \param jacDfunc
    //! \param tin
    //! \param y_
    //! \param yp_
    //! \param job         control code for IDA options // todo: improve scheme
    //! \param reltol      Real tolerance
    //! \param abstol      Absolute tolerance
    //! \param userData    Optional user data
    //! \param pEventTime  Event function times
    //! \param pEventFnVal Event function values
    //! \param pEventIndx  Event function indices
    hwIdaWrap(IDAResFn_client      sysfunc, 
              IDARootFn_client     rootfunc,
              int                  nrtfn,
              IDADenseJacFn_client jacDfunc,
              double               tin,
              const hwMatrix&      y_, 
              const hwMatrix&      yp_, 
              const char*          job,
              double               reltol   = 0.001,
              const hwMatrix*      abstol   = nullptr,
              double               maxstep  = -999.0,
              const hwMatrix*      userData = nullptr,
              hwMatrix*            pEventTime  = nullptr,
              hwMatrix*            pEventFnVal = nullptr,
              hwMatrix*            pEventIndx  = nullptr);

    //!
    //! Destructor
    //!
    virtual ~hwIdaWrap();

    //!
    //! Performs an integration step
    //! \param t
    //! \param tout
    //! \param flag
    //!
    void TakeStep(double& t, 
                  double  tout, 
                  int&    flag);
    //!
    //! Returns true if successful
    //! \param flag
    //!
    bool Success(int flag);
    //!
    //! Returns true if execution can continue
    //! \param flag
    //!
    bool Continue(int flag);
    //!
    //! Respond to events that have returned zeros
    //! \param flag
    //!
    int ManageEvents(double t, bool init);

private:
    hwMatrix        m_yp;          //!< output derivative function vector
    void*           ida_mem;       //!< SUNDIALS IDA internal memory
    N_Vector        y;             //!< SUNDIALS DAE output vector
    N_Vector        yp;            //!< SUNDIALS derivative of DAE output
    SUNMatrix       A;             //!< SUNDIALS dense matrix
    SUNLinearSolver LS;            //!< SUNDIALS linear solver memory
    int             m_nrtfn;       //!< SUNDIALS number of events
    hwMatrix*       m_pEventTime;  //!< Event function times
    hwMatrix*       m_pEventFnVal; //!< Event function values
    hwMatrix*       m_pEventIndx;  //!< Event function indices

    //!
    //! Set stop time for one step mode
    //! \param tstop Stop time
    //!
    void SetStopTime(double tstop);
    //!
    //! Check IDA flag
    //! \param flagvalue
    //! \param funcname  Function name
    //! \param opt       Function return value
    //!
int Check_flag(void*       flagvalue, 
               const char* funcname, 
               int         opt);
};

#endif // _DiffEq_IDA_h
