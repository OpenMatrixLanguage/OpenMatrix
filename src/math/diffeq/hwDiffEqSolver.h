/**
* @file hwDiffEqSolver.h
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
#ifndef _DiffEq_DiffEqSolver_h
#define _DiffEq_DiffEqSolver_h

#include "DiffEqExports.h"
#include "hwMatrix.h"

//------------------------------------------------------------------------------
//! 
//! \class hwDiffEqSolver
//! \brief Differential equation solver base class
//!
//------------------------------------------------------------------------------
class hwDiffEqSolver
{
public:
    //!
    //! Constructor
    //! \param y
    //!
    hwDiffEqSolver(const hwMatrix& y);
    //!
    //! Destructor
    //!
    virtual ~hwDiffEqSolver();

    //!
    //! Gets status
    //!
    const hwMathStatus& GetStatus() const { return m_status; }
    //!
    //! Performs an integration step
    //! \param t
    //! \param tout
    //! \param flag
    //!
    virtual void TakeStep(double& t, 
                          double  tout, 
                          int&    flag) = 0;
    //! 
    //! Returns true if successful
    //! \param flag
    //!
    virtual bool Success(int flag) = 0;
    //! 
    //! Returns true if execution can continue
    //! \param flag
    virtual bool Continue(int flag) = 0;
    //!
    //! Fills output matrix with results from each step
    //! \param tStart
    //! \param tStop
    //! \param numTimes
    //! \param ySolution Output matrix
    //!
    hwMathStatus FillMatrix(double    tStart, 
                            double    tStop,
                            int       numTimes,
                            hwMatrix& ySolution);
    //!
    //! Fills output matrix with results from each step
    //! \param time
    //! \param timeSolution
    //! \param ySolution    Output matrix
    //!
    hwMathStatus FillMatrix(const hwMatrix& time, 
                            hwMatrix*       timeSolution,
                            hwMatrix&       ySolution);
protected:
    hwMathStatus m_status;      //!< Math status
    bool         OneStepMode;   //!< Flag indicating one step mode
    hwMatrix     m_y;           //!< Solution function vector

private:
    //! 
    //! Sets stop time
    //! \param tstop
    //!
    virtual void SetStopTime(double tstop) {};
};

#endif // _DiffEq_DiffEqSolver_h
