/**
* @file OmlOptimizationTboxFuncs.h
* @date September 2017
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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

#ifndef __OMLOPTIMIZATIONTBOXFUNCS_H__
#define __OMLOPTIMIZATIONTBOXFUNCS_H__

#pragma warning(disable: 4251)

#include "EvaluatorInt.h"

#include "OmlOptimizationTboxDefs.h"

//------------------------------------------------------------------------------
//!
//! \brief oml optimization functions
//!
//------------------------------------------------------------------------------

//!
//! Entry point which registers optimization functions with oml
//! \param eval Evaluator interface
//!
extern "C" OMLOPTIMIZATIONTBOX_DECLS int InitDll(EvaluatorInterface eval);
//!
//! Returns toolbox version
//! \param eval Evaluator interface
//!
extern "C" OMLOPTIMIZATIONTBOX_DECLS double GetToolboxVersion(EvaluatorInterface eval);

//!
//! Finds the minimum of a univariate real function within an interval [fminbnd]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFminbnd(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs);
//!
//! Finds the unconstrained minimum of a real function using the Nelder-Mead 
//! simplex algorithm [fminsearch]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFminsearch(EvaluatorInterface           eval, 
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs);
//!
//! Finds the unconstrained minimum of a real function [fminunc]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFminunc(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs);
//!
//! Finds a solution of a system of real, nonlinear equations [fsolve]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFsolve(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs);
//!
//! Finds root of an univariate real function within an interval [fzero]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFzero(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs);
//!
//! Finds the equation parameters that produce the least squares best fit to a 
//! data set [lsqcurvefit]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlLsqcurvefit(EvaluatorInterface           eval, 
                    const std::vector<Currency>& inputs, 
                    std::vector<Currency>&       outputs);
//!
//! Sets options for optimization functions [optimset]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlOptimset(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs);
//!
//! Sets options for optimization functions [optimoptions]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlOptimoptions(EvaluatorInterface           eval, 
                     const std::vector<Currency>& inputs, 
                     std::vector<Currency>&       outputs);


#endif // __OMLOPTIMIZATIONTBOXFUNCS_H__
