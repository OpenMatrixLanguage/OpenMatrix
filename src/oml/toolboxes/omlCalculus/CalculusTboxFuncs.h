/**
* @file CalculusTboxFuncs.h
* @date January 2015
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

#ifndef __CALCULUSTBOXFUNCS_OML_H__
#define __CALCULUSTBOXFUNCS_OML_H__

#pragma warning(disable: 4251)

#include "EvaluatorInt.h"
#include "CalculusTboxDefs.h"

//------------------------------------------------------------------------------
//!
//! \brief oml Calculus functions
//!
//------------------------------------------------------------------------------

extern "C" 
{
    //!
    //! Entry point which registers oml Calculus functions with oml
    //! \param eval Evaluator interface
    //!
    CALCULUSOMLTBOX_DECLS int InitDll(EvaluatorInterface eval);
}

//! 
//! Returns true after performing numerical integration (quadrature)
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlQuad(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//! 
//! Returns true after performing numerical integration using Simpson's rule
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlQuadv(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs);
//! 
//! Performs numerical integration  of discrete data using the trapezoid rule
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlTrapz(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs);
//! 
//! Cumulative numerical integration  of discrete data using the trapezoid rule
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlCumtrapz(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs);

#endif // __CALCULUSTBOXFUNCS_OML_H__