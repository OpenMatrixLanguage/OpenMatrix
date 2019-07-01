/**
* @file PolynomialTboxFuncs.h
* @date January 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#ifndef __POLYNOMIALTBOXFUNCS_H__
#define __POLYNOMIALTBOXFUNCS_H__

#pragma warning(disable: 4251)

#include "EvaluatorInt.h"
#include "PolynomialTboxDefs.h"

//------------------------------------------------------------------------------
//!
//! \brief oml Polynomial functions
//!
//------------------------------------------------------------------------------

extern "C" 
{
    //!
    //! Entry point which registers polynomial functions with oml
    //! \param eval Evaluator interface
    //!
    POLYNOMIALOMLTBOX_DECLS int InitDll(EvaluatorInterface eval);
    //!
    //! Returns toolbox version
    //! \param eval Evaluator interface
    //!
    POLYNOMIALOMLTBOX_DECLS double GetToolboxVersion(EvaluatorInterface eval);
}

//!
//! Returns true after computing the roots of a polynomial
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlRoots(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs);
//!
//! Returns true after interpolating (x,y) data with a cubic spline
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlSpline(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs);
//!
//! Returns true after executing one-dimensional interpolation
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlInterp1(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs);
//!
//! Returns true after executing two-dimensional interpolation
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlInterp2(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs);
//!
//! Returns true after executing two-dimensional interpolation
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlIsoNormals(EvaluatorInterface           eval,
                   const std::vector<Currency>& inputs,
                   std::vector<Currency>&       outputs);
//!
//! Returns true and the derivative of the polynomial
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlDeconv(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs);
//!
//! Returns true and the derivative of the polynomial
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlPolyder(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs);
//!
//! Returns true and the integral of a polynomial 
//! coefficients are given by the input vector
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlPolyint(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs);

//!
//! Helper method for interp1 command
//! \param eval      Evaluator interface
//! \param input     Inputs
//! \param extrap    True if extrap is set
//! \param method    linear/spline methods
//! \param setExtrap True if extrap needs to be set
//! \param setMethod True if method needs to be set
//!
void interpOptionsHelper(EvaluatorInterface& eval, 
                         const Currency&     input, 
                         bool&               extrap, 
                         std::string&        method, 
                         bool&               setExtrap,
                         bool&               setMethod);

#endif // __POLYNOMIALTBOXFUNCS_H__