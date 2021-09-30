/**
* @file MathUtilsTboxFuncs.h
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

#ifndef __MATHUTILSTBOXFUNCS_H__
#define __MATHUTILSTBOXFUNCS_H__

#pragma warning(disable: 4251)

#include "EvaluatorInt.h"
#include "MathUtilsTboxDefs.h"
#include <map>

//------------------------------------------------------------------------------
//!
//! \brief omlMathUtils toolbox functions
//!
//------------------------------------------------------------------------------

extern "C" 
{
    //!
    //! Entry point which registers toolbox with oml
    //! \param eval Evaluator interface
    //!
    MATHUTILSOMLTBOX_DECLS int InitDll(EvaluatorInterface eval);
    //!
    //! Returns toolbox version
    //! \param eval Evaluator interface
    //!
    MATHUTILSOMLTBOX_DECLS double GetToolboxVersion(EvaluatorInterface eval);
}

//!
//! Returns true after executing beta function
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlBeta(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Returns true after executing log beta function
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlBetaLn(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs);
//!
//! Returns true after executing gamma function
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlGamma(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs);
//!
//! Returns true after executing log gamma function
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlGammaLn(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs);
//!
//! Returns true after executing gamma function
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlFactorial(EvaluatorInterface           eval,
                  const std::vector<Currency>& inputs,
                  std::vector<Currency>&       outputs);
//!
//! Returns true after dividing range of data into given number of equal bins
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlBins(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Returns true after computing rational fraction approximation
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlRat(EvaluatorInterface           eval, 
            const std::vector<Currency>& inputs, 
            std::vector<Currency>&       outputs);

#endif // __MATHUTILSTBOXFUNCS_H__