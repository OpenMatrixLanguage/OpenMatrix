/**
* @file MatioTboxFuncs.h
* @date November 2015
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

#ifndef __MATIOTBOXFUNCS_H__
#define __MATIOTBOXFUNCS_H__

#pragma warning(disable: 4251)

#include "EvaluatorInt.h"
#include "MatioTboxDefs.h"

//------------------------------------------------------------------------------
//!
//! \brief oml load/save functions which use the MATIO library
//!
//------------------------------------------------------------------------------

extern "C" 
{
    //!
    //! Entry point which registers load/save functions with oml
    //! \param eval Evaluator interface
    //!
    MATIOHML2TBOX_DECLS int InitDll(EvaluatorInterface eval);
}

//!
//! Returns true after loading the given file using MATIO library [load command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlLoad(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Returns true after saving the given file using MATIO library [save command]
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlSave(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);

#endif // __MATIOTBOXFUNCS_H__