/**
* @file BetterCalc.h
* @date June 2015
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

#ifndef __BETTERCALC_H__
#define __BETTERCALC_H__

class ConsoleWrapper;

//!
//! Handles user interrupt returns true if application needs to quit
//! \param omlWrapper Oml interpreter wrapper
//!
bool HandleUserInterrupt(ConsoleWrapper* omlWrapper);
//!
//! Sets user interrupt
//!
void SetUserInterrupt();
//!
//! Sets user interrupt handler - hooks to control C
//!
void SetUserInterruptHandler();
//!
//! Prints banner
//!
void PrintBanner();
//!
//! Gets input command
//! \param interpwrapper Oml interp wrapper
//!
std::string GetInputCommand(ConsoleWrapper* interpwrapper);
//!
//! Clears command window (clc command)
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool hml_clc(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs);
//!
//! Encrypts given file (encrypt command)
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool oml_encryptfile(EvaluatorInterface           eval, 
                     const std::vector<Currency>& inputs, 
                     std::vector<Currency>&       outputs);

//!
//! Prints help
//!
void OML_help();
//!
//! Gets the version string of the application
//! \param appdir Application directory
//!
std::string GetVersion(const std::string& appdir);
//! Returns true if successful in getting the version string (version command)
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlVersion(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs);
//! 
//! Registers oml built in functions
//! 
void RegisterBuiltInFuncs();
//!
//! Gets argc (getargc command)
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlGetArgC(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs);
//!
//! Gets argv at the given index (getargv command)
//! \param eval    Evaluator interface
//! \param inputs  Vector of inputs
//! \param outputs Vector of outputs
//!
bool OmlGetArgV(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs);

#endif

// End of file:

