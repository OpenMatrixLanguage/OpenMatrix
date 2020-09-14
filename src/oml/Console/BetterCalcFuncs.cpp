/**
* @file BetterCalcFuncs.cpp
* @date June 2015
* Copyright (C) 2015-2020 Altair Engineering, Inc.  
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

// Begin defines/includes

#include "../Runtime/BuiltInFuncsCore.h"
#include "../Runtime/BuiltInFuncsUtils.h"
#include "../Runtime/CurrencyDisplay.h"
#include "../Runtime/Interpreter.h"
#include "../Runtime/StructData.h"
#include "../Runtime/ErrorInfo.h"
#include "../Runtime/OML_Error.h"

#include "BetterCalc.h"
#include "ConsoleWrapper.h"

#include <cassert>

extern Interpreter*    interp;
extern ConsoleWrapper* wrapper;

#define OML_PRODUCT "OpenMatrix "
#include "OpenMatrix_Version.h"
// End defines/includes

//------------------------------------------------------------------------------
// Prints help
//------------------------------------------------------------------------------
void OML_help()
{
	std::cout << OML_PRODUCT << " Console Help:" << std::endl;
	std::cout << "-e \"any valid OML command\" --> execute the OML command in " << OML_PRODUCT << " and then quit." << std::endl;
	std::cout << "-f foo.oml --> load and execute the OML script called foo.oml in " << OML_PRODUCT << " and then quit" << std::endl;
	std::cout << "-help      --> give the list of supported optional arguments" << std::endl;
	std::cout << "-version   --> give the version of " << OML_PRODUCT << std::endl;
}
//------------------------------------------------------------------------------
//! Gets the version string of the application
//------------------------------------------------------------------------------
std::string GetVersion(const std::string& appdir)
{
    std::string version (OML_PRODUCT);
    version += " Version ";
    version += OML_VERSION;

    return version;
}
//------------------------------------------------------------------------------
//! Returns true if successful in getting the version string
//! \param[in]  eval    Evaluator interface
//! \param[in]  inputs  Vector of inputs
//! \param[out] outputs Vector of outputs
//------------------------------------------------------------------------------
bool OmlVersion(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
	std::string version = GetVersion(eval.GetApplicationDir());
	Currency out (version);
    out.DispOutput();

    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Gets argc (getargc command)
//------------------------------------------------------------------------------
bool OmlGetArgC(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    assert(wrapper);

    outputs.push_back(wrapper->GetArgc());
    return true;
}
//------------------------------------------------------------------------------
// Gets argv at the given index (getargv command)
//------------------------------------------------------------------------------
bool OmlGetArgV(EvaluatorInterface           eval,
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    assert(wrapper);

    if (inputs.size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    if (!inputs[0].IsPositiveInteger())
    {
        throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_VALUE);
    }
    int idx  = static_cast<int>(inputs[0].Scalar()) - 1;
    int argc = wrapper->GetArgc();
    if (idx > argc)
    {
        std::string msg ("Error: invalid input in argument 1; value must be ");
        if (argc == 1)
        {
            msg += "1";
        }
        else
        {
            msg += "in the range of 1 and " + std::to_string(static_cast<long long>(argc));
        }
        throw OML_Error(msg);
    }
    outputs.push_back(wrapper->GetArgv(idx));
    return true;
}

