/**
* @file OmlDiffEqTboxFuncs.cxx
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

#include "OmlDiffEqTboxFuncs.h"

#include <cassert>

#include "OML_Error.h"
#include "StructData.h"

#define DIFF "DifferentialEquations"

//------------------------------------------------------------------------------
// Entry point which registers differential equations functions with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("odeset", &OmlOdeset, 
                                 FunctionMetaData(-1,-1, DIFF));
    eval.RegisterBuiltInFunction("ode45", &OmlOde45, 
                                 FunctionMetaData(4,2, DIFF));
    eval.RegisterBuiltInFunction("ode113", &OmlOde113,
                                 FunctionMetaData(4,2, DIFF));
    eval.RegisterBuiltInFunction("ode15s", &OmlOde15s,
                                 FunctionMetaData(3,2, DIFF));
    eval.RegisterBuiltInFunction("ode15i", &OmlOde15i,
                                 FunctionMetaData(3,2, DIFF));
    return 1;
}
//------------------------------------------------------------------------------
// Sets ordinary differential equation options [odeset]
//------------------------------------------------------------------------------
bool OmlOdeset(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    size_t nargin = (inputs.empty()) ? 0 : inputs.size();

    if (nargin == 0)
    {
        std::string msg = "Available ODE options:\n";
        msg += "RelTol\n";
        msg += "AbsTol\n";
        msg += "Jacobian\n";

        Currency cmsg (msg);
        cmsg.DispOutput();
        eval.PrintResult(cmsg);

        if (eval.GetNargoutValue() > 0)
        { 
            Currency out (EvaluatorInterface::allocateStruct());
            assert(out.Struct());
            out.Struct()->SetValue(0, -1, "RelTol",   1.0e-3);
            out.Struct()->SetValue(0, -1, "AbsTol",   1.0e-6);
            out.Struct()->SetValue(0, -1, "Jacobian", Currency());
            outputs.push_back(out);
        }
        return true;
    }

    if (nargin % 2 != 0)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    Currency out = EvaluatorInterface::allocateStruct();

    for (int i = 0; i < nargin; i += 2)
    {
        Currency cur = inputs[i];
        if (!cur.IsString())
        {
            throw OML_Error(OML_ERR_STRING, i + 1);
        }

        std::string val (cur.StringVal());
        if (val == "RelTol")
        {
            if (!inputs[i+1].IsScalar())
            {
                throw OML_Error(OML_ERR_SCALAR, i + 2, OML_VAR_RELTOL);
            }
        }
        else if (val == "AbsTol")
        {
            if (!inputs[i+1].IsScalar() && !inputs[i+1].IsMatrix())
            {
                throw OML_Error(OML_ERR_SCALARVECTOR, i + 2, OML_VAR_ABSTOL);
            }
        }
        else if (val == "Jacobian")
        {
            if (!inputs[i+1].IsFunctionHandle())
            {
                throw OML_Error(OML_ERR_HANDLE, i + 2, OML_VAR_JACOBIAN);
            }
        }
        else
        {
            std::string strnum  = std::to_string(static_cast<long long>(i + 1));
            std::string message = "Error: unsupported parameter: " + val +
                                  " [" + strnum + "]";
            throw OML_Error(message);
        }

        out.Struct()->SetValue(0, -1, val, inputs[i+1]);
    }

    outputs.push_back(out);
    return true;
}
