/**
* @file OmlOptimizationTboxFuncs.cxx
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

#include "OmlOptimizationTboxFuncs.h"

#include <cassert>

#include "hwMatrix.h"

#include "BuiltInFuncsUtils.h"
#include "OML_Error.h"
#include "StructData.h"

#define OPTIM "Optimization"

//------------------------------------------------------------------------------
// Entry point which registers optimization functions with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("fminbnd",     &OmlFminbnd,     FunctionMetaData(-3, 2, OPTIM));
    eval.RegisterBuiltInFunction("fminsearch",  &OmlFminsearch,  FunctionMetaData(-3, 2, OPTIM));
    eval.RegisterBuiltInFunction("fminunc",     &OmlFminunc,     FunctionMetaData(-3, 2, OPTIM));
    eval.RegisterBuiltInFunction("fsolve",      &OmlFsolve,      FunctionMetaData(-3, 2, OPTIM));
    eval.RegisterBuiltInFunction("fzero",       &OmlFzero,       FunctionMetaData(-3, 2, OPTIM));
    eval.RegisterBuiltInFunction("lsqcurvefit", &OmlLsqcurvefit, FunctionMetaData(-5, 2, OPTIM));
    eval.RegisterBuiltInFunction("optimset",    &OmlOptimset,    FunctionMetaData(-1, 1, OPTIM));

    //\todo: Uncomment when complete
    //eval.RegisterBuiltInFunction("optimoptions", &OmlOptimoptions, FunctionMetaData(-1, -1, OPTIM));
    return 1;
}
//------------------------------------------------------------------------------
// Sets options for optimization functions
//------------------------------------------------------------------------------
bool OmlOptimset(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)
{
    int nargin = static_cast<int>(inputs.size());
    
    if (nargin == 0)  // Just prints options
    {
        std::string msg ("Available optimization options:\n\n");
        msg += "Display\n";
        msg += "GradConstr\n";
        msg += "GradObj\n";
        msg += "Jacobian\n";
        msg += "MaxFunEvals\n";
        msg += "MaxIter\n";
        msg += "TolCon\n";
        msg += "TolFun\n";
        msg += "TolKKT\n";
        msg += "TolX";

        Currency tmp(msg);
        tmp.DispOutput();
        eval.PrintResult(tmp);

        if (eval.GetNargoutValue() <= 0)
            return true;

        Currency out = EvaluatorInterface::allocateStruct();
        assert(out.Struct());
        out.Struct()->SetValue(0, -1, "Display",     "off");
        out.Struct()->SetValue(0, -1, "GradConstr",  "off");
        out.Struct()->SetValue(0, -1, "GradObj",     "off");
        out.Struct()->SetValue(0, -1, "Jacobian",    "off");
        out.Struct()->SetValue(0, -1, "MaxFunEvals", 400);
        out.Struct()->SetValue(0, -1, "MaxIter",     100);
        out.Struct()->SetValue(0, -1, "TolCon",      0.5);
        out.Struct()->SetValue(0, -1, "TolFun",      1.0e-7);
        out.Struct()->SetValue(0, -1, "TolKKT",      1.0e-4);
        out.Struct()->SetValue(0, -1, "TolX",        EvaluatorInterface::allocateMatrix());
        outputs.push_back(out);
        return true;
    }

    if (nargin % 2 != 0)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency outStruct = EvaluatorInterface::allocateStruct();

    for (int i = 0; i < nargin; i += 2)
    {
        if (!inputs[i].IsString())
            throw OML_Error(OML_ERR_STRING, i + 1, OML_VAR_PARAMETER);

        std::string opt (inputs[i].StringVal());

        assert(i + 1 < nargin);
        Currency cur (inputs[i + 1]); 

        if (opt == "MaxIter")  //\todo: Should this not be integer?
        {
            if (!cur.IsScalar())
                throw OML_Error(OML_ERR_SCALAR, i+2, OML_VAR_MAXITER);
        }
        else if (opt == "MaxFunEvals") //\todo: Should this not be integer?
        {
            if (!cur.IsScalar())
                throw OML_Error(OML_ERR_SCALAR, i+2, OML_VAR_MAXFUNEVALS);
        }
        else if (opt == "TolFun")
        {
            if (!cur.IsScalar())
                throw OML_Error(OML_ERR_SCALAR, i+2, OML_VAR_TOLFUN);
        }
        else if (opt == "TolX")
        {
            if (!cur.IsScalar())
                throw OML_Error(OML_ERR_SCALAR, i+2, OML_VAR_TOLFUN);
        }
        else if (opt == "TolCon")
        {
            if (!cur.IsScalar())
                throw OML_Error(OML_ERR_SCALAR, i+2, OML_VAR_TOLCON);
        }
        else if (opt == "TolKKT")
        {
            if (!cur.IsScalar())
                throw OML_Error(OML_ERR_SCALAR, i+2, OML_VAR_TOLKKT);
        }
        else if (opt == "GradObj" || opt == "GradConstr" || opt == "Jacobian" ||
                 opt == "Hessian")
        {
            if (!cur.IsString())
                throw OML_Error(OML_ERR_STRING, i);

            std::string val (cur.StringVal());
            if (val != "on" && val != "off")
            {
                std::string message = "Error: argument " + 
                    std::to_string(static_cast<long long>(i + 2)) + " must be \"on\" or \"off\"";
                throw OML_Error(message);
            }
        }
        else if (opt == "Display")
        {
            if (!cur.IsString())
                throw OML_Error(OML_ERR_STRING, i);

            std::string val (cur.StringVal());
            if (val != "iter" && val != "off")
            {
                std::string msg = "Error: argument " + 
                    std::to_string(static_cast<long long>(i + 2)) + 
                    " must be \"iter\" or \"off\"";
                throw OML_Error(msg);
            }
        }
        else
        {
            std::string msg = "Error: parameter " + opt + "[";
            msg += std::to_string(static_cast<long long>(i + 2)) + "] is unsupported";
            throw OML_Error(msg);
        }

        outStruct.Struct()->SetValue(0, -1, opt, cur);
    }

    outputs.push_back(outStruct);
    return true;
}
//------------------------------------------------------------------------------
// Sets options [optimoptions] - //\todo: To be completed and registered
//------------------------------------------------------------------------------
bool OmlOptimoptions(EvaluatorInterface           eval, 
                     const std::vector<Currency>& inputs, 
                     std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    if (nargin == 0)
    {
        std::string msg ("Available optimization options:\n");
        msg += "TolCon\n";
        msg += "TolFun\n";        
        msg += "Display\n";
        eval.PrintResult(msg);

        if (eval.GetNargoutValue() > 0)
        {
            outputs.push_back(EvaluatorInterface::allocateMatrix());
        }
        return true;
    }

    if (nargin % 2 != 0)
        throw OML_Error(OML_ERR_NUMARGIN);

    //\todo: Complete implementation
    return true;
}
