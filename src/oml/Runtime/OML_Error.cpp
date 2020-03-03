/**
* @file OML_Error.cpp
* @date November 2015
* Copyright (C) 2015-2019 Altair Engineering, Inc.  
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

#include "OML_Error.h"

#include "EvaluatorInt.h"

// Note: each OML_MSG_ITEM is paired with OML_ERR_ITEM in omlMathErrCode below. Likewise,
// each  OML_STR_QUANTITY is paired with OML_VAR_QUANTITY in omlMathVarCode. Here
// are some options for how to throw an error.
// 1. throw OML_Error(OML_ERR_ITEM);
// 2. throw OML_Error(OML_ERR_ITEM, arg_num);
// 3. throw OML_Error(OML_ERR_ITEM, arg_num, OML_VAR_QUANTITY);
// See the OML_Error constructors for all options. When an argument number is to be included
// in the message, OML_MSG_ITEM should be written in the form "Error: problem; solution"
// Error message definitions
#define OML_MSG_NUMARGIN      "Error: invalid function call; incorrect number of input arguments"
#define OML_MSG_NUMARGOUT     "Error: invalid function call; incorrect number of output arguments"
#define OML_MSG_NUMARGINOUT   "Error: invalid function call; incorrect combination of input/output arguments"
#define OML_MSG_CELL          "Error: invalid input; must be cell"
#define OML_MSG_CELLARRAY     "Error: invalid input; must be cell array"
#define OML_MSG_STRUCT        "Error: invalid input; must be struct"
#define OML_MSG_STRING        "Error: invalid input; must be string"
#define OML_MSG_INTSTRING     "Error: invalid input; must be integer or string"
#define OML_MSG_SCALARSTRING  "Error: invalid input; must be scalar or string"
#define OML_MSG_BAD_STRING    "Error: unsupported option; see help for valid options"
#define OML_MSG_NUMERIC       "Error: invalid input; must be numeric"
#define OML_MSG_SCALAR        "Error: invalid input; must be a scalar"
#define OML_MSG_VECTOR        "Error: invalid input; must be a vector"
#define OML_MSG_VECTOR2       "Error: invalid input; must be 2 element vector"
#define OML_MSG_SCALARVECTOR  "Error: invalid input; must be a scalar or vector"
#define OML_MSG_SCALARMATRIX  "Error: invalid input; must be a scalar or matrix"
#define OML_MSG_SCALARCOMPLEXMTX "Error: invalid input; must be a scalar, complex or matrix"
#define OML_MSG_STRINGVECTOR  "Error: invalid input; must be a string or vector"
#define OML_MSG_INTVECSTR     "Error: invalid input; must be an integer, vector, or string"
#define OML_MSG_REALVECTOR    "Error: invalid input; must be a real vector"
#define OML_MSG_NNINTVECTOR   "Error: invalid input; must be nonnegative integer vector"
#define OML_MSG_POSINTVECTOR  "Error: invalid input; must be positive integer vector"
#define OML_MSG_MATRIX        "Error: invalid input; must be a matrix"
#define OML_MSG_REALMATRIX    "Error: invalid input; must be a real matrix"
#define OML_MSG_EMPTYMATRIX   "Error: invalid input; must be empty [] matrix"
#define OML_MSG_HANDLE        "Error: invalid input; must be function handle"
#define OML_MSG_HANDLE_EMPTY  "Error: invalid input; must be function handle or []"
#define OML_MSG_FUNCNAME      "Error: invalid input; function not found"
#define OML_MSG_ACCUMFUNC     "Error: invalid accumulator; must have vector input and return scalar or complex"
#define OML_MSG_REAL          "Error: invalid input; must be real"
#define OML_MSG_INTEGER       "Error: invalid input; must be integer"
#define OML_MSG_NATURALNUM    "Error: invalid input; must be nonnegative integer"
#define OML_MSG_POSINTEGER    "Error: invalid input; must be positive integer"
#define OML_MSG_FINITE        "Error: invalid value; must be finite"
#define OML_MSG_VECLENDIM     "Error: invalid vector length in specified dimension"
#define OML_MSG_ARRAYSIZE     "Error: incompatible matrices; dimensions must be consistent"
#define OML_MSG_ARRAYCATDIM   "Error: incompatible array sizes; non-concatenated dimensions must match"
#define OML_MSG_CELLSIZE      "Error: incompatible cell sizes; must match"
#define OML_MSG_OPTION        "Error: invalid option argument"
#define OML_MSG_OPTIONVAL     "Error: invalid option; is incorrectly specified"
#define OML_MSG_FUNCSWITCH    "Error: invalid option; must be either 'on' or 'off'"
#define OML_MSG_NOBUILTIN     "Error: built in function not supported in this context"
#define OML_MSG_CONSTRARG2    "Error: invalid constraint function; argument 2 must be []"
#define OML_MSG_CONSTRRET2    "Error: invalid constraint function; can have at most 2 returns"
#define OML_MSG_CONSTRRET4    "Error: invalid constraint function; can have at most 4 returns"
#define OML_MSG_ANALYTICGRADS "Error: invalid options; GradObj and GradConstr be set the same"
#define OML_MSG_UNSUPPORTDIM  "Error: unsupported dimension; matrices with more than 2 dimensions are not currently allowed"
#define OML_MSG_FLAG_01       "Error: invalid input; must be 0 or 1"
#define OML_MSG_FORMAT        "Error: invalid input; format(s) cannot be applied"
#define OML_MSG_PNORM         "Error: invalid input; p must be positive"
#define OML_MSG_NORM_STRING3  "Error: invalid string; options are 'rows', 'cols', or 'columns'"
#define OML_MSG_NOCOMPLEX     "Error: invalid input; cannot be a complex number"
#define OML_MSG_NATURALNUM_MATRIX_CELLARRAY "Error: invalid input; must be a nonnegative integer, matrix or cell array"
#define OML_MSG_STRING_MATRIX_CELLARRAY     "Error: invalid input; must be a string, matrix or cell array"
#define OML_MSG_STRING_ONEDIMENSION         "Error: invalid input; string input must be one-dimensional"
#define OML_MSG_STRSCALARCOMPLEXMTX         "Error: invalid input; must be a string, scalar, complex or matrix"
#define OML_MSG_FINITE_NATURALNUM           "Error: invalid input; must be a finite, nonnegative integer"
#define OML_MSG_STRING_NATURALNUM           "Error: invalid input; must be a string or a nonnegative integer"
#define OML_MSG_POSITIVE_SCALAR             "Error: invalid input; must be a finite, positive scalar"
#define OML_MSG_SCALAR_COMPLEX              "Error: invalid input; must be a scalar or a complex number"
#define OML_MSG_STRING_STRINGCELL           "Error: invalid input; must be a string or a cell array of strings"
#define OML_MSG_INVALID_INDEX               "Error: invalid index; must be a positive integer"
#define OML_MSG_INVALID_RANGE               "Error: invalid input; must be in valid range"
#define OML_MSG_INVALID_BASE                "Error: invalid input; base must be an integer >= 2"
#define OML_MSG_INVALID_DLL                 "Error: invalid dynamic library"
#define OML_MSG_INVALID_INITDLL             "Error: invalid initialization function in dynamic library"
#define OML_MSG_GUI_CMDEXEC_FAIL            "Error: command execution failed in GUI";
#define OML_MSG_POS_INTEGER_VEC_MTX         "Error: invalid index; must be a positive integer, vector or matrix of positive integers";
#define OML_MSG_STRING_INTEGER              "Error: invalid input; must be a string or an integer";
#define OML_MSG_SCALAR_VECTOR_STRING        "Error: invalid input; must be a scalar, vector or string";
#define OML_MSG_POS_INTEGER_MTX_INF         "Error: invalid input; must be a positive integer, matrix of positive integers or infinity"
#define OML_MSG_POS_INTEGER_MTX_SIZE2       "Error: invalid input; matrix of positive integers must have at least two elements"
#define OML_MSG_CELLMTXSTRUCT               "Error: invalid input; must be a cell array, matrix or struct"
#define OML_MSG_CELLSTRING                  "Error: invalid input; must be a cell array or string"
#define OML_MSG_SCALAROUTOFRANGE            "Error: invalid input; scalar is out of character range for strings"
#define OML_MSG_INVALIDSTRUCTINDEX          "Error: invalid input; struct cannot be indexed"
#define OML_MSG_INVALIDINDEX                "Error: invalid index; must be within allowable bounds of input"
#define OML_MSG_TRIANGMATTYPE               "Error: invalid input; must be 'lower' or 'upper'"
#define OML_MSG_MTXSTRING                   "Error: invalid input; must be a matrix or string"
#define OML_MSG_POSINTEGER_VEC              "Error: invalid input; must be a positive integer or vector"
#define OML_MSG_NONEMPTY_STR                "Error: invalid input; must be a non-empty string"
#define OML_MSG_ONEROW_STRING               "Error: invalid input; must be a string with one row"
#define OML_MSG_SCALAR_REALMTX              "Error: invalid input; must be a scalar or real matrix"
#define OML_MSG_INTEGER_INTMTX              "Error: invalid input; must be an integer or a matrix of integers"
#define OML_MSG_LOGICAL                     "Error: invalid input; must be true or false"
#define OML_MSG_INVALID_VERSION             "Error: invalid version"
#define OML_MSG_STRING_FILESTREAM           "Error: invalid input; must be a string or a valid file stream"
#define OML_MSG_FILE_FILESTREAM             "Error: invalid input; must be a file name or a valid file stream"
#define OML_MSG_FILE_NOTFOUND               "Error: invalid input; cannot find file"
#define OML_MSG_FILE_CANNOTOPEN             "Error: invalid input; cannot open file"
#define OML_MSG_OPT_UNSUPPORTED             "Error: invalid input; option is not supported"
#define OML_MSG_INVALIDTIMERANGE            "Error: invalid input; time index out of range;"
#define OML_MSG_INVALIDTIMESTEP             "Error: invalid input; time step is not valid"
#define OML_MSG_NONNEGATIVE_SCALAR          "Error: invalid input; must be a non-negative scalar"
#define OML_MSG_FILE_CANNOTREAD             "Error: invalid input; cannot read file"
#define OML_MSG_STRINGSTRUCT                "Error: invalid input; must be a string or struct"
#define OML_MSG_SCAL_COMP_MTX_STR_CELL      "Error: invalid input; must be a scalar, complex, matrix, string or cell"
#define OML_MSG_INVALID_DATE_FMT            "Error: invalid format; check help for valid date formats"
#define OML_MSG_CANNOTAPPLY_DATE_FMT        "Error: cannot apply format; check help for applicable date formats"
#define OML_MSG_INPUT_EMTPY                 "Error: invalid input; cannot be empty"
#define OML_MSG_MULTILINE_STRING            "Error: invalid input; cannot be a multiline string"
#define OML_MSG_INVALID_FILE_MODE           "Error: invalid file mode; check help for valid file modes"
#define OML_MSG_FILE_CANNOTWRITE            "Error: invalid input; cannot write to file"
#define OML_MSG_CELL_MTX                    "Error: invalid input; must be a cell or matrix"
#define OML_MSG_INVALIDTIMEVAL              "Error: invalid input: value is out of range for time function"
#define OML_MSG_HANDLE_STRING_CELL          "Error: invalid input: must be a function handle, string or cell with function details"
#define OML_MSG_INTERNAL                    "Error: internal error"
#define OML_MSG_AUTHENTICATE                "Error: authentication failure"

// plot messages
#define OML_MSG_PLOT_DIM_NOT_MATCH              "Error: invalid inputs; data dimensions do not match"
#define OML_MSG_PLOT_MISSING_VALUE				"Error: invalid input; missing value for property"
#define OML_MSG_PLOT_UNMATCHED_AXES_TYPE		"Error: axes types are mismatched; turn hold off"
#define OML_MSG_PLOT_INVALID_PROPERTY			"Error: invalid input; cannot find property"
#define OML_MSG_PLOT_INVALID_OBJECT_HANDLE		"Error: invalid input; cannot find object handle"
#define OML_MSG_PLOT_INVALID_FIGURE_HANDLE		"Error: invalid input; cannot find figure handle"
#define OML_MSG_PLOT_INVALID_AXES_HANDLE		"Error: invalid input; cannot find axes handle"
#define OML_MSG_PLOT_INVALID_CURVE_HANDLE		"Error: invalid input; cannot find curve handle"
#define OML_MSG_PLOT_INVALID_COLOR       		"Error: invalid color option"
#define OML_MSG_PLOT_NOT_SUPPORTED				"Error: command not supported"
#define OML_MSG_PLOT_NOT_SUPPORTED_FOR_2D       "Error: command supported only for 3D plots"
#define OML_MSG_PLOT_UNKNOWN_ERROR				"Error: internal error; operation cannot be completed"
#define OML_MSG_PLOT_ZERORANGE                  "Error: invalid data; has a range of zero"
#define OML_MSG_PLOT_ZIN2D                      "Error: invalid input; z axis is not applicable for 2D plots"
#define OML_MSG_PLOT_UNSUPPORTED_FORMAT         "Error: unsupported file format"
#define OML_MSG_PLOT_EMPTY_PROPERTY             "Error: invalid operation; property name cannot be empty"
#define OML_MSG_PLOT_MATXY_NOT_MATCH            "Error: size of x and y must match"
#define OML_MSG_PLOT_MATXZ_NOT_MATCH            "Error: size of x and z must match"
#define OML_MSG_PLOT_MATYZ_NOT_MATCH            "Error: size of y and z must match"
#define OML_MSG_PLOT_X_Z2_NOT_MATCH             "Error: length of x must match the number of columns of z"
#define OML_MSG_PLOT_Y_Z1_NOT_MATCH             "Error: length of y must match the number of rows of z"
#define OML_MSG_PLOT_XZ_NOT_MATCH               "Error: length of x and z must match"
#define OML_MSG_PLOT_YZ_NOT_MATCH               "Error: length of y and z must match"
#define OML_MSG_PLOT_CONST_PROPERTY             "Error: property is read only; cannot be updated"
#define OML_MSG_PLOT_CANNOT_OPEN_IMAGE          "Error: invalid path; cannot load image"
#define OML_MSG_PLOT_NEED_NORM_DATA             "Error: invalid operation; range of normalized value is [0 1]"
#define OML_MSG_PLOT_NEED_PIXEL_DATA            "Error: invalid operation; pixel value should larger than 1"

// \todo: Used only in open matrix
#define OML_MSG_PLOT_AMBIGUOUS_PROPERTY "Error: ambiguous property"

// ABF toolbox messages
#define OML_MSG_ABF_CREATE_FAILED      "Error: failed to create ABF writer; check input"
#define OML_MSG_ABF_WRITE_FAILED       "Error: file not opened for writing; open file"
#define OML_MSG_ABF_SUBCASE_EMPTY      "Error: subcase name is missing; check input"
#define OML_MSG_ABF_EXPORT_DONE        "Error: data is already exported"
#define OML_MSG_ABF_WRITE_IN_PROGRESS  "Error: data write is in progress"

// HW reader messages
#define OML_MSG_HWREADER_TIMECHANNELS_COMPARE    "Time channels does not match"

// Variable type definitions
#define OML_STR_MATRIX          "matrix"
#define OML_STR_VECTOR          "vector"
#define OML_STR_STRUCT          "struct"
#define OML_STR_CELL            "cell"
#define OML_STR_STRING          "string"
#define OML_STR_INDEX           "index"
#define OML_STR_ORDER           "order"
#define OML_STR_DIM             "dimension"
#define OML_STR_DIMS            "dimensions"
#define OML_STR_TYPE            "type"
#define OML_STR_VALUE           "value"
#define OML_STR_VARIABLE        "variable"
#define OML_STR_DATA            "data"
#define OML_STR_FUNC            "function"
#define OML_STR_INPUT           "input"
#define OML_STR_PARAMETER       "parameter"
#define OML_STR_CONTEXT         "context"
#define OML_STR_TEMPLATE        "template"
#define OML_STR_JACOBIAN        "Jacobian"
#define OML_STR_GRADOBJ         "GradObj"
#define OML_STR_GRADCONSTR      "GradConstr"
#define OML_STR_ABSTOL          "AbsTol"
#define OML_STR_RELTOL          "RelTol"
#define OML_STR_MAXSTEP         "MaxStep"
#define OML_STR_TOLX            "TolX"
#define OML_STR_TOLFUN          "TolFun"
#define OML_STR_TOLFUNABS       "TolFunAbs"
#define OML_STR_TOLFUNREL       "TolFunRel"
#define OML_STR_TOLCON          "TolCon"
#define OML_STR_CONRET          "Constrant Retention"
#define OML_STR_MOVE            "Move Limit Fraction"
#define OML_STR_MPERT           "Perturbation Method"
#define OML_STR_PERT            "Initial Perturbation Value"
#define OML_STR_TOLKKT          "TolKKT"
#define OML_STR_MAXFUNEVALS     "MaxFunEvals"
#define OML_STR_MAXITER         "MaxIter"
#define OML_STR_DISPLAY         "Display"
#define OML_STR_SKIPVAL         "skip value"
#define OML_STR_ORIGIN          "orgin"
#define OML_STR_FILEID          "file ID"
#define OML_STR_OFFSET          "offset"
#define OML_STR_LENGTH          "length"
#define OML_STR_OPTION          "option"

//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
OML_Error::OML_Error(omlMathErrCode errCode, int arg1, int arg2)
    : m_formatMsg (true)
	, m_errCode   (errCode)
	, m_arg1      (arg1)
	, m_arg2      (arg2)
	, m_varCode   (OML_VAR_NONE)
{
}
//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
OML_Error::OML_Error(omlMathErrCode errCode, int arg1, omlMathVarCode varCode)
    : m_formatMsg (true)
	, m_errCode   (errCode)
	, m_arg1      (arg1)
	, m_arg2      (-1)
	, m_varCode   (varCode)
{
}
//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
OML_Error::OML_Error(omlMathErrCode errCode, 
	                 int            arg1, 
	                 int            arg2, 
	                 omlMathVarCode varCode)
    : m_formatMsg (true)
	, m_errCode   (errCode)
	, m_arg1      (arg1)
	, m_arg2      (arg2)
	, m_varCode   (varCode)
{
}
//-----------------------------------------------------------------------------
// Constructor - non standard
//-----------------------------------------------------------------------------
OML_Error::OML_Error(const std::string& message)
    : m_formatMsg (true)
	, m_errCode   (OML_ERR_NONE)
	, m_arg1      (-1)
	, m_arg2      (-1)
	, m_varCode   (OML_VAR_NONE)
	, m_message   (message)
{
    // Please use this constructor sparingly. It should be primarliy for special cases.
    // The use of other constructors ensures better standardization.
}
//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
OML_Error::OML_Error(const hwMathStatus& status)
    : m_status    (status)
    , m_formatMsg (true)
	, m_errCode   (OML_ERR_NONE)
	, m_arg1      (-1)
	, m_arg2      (-1)
	, m_varCode   (OML_VAR_NONE)
{
}
//------------------------------------------------------------------------------
// Non standard constructor - use sparingly
//------------------------------------------------------------------------------
OML_Error::OML_Error(const std::string& message, bool formatMsg)
    : m_errCode   (OML_ERR_NONE)
    , m_arg1      (-1)
    , m_arg2      (-1)
    , m_varCode   (OML_VAR_NONE)
    , m_message   (message)
    , m_formatMsg (formatMsg) 
{
    // This should be used mainly for special cases. Other CTORs ensure better
    // standardization of errors
}
//------------------------------------------------------------------------------
// Returns error message
//------------------------------------------------------------------------------
std::string OML_Error::GetErrorMessage() const
{
	if (!m_message.empty())
	{
		return m_message;
	}

    std::string msgStr;

    if (m_errCode != OML_ERR_NONE)
    {
        if (m_arg1 != -1 && m_arg2 == -1)       // one argument
        {
            std::string message = GetOmlErrorMessage(m_errCode);
            size_t pos = message.find(';');
            size_t len = message.length();

            std::string argChar = std::to_string(static_cast<long long>(m_arg1));

            msgStr = (pos == std::string::npos) ? message : 
                                                  message.substr(0, pos);
                     
            // get message up to semicolon, insert argument number
            msgStr += " in argument " + argChar;

            // insert variable name if exists, append rest of message
            if (m_varCode != OML_VAR_NONE)
            {
                if (pos != std::string::npos && pos < len - 1)
                {
                    msgStr += "; ";
                }
                else 
                {
                    msgStr += " ";
                }
                msgStr += GetOmlVarStr(m_varCode);
                if (pos != std::string::npos)
                {
                    msgStr += message.substr(pos+1, len-pos);
                }
            }
            else if (pos != std::string::npos && pos < len - 1)
            {
                msgStr += "; " + message.substr(pos+2, len-pos);
            }
        }
        else if (m_arg1 != -1 && m_arg2 != -1)  // two arguments
        {
            std::string message = GetOmlErrorMessage(m_errCode);
            size_t pos = message.find(';');
            size_t len = message.length();

			std::string strarg1 (std::to_string(m_arg1));
			std::string strarg2 (std::to_string(m_arg2));

            // get message up to semicolon, insert argument number
            msgStr = message.substr(0, pos) + " in arguments " + strarg1 + "," + strarg2 + "; ";

            // insert variable name if exists, append rest of message
			if (m_varCode != OML_VAR_NONE)
			{
				msgStr += GetOmlVarStr(m_varCode) + message.substr(pos + 1, len - pos);
			}
			else
			{
				msgStr += message.substr(pos + 2, len - pos);
			}
        }
        else    // no arguments
        {
            std::string message = GetOmlErrorMessage(m_errCode);
            size_t pos = message.find(';');

            if (pos != std::string::npos)
            {
                // get message up to semicolon
                msgStr = message.substr(0, pos) + "; ";

                // insert variable name if exists, append rest of message
                size_t len = message.length();

				if (m_varCode != OML_VAR_NONE)
				{
					msgStr += GetOmlVarStr(m_varCode) + message.substr(pos + 1, len - pos);
				}
				else
				{
					msgStr += message.substr(pos + 2, len - pos);
				}
            }
			else
			{
				msgStr = message;
			}
        }
    }
    else if (!m_status.IsOk())
    {
		if (m_status.GetArg1() == 0)
		{
			m_status.ResetArgs();
		}
        msgStr = m_status.GetMessage();
    }

    return msgStr;
}
//-----------------------------------------------------------------------------
// Returns error message for given error code
//-----------------------------------------------------------------------------
std::string OML_Error::GetOmlErrorMessage(omlMathErrCode errCode) const
{
    std::string msgStr;

    switch (errCode)
    {
    case OML_ERR_NUMARGIN:                      msgStr = OML_MSG_NUMARGIN;                      break;
    case OML_ERR_NUMARGOUT:                     msgStr = OML_MSG_NUMARGOUT;                     break;
    case OML_ERR_NUMARGINOUT:                   msgStr = OML_MSG_NUMARGINOUT;                   break;
    case OML_ERR_CELL:                          msgStr = OML_MSG_CELL;                          break;
    case OML_ERR_CELLARRAY:                     msgStr = OML_MSG_CELLARRAY;                     break;
    case OML_ERR_STRUCT:                        msgStr = OML_MSG_STRUCT;                        break;
    case OML_ERR_STRING:                        msgStr = OML_MSG_STRING;                        break;
    case OML_ERR_INTSTRING:                     msgStr = OML_MSG_INTSTRING;                     break;
    case OML_ERR_SCALARSTRING:                  msgStr = OML_MSG_SCALARSTRING;                  break;
    case OML_ERR_BAD_STRING:                    msgStr = OML_MSG_BAD_STRING;                    break;
    case OML_ERR_NUMERIC:                       msgStr = OML_MSG_NUMERIC;                       break;
    case OML_ERR_SCALAR:                        msgStr = OML_MSG_SCALAR;                        break;
    case OML_ERR_VECTOR:                        msgStr = OML_MSG_VECTOR;                        break;
    case OML_ERR_VECTOR2:                       msgStr = OML_MSG_VECTOR2;                       break;
    case OML_ERR_SCALARVECTOR:                  msgStr = OML_MSG_SCALARVECTOR;                  break;
    case OML_ERR_SCALARMATRIX:                  msgStr = OML_MSG_SCALARMATRIX;                  break;
    case OML_ERR_SCALARCOMPLEXMTX:              msgStr = OML_MSG_SCALARCOMPLEXMTX;              break;
    case OML_ERR_STRINGVECTOR:                  msgStr = OML_MSG_STRINGVECTOR;                  break;
    case OML_ERR_INTVECSTR:                     msgStr = OML_MSG_INTVECSTR;                     break;
    case OML_ERR_REALVECTOR:                    msgStr = OML_MSG_REALVECTOR;                    break;
    case OML_ERR_NNINTVECTOR:                   msgStr = OML_MSG_NNINTVECTOR;                   break;
    case OML_ERR_POSINTVECTOR:                  msgStr = OML_MSG_POSINTVECTOR;                  break;
    case OML_ERR_MATRIX:                        msgStr = OML_MSG_MATRIX;                        break;
    case OML_ERR_REALMATRIX:                    msgStr = OML_MSG_REALMATRIX;                    break;
    case OML_ERR_EMPTYMATRIX:                   msgStr = OML_MSG_EMPTYMATRIX;                   break;
    case OML_ERR_HANDLE:                        msgStr = OML_MSG_HANDLE;                        break;
    case OML_ERR_HANDLE_EMPTY:                  msgStr = OML_MSG_HANDLE_EMPTY;                  break;
    case OML_ERR_FUNCNAME:                      msgStr = OML_MSG_FUNCNAME;                      break;
    case OML_ERR_ACCUMFUNC:                     msgStr = OML_MSG_ACCUMFUNC;                     break;
    case OML_ERR_REAL:                          msgStr = OML_MSG_REAL;                          break;
    case OML_ERR_INTEGER:                       msgStr = OML_MSG_INTEGER;                       break;
    case OML_ERR_NATURALNUM:                    msgStr = OML_MSG_NATURALNUM;                    break;
    case OML_ERR_POSINTEGER:                    msgStr = OML_MSG_POSINTEGER;                    break;
    case OML_ERR_FINITE:                        msgStr = OML_MSG_FINITE;                        break;
    case OML_ERR_VECLENDIM:                     msgStr = OML_MSG_VECLENDIM;                     break;
    case OML_ERR_ARRAYSIZE:                     msgStr = OML_MSG_ARRAYSIZE;                     break;
    case OML_ERR_ARRAYCATDIM:                   msgStr = OML_MSG_ARRAYCATDIM;                   break;
    case OML_ERR_CELLSIZE:                      msgStr = OML_MSG_CELLSIZE;                      break;
    case OML_ERR_OPTION:                        msgStr = OML_MSG_OPTION;                        break;
    case OML_ERR_OPTIONVAL:                     msgStr = OML_MSG_OPTIONVAL;                     break;
    case OML_ERR_FUNCSWITCH:                    msgStr = OML_MSG_FUNCSWITCH;                    break;
    case OML_ERR_NOBUILTIN:                     msgStr = OML_MSG_NOBUILTIN;                     break;
    case OML_ERR_CONSTRARG2:                    msgStr = OML_MSG_CONSTRARG2;                    break;
    case OML_ERR_CONSTRRET2:                    msgStr = OML_MSG_CONSTRRET2;                    break;
    case OML_ERR_CONSTRRET4:                    msgStr = OML_MSG_CONSTRRET4;                    break;
    case OML_ERR_ANALYTICGRADS:                 msgStr = OML_MSG_ANALYTICGRADS;                 break;
    case OML_ERR_UNSUPPORTDIM:                  msgStr = OML_MSG_UNSUPPORTDIM;                  break;
    case OML_ERR_FLAG_01:                       msgStr = OML_MSG_FLAG_01;                       break;
    case OML_ERR_FORMAT:                        msgStr = OML_MSG_FORMAT;                        break;
    case OML_ERR_PNORM:                         msgStr = OML_MSG_PNORM;                         break;
    case OML_ERR_NORM_STRING3:                  msgStr = OML_MSG_NORM_STRING3;                  break;
    case OML_ERR_NOCOMPLEX:                     msgStr = OML_MSG_NOCOMPLEX;                     break;
    case OML_ERR_NATURALNUM_MATRIX_CELLARRAY:   msgStr = OML_MSG_NATURALNUM_MATRIX_CELLARRAY;   break;
    case OML_ERR_STRING_MATRIX_CELLARRAY:       msgStr = OML_MSG_STRING_MATRIX_CELLARRAY;       break;    
    case OML_ERR_STRING_ONEDIMENSION:           msgStr = OML_MSG_STRING_ONEDIMENSION;           break;
    case OML_ERR_STRSCALARCOMPLEXMTX:           msgStr = OML_MSG_STRSCALARCOMPLEXMTX;           break;
    case OML_ERR_FINITE_NATURALNUM:             msgStr = OML_MSG_FINITE_NATURALNUM;             break;
    case OML_ERR_STRING_NATURALNUM:             msgStr = OML_MSG_STRING_NATURALNUM;             break;
    case OML_ERR_POSITIVE_SCALAR:               msgStr = OML_MSG_POSITIVE_SCALAR;               break;
    case OML_ERR_SCALAR_COMPLEX:                msgStr = OML_MSG_SCALAR_COMPLEX;                break;
    case OML_ERR_STRING_STRINGCELL:             msgStr = OML_MSG_STRING_STRINGCELL;             break;
    case OML_ERR_INVALID_INDEX:                 msgStr = OML_MSG_INVALID_INDEX;                 break;
    case OML_ERR_INVALID_RANGE:                 msgStr = OML_MSG_INVALID_RANGE;                 break;
    case OML_ERR_INVALID_BASE:                  msgStr = OML_MSG_INVALID_BASE;                  break;
    case OML_ERR_INVALID_DLL:                   msgStr = OML_MSG_INVALID_DLL;                   break;
    case OML_ERR_INVALID_INITDLL:               msgStr = OML_MSG_INVALID_INITDLL;               break;
    case OML_ERR_GUI_CMDEXEC_FAIL:              msgStr = OML_MSG_GUI_CMDEXEC_FAIL;              break;  
    case OML_ERR_POS_INTEGER_VEC_MTX:           msgStr = OML_MSG_POS_INTEGER_VEC_MTX;           break;
    case OML_ERR_STRING_INTEGER:                msgStr = OML_MSG_STRING_INTEGER;                break;
    case OML_ERR_SCALAR_VECTOR_STRING:          msgStr = OML_MSG_SCALAR_VECTOR_STRING;          break;
    case OML_ERR_POS_INTEGER_MTX_INF:           msgStr = OML_MSG_POS_INTEGER_MTX_INF;           break;
    case OML_ERR_POS_INTEGER_MTX_SIZE2:         msgStr = OML_MSG_POS_INTEGER_MTX_SIZE2;         break;
    case OML_ERR_CELLMTXSTRUCT:                 msgStr = OML_MSG_CELLMTXSTRUCT;                 break;
    case OML_ERR_CELLSTRING:                    msgStr = OML_MSG_CELLSTRING;                    break;
    case OML_ERR_SCALAROUTOFRANGE:              msgStr = OML_MSG_SCALAROUTOFRANGE;              break;
    case OML_ERR_INVALIDSTRUCTINDEX:            msgStr = OML_MSG_INVALIDSTRUCTINDEX;            break;
    case OML_ERR_INVALIDINDEX:                  msgStr = OML_MSG_INVALIDINDEX;                  break;
    case OML_ERR_TRIANGMATTYPE:                 msgStr = OML_MSG_TRIANGMATTYPE;                 break;
    case OML_ERR_MTXSTRING:                     msgStr = OML_MSG_MTXSTRING;                     break;
    case OML_ERR_POSINTEGER_VEC:                msgStr = OML_MSG_POSINTEGER_VEC;                break;
    case OML_ERR_NONEMPTY_STR:                  msgStr = OML_MSG_NONEMPTY_STR;                 break;           
    case OML_ERR_ONEROW_STRING:                 msgStr = OML_MSG_ONEROW_STRING;                 break;
    case OML_ERR_SCALAR_REALMTX:                msgStr = OML_MSG_SCALAR_REALMTX; break;
    case OML_ERR_INTEGER_INTMTX:                msgStr = OML_MSG_INTEGER_INTMTX; break;
    case OML_ERR_LOGICAL:                       msgStr = OML_MSG_LOGICAL; break;
    case OML_ERR_INVALID_VERSION:               msgStr = OML_MSG_INVALID_VERSION; break;
    case OML_ERR_STRING_FILESTREAM:             msgStr = OML_MSG_STRING_FILESTREAM;             break;
    case OML_ERR_FILE_FILESTREAM:               msgStr = OML_MSG_FILE_FILESTREAM; break;
	case OML_ERR_FILE_NOTFOUND:                 msgStr = OML_MSG_FILE_NOTFOUND; break;
	case OML_ERR_FILE_CANNOTOPEN:               msgStr = OML_MSG_FILE_CANNOTOPEN; break;
	case OML_ERR_OPT_UNSUPPORTED:               msgStr = OML_MSG_OPT_UNSUPPORTED; break;
	case OML_ERR_INVALIDTIMERANGE:              msgStr = OML_MSG_INVALIDTIMERANGE; break;
	case OML_ERR_INVALIDTIMESTEP:               msgStr = OML_MSG_INVALIDTIMESTEP; break;
    case OML_ERR_NONNEGATIVE_SCALAR:            msgStr = OML_MSG_NONNEGATIVE_SCALAR; break;
    case OML_ERR_FILE_CANNOTREAD:               msgStr = OML_MSG_FILE_CANNOTREAD; break;
    case OML_ERR_STRINGSTRUCT:                  msgStr = OML_MSG_STRINGSTRUCT; break;
    case OML_ERR_SCAL_COMP_MTX_STR_CELL:        msgStr = OML_MSG_SCAL_COMP_MTX_STR_CELL; break;
    case OML_ERR_INVALID_DATE_FMT:              msgStr = OML_MSG_INVALID_DATE_FMT; break;
    case OML_ERR_CANNOTAPPLY_DATE_FMT:          msgStr = OML_MSG_CANNOTAPPLY_DATE_FMT; break;
	case OML_ERR_INPUT_EMPTY:                   msgStr = OML_MSG_INPUT_EMTPY; break;
    case OML_ERR_MULTILINE_STRING:              msgStr = OML_MSG_MULTILINE_STRING; break;
    case OML_ERR_INVALID_FILE_MODE:             msgStr = OML_MSG_INVALID_FILE_MODE; break;
    case OML_ERR_FILE_CANNOTWRITE:              msgStr = OML_MSG_FILE_CANNOTWRITE; break;
    case OML_ERR_CELL_MTX:                      msgStr = OML_MSG_CELL_MTX; break;
    case OML_ERR_INVALIDTIMEVAL:                msgStr = OML_MSG_INVALIDTIMEVAL; break;
    case OML_ERR_HANDLE_STRING_CELL:            msgStr = OML_MSG_HANDLE_STRING_CELL; break;
    case OML_ERR_INTERNAL:                      msgStr = OML_MSG_INTERNAL; break;
    case OML_ERR_AUTHENTICATE:                  msgStr = OML_MSG_AUTHENTICATE; break;

	// plot error messages:
    case OML_ERR_PLOT_DIM_NOT_MATCH:            msgStr = OML_MSG_PLOT_DIM_NOT_MATCH;            break;
    case OML_ERR_PLOT_MISSING_VALUE:            msgStr = OML_MSG_PLOT_MISSING_VALUE;            break;
    case OML_ERR_PLOT_UNMATCHED_AXES_TYPE:      msgStr = OML_MSG_PLOT_UNMATCHED_AXES_TYPE;      break;
    case OML_ERR_PLOT_INVALID_PROPERTY:         msgStr = OML_MSG_PLOT_INVALID_PROPERTY;         break;
    case OML_ERR_PLOT_INVALID_OBJECT_HANDLE:    msgStr = OML_MSG_PLOT_INVALID_OBJECT_HANDLE;    break;
    case OML_ERR_PLOT_INVALID_FIGURE_HANDLE:    msgStr = OML_MSG_PLOT_INVALID_FIGURE_HANDLE;    break;
    case OML_ERR_PLOT_INVALID_AXES_HANDLE:      msgStr = OML_MSG_PLOT_INVALID_AXES_HANDLE;      break;
    case OML_ERR_PLOT_INVALID_CURVE_HANDLE:     msgStr = OML_MSG_PLOT_INVALID_CURVE_HANDLE;     break;
    case OML_ERR_PLOT_INVALID_COLOR:            msgStr = OML_MSG_PLOT_INVALID_COLOR;            break;
    case OML_ERR_PLOT_NOT_SUPPORTED:            msgStr = OML_MSG_PLOT_NOT_SUPPORTED;            break;
    case OML_ERR_PLOT_NOT_SUPPORTED_FOR_2D:     msgStr = OML_MSG_PLOT_NOT_SUPPORTED_FOR_2D;     break;
    case OML_ERR_PLOT_UNKNOWN_ERROR:            msgStr = OML_MSG_PLOT_UNKNOWN_ERROR;            break;
    case OML_ERR_PLOT_ZERORANGE:                msgStr = OML_MSG_PLOT_ZERORANGE;                break;
    case OML_ERR_PLOT_ZIN2D:                    msgStr = OML_MSG_PLOT_ZIN2D;                    break;
    case OML_ERR_PLOT_UNSUPPORTED_FORMAT:       msgStr = OML_MSG_PLOT_UNSUPPORTED_FORMAT;       break;
    case OML_ERR_PLOT_EMPTY_PROPERTY:           msgStr = OML_MSG_PLOT_EMPTY_PROPERTY;           break;
    case OML_ERR_PLOT_MATXY_NOT_MATCH:          msgStr = OML_MSG_PLOT_MATXY_NOT_MATCH;          break;
    case OML_ERR_PLOT_MATXZ_NOT_MATCH:          msgStr = OML_MSG_PLOT_MATXZ_NOT_MATCH;          break;
    case OML_ERR_PLOT_MATYZ_NOT_MATCH:          msgStr = OML_MSG_PLOT_MATYZ_NOT_MATCH;          break;
    case OML_ERR_PLOT_XZ_NOT_MATCH:             msgStr = OML_MSG_PLOT_XZ_NOT_MATCH;             break;
    case OML_ERR_PLOT_YZ_NOT_MATCH:             msgStr = OML_MSG_PLOT_YZ_NOT_MATCH;             break;
    case OML_ERR_PLOT_X_Z2_NOT_MATCH:           msgStr = OML_MSG_PLOT_X_Z2_NOT_MATCH;           break;
    case OML_ERR_PLOT_Y_Z1_NOT_MATCH:           msgStr = OML_MSG_PLOT_Y_Z1_NOT_MATCH;           break;
    case OML_ERR_PLOT_CONST_PROPERTY:           msgStr = OML_MSG_PLOT_CONST_PROPERTY;           break;
    case OML_ERR_PLOT_CANNOT_OPEN_IMAGE:        msgStr = OML_MSG_PLOT_CANNOT_OPEN_IMAGE;        break;
    case OML_ERR_PLOT_NEED_NORM_DATA:           msgStr = OML_MSG_PLOT_NEED_NORM_DATA;           break;
    case OML_ERR_PLOT_NEED_PIXEL_DATA:          msgStr = OML_MSG_PLOT_NEED_PIXEL_DATA;          break;
    case OML_ERR_PLOT_AMBIGUOUS_PROPERTY:       msgStr = OML_MSG_PLOT_AMBIGUOUS_PROPERTY;       break;

	// ABF error messages:
	case OML_ERR_ABF_CREATE_FAILED:             msgStr = OML_MSG_ABF_CREATE_FAILED;             break;
	case OML_ERR_ABF_WRITE_FAILED:              msgStr = OML_MSG_ABF_WRITE_FAILED;              break;
	case OML_ERR_ABF_SUBCASE_EMPTY:             msgStr = OML_MSG_ABF_SUBCASE_EMPTY;             break;
	case OML_ERR_ABF_EXPORT_DONE:               msgStr = OML_MSG_ABF_EXPORT_DONE;               break;

    // HW reader error messages:
    case OML_ERR_HWREADER_TIMECHANNELS_COMPARE: msgStr = OML_MSG_HWREADER_TIMECHANNELS_COMPARE; break;
    default: break;
    }

    return msgStr;
}
//-----------------------------------------------------------------------------
// Returns variable type string for given code
//-----------------------------------------------------------------------------
std::string OML_Error::GetOmlVarStr(omlMathVarCode varCode) const
{
    std::string varStr;

    switch (varCode)
    {
    case OML_VAR_MATRIX:       varStr = OML_STR_MATRIX;       break;
    case OML_VAR_VECTOR:       varStr = OML_STR_VECTOR;       break;
    case OML_VAR_STRUCT:       varStr = OML_STR_STRUCT;       break;
    case OML_VAR_CELL:         varStr = OML_STR_CELL;         break;
    case OML_VAR_STRING:       varStr = OML_STR_STRING;       break;
    case OML_VAR_INDEX:        varStr = OML_STR_INDEX;        break;
    case OML_VAR_ORDER:        varStr = OML_STR_ORDER;        break;
    case OML_VAR_DIM:          varStr = OML_STR_DIM;          break;
    case OML_VAR_DIMS:         varStr = OML_STR_DIMS;         break;
    case OML_VAR_TYPE:         varStr = OML_STR_TYPE;         break;
    case OML_VAR_VALUE:        varStr = OML_STR_VALUE;        break;
    case OML_VAR_VARIABLE:     varStr = OML_STR_VARIABLE;     break;
    case OML_VAR_DATA:         varStr = OML_STR_DATA;         break;
    case OML_VAR_FUNC:         varStr = OML_STR_FUNC;         break;
    case OML_VAR_INPUT:        varStr = OML_STR_INPUT;        break;
    case OML_VAR_PARAMETER:    varStr = OML_STR_PARAMETER;    break;
    case OML_VAR_CONTEXT:      varStr = OML_STR_CONTEXT;      break;
    case OML_VAR_TEMPLATE:     varStr = OML_STR_TEMPLATE;     break;
    case OML_VAR_JACOBIAN:     varStr = OML_STR_JACOBIAN;     break;
    case OML_VAR_GRADOBJ:      varStr = OML_STR_GRADOBJ;      break;
    case OML_VAR_GRADCONSTR:   varStr = OML_STR_GRADCONSTR;   break;
    case OML_VAR_ABSTOL:       varStr = OML_STR_ABSTOL;       break;
    case OML_VAR_RELTOL:       varStr = OML_STR_RELTOL;       break;
    case OML_VAR_MAXSTEP:      varStr = OML_STR_MAXSTEP;      break;
    case OML_VAR_TOLX:         varStr = OML_STR_TOLX;         break;
    case OML_VAR_TOLFUN:       varStr = OML_STR_TOLFUN;       break;
    case OML_VAR_TOLFUNABS:    varStr = OML_STR_TOLFUNABS;    break;
    case OML_VAR_TOLFUNREL:    varStr = OML_STR_TOLFUNREL;    break;
    case OML_VAR_TOLCON:       varStr = OML_STR_TOLCON;       break;
    case OML_VAR_CONRET:       varStr = OML_STR_CONRET;       break;
    case OML_VAR_MOVE:         varStr = OML_STR_MOVE;         break;
    case OML_VAR_MPERT:        varStr = OML_STR_MPERT;        break;
    case OML_VAR_PERT:         varStr = OML_STR_PERT;         break;
    case OML_VAR_TOLKKT:       varStr = OML_STR_TOLKKT;       break;
    case OML_VAR_MAXFUNEVALS:  varStr = OML_STR_MAXFUNEVALS;  break;
    case OML_VAR_MAXITER:      varStr = OML_STR_MAXITER;      break;
    case OML_VAR_DISPLAY:      varStr = OML_STR_DISPLAY;      break;
    case OML_VAR_SKIPVAL:      varStr = OML_STR_SKIPVAL;      break;
    case OML_VAR_ORIGIN:       varStr = OML_STR_ORIGIN;       break;
    case OML_VAR_FILEID:       varStr = OML_STR_FILEID;       break;
    case OML_VAR_OFFSET:       varStr = OML_STR_OFFSET;       break;
    case OML_VAR_LENGTH:       varStr = OML_STR_LENGTH;       break;
    case OML_VAR_OPTION:       varStr = OML_STR_OPTION;       break;
    default: break;
    }

    return varStr;
}
