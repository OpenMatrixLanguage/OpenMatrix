/**
* @file OML_Error.cpp
* @date November 2015
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

#include "OML_Error.h"
#include "EvaluatorInt.h"

#if defined(_DARWIN) || defined(LINUX)
  #include <stdio.h>
  #define _itoa_s(value, str, size, base) ::sprintf((str), "%d", value)
#endif

OML_Error::OML_Error(omlMathErrCode errCode, int arg1, int arg2)
    : m_formatMsg (true)
{
    m_errCode = errCode;
    m_arg1 = arg1;
    m_arg2 = arg2;
    m_varCode = OML_VAR_NONE;
}

OML_Error::OML_Error(omlMathErrCode errCode, int arg1, omlMathVarCode varCode)
    : m_formatMsg (true)
{
    m_errCode = errCode;
    m_arg1 = arg1;
    m_arg2 = -1;
    m_varCode = varCode;
}

OML_Error::OML_Error(omlMathErrCode errCode, int arg1, int arg2, omlMathVarCode varCode)
    : m_formatMsg (true)
{
    m_errCode = errCode;
    m_arg1 = arg1;
    m_arg2 = arg2;
    m_varCode = varCode;
}

OML_Error::OML_Error(const std::string& message)
    : m_formatMsg (true)
{
    // Please use this constructor sparingly. It should be primarliy for special cases.
    // The use of other constructors ensures better standardization.
    m_errCode = OML_ERR_NONE;
    m_arg1 = -1;
    m_arg2 = -1;
    m_varCode = OML_VAR_NONE;
    m_message = message;
}

OML_Error::OML_Error(const hwMathStatus& status)
    : m_status(status)
    , m_formatMsg (true)
{
    m_errCode = OML_ERR_NONE;
    m_arg1 = -1;
    m_arg2 = -1;
    m_varCode = OML_VAR_NONE;
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
std::string OML_Error::GetErrorMessage() const
{
    if (!m_message.empty())
        return m_message;

    std::string msgStr;

    if (m_errCode != OML_ERR_NONE)
    {
        if (m_arg1 != -1 && m_arg2 == -1)       // one argument
        {
            std::string message = GetComposeErrMsg(m_errCode);
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
                msgStr += GetComposeVarStr(m_varCode);
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
            std::string message = GetComposeErrMsg(m_errCode);
            size_t pos = message.find(';');
            size_t len = message.length();
            char argChar1[8];
            char argChar2[8];
            _itoa_s(m_arg1, argChar1, 8, 10);
            _itoa_s(m_arg2, argChar2, 8, 10);

            // get message up to semicolon, insert argument number
            msgStr = message.substr(0, pos) + " in arguments " + (std::string(argChar1) + "," + std::string(argChar2) + "; ");

            // insert variable name if exists, append rest of message
            if (m_varCode != OML_VAR_NONE)
                msgStr += GetComposeVarStr(m_varCode) + message.substr(pos+1, len-pos);
            else
                msgStr += message.substr(pos+2, len-pos);
        }
        else    // no arguments
        {
            std::string message = GetComposeErrMsg(m_errCode);
            size_t pos = message.find(';');

            if (pos != std::string::npos)
            {
                // get message up to semicolon
                msgStr = message.substr(0, pos) + "; ";

                // insert variable name if exists, append rest of message
                size_t len = message.length();

                if (m_varCode != OML_VAR_NONE)
                    msgStr += GetComposeVarStr(m_varCode) + message.substr(pos+1, len-pos);
                else
                    msgStr += message.substr(pos+2, len-pos);
            }
            else
                msgStr = message;
        }
    }
    else if (!m_status.IsOk())
    {
        if (m_status.GetArg1() == 0)
            m_status.ResetArgs();

        msgStr = m_status.GetMessage();
    }

    return msgStr;
}

std::string GetComposeErrMsg(omlMathErrCode errCode)
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
    case OML_ERR_STRING_FILESTREAM:             msgStr = OML_MSG_STRING_FILESTREAM;             break;
    case OML_ERR_SCALAR_COMPLEX:                msgStr = OML_MSG_SCALAR_COMPLEX;                break;
    case OML_ERR_STRING_STRINGCELL:             msgStr = OML_MSG_STRING_STRINGCELL;             break;
    case OML_ERR_INVALID_INDEX:                 msgStr = OML_MSG_INVALID_INDEX;                 break;
    case OML_ERR_INVALID_RANGE:                 msgStr = OML_MSG_INVALID_RANGE;                 break;
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
    case OML_ERR_FILE_FILESTREAM:               msgStr = OML_MSG_FILE_FILESTREAM; break;

    // plot error messages:
    case OML_ERR_PLOT_OUT_OF_RANGE:             msgStr = OML_MSG_PLOT_OUT_OF_RANGE;             break;
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
    case OML_ERR_PLOT_NOT_SUPPORTED_OPTION:     msgStr = OML_MSG_PLOT_NOT_SUPPORTED_OPTION;     break;
    case OML_ERR_PLOT_UNKNOWN_ERROR:            msgStr = OML_MSG_PLOT_UNKNOWN_ERROR;            break;
    case OML_ERR_PLOT_ZERORANGE:                msgStr = OML_MSG_PLOT_ZERORANGE;                break;
    case OML_ERR_PLOT_ZIN2D:                    msgStr = OML_MSG_PLOT_ZIN2D;                    break;
    case OML_ERR_PLOT_LIMDATA:                  msgStr = OML_MSG_PLOT_LIMDATA;                  break;
    case OML_ERR_PLOT_UNSUPPORTED_FORMAT:       msgStr = OML_MSG_PLOT_UNSUPPORTED_FORMAT;       break;
    case OML_ERR_PLOT_AMBIGUOUS_PROPERTY:       msgStr = OML_MSG_PLOT_AMBIGUOUS_PROPERTY;       break;
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
    default: break;
    }

    return msgStr;
}

std::string GetComposeVarStr(omlMathVarCode varCode)
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
    case OML_VAR_TOLX:         varStr = OML_STR_TOLX;         break;
    case OML_VAR_TOLFUN:       varStr = OML_STR_TOLFUN;       break;
    case OML_VAR_TOLCON:       varStr = OML_STR_TOLCON;       break;
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
