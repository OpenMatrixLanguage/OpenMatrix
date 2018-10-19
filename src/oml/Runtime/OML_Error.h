/**
* @file OML_Error.h
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

#ifndef _OML_Math_Error_h
#define _OML_Math_Error_h

// Note: each OML_MSG_ITEM is paired with OML_ERR_ITEM in omlMathErrCode below. Likewise,
// each  OML_STR_QUANTITY is paired with OML_VAR_QUANTITY in omlMathVarCode. Here
// are some options for how to throw an error.
// 1. throw OML_Error(OML_ERR_ITEM);
// 2. throw OML_Error(OML_ERR_ITEM, arg_num);
// 3. throw OML_Error(OML_ERR_ITEM, arg_num, OML_VAR_QUANTITY);
// See the OML_Error constructors for all options. When an argument number is to be included
// in the message, OML_MSG_ITEM should be written in the form "Error: problem; solution"
#define OML_MSG_NUMARGIN      "Error: invalid function call; incorrect number of input arguments"
#define OML_MSG_NUMARGOUT     "Error: invalid function call; incorrect number of output arguments"
#define OML_MSG_NUMARGINOUT   "Error: invalid function call; incorrect combination of input/output arguments"
#define OML_MSG_CELL          "Error: invalid input; must be cell"
#define OML_MSG_CELLARRAY     "Error: invalid input; must be cell array"
#define OML_MSG_STRUCT        "Error: invalid input; must be struct"
#define OML_MSG_STRING        "Error: invalid input; must be string"
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
#define OML_MSG_EMPTYMATRIX   "Error: invalid input; must be empty matrix"
#define OML_MSG_HANDLE        "Error: invalid input; must be function handle"
#define OML_MSG_HANDLE_EMPTY  "Error: invalid input; must be function handle or []"
#define OML_MSG_FUNCNAME      "Error: invalid input; function not found"
#define OML_MSG_REAL          "Error: invalid input; must be real"
#define OML_MSG_INTEGER       "Error: invalid input; must be integer"
#define OML_MSG_NATURALNUM    "Error: invalid input; must be nonnegative integer"
#define OML_MSG_POSINTEGER    "Error: invalid input; must be positive integer"
#define OML_MSG_FINITE        "Error: invalid value; must be finite"
#define OML_MSG_VECLENDIM     "Error: invalid vector length in specified dimension"
#define OML_MSG_ARRAYSIZE     "Error: incompatible array sizes; must match"
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
#define OML_MSG_FORMAT        "Error: invalid format specifier"
#define OML_MSG_PNORM         "Error: invalid input; p must be positive"
#define OML_MSG_NOCOMPLEX     "Error: invalid input; cannot be a complex number"
#define OML_MSG_NATURALNUM_MATRIX_CELLARRAY "Error: invalid input; must be a nonnegative integer, matrix or cell array"
#define OML_MSG_STRING_MATRIX_CELLARRAY     "Error: invalid input; must be a string, matrix or cell array"
#define OML_MSG_STRING_ONEDIMENSION         "Error: invalid input; string input must be one-dimensional"
#define OML_MSG_STRSCALARCOMPLEXMTX         "Error: invalid input; must be a string, scalar, complex or matrix"
#define OML_MSG_FINITE_NATURALNUM           "Error: invalid input; must be a finite, nonnegative integer"
#define OML_MSG_STRING_NATURALNUM           "Error: invalid input; must be a string or a nonnegative integer"
#define OML_MSG_POSITIVE_SCALAR             "Error: invalid input; must be a finite, positive scalar"
#define OML_MSG_STRING_FILESTREAM           "Error: invalid input; must be a string or a valid file stream"
#define OML_MSG_SCALAR_COMPLEX              "Error: invalid input; must be a scalar or a complex number"
#define OML_MSG_STRING_STRINGCELL           "Error: invalid input; must be a string or a cell array of strings"
#define OML_MSG_INVALID_INDEX               "Error: invalid index; must be a positive integer"
#define OML_MSG_INVALID_RANGE               "Error: invalid input; must be in valid range"
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
#define OML_MSG_INVALID_TIME_RANGE          "Error: time index out of range; check input"
#define OML_MSG_INVALID_TIME_STEP           "Error: invalid time step; check input"
#define OML_MSG_MTXSTRING                   "Error: invalid input; must be a matrix or string"
#define OML_MSG_POSINTEGER_VEC              "Error: invalid input; must be a positive integer or vector"
#define OML_MSG_NONEMPTY_STR                "Error: invalid input; must be a non-empty string"
#define OML_MSG_ONEROW_STRING               "Error: invalid input; must be a string with one row"
#define OML_MSG_SCALAR_REALMTX              "Error: invalid input; must be a scalar or real matrix"
#define OML_MSG_INTEGER_INTMTX              "Error: invalid input; must be an integer or a matrix of integers"
#define OML_MSG_LOGICAL                     "Error: invalid input; must be true or false"
#define OML_MSG_INVALID_DLL                 "Error: invalid dynamic library"
#define OML_MSG_INVALID_VERSION        "Error: invalid version"

// plot messages
#define OML_MSG_PLOT_OUT_OF_RANGE				"Error: index out of range; check input"
#define OML_MSG_PLOT_DIM_NOT_MATCH              "Error: data dimension do not match; check input"
#define OML_MSG_PLOT_MISSING_VALUE				"Error: invalid input; property must be followd by a value"
#define OML_MSG_PLOT_UNMATCHED_AXES_TYPE		"Error: axes type not matched; turn hold off"
#define OML_MSG_PLOT_INVALID_PROPERTY			"Error: invalid property; check input"
#define OML_MSG_PLOT_INVALID_OBJECT_HANDLE		"Error: invalid object handle; check input"
#define OML_MSG_PLOT_INVALID_FIGURE_HANDLE		"Error: invalid figure handle; check input"
#define OML_MSG_PLOT_INVALID_AXES_HANDLE		"Error: invalid axes handle; check input"
#define OML_MSG_PLOT_INVALID_CURVE_HANDLE		"Error: invalid curve handle; check input"
#define OML_MSG_PLOT_INVALID_COLOR       		"Error: invalid color option; check input"
#define OML_MSG_PLOT_NOT_SUPPORTED				"Error: command not supported"
#define OML_MSG_PLOT_NOT_SUPPORTED_FOR_2D       "Error: command supported for 3d plot only"
#define OML_MSG_PLOT_NOT_SUPPORTED_OPTION		"Error: invalid argument; option not supported"
#define OML_MSG_PLOT_UNKNOWN_ERROR				"Error: unknown error"
#define OML_MSG_PLOT_ZERORANGE                  "Error: invalid data; has a range of zero"
#define OML_MSG_PLOT_ZIN2D                      "Error: no z axis in 2D plot"
#define OML_MSG_PLOT_LIMDATA                    "Error: invalid vector; must have 2 elements"
#define OML_MSG_PLOT_UNSUPPORTED_FORMAT         "Error: unsupported file format"
#define OML_MSG_PLOT_AMBIGUOUS_PROPERTY         "Error: ambiguous property"
#define OML_MSG_PLOT_MATXY_NOT_MATCH            "Error: size of x and y must match"
#define OML_MSG_PLOT_MATXZ_NOT_MATCH            "Error: size of x and z must match"
#define OML_MSG_PLOT_MATYZ_NOT_MATCH            "Error: size of y and z must match"
#define OML_MSG_PLOT_X_Z2_NOT_MATCH             "Error: length of x must match the number of columns of z"
#define OML_MSG_PLOT_Y_Z1_NOT_MATCH             "Error: length of y must match the number of rows of z"
#define OML_MSG_PLOT_XZ_NOT_MATCH               "Error: length of x and z must match"
#define OML_MSG_PLOT_YZ_NOT_MATCH               "Error: length of y and z must match"
#define OML_MSG_PLOT_CONST_PROPERTY             "Error: property is read only; cannot be changed"
#define OML_MSG_PLOT_CANNOT_OPEN_IMAGE          "Error: cannot load image file. Check file path"
#define OML_MSG_PLOT_NEED_NORM_DATA             "Error: the range of normalized value is [0 1]"
#define OML_MSG_PLOT_NEED_PIXEL_DATA            "Error: pixel value should larger than 1"
#define OML_MSG_PLOT_GNUPLOT_NOT_FOUND          "Error: gnuplot not found"

enum omlMathErrCode
{
    OML_ERR_NONE = 0,
    OML_ERR_NUMARGIN,
    OML_ERR_NUMARGOUT,
    OML_ERR_NUMARGINOUT,
    OML_ERR_CELL,
    OML_ERR_CELLARRAY,
    OML_ERR_STRUCT,
    OML_ERR_STRING,
    OML_ERR_NUMERIC,
    OML_ERR_SCALAR,
    OML_ERR_VECTOR,
    OML_ERR_VECTOR2,
    OML_ERR_SCALARVECTOR,
    OML_ERR_SCALARMATRIX,
    OML_ERR_SCALARCOMPLEXMTX,
    OML_ERR_STRINGVECTOR,
    OML_ERR_INTVECSTR,
    OML_ERR_REALVECTOR,
    OML_ERR_NNINTVECTOR,
    OML_ERR_POSINTVECTOR,
    OML_ERR_MATRIX,
    OML_ERR_REALMATRIX,
    OML_ERR_EMPTYMATRIX,
    OML_ERR_HANDLE,
    OML_ERR_HANDLE_EMPTY,
    OML_ERR_FUNCNAME,
    OML_ERR_REAL,
    OML_ERR_INTEGER,
    OML_ERR_NATURALNUM,
    OML_ERR_POSINTEGER,
    OML_ERR_FINITE,
    OML_ERR_VECLENDIM,
    OML_ERR_ARRAYSIZE,
    OML_ERR_ARRAYCATDIM,
    OML_ERR_CELLSIZE,
    OML_ERR_OPTION,
    OML_ERR_OPTIONVAL,
    OML_ERR_FUNCSWITCH,
    OML_ERR_NOBUILTIN,
    OML_ERR_CONSTRARG2,
    OML_ERR_CONSTRRET2,
    OML_ERR_CONSTRRET4,
    OML_ERR_ANALYTICGRADS,
    OML_ERR_UNSUPPORTDIM,
    OML_ERR_FLAG_01,
	OML_ERR_FORMAT,
    OML_ERR_PNORM,
    OML_ERR_NOCOMPLEX,
    OML_ERR_NATURALNUM_MATRIX_CELLARRAY,
    OML_ERR_STRING_MATRIX_CELLARRAY,
    OML_ERR_STRING_ONEDIMENSION,
    OML_ERR_STRSCALARCOMPLEXMTX,
    OML_ERR_FINITE_NATURALNUM,
    OML_ERR_STRING_NATURALNUM,
    OML_ERR_POSITIVE_SCALAR,
    OML_ERR_POS_INTEGER_MTX_INF,
    OML_ERR_POS_INTEGER_MTX_SIZE2,
    OML_ERR_STRING_FILESTREAM,
    OML_ERR_SCALAR_COMPLEX,
    OML_ERR_STRING_STRINGCELL,
    OML_ERR_INVALID_INDEX,
    OML_ERR_INVALID_RANGE,
    OML_ERR_INVALID_DLL,
    OML_ERR_INVALID_INITDLL,
    OML_ERR_GUI_CMDEXEC_FAIL,
    OML_ERR_POS_INTEGER_VEC_MTX,
    OML_ERR_STRING_INTEGER,
    OML_ERR_SCALAR_VECTOR_STRING,
    OML_ERR_CELLMTXSTRUCT,
    OML_ERR_CELLSTRING,
    OML_ERR_SCALAROUTOFRANGE,
    OML_ERR_INVALIDSTRUCTINDEX,
    OML_ERR_INVALIDINDEX,
    OML_ERR_TRIANGMATTYPE,
    OML_ERR_MTXSTRING,
    OML_ERR_POSINTEGER_VEC,
    OML_ERR_NONEMPTY_STR,
    OML_ERR_ONEROW_STRING,
    OML_ERR_SCALAR_REALMTX,
    OML_ERR_INTEGER_INTMTX,
    OML_ERR_LOGICAL,
    OML_ERR_INVALID_VERSION,

    // plot codes
    OML_ERR_PLOT_OUT_OF_RANGE,
    OML_ERR_PLOT_DIM_NOT_MATCH,
    OML_ERR_PLOT_MISSING_VALUE,
    OML_ERR_PLOT_UNMATCHED_AXES_TYPE,
    OML_ERR_PLOT_INVALID_PROPERTY,
    OML_ERR_PLOT_INVALID_OBJECT_HANDLE,
    OML_ERR_PLOT_INVALID_FIGURE_HANDLE,
    OML_ERR_PLOT_INVALID_AXES_HANDLE,
    OML_ERR_PLOT_INVALID_CURVE_HANDLE,
    OML_ERR_PLOT_INVALID_COLOR,
    OML_ERR_PLOT_NOT_SUPPORTED,
    OML_ERR_PLOT_NOT_SUPPORTED_FOR_2D,
    OML_ERR_PLOT_NOT_SUPPORTED_OPTION,
    OML_ERR_PLOT_UNKNOWN_ERROR,
    OML_ERR_PLOT_ZERORANGE,
    OML_ERR_PLOT_ZIN2D,
    OML_ERR_PLOT_LIMDATA,
    OML_ERR_PLOT_UNSUPPORTED_FORMAT,
    OML_ERR_PLOT_AMBIGUOUS_PROPERTY,
    OML_ERR_PLOT_MATXY_NOT_MATCH,
    OML_ERR_PLOT_MATXZ_NOT_MATCH,
    OML_ERR_PLOT_MATYZ_NOT_MATCH,
    OML_ERR_PLOT_XZ_NOT_MATCH,
    OML_ERR_PLOT_YZ_NOT_MATCH,
    OML_ERR_PLOT_X_Z2_NOT_MATCH,
    OML_ERR_PLOT_Y_Z1_NOT_MATCH,
	OML_ERR_PLOT_CONST_PROPERTY,
    OML_ERR_PLOT_CANNOT_OPEN_IMAGE,
    OML_ERR_PLOT_NEED_NORM_DATA,
    OML_ERR_PLOT_NEED_PIXEL_DATA,
    OML_ERR_PLOT_GNUPLOT_NOT_FOUND,

    OML_ERR_END
};

#define OML_STR_MATRIX          "matrix"
#define OML_STR_VECTOR          "vector"
#define OML_STR_STRUCT          "struct"
#define OML_STR_STRING          "string"
#define OML_STR_INDEX           "index"
#define OML_STR_ORDER           "order"
#define OML_STR_DIM             "dimension"
#define OML_STR_DIMS            "dimensions"
#define OML_STR_TYPE            "type"
#define OML_STR_VALUE           "value"
#define OML_STR_VARIABLE        "variable"
#define OML_STR_DATA            "data"
#define OML_STR_INPUT           "input"
#define OML_STR_PARAMETER       "parameter"
#define OML_STR_CONTEXT         "context"
#define OML_STR_TEMPLATE        "template"
#define OML_STR_JACOBIAN        "Jacobian"
#define OML_STR_GRADOBJ         "GradObj"
#define OML_STR_GRADCONSTR      "GradConstr"
#define OML_STR_ABSTOL          "AbsTol"
#define OML_STR_RELTOL          "RelTol"
#define OML_STR_TOLX            "TolX"
#define OML_STR_TOLFUN          "TolFun"
#define OML_STR_TOLCON          "TolCon"
#define OML_STR_TOLKKT          "TolKKT"
#define OML_STR_MAXFUNEVALS     "MaxFunEvals"
#define OML_STR_MAXITER         "MaxIter"
#define OML_STR_DISPLAY         "Display"
#define OML_STR_SKIPVAL         "skip value"
#define OML_STR_ORIGIN          "orgin"
#define OML_STR_FILEID          "file ID"
#define OML_STR_OFFSET          "offset"
#define OML_STR_LENGTH          "length"

enum omlMathVarCode
{
    OML_VAR_NONE = 0,
    OML_VAR_MATRIX,
    OML_VAR_VECTOR,
    OML_VAR_STRUCT,
    OML_VAR_STRING,
    OML_VAR_INDEX,
    OML_VAR_ORDER,
    OML_VAR_DIM,
    OML_VAR_DIMS,
    OML_VAR_TYPE,
    OML_VAR_VALUE,
    OML_VAR_VARIABLE,
    OML_VAR_DATA,
    OML_VAR_INPUT,
    OML_VAR_PARAMETER,
    OML_VAR_CONTEXT,
    OML_VAR_TEMPLATE,
    OML_VAR_JACOBIAN,
    OML_VAR_GRADOBJ,
    OML_VAR_GRADCONSTR,
    OML_VAR_ABSTOL,
    OML_VAR_RELTOL,
    OML_VAR_TOLX,
    OML_VAR_TOLFUN,
    OML_VAR_TOLCON,
    OML_VAR_TOLKKT,
    OML_VAR_MAXFUNEVALS,
    OML_VAR_MAXITER,
    OML_VAR_DISPLAY,
    OML_VAR_SKIPVAL,
    OML_VAR_ORIGIN,
    OML_VAR_FILEID,
    OML_VAR_OFFSET,
    OML_VAR_LENGTH,
    OML_VAR_END
};

#include <string>
#include "Hml2Dll.h"
#include "hwMathStatus.h"

using std::string;

class HML2DLL_DECLS OML_Error
{
public:
	OML_Error(omlMathErrCode errCode, int arg1 = -1, int arg2 = -1);
	OML_Error(omlMathErrCode errCode, int arg1, omlMathVarCode varCode);
	OML_Error(omlMathErrCode errCode, int arg1, int arg2, omlMathVarCode varCode);
	explicit OML_Error(const std::string& message);   // non-standard format, please use sparingly
	explicit OML_Error(const hwMathStatus& status);

    //!
    //! Constructor - non-standard, please use sparingly
    //! \param message   Message to display
    //! \param formatMsg True if line info needs to be appended to message
    //!
	OML_Error(const std::string& message,
              bool               formatMsg);   

	std::string GetErrorMessage() const;
    omlMathErrCode ErrCode() const { return m_errCode; }
    int Arg1() const { return m_arg1; }
    void Arg1(int arg) { m_arg1 = arg; }
    int Arg2() const { return m_arg2; }
    void Arg2(int arg) { m_arg2 = arg; }

    //!
    //! Returns true if error needs to be formatted with line/file info
    //!
    bool GetFormatMessage() const { return m_formatMsg; }
private:
    omlMathErrCode m_errCode;
    int m_arg1;
    int m_arg2;
    omlMathVarCode m_varCode;
    std::string m_message;
    mutable hwMathStatus m_status;
    bool m_formatMsg;          //!< True if format (line info) needs to be added
};

std::string GetComposeErrMsg(omlMathErrCode errCode);
std::string GetComposeVarStr(omlMathVarCode varCode);

// NOTICE: DEPRECATED #define list. PLEASE DO NOT ADD TO IT.

// The #defines in this list should be consolidated with the
// #define OML_MSG_XXXXXX list above, and eventually with those
// in /hwcommon/math/core/Globals.h as well.

#define HW_ERROR_NOTACCEPTINPAFTERSEEDSTATE "Error: Cannot accept inputs after seed/state"
#define HW_ERROR_ALLINPMATCHSPECSIZE "Error: all inputs must match the specified size"
#define HW_ERROR_ALREADYSETDIM "Error: already set dimension"
#define HW_ERROR_ALREADYSETALPHA "Error: already set alpha value"

#define HW_ERROR_BOUNDNOTYETSUPP5THPAR "Error: bounds are not yet supported, fifth input must be []"
#define HW_ERROR_BOUNDNOTYETSUPP6THPAR "Error: bounds are not yet supported, sixth input must be []"
#define HW_ERROR_POPINITRANGEINCSPEC "Error: 'PopInitRange' is incorrectly specified"
#define HW_ERROR_INVINITPOPRANGE "Error: invalid initial population range"
#define HW_ERROR_INITPOPINCSPEC "Error: 'InitialPopulation' is incorrectly specified"
#define HW_ERROR_POPSIZEINCSPEC "Error: 'PopulationSize' is incorrectly specified"
#define HW_ERROR_DESVARNOTINIT "Error: design variables could not be initialized"
#define HW_ERROR_INITSTEMPTOSPECDIM "Error: initial state must be empty to specify dimension"
#define HW_ERROR_1STAND2NDINPBOTHSCALORBOTHVEC "Error: first and second inputs must be both scalars or both vectors"
#define HW_ERROR_NOMORE1OUTALLOWEDW1INP "Error: no more than one output is allowed with 1 input"

#define HW_ERROR_UNSUPFORMTYPE "Error: unsupported format type"
#define HW_ERROR_INVALIDFORMAT "Error: invalid format"

#define HW_ERROR_UNSUPP2DIM "Error: more than 2 dimensions is not supported yet"

#define HW_ERROR_INPSTRONLY1D "Error: input string can only be one-dimensional"
#define HW_ERROR_STRINPMUST1DIM "Error: string inputs must be 1 dimensional; use a cell array for multiple strings"
#define HW_ERROR_TYPE1DSTR "Error: type must be a 1D string"
#define HW_ERROR_STRDIM "Error: all string inputs must be of the same dimensions"

#define HW_ERROR_OUTMEM "Error: out of memory"

#define HW_ERROR_MATSIZE "Error: matrix sizes must agree"
#define HW_ERROR_SIZENOMATCH "Error: sizes do not match"
#define HW_ERROR_MATRIXDIM "Error: matrices must have the same dimensions"
#define HW_ERROR_INCOMPDIM "Error: incompatible dimensions"
#define HW_ERROR_INPMATSAMESIZEALONGSPECDIMEN "Error: input matrices must have the same size along the specified dimension"

#define HW_ERROR_INV1STIND "Error: invalid first index"
#define HW_ERROR_INV2NDIND "Error: invalid second index"
#define HW_ERROR_INVCELLIND "Error: invalid cell index"
#define HW_ERROR_INVCELLINDTYPE "Error: invalid cell index type"
#define HW_ERROR_CELLINDEXRANGE "Error: cell index out of range"
#define HW_ERROR_INDEXRANGE "Error: index out of range"
#define HW_ERROR_INVINDTYPE "Error: invalid index type"
#define HW_ERROR_INDEXPOSINT "Error: index must be a positive integer"
#define	HW_ERROR_INDGRDIM "Error: index exceeds dimension"
#define HW_ERROR_INVIND "Error: invalid index"
#define HW_ERROR_INVINDOP "Error: invalid index operation"

#define HW_ERROR_UNSUPINDOP "Error: unsupported indexing operation"

#define HW_ERROR_DIMFINITE "Error: dimensions must be finite"
#define HW_ERROR_VARIABLENAMETOOLONG "Error: variable name too long"

#define HW_ERROR_MAXFUNCDEPTH "Error: maximum function depth reached"
#define HW_ERROR_INVPERSISKEY "Error: invalid use of persistent keyword"
#define HW_ERROR_STRSAMENUMOFCOLWCOMPROW "Error: strings must have same number of columns when comparing rows"
#define HW_ERROR_SAMENUMOFCOLWCOMPROW "Error: must have same number of columns when comparing rows"
#define HW_ERROR_NOTUSECELLSEARCHFORROW "Error: cannot use cell arrays when searching for rows"
#define HW_ERROR_READPIPEENDFAIL "Error: failed to read the pipe to the end"
#define HW_ERROR_PROBOPENPIPE "Error: problem opening pipe"
#define HW_ERROR_INVTOL "Error: invalid tolerance"
#define HW_ERROR_CONVABSPATH "Error: problem converting to absolute path"
#define HW_ERROR_INPSTRMUSTFILEDIR "Error: input must be a string name referring to a file or directory"
#define HW_ERROR_NOTFINDCURWORKDIR "Error: could not find current working directory"
#define HW_ERROR_NOTLOCALLSUBEXP "Error: problem locating all sub-expressions"
#define HW_ERROR_NOTEMPSUBEXPNAME "Error: cannot have an empty sub-expression name"
#define HW_ERROR_SUBEXPNAMNOCLOSINGGR "Error: sub-expression name in regular expression pattern did not have a closing '>'"
#define HW_ERROR_NOTCONVINPTODOUBLE "Error: cannot convert input to double"
#define HW_ERROR_LENNOTNEG "Error: length cannot be negative"
#define HW_ERROR_SCALOUTCHARRANGE "Error: scalar outside of character range"
#define HW_ERROR_INVINPSTRUCTANDCOMPNOTCONVSTR "Error: invalid input type: structs and complex values cannot be converted to strings"

#define HW_ERROR_MONTHMUSTINT "Error: number of months must be an integer"
#define HW_ERROR_DATENOTCOMP "Error: date cannot be complex"
#define HW_ERROR_PROBCREATTEMPF "Error: problem creating temporary file"
#define HW_ERROR_OPTNOTSING "Error: option is not singular"
#define HW_ERROR_MATMUST3COL "Error: matrix must have 3 columns"
#define HW_ERROR_MATMUST2OR3COL "Error: matrix must have 2 or 3 columns"
#define HW_ERROR_1INPMUSTMAT "Error: single input must be a matrix"
#define HW_ERROR_WHENUNIONINDROWVECSAMECOL "Error: when finding union of individual rows, each vector must have the same number of columns"
#define HW_ERROR_CELLINPMUSTSTR "Error: cell inputs must only contain strings"
#define HW_ERROR_NOUNIONBYROWCELLINP "Error: cannot find union by rows for cell input"
#define HW_ERROR_INVTYPESETELE "Error: invalid type to set all elements to"
#define HW_ERROR_READ "Error: read error"
#define HW_ERROR_NOTCONVNANTOLOG "Error: cannot convert NaN to logical"
#define HW_ERROR_NOTCONVCOMPTOLOG "Error: cannot convert complex to logical"
#define HW_ERROR_MATINPMUSTVEC "Error: all matrix inputs must be vectors"
#define HW_ERROR_3INP1STVEC "Error: when using 3 arguments, the first argument must be a vector"
#define HW_ERROR_INPMUSTSAMESIZE "Error: inputs must be the same size"
#define HW_ERROR_INVBALFLAG "Error: invalid balance flag"

#define HW_ERROR_INVCLASSNAME "Error: invalid class name"
#define HW_ERROR_INVSTRVAL "Error: invalid string value"
#define HW_ERROR_CONVFILEMODTIME "Error: problem converting file modified time"
#define HW_ERROR_TRYSTR "Error: try expression must be a string"
#define HW_ERROR_DIRECT1STLAST "Error: direction must be 'first' or 'last'"
#define HW_ERROR_NUMVALFINDPOSINT "Error: number of values to find must be a positive integer"
#define HW_ERROR_INVINP "Invalid inputs"
#define HW_ERROR_INVINPVAL "Error: invalid input value"
#define HW_ERROR_INVINPTYPE "Error: invalid input type"
#define HW_ERROR_INVINP5THARG "Error: invalid input type for 5th argument"
#define HW_ERROR_INVSTRVALINFNEGINFFRO "Error: invalid string value; must be 'inf', '-inf', or 'fro'"
#define HW_ERROR_INPSQUAREMATVEC "Error: input must be square matrix or vector"
#define HW_ERROR_3RDINPMUSTEMPMATSPMATUNSUPP "Error: third input must be an empty matrix; sparse matrices are not yet supported"

#define HW_ERROR_NOTSETDIRECTTWICE "Error: cannot set direction twice"
#define HW_ERROR_INPONLY1STLASTORROWS "Error: input can only be 'first', 'last', or 'rows'"
#define HW_ERROR_NOTCONVINPTOSTR "Error: error converting input to string"
#define HW_ERROR_DELCELLINVAMTOFELE "Error: delimiter cell array had an invalid amount of elements"
#define HW_ERROR_DELTYPESIMPORREGEX "Error: value for delimitertype must either be 'simple' or 'regularexpression'"
#define HW_ERROR_ "Error: value for collapsedelimiters must be scalar"
#define HW_ERROR_DELTYPEALLSET "Error: the value for delimitertype was already set"
#define HW_ERROR_COLAPSDELALLSET "Error: the value for collapsedelimiters was already set"
#define HW_ERROR_DELONLYSTR "Error: delimiters can only be strings"
#define HW_ERROR_INVMODESTR "Error: invalid mode string"
#define HW_ERROR_NOTFTELLONSTDINOUTERR "Error: cannot call ftell on stdin, stdout, or stderr"
#define HW_ERROR_NOTFSEEKONSTDINOUTERR "Error: cannot call fseek on stdin, stdout, or stderr"
#define HW_ERROR_FRSTINPSTRORFILESTR "Error: first input should be a string of a file name"
#define HW_ERROR_INVINPTYPESIZEINP "Error: invalid input type for size input"
#define HW_ERROR_NOTSKIPBUILTINFUNC "Error: cannot use skip on built-in function streams"
#define HW_ERROR_INVPRECTYPE "Error: invalid precision type"
#define HW_ERROR_PRECBLOCKPOSINT "Error: precision block size must be a positive integer"
#define HW_ERROR_INPNOTSTRUCT "Error: input cannot contain structs"
#define HW_ERROR_INP1DINPMAT3ELE "Error: input at least one dimension of input matrices must have 3 elements"
#define HW_ERROR_BADFUNCHANDLE "Error: bad function handle; was the function cleared or created in a different scope?"
#define HW_ERROR_FUNCNAMENOTSPACE "Error: function name cannot have spaces"
#define HW_ERROR_MISOUTPUTTEMPLATE "Error: missing output template"

#define HW_ERROR_INVOPTMUSTROW "Error: invalid option; must be 'rows'"
#define HW_ERROR_INVOPTMUSTROWORCOL "Error: invalid input string; must be 'rows' or 'cols'"
#define HW_ERROR_INPSMUSTCELLORSTR "Error: inputs must be cell arrays or strings"
#define HW_ERROR_VALNOTSTRUCT "Error: value is not a struct"
#define HW_ERROR_INPVECSORTROW "Error: input must be a vector if not checking sort by rows"
#define HW_ERROR_NOTSETOPTROWTWICE "Error: cannot set option 'rows' twice"
#define HW_ERROR_OUTNOTUNI "Error: outputs were not uniform"
#define HW_ERROR_MISSVALUNIFOUTOPT "Error: missing value for UniformOutput option"
#define HW_ERROR_PATHSEP1CHAR "Error: path separator can only be one character"
#define HW_ERROR_NOTCONCINPTYPE "Error: cannot concatenate input types"
#define HW_ERROR_CONCATSTRUCTWSTRUCT "Error: can only concatenate structs with other structs"
#define HW_ERROR_STRSAMENUMROWCONCAT "Error: strings must have the same number of rows to be concatenated"
#define HW_ERROR_CONCATMATWSCALCOMPMATORCELL "Error: can only concatenate matrices with scalars, complex values, other matrices, or cells"
#define HW_ERROR_INVSHAPEFULLSAMEVALID "Error: invalid shape; must be 'full', 'same', or 'valid'"
#define HW_ERROR_IFCONVMATINPMUSTVEC "Error: if using convolution matrix, other inputs must be vectors"
#define HW_ERROR_DUPFIELD "Error: duplicate field"
#define HW_ERROR_FIELDNAMECELLSTR "Error: field names cell array must only contain strings"
#define HW_ERROR_FIELDNAMEDIMINPCELL "Error: field names did not match dimension of input cell"
#define HW_ERROR_INDINPSAMESIZE "Error: all index inputs must be the same size"
#define HW_ERROR_STRUCTMUSTSAMEFIELDNAME "Error: struct must have the same field names to add to struct array"
#define HW_ERROR_NOTINDINPTYPE "Error: cannot index input type"
#define HW_ERROR_NOTINDTYPEOBJ "Error: cannot index this type of object"
#define HW_ERROR_NOTCELLINDNONCELL "Error: cannot use cell indexing on a non-cell"
#define HW_ERROR_NOTCOMPTOSTR "Error: cannot convert complex value to string"
#define HW_ERROR_NOTSETEMPTYSTRUCT "Error: cannot set the value of an empty struct"

#define HW_ERROR_SETDIMMOREONCE "Error: cannot set dimension more than once"
#define HW_ERROR_SETMODEMOREONCE "Error: cannot set mode more than once"
#define HW_ERROR_NOTSETEXTRAPMOREONCE "Error: cannot set extrap more than once"
#define HW_ERROR_NOTSETMETHODMOREONCE "Error: cannot set method more than once"
#define HW_ERROR_SEEDVALMORETHANONCE "Error: seed value was set more than once"
#define HW_ERROR_STATESETMORETHANONCE "Error: state was set more than once"

#define HW_ERROR_ENVVARFAIL "Error: problem setting environment variable"
#define HW_ERROR_OPTSTRUCTSIZE1 "Error: opts struct must have a size of 1"
#define HW_ERROR_NOTCALLRECURSEEVALIN "Error: cannot call evalin inside a call to evalin"
#define HW_ERROR_CONTBASEORCALL "Error: context must be 'base' or 'caller'"
#define HW_ERROR_2NDINPFINITEVAL "Error: second input must contain finite values"

#define HW_ERROR_TEMPLATESTR "Error: template must be a string"
#define HW_ERROR_INPUTSTRING "Error: input must be a string"
#define HW_ERROR_INPUTALLSTRINGS "Error: inputs must all be strings"
#define HW_ERROR_FUNCNAMESTR "Error: function name must be a string"
#define HW_ERROR_FIELDNAMESTR "Error: field names must be strings"
#define HW_ERROR_FIELDNAMENOTEMPTYSTR "Error: field name cannot be an empty string"
#define HW_ERROR_EVERYVALFIELD "Error: must have a value for every field"

#define HW_ERROR_4THINPBOOL "Error: 4th argument must be true or false"
#define HW_ERROR_2NDINPCELLSTR "Error: second input must be a cell array or string"
#define HW_ERROR_FUNCHANDLNAMEINPUT "Error: input must be a function handle or function name"
#define HW_ERROR_INVTYPEVALDOTPARCURL "Error: invalid type value; must be '.', '()', or '{}'"
#define HW_ERROR_SUBMUSTCELL "Error: subs must be ':' or a cell array"
#define HW_ERROR_NOTINDSTRUCT1X1 "Error: cannot index struct array; must be 1x1"
#define HW_NONSTRUCTFIELD "Error: can't get field value of non-struct value"
#define HW_ERROR_UNSUPINDCHAIN "Error: no support for index chaining"
#define HW_ERROR_INVSTRUCTFIELD "Error: invalid struct field"
#define HW_ERROR_INVSTRUCTFIELDTYPESUB "Error: invalid struct fields; must be 'type' and 'subs'"
#define HW_ERROR_ASSERTFAIL "Error: assert failed"
#define HW_ERROR_INVFUNCNARGIN "Error: invalid function to call nargin on"
#define HW_ERROR_INVFUNCNAME "Error: invalid function name"
#define HW_ERROR_INVDIR "Error: invalid directory name"
#define HW_ERROR_INPROWVECT "Error: inputs must be row vectors"
#define HW_ERROR_PROBCHANGCURDIR "Error: problem changing current directory"
#define HW_ERROR_NOFMATCH "Error: no files matching the pattern were found"
#define HW_ERROR_INVFNAMEPAT "Error: invalid file name pattern"
#define HW_ERROR_2NDARGONLYCLR "Error: second argument can only be 'clear'"
#define HW_ERROR_PROBLOADTREE "Error: problem loading tree"
#define HW_ERROR_PROBSERTREE "Error: problem serializing tree"
#define HW_ERROR_INVCELLEXT "Error: invalid cell extraction"
#define HW_ERROR_INVAST "Error: Invalid AST"
#define HW_ERROR_NOCONTEXTENDFUNC "Error: no available context for end function"
#define HW_ERROR_INVFUNCCALLNARGOUT "Error: invalid function to call nargout on"
#define HW_ERROR_INVCALLNARGOUT "Error: invalid call to nargout"
#define HW_ERROR_INVCALLNARGIN "Error: invalid call to nargin"
#define HW_ERROR_INVCONTNARGOUT "Error: invalid context for nargout"
#define HW_ERROR_INVCONTNARGIN "Error: invalid context for nargin"
#define HW_ERROR_CANNOTCHECKLOCKNOTFUNC "Error: not in a function; cannot check if locked"
#define HW_ERROR_UNLOCKFAILNOTFUNC "Error: not in a function; unlock failed"
#define HW_ERROR_LOCKFAILNOTFUNC "Error: not in a function; lock failed"
#define HW_ERROR_ILLASSIGN "Error: illegal assignment"
#define HW_ERROR_FILENOTEXEC "Error: OML file cannot have executable statements"
#define HW_ERROR_INVFIELDNAME "Error: invalid field name"
#define HW_ERROR_INVCELLRESIZE "Error: invalid cell resize"
#define HW_ERROR_UNKNOWNLHS "Error: unknown LHS"
#define HW_ERROR_UNKNOWNVAR "Error: unknown variable"
#define HW_ERROR_INVSTARTVAL "Error: invalid start value"
#define HW_ERROR_INVENDVAL "Error: invalid end value"
#define HW_ERROR_INVCELLDEL "Error: invalid cell deletion"
#define HM_ERROR_INDSAMESIZE "Error: index must be same size as RHS"
#define HW_ERROR_INVRHS "Error: invalid RHS"
#define HW_ERROR_RHSSIZEIND "Error: index and RHS must be the same size"
#define HW_ERROR_INVRHSSIZE "Error: invalid RHS size"
#define HW_ERROR_INVMATRESIZE "Error: invalid matrix resize"
#define HW_ERROR_INVRIGHT "Error: invalid right hand side"
#define HW_ERROR_UNASSIGNEMPTRIGHT "Error: unable to assign empty right hand side"
#define HW_ERROR_NOTUSEKEY "Error: cannot use reserved keyword"
#define HW_ERROR_NOTMIXCOMPSTR "Error: cannot mix complex and strings"
#define HW_ERROR_NOTCREATEMATRIX "Error: could not create matrix"
#define HW_ERROR_NOTCREATECELLARRAY "Error: could not create cell array"
#define HW_ERROR_INVMATRIXINP "Error: invalid matrix input"
#define HW_ERROR_UNKOWNANTLR "Error: unknown ANTLR error"
#define HW_ERROR_INVSWITCH "Error: invalid switchcase input"
#define HW_ERROR_UNKOWNFUN "Error: unknown function"
#define HW_ERROR_UNKOWNTYPE "Error: unknown type"
#define HW_ERROR_INVFUNCALL "Error: invalid function call"
#define HW_ERROR_LOGOPCOMP "Error: cannot perform logical operation on a complex value"
#define HW_ERROR_COMPMATUNSUP "Error: complex matrices not supported yet"
#define HW_ERROR_MISSRETURNS "Error: missing return values"
#define HW_ERROR_MISSCELLVAL "Error: missing cell values"
#define HW_ERROR_NOTIMP "Error: not implemented yet"
#define HW_ERROR_UNSUPRANGEOP "Error: unsupported range operation"
#define HW_ERROR_UNSUPOP "Error: unsupported operation"
#define HW_ERROR_UNSUPCOMP "Error: unsupported comparision"

#define HW_ERROR_INPUTSTRUCT  "Error: invalid input type; must be struct"
#define HW_ERROR_INPUTSCALARCOMPLEXMATRIX  "Error: invalid input type; must be scalar, complex, or a matrix"
#define HW_ERROR_INPUTSCALARMATRIX  "Error: invalid input type; must be scalar or a matrix"
#define HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING  "Error: invalid input type; must be scalar, complex, a matrix, or a string"
#define HW_ERROR_INPUTSCALARMTXSTRING  "Error: invalid input type; must be scalar, a matrix, or a string"
#define HW_ERROR_INPUTSCALARMTXSTRINGCELL  "Error: invalid input type; must be a scalar, matrix, cell array, or a string"
#define HW_ERROR_INPUTSTRCELLMTX "Error: input must be a scalar, complex, matrix, cell array, or string"
#define HW_ERROR_INPUTCELLARRAY "Error: input must be a cell array"
#define HW_ERROR_INPUTHCATDIM "Error: inputs must have the same number of rows to be concatenated horizontally"
#define HW_ERROR_INPUTVCATDIM "Error: inputs must have the same number of columns to be concatenated vertically"
#define HW_ERROR_INPUTCATSTRUCT "Error: structs must have the same field names to be concatenated"
#define HW_ERROR_INPUTSTRINGCELLARRAY "Error: all inputs must be strings or cell arrays of strings"
#define HW_ERROR_UNEVENDIMENSIONS "Error: input did not make even dimensions; matrix cannot be constructed"
#define HW_ERROR_MIXEDCELLELEMS "Error: cell array contains mixed element types; cannot concatenate them"
#define HW_ERROR_DIM3ELEM "Error: matrices must have 3 elements in specified direction"
#define HW_ERROR_ADDPATHLOC "Error: append location must be 0, 1, '-begin', or '-end'"
#define HW_ERROR_INPUTNOCONJ "Error: not all inputs had matching conjugates"

#define HW_ERROR_INPUTALLCELL "Error: all inputs must be cell arrays"
#define HW_ERROR_INPUTREALSTR "Error: input must be a scalar, real matrix, or string"
#define HW_ERROR_INPSCSTR "Error: input must be a scalar or string"
#define HW_ERROR_INPUTSTRINGFUNC "Error: function name must be a string or function handle"
#define HW_ERROR_CELLELEMSTR "Error: cell array elements must be strings"
#define HW_ERROR_INPUTISLOGICAL "Error: input cannot be logical"

#define HW_ERROR_OPTIONSTRING "Error: option must be a string"
#define HW_ERROR_CELLINPSAMESASIZE "Error: all cell array inputs must be the same size"
#define HW_ERROR_CELLARRAYSSAMESIZE "Error: cell arrays must be the same size"
#define HW_ERROR_INDEXSTRUCTFIELDS "Error: cannot index list of struct fields"
#define HW_ERROR_DOTSTRUCTPARENTHESIS "Error: invalid type value; must be '.', '()', or '{}'"
#define HW_ERROR_REQCOMPSTARTINDEXLESS "Error: Comp and Req start index must be less than their respective end indexes"
#define HW_ERROR_INPUTONOROFF "Error: Input must be 'on', 'off', or a text file name"
#define HW_ERROR_ILLEGALCHARFORFILE "Error: file name input contains an illegal character"
#define HW_ERROR_DYNFIELD "Error: dynamic field must be a string"
// End deprecated #define list

#endif // _OML_Math_Error_h
