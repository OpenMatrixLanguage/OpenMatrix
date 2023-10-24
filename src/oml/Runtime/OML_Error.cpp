/**
* @file OML_Error.cpp
* @date November 2015
* Copyright (C) 2015-2022 Altair Engineering, Inc.  
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
#define OML_MSG_SCALARSTRING  "Error: invalid input; must be scalar or string"
#define OML_MSG_INTSTRING     "Error: invalid input; must be integer or string"
#define OML_MSG_POSINTALL     "Error: invalid input; must be positive integer or 'all'"
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
#define OML_MSG_POSINTMATRIX  "Error: invalid input; must be positive integer matrix"
#define OML_MSG_MATRIX        "Error: invalid input; must be a matrix"
#define OML_MSG_REALMATRIX    "Error: invalid input; must be a real matrix"
#define OML_MSG_EMPTYMATRIX   "Error: invalid input; must be empty [] matrix"
#define OML_MSG_VEC_2DMAT     "Error: invalid input; must be a vector or 2D matrix"
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
#define OML_MSG_ARGDIMENSIONS "Error: incompatible argument dimensions"
#define OML_MSG_ARRAYSIZE     "Error: incompatible matrices; dimensions must be consistent"
#define OML_MSG_ARRAYCATDIM   "Error: incompatible array sizes; non-concatenated dimensions must match"
#define OML_MSG_CELLSIZE      "Error: incompatible cell sizes; must match"
#define OML_MSG_OPTION        "Error: invalid option argument"
#define OML_MSG_OPTIONVAL     "Error: invalid option; is incorrectly specified"
#define OML_MSG_FUNCSWITCH    "Error: invalid option; must be either 'on' or 'off'"
#define OML_MSG_NOBUILTIN     "Error: built in function not supported in this context"
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
#define OML_MSG_ALLOCFAILED                 "Error: allocation failure"
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
#define OML_MSG_SCALAR_REALVEC              "Error: invalid input; must be a scalar or real vector"
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
#define OML_MSG_ERR_OBJ_HANDLE              "Error: invalid input: must be an object handle"
#define OML_MSG_ERR_HANDLE_STRING           "Error: invalid input: must be an object handle or string"
#define OML_MSG_INTERNAL                    "Error: internal error"
#define OML_MSG_AUTHENTICATE                "Error: authentication failure"
#define OML_MSG_UNICODE_FILENAME            "Error: invalid input: file name cannot have Unicode characters"
#define OML_MSG_STRING_POSINTEGER           "Error: invalid input: must be a string or positive integer"

// optimization messages
#define OML_MSG_OBJSTRRET1                  "Error: invalid objective function; must have exactly 1 return"
#define OML_MSG_CONSTRARG2                  "Error: invalid constraint function; argument 2 must be []"
#define OML_MSG_CONSTRRET2                  "Error: invalid constraint function; can have at most 2 returns"
#define OML_MSG_CONSTRRET4                  "Error: invalid constraint function; can have at most 4 returns"
#define OML_MSG_ANALYTICGRADS               "Error: invalid options; GradObj and GradConstr must be used together"

// signal processing messages
#define OML_MSG_ESTIMATOR_STR               "Error: invalid argument; 'Estimator' is expected"
#define OML_MSG_ESTIMATOR_TYPE              "Error: invalid estimator; 'H1', 'H2' or 'H3' is expected"

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
#define OML_MSG_PLOT_CONST_PROPERTY             "Error: property is read-only; cannot update"
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

// HDF5 reader messages
#define OML_MSG_HDF5_INVALID_FILE                "Error: not a valid hdf5 file"
#define OML_MSG_HDF5_FILE_READ_FAILED            "Error: failed to read hdf5 file"
#define OML_MSG_HDF5_INVALID_GROUP               "Error: input is not a group"
#define OML_MSG_HDF5_GROUP_READ_FAILED           "Error: failed to read group"
#define OML_MSG_HDF5_INVALID_DATASET             "Error: input is not a dataset"
#define OML_MSG_HDF5_DATASET_READ_FAILED         "Error: failed to read dataset"
#define OML_MSG_HDF5_NEITHER_DATASET_NOR_GROUP   "Error: location is invalid; location must be a dataset or a group "
#define OML_MSG_HDF5_ATTRIBUTE_READ_FAILED       "Error: failed to read attributes"
#define OML_MSG_HDF5_DATATYPE_READ_FAILED        "Error: failed to read datatype"
#define OML_MSG_HDF5_DATASPACE_READ_FAILED       "Error: failed to read dataspace"
#define OML_MSG_HDF5_POINTS_MATRIXDIM_INVALID    "Error: number of columns in points selection matrix does not match number of dimensions in data"
#define OML_MSG_HDF5_POINTS_SELECTION_INVALID    "Error: selected points are out of range"
#define OML_MSG_HDF5_UNSUPPORTDIM                "Error: data with more than seven dimentions are not supported"
#define OML_MSG_HDF5_UNSUPPORT_DATA              "Error: unsupported data"
#define OML_MSG_HDF5_FILE_CREATION_FAILED        "Error: failed to create HDF5 file"
#define OML_MSG_HDF5_GROUP_EXIST                 "Error: group already exists"
#define OML_MSG_HDF5_GROUP_CREATION_FAILED       "Error: failed to create group"
#define OML_MSG_HDF5_GROUP_RENAME_FAILED         "Error: failed to rename group"
#define OML_MSG_HDF5_GROUP_REMOVE_FAILED         "Error: failed to remove group"
#define OML_MSG_HDF5_DATASET_CREATION_FAILED     "Error: failed to create dataset"
#define OML_MSG_HDF5_WRITE_FAILED                "Error: failed to write"
#define OML_MSG_HDF5_DATASET_RENAME_FAILED       "Error: failed to rename dataset"
#define OML_MSG_HDF5_DATASET_REMOVE_FAILED       "Error: failed to remove dataset"
#define OML_MSG_HDF5_DATASET_EXIST               "Error: dataset already exists"
#define OML_MSG_HDF5_CHUNK_MATRIXDIM_INVALID     "Error: number of columns in chunk size matrix must match the number of dimensions in data"
#define OML_MSG_HDF5_CHUNK_MATRIX_INVALID_ROW    "Error: number of rows in chunk size matrix must be 1"
#define OML_MSG_HDF5_EMPTYMATRIX                 "Error: must not be empty matrix"
#define OML_MSG_HDF5_EMPTYCELL                   "Error: must not be empty cell"
#define OML_MSG_HDF5_INVALID_CELL_MEMBERS        "Error: all cell members must have same dimensions"
#define OML_MSG_HDF5_STRUCT_INVALID_DIMS         "Error: invalid struct dimensions; only single element structs are supported"
#define OML_MSG_HDF5_ATTRIBUTE_CREATION_FAILED   "Error: failed to create attribute"
#define OML_MSG_HDF5_INVALID_LOCATION            "Error: invalid location"
#define OML_MSG_HDF5_ATTRIBUTE_EXIST             "Error: attribute already exists"
#define OML_MSG_HDF5_INVALID_ATTRIBUTE           "Error: invalid attribute"
#define OML_MSG_HDF5_INVALID_STRUCT_MEMBERS      "Error: all values in struct must have same dimensions"
#define OML_MSG_HDF5_INVALID_ATTRIBUTE_DATA      "Error: attribute value must be a scalar,string or real matrix."

// oml menu apis codes
#define OML_MSG_MENUAPI_SYS_NAME                 "Error: cannot use system ribbon name."
#define OML_MSG_MENUAPI_SYS_DEL                  "Error: cannot delete system ribbons."
#define OML_MSG_MENUAPI_SYS_MOD                  "Error: cannot modify system ribbons."
#define OML_MSG_MENUAPI_INVALID_HANDLE           "Error: invalid handle"
#define OML_MSG_MENUAPI_INVALID_HANDLE_TYPE      "Error: invalid handle type"

// Image codes
#define OML_MSG_IMAGE_ALPHA                      "Error: invalid handle; cannot retrieve alpha data"
#define OML_MSG_IMAGE_GRAYSCALE                  "Error: invalid handle; cannot retrieve grayscale data"
#define OML_MSG_IMAGE_RGB                        "Error: invalid handle; cannot retrieve RGB data"
#define OML_MSG_IMAGE_HANDLE                     "Error: invalid handle; cannot retrieve image data"

#define OML_MSG_REALVECTOR2                      "Error: invalid input; must be a  real vector with 2 elements"
#define OML_MSG_REALVECTOR3                      "Error: invalid input; must be a  real vector with 3 elements"

// RPC Writer
#define OML_MSG_RPC_INVALID_TIME_MAT             "Error: time data should be 1xN"
#define OML_MSG_RPC_INVALID_DATA_MAT             "Error: at least one row should be in the data matrix"
#define OML_MSG_RPC_INVALID_TIME_DATA_MAT        "Error: number of columns in time and data matrix should match"
#define OML_MSG_RPC_INVALID_HEADER               "Error: invalid header"
#define OML_MSG_RPC_WRITE_FAILED                 "Error: failed to write"
#define OML_MSG_RPC_INVALID_CURVENAME            "Error: curve name should be a string"
#define OML_MSG_RPC_LOAD_FAILED                  "Error: faild to load writer"
#define OML_MSG_RPC_XVECT_NUM_POINTS             "Error: X vector should have at least two data points"
#define OML_MSG_RPC_XVECT_SPREAD                 "Error: X vector is not evenly distributed"

//ReadCAE Builder && HW reader messages
#define OML_MSG_HWREADER_INVALID_TYPE            "Error: invalid datatype"
#define OML_MSG_HWREADER_INVALID_REQUEST         "Error: invalid requests"
#define OML_MSG_HWREADER_MISSING_COMPONENT       "Error: missing component"
#define OML_MSG_HWREADER_INVALID_COMPONENT       "Error: invalid component"
#define OML_MSG_HWREADER_MISSING_TIME            "Error: missing time"
#define OML_MSG_HWREADER_FAIL_READER_INIT        "Error: reader initialization failed"
#define OML_MSG_HWREADER_READING_FAILED          "Error: problem reading file"
#define OML_MSG_HWREADER_INVALID_FILE_ENTRY      "Error: filename cannot be an empty string. Specify a valid filename"
#define OML_MSG_HWREADER_TIMECHANNELS_COMPARE    "Error: time channels does not match"
#define OML_MSG_HWREADER_SUBCASE_INVALID_RANGE   "Error: invalid input; subcase index out of range;"
#define OML_MSG_HWREADER_INVALID_TIME            "Error: invalid time"
#define OML_MSG_HWREADER_NEED_SUBCASE            "Error: must specify a subcase for files with subcases"
#define OML_MSG_HWREADER_NO_SUBCASE              "Error: file has no subcases"
#define OML_MSG_HWREADER_STRING_POSINTEGER       "Error: invalid input; must be a string or an integer"
#define OML_MSG_HWREADER_INVALID_TRC             "Error: invalid type, request, and component combination"
#define OML_MSG_HWREADER_SUBCASE_READ_FAIL       "Error: problem getting subcases from file"
#define OML_MSG_HWREADER_TYPE_READ_FAIL          "Error: problem getting data types from file"
#define OML_MSG_HWREADER_REQUEST_READ_FAIL       "Error: problem getting data requests from file"
#define OML_MSG_HWREADER_COMP_READ_FAIL          "Error: problem getting components from file"
#define OML_MSG_HWREADER_FILTER_COMP_READ_FAIL   "Error: problem getting filtered components from file"
#define OML_MSG_HWREADER_INVALID_TYPE_INDX       "Error: invalid datatype index"
#define OML_MSG_HWREADER_INVALID_REQUEST_INDX    "Error: invalid request index"
#define OML_MSG_HWREADER_INVALID_COMPONENT_INDX  "Error: invalid component index"
#define OML_MSG_HWREADER_CELL_INDX               "Error: index must be less than the number of rows in the input cell"
#define OML_MSG_HWREADER_INDX_POSINTEGER         "Error: index must be a positive integer"
#define OML_MSG_HWREADER_EXTRACT_STR_SCALAR      "Error: extract fields must either be strings or indices"
#define OML_MSG_HWREADER_DIMENSIONS_MATCH        "Error: dimensions must match"
#define OML_MSG_HWREADER_INDICES_REAL            "Error: indices must be real"
#define OML_MSG_HWREADER_TRC_DATA_TYPE           "Error: last three inputs must be Type, Req, and Comp strings or integers"
#define OML_MSG_HWREADER_STRC_DATA_TYPE          "Error: last four inputs must be Subcase, Type, Req, and Comp strings or integers"
#define OML_MSG_HWREADER_CELL_4_5_COLS           "Error: cell array input must have 4 or 5 columns"
#define OML_MSG_HWREADER_CELL_3_4_COLS           "Error: cell array input must have 3 or 4 columns"
#define OML_MSG_HWREADER_REQ_CELL_1_N            "Error: request input cell should be 1,n size"
#define OML_MSG_HWREADER_COMP_CELL_1_N           "Error: component input cell should be 1,n size"
#define OML_MSG_HWREADER_TIME_CELL_1_N           "Error: time input cell should be 1,n size"
#define OML_MSG_HWREADER_OUTPUT_CELL_MAT         "Error: invalid output type requested. Output can be cell or matrix. Valid input value 0 or 1"
#define OML_MSG_HWREADER_INVALID_UNIT_TYPE       "Error: invalid unit type"

#define OML_MSG_SPREADSHEET_RANGE                "Error: invalid input; must be a spreadsheet style range"
#define OML_MSG_WINDOWSONLY_OP                   "Error: invalid operation; operation is only valid on Windows"
#define OML_MSG_FILEALREADYOPEN_OP               "Error: invalid operation; file is already open in other application(s)"
#define OML_MSG_OMLINTERFACE_OP                  "Error: invalid operation; cannot complete operation with 'oml'(default) interface"

#define OML_MSG_MQTT_INVALID_OPTION              "Error: invalid option"
#define OML_MSG_MQTT_CLIENTID_INUSE              "Error: client id is in use"
#define OML_MSG_MQTT_CLIENTID_INVALID            "Error: no client exists with given client id"
#define OML_MSG_MQTT_CLIENTCREATION_FAIL         "Error: failed to create client"
#define OML_MSG_MQTT_CLIENTID_SIZE               "Error: invalid client id; client id length must be less than 24"
#define OML_MSG_MQTT_MIN_IDLETIME                "Error: invalid idle time; idle time must be at least 0.001 seconds"


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
#define OML_STR_RANDSEED        "Seed"
#define OML_STR_TOLX            "TolX"
#define OML_STR_TOLFUN          "TolFun"
#define OML_STR_TOLFUNABS       "TolFunAbs"
#define OML_STR_TOLFUNREL       "TolFunRel"
#define OML_STR_TOLCON          "TolCon"
#define OML_STR_CONRET          "Constrant Retention"
#define OML_STR_MOVELIM         "Move Limit Fraction"
#define OML_STR_PERTM           "Perturbation Method"
#define OML_STR_PERTV           "Initial Perturbation Value"
#define OML_STR_POPSIZE         "Population Size"
#define OML_STR_INITPOP         "Initial Population"
#define OML_STR_CRODIST         "Crowding Distance"
#define OML_STR_INITSAMPNTS     "Number of Initial Sample Points"
#define OML_STR_MAXFAIL         "Maximum Failed Iterations"
#define OML_STR_METHOD          "Method"
#define OML_STR_PNTSPERITER     "Points Per Iteration"
#define OML_STR_STOPNOIMPR      "Iterations with no Improvement"
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
std::string OML_Error::GetOmlErrorMessage(omlMathErrCode code) const
{
    switch (code)
    {
        case OML_ERR_NUMARGIN:                    return OML_MSG_NUMARGIN;
        case OML_ERR_NUMARGOUT:                   return OML_MSG_NUMARGOUT;
        case OML_ERR_NUMARGINOUT:                 return OML_MSG_NUMARGINOUT;
        case OML_ERR_CELL:                        return OML_MSG_CELL;
        case OML_ERR_CELLARRAY:                   return OML_MSG_CELLARRAY;
        case OML_ERR_STRUCT:                      return OML_MSG_STRUCT;
        case OML_ERR_STRING:                      return OML_MSG_STRING;
        case OML_ERR_SCALARSTRING:                return OML_MSG_SCALARSTRING;
        case OML_ERR_INTSTRING:                   return OML_MSG_INTSTRING;
        case OML_ERR_POSINTALL:                   return OML_MSG_POSINTALL;
        case OML_ERR_BAD_STRING:                  return OML_MSG_BAD_STRING;
        case OML_ERR_NUMERIC:                     return OML_MSG_NUMERIC;
        case OML_ERR_SCALAR:                      return OML_MSG_SCALAR;
        case OML_ERR_VECTOR:                      return OML_MSG_VECTOR;
        case OML_ERR_VECTOR2:                     return OML_MSG_VECTOR2;
        case OML_ERR_SCALARVECTOR:                return OML_MSG_SCALARVECTOR;
        case OML_ERR_SCALARMATRIX:                return OML_MSG_SCALARMATRIX;
        case OML_ERR_SCALARCOMPLEXMTX:            return OML_MSG_SCALARCOMPLEXMTX;
        case OML_ERR_STRINGVECTOR:                return OML_MSG_STRINGVECTOR;
        case OML_ERR_INTVECSTR:                   return OML_MSG_INTVECSTR;
        case OML_ERR_REALVECTOR:                  return OML_MSG_REALVECTOR;
        case OML_ERR_NNINTVECTOR:                 return OML_MSG_NNINTVECTOR;
        case OML_ERR_POSINTVECTOR:                return OML_MSG_POSINTVECTOR;
        case OML_ERR_POSINTMATRIX:                return OML_MSG_POSINTMATRIX;
        case OML_ERR_MATRIX:                      return OML_MSG_MATRIX;
        case OML_ERR_REALMATRIX:                  return OML_MSG_REALMATRIX;
        case OML_ERR_EMPTYMATRIX:                 return OML_MSG_EMPTYMATRIX;
        case OML_ERR_VEC_2DMAT:                   return OML_MSG_VEC_2DMAT;
        case OML_ERR_HANDLE:                      return OML_MSG_HANDLE;
        case OML_ERR_HANDLE_EMPTY:                return OML_MSG_HANDLE_EMPTY;
        case OML_ERR_FUNCNAME:                    return OML_MSG_FUNCNAME;
        case OML_ERR_ACCUMFUNC:                   return OML_MSG_ACCUMFUNC;
        case OML_ERR_REAL:                        return OML_MSG_REAL;
        case OML_ERR_INTEGER:                     return OML_MSG_INTEGER;
        case OML_ERR_NATURALNUM:                  return OML_MSG_NATURALNUM;
        case OML_ERR_POSINTEGER:                  return OML_MSG_POSINTEGER;
        case OML_ERR_FINITE:                      return OML_MSG_FINITE;
        case OML_ERR_VECLENDIM:                   return OML_MSG_VECLENDIM;
        case OML_ERR_ARGDIMENSIONS:               return OML_MSG_ARGDIMENSIONS;
        case OML_ERR_ARRAYSIZE:                   return OML_MSG_ARRAYSIZE;
        case OML_ERR_ARRAYCATDIM:                 return OML_MSG_ARRAYCATDIM;
        case OML_ERR_CELLSIZE:                    return OML_MSG_CELLSIZE;
        case OML_ERR_OPTION:                      return OML_MSG_OPTION;
        case OML_ERR_OPTIONVAL:                   return OML_MSG_OPTIONVAL;
        case OML_ERR_FUNCSWITCH:                  return OML_MSG_FUNCSWITCH;
        case OML_ERR_NOBUILTIN:                   return OML_MSG_NOBUILTIN;
        case OML_ERR_UNSUPPORTDIM:                return OML_MSG_UNSUPPORTDIM;
        case OML_ERR_FLAG_01:                     return OML_MSG_FLAG_01;
        case OML_ERR_FORMAT:                      return OML_MSG_FORMAT;
        case OML_ERR_PNORM:                       return OML_MSG_PNORM;
        case OML_ERR_NORM_STRING3:                return OML_MSG_NORM_STRING3;
        case OML_ERR_NOCOMPLEX:                   return OML_MSG_NOCOMPLEX;
        case OML_ERR_NATURALNUM_MATRIX_CELLARRAY: return OML_MSG_NATURALNUM_MATRIX_CELLARRAY;
        case OML_ERR_STRING_MATRIX_CELLARRAY:     return OML_MSG_STRING_MATRIX_CELLARRAY;
        case OML_ERR_STRING_ONEDIMENSION:         return OML_MSG_STRING_ONEDIMENSION;
        case OML_ERR_STRSCALARCOMPLEXMTX:         return OML_MSG_STRSCALARCOMPLEXMTX;
        case OML_ERR_FINITE_NATURALNUM:           return OML_MSG_FINITE_NATURALNUM;
        case OML_ERR_STRING_NATURALNUM:           return OML_MSG_STRING_NATURALNUM;
        case OML_ERR_POSITIVE_SCALAR:             return OML_MSG_POSITIVE_SCALAR;
        case OML_ERR_SCALAR_COMPLEX:              return OML_MSG_SCALAR_COMPLEX;
        case OML_ERR_STRING_STRINGCELL:           return OML_MSG_STRING_STRINGCELL;
        case OML_ERR_ALLOCFAILED:                 return OML_MSG_ALLOCFAILED;
        case OML_ERR_INVALID_INDEX:               return OML_MSG_INVALID_INDEX;
        case OML_ERR_INVALID_RANGE:               return OML_MSG_INVALID_RANGE;
        case OML_ERR_INVALID_BASE:                return OML_MSG_INVALID_BASE;
        case OML_ERR_INVALID_DLL:                 return OML_MSG_INVALID_DLL;
        case OML_ERR_INVALID_INITDLL:             return OML_MSG_INVALID_INITDLL;
        case OML_ERR_GUI_CMDEXEC_FAIL:            return OML_MSG_GUI_CMDEXEC_FAIL;
        case OML_ERR_POS_INTEGER_VEC_MTX:         return OML_MSG_POS_INTEGER_VEC_MTX;
        case OML_ERR_STRING_INTEGER:              return OML_MSG_STRING_INTEGER;
        case OML_ERR_SCALAR_VECTOR_STRING:        return OML_MSG_SCALAR_VECTOR_STRING;
        case OML_ERR_POS_INTEGER_MTX_INF:         return OML_MSG_POS_INTEGER_MTX_INF;
        case OML_ERR_POS_INTEGER_MTX_SIZE2:       return OML_MSG_POS_INTEGER_MTX_SIZE2;
        case OML_ERR_CELLMTXSTRUCT:               return OML_MSG_CELLMTXSTRUCT;
        case OML_ERR_CELLSTRING:                  return OML_MSG_CELLSTRING;
        case OML_ERR_SCALAROUTOFRANGE:            return OML_MSG_SCALAROUTOFRANGE;
        case OML_ERR_INVALIDSTRUCTINDEX:          return OML_MSG_INVALIDSTRUCTINDEX;
        case OML_ERR_INVALIDINDEX:                return OML_MSG_INVALIDINDEX;
        case OML_ERR_TRIANGMATTYPE:               return OML_MSG_TRIANGMATTYPE;
        case OML_ERR_MTXSTRING:                   return OML_MSG_MTXSTRING;
        case OML_ERR_POSINTEGER_VEC:              return OML_MSG_POSINTEGER_VEC;
        case OML_ERR_NONEMPTY_STR:                return OML_MSG_NONEMPTY_STR;
        case OML_ERR_ONEROW_STRING:               return OML_MSG_ONEROW_STRING;
        case OML_ERR_SCALAR_REALVEC:              return OML_MSG_SCALAR_REALVEC;
        case OML_ERR_SCALAR_REALMTX:              return OML_MSG_SCALAR_REALMTX;
        case OML_ERR_INTEGER_INTMTX:              return OML_MSG_INTEGER_INTMTX;
        case OML_ERR_LOGICAL:                     return OML_MSG_LOGICAL;
        case OML_ERR_INVALID_VERSION:             return OML_MSG_INVALID_VERSION;
        case OML_ERR_STRING_FILESTREAM:           return OML_MSG_STRING_FILESTREAM;
        case OML_ERR_FILE_FILESTREAM:             return OML_MSG_FILE_FILESTREAM;
	    case OML_ERR_FILE_NOTFOUND:               return OML_MSG_FILE_NOTFOUND;
	    case OML_ERR_FILE_CANNOTOPEN:             return OML_MSG_FILE_CANNOTOPEN;
	    case OML_ERR_OPT_UNSUPPORTED:             return OML_MSG_OPT_UNSUPPORTED;
	    case OML_ERR_INVALIDTIMERANGE:            return OML_MSG_INVALIDTIMERANGE;
	    case OML_ERR_INVALIDTIMESTEP:             return OML_MSG_INVALIDTIMESTEP;
        case OML_ERR_NONNEGATIVE_SCALAR:          return OML_MSG_NONNEGATIVE_SCALAR;
        case OML_ERR_FILE_CANNOTREAD:             return OML_MSG_FILE_CANNOTREAD;
        case OML_ERR_STRINGSTRUCT:                return OML_MSG_STRINGSTRUCT;
        case OML_ERR_SCAL_COMP_MTX_STR_CELL:      return OML_MSG_SCAL_COMP_MTX_STR_CELL;
        case OML_ERR_INVALID_DATE_FMT:            return OML_MSG_INVALID_DATE_FMT;
        case OML_ERR_CANNOTAPPLY_DATE_FMT:        return OML_MSG_CANNOTAPPLY_DATE_FMT;
	    case OML_ERR_INPUT_EMPTY:                 return OML_MSG_INPUT_EMTPY;
        case OML_ERR_MULTILINE_STRING:            return OML_MSG_MULTILINE_STRING;
        case OML_ERR_INVALID_FILE_MODE:           return OML_MSG_INVALID_FILE_MODE;
        case OML_ERR_FILE_CANNOTWRITE:            return OML_MSG_FILE_CANNOTWRITE;
        case OML_ERR_CELL_MTX:                    return OML_MSG_CELL_MTX;
        case OML_ERR_INVALIDTIMEVAL:              return OML_MSG_INVALIDTIMEVAL;
        case OML_ERR_HANDLE_STRING_CELL:          return OML_MSG_HANDLE_STRING_CELL;
        case OML_ERR_OBJ_HANDLE:                  return OML_MSG_ERR_OBJ_HANDLE;
        case OML_ERR_HANDLE_STRING:               return OML_MSG_ERR_HANDLE_STRING;
        case OML_ERR_INTERNAL:                    return OML_MSG_INTERNAL;
        case OML_ERR_AUTHENTICATE:                return OML_MSG_AUTHENTICATE;
        case OML_ERR_UNICODE_FILENAME:            return OML_MSG_UNICODE_FILENAME;
        case OML_ERR_STRING_POSINTEGER:           return OML_MSG_STRING_POSINTEGER;

        // optimization error messages:
        case OML_ERR_OBJSTRRET1:                  return OML_MSG_OBJSTRRET1;
        case OML_ERR_CONSTRARG2:                  return OML_MSG_CONSTRARG2;
        case OML_ERR_CONSTRRET2:                  return OML_MSG_CONSTRRET2;
        case OML_ERR_CONSTRRET4:                  return OML_MSG_CONSTRRET4;
        case OML_ERR_ANALYTICGRADS:               return OML_MSG_ANALYTICGRADS;

        // signal processing error messages:
        case OML_ERR_ESTIMATOR_STR:               return OML_MSG_ESTIMATOR_STR;
        case OML_ERR_ESTIMATOR_TYPE:              return OML_MSG_ESTIMATOR_TYPE;
        
        // plot error messages:
        case OML_ERR_PLOT_DIM_NOT_MATCH:          return OML_MSG_PLOT_DIM_NOT_MATCH;
        case OML_ERR_PLOT_MISSING_VALUE:          return OML_MSG_PLOT_MISSING_VALUE;
        case OML_ERR_PLOT_UNMATCHED_AXES_TYPE:    return OML_MSG_PLOT_UNMATCHED_AXES_TYPE;
        case OML_ERR_PLOT_INVALID_PROPERTY:       return OML_MSG_PLOT_INVALID_PROPERTY;
        case OML_ERR_PLOT_INVALID_OBJECT_HANDLE:  return OML_MSG_PLOT_INVALID_OBJECT_HANDLE;
        case OML_ERR_PLOT_INVALID_FIGURE_HANDLE:  return OML_MSG_PLOT_INVALID_FIGURE_HANDLE;
        case OML_ERR_PLOT_INVALID_AXES_HANDLE:    return OML_MSG_PLOT_INVALID_AXES_HANDLE;
        case OML_ERR_PLOT_INVALID_CURVE_HANDLE:   return OML_MSG_PLOT_INVALID_CURVE_HANDLE;
        case OML_ERR_PLOT_INVALID_COLOR:          return OML_MSG_PLOT_INVALID_COLOR;
        case OML_ERR_PLOT_NOT_SUPPORTED:          return OML_MSG_PLOT_NOT_SUPPORTED;
        case OML_ERR_PLOT_NOT_SUPPORTED_FOR_2D:   return OML_MSG_PLOT_NOT_SUPPORTED_FOR_2D;
        case OML_ERR_PLOT_UNKNOWN_ERROR:          return OML_MSG_PLOT_UNKNOWN_ERROR;
        case OML_ERR_PLOT_ZERORANGE:              return OML_MSG_PLOT_ZERORANGE;
        case OML_ERR_PLOT_ZIN2D:                  return OML_MSG_PLOT_ZIN2D;
        case OML_ERR_PLOT_UNSUPPORTED_FORMAT:     return OML_MSG_PLOT_UNSUPPORTED_FORMAT;
        case OML_ERR_PLOT_EMPTY_PROPERTY:         return OML_MSG_PLOT_EMPTY_PROPERTY;
        case OML_ERR_PLOT_MATXY_NOT_MATCH:        return OML_MSG_PLOT_MATXY_NOT_MATCH;
        case OML_ERR_PLOT_MATXZ_NOT_MATCH:        return OML_MSG_PLOT_MATXZ_NOT_MATCH;
        case OML_ERR_PLOT_MATYZ_NOT_MATCH:        return OML_MSG_PLOT_MATYZ_NOT_MATCH;
        case OML_ERR_PLOT_XZ_NOT_MATCH:           return OML_MSG_PLOT_XZ_NOT_MATCH;
        case OML_ERR_PLOT_YZ_NOT_MATCH:           return OML_MSG_PLOT_YZ_NOT_MATCH;
        case OML_ERR_PLOT_X_Z2_NOT_MATCH:         return OML_MSG_PLOT_X_Z2_NOT_MATCH;
        case OML_ERR_PLOT_Y_Z1_NOT_MATCH:         return OML_MSG_PLOT_Y_Z1_NOT_MATCH;
        case OML_ERR_PLOT_CONST_PROPERTY:         return OML_MSG_PLOT_CONST_PROPERTY;
        case OML_ERR_PLOT_CANNOT_OPEN_IMAGE:      return OML_MSG_PLOT_CANNOT_OPEN_IMAGE;
        case OML_ERR_PLOT_NEED_NORM_DATA:         return OML_MSG_PLOT_NEED_NORM_DATA;
        case OML_ERR_PLOT_NEED_PIXEL_DATA:        return OML_MSG_PLOT_NEED_PIXEL_DATA;
        case OML_ERR_PLOT_AMBIGUOUS_PROPERTY:     return OML_MSG_PLOT_AMBIGUOUS_PROPERTY;

	    // ABF error messages:
	    case OML_ERR_ABF_CREATE_FAILED:  return OML_MSG_ABF_CREATE_FAILED;
	    case OML_ERR_ABF_WRITE_FAILED:   return OML_MSG_ABF_WRITE_FAILED;
	    case OML_ERR_ABF_SUBCASE_EMPTY:  return OML_MSG_ABF_SUBCASE_EMPTY;
	    case OML_ERR_ABF_EXPORT_DONE:    return OML_MSG_ABF_EXPORT_DONE;

        // HDF5 reader/writer error messages:
        case OML_ERR_HDF5_INVALID_FILE:              return OML_MSG_HDF5_INVALID_FILE;
        case OML_ERR_HDF5_FILE_READ_FAILED:          return OML_MSG_HDF5_FILE_READ_FAILED;
        case OML_ERR_HDF5_INVALID_GROUP:             return OML_MSG_HDF5_INVALID_GROUP;
        case OML_ERR_HDF5_GROUP_READ_FAILED:         return OML_MSG_HDF5_GROUP_READ_FAILED;
        case OML_ERR_HDF5_INVALID_DATASET:           return OML_MSG_HDF5_INVALID_DATASET;
        case OML_ERR_HDF5_DATASET_READ_FAILED:       return OML_MSG_HDF5_DATASET_READ_FAILED;
        case OML_ERR_HDF5_NEITHER_DATASET_NOR_GROUP: return OML_MSG_HDF5_NEITHER_DATASET_NOR_GROUP;
        case OML_ERR_HDF5_ATTRIBUTE_READ_FAILED:     return OML_MSG_HDF5_ATTRIBUTE_READ_FAILED;
        case OML_ERR_HDF5_DATATYPE_READ_FAILED:      return OML_MSG_HDF5_DATATYPE_READ_FAILED;
        case OML_ERR_HDF5_DATASPACE_READ_FAILED:     return OML_MSG_HDF5_DATASPACE_READ_FAILED;
        case OML_ERR_HDF5_POINTS_MATRIXDIM_INVALID:  return OML_MSG_HDF5_POINTS_MATRIXDIM_INVALID;
        case OML_ERR_HDF5_POINTS_SELECTION_INVALID:  return OML_MSG_HDF5_POINTS_SELECTION_INVALID;
        case OML_ERR_HDF5_UNSUPPORTDIM:              return OML_MSG_HDF5_UNSUPPORTDIM;
        case OML_ERR_HDF5_UNSUPPORT_DATA:            return OML_MSG_HDF5_UNSUPPORT_DATA;
        case OML_ERR_HDF5_FILE_CREATION_FAILED:      return OML_MSG_HDF5_FILE_CREATION_FAILED;
        case OML_ERR_HDF5_GROUP_EXIST:               return OML_MSG_HDF5_GROUP_EXIST;
        case OML_ERR_HDF5_GROUP_CREATION_FAILED:     return OML_MSG_HDF5_GROUP_CREATION_FAILED;
        case OML_ERR_HDF5_GROUP_RENAME_FAILED:       return OML_MSG_HDF5_GROUP_RENAME_FAILED;
        case OML_ERR_HDF5_GROUP_REMOVE_FAILED:       return OML_MSG_HDF5_GROUP_REMOVE_FAILED;
        case OML_ERR_HDF5_DATASET_CREATION_FAILED:   return OML_MSG_HDF5_DATASET_CREATION_FAILED;
        case OML_ERR_HDF5_WRITE_FAILED:              return OML_MSG_HDF5_WRITE_FAILED;
        case OML_ERR_HDF5_DATASET_RENAME_FAILED:     return OML_MSG_HDF5_DATASET_RENAME_FAILED;
        case OML_ERR_HDF5_DATASET_REMOVE_FAILED:     return OML_MSG_HDF5_DATASET_REMOVE_FAILED;
        case OML_ERR_HDF5_DATASET_EXIST:             return OML_MSG_HDF5_DATASET_EXIST;
        case OML_ERR_HDF5_CHUNK_MATRIXDIM_INVALID:   return OML_MSG_HDF5_CHUNK_MATRIXDIM_INVALID;
        case OML_ERR_HDF5_CHUNK_MATRIX_INVALID_ROW:  return OML_MSG_HDF5_CHUNK_MATRIX_INVALID_ROW;
        case OML_ERR_HDF5_EMPTYMATRIX:               return OML_MSG_HDF5_EMPTYMATRIX;
        case OML_ERR_HDF5_EMPTYCELL:                 return OML_MSG_HDF5_EMPTYCELL;
        case OML_ERR_HDF5_INVALID_CELL_MEMBERS:      return OML_MSG_HDF5_INVALID_CELL_MEMBERS;
        case OML_ERR_HDF5_STRUCT_INVALID_DIMS:       return OML_MSG_HDF5_STRUCT_INVALID_DIMS;
        case OML_ERR_HDF5_ATTRIBUTE_CREATION_FAILED: return OML_MSG_HDF5_ATTRIBUTE_CREATION_FAILED;
        case OML_ERR_HDF5_INVALID_LOCATION:          return OML_MSG_HDF5_INVALID_LOCATION;
        case OML_ERR_HDF5_ATTRIBUTE_EXIST:           return OML_MSG_HDF5_ATTRIBUTE_EXIST;
        case OML_ERR_HDF5_INVALID_ATTRIBUTE:         return OML_MSG_HDF5_INVALID_ATTRIBUTE;
        case OML_ERR_HDF5_INVALID_STRUCT_MEMBERS:    return OML_MSG_HDF5_INVALID_STRUCT_MEMBERS;
        case OML_ERR_HDF5_INVALID_ATTRIBUTE_DATA:    return OML_MSG_HDF5_INVALID_ATTRIBUTE_DATA;

        //oml menu apis code
        case OML_ERR_MENUAPI_SYS_NAME:            return OML_MSG_MENUAPI_SYS_NAME;
        case OML_ERR_MENUAPI_SYS_DEL:             return OML_MSG_MENUAPI_SYS_DEL;
        case OML_ERR_MENUAPI_SYS_MOD:             return OML_MSG_MENUAPI_SYS_MOD;
        case OML_ERR_MENUAPI_INVALID_HANDLE:      return OML_MSG_MENUAPI_INVALID_HANDLE;
        case OML_ERR_MENUAPI_INVALID_HANDLE_TYPE: return OML_MSG_MENUAPI_INVALID_HANDLE_TYPE;

        case OML_ERR_IMAGE_ALPHA:     return OML_MSG_IMAGE_ALPHA;
        case OML_ERR_IMAGE_GRAYSCALE: return OML_MSG_IMAGE_GRAYSCALE;
        case OML_ERR_IMAGE_RGB:       return OML_MSG_IMAGE_RGB;
        case OML_ERR_IMAGE_HANDLE:    return OML_MSG_IMAGE_HANDLE;
        case OML_ERR_REALVECTOR2:     return OML_MSG_REALVECTOR2;
        case OML_ERR_REALVECTOR3:     return OML_MSG_REALVECTOR3;

        // RPC Writer
        case OML_ERR_RPC_INVALID_TIME_MAT:      return OML_MSG_RPC_INVALID_TIME_MAT;
        case OML_ERR_RPC_INVALID_DATA_MAT:      return OML_MSG_RPC_INVALID_DATA_MAT;
        case OML_ERR_RPC_INVALID_TIME_DATA_MAT: return OML_MSG_RPC_INVALID_TIME_DATA_MAT;
        case OML_ERR_RPC_INVALID_HEADER:        return OML_MSG_RPC_INVALID_HEADER;
        case OML_ERR_RPC_WRITE_FAILED:          return OML_MSG_RPC_WRITE_FAILED;
        case OML_ERR_RPC_INVALID_CURVENAME:     return OML_MSG_RPC_INVALID_CURVENAME;
        case OML_ERR_RPC_LOAD_FAILED:           return OML_MSG_RPC_LOAD_FAILED;
        case OML_ERR_RPC_XVECT_NUM_POINTS:      return OML_MSG_RPC_XVECT_NUM_POINTS;
        case OML_ERR_RPC_XVECT_SPREAD:          return OML_MSG_RPC_XVECT_SPREAD;

        // ReadCAE Builder && HW reader error messages:
        case OML_ERR_HWREADER_INVALID_TYPE:           return OML_MSG_HWREADER_INVALID_TYPE;
        case OML_ERR_HWREADER_INVALID_REQUEST:        return OML_MSG_HWREADER_INVALID_REQUEST;
        case OML_ERR_HWREADER_MISSING_COMPONENT:      return OML_MSG_HWREADER_MISSING_COMPONENT;
        case OML_ERR_HWREADER_INVALID_COMPONENT:      return OML_MSG_HWREADER_INVALID_COMPONENT;
        case OML_ERR_HWREADER_MISSING_TIME:           return OML_MSG_HWREADER_MISSING_TIME;
        case OML_ERR_HWREADER_FAIL_READER_INIT:       return OML_MSG_HWREADER_FAIL_READER_INIT;
        case OML_ERR_HWREADER_READING_FAILED:         return OML_MSG_HWREADER_READING_FAILED;
        case OML_ERR_HWREADER_INVALID_FILE_ENTRY:     return OML_MSG_HWREADER_INVALID_FILE_ENTRY;
        case OML_ERR_HWREADER_TIMECHANNELS_COMPARE:   return OML_MSG_HWREADER_TIMECHANNELS_COMPARE;
        case OML_ERR_HWREADER_SUBCASE_INVALID_RANGE:  return OML_MSG_HWREADER_SUBCASE_INVALID_RANGE;
        case OML_ERR_HWREADER_INVALID_TIME:           return OML_MSG_HWREADER_INVALID_TIME;
        case OML_ERR_HWREADER_NEED_SUBCASE:           return OML_MSG_HWREADER_NEED_SUBCASE;
        case OML_ERR_HWREADER_NO_SUBCASE:             return OML_MSG_HWREADER_NO_SUBCASE;
        case OML_ERR_HWREADER_STRING_POSINTEGER:      return OML_MSG_HWREADER_STRING_POSINTEGER;
        case OML_ERR_HWREADER_INVALID_TRC:            return OML_MSG_HWREADER_INVALID_TRC;
        case OML_ERR_HWREADER_SUBCASE_READ_FAIL:      return OML_MSG_HWREADER_SUBCASE_READ_FAIL;
        case OML_ERR_HWREADER_TYPE_READ_FAIL:         return OML_MSG_HWREADER_TYPE_READ_FAIL;
        case OML_ERR_HWREADER_REQUEST_READ_FAIL:      return OML_MSG_HWREADER_REQUEST_READ_FAIL;
        case OML_ERR_HWREADER_COMP_READ_FAIL:         return OML_MSG_HWREADER_COMP_READ_FAIL;
        case OML_ERR_HWREADER_FILTER_COMP_READ_FAIL:  return OML_MSG_HWREADER_FILTER_COMP_READ_FAIL;
        case OML_ERR_HWREADER_INVALID_TYPE_INDX:      return OML_MSG_HWREADER_INVALID_TYPE_INDX;
        case OML_ERR_HWREADER_INVALID_REQUEST_INDX:   return OML_MSG_HWREADER_INVALID_REQUEST_INDX;
        case OML_ERR_HWREADER_INVALID_COMPONENT_INDX: return OML_MSG_HWREADER_INVALID_COMPONENT_INDX;
        case OML_ERR_HWREADER_CELL_INDX:              return OML_MSG_HWREADER_CELL_INDX;
        case OML_ERR_HWREADER_INDX_POSINTEGER:        return OML_MSG_HWREADER_INDX_POSINTEGER;
        case OML_ERR_HWREADER_EXTRACT_STR_SCALAR:     return OML_MSG_HWREADER_EXTRACT_STR_SCALAR;
        case OML_ERR_HWREADER_DIMENSIONS_MATCH:       return OML_MSG_HWREADER_DIMENSIONS_MATCH;
        case OML_ERR_HWREADER_INDICES_REAL:           return OML_MSG_HWREADER_INDICES_REAL;
        case OML_ERR_HWREADER_TRC_DATA_TYPE:          return OML_MSG_HWREADER_TRC_DATA_TYPE;
        case OML_ERR_HWREADER_STRC_DATA_TYPE:         return OML_MSG_HWREADER_STRC_DATA_TYPE;
        case OML_ERR_HWREADER_CELL_4_5_COLS:          return OML_MSG_HWREADER_CELL_4_5_COLS;
        case OML_ERR_HWREADER_CELL_3_4_COLS:          return OML_MSG_HWREADER_CELL_3_4_COLS;
        case OML_ERR_HWREADER_REQ_CELL_1_N:           return OML_MSG_HWREADER_REQ_CELL_1_N;
        case OML_ERR_HWREADER_COMP_CELL_1_N:          return OML_MSG_HWREADER_COMP_CELL_1_N;
        case OML_ERR_HWREADER_TIME_CELL_1_N:          return OML_MSG_HWREADER_TIME_CELL_1_N;
        case OML_ERR_HWREADER_OUTPUT_CELL_MAT:        return OML_MSG_HWREADER_OUTPUT_CELL_MAT;
        case OML_ERR_HWREADER_INVALID_UNIT_TYPE:      return OML_MSG_HWREADER_INVALID_UNIT_TYPE;

        case OML_ERR_SPREADSHEET_RANGE:  return OML_MSG_SPREADSHEET_RANGE;
        case OML_ERR_WINDOWSONLY_OP:     return OML_MSG_WINDOWSONLY_OP;
        case OML_ERR_FILEALREADYOPEN_OP: return OML_MSG_FILEALREADYOPEN_OP;
        case OML_ERR_OMLINTERFACE_OP:    return OML_MSG_OMLINTERFACE_OP;

        case OML_ERR_MQTT_INVALID_OPTION:             return OML_MSG_MQTT_INVALID_OPTION;
        case OML_ERR_MQTT_CLIENTID_INUSE:             return OML_MSG_MQTT_CLIENTID_INUSE;
        case OML_ERR_MQTT_CLIENTID_INVALID:           return OML_MSG_MQTT_CLIENTID_INVALID;
        case OML_ERR_MQTT_CLIENTCREATION_FAIL:        return OML_MSG_MQTT_CLIENTCREATION_FAIL;
        case OML_ERR_MQTT_CLIENTID_SIZE:              return OML_MSG_MQTT_CLIENTID_SIZE;
        case OML_ERR_MQTT_MIN_IDLETIME:               return OML_MSG_MQTT_MIN_IDLETIME;

        default: break;
    }

    return std::string();
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
    case OML_VAR_RANDSEED:     varStr = OML_STR_RANDSEED;     break;
    case OML_VAR_TOLX:         varStr = OML_STR_TOLX;         break;
    case OML_VAR_TOLFUN:       varStr = OML_STR_TOLFUN;       break;
    case OML_VAR_TOLFUNABS:    varStr = OML_STR_TOLFUNABS;    break;
    case OML_VAR_TOLFUNREL:    varStr = OML_STR_TOLFUNREL;    break;
    case OML_VAR_TOLCON:       varStr = OML_STR_TOLCON;       break;
    case OML_VAR_CONRET:       varStr = OML_STR_CONRET;       break;
    case OML_VAR_MOVELIM:      varStr = OML_STR_MOVELIM;      break;
    case OML_VAR_PERTM:        varStr = OML_STR_PERTM;        break;
    case OML_VAR_PERTV:        varStr = OML_STR_PERTV;        break;
    case OML_VAR_POPSIZE:      varStr = OML_STR_POPSIZE;      break;
    case OML_VAR_INITPOP:      varStr = OML_STR_INITPOP;      break;
    case OML_VAR_CRODIST:      varStr = OML_STR_CRODIST;      break;
    case OML_VAR_INITSAMPNTS:  varStr = OML_STR_INITSAMPNTS;  break;
    case OML_VAR_MAXFAIL:      varStr = OML_STR_MAXFAIL;      break;
    case OML_VAR_METHOD:       varStr = OML_STR_METHOD;       break;
    case OML_VAR_PNTSPERITER:  varStr = OML_STR_PNTSPERITER;  break;
    case OML_VAR_STOPNOIMPR :  varStr = OML_STR_STOPNOIMPR;   break;
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
