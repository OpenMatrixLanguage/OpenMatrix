/**
* @file BuiltInFuncs.cpp
* @date October 2013
* Copyright (C) 2013-2018 Altair Engineering, Inc.  
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

// Begin defines/includes
#if OS_WIN
#	if _MSC_VER == 1900
#define _CRT_NO_VA_START_VALIDATION
#   endif
#endif

#ifdef OS_WIN
#    define NOMINMAX
#endif

#include "BuiltInFuncs.h"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <memory>

#include "Currency.h"
#include "BuiltInFuncsConvert.h"
#include "BuiltInFuncsCore.h"
#include "BuiltInFuncsData.h"
#include "BuiltInFuncsElemMath.h"
#include "BuiltInFuncsFile.h"
#include "BuiltInFuncsString.h"
#include "BuiltInFuncsSystem.h"
#include "BuiltInFuncsUtils.h"
#include "MemoryScope.h"
#include "Interpreter.h"
#include "FunctionInfo.h"
#include "ErrorInfo.h"
#include "OML_Error.h"
#include "hwComplex.h"
#include "StructData.h"
#include "MatrixNUtils.h"

#include <cmath>
#include <time.h>
#include <stdio.h>
#include <cctype>
#include <cerrno>
#include <cstdint>
#include <regex>
#include <deque>
#include <functional>
#include <sys/stat.h>
#include "ExprCppTreeParser.h"
#include <fstream>
#include <cfloat>
#ifdef OS_WIN
#include <Windows.h>
#include <WinBase.h>
#define S_ISDIR(mode) (mode & S_IFDIR)
#else
#include <unistd.h>
#include <glob.h>
#include <libgen.h>
#include <dirent.h>
#include <sys/resource.h>
#include <sys/times.h>
#include <time.h>
#include <dlfcn.h>
#endif

#define hconcat concat<true>
#define vconcat concat<false>
#define HW_MAKE_BOOL_CURRENCY(val) val ? getTrue() : getFalse()
#define HW_VALID_FORMATCHARS "diouxXfeEgGcs%"
#define HW_BAD_VARNAME_CHARS "\"'!#%^&*()-+=./`~<>[]{},:;\\|\? \0"
#define HW_TOLERANCE 1e-15
#define SEEK_OFFSET -1

#undef GetMessage

static clock_t tictime = -1;

enum DataType
{
    Int,
    Short,
    Long,
    Char,
    Float,
    Double,
    LongLong,
    Int8,
    Int16,
    Int32,
    Int64
};

struct Precision
{
    bool sign;
    size_t numBytes;
    int blockSize;
    DataType dtype;
    Precision (bool issigned, size_t numOfBytes, int mult, DataType type) : 
              sign(issigned), numBytes(numOfBytes), blockSize(mult), dtype(type) {} 
};

enum FormatType
{
    Percent,
    String,
    Truncate,
    NoTruncate
};

enum FileStreamType
{
    Stdout,
    Stderr,
    Other
};

#ifdef OS_WIN
static char pathsep = ';';
#else
static char pathsep = ':';
#endif

#define CORE   "CoreMinimalInterpreter"
#define DATA   "DataStructures"
#define DIFF   "DifferentialEquations"
#define ELEM   "ElementaryMath"
#define FILEIO "FileIO"
#define READER "HWReaders"
#define LINA   "LinearAlgebra"
#define LOGI   "Logical"
#define OPTI   "Optimization"
#define PLOT   "Plotting"
#define POLY   "PolynomialMath"
#define SIGN   "SignalProcessing"
#define STATAN "StatisticalAnalysis"
#define STNG   "StringOperations"
#define SYSTEM "System"
#define TIME   "Time"
#define TRIG   "Trigonometry"
#define GUI    "Gui"

// End defines/includes

//------------------------------------------------------------------------------
// Dummy built in function
//------------------------------------------------------------------------------
bool DummyVoid(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    return true;
}
//------------------------------------------------------------------------------
// Maps the built in functions
//------------------------------------------------------------------------------
void mapBuiltInFuncs(std::map<std::string, BuiltinFunc>* std_functions)
{
    (*std_functions)["nargout"]            = BuiltinFunc(oml_nargout, FunctionMetaData(1, 1, CORE));
    (*std_functions)["nargin"]             = BuiltinFunc(oml_nargin, FunctionMetaData(1, 1, CORE));
    (*std_functions)["bitand"]             = BuiltinFunc(oml_bitand, FunctionMetaData(2, 1, LOGI));
    (*std_functions)["bitor"]              = BuiltinFunc(oml_bitor, FunctionMetaData(2, 1, LOGI));
    (*std_functions)["bitxor"]             = BuiltinFunc(oml_bitxor, FunctionMetaData(2, 1, LOGI));
    (*std_functions)["isdir"]              = BuiltinFunc(oml_isdir, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["rmfield"]            = BuiltinFunc(oml_rmfield, FunctionMetaData(2, 1, DATA));
    (*std_functions)["ispc"]               = BuiltinFunc(oml_ispc, FunctionMetaData(0, 1, FILEIO));
    (*std_functions)["ismac"]              = BuiltinFunc(oml_ismac, FunctionMetaData(0, 1, FILEIO));
    (*std_functions)["isunix"]             = BuiltinFunc(oml_isunix, FunctionMetaData(0, 1, FILEIO));
    (*std_functions)["assignin"]           = BuiltinFunc(oml_assignin, FunctionMetaData(3, 0, CORE));
    (*std_functions)["evalin"]             = BuiltinFunc(oml_evalin, FunctionMetaData(2, -1, CORE));
    (*std_functions)["struct"]             = BuiltinFunc(oml_struct, FunctionMetaData(-1, 1, DATA));
    (*std_functions)["iscell"]             = BuiltinFunc(oml_iscell, FunctionMetaData(1, 1, DATA));
    (*std_functions)["isfinite"]           = BuiltinFunc(oml_isfinite, FunctionMetaData(1, 1, DATA));
    (*std_functions)["and"]                = BuiltinFunc(oml_and, FunctionMetaData(-3, 1, CORE));
    (*std_functions)["or"]                 = BuiltinFunc(oml_or, FunctionMetaData(-3, 1, CORE));
    (*std_functions)["uplus"]              = BuiltinFunc(oml_uplus, FunctionMetaData(1, 1, CORE));
    (*std_functions)["not"]                = BuiltinFunc(oml_not, FunctionMetaData(1, 1, CORE));
    (*std_functions)["uminus"]             = BuiltinFunc(oml_uminus, FunctionMetaData(1, 1, CORE));
    (*std_functions)["rdivide"]            = BuiltinFunc(oml_rdivide, FunctionMetaData(2, 1, CORE));
    (*std_functions)["ldivide"]            = BuiltinFunc(oml_ldivide, FunctionMetaData(2, 1, CORE));
    (*std_functions)["mrdivide"]           = BuiltinFunc(oml_mrdivide, FunctionMetaData(2, 1, CORE));
    (*std_functions)["mldivide"]           = BuiltinFunc(oml_mldivide, FunctionMetaData(2, 1, CORE));
    (*std_functions)["plus"]               = BuiltinFunc(oml_plus, FunctionMetaData(-3, 1, CORE));
    (*std_functions)["times"]              = BuiltinFunc(oml_times, FunctionMetaData(-3, 1, CORE));
    (*std_functions)["mtimes"]             = BuiltinFunc(oml_mtimes, FunctionMetaData(-3, 1, CORE));
    (*std_functions)["minus"]              = BuiltinFunc(oml_minus, FunctionMetaData(2, 1, CORE));
    (*std_functions)["power"]              = BuiltinFunc(oml_power, FunctionMetaData(2, 1, CORE));
    (*std_functions)["mpower"]             = BuiltinFunc(oml_mpower, FunctionMetaData(2, 1, CORE));
    (*std_functions)["rmdir"]              = BuiltinFunc(oml_rmdir, FunctionMetaData(1, 3, SYSTEM));
    (*std_functions)["date"]               = BuiltinFunc(oml_date, FunctionMetaData(0, 1, SYSTEM));
    (*std_functions)["restoredefaultpath"] = BuiltinFunc(oml_restorepath, FunctionMetaData(0, 1, CORE));
    (*std_functions)["lt"]                 = BuiltinFunc(oml_lt, FunctionMetaData(2, 1, CORE));
    (*std_functions)["gt"]                 = BuiltinFunc(oml_gt, FunctionMetaData(2, 1, CORE));
    (*std_functions)["ge"]                 = BuiltinFunc(oml_geq, FunctionMetaData(2, 1, CORE));
    (*std_functions)["le"]                 = BuiltinFunc(oml_leq, FunctionMetaData(2, 1, CORE));
    (*std_functions)["eq"]                 = BuiltinFunc(oml_eq, FunctionMetaData(2, 1, CORE));
    (*std_functions)["ne"]                 = BuiltinFunc(oml_ne, FunctionMetaData(2, 1, CORE));
    (*std_functions)["linsolve"]           = BuiltinFunc(oml_linsolve, FunctionMetaData(3, 1, LINA));
    (*std_functions)["pinv"]               = BuiltinFunc(oml_pinv, FunctionMetaData(1, 1, LINA));
    (*std_functions)["normalize"]          = BuiltinFunc(oml_normalize, FunctionMetaData(1, 1, CORE));
    (*std_functions)["expm"]               = BuiltinFunc(oml_expm, FunctionMetaData(1, 1, LINA));
    (*std_functions)["mkdir"]              = BuiltinFunc(oml_mkdir, FunctionMetaData(1, 3, SYSTEM));
    (*std_functions)["cond"]               = BuiltinFunc(oml_cond, FunctionMetaData(1, 1, LINA));
    (*std_functions)["sort"]               = BuiltinFunc(oml_sort, FunctionMetaData(3, 2, CORE));
    (*std_functions)["isnumeric"]          = BuiltinFunc(oml_isnumeric, FunctionMetaData(1, 1, CORE));
    (*std_functions)["all"]                = BuiltinFunc(oml_all, FunctionMetaData(2, 1, CORE));
    (*std_functions)["any"]                = BuiltinFunc(oml_any, FunctionMetaData(2, 1, CORE));
    (*std_functions)["error"]              = BuiltinFunc(oml_error, FunctionMetaData(-2, 0, CORE));
    (*std_functions)["warning"]            = BuiltinFunc(oml_warning, FunctionMetaData(-2, 0, CORE));
    (*std_functions)["isinf"]              = BuiltinFunc(oml_isinf, FunctionMetaData(1, 1, CORE));
    (*std_functions)["isnan"]              = BuiltinFunc(oml_isnan, FunctionMetaData(1, 1, CORE));
    (*std_functions)["isscalar"]           = BuiltinFunc(oml_isscalar, FunctionMetaData(1, 1, CORE));
    (*std_functions)["isstr"]              = BuiltinFunc(oml_isstr, FunctionMetaData(1, 1, STNG));
    (*std_functions)["ischar"]             = BuiltinFunc(oml_isstr, FunctionMetaData(1, 1, STNG));
    (*std_functions)["ismatrix"]           = BuiltinFunc(oml_ismatrix, FunctionMetaData(1, 1, CORE));
    (*std_functions)["issymmetric"]        = BuiltinFunc(oml_issymmetric, FunctionMetaData(-2, 1, LINA));
	(*std_functions)["issquare"]           = BuiltinFunc(oml_issquare, FunctionMetaData(1, 1, LINA));
    (*std_functions)["ishermitian"]        = BuiltinFunc(oml_ishermitian, FunctionMetaData(-2, 1, LINA));
    (*std_functions)["genvarname"]         = BuiltinFunc(oml_genvarname, FunctionMetaData(2, 1, CORE));
    (*std_functions)["sub2ind"]            = BuiltinFunc(oml_sub2ind, FunctionMetaData(-2, 1, ELEM));
    (*std_functions)["setdiff"]            = BuiltinFunc(oml_setdiff, FunctionMetaData(-1, -1, ELEM));
    (*std_functions)["intersect"]          = BuiltinFunc(oml_intersect, FunctionMetaData(-1, -1, ELEM));
    (*std_functions)["setxor"]             = BuiltinFunc(oml_setxor, FunctionMetaData(-1, -1, ELEM));
    (*std_functions)["double"]             = BuiltinFunc(oml_double, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["regexp"]             = BuiltinFunc(oml_regexp, FunctionMetaData(-1, -1, STNG));
    (*std_functions)["cell2struct"]        = BuiltinFunc(oml_cell2struct, FunctionMetaData(3, 1, DATA));
    (*std_functions)["cplxpair"]           = BuiltinFunc(oml_cplxpair, FunctionMetaData(3, 1, ELEM));
    (*std_functions)["ind2sub"]            = BuiltinFunc(oml_ind2sub, FunctionMetaData(2, -1, ELEM));
    (*std_functions)["isprime"]            = BuiltinFunc(oml_isprime, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["conv2"]              = BuiltinFunc(oml_conv2, FunctionMetaData(4, 1, LINA));
    (*std_functions)["vertcat"]            = BuiltinFunc(oml_vertcat, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["isletter"]           = BuiltinFunc(oml_isletter, FunctionMetaData(1, 1, STNG));
    (*std_functions)["isspace"]            = BuiltinFunc(oml_isspace, FunctionMetaData(1, 1, STNG));
    (*std_functions)["horzcat"]            = BuiltinFunc(oml_horzcat, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["celldisp"]           = BuiltinFunc(oml_celldisp, FunctionMetaData(2, 0, DATA));
    (*std_functions)["gcd"]                = BuiltinFunc(oml_gcd, FunctionMetaData(-1, -2, ELEM));
    (*std_functions)["pathsep"]            = BuiltinFunc(oml_pathsep, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["func2str"]           = BuiltinFunc(oml_func2str, FunctionMetaData(1, 1, STNG));
    (*std_functions)["str2func"]           = BuiltinFunc(oml_str2func, FunctionMetaData(1, 1, STNG));
    (*std_functions)["strncmpi"]           = BuiltinFunc(oml_strncmpi, FunctionMetaData(3, 1, STNG));
    (*std_functions)["strcmpi"]            = BuiltinFunc(oml_strcmpi, FunctionMetaData(2, 1, STNG));
    (*std_functions)["cellfun"]            = BuiltinFunc(oml_cellfun, FunctionMetaData(-1, -1, DATA));
    (*std_functions)["datenum"]            = BuiltinFunc(oml_datenum, FunctionMetaData(6, 1, TIME));
    (*std_functions)["issorted"]           = BuiltinFunc(oml_issorted, FunctionMetaData(3, 1, ELEM));
    (*std_functions)["getfield"]           = BuiltinFunc(oml_getfield, FunctionMetaData(-2, -1, DATA));
    (*std_functions)["ismember"]           = BuiltinFunc(oml_ismember, FunctionMetaData(3, 2, ELEM));
    (*std_functions)["arg"]                = BuiltinFunc(oml_angle, FunctionMetaData(1, 1, CORE));
    (*std_functions)["angle"]              = BuiltinFunc(oml_angle, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["isfield"]            = BuiltinFunc(oml_isfield, FunctionMetaData(2, 1, DATA));
    (*std_functions)["which"]              = BuiltinFunc(oml_which, FunctionMetaData(-1, -1, CORE));
    (*std_functions)["chdir"]              = BuiltinFunc(oml_cd, FunctionMetaData(1, 1, SYSTEM));
    (*std_functions)["cd"]                 = BuiltinFunc(oml_cd, FunctionMetaData(1, 1, SYSTEM));
    (*std_functions)["linspace"]           = BuiltinFunc(oml_linspace, FunctionMetaData(3, 1, CORE));
    (*std_functions)["isglobal"]           = BuiltinFunc(oml_isglobal, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["printf"]             = BuiltinFunc(oml_printf, FunctionMetaData(-2, 1, FILEIO));
    (*std_functions)["fprintf"]            = BuiltinFunc(oml_fprintf, FunctionMetaData(-3, 1, FILEIO));
    (*std_functions)["str2num"]            = BuiltinFunc(oml_str2num, FunctionMetaData(1, 2, STNG));
    (*std_functions)["str2double"]         = BuiltinFunc(oml_str2double, FunctionMetaData(1, 2, STNG));
    (*std_functions)["fieldnames"]         = BuiltinFunc(oml_fieldnames, FunctionMetaData(1, 1, DATA));
    (*std_functions)["sprintf"]            = BuiltinFunc(oml_sprintf, FunctionMetaData(-2, 2, STNG));
    (*std_functions)["feval"]              = BuiltinFunc(oml_feval, FunctionMetaData(-2, -1, CORE));
    (*std_functions)["rmpath"]             = BuiltinFunc(oml_rmpath, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["addpath"]            = BuiltinFunc(oml_addpath, FunctionMetaData(-1, 1, CORE));
	(*std_functions)["registerpath"]           = BuiltinFunc(oml_addpath2, FunctionMetaData(2, 1, CORE));
    (*std_functions)["path"]               = BuiltinFunc(oml_path, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["repmat"]             = BuiltinFunc(oml_repmat, FunctionMetaData(-2, 1, ELEM));
    (*std_functions)["asind"]              = BuiltinFunc(oml_asind, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["deg2rad"]            = BuiltinFunc(oml_deg2rad, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["rad2deg"]            = BuiltinFunc(oml_rad2deg, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["acosd"]              = BuiltinFunc(oml_acosd, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["cosd"]               = BuiltinFunc(oml_cosd, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["sind"]               = BuiltinFunc(oml_sind, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["cross"]              = BuiltinFunc(oml_cross, FunctionMetaData(3, 1, LINA));
    (*std_functions)["isstruct"]           = BuiltinFunc(oml_isstruct, FunctionMetaData(1, 1, DATA));
    (*std_functions)["cell2mat"]           = BuiltinFunc(oml_cell2mat, FunctionMetaData(1, 1, DATA));
    (*std_functions)["fgetl"]              = BuiltinFunc(oml_fgetl, FunctionMetaData(2, 1, FILEIO));
    (*std_functions)["realmax"]            = BuiltinFunc(oml_realmax, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["realmin"]            = BuiltinFunc(oml_realmin, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["run"]                = BuiltinFunc(oml_run, FunctionMetaData(1, 0, CORE));
	(*std_functions)["run2"]               = BuiltinFunc(oml_run2, FunctionMetaData(1, 0, CORE));
    (*std_functions)["fflush"]             = BuiltinFunc(oml_fflush, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["frewind"]            = BuiltinFunc(oml_frewind, FunctionMetaData(1, 0, FILEIO));
    (*std_functions)["fwrite"]             = BuiltinFunc(oml_fwrite, FunctionMetaData(5, 1, FILEIO));
    (*std_functions)["nan"]                = BuiltinFunc(oml_nan, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["NaN"]                = BuiltinFunc(oml_nan, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["Inf"]                = BuiltinFunc(oml_inf, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["inf"]                = BuiltinFunc(oml_inf, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["fread"]              = BuiltinFunc(oml_fread, FunctionMetaData(5, 2, FILEIO));
    (*std_functions)["stdin"]              = BuiltinFunc(oml_stdin, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["stdout"]             = BuiltinFunc(oml_stdout, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["stderr"]             = BuiltinFunc(oml_stderr, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["SEEK_END"]           = BuiltinFunc(oml_seek_end, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["SEEK_CUR"]           = BuiltinFunc(oml_seek_cur, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["SEEK_SET"]           = BuiltinFunc(oml_seek_set, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["fseek"]              = BuiltinFunc(oml_fseek, FunctionMetaData(3, 1, FILEIO));
    (*std_functions)["feof"]               = BuiltinFunc(oml_feof, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["isbool"]             = BuiltinFunc(oml_islogical, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["islogical"]          = BuiltinFunc(oml_islogical, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["ftell"]              = BuiltinFunc(oml_ftell, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["fgets"]              = BuiltinFunc(oml_fgets, FunctionMetaData(2, 1, FILEIO));
    (*std_functions)["fclose"]             = BuiltinFunc(oml_fclose, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["fopen"]              = BuiltinFunc(oml_fopen, FunctionMetaData(3, 2, FILEIO));
    (*std_functions)["char"]               = BuiltinFunc(oml_char, FunctionMetaData(-1, 1, STNG));
    (*std_functions)["strsplit"]           = BuiltinFunc(oml_strsplit, FunctionMetaData(-3, 2, STNG));
    (*std_functions)["strjoin"]            = BuiltinFunc(oml_strjoin, FunctionMetaData(2, 1, STNG));
    (*std_functions)["strrep"]             = BuiltinFunc(oml_strrep, FunctionMetaData(5, 1, STNG));
    (*std_functions)["strncmp"]            = BuiltinFunc(oml_strncmp, FunctionMetaData(3, 1, STNG));
    (*std_functions)["strcmp"]             = BuiltinFunc(oml_strcmp, FunctionMetaData(2, 1, STNG));
    (*std_functions)["strcat"]             = BuiltinFunc(oml_strcat, FunctionMetaData(-1, 1, STNG));
    (*std_functions)["tolower"]            = BuiltinFunc(oml_lower, FunctionMetaData(1, 1, STNG));
    (*std_functions)["lower"]              = BuiltinFunc(oml_lower, FunctionMetaData(1, 1, STNG));
    (*std_functions)["toupper"]            = BuiltinFunc(oml_upper, FunctionMetaData(1, 1, STNG));
    (*std_functions)["upper"]              = BuiltinFunc(oml_upper, FunctionMetaData(1, 1, STNG));
    (*std_functions)["strfind"]            = BuiltinFunc(oml_strfind, FunctionMetaData(2, 1, STNG));
    (*std_functions)["isvarname"]          = BuiltinFunc(oml_isvarname, FunctionMetaData(1, 1, CORE));
    (*std_functions)["unique"]             = BuiltinFunc(oml_unique, FunctionMetaData(-2, -1, ELEM));
    (*std_functions)["cell"]               = BuiltinFunc(oml_cell, FunctionMetaData(-1, 1, DATA));
    (*std_functions)["toc"]                = BuiltinFunc(oml_toc, FunctionMetaData(1, 1, TIME));
    (*std_functions)["tic"]                = BuiltinFunc(oml_tic, FunctionMetaData(0, 1, TIME));
    (*std_functions)["factor"]             = BuiltinFunc(oml_factor, FunctionMetaData(1, 2, ELEM));
    (*std_functions)["polyval"]            = BuiltinFunc(oml_polyval, FunctionMetaData(4, 1, POLY));
    (*std_functions)["triu"]               = BuiltinFunc(oml_triu, FunctionMetaData(2, 1, LINA));
    (*std_functions)["tril"]               = BuiltinFunc(oml_tril, FunctionMetaData(2, 1, LINA));
    (*std_functions)["iscellstr"]          = BuiltinFunc(oml_iscellstr, FunctionMetaData(1, 1, DATA));   
    (*std_functions)["poly"]               = BuiltinFunc(oml_poly, FunctionMetaData(1, 1, POLY));
    (*std_functions)["fix"]                = BuiltinFunc(oml_fix, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["rem"]                = BuiltinFunc(oml_rem, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["reshape"]            = BuiltinFunc(oml_reshape, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["permute"]            = BuiltinFunc(oml_permute, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["squeeze"]            = BuiltinFunc(oml_squeeze, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["norm"]               = BuiltinFunc(oml_norm, FunctionMetaData(2, 1, LINA));
    (*std_functions)["J"]                  = BuiltinFunc(oml_i, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["I"]                  = BuiltinFunc(oml_i, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["j"]                  = BuiltinFunc(oml_i, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["i"]                  = BuiltinFunc(oml_i, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["find"]               = BuiltinFunc(oml_find, FunctionMetaData(3, 3, ELEM));
    (*std_functions)["eval"]               = BuiltinFunc(oml_eval, FunctionMetaData(2, -1, CORE));   
    (*std_functions)["length"]             = BuiltinFunc(oml_length, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["ndims"]              = BuiltinFunc(oml_ndims, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["size"]               = BuiltinFunc(oml_size, FunctionMetaData(2, -1, ELEM));
    (*std_functions)["isreal"]             = BuiltinFunc(oml_isreal, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["iscomplex"]          = BuiltinFunc(oml_iscomplex, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["isequal"]            = BuiltinFunc(oml_isequal, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["isempty"]            = BuiltinFunc(oml_isempty, FunctionMetaData(1, 1, ELEM));    
    (*std_functions)["clock"]              = BuiltinFunc(oml_clock, FunctionMetaData(0, 1, TIME));
    (*std_functions)["pwd"]                = BuiltinFunc(oml_pwd, FunctionMetaData(0, 1, FILEIO));
    (*std_functions)["balance"]            = BuiltinFunc(oml_balance, FunctionMetaData(3, 4, LINA));
    (*std_functions)["union"]              = BuiltinFunc(oml_union, FunctionMetaData(3, 3, ELEM));
    (*std_functions)["rank"]               = BuiltinFunc(oml_rank, FunctionMetaData(-2, 1, LINA));
    (*std_functions)["primes"]             = BuiltinFunc(oml_primes, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["hypot"]              = BuiltinFunc(oml_hypot, FunctionMetaData(2, 1, TRIG));
    (*std_functions)["qr"]                 = BuiltinFunc(oml_qr, FunctionMetaData(1, 2, LINA));
    (*std_functions)["lu"]                 = BuiltinFunc(oml_lu, FunctionMetaData(1, 3, LINA));
    (*std_functions)["schur"]              = BuiltinFunc(oml_schur, FunctionMetaData(1, 2, LINA));
    (*std_functions)["svd"]                = BuiltinFunc(oml_svd, FunctionMetaData(2, 3, LINA));
    (*std_functions)["eps"]                = BuiltinFunc(oml_eps, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["eig"]                = BuiltinFunc(oml_eig, FunctionMetaData(2, 2, LINA));
    (*std_functions)["chol"]               = BuiltinFunc(oml_chol, FunctionMetaData(2, 2, LINA));
    (*std_functions)["complex"]            = BuiltinFunc(oml_complex, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["transpose"]          = BuiltinFunc(oml_transpose, FunctionMetaData(1, 1, LINA));
    (*std_functions)["eye"]                = BuiltinFunc(oml_eye, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["sign"]               = BuiltinFunc(oml_sign, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["pow2"]               = BuiltinFunc(oml_pow2, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["dot"]                = BuiltinFunc(oml_dot, FunctionMetaData(3, 1, LINA));
    (*std_functions)["kron"]               = BuiltinFunc(oml_kron, FunctionMetaData(2, 1, LINA));
    (*std_functions)["diff"]               = BuiltinFunc(oml_diff, FunctionMetaData(3, 1, ELEM));
    (*std_functions)["diag"]               = BuiltinFunc(oml_diag, FunctionMetaData(-2, 1, ELEM));
    (*std_functions)["conv"]               = BuiltinFunc(oml_conv, FunctionMetaData(3, 1, LINA));
    (*std_functions)["e"]                  = BuiltinFunc(oml_e, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["pi"]                 = BuiltinFunc(oml_pi, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["ones"]               = BuiltinFunc(oml_ones, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["zeros"]              = BuiltinFunc(oml_zeros, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["trace"]              = BuiltinFunc(oml_trace, FunctionMetaData(1, 2, LINA));
    (*std_functions)["det"]                = BuiltinFunc(oml_det, FunctionMetaData(1, 2, LINA));
    (*std_functions)["rcond"]              = BuiltinFunc(oml_rcond, FunctionMetaData(1, 1, LINA));
    (*std_functions)["mod"]                = BuiltinFunc(oml_mod, FunctionMetaData(2, 2, ELEM));
    (*std_functions)["exp"]                = BuiltinFunc(oml_exp, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["real"]               = BuiltinFunc(oml_real, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["imag"]               = BuiltinFunc(oml_imag, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["conj"]               = BuiltinFunc(oml_conj, FunctionMetaData(1, 1, LINA));
    (*std_functions)["sum"]                = BuiltinFunc(oml_sum, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["cumsum"]             = BuiltinFunc(oml_cumsum, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["prod"]               = BuiltinFunc(oml_prod, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["cumprod"]            = BuiltinFunc(oml_cumprod, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["log10"]              = BuiltinFunc(oml_log10, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["log2"]               = BuiltinFunc(oml_log2, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["log"]                = BuiltinFunc(oml_log, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["abs"]                = BuiltinFunc(oml_abs, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["cosh"]               = BuiltinFunc(oml_cosh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["acosh"]              = BuiltinFunc(oml_acosh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["acos"]               = BuiltinFunc(oml_acos, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["cos"]                = BuiltinFunc(oml_cos, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["sin"]                = BuiltinFunc(oml_sin, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["sinh"]               = BuiltinFunc(oml_sinh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["asinh"]              = BuiltinFunc(oml_asinh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["asin"]               = BuiltinFunc(oml_asin, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["tan"]                = BuiltinFunc(oml_tan, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["tanh"]               = BuiltinFunc(oml_tanh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["atanh"]              = BuiltinFunc(oml_atanh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["atan"]               = BuiltinFunc(oml_atan, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["sqrt"]               = BuiltinFunc(oml_sqrt, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["disp"]               = BuiltinFunc(oml_print, FunctionMetaData(1, 0, CORE));
    (*std_functions)["ceil"]               = BuiltinFunc(oml_ceil, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["floor"]              = BuiltinFunc(oml_floor, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["round"]              = BuiltinFunc(oml_round, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["max"]                = BuiltinFunc(oml_max, FunctionMetaData(3, 2, ELEM));
    (*std_functions)["min"]                = BuiltinFunc(oml_min, FunctionMetaData(3, 2, ELEM));
    (*std_functions)["inv"]                = BuiltinFunc(oml_inv, FunctionMetaData(1, 2, LINA));
    (*std_functions)["logical"]            = BuiltinFunc(oml_logical, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["clear"]              = BuiltinFunc(oml_clear, FunctionMetaData(-1, 0, CORE));
    (*std_functions)["end"]                = BuiltinFunc(oml_end, FunctionMetaData(OML_NO_NARG, OML_NO_NARG, CORE));
    (*std_functions)["nargout"]            = BuiltinFunc(oml_nargout, FunctionMetaData(1, 1, CORE));
    (*std_functions)["nargin"]             = BuiltinFunc(oml_nargin, FunctionMetaData(1, 1, CORE));
    (*std_functions)["who"]                = BuiltinFunc(oml_who, FunctionMetaData(-1, -1, CORE));
    (*std_functions)["true"]               = BuiltinFunc(oml_true, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["false"]              = BuiltinFunc(oml_false, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["refcnt"]             = BuiltinFunc(oml_refcnt,FunctionMetaData( OML_NO_NARG, OML_NO_NARG, CORE));
    (*std_functions)["clc"]                = BuiltinFunc(DummyVoid, FunctionMetaData(0, 0, CORE));
    (*std_functions)["lasterr"]            = BuiltinFunc(oml_lasterr, FunctionMetaData(1, 1, CORE));
    (*std_functions)["lastwarn"]           = BuiltinFunc(oml_lastwarn, FunctionMetaData(1, 1, CORE));
    (*std_functions)["continue"]           = BuiltinFunc(oml_continue, FunctionMetaData(OML_NO_NARG, OML_NO_NARG, CORE));
    (*std_functions)["break"]              = BuiltinFunc(oml_break, FunctionMetaData(OML_NO_NARG, OML_NO_NARG, CORE));
    (*std_functions)["pol2cart"]           = BuiltinFunc(oml_pol2cart, FunctionMetaData(3, -1, ELEM));
    (*std_functions)["deblank"]            = BuiltinFunc(oml_deblank, FunctionMetaData(1, 1, STNG));
    (*std_functions)["isvector"]           = BuiltinFunc(oml_isvector, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["cart2pol"]           = BuiltinFunc(oml_cart2pol, FunctionMetaData(3, -1, ELEM));
    (*std_functions)["sph2cart"]           = BuiltinFunc(oml_sph2cart, FunctionMetaData(3, -1, ELEM));
    (*std_functions)["cart2sph"]           = BuiltinFunc(oml_cart2sph, FunctionMetaData(3, -1, ELEM));
    (*std_functions)["struct2cell"]        = BuiltinFunc(oml_struct2cell, FunctionMetaData(1, 1, DATA));
    (*std_functions)["fullfile"]           = BuiltinFunc(oml_fullfile, FunctionMetaData(-2, 1, FILEIO));
    (*std_functions)["atan2"]              = BuiltinFunc(oml_atan2, FunctionMetaData(2, 1, TRIG));
    (*std_functions)["builtin"]            = BuiltinFunc(oml_builtin, FunctionMetaData(-2, -1, CORE));
    (*std_functions)["assert"]             = BuiltinFunc(oml_assert, FunctionMetaData(-1, 0, CORE));
    (*std_functions)["narginchk"]          = BuiltinFunc(oml_narginchk, FunctionMetaData(2, 0, CORE));
    (*std_functions)["nargoutchk"]         = BuiltinFunc(oml_nargoutchk, FunctionMetaData(2, 0, CORE));
    (*std_functions)["fileparts"]          = BuiltinFunc(oml_fileparts, FunctionMetaData(1, 3, FILEIO));
    (*std_functions)["subsref"]            = BuiltinFunc(oml_subsref, FunctionMetaData(2, 1, CORE));
    (*std_functions)["subsasgn"]           = BuiltinFunc(oml_subsasgn, FunctionMetaData(3, 1, DATA));
    (*std_functions)["flintmax"]           = BuiltinFunc(oml_flintmax, FunctionMetaData(1, 1, CORE));
    (*std_functions)["mislocked"]          = BuiltinFunc(oml_mislocked, FunctionMetaData(1, 1, CORE));
    (*std_functions)["mlock"]              = BuiltinFunc(oml_mlock, FunctionMetaData(0, 0, CORE));
    (*std_functions)["munlock"]            = BuiltinFunc(oml_munlock, FunctionMetaData(1, 0, CORE));
    (*std_functions)["ferror"]             = BuiltinFunc(oml_ferror, FunctionMetaData(2, 2, FILEIO));
    (*std_functions)["meshgrid"]           = BuiltinFunc(oml_meshgrid, FunctionMetaData(-2, -1, PLOT));
    (*std_functions)["ndgrid"]             = BuiltinFunc(oml_ndgrid, FunctionMetaData(-2, -1, CORE));
	(*std_functions)["exit"]               = BuiltinFunc(oml_exit, FunctionMetaData(0, 0, CORE));
	(*std_functions)["quit"]               = BuiltinFunc(oml_exit, FunctionMetaData(0, 0, CORE));
	(*std_functions)["cputime"]            = BuiltinFunc(oml_cputime, FunctionMetaData(0, 1, SYSTEM));
	(*std_functions)["diary"]              = BuiltinFunc(oml_diary, FunctionMetaData(1, 0, SYSTEM));
	(*std_functions)["numel"]              = BuiltinFunc(oml_numel, FunctionMetaData(-2, 1, ELEM));
    (*std_functions)["rot90"]              = BuiltinFunc(oml_rot90, FunctionMetaData(-2, 1, ELEM));
	(*std_functions)["help"]			   = BuiltinFunc(oml_help, FunctionMetaData(1,1, CORE));

    // Elementary math functions
    (*std_functions)["dec2hex"]            = BuiltinFunc(BuiltInFuncsConvert::hml_dec2hex,
                                                         FunctionMetaData(2, 1, ELEM));
    (*std_functions)["hex2dec"]            = BuiltinFunc(BuiltInFuncsConvert::hml_hex2dec,
                                                         FunctionMetaData(1, 1, ELEM));
    (*std_functions)["dec2bin"]            = BuiltinFunc(BuiltInFuncsConvert::hml_dec2bin,
                                                         FunctionMetaData(2, 1, ELEM));
    (*std_functions)["bin2dec"]            = BuiltinFunc(BuiltInFuncsConvert::hml_bin2dec,
                                                         FunctionMetaData(1, 1, ELEM));
    (*std_functions)["flip"]      = BuiltinFunc(BuiltInFuncsElemMath::Flip,   FunctionMetaData(-2, 1, ELEM));
    (*std_functions)["fliplr"]    = BuiltinFunc(BuiltInFuncsElemMath::Fliplr, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["flipud"]    = BuiltinFunc(BuiltInFuncsElemMath::Flipud, FunctionMetaData(1, 1, ELEM));

    // Core functions
    (*std_functions)["class"]     = BuiltinFunc(BuiltInFuncsCore::Class, FunctionMetaData(1, 1, CORE));
    (*std_functions)["isa"]       = BuiltinFunc(BuiltInFuncsCore::IsA, FunctionMetaData(2, 1, CORE));
    (*std_functions)["logspace"]  = BuiltinFunc(oml_logspace,               FunctionMetaData(3, 1, CORE));
	(*std_functions)["format"]    = BuiltinFunc(BuiltInFuncsCore::Format,   FunctionMetaData(1, 0, CORE));
    (*std_functions)["input"]     = BuiltinFunc(BuiltInFuncsCore::Input,    FunctionMetaData(2, 1, CORE));
    (*std_functions)["pause"]     = BuiltinFunc(BuiltInFuncsCore::Pause,    FunctionMetaData(1, 0, CORE));
    (*std_functions)["getenv"]    = BuiltinFunc(oml_getenv, FunctionMetaData(1, 1, CORE));
    (*std_functions)["putenv"]    = BuiltinFunc(oml_setenv, FunctionMetaData(2, 0, CORE));
    (*std_functions)["setenv"]    = BuiltinFunc(oml_setenv, FunctionMetaData(2, 0, CORE));
    (*std_functions)["whos"]      = BuiltinFunc(BuiltInFuncsCore::Whos,     FunctionMetaData(-1, 1, CORE));
    (*std_functions)["addtoolbox"] = BuiltinFunc(BuiltInFuncsCore::AddToolbox,  FunctionMetaData(1, 0, CORE));
    (*std_functions)["addToolbox"] = BuiltinFunc(BuiltInFuncsCore::AddToolbox,  FunctionMetaData(1, 0, CORE));
    (*std_functions)["funccount"]  = BuiltinFunc(BuiltInFuncsCore::Funccount,   FunctionMetaData(0, 1, CORE));
    (*std_functions)["funclist"]   = BuiltinFunc(BuiltInFuncsCore::Funclist,  FunctionMetaData(-1, 0, CORE));
    (*std_functions)["varlist"]    = BuiltinFunc(BuiltInFuncsCore::Varlist,   FunctionMetaData(0, 1, CORE));
	(*std_functions)["omlfilename"] = BuiltinFunc(BuiltInFuncsCore::Omlfilename, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["exist"]       = BuiltinFunc(BuiltInFuncsCore::Exist, FunctionMetaData(2, 1, CORE));
    (*std_functions)["type"]        = BuiltinFunc(BuiltInFuncsCore::Type,  FunctionMetaData(-1, 1, CORE));
    (*std_functions)["omlpaginate"] = BuiltinFunc(BuiltInFuncsCore::OmlPaginate, FunctionMetaData(1, 0, CORE));
    (*std_functions)["sleep"]       = BuiltinFunc(BuiltInFuncsCore::OmlSleep, FunctionMetaData(1, 0, TIME));
    (*std_functions)["errormsgonly"] = BuiltinFunc(oml_erromsgonly, FunctionMetaData(-2, 0, CORE));

    // Data structures
    (*std_functions)["cat"]       = BuiltinFunc(oml_cat, FunctionMetaData(-2, 1, DATA));
    (*std_functions)["setfield"]  = BuiltinFunc(BuiltInFuncsData::Setfield, FunctionMetaData(-3, 1, DATA));
    (*std_functions)["mat2cell"]  = BuiltinFunc(BuiltInFuncsData::Mat2Cell, FunctionMetaData(-2, 1, DATA));

    // File I/O
	(*std_functions)["textread"]  = BuiltinFunc(BuiltInFuncsFile::Textread, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["dlmwrite"]  = BuiltinFunc(BuiltInFuncsFile::Dlmwrite, FunctionMetaData(-5, 1, FILEIO));
    (*std_functions)["copyfile"]  = BuiltinFunc(BuiltInFuncsFile::Copyfile, FunctionMetaData(3, 3, FILEIO));
    (*std_functions)["importdata"] = BuiltinFunc(BuiltInFuncsFile::Importdata, FunctionMetaData(-1, -1, FILEIO));
	(*std_functions)["textscan"]   = BuiltinFunc(BuiltInFuncsFile::Textscan, FunctionMetaData(-1, 1, FILEIO));
    (*std_functions)["fscanf"]     = BuiltinFunc(BuiltInFuncsFile::Fscanf, FunctionMetaData(-3, 1, FILEIO));
    (*std_functions)["rename"]     = BuiltinFunc(BuiltInFuncsFile::Rename, FunctionMetaData(2, 2, FILEIO));

    // String functions
    (*std_functions)["blanks"]   = BuiltinFunc(BuiltInFuncsString::hml_blanks, 
                                       FunctionMetaData(1, 1, STNG));
    (*std_functions)["strtok"]   = BuiltinFunc(BuiltInFuncsString::hml_strtok, 
                                       FunctionMetaData(2, 2, STNG));
    (*std_functions)["strvcat"]  = BuiltinFunc(BuiltInFuncsString::hml_strvcat,
                                       FunctionMetaData(1, 1, STNG));
    (*std_functions)["sscanf"]   = BuiltinFunc(BuiltInFuncsString::Sscanf,
                                       FunctionMetaData(-3, 1, STNG));
    (*std_functions)["toascii"]  = BuiltinFunc(BuiltInFuncsString::ToAscii, FunctionMetaData(1, 1, STNG));
    (*std_functions)["strtrim"]  = BuiltinFunc(BuiltInFuncsString::Strtrim, FunctionMetaData(1, 1, STNG));
    (*std_functions)["mat2str"]  = BuiltinFunc(BuiltInFuncsString::Mat2Str, FunctionMetaData(-2, 1, STNG));
    (*std_functions)["regexprep"] = BuiltinFunc(BuiltInFuncsString::Regexprep, FunctionMetaData(1, -4, STNG));
    (*std_functions)["circshift"] = BuiltinFunc(oml_circshift, FunctionMetaData(1, -4, DATA));
    (*std_functions)["checksyntax"] = BuiltinFunc(oml_checksyntax, FunctionMetaData(1, 1, CORE));
    (*std_functions)["ast"]         = BuiltinFunc(oml_ast, FunctionMetaData(1, 1, CORE));
    (*std_functions)["num2str"]     = BuiltinFunc(BuiltInFuncsString::Num2Str, FunctionMetaData(2, 1, STNG));
    (*std_functions)["str2mat"]     = BuiltinFunc(BuiltInFuncsString::Str2mat, FunctionMetaData(1, -1, STNG));

    // System functions
    (*std_functions)["dir"]      = BuiltinFunc(BuiltInFuncsSystem::Dir, FunctionMetaData(1, 1, SYSTEM));
    (*std_functions)["ls"]       = BuiltinFunc(BuiltInFuncsSystem::Ls, FunctionMetaData(-1,-1, SYSTEM));
    (*std_functions)["startupscript_edit"] = BuiltinFunc(DummyVoid, FunctionMetaData(0, 0, SYSTEM));
    (*std_functions)["system"]   = BuiltinFunc(BuiltInFuncsSystem::System, FunctionMetaData(-2, -2, SYSTEM));
    (*std_functions)["dos"]      = BuiltinFunc(BuiltInFuncsSystem::System, FunctionMetaData(-2, -2, SYSTEM));
    (*std_functions)["unix"]     = BuiltinFunc(BuiltInFuncsSystem::Unix,   FunctionMetaData(2, 2, SYSTEM));
    (*std_functions)["delete"]   = BuiltinFunc(BuiltInFuncsSystem::Delete, FunctionMetaData(-1, 0, SYSTEM));

    // Client specific environment related functions
    (*std_functions)["getbaseenv"]    = BuiltinFunc(BuiltInFuncsSystem::GetBaseEnv, 
                                        FunctionMetaData(OML_NO_NARG, OML_NO_NARG));
    (*std_functions)["getcurrentenv"] = BuiltinFunc(BuiltInFuncsSystem::GetCurrentEnv, 
                                        FunctionMetaData(OML_NO_NARG, OML_NO_NARG));
    (*std_functions)["getnewenv"]     = BuiltinFunc(BuiltInFuncsSystem::GetNewEnv, 
                                        FunctionMetaData(OML_NO_NARG, OML_NO_NARG));
    (*std_functions)["getenvvalue"]   = BuiltinFunc(BuiltInFuncsSystem::GetEnvVal, 
                                        FunctionMetaData(OML_NO_NARG, OML_NO_NARG));
    (*std_functions)["setenvvalue"]   = BuiltinFunc(BuiltInFuncsSystem::SetEnvVal, 
                                        FunctionMetaData(OML_NO_NARG, OML_NO_NARG));
    (*std_functions)["clearenvvalue"] = BuiltinFunc(BuiltInFuncsSystem::ClearEnvVal, 
                                        FunctionMetaData(OML_NO_NARG, OML_NO_NARG));
    (*std_functions)["importenv"]     = BuiltinFunc(BuiltInFuncsSystem::ImportEnv, 
                                        FunctionMetaData(OML_NO_NARG, OML_NO_NARG));
    (*std_functions)["importenvin"]   = BuiltinFunc(BuiltInFuncsSystem::ImportEnvIn, 
                                        FunctionMetaData(OML_NO_NARG, OML_NO_NARG));
    (*std_functions)["cloneenv"]      = BuiltinFunc(BuiltInFuncsSystem::CloneEnv, 
                                        FunctionMetaData(OML_NO_NARG, OML_NO_NARG));
    // End of client specific environment related functions

    // GUI
    (*std_functions)["msgbox"]     = BuiltinFunc(dummyint,    FunctionMetaData(-2, -1, GUI));
    (*std_functions)["errordlg"]   = BuiltinFunc(dummyint,    FunctionMetaData(-2, -1, GUI));
    (*std_functions)["warndlg"]    = BuiltinFunc(dummyint,    FunctionMetaData(-2, -1, GUI));
    (*std_functions)["uigetfile"]  = BuiltinFunc(dummystring, FunctionMetaData(-3, -3, GUI));
    (*std_functions)["inputdlg"]   = BuiltinFunc(dummycell,   FunctionMetaData(-4, 1, GUI));
    (*std_functions)["uigetdir"]   = BuiltinFunc(dummystring, FunctionMetaData(-2, 1, GUI));
    (*std_functions)["questdlg"]   = BuiltinFunc(dummystring, FunctionMetaData(-6, 1, GUI));
    (*std_functions)["uiputfile"]  = BuiltinFunc(dummystring, FunctionMetaData(-3, -3, GUI));
    (*std_functions)["edit"]       = BuiltinFunc(DummyVoid,   FunctionMetaData(0, 1, GUI));
	(*std_functions)["memoryuse"]  = BuiltinFunc(oml_memoryuse, FunctionMetaData(0, 1, CORE));
	(*std_functions)["analyze"]    = BuiltinFunc(oml_analyze, FunctionMetaData(1, 0, CORE));

	bool experimental = true;

#if !defined(_DEBUG)
	experimental = Currency::GetExperimental();
#endif

	if (experimental)
	{
		(*std_functions)["writepfile"] = BuiltinFunc(oml_writepfile, FunctionMetaData(2, 0, CORE));

		// helptest function is restricted to experimental mode
	    (*std_functions)["helptest"]      = BuiltinFunc(oml_helptest, FunctionMetaData(1,1, CORE));
	    (*std_functions)["helpmodule"]      = BuiltinFunc(oml_helpmodule, FunctionMetaData(1,1, CORE));
	}
}
//------------------------------------------------------------------------------
// Display the error status of fileID [ferror]
//------------------------------------------------------------------------------
bool oml_ferror(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs, 
                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin > 1 && readOption(eval, inputs[1]) != "clear")
        throw OML_Error(HW_ERROR_2NDARGONLYCLR);

    int fileID = getFileFromInput(eval, inputs[0]);
    checkFileIndex(eval, fileID, true);
    FILE* f = eval.GetFile(fileID);
    
    int err = ferror(f);
    if (eval.GetNargoutValue() > 1)
        outputs.push_back(err ? 1.0 : 0.0);
    outputs.push_back(err ? std::string(strerror(ferror(f))) : std::string());

    if (nargin > 1)
        clearerr(f);

    return true;
}
//------------------------------------------------------------------------------
// Returns true and locks a function into memory [mlock]
//------------------------------------------------------------------------------
bool oml_mlock(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    eval.LockCurrent();
    return true;
}
//------------------------------------------------------------------------------
// Returns true and unlocks current/given function [munlock]
//------------------------------------------------------------------------------
bool oml_munlock(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin > 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin)
    {
        if (inputs[0].IsString())
        {
            std::string fname = readString(inputs[0]);
            eval.Unlock(fname);
        }
        else
            throw OML_Error(HW_ERROR_INPUTSTRING);
    }
    else
    {
        eval.UnlockCurrent();
    }

    return true;
}
bool oml_mislocked(EvaluatorInterface           eval, 
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin > 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin)
    {
        if (inputs[0].IsString())
        {
            std::string fname = readString(inputs[0]);
            outputs.push_back(HW_MAKE_BOOL_CURRENCY(eval.IsLocked(fname)));
        }
        else
            throw OML_Error(HW_ERROR_INPUTSTRING);
    }
    else
    {
        outputs.push_back(HW_MAKE_BOOL_CURRENCY(eval.IsCurrentLocked()));
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true and if given function is locked [mislocked]
//------------------------------------------------------------------------------
bool oml_flintmax(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin > 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    bool single = false;

    if (nargin)
    {
        if (inputs[0].IsString())
        {
            std::string opt = readOption(eval, inputs[0]);
            if (opt == "single")
                single = true;
            else if (opt == "double")
                single = false;
            else
                throw OML_Error(HW_ERROR_INVALIDOPTION(opt));
        }
        else
            throw OML_Error(HW_ERROR_INPUTSTRING);
    }

    if (single)
        outputs.push_back(std::pow(2.0f, FLT_MANT_DIG));
    else
        outputs.push_back(std::pow(2.0, DBL_MANT_DIG));

    return true;
}
//------------------------------------------------------------------------------
// Verifies field names
//------------------------------------------------------------------------------
static void VerifySubsStructFields(const StructData* sd)
{
    const std::map<std::string,int>& fields = sd->GetFieldNames();
    if (fields.size() != 2 || !fields.count("type") || !fields.count("subs"))
        throw OML_Error(HW_ERROR_INVSTRUCTFIELDTYPESUB);
}
//------------------------------------------------------------------------------
// Assigns a value to an element of a cell or matrix [subsasgn]
//------------------------------------------------------------------------------
bool oml_subsasgn(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency* index = nullptr;

    if (!inputs[1].IsStruct())
        throw OML_Error(HW_ERROR_INPUTSTRUCT);

    const StructData* sd = inputs[1].Struct();
    VerifySubsStructFields(sd);
    int size = sd->M()*sd->N();
    
    if (size == 0)
    {
        outputs.push_back(inputs[2]);
        return true;
    }

    if (size > 1)
        throw OML_Error(HW_ERROR_UNSUPINDCHAIN);

    Currency val = inputs[0];

    for (int i = 0; i < size; ++i)
    {
        const Currency& type = sd->GetValue(i, -1, std::string("type"));
        const Currency& subs = sd->GetValue(i, -1, std::string("subs"));
        if (!type.IsString() || type.Matrix()->M() != 1)
            throw OML_Error(HW_ERROR_TYPE1DSTR);

        std::string str_type = type.StringVal();
        if (str_type == ".")
        {
            if (!val.IsStruct())
                throw OML_Error(HW_NONSTRUCTFIELD);

            if (!subs.IsString() || subs.Matrix()->M() != 1)
                throw OML_Error(HW_ERROR_INVSTRUCTFIELD);

            StructData* sd_val = val.Struct();

            if (sd_val->M() != 1 || sd_val->N() != 1)
                throw OML_Error(HW_ERROR_NOTINDSTRUCT1X1);

            //if (i < size-1)
                //val = (Currency*) &sd_val->GetValue(0, 0, subs.StringVal());
            //else
                sd_val->SetValue(0, 0, subs.StringVal(), inputs[2]);
            continue;
        }
        else
        {
            std::vector<Currency> index_params;
            if (subs.IsCharacter() && subs.StringVal() == ":")
            {
                index_params.push_back(Currency(0.0, Currency::TYPE_COLON));
            }
            else if (subs.IsCellArray())
            {
                HML_CELLARRAY* cell = subs.CellArray();
                for (int i = 0; i < cell->Size(); ++i)
                {
                    const Currency& idx = (*cell)(i);
                    if (idx.IsString() && idx.StringVal() == ":")
                    {
                        index_params.push_back(Currency(0.0, Currency::TYPE_COLON));
                    }
                    else
                    {
                        index_params.push_back(idx);
                    }
                }
            }
            else
                throw OML_Error(HW_ERROR_SUBMUSTCELL);

            if (str_type == "()")
            {
                 outputs.push_back(eval.Assignment(val, index_params, inputs[2]));
				 return true;
            }
            else if (str_type == "{}")
            {

                if (!val.IsCellArray())
                    throw OML_Error(HW_ERROR_INVINDTYPE);

                outputs.push_back(eval.CellAssignment(val, index_params, inputs[2]));
				return true;
            }
            else
                throw OML_Error(HW_ERROR_DOTSTRUCTPARENTHESIS);
        }
    }

    outputs.push_back(val);
    return true;
}
//------------------------------------------------------------------------------
// Extracts a subset of a collection with given indexing method/range [subsref]
//------------------------------------------------------------------------------
bool oml_subsref(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency val = inputs[0];
    if (!inputs[1].IsStruct())
        throw OML_Error(HW_ERROR_INPUTSTRUCT);

    const StructData* sd = inputs[1].Struct();
    VerifySubsStructFields(sd);
    int size = sd->M() * sd->N();
    for (int i = 0 ; i < size; ++i)
    {
        const Currency& type = sd->GetValue(i, -1, std::string("type"));
        const Currency& subs = sd->GetValue(i, -1, std::string("subs"));
        if (!type.IsString() || type.Matrix()->M() != 1)
            throw OML_Error(HW_ERROR_TYPE1DSTR);

        std::string str_type = type.StringVal();
        if (str_type == ".")
        {
            if (!val.IsStruct())
                throw OML_Error(HW_NONSTRUCTFIELD);

            if (!subs.IsString() || subs.Matrix()->M() != 1)
                throw OML_Error(HW_ERROR_INVSTRUCTFIELD);

            const StructData* sd_val = val.Struct();

            if (sd_val->M() != 1 || sd_val->N() != 1)
                throw OML_Error(HW_ERROR_NOTINDSTRUCT1X1);

            val = sd_val->GetValue(0, 0, subs.StringVal());

            continue;
        }
        else
        {
            std::vector<Currency> index_params;
            if (subs.IsCharacter() && subs.StringVal() == ":")
            {
                index_params.push_back(Currency(0.0, Currency::TYPE_COLON));
            }
            else if (subs.IsCellArray())
            {
                HML_CELLARRAY* cell = subs.CellArray();
                for (int i = 0; i < cell->Size(); ++i)
                {
                    const Currency& idx = (*cell)(i);
                    if (idx.IsString() && idx.StringVal() == ":")
                    {
                        index_params.push_back(Currency(0.0, Currency::TYPE_COLON));
                    }
                    else
                    {
                        index_params.push_back(idx);
                    }
                }
            }
            else
                throw OML_Error(HW_ERROR_SUBMUSTCELL);

            if (str_type == "()")
            {
                val = eval.VariableIndex(val, index_params);
            }
            else if (str_type == "{}")
            {
                if (!val.IsCellArray())
                    throw OML_Error(HW_ERROR_INVINDTYPE);
                val = eval.CellIndex(val, index_params);
            }
            else
                throw OML_Error(HW_ERROR_INVTYPEVALDOTPARCURL);
        }
    }
    
    outputs.push_back(val);
    return true;
}
//------------------------------------------------------------------------------
// Gets directory, name and extension components given file name [fileparts]
//------------------------------------------------------------------------------
bool oml_fileparts(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsString())
        throw OML_Error(HW_ERROR_INPUTSTRING);

    if (inputs[0].Matrix()->M() > 1)
        throw OML_Error(HW_ERROR_TYPE1DSTR);

    std::string filename = inputs[0].StringVal();

#if OS_WIN
    size_t lastslash = filename.find_last_of("/\\");
#else
    size_t lastslash = filename.find_last_of('/');
#endif

    if (lastslash == std::string::npos)
    {
        outputs.push_back(std::string());
        lastslash = 0;
    }
    else
    {
        outputs.push_back(filename.substr(0, lastslash));
        ++lastslash;
    }

    size_t lastdot = filename.find_last_of('.');
    if (lastdot < lastslash || lastdot == std::string::npos)
    {
        outputs.push_back(filename.substr(lastslash));
        outputs.push_back(std::string());
    }
    else
    {
        outputs.push_back(filename.substr(lastslash, lastdot-lastslash));
        outputs.push_back(filename.substr(lastdot));
    }

    return true;
}
//------------------------------------------------------------------------------
// Checks for correct number of inputs [narginchk]
//------------------------------------------------------------------------------
bool oml_narginchk(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsInteger())
        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_VALUE);

    if (inputs[0].Scalar() == 0.0)
        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_VALUE);

    if (!inputs[1].IsInteger())
        throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_VALUE);

    if (inputs[1].Scalar() == 0.0)
        throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_VALUE);

    int min = (int) inputs[0].Scalar();
    int max = (int) inputs[1].Scalar();

    EvaluatorInterface caller(eval, false);
    int actual = caller.GetNarginValue();

    if (actual < min)
        throw OML_Error(OML_ERR_NUMARGIN);
    if (actual > max)
        throw OML_Error(OML_ERR_NUMARGIN);
    return true;
}
//------------------------------------------------------------------------------
// Checks for correct number of outputs [nargoutchk]
//------------------------------------------------------------------------------
bool oml_nargoutchk(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsInteger())
        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_VALUE);

    if (inputs[0].Scalar() == 0.0)
        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_VALUE);

    if (!inputs[1].IsInteger())
        throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_VALUE);

    if (inputs[1].Scalar() == 0.0)
        throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_VALUE);

    int min = (int) inputs[0].Scalar();
    int max = (int) inputs[1].Scalar();
    EvaluatorInterface caller(eval, false);
    int actual = caller.GetNargoutValue();

    if (actual < min)
        throw OML_Error(OML_ERR_NUMARGOUT);
    if (actual > max)
        throw OML_Error(OML_ERR_NUMARGOUT);
    return true;
}
//------------------------------------------------------------------------------
// Returns an error if an assertion fails [assert]
//------------------------------------------------------------------------------
bool oml_assert(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin == 1 || (inputs[0].IsLogical() && inputs[1].IsString()))
    {
        bool retval = true;
        const Currency& cur = inputs[0];

        if (cur.IsScalar())
            retval = iszero(cur.Scalar());
        else if (cur.IsComplex())
            retval = iszero(cur.Complex());
        else if (cur.IsMatrix())
        {
            const hwMatrix* mtx = cur.Matrix();
            if (!mtx || mtx->Size() == 0)
                retval = true;
			else
				retval = (static_cast<int>(all(eval, cur.Matrix()))) ? false : true;
        }
        else if (cur.IsNDMatrix())
        {
            return oml_MatrixNUtil1(eval, inputs, outputs, oml_assert);
        }

        if (retval)
        {
            if (nargin == 1)
            {
                throw OML_Error(HW_ERROR_ASSERTFAIL);
            }
            else
            {
                _OML_Error(eval, inputs.cbegin() + 1, inputs.cend());
            }
        }
    }
    else // compare values
    {
        if (nargin == 2) // no tolerance
        {
            if (!isequal(inputs[0], inputs[1]))
                throw OML_Error(HW_ERROR_ASSERTFAIL);
        }
        else
        {
            if (!isequal(inputs[0], inputs[1], inputs[2]))
                throw OML_Error(HW_ERROR_ASSERTFAIL);
        }
    }
    return true;
}
//------------------------------------------------------------------------------
// Calls given builtin function [builtin]
//------------------------------------------------------------------------------
bool oml_builtin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() < 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsString())
        throw OML_Error(HW_ERROR_INPUTSTRING);

    std::string fname = readString(inputs[0]);
    FUNCPTR fptr = eval.GetStdFunction(fname);
    if (!fptr)
        throw OML_Error("Error: invalid function name: " + fname);

    return (*fptr)(eval, std::vector<Currency>(inputs.cbegin() + 1, inputs.cend()), outputs);
}
//------------------------------------------------------------------------------
// Returns the inverse tangent in radians [atan2]
//------------------------------------------------------------------------------
bool oml_atan2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsScalar())
    {
        double y = inputs[0].Scalar();

        if (inputs[1].IsScalar())
        {
            outputs.push_back(std::atan2(y, inputs[1].Scalar()));
        }
        else if (inputs[1].IsMatrix() && inputs[1].Matrix()->IsRealData())
        {
            const hwMatrix* x = inputs[1].Matrix();
            hwMatrix* result = EvaluatorInterface::allocateMatrix(x->M(), x->N(), hwMatrix::REAL);

            if (x->IsReal())
            {
                for (int i = 0; i < x->Size(); ++i)
                {
                    (*result)(i) = std::atan2(y, (*x)(i));
                }
            }
            else
            {
                for (int i = 0; i < x->Size(); ++i)
                {
                    (*result)(i) = std::atan2(y, x->z(i).Real());
                }
            }
            outputs.push_back(result);
        }
        else if (inputs[1].IsNDMatrix())
        {
            return oml_MatrixNUtil2(eval, inputs, outputs, oml_atan2);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
    }
    else if (inputs[0].IsMatrix() && inputs[0].Matrix()->IsRealData())
    {
        const hwMatrix* y = inputs[0].Matrix();

        if (inputs[1].IsScalar())
        {
            hwMatrix* result = EvaluatorInterface::allocateMatrix(y->M(), y->N(), hwMatrix::REAL);

            if (y->IsReal())
            {
                for (int i = 0; i < y->Size(); ++i)
                {
                    (*result)(i) = std::atan2((*y)(i), inputs[1].Scalar());
                }
            }
            else
            {
                for (int i = 0; i < y->Size(); ++i)
                {
                    (*result)(i) = std::atan2(y->z(i).Real(), inputs[1].Scalar());
                }
            }
            outputs.push_back(result);
        }
        else if (inputs[1].IsMatrix() && inputs[1].Matrix()->IsRealData())
        {
            const hwMatrix* x = inputs[1].Matrix();

            if (!sameSize(y,x))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            hwMatrix* result = EvaluatorInterface::allocateMatrix(y->M(), y->N(), hwMatrix::REAL);

            for (int i = 0; i < y->Size(); ++i)
            {
                (*result)(i) = std::atan2(realval(y,i), realval(x,i));
            }
            outputs.push_back(result);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_DIM);
    }
    else if (inputs[0].IsNDMatrix() || inputs[1].IsNDMatrix())
    {
        return oml_MatrixNUtil2(eval, inputs, outputs, oml_atan2);
    }
    else
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DIM);

    return true;
}
//------------------------------------------------------------------------------
// Returns true and joins path names to build complete filename [fullfile]
//------------------------------------------------------------------------------
bool oml_fullfile(EvaluatorInterface           eval, 
                  const std::vector<Currency>& inputs, 
                  std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.empty() ? 0 : inputs.size();

    if (nargin < 1) throw OML_Error(OML_ERR_NUMARGIN);

    bool hasslash = false;
    std::stringstream ss;
    for (size_t i = 0; i < nargin; ++i)
    {
        const Currency& cur = inputs[i];

        if (!cur.IsString()) 
            throw OML_Error(OML_ERR_STRING, static_cast<int>(i)+1, OML_VAR_TYPE);
        
        if (cur.Matrix()->IsEmpty()) continue;
        
        if (cur.Matrix()->M() != 1)  throw OML_Error(HW_ERROR_INPROWVECT);

        std::string s = cur.StringVal();
        if (s.length())
        {
            if (i && !hasslash)
                ss << DIRECTORY_DELIM;

            ss << s;
            hasslash = s.back() == DIRECTORY_DELIM;
        }
    }

    outputs.push_back(ss.str());
    return true;
}
//------------------------------------------------------------------------------
// Creates cell array from struct [struct2cell]
//------------------------------------------------------------------------------
bool oml_struct2cell(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsStruct())
        throw OML_Error(HW_ERROR_INPUTSTRUCT);

    const StructData* sd = inputs[0].Struct();

    int m = sd->M();
    int n = sd->N();

    if (m && n)
    {
        if (n > 1)
            throw OML_Error(OML_ERR_UNSUPPORTDIM, 1);

        const std::map<std::string, int> fields = sd->GetFieldNames();
        HML_CELLARRAY* cell = EvaluatorInterface::allocateCellArray((int)fields.size(), m);
        std::map<std::string, int>::const_iterator iter;

        for (int i = 0; i < m; ++i)
        {
            int f = 0;
            for (iter = fields.cbegin(); iter != fields.cend(); ++iter)
            {
                (*cell)(f++, i) = sd->GetValue(i,0,iter->first);
            }
        }
        outputs.push_back(cell);
    }
    else // struct is size 0
    {
        outputs.push_back(EvaluatorInterface::allocateCellArray());
    }
    return true;
}
//------------------------------------------------------------------------------
// Converts Cartesian coordinates to spherical [cart2sph]
//------------------------------------------------------------------------------
bool oml_cart2sph(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    int nargout = getNumOutputs(eval);

    if (nargin != 1 && nargin != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin == 1)
    {
        if (!inputs[0].IsMatrix())
            throw OML_Error(HW_ERROR_1INPMUSTMAT);

        const hwMatrix* mtx = inputs[0].Matrix();
        if (mtx->N() != 3)
            throw OML_Error(HW_ERROR_MATMUST3COL);

        bool usez = mtx->N() > 2;
        if (nargout > 1)
        {
            hwMatrix* x = EvaluatorInterface::allocateMatrix(mtx->M(), 1, mtx->Type());
            hwMatrix* y = EvaluatorInterface::allocateMatrix(mtx->M(), 1, mtx->Type());
            hwMatrix* z = EvaluatorInterface::allocateMatrix(mtx->M(), 1, mtx->Type());

            if (mtx->IsReal())
            {
                for (int i = 0; i < mtx->M(); ++i)
                    cart2sphHelper((*mtx)(i,0), (*mtx)(i,1), (*mtx)(i,2), (*x)(i), (*y)(i), (*z)(i));
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_MATRIX);
            }

            outputs.push_back(x);
            outputs.push_back(y);
            outputs.push_back(z);
        }

        else
        {
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

            if (mtx->IsReal())
            {
                for (int i = 0; i < mtx->M(); ++i)
                    cart2sphHelper((*mtx)(i,0), (*mtx)(i,1), (*mtx)(i,2), (*result)(i,0), (*result)(i,1), (*result)(i,2));
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_MATRIX);
            }

            outputs.push_back(result);
        }
    }
    else if (nargin > 1)
    {
        if (nargout > 1)
        {
            outputs = mtxFun(eval, inputs, (int)nargin, &cart2sphHelper);
        }
        else
        {
            std::vector<Currency> results = mtxFun(eval, inputs, (int)nargin, &cart2sphHelper);
            const hwMatrix* x = results[0].Matrix();
            const hwMatrix* y = results[1].Matrix();
            const hwMatrix* z = results[2].Matrix();
            hwMatrix* outmtx = EvaluatorInterface::allocateMatrix(x->Size(), (int)nargin, x->Type());
            
            if (x->IsReal())
            {
                for (int i = 0; i < outmtx->M(); ++i)
                {
                    (*outmtx)(i,0) = (*x)(i);
                    (*outmtx)(i,1) = (*y)(i);
                    (*outmtx)(i,2) = (*z)(i);
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL);
            }

            outputs.push_back(outmtx);
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Converts spherical coordinates to Cartesian [sph2cart]
//------------------------------------------------------------------------------
bool oml_sph2cart(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    int nargout = getNumOutputs(eval);

    if (nargin != 1 && nargin != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin == 1)
    {
        if (!inputs[0].IsMatrix())
            throw OML_Error(HW_ERROR_1INPMUSTMAT);

        const hwMatrix* mtx = inputs[0].Matrix();
        if (mtx->N() != 3)
            throw OML_Error(HW_ERROR_MATMUST3COL);

        bool usez = mtx->N() > 2;
        if (nargout > 1)
        {
            hwMatrix* x = EvaluatorInterface::allocateMatrix(mtx->M(), 1, mtx->Type());
            hwMatrix* y = EvaluatorInterface::allocateMatrix(mtx->M(), 1, mtx->Type());
            hwMatrix* z = EvaluatorInterface::allocateMatrix(mtx->M(), 1, mtx->Type());

            if (mtx->IsReal())
            {
                for (int i = 0; i < mtx->M(); ++i)
                    sph2cartHelper((*mtx)(i,0), (*mtx)(i,1), (*mtx)(i,2), (*x)(i), (*y)(i), (*z)(i));
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_MATRIX);
            }

            outputs.push_back(x);
            outputs.push_back(y);
            outputs.push_back(z);
        }

        else
        {
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

            if (mtx->IsReal())
            {
                for (int i = 0; i < mtx->M(); ++i)
                    sph2cartHelper((*mtx)(i,0), (*mtx)(i,1), (*mtx)(i,2), (*result)(i,0), (*result)(i,1), (*result)(i,2));
            }
            else
            {
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_MATRIX);
            }

            outputs.push_back(result);
        }
    }
    else if (nargin > 1)
    {
        if (nargout > 1)
        {
            outputs = mtxFun(eval, inputs, (int)nargin, &sph2cartHelper);
        }
        else
        {
            std::vector<Currency> results = mtxFun(eval, inputs, (int)nargin, &sph2cartHelper);
            const hwMatrix* x = results[0].Matrix();
            const hwMatrix* y = results[1].Matrix();
            const hwMatrix* z = results[2].Matrix();
            hwMatrix* outmtx = EvaluatorInterface::allocateMatrix(x->Size(), (int)nargin, x->Type());
            
            if (x->IsReal())
            {
                for (int i = 0; i < outmtx->M(); ++i)
                {
                    (*outmtx)(i,0) = (*x)(i);
                    (*outmtx)(i,1) = (*y)(i);
                    (*outmtx)(i,2) = (*z)(i);
                }
            }
            else
            {
                throw OML_Error(OML_ERR_REAL);
            }

            outputs.push_back(outmtx);
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Converts Cartesian coordinates to polar or cylindrical [cart2pol]
//------------------------------------------------------------------------------
bool oml_cart2pol(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    int nargout = getNumOutputs(eval);

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin == 1)
    {
        if (!inputs[0].IsMatrix())
            throw OML_Error(HW_ERROR_1INPMUSTMAT);

        const hwMatrix* mtx = inputs[0].Matrix();
        if (mtx->N() != 2 && mtx->N() != 3)
            throw OML_Error(HW_ERROR_MATMUST2OR3COL);

        bool usez = mtx->N() > 2;
        if (mtx->IsReal())
        {
            if (nargout > 1)
            {
                hwMatrix* x = EvaluatorInterface::allocateMatrix(mtx->M(), 1, hwMatrix::REAL);
                hwMatrix* y = EvaluatorInterface::allocateMatrix(mtx->M(), 1, hwMatrix::REAL);
                hwMatrix* z;
                if (usez)
                    z = EvaluatorInterface::allocateMatrix(mtx->M(), 1, hwMatrix::REAL);

                for (int i = 0; i < mtx->M(); ++i)
                {
                    std::pair<double,double> singleResult = cart2polHelper((*mtx)(i,0), (*mtx)(i,1));

                    (*x)(i) = singleResult.first;
                    (*y)(i) = singleResult.second;
                    if (usez)
                        (*z)(i) = (*mtx)(i,2);
                }
                outputs.push_back(x);
                outputs.push_back(y);
                if (usez)
                    outputs.push_back(z);
                else
                    outputs.push_back(EvaluatorInterface::allocateMatrix());
            }
            else
            {
                hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

                for (int i = 0; i < mtx->M(); ++i)
                {
                    std::pair<double,double> singleResult = cart2polHelper((*mtx)(i,0), (*mtx)(i,1));

                    (*result)(i,0) = singleResult.first;
                    (*result)(i,1) = singleResult.second;
                    if (usez)
                        (*result)(i,2) = (*mtx)(i,2);
                }
                outputs.push_back(result);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_MATRIX);
        }

        return true;
    }
    else if (nargin > 1)
    {
        bool ignorez = nargin > 2 && inputs[2].IsEmpty();

        if (nargout > 1)
        {
            if (ignorez)
                outputs = mtxFun(eval, std::vector<Currency>(inputs.cbegin(), inputs.cend() - 1), (int)(nargin - 1), &cart2polHelper);
            else
                outputs = mtxFun(eval, inputs, (int)nargin, &cart2polHelper);

            if (ignorez)
                outputs.push_back(EvaluatorInterface::allocateMatrix());
        }
        else
        {
            bool usez = nargin > 2 && !ignorez;
            std::vector<Currency> results;
            if (ignorez)
                results = mtxFun(eval, std::vector<Currency>(inputs.cbegin(), inputs.cend() - 1), (int)(nargin - 1), &cart2polHelper);
            else
                results = mtxFun(eval, inputs, (int)nargin, &cart2polHelper);

            const hwMatrix* x = results[0].Matrix();
            const hwMatrix* y = results[1].Matrix();
            const hwMatrix* z;
            if (usez)
                z = results[2].Matrix();
            hwMatrix* outmtx = EvaluatorInterface::allocateMatrix(x->Size(), ignorez ? (int)(nargin - 1) : (int)nargin, x->Type());

            if (x->IsReal())
            {
                for (int i = 0; i < outmtx->M(); ++i)
                {
                    (*outmtx)(i,0) = (*x)(i);
                    (*outmtx)(i,1) = (*y)(i);
                    if (usez)
                        (*outmtx)(i,2) = (*z)(i);
                }
            }
            else
            {
                throw OML_Error(OML_ERR_NUMARGIN);
            }

            outputs.push_back(outmtx);
        }
    }
    else
        throw OML_Error(OML_ERR_NUMARGIN);

    return true;
}
//------------------------------------------------------------------------------
// Returns true if given input is a vector [isvector]
//------------------------------------------------------------------------------
bool oml_isvector(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& in1 = inputs[0];
    if (in1.IsScalar() || in1.IsComplex() || in1.IsFunctionHandle())
    {
        outputs.push_back(getTrue());
    }
    else if (in1.IsMatrix() || in1.IsString())
    {
        outputs.push_back(HW_MAKE_BOOL_CURRENCY(in1.Matrix()->IsVector()));
    }
    else if (in1.CellArray())
    {
        outputs.push_back(HW_MAKE_BOOL_CURRENCY(in1.CellArray()->IsVector()));
    }
    else if (in1.IsStruct())
    {
        StructData* sd = in1.Struct();
        outputs.push_back(HW_MAKE_BOOL_CURRENCY(sd->M() == 1 || sd->N() == 1));
    }
    else
    {
        outputs.push_back(getFalse());
    }
    return true;
}
//------------------------------------------------------------------------------
// Trims trailing whitespaces from strings, matrices and cell arrays [deblank]
//------------------------------------------------------------------------------
bool oml_deblank(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(deblankHelper(eval, inputs[0]));
    return true;
}
//------------------------------------------------------------------------------
// Converts polar or cylindrical coordinates to Cartesian [pol2cart]
//------------------------------------------------------------------------------
bool oml_pol2cart(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    int nargout = getNumOutputs(eval);

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin == 1)
    {
        if (!inputs[0].IsMatrix())
            throw OML_Error(HW_ERROR_1INPMUSTMAT);

        const hwMatrix* mtx = inputs[0].Matrix();
        if (mtx->N() != 2 && mtx->N() != 3)
            throw OML_Error(HW_ERROR_MATMUST2OR3COL);

        bool usez = mtx->N() > 2;
        if (mtx->IsReal())
        {
            if (nargout > 1)
            {
                hwMatrix* x = EvaluatorInterface::allocateMatrix(mtx->M(), 1, hwMatrix::REAL);
                hwMatrix* y = EvaluatorInterface::allocateMatrix(mtx->M(), 1, hwMatrix::REAL);
                hwMatrix* z;
                if (usez)
                    z = EvaluatorInterface::allocateMatrix(mtx->M(), 1, hwMatrix::REAL);

                for (int i = 0; i < mtx->M(); ++i)
                {
                    std::pair<double,double> singleResult = pol2cartHelper((*mtx)(i,0), (*mtx)(i,1));

                    (*x)(i) = singleResult.first;
                    (*y)(i) = singleResult.second;
                    if (usez)
                        (*z)(i) = (*mtx)(i,2);
                }
                outputs.push_back(x);
                outputs.push_back(y);
                if (usez)
                    outputs.push_back(z);
                else
                    outputs.push_back(EvaluatorInterface::allocateMatrix());
            }
            else
            {
                hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

                for (int i = 0; i < mtx->M(); ++i)
                {
                    std::pair<double,double> singleResult = pol2cartHelper((*mtx)(i,0), (*mtx)(i,1));

                    (*result)(i,0) = singleResult.first;
                    (*result)(i,1) = singleResult.second;
                    if (usez)
                        (*result)(i,2) = (*mtx)(i,2);
                }
                outputs.push_back(result);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_MATRIX);
        }

        return true;
    }
    else if (nargin > 1)
    {
        bool ignorez = nargin > 2 && inputs[2].IsEmpty();

        if (inputs[0].IsMatrix())
        {
            if (!inputs[0].Matrix()->IsRealData())
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);
        }

        if (inputs[1].IsMatrix())
        {
            if (!inputs[1].Matrix()->IsRealData())
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
        }

        if (nargin > 2 && inputs[2].IsMatrix())
        {
            if (!inputs[2].Matrix()->IsRealData())
                throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);
        }

        if (nargout > 1)
        {
            if (ignorez)
                outputs = mtxFun(eval, std::vector<Currency>(inputs.cbegin(), inputs.cend() - 1), (int)(nargin - 1), &pol2cartHelper);
            else
                outputs = mtxFun(eval, inputs, (int)nargin, &pol2cartHelper);

            if (ignorez)
                outputs.push_back(EvaluatorInterface::allocateMatrix());
        }
        else
        {
            bool usez = nargin > 2 && !ignorez;
            std::vector<Currency> results;
            if (ignorez)
                results = mtxFun(eval, std::vector<Currency>(inputs.cbegin(), inputs.cend() - 1), (int)(nargin - 1), &pol2cartHelper);
            else
                results = mtxFun(eval, inputs, (int)nargin, &pol2cartHelper);

            const hwMatrix* x = results[0].Matrix();
            const hwMatrix* y = results[1].Matrix();
            const hwMatrix* z;

            if (usez)
                z = results[2].Matrix();

            hwMatrix* outmtx = EvaluatorInterface::allocateMatrix((int)x->Size(), ignorez ? (int)(nargin - 1) : (int)nargin, x->Type());

            if (!x->IsReal())
               throw OML_Error(OML_ERR_REAL, 1, OML_VAR_MATRIX);

            if (!y->IsReal())
               throw OML_Error(OML_ERR_REAL, 2, OML_VAR_MATRIX);

            for (int i = 0; i < outmtx->M(); ++i)
            {
                (*outmtx)(i,0) = (*x)(i);
                (*outmtx)(i,1) = (*y)(i);
                if (usez)
                    (*outmtx)(i,2) = (*z)(i);
            }

            outputs.push_back(outmtx);
        }
    }
    else
        throw OML_Error(OML_ERR_NUMARGIN);

    return true;
}
//------------------------------------------------------------------------------
// Gets number of inputs passed to a function [nargin]
//------------------------------------------------------------------------------
bool oml_nargin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() == 0) // inside a function
    {
        int nargin = eval.GetNarginValue();

        if (nargin == -1)
            throw OML_Error(HW_ERROR_INVCONTNARGIN);

        outputs.push_back(nargin);
        return true;
    }
    else if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    FunctionInfo* fi = nullptr;
    if (inputs[0].IsString())
    {
        FUNCPTR fptr = nullptr;
        std::string fname = readString(inputs[0]);
        if (!eval.FindFunctionByName(fname, &fi, &fptr))
            throw OML_Error(HW_ERROR_INVFUNCNAME);

        if (fptr)
        {
            if (eval.HasBuiltin(fname))
            {
                int n = eval.NarginFor(fname);
                if (n != OML_NO_NARG)
                {
                    outputs.push_back(static_cast<double>(n));
                    return true;
                }
            }
			throw OML_Error(HW_ERROR_INVFUNCNARGIN);
        }
    }
    else if (inputs[0].IsFunctionHandle())
    {
        fi = inputs[0].FunctionHandle();

		if (fi && fi->IsBuiltIn())
		{
            if (eval.HasBuiltin(fi->FunctionName()))
            {
                int n = eval.NarginFor(fi->FunctionName());
                if (n != OML_NO_NARG)
                {
                    outputs.push_back(static_cast<double>(n));
                    return true;
                }
            }
            throw OML_Error(HW_ERROR_INVFUNCNARGIN);
		}
    }
    else
	{
        throw OML_Error(HW_ERROR_INPUTSTRINGFUNC);
	}

    outputs.push_back(static_cast<double>(fi->Nargin()));
    return true;
}
//------------------------------------------------------------------------------
// Gets number of inputs returned from a function [nargout]
//------------------------------------------------------------------------------
bool oml_nargout(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() == 0)
    {
        int nargout = eval.GetNargoutValue();

        if (nargout == -1)
            throw OML_Error(HW_ERROR_INVCONTNARGOUT);

        outputs.push_back(nargout);
        return true;
    }
    else if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    FunctionInfo* fi = nullptr;
    if (inputs[0].IsString())
    {
        FUNCPTR fptr = nullptr;
        std::string fname = readString(inputs[0]);
        if (!eval.FindFunctionByName(fname, &fi, &fptr))
            throw OML_Error(HW_ERROR_INVFUNCNAME);

        if (fptr)
        {
            if (eval.HasBuiltin(fname))
            {
                int n = eval.NargoutFor(fname);
                if (n != OML_NO_NARG)
                {
                    outputs.push_back(static_cast<double>(n));
                    return true;
                }
            }
            throw OML_Error(HW_ERROR_INVFUNCCALLNARGOUT);
        }
    }
    else if (inputs[0].IsFunctionHandle())
    {
        fi = inputs[0].FunctionHandle();
    }
    else
        throw OML_Error(HW_ERROR_FUNCHANDLNAMEINPUT);

    outputs.push_back(static_cast<double>(fi->Nargout()));
    return true;
}
//------------------------------------------------------------------------------
// Performs the specified bitwise operation
//------------------------------------------------------------------------------
static bool bitOperatorFunc(EvaluatorInterface&          eval, 
                            const std::vector<Currency>& inputs, 
                            std::vector<Currency>&       outputs, 
                            int                          (*bitop)(int, int))
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_NUMERIC, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_NUMERIC, 2, OML_VAR_DATA);

    const hwMatrix* x = inputs[0].ConvertToMatrix();
    const hwMatrix* y = inputs[1].ConvertToMatrix();

    if (!x->IsRealData())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

    if (!y->IsRealData())
        throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

    int m = 0;
    int n = 0;
    getMatchingSizes(x, y, &m, &n);

    if (m == -1)
    {
		if (IsNaN_T((*x)(0)) || IsNaN_T((*y)(0)))
		{
			outputs.push_back(0.0);
		}
		else
		{
			outputs.push_back(static_cast<double>((*bitop)(
				(int)realval(x, 0), (int)realval(y, 0))));
		}
    }
    else
    {
        hwMatrix* result = eval.allocateMatrix(m, n, hwMatrix::REAL);
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                double xx = realvalorscalar(x,i,j);
                double yy = realvalorscalar(y,i,j);
                if (isnonnegint(xx) && isnonnegint(yy))
                    (*result)(i,j) = static_cast<double>((*bitop)((int)xx,(int)yy));
                else
                    (*result)(i,j) = 0.0;  // Handle negative integers
            }
        }
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Performs bitwise AND operation [bitand]
//------------------------------------------------------------------------------
bool oml_bitand(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return bitOperatorFunc(eval, inputs, outputs, &bitAnd);
}
//------------------------------------------------------------------------------
// Performs bitwise OR operation [bitor]
//------------------------------------------------------------------------------
bool oml_bitor(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return bitOperatorFunc(eval, inputs, outputs, &bitOr);
}
//------------------------------------------------------------------------------
// Performs bitwise XOR operation [bitxor]
//------------------------------------------------------------------------------
bool oml_bitxor(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return bitOperatorFunc(eval, inputs, outputs, &bitXor);
}
//------------------------------------------------------------------------------
// Returns true if input is a directory [isdir]
//------------------------------------------------------------------------------
bool oml_isdir(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::string path = readString(inputs[0]);
    outputs.push_back(HW_MAKE_BOOL_CURRENCY(isDirectory(path)));
    return true;
}
//------------------------------------------------------------------------------
// Removes given field in a structure [rmfield]
//------------------------------------------------------------------------------
bool oml_rmfield(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsStruct())
        throw OML_Error(OML_ERR_STRUCT, 1, OML_VAR_TYPE);
    
    StructData* strct = eval.allocateStruct(inputs[0].Struct());
    Currency out(strct);
    
    if (inputs[1].IsString())
    {
        removeFields(eval, strct, inputs[1].Matrix());
    }
    else if (inputs[1].IsCellArray())
    {
        HML_CELLARRAY* cell = inputs[1].CellArray();
        if (!isstr(cell))
            throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

        for (int i = 0; i < cell->Size(); ++i)
            removeFields(eval, strct, (*cell)(i).Matrix());
    }
    else
        throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Returns true if running on windows [ispc]
//------------------------------------------------------------------------------
bool oml_ispc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
#ifdef OS_WIN
    outputs.push_back(getTrue());
#else
    outputs.push_back(getFalse());
#endif
    return true;
}
//------------------------------------------------------------------------------
// Returns true if running on OSX [ismac]
//------------------------------------------------------------------------------
bool oml_ismac(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
#if defined(__APPLE__) && defined(__MACH__) // true for OSX, false otherwise
    outputs.push_back(getTrue());
#else
    outputs.push_back(getFalse());
#endif
    return true;
}
//------------------------------------------------------------------------------
// Returns true if running on Linux [isunix]
//------------------------------------------------------------------------------
bool oml_isunix(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
#if !defined(OS_WIN) && (!defined(__APPLE__) || defined(__MACH__)) // if apple, must be OSX
    outputs.push_back(getTrue());
#else
    outputs.push_back(getFalse());
#endif
    return true;
}
//------------------------------------------------------------------------------
// Assigns value to the given variable in the given context [assignin]
//------------------------------------------------------------------------------
bool oml_assignin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsString())
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_CONTEXT);

    if (!inputs[1].IsString())
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_VARIABLE);

    std::string context = readOption(eval, inputs[0]);
    std::string varname = readString(inputs[1]);
    
    bool base_context;
    if (context == "base")
        base_context = true;
    else if (context == "caller")
        base_context = false;
    else
        throw OML_Error(HW_ERROR_CONTBASEORCALL);
    
    if (!isvarname(varname))
        throw OML_Error("Error: invalid variable name: " + varname);

    EvaluatorInterface context_eval(eval, base_context);
    context_eval.SetValue(varname, inputs[2]);
    return true;
}
//------------------------------------------------------------------------------
// Evaluates variable in context [evalin]
//------------------------------------------------------------------------------
bool oml_evalin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsString())
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

    if (!inputs[1].IsString())
        throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);

    if (nargin > 2 && !inputs[2].IsString())
        throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);

    if (eval.isUsedForEvalin())
        throw OML_Error(HW_ERROR_NOTCALLRECURSEEVALIN);

    std::string context = readOption(eval, inputs[0]);
    std::string tryeval = readString(inputs[1]);
    std::string catcheval;
    if (nargin > 2)
        catcheval = readString(inputs[2]);

    std::vector<Currency> new_inputs(inputs.begin()+1, inputs.end());

    bool base_context;
    if (context == "base")
        base_context = true;
    else if (context == "caller")
        base_context = false;
    else
        throw OML_Error(HW_ERROR_CONTBASEORCALL);

    EvaluatorInterface new_eval(eval, base_context, false);
    oml_eval(new_eval, new_inputs, outputs);
    eval.transferResultsFrom(new_eval);
    return true;
}
//------------------------------------------------------------------------------
// Creates a struct [struct]
//------------------------------------------------------------------------------
bool oml_struct(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin == 1 && inputs[0].IsEmpty())
    {
        StructData*     sd  = EvaluatorInterface::allocateStruct();
		const hwMatrix* mtx = inputs[0].Matrix();

        sd->DimensionNew(mtx->M(), mtx->N()); // retain the size of the empty object passed in
        outputs.push_back(sd);
        return true;
    }
    else if (nargin % 2)
        throw OML_Error(HW_ERROR_EVERYVALFIELD);

    std::unique_ptr<StructData> strct(EvaluatorInterface::allocateStruct());

    if (nargin)
    {
        int m = -1;
        int n = -1;
        bool checkingField = true;
        std::vector<Currency>::const_iterator iter = inputs.cbegin();

        while (iter != inputs.cend())
        {
            if (checkingField)
            {
                if (!iter->IsString())
                    throw OML_Error(HW_ERROR_FIELDNAMESTR);
            }
            else
            {
                if (iter->IsCellArray())
                {
                    HML_CELLARRAY* cell = iter->CellArray();
                    if (cell->Size() != 1)
                    {
                        if (m == -1)
                        {
                            m = cell->M();
                            n = cell->N();
                        }
                        else if (cell->M() != m || cell->N() != n)
                        {
                            throw OML_Error(HW_ERROR_CELLARRAYSSAMESIZE);
                        }
                    }
                }
            }

            checkingField = !checkingField;
            ++iter;
        }

        if (m == -1)
            strct->DimensionNew(1,1); // 1x1
        else
            strct->DimensionNew(m, n);

        iter = inputs.cbegin();
        if (m == 0 || n == 0)
        {
            // add fields to empty struct
            for (; iter != inputs.cend(); iter += 2)
            {
                std::string field = readString(*iter);
                strct->addField(field);
            }
        }
        else
        {
            while (iter != inputs.cend())
            {
                std::string field = readString(*iter++);
                if (iter->IsCellArray())
                {
                    HML_CELLARRAY* cell = iter->CellArray();
                    if (m == -1)
                    {
                        strct->SetValue(0, 0, field, (*cell)(0));
                    }
                    else if (cell->Size() == 1)
                    {
                        const Currency& elem = (*cell)(0);
                        for (int i = 0; i < m; ++i)
                        {
                            for (int j = 0; j < n; ++j)
                            {
                                strct->SetValue(i, j, field, elem);
                            }
                        }
                    }
                    else
                    {
                        for (int i = 0; i < m; ++i)
                        {
                            for (int j = 0; j < n; ++j)
                            {
                                strct->SetValue(i, j, field, (*cell)(i,j));
                            }
                        }
                    }
                }
                else
                {
                    if (m == -1)
                    {
                        strct->SetValue(0, 0, field, *iter);
                    }
                    else
                    {
                        for (int i = 0; i < m; ++i)
                        {
                            for (int j = 0; j < n; ++j)
                            {
                                strct->SetValue(i, j, field, *iter);
                            }
                        }
                    }
                }
                ++iter;
            }
        }
    }
    else
        strct->DimensionNew(1,1);

    outputs.push_back(strct.release());
    return true;
}
//------------------------------------------------------------------------------
// Returns true if given input is a cell [iscell]
//------------------------------------------------------------------------------
bool oml_iscell(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(HW_MAKE_BOOL_CURRENCY(inputs[0].IsCellArray()));
    return true;
}
//------------------------------------------------------------------------------
// Returns true if given input is finite [isfinite]
//------------------------------------------------------------------------------
bool oml_isfinite(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return allowComplexToLogical(inputs, outputs, &::isfinite, true);
}
//------------------------------------------------------------------------------
// Unary plus sign [uplus]
//------------------------------------------------------------------------------
bool oml_uplus(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(inputs[0]);
    return true;
}
//------------------------------------------------------------------------------
// Performs logical negation [not]
//------------------------------------------------------------------------------
bool oml_not(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return unOperatorFunc(eval, inputs, outputs, NEGATE);
}
//------------------------------------------------------------------------------
// Unary minus sign [uminus]
//------------------------------------------------------------------------------
bool oml_uminus(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return unOperatorFunc(eval, inputs, outputs, UMINUS);
}
//------------------------------------------------------------------------------
// Perform addition on a list of arguments from left to right [plus]
//------------------------------------------------------------------------------
bool oml_plus(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return binOperatorFuncVararg(eval, inputs, outputs, PLUS);
}
//------------------------------------------------------------------------------
// Perform multiplication on a list of arguments from left to right [times]
//------------------------------------------------------------------------------
bool oml_times(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return binOperatorFuncVararg(eval, inputs, outputs, ETIMES);
}
//------------------------------------------------------------------------------
// Performs matrix multiplication [mtimes]
//------------------------------------------------------------------------------
bool oml_mtimes(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return binOperatorFuncVararg(eval, inputs, outputs, TIMES);
}
//------------------------------------------------------------------------------
// Performs subtraction [minus]
//------------------------------------------------------------------------------
bool oml_minus(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return binOperatorFunc(eval, inputs, outputs, MINUS);
}
//------------------------------------------------------------------------------
// Perform element-wise matrix right division [rdivide]
//------------------------------------------------------------------------------
bool oml_rdivide(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return binOperatorFunc(eval, inputs, outputs, EDIV);
}
//------------------------------------------------------------------------------
// Perform element-wise matrix left division [ldivide]
//------------------------------------------------------------------------------
bool oml_ldivide(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return binOperatorFunc(eval, inputs, outputs, ELDIV);
}
//------------------------------------------------------------------------------
// Perform left division of scalars and/or matrices [mldivide]
//------------------------------------------------------------------------------
bool oml_mldivide(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return binOperatorFunc(eval, inputs, outputs, LDIV);
}
//------------------------------------------------------------------------------
// Perform right division of scalars and/or matrices [mrdivide]
//------------------------------------------------------------------------------
bool oml_mrdivide(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return binOperatorFunc(eval, inputs, outputs, DIV);
}
//------------------------------------------------------------------------------
// Performs power operation [power]
//------------------------------------------------------------------------------
bool oml_power(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return binOperatorFunc(eval, inputs, outputs, DOTPOW);
}
//------------------------------------------------------------------------------
// Performs matrix power operation [mpower]
//------------------------------------------------------------------------------
bool oml_mpower(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return binOperatorFunc(eval, inputs, outputs, POW);
}
//------------------------------------------------------------------------------
// Performs logical conjunction equivalent to the & operator [and]
//------------------------------------------------------------------------------
bool oml_and(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return logOperatorFunc(eval, inputs, outputs, AND);
}
//------------------------------------------------------------------------------
// Performs logical disjunction, equivalent to the | operator [or]
//------------------------------------------------------------------------------
bool oml_or(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return logOperatorFunc(eval, inputs, outputs, OR);
}
//------------------------------------------------------------------------------
// Removes the given directory [rmdir]
//------------------------------------------------------------------------------
bool oml_rmdir(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::string filename = readString(toCurrencyStr(eval, inputs[0], false, false));
    filename = BuiltInFuncsUtils::GetAbsolutePath(filename);
#ifdef OS_WIN
    BOOL result = RemoveDirectory((LPCSTR) filename.c_str());
    if (!result)
    {
        DWORD err = GetLastError();
        outputs.push_back(getFalse());
        if (err == ERROR_PATH_NOT_FOUND)
            outputs.push_back(std::string("No such file or directory"));
        else
            outputs.push_back(std::string("Problem removing directory"));
        outputs.push_back(std::string("rmdir"));
    }
#else
    int result = rmdir(filename.c_str());
    if (result)
    {
        outputs.push_back(getFalse());
        outputs.push_back(std::string(strerror(errno)));
        outputs.push_back(std::string("rmdir"));
    }
#endif
    else
    {
        outputs.push_back(getTrue());
        outputs.push_back(std::string());
        outputs.push_back(std::string());

		eval.OnRefreshDirs();
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns the current date as a string in the format DD-MMM-YYYY [date]
//------------------------------------------------------------------------------
bool oml_date(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    time_t t = time(0);
    struct tm *currentTime = localtime(&t);
    char day[4];
    sprintf(day, "%.2d-", (int) currentTime->tm_mday);
    std::string date(day);

    switch (currentTime->tm_mon)
    {
        case 0:
            date += "Jan-";
            break;
        case 1:
            date += "Feb-";
            break;
        case 2:
            date += "Mar-";
            break;
        case 3:
            date += "Apr-";
            break;
        case 4:
            date += "May-";
            break;
        case 5:
            date += "Jun-";
            break;
        case 6:
            date += "Jul-";
            break;
        case 7:
            date += "Aug-";
            break;
        case 8:
            date += "Sep-";
            break;
        case 9:
            date += "Oct-";
            break;
        case 10:
            date += "Nov-";
            break;
        case 11:
            date += "Dec-";
            break;
    }

    date += std::to_string(static_cast<long long>(currentTime->tm_year + 1900));
    outputs.push_back(date);
    return true;
}
//------------------------------------------------------------------------------
// Resets search path for functions [restoredefaultpath]
//------------------------------------------------------------------------------
bool oml_restorepath(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    eval.RestorePath();
    if (getNumOutputs(eval))
        outputs.push_back(getPathString(eval, pathsep));

    return true;
}
//------------------------------------------------------------------------------
// Performs less than comparison [lt]
//------------------------------------------------------------------------------
bool oml_lt(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return ineqOperatorFunc(eval, inputs, outputs, LTHAN);
}
//------------------------------------------------------------------------------
// Performs greater than comparison [gt]
//------------------------------------------------------------------------------
bool oml_gt(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return ineqOperatorFunc(eval, inputs, outputs, GTHAN);
}
//------------------------------------------------------------------------------
// Performs greater than or equal comparison [ge]
//------------------------------------------------------------------------------
bool oml_geq(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return ineqOperatorFunc(eval, inputs, outputs, GEQ);
}
//------------------------------------------------------------------------------
// Performs less than or equal comparison [le]
//------------------------------------------------------------------------------
bool oml_leq(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return ineqOperatorFunc(eval, inputs, outputs, LEQ);
}
//------------------------------------------------------------------------------
// Performs equality comparison [eq]
//------------------------------------------------------------------------------
bool oml_eq(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return eqOperatorFunc(eval, inputs, outputs, EQUAL);
}
//------------------------------------------------------------------------------
// Performs inequality comparison [ne]
//------------------------------------------------------------------------------
bool oml_ne(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return eqOperatorFunc(eval, inputs, outputs, NEQUAL);
}
#if 0 // Unused code
bool oml_qrsolve(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return twoMatrixCaller(eval, inputs, outputs, &hwMatrix::QRSolve);
}
#endif
//------------------------------------------------------------------------------
// Solves a linear system of equations [linsolve]
//------------------------------------------------------------------------------
bool oml_linsolve(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin == 2)
    {
        return twoMatrixCaller(eval, inputs, outputs, &hwMatrix::LSolve);
    }
    else // nargin == 3
    {
        if (!inputs[2].IsStruct())
            throw OML_Error(OML_ERR_OPTION, 3, OML_VAR_STRUCT);

        std::vector<Currency> newinputs(inputs.cbegin(), inputs.cbegin() + 2);
        const StructData* sd = inputs[2].Struct();
        bool lt, ut, uhess, sym, posdef, rect, transa;
        lt = ut = uhess = sym = posdef = rect = transa = false;

        if (sd->M() != 1 || sd->N() != 1)
            throw OML_Error(HW_ERROR_OPTSTRUCTSIZE1);

        const std::map<std::string, int>& fields = sd->GetFieldNames();
        std::map<std::string, int>::const_iterator iter;
        for (iter = fields.cbegin(); iter != fields.cend(); ++iter)
        {
            std::string op = iter->first;
            bool val = boolFromCur(sd->GetValue(0, 0, op));
            if (op == "LT")
                lt = val;
            else if (op == "UT")
                ut = val;
            else if (op == "UHESS")
                uhess = val;
            else if (op == "SYM")
                sym = val;
            else if (op == "POSDEF")
                posdef = val;
            else if (op == "RECT")
                rect = val;
            else if (op == "TRANSA")
                transa = val;
            else
                BuiltInFuncsUtils::SetWarning(eval, "Warning: ignoring unknown option " + op);
        }

        if (posdef)
            return twoMatrixCaller(eval, newinputs, outputs, &hwMatrix::LSolveSPD);
        else
            return twoMatrixCaller(eval, newinputs, outputs, &hwMatrix::LSolve);
    }
}
//------------------------------------------------------------------------------
// Computes the matrix pseudo-inverse [pinv]
//------------------------------------------------------------------------------
bool oml_pinv(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return oneMatrixCaller(eval, inputs, outputs, &hwMatrix::Pinv);
}
//------------------------------------------------------------------------------
// Normalizes vectors [normalize]
//------------------------------------------------------------------------------
bool oml_normalize(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    hwMatrix* mtx = matrixCopyFromInput(inputs[0], false);

    if (!mtx)
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    Currency mtxcur(mtx);
    if (mtx->IsVector())
    {
        if (mtx->IsRealData())
        {
            if (mtx->IsReal())
                BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Normalize());
            else
            {
                std::unique_ptr<hwMatrix> mtxcopy(makeRealCopy(mtx));
                BuiltInFuncsUtils::CheckMathStatus(eval, mtxcopy->Normalize());
            }
        }
        else
            throw OML_Error(HW_MATH_MSG_COMPLEX);
    }
    else
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);
    outputs.push_back(mtxcur);
    return true;
}
//------------------------------------------------------------------------------
// Computes the matrix exponential [expm]
//------------------------------------------------------------------------------
bool oml_expm(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& in1 = inputs[0];
    if (in1.IsScalar() || in1.IsComplex())
    {
        return oml_exp(eval, inputs, outputs);
    }
    else if (in1.IsMatrix() || in1.IsString())
    {
        const hwMatrix* mtx = in1.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix();
        Currency out(result);
        BuiltInFuncsUtils::CheckMathStatus(eval, result->MatExp(*mtx));
        outputs.push_back(out);
        return true;
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
}
//------------------------------------------------------------------------------
// Returns true after making directories [mkdir]
//------------------------------------------------------------------------------
bool oml_mkdir(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.empty() ? 0: inputs.size();
    if (nargin < 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    Currency cur1 = inputs[0];
    if (!cur1.IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }

    std::string parentdir;
    std::string newdir;

    if (nargin == 1)
    {
        newdir = cur1.StringVal();
    }
    else 
    {
        Currency cur2 = inputs[1];
        if (!cur2.IsString())
        {
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
        }
        newdir    = cur2.StringVal();
        parentdir = cur1.StringVal();
    }

#ifdef OS_WIN
    std::string fullpath (newdir);
	if (!parentdir.empty())
	{
		BuiltInFuncsUtils::StripTrailingSlash(parentdir);
        fullpath = parentdir + "/" + newdir;
	}
    fullpath = BuiltInFuncsUtils::Normpath(fullpath);

    BOOL result = CreateDirectory((LPCSTR) fullpath.c_str(), nullptr);
    if (!result)
    {
        DWORD err = GetLastError();
        outputs.push_back(getFalse());
        if (err == ERROR_ALREADY_EXISTS)
            outputs.push_back(std::string("directory exists"));
        else if (err == ERROR_PATH_NOT_FOUND)
            outputs.push_back(std::string("No such file or directory"));
        else
            outputs.push_back(std::string("Could not create directory"));
        outputs.push_back(std::string("mkdir"));
    }
#else
    if (parentdir.empty() && !BuiltInFuncsUtils::IsAbsolutePath(newdir))
    {
        parentdir = BuiltInFuncsUtils::GetCurrentWorkingDir();
    }
    else if (!parentdir.empty())
    {
        parentdir = BuiltInFuncsUtils::GetAbsolutePath(parentdir);
    }

    std::string fullpath = (!parentdir.empty()) ?
        BuiltInFuncsUtils::Normpath(parentdir + "/" + newdir) : newdir;

    // Use the system command to make recursive directories
    std::string strcmd ("mkdir -p " + fullpath + " > /dev/null");
    int result = system(strcmd.c_str());
    if (result != 0)
    {
        std::string err (strerror(errno));
        std::string msg = (!err.empty() && err != "No error") ? err : "";

        outputs.push_back(false);
        outputs.push_back(msg);
        outputs.push_back("mkdir");
    }
#endif
    else
    {
        outputs.push_back(true);
        outputs.push_back(std::string());
        outputs.push_back(std::string());

        eval.OnRefreshDirs();
    }

    return true;
}
//------------------------------------------------------------------------------
// Sets a local environment variable [setenv]
//------------------------------------------------------------------------------
bool oml_setenv(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::string varname = readString(toCurrencyStr(eval, inputs[0], false, false));
    std::string value   = readString(toCurrencyStr(eval, inputs[1], false, false));
#ifdef OS_WIN
    if (_putenv_s(varname.c_str(), value.c_str()))
#else
    if (setenv(varname.c_str(), value.c_str(), 1))
#endif
        throw OML_Error(HW_ERROR_ENVVARFAIL);
    return true;
}
//------------------------------------------------------------------------------
// Computes the condition number of a matrix [cond]
//------------------------------------------------------------------------------
bool oml_cond(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsLogical())
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_VARIABLE);

    double cond;
    const hwMatrix* mtx = inputs[0].ConvertToMatrix();
    BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Cond(cond));
    outputs.push_back(cond);
    return true;
}
//------------------------------------------------------------------------------
// Sorts elements of the given input in ascending order by default [sort]
//------------------------------------------------------------------------------
bool oml_sort(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    bool ascend = true;
    int dim = 0;

    if (nargin > 1)
    {
        bool modeset = false;
        bool dimset  = false;

        for (int i = 1; i < nargin; ++i)
        {
            const Currency& cur = inputs[i];
            if (cur.IsString())
            {
                if (modeset)
                    throw OML_Error(HW_ERROR_SETMODEMOREONCE);
                std::string str = readOption(eval, cur);
                if (str == "descend")
                    ascend = false;
                else if (str != "ascend")
                    throw OML_Error(HW_ERROR_INVALIDOPTION(str));
                modeset = true;
            }
            else
            {
                if (dimset)
                    throw OML_Error(HW_ERROR_SETDIMMOREONCE);

                if (!cur.IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, i+1, OML_VAR_DIM);

                dim = (int) cur.Scalar();
                dimset = true;
            }
        }
    }

    int temp = 0;
    const Currency& in1 = inputs[0];
    if (in1.IsCellArray())
    {
        if (!isstr(in1.CellArray()))
        {
            throw OML_Error(HW_ERROR_CELLELEMSTR);
        }

        HML_CELLARRAY* cell = in1.CellArray();

        if (!dim)
        {
            dim = cell->M() == 1 ? 2 : 1;
        }
        else if (dim > 2)
        {
            outputs.push_back(in1);
            return true;
        }

        std::unique_ptr<hwMatrix> indices(EvaluatorInterface::allocateMatrix(cell->M(), cell->N(), hwMatrix::REAL));
        std::pair<int&, hwMatrix*> index_data(temp, indices.get());

        Currency out;
        if (ascend)
            out = oml_Matrix_sort(eval, cell, dim, &sort<true>, &index_data);
        else
            out = oml_Matrix_sort(eval, cell, dim, &sort<false>, &index_data);
        cell = out.CellArray();
        for (int i = 0; i < cell->Size(); ++i)
            (*cell)(i).SetMask(Currency::MASK_STRING);
        outputs.push_back(out);
        outputs.push_back(indices.release());
    }
    else if (in1.IsMatrix() || in1.IsScalar() || in1.IsComplex() || in1.IsString())
    {
        const hwMatrix* mtx = in1.ConvertToMatrix();

        if (!dim)
        {
            dim = mtx->M() == 1 ? 2 : 1;
        }
        else if (dim > 2)
        {
            outputs.push_back(in1);
            return true;
        }
        
        std::unique_ptr<hwMatrix> indices(EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL));
        std::pair<int&, hwMatrix*> index_data(temp, indices.get());

        Currency out;
        if (ascend)
            out = oml_Matrix_sort(eval, mtx, dim, &sort<true>, &index_data);
        else
            out = oml_Matrix_sort(eval, mtx, dim, &sort<false>, &index_data);

        out.SetMask(in1.GetMask());
        outputs.push_back(out);
        outputs.push_back(indices.release());
    }
    else if (in1.IsNDMatrix())
    {
        if (nargin == 1)
        {
            return oml_MatrixNUtil4(eval, inputs, outputs, oml_sort);
        }
        else if (nargin == 2)
        {
            if (inputs[1].IsString())       // sort(ND, mode)
                return oml_MatrixNUtil4(eval, inputs, outputs, oml_sort);
            else                            // sort(ND, dim)
                return oml_MatrixNUtil4(eval, inputs, outputs, oml_sort, 2);
        }
        else if (nargin == 3)
        {
            // sort(ND, dim, mode)
            return oml_MatrixNUtil4(eval, inputs, outputs, oml_sort, 2);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is numeric [isnumeric]
//------------------------------------------------------------------------------
bool oml_isnumeric(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& in1 = inputs[0];
    outputs.push_back(HW_MAKE_BOOL_CURRENCY(!in1.IsLogical() && (in1.IsScalar() || in1.IsComplex() ||
                                                                 in1.IsMatrix() || in1.IsNDMatrix())));
    return true;
}
//------------------------------------------------------------------------------
// Returns true if all elements in input are nonzero [all]
//------------------------------------------------------------------------------
bool oml_all(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin == 1 && inputs[0].IsMatrix())
	{
		const hwMatrix* matrix = inputs[0].Matrix();

        if (matrix->IsVector())
        {
            double value = all(eval, matrix);

            if (value)
            {
                outputs.push_back(true);
            }
            else
            {
                outputs.push_back(false);
            }

            return true;
        }
	}

    if (!anyall(eval, inputs, outputs, &all, true))
    {
        if (nargin == 1)
        {
            oml_MatrixNUtil3(eval, inputs, outputs, oml_all);
        }
        else
        {
            oml_MatrixNUtil3(eval, inputs, outputs, oml_all, 2);
        }

        outputs[0].SetMask(Currency::MASK_LOGICAL);
        return true;
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true if any element in input is nonzero [any]
//------------------------------------------------------------------------------
bool oml_any(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin == 1 && inputs[0].IsMatrix())
	{
		const hwMatrix* matrix = inputs[0].Matrix();

        if (matrix->IsVector())
        {
            double value = any(eval, matrix);

            if (value)
            {
                outputs.push_back(true);
            }
            else
            {
                outputs.push_back(false);
            }

            return true;
        }
	}

    if (!anyall(eval, inputs, outputs, &any, false))
    {
        if (nargin == 1)
        {
            oml_MatrixNUtil3(eval, inputs, outputs, oml_any);
        }
        else
        {
            oml_MatrixNUtil3(eval, inputs, outputs, oml_any, 2);
        }

        outputs[0].SetMask(Currency::MASK_LOGICAL);
        return true;
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input has a size of 1
//------------------------------------------------------------------------------
static bool is1x1(const Currency& cur)
{
    if (cur.IsScalar() || cur.IsComplex() || cur.IsFunctionHandle())
        return true;
    else if (cur.IsString() || cur.IsMatrix())
        return cur.Matrix()->Size() == 1;
    else if (cur.IsCellArray())
        return cur.CellArray()->Size() == 1;
    else if (cur.IsStruct())
        return cur.Struct()->M() * cur.Struct()->N() == 1;
    return false;
}
//------------------------------------------------------------------------------
// Returns true if input is scalar [isscalar]
//------------------------------------------------------------------------------
bool oml_isscalar(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(HW_MAKE_BOOL_CURRENCY(is1x1(inputs[0])));
    return true;
}
//------------------------------------------------------------------------------
// Helper method to throw error
//------------------------------------------------------------------------------
//#undef FormatMessage
void _OML_Error(EvaluatorInterface& eval, std::vector<Currency>::const_iterator start, std::vector<Currency>::const_iterator end)
{
    if (start != end)
    {
        if (!start->IsString())
        {
            throw OML_Error(HW_ERROR_TEMPLATESTR);
        }
        std::string msg = (start + 1 == end) ? (*start).StringVal() : 
                                               sprintf(eval, start, end);

        if (msg.empty())
        {
            return;
        }
        std::string lower (msg);
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        if (lower.find("error: ", 0, 7) == std::string::npos)
        {
            msg.insert(0, "Error: ");
        }
        throw OML_Error(msg);
    }
}
//------------------------------------------------------------------------------
// Throws an error with the given message/template [error]
//------------------------------------------------------------------------------
bool oml_error(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    int nargin = (inputs.empty()) ? 0 : static_cast<int>(inputs.size());
    if (nargin < 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    Currency cur = inputs[0];
    if (!cur.IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }
    std::string msg = (nargin == 1) ? cur.StringVal() : 
                      sprintf(eval, inputs.cbegin(), inputs.cend());
    if (msg.empty())
    {
        return true;
    }

    std::string lower (msg);
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower.find("error: ", 0, 7) == std::string::npos)
    {
        msg.insert(0, "Error: ");
    }

    throw OML_Error(msg);
    return true;
}
//------------------------------------------------------------------------------
// Sets warning [warning]
//------------------------------------------------------------------------------
bool oml_warning(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (!nargin)
        throw OML_Error(OML_ERR_NUMARGIN);
    
    if (!inputs[0].IsString())
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

    std::string msg = nargin > 1 ? sprintf(eval, inputs.cbegin(), inputs.cend()) : readString(inputs[0]);
    if (msg.length())
        BuiltInFuncsUtils::SetWarning(eval, "Warning: " + msg);
    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is infinite [isinf]
//------------------------------------------------------------------------------
bool oml_isinf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return allowComplexToLogical(inputs, outputs, &isinfinity, false);
}
//------------------------------------------------------------------------------
// Returns true if input is not a number [isnan]
//------------------------------------------------------------------------------
bool oml_isnan(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return allowComplexToLogical(inputs, outputs, &IsNaN_T<double>, false);
}
//------------------------------------------------------------------------------
// Returns true if input is a string [isstr]
//------------------------------------------------------------------------------
bool oml_isstr(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];
    outputs.push_back(HW_MAKE_BOOL_CURRENCY(cur.IsString()));
    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is a matrix [ismatrix]
//------------------------------------------------------------------------------
bool oml_ismatrix(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];
    outputs.push_back(HW_MAKE_BOOL_CURRENCY(cur.IsScalar() || cur.IsComplex() || cur.IsMatrix() || cur.IsString() || cur.IsNDMatrix()));
    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is a symmetric matrix [issymmetric]
//------------------------------------------------------------------------------
bool oml_issymmetric(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];

    if (input1.IsScalar())
    {
        outputs.push_back(true);
        return true;
    }

    if (!input1.IsMatrix())
    {
        outputs.push_back(false);
        return true;
    }

    const hwMatrix* matrix = input1.Matrix();

    if (!matrix->IsSquare())
    {
        outputs.push_back(false);
        return true;
    }

    if (nargin == 1)
    {
        outputs.push_back(matrix->IsSymmetric());
    }
    else // (nargin == 2)
    {
        const Currency& input2 = inputs[1];

        if (!input2.IsScalar())
            throw OML_Error(OML_ERR_SCALAR, 2, OML_VAR_PARAMETER);

        double tol =  input2.Scalar();
        double norm1, norm2;
        hwMatrix trans;
        hwMatrix diff;
        hwMathStatus status;

        status = trans.Transpose(*matrix);
        diff = *matrix - trans;

        status = diff.Norm(norm1, "inf");
        status = matrix->Norm(norm2, "inf");

        if (norm2 == 0 || norm1 <= norm2 * tol)
            outputs.push_back(true);
        else
            outputs.push_back(false);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is a square matrix [issquare]
//------------------------------------------------------------------------------
bool oml_issquare(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];

    if (input1.IsScalar() || input1.IsComplex())
    {
        outputs.push_back(true);
        return true;
    }

    if (!input1.IsMatrixOrString())
    {
        outputs.push_back(false);
        return true;
    }

    const hwMatrix* matrix = input1.Matrix();

    if (!matrix->IsSquare())
        outputs.push_back(false);
	else
        outputs.push_back(true);

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is a hermitian matrix [ishermitian]
//------------------------------------------------------------------------------
bool oml_ishermitian(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];

    if (!input1.IsMatrix())
    {
        outputs.push_back(false);
        return true;
    }

    const hwMatrix* matrix = input1.Matrix();

    if (matrix->IsReal())
        return oml_issymmetric(eval, inputs, outputs);

    if (!matrix->IsSquare())
    {
        outputs.push_back(false);
        return true;
    }

    if (nargin == 1)
    {
        outputs.push_back(matrix->IsHermitian());
    }
    else // (nargin == 2)
    {
        const Currency& input2 = inputs[1];

        if (!input2.IsScalar())
            throw OML_Error(OML_ERR_SCALAR, 2, OML_VAR_PARAMETER);

        double tol = input2.Scalar();
        double norm1, norm2;
        hwMatrix trans;
        hwMatrix herm;
        hwMatrix diff;
        hwMathStatus status;

        status = trans.Transpose(*matrix);
        herm.Conjugate(trans);
        diff = *matrix - herm;

        status = diff.Norm(norm1, "inf");
        status = matrix->Norm(norm2, "inf");

        if (norm2 == 0 || norm1 < norm2 * tol)
            outputs.push_back(true);
        else
            outputs.push_back(false);
    }

    return true;
}
//------------------------------------------------------------------------------
// Generates unique variable names, based on the given input(s) [genvarname]
//------------------------------------------------------------------------------
bool oml_genvarname(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    // If either input is a cell array, the output will be a cell array
    // Otherwise the output will be a string
    // Every string in the input must be 1-dimensional
    //
    // If a variable name is used, it will be appended with a number (e.g. x --> x1, x2, etc)
    // If the variable name is a number or keyword, it will be prepended with an underscore
    // If the variable name ends with a number, it will be appended with an underscore before
    // appending with a number
    // Every character the variable name contains that cannot be in a variable name is replaced
    // with an underscore

    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<std::string> varnameBases, varnames, exceptions, kwds = eval.GetKeywords();
    const Currency& in1 = inputs[0];
    std::sort(kwds.begin(), kwds.end());

    get1DStringsFromInput(in1, varnameBases);

    if (nargin > 1)
    {
        get1DStringsFromInput(inputs[1], exceptions);
    }

    for (size_t i = 0; i < varnameBases.size(); ++i)
    {
        std::sort(exceptions.begin(), exceptions.end());
        std::string &base = varnameBases[i];

        // replace invalid character with an underscore
        size_t loc = base.find_first_of(HW_BAD_VARNAME_CHARS);
        while (loc != std::string::npos)
        {
            base[loc] = '_';
            loc = base.find_first_of(HW_BAD_VARNAME_CHARS);
        }

        // loop in case we have a keyword that begins with an underscore
        while (std::count_if(base.begin(), base.end(), static_cast<int (*)(int)>(&isdigit)) == base.length() || std::binary_search(kwds.begin(), kwds.end(), base))
            base = "_" + base;

        // don't add appending numbers if we don't need to
        std::string varname;
        if (std::binary_search(exceptions.begin(), exceptions.end(), base))
        {
            if (isdigit(base.back()))
                base.push_back('_');

            int sup = 0;
            while (++sup)
            {
                varname = base + std::to_string((long long) sup);
                if (!std::binary_search(exceptions.begin(), exceptions.end(), varname))
                    break;
            }
        }
        else
            varname = base;

        varnames.push_back(varname);
        exceptions.push_back(varname);
    }

    if (in1.IsCellArray())
    {
        HML_CELLARRAY* incell = in1.CellArray();
        HML_CELLARRAY* out = EvaluatorInterface::allocateCellArray(incell->M(), incell->N());
        for (size_t i = 0; i < varnames.size(); ++i)
        {
            (*out)((int) i) = varnames[i];
        }
        outputs.push_back(out);
    }
    else
    {
        outputs.push_back(varnames[0]);
    }
    return true;
}
//------------------------------------------------------------------------------
// Converts subscript to a linear index [sub2ind]
//------------------------------------------------------------------------------
bool oml_sub2ind(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    int prod = 1;
    const Currency &in1 = inputs[0];
    std::vector<int> dims, cumprod;
    std::vector<std::vector<int> > indices;

    // get matrix dimensions from input
    if (in1.IsScalar())
    {
        int i = posIntFromDouble(in1.Scalar());
        if (!i)
            throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DIM);

        dims.push_back(i);
        cumprod.push_back(prod);
        prod *= i;
    }
    else if (in1.IsComplex())
    {
        throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DIM);
    }
    else if (in1.IsMatrix())
    {
        const hwMatrix *mtx = in1.Matrix();
        int size = mtx->Size();
        if (!size)
            throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DIMS);

        if (!mtx->IsRealData())
            throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DIM);

        for (int i = 0; i < size; ++i)
        {
            int ii = posIntFromDouble(realval(mtx,i));
            if (!ii)
                throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DIM);
            dims.push_back(ii);
            cumprod.push_back(prod);
            prod *= ii;
        }
    }
    else
        throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DIM);

    // get specified indices from input
    int m = -1;
    int n = -1;

    for (int i = 1; i < static_cast<int>(nargin); ++i)
    {
        std::vector<int> index;
        const Currency input = inputs[i];

        if (input.IsLogical())
            throw OML_Error(HW_ERROR_INPUTISLOGICAL);

        if (input.IsScalar())
        {
            if (m == -1 || (m == 1 && n == 1))
            {
                m = n = 1;
                int ii = posIntFromDouble(input.Scalar());
                if (!ii)
                    throw OML_Error(OML_ERR_POSINTEGER, i+1, OML_VAR_INDEX);
                index.push_back(ii);
            }
            else
                throw OML_Error(HW_ERROR_INDINPSAMESIZE);
        }
        else if (input.IsComplex())
        {
            throw OML_Error(OML_ERR_POSINTEGER, i+1, OML_VAR_INDEX);
        }
        else if (input.IsMatrix())
        {
            const hwMatrix *mtx = input.Matrix();
            int mtxsize = mtx->Size();

            if (!mtx->IsRealData())
                throw OML_Error(OML_ERR_POSINTEGER, i+1, OML_VAR_INDEX);

            if (m == -1)
            {
                m = mtx->M();
                n = mtx->N();
            }
            else if (!(m == mtx->M() && n == mtx->N()))
                throw OML_Error(HW_ERROR_INDINPSAMESIZE);

            for (int j = 0; j < mtxsize; ++j)
            {
                int ii = posIntFromDouble(realval(mtx, j));
                if (!ii)
                    throw OML_Error(OML_ERR_POSINTEGER, i+1, OML_VAR_INDEX);
                index.push_back(ii);
            }
        }
        else
            throw OML_Error(OML_ERR_POSINTEGER, i+1, OML_VAR_INDEX);

        indices.push_back(index);
    }

    // determine linear index
    hwMatrix *outmtx = EvaluatorInterface::allocateMatrix(max(m, 0), max(n, 0), 1.0);
    Currency out(outmtx);

    if (nargin == 2)
    {
        std::vector<int> index = indices[0];
        for (size_t i = 0; i < outmtx->Size(); ++i)
        {
            int temp = index[i];
            if (temp > prod)
                throw OML_Error(HW_ERROR_INDGRDIM);
            (*outmtx)((int) i) = temp;
        }
        outputs.push_back(out);
        return true;
    }

    for (size_t i = 0; i < indices.size(); ++i)
    {
        std::vector<int> index = indices[i];
        int thisdim, thisprod;
        if (i >= dims.size())
        {
            thisdim = 1;
            thisprod = cumprod.back();
        }
        else
        {
            thisdim = dims[i];
            thisprod = cumprod[i];
        }

        for (size_t j = 0; j < index.size(); ++j)
        {
            int currentIndex = index[j];
            if (currentIndex > thisdim)
                throw OML_Error("Error: index (" + std::to_string((long long) currentIndex) + ") exceeds dimension (" + std::to_string((long long) thisdim) + ")");

            (*outmtx)((int) j) += thisprod * (currentIndex - 1);
        }
    }

    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Returns a matrix with elements that are not in another specified matrix [setdiff]
//------------------------------------------------------------------------------
bool oml_setdiff(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    // don't output third index because it will always be an empty matrix -- it's useless
    std::vector<Currency> temp;
    sortBasedOperation(eval, inputs, temp, false, &dosetdiff, &dosetdiff, &dosetdiff, &dosetdiff);
    for (size_t i = 0; i < 2 && i < temp.size(); ++i)
        outputs.push_back(temp[i]);
    return true;
}
//------------------------------------------------------------------------------
// Returns the elements in both a and b [intersect]
//------------------------------------------------------------------------------
bool oml_intersect(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return sortBasedOperation(eval, inputs, outputs, true, &dointersect, &dointersect, &dointersect, &dointersect);
}
//------------------------------------------------------------------------------
// Returns the elements exclusive to a or b [setxor]
//------------------------------------------------------------------------------
bool oml_setxor(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return sortBasedOperation(eval, inputs, outputs, false, &dosetxor, &dosetxor, &dosetxor, &dosetxor);
}
//------------------------------------------------------------------------------
// Converts input into double precision type [double]
//------------------------------------------------------------------------------
bool oml_double(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

	Currency ret = inputs[0];

    if (ret.IsScalar() || ret.IsComplex() || ret.IsMatrix() || ret.IsString())
	{
		ret.SetMask(Currency::MASK_DOUBLE);
        outputs.push_back(ret);
	}
    else
	{
        throw OML_Error(HW_ERROR_NOTCONVINPTODOUBLE);
	}

    return true;
}
//------------------------------------------------------------------------------
// Returns true if successul in string matching for regular expressions [regexp]
//------------------------------------------------------------------------------
bool oml_regexp(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.empty() || inputs.size() < 2) throw OML_Error(OML_ERR_NUMARGIN);

    // Get output options, in user specified order if applicable
    std::vector<std::string> options (GetRegExpOptions(eval, inputs));

    bool result = GetRegExpOutput(eval, inputs[0], inputs[1], options, outputs);
    return result;
}
//------------------------------------------------------------------------------
// Converts cell to struct [cell2struct]
//------------------------------------------------------------------------------
bool oml_cell2struct(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    HML_CELLARRAY *cell;
    if (inputs[0].IsCellArray())
        cell = inputs[0].CellArray();
    else
        throw OML_Error(OML_ERR_CELL, 1, OML_VAR_TYPE);

    int dim;
    if (nargin > 2)
    {
        if (inputs[2].IsScalar())
            dim = (int) inputs[2].Scalar();
        else
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);

        if (dim < 1)
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);

        if (dim > 2)
            throw OML_Error(OML_ERR_UNSUPPORTDIM, 3);
    }
    else
    {
        dim = cell->M() == 1 ? 2 : 1;
    }

    std::vector<std::string> fields;
    const Currency &fieldcur = inputs[1];
    int size_along_dim = dim == 1 ? cell->M() : cell->N();
    if (fieldcur.IsString())
    {
        int m = fieldcur.Matrix()->M();
        if (m != size_along_dim)
            throw OML_Error(HW_ERROR_FIELDNAMEDIMINPCELL);
        for (int i = 0; i < m; i++)
        {
            std::string f = readString(fieldcur, i);
            if (f.empty())
                throw OML_Error(HW_ERROR_FIELDNAMENOTEMPTYSTR);
            fields.push_back(f);
        }
    }
    else if (fieldcur.IsCellArray())
    {
        HML_CELLARRAY *fcell = fieldcur.CellArray();
        int size = fcell->Size();
        if (size != size_along_dim)
            throw OML_Error(HW_ERROR_FIELDNAMEDIMINPCELL);

        if (!isstr(fcell))
            throw OML_Error(HW_ERROR_FIELDNAMECELLSTR);

        for (int i = 0; i < size; i++)
        {
            std::string f = readString((*fcell)(i));
            if (f.empty())
                throw OML_Error(HW_ERROR_FIELDNAMENOTEMPTYSTR);
            fields.push_back(f);
        }
    }
    else
        throw OML_Error(HW_ERROR_2NDINPCELLSTR);

    // verify there are no duplicates in field names
    std::set<std::string> fieldset(fields.begin(), fields.end());
    if (fieldset.size() != fields.size())
        throw OML_Error(HW_ERROR_DUPFIELD);

    StructData *sd = new StructData();
    Currency outstruct(sd);

    // this will have to be changed for ND cell inputs
    sd->DimensionNew((dim == 1 ? cell->N() : cell->M()), 1);

    int struct_size = sd->M() * sd->N();
    for (size_t i = 0; i < fields.size(); i++) // fields.size() == size_along_dim
    {
        for (int j = 0; j < struct_size; j++)
        {
            sd->SetValue(j, -1, fields[i], dim == 1 ? (*cell)((int) i, j) : (*cell)(j, (const int)i));
        }
    }

    outputs.push_back(outstruct);
    return true;
}
//------------------------------------------------------------------------------
// Sorts input into its ascending pairs by real part [cplxpair]
//------------------------------------------------------------------------------
bool oml_cplxpair(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    double tol = 100 * MachPrecision(1.0); // default tolerance
    int dim;

    const hwMatrix *mtx = NULL;
    std::deque<hwComplex> cplxs;
    std::vector<double> reals;

    const Currency &cur = inputs[0];
    if (cur.IsScalar())
    {
        outputs.push_back(cur);
        return true;
    }
    else if (cur.IsComplex())
    {
        if (cur.Complex().IsReal(tol))
        {
            outputs.push_back(cur);
            return true;
        }
        throw OML_Error(HW_ERROR_INPUTNOCONJ);
    }
    else if (cur.IsMatrix() || cur.IsString())
    {
        mtx = cur.Matrix();
        if (mtx->IsEmpty())
        {
            outputs.push_back(cur);
            return true;
        }
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

    if (nargin > 1)
    {
        if (!inputs[1].IsScalar())
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

        tol = inputs[1].Scalar();
    }

    if (nargin > 2)
    {
        dim = (int) inputs[2].Scalar();

        if (!IsInteger(inputs[2].Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 3, OML_VAR_DIM);
    }
    else
    {
        dim = mtx->M() == 1 ? 2 : 1;
    }

    if (mtx->IsReal())
    {
		int loop_limit = dim == 1 ? mtx->M() : mtx->N();
        for (int i = 0; i < loop_limit; i++)
            reals.push_back((*mtx)(i));
    }
    else
    {
        for (int i = 0; i < mtx->Size(); i++)
        {
            hwComplex c = mtx->z(i);
            if (c.IsReal(tol))
                reals.push_back(c.Real());
            else
                cplxs.push_back(c);
        }
    }

    std::sort(reals.begin(), reals.end());
    std::sort(cplxs.begin(), cplxs.end(), std::bind(complexLessThanTol, std::placeholders::_1, std::placeholders::_2, tol));
    std::vector<hwComplex> retvec;
    while (cplxs.size())
    {
        const hwComplex &c = cplxs[0];
        hwComplex conj = c.Conjugate();
        size_t i = (size_t) indexOfTol(cplxs, conj, tol);
        if (i == -1)
            throw OML_Error(HW_ERROR_INPUTNOCONJ);

        conj = cplxs[i];
        if (c.Imag() < 0)
        {
            retvec.push_back(c);
            retvec.push_back(conj);
        }
        else
        {
            retvec.push_back(conj);
            retvec.push_back(c);
        }

        cplxs.erase(cplxs.begin() + i);
        cplxs.erase(cplxs.begin());
    }

    for (std::vector<double>::iterator i = reals.begin(); i != reals.end(); i++)
        retvec.push_back(hwComplex(*i, 0.0));

    outputs.push_back(containerToMatrix(retvec));
    return true;
}
//------------------------------------------------------------------------------
// Returns subscript from the linear index (input) [ind2sub]
//------------------------------------------------------------------------------
bool oml_ind2sub(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &in1 = inputs[0];
    const Currency &in2 = inputs[1];

    int nargout = max(getNumOutputs(eval), 1);
    int max_size = 1;
    std::vector<int> dimlist;
    if (in1.IsScalar() || in1.IsComplex())
    {
        double d = in1.IsScalar() ? in1.Scalar() : in1.Complex().Real();
        if (!isposint(d))
            throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DIM);
        max_size = (int) d;
        dimlist.push_back(max_size);
    }
    else if (in1.IsMatrix())
    {
        const hwMatrix *mtx = in1.Matrix();
        if (!mtx->Size())
            throw OML_Error(OML_ERR_POSINTEGER, 1);

        for (int i = 0 ; i < mtx->Size(); i++)
        {
            double d = realval(mtx, i);
            if (!isposint(d))
                throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DIM);
            max_size *= (int) d;
            dimlist.push_back((int) d);
        }
    }
    else
        throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DIM);

    // determine the indices to convert from second input
    std::pair<int, int> dim_out;
    std::vector<int> indices;

    if (in2.IsScalar())
    {
        double d = in2.Scalar();
        if (!isposint(d))
            OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_INDEX);
        if ((int) d > max_size)
            throw OML_Error(HW_ERROR_INDGRDIM);
        dim_out.first = dim_out.second = 1;
        indices.push_back((int) d);

    }
    else if (in2.IsComplex())
    {
        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_INDEX);
    }
    else if (in2.IsMatrix() || in2.IsString())
    {
        const hwMatrix *mtx = in2.Matrix();
        dim_out.first = mtx->M();
        dim_out.second = mtx->N();
        if (mtx->IsRealData())
        {
            for (int i = 0; i < mtx->Size(); i++)
            {
                double d = realval(mtx, i);
                if (!isposint(d))
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_INDEX);
                if ((int) d > max_size)
                    throw OML_Error(HW_ERROR_INDGRDIM);
                indices.push_back((int) d);
            }
        }
        else
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_INDEX);
    }
    else
        throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_INDEX);

    // prepare output matrices
    std::vector<hwMatrix*> index_matrices;
    std::vector<Currency> out_curs;
    for (int i = 0; i < nargout; i++)
    {
        hwMatrix *mtx = EvaluatorInterface::allocateMatrix(dim_out.first, dim_out.second, 1.0);
        index_matrices.push_back(mtx);
        out_curs.push_back(mtx);
    }

    // do the calculation
    for (size_t i = 0; i < indices.size(); i++)
    {
        int lindex = indices[i] - 1;
        for (size_t j = 0; j < dimlist.size() && j < nargout; j++)
        {
            int dim = dimlist[j];
            hwMatrix *m = index_matrices[j];
            (*m)((const int)i) = (double) (lindex % dim) + 1.0;
            lindex /= dim;
        }
    }

    // output values
    for (size_t i = 0; i < index_matrices.size(); i++)
    {
        outputs.push_back(out_curs[i]);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is prime [isprime]
//------------------------------------------------------------------------------
bool oml_isprime(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

	if (inputs[0].IsScalar())
    {
		double d = inputs[0].Scalar();
        if (!isint(d) || d < 0.0)
            throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_VALUE);
        int p = (int) d;

		Currency out((double) isprime(p));
		out.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(out);
         
	}
	else if (inputs[0].IsMatrix())
	{
		const hwMatrix* mtx = (hwMatrix*)inputs[0].Matrix();

		if (!mtx->IsReal())
			throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_VALUE);

        hwMatrix *out = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

        for (int i = 0; i < mtx->Size(); i++)
        {
            Currency k = (*mtx)(i);

			if (!k.IsScalar())
				throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_VALUE);

			double p = k.Scalar();

			if (!isint(p) || p < 0.0)
				throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_VALUE);
				
            (*out)(i) = (double) isprime((int) p);
        }

        Currency out_cur(out);
        out_cur.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(out_cur);
	}
    else if (inputs[0].IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_isprime);
    }
	else
	{
		throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_VALUE);
	}
    
    return true;
}
//------------------------------------------------------------------------------
// Computes 2D convolution [conv2]
//------------------------------------------------------------------------------
bool oml_conv2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    // note -- this produces different results than in hml

    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar() && !inputs[1].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);

    hwMatrix a(*inputs[0].ConvertToMatrix());
    hwMatrix b(*inputs[1].ConvertToMatrix());
    std::unique_ptr<hwMatrix> mtx;

    hwMatrix *result = EvaluatorInterface::allocateMatrix();

    Currency rcur(result);
    std::string shape("full");

    enum
    {
        Full,
        Same,
        Valid
    } shapeCode;

    if (nargin == 3)
    {
        const Currency &in3 = inputs[2];
        if (in3.IsString())
        {
            shape = readOption(eval, in3);
        }
        else
        {
            mtx.reset(safeMatrixCopyFromInput(in3, false));
        }
    }
    else if (nargin == 4)
    {
        mtx.reset(matrixCopyFromInput(inputs[2], false));
        shape = readOption(eval, inputs[3]);
    }

    if (mtx.get() && !(a.IsVector() && b.IsVector()))
        throw OML_Error(HW_ERROR_IFCONVMATINPMUSTVEC);

    if (shape == "full")
        shapeCode = Full;
    else if (shape == "same")
        shapeCode = Same;
    else if (shape == "valid")
        shapeCode = Valid;
    else
        throw OML_Error(HW_ERROR_INVSHAPEFULLSAMEVALID);

    if (!(b.Size() && a.Size()))
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        return true;
    }

    int top, left, m, n;

    if (mtx.get())
    {
        int bs = b.Size();
        int as = a.Size();
        int mm = mtx->M();
        int nn = mtx->N();

        if (mm * nn == 0)
        {
            outputs.push_back(EvaluatorInterface::allocateMatrix());
            return true;
        }

        m = max(as, mm);
        n = max(bs, nn);
        if (a.N() != 1)
            safeResizeMatrix(eval, 1, m, a, true);
        else
            safeResizeMatrix(eval, m, 1, a, true);

        if (b.M() != 1)
            safeResizeMatrix(eval, n, 1, b, true);
        else
            safeResizeMatrix(eval, 1, n, b, true);

        safeResizeMatrix(eval, m, n, *mtx, true);
        BuiltInFuncsUtils::CheckMathStatus(eval, result->Conv2D(a, b, *mtx));

        if (shapeCode == Same)
        {
            top = as / 2;
            left = bs / 2;
            m = mm;
            n = nn;
        }
        else if (shapeCode == Valid)
        {
            top = m / 2 - 1;
            left = n / 2 - 1;
            m = mm - as + 1;
            n = nn - bs + 1;
        }
        else
        {
            top = left = 0;
            m = mm + as - 1;
            n = nn + bs - 1;
        }
    }
    else
    {
        int am = a.M();
        int an = a.N();
        int bm = b.M();
        int bn = b.N();
        m = max(am, bm);
        n = max(an, bn);
        safeResizeMatrix(eval, m, n, a, true);
        safeResizeMatrix(eval, m, n, b, true);
        BuiltInFuncsUtils::CheckMathStatus(eval, result->Conv2D(a, b));

        top = bm / 2;
        left = bn / 2;
        if (shapeCode == Same)
        {
            m = am;
            n = an;
        }
        else if (shapeCode == Valid)
        {
            m = am - bm + 1;
            n = an - bn + 1;
        }
        else
        {
            top = left = 0;
            m = am + bm - 1;
            n = an + bn - 1;
        }
    }

    hwMatrix *temp = getInnerMatrix(result, top, top + m, left, left + n);
    result = temp;
    rcur = Currency(result);
    outputs.push_back(rcur);
    return true;
}
//------------------------------------------------------------------------------
// Helper method for matrix concatenation - assumes both inputs of the same type
// real or complex
//------------------------------------------------------------------------------
template<bool HORIZ, typename T1, typename T2>
static hwTMatrix<T1, T2>* concat(const hwTMatrix<T1, T2> *m1, const hwTMatrix<T1, T2> *m2)
{
    if (m1->M() != m2->M())
    {
        if ((m1->M() == 0) && (m1->N() == 0))
        {
            return new hwTMatrix<T1, T2>(*m2);
        }
        else if ((m2->M() == 0) && (m2->N() == 0))
        {
            return new hwTMatrix<T1, T2>(*m1);
        }
        else
            throw OML_Error(HORIZ ? HW_ERROR_INPUTHCATDIM : HW_ERROR_INPUTVCATDIM);
    }

    if (m1->IsReal())
    {
        hwTMatrix<T1, T2> *newmtx = new hwTMatrix<T1, T2>(m1->M(), m1->N() + m2->N(),
            hwTMatrix<T1, T2>::REAL);

        int i;
        for (i = 0; i < m1->Size(); i++)
            (*newmtx)(i) = (*m1)(i);

        for (int j = i; j < newmtx->Size(); j++)
            (*newmtx)(j) = (*m2)(j - i);

        return newmtx;
    }

    hwTMatrix<T1, T2> *newmtx = new hwTMatrix<T1, T2>(m1->M(), m1->N() + m2->N(),
        hwTMatrix<T1, T2>::COMPLEX);

    int i;
    for (i = 0; i < m1->Size(); i++)
        newmtx->z(i) = m1->z(i);

    for (int j = i; j < newmtx->Size(); j++)
        newmtx->z(j) = m2->z(j - i);

    return newmtx;
}
//------------------------------------------------------------------------------
// Helper method for cell array concatenation
//------------------------------------------------------------------------------
template <bool HORIZ>
static HML_CELLARRAY* concat(EvaluatorInterface& eval, const HML_CELLARRAY *cell, Currency right)
{
    if (right.IsCellArray())
        return concat<HORIZ>(cell, right.CellArray());
    if (!HORIZ)
        right = transpose(eval, right);

    if (cell->M() != 1)
    {
        if (!cell->Size())
        {
            HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(1, 1);
            (*newcell)(0) = right;
            return newcell;
        }
        else
            throw OML_Error(HORIZ ? HW_ERROR_INPUTHCATDIM : HW_ERROR_INPUTVCATDIM);
    }

    HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(cell->M(), cell->N() + 1);

    int i;
    for (i = 0; i < cell->Size(); i++)
        (*newcell)(i) = (*cell)(i);

    (*newcell)(i) = right;

    return newcell;
}
//------------------------------------------------------------------------------
// Helper method for cell concatenation
//------------------------------------------------------------------------------
template <bool HORIZ>
static HML_CELLARRAY* concat(EvaluatorInterface& eval, Currency left, const HML_CELLARRAY *cell)
{
    if (left.IsCellArray())
        return concat<HORIZ>(left.CellArray(), cell);
    if (!HORIZ)
        left = transpose(eval, left);

    if (cell->Size())
    {
        if (cell->M() != 1)
            throw OML_Error(HORIZ ? HW_ERROR_INPUTHCATDIM : HW_ERROR_INPUTVCATDIM);

        HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(cell->M(), cell->N() + 1);

        (*newcell)(0) = left;

        for (int i = 0; i < cell->Size(); i++)
            (*newcell)(i + 1) = (*cell)(i);

        return newcell;
    }
    else
    {
        HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(1, 1);
        (*newcell)(0) = left;
        return newcell;
    }
}
//------------------------------------------------------------------------------
// Helper method for concatenation
//------------------------------------------------------------------------------
template <bool HORIZ>
static Currency concat(EvaluatorInterface& eval, const Currency &left, const Currency &right)
{
    if (left.IsCellArray())
        return concat<HORIZ>(eval, left.CellArray(), right);
    else if (right.IsCellArray())
        return concat<HORIZ>(eval, left, right.CellArray());

    if (left.IsScalar() || left.IsComplex() || left.IsMatrix() || left.IsString())
    {
        hwMatrix *lm = matrixCopyFromInput(left, true);
        Currency lc(lm);
        hwMatrix *rm = matrixCopyFromInput(right, true);
        if (!rm)
            throw OML_Error(HW_ERROR_CONCATMATWSCALCOMPMATORCELL);
        Currency rc(rm);

        bool complex = checkMakeComplex(eval, lm, rm);

        if (left.IsString() || right.IsString())
        {
            if (complex)
            {
                lm = tostr(eval, lm, false);
                rm = tostr(eval, rm, false);
                lc = Currency(lm);
            }
            return addStringMask(concat<HORIZ>(lm, rm));
        }

        return concat<HORIZ>(lm, rm);
    }
    else if (left.IsStruct())
    {
        if (right.IsStruct())
        {
            StructData *l = left.Struct();
            StructData *r = right.Struct();
            int m = l->M();
            int ln = l->N();
            int rn = r->N();

            if (m != r->M())
            {
                if (m * ln == 0)
                    return new StructData(*r);
                else if (m * rn == 0)
                    return new StructData(*l);
                else
                    throw OML_Error(HW_ERROR_INPUTHCATDIM);
            }

            const std::map<std::string, int> lfields = l->GetFieldNames(), rfields = r->GetFieldNames();
            if (lfields.size() != rfields.size())
                throw OML_Error(HW_ERROR_INPUTCATSTRUCT);

            std::map<std::string, int>::const_iterator iter;
            for (iter = lfields.cbegin(); iter != lfields.cend(); iter++)
            {
                if (!rfields.count(iter->first))
                    throw OML_Error(HW_ERROR_INPUTCATSTRUCT);
            }

            StructData *result = new StructData();
            Currency rescur(result);

            result->DimensionNew(m, ln + rn);

            int i, size = m * ln;
            for (i = 0; i < size; i++)
            {
                for (iter = lfields.cbegin(); iter != lfields.cend(); iter++)
                {
                    result->SetValue(i, -1, iter->first, l->GetValue(i, -1, iter->first));
                }
            }

            size = m * (ln + rn);
            for (int j = i; j < size; j++)
            {
                for (iter = lfields.cbegin(); iter != lfields.cend(); iter++)
                {
                    result->SetValue(j, -1, iter->first, r->GetValue(j - i, -1, iter->first));
                }
            }

            return rescur;
        }
        else
            throw OML_Error(HW_ERROR_CONCATSTRUCTWSTRUCT);
    }
    else
    {
        throw OML_Error(HW_ERROR_NOTCONCINPTYPE);
    }
}
//------------------------------------------------------------------------------
// Vertically concatenates a1 through aN along the first dimension [vertcat]
//------------------------------------------------------------------------------
bool oml_vertcat(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    // transposes each input, then calls horzcat - the final output is again transposed
    // although less efficient, this approach saves a lot of duplicated code

    size_t nargin = inputs.size();
    bool any_ND = false;

    for (int i = 0; i < nargin; ++i)
    {
        const Currency& cur = inputs[i];

        if (cur.IsNDMatrix())
        {
            any_ND = true;
            break;
        }
    }

    if (nargin && any_ND)
    {
        std::vector<Currency> new_inputs;

        new_inputs.push_back(1);
        
        for (std::vector<Currency>::const_iterator it = inputs.begin(); it != inputs.end(); it++)
            new_inputs.push_back(*it);

        try
        {
            return oml_cat(eval, new_inputs, outputs);
        }
        catch (OML_Error& e)
        {
            if (e.Arg1() != -1)
                e.Arg1(e.Arg1() - 1);

            if (e.Arg2() != -1)
                e.Arg2(e.Arg2() - 1);

            throw e;
        }
    }

    if (nargin)
    {
        Currency current = getFirstInputOrCell(inputs);
        current = transpose(eval, current);
        for (size_t i = current.IsCellArray() ? 0 : 1; i < nargin; i++)
        {
            current = vconcat(eval, current, transpose(eval, inputs[i]));
        }
        outputs.push_back(transpose(eval, current));
    }
    else
        outputs.push_back(EvaluatorInterface::allocateMatrix());
    return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template <int (*F)(int)>
static double wrapLogical(double d)
{
    static const int modby = (const int)pow(2.0, (int) sizeof(char) * 8);
    int i = ((int) d) % modby;

    if (i < 0)
        i += modby;

    return (*F)(i) ? 1.0 : 0.0;
}
//------------------------------------------------------------------------------
// Returns true if input or elements of input are letters [isletter]
//------------------------------------------------------------------------------
bool oml_isletter(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return scalarToLogicalFunc(inputs, outputs, &wrapLogical<&isalpha>);
}
//------------------------------------------------------------------------------
// Returns true if input or elements of input are spaces [isspace]
//------------------------------------------------------------------------------
bool oml_isspace(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return scalarToLogicalFunc(inputs, outputs, &wrapLogical<&isspace>);
}
//------------------------------------------------------------------------------
// Horizontally concatenates inputs [horzcat]
//------------------------------------------------------------------------------
bool oml_horzcat(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    bool any_ND = false;

    for (int i = 0; i < nargin; ++i)
    {
        const Currency& cur = inputs[i];

        if (cur.IsNDMatrix())
        {
            any_ND = true;
            break;
        }
    }

    if (nargin && any_ND)
    {
        std::vector<Currency> new_inputs;

        new_inputs.push_back(2);
        
        for (std::vector<Currency>::const_iterator it = inputs.begin(); it != inputs.end(); it++)
            new_inputs.push_back(*it);

        try
        {
            return oml_cat(eval, new_inputs, outputs);
        }
        catch (OML_Error& e)
        {
            if (e.Arg1() != -1)
                e.Arg1(e.Arg1() - 1);

            if (e.Arg2() != -1)
                e.Arg2(e.Arg2() - 1);

            throw e;
        }
    }

    if (nargin)
    {
        Currency current = getFirstInputOrCell(inputs);
        for (size_t i = current.IsCellArray() ? 0 : 1; i < nargin; i++)
        {
            current = hconcat(eval, current, inputs[i]);
        }
        outputs.push_back(current);
    }
    else
        outputs.push_back(EvaluatorInterface::allocateMatrix());
    return true;
}
//------------------------------------------------------------------------------
// Recursively display contents of cell array [celldisp]
//------------------------------------------------------------------------------
bool oml_celldisp(EvaluatorInterface           eval, 
                  const std::vector<Currency>& inputs, 
                  std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2) throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];
    if (!input.IsCellArray()) throw OML_Error(OML_ERR_CELLARRAY, 1, OML_VAR_TYPE);

    std::string cellname (input.GetOutputName());

    if (nargin == 2)
    {
        const Currency& cur = inputs[1];
        if (!cur.IsString()) throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);

        cellname = orderedStringVal(toCurrencyStr(eval, cur, false, false));
    }

    celldisp(eval, input.CellArray(), cellname);
    return true;
}
//------------------------------------------------------------------------------
// Returns the greatest common deviser of the input [gcd]
//------------------------------------------------------------------------------
bool oml_gcd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() < 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    bool complex = false;
    int size;
    int m = -1, n = -1;

    // verify input and get output size
    for (size_t i = 0; i < inputs.size(); i++)
    {
        const Currency &cur = inputs[i];
        if (cur.IsMatrix())
        {
            const hwMatrix *mtx = cur.Matrix();
            if (!isint(mtx))
                throw OML_Error(OML_ERR_INTEGER, (int)i+1, OML_VAR_VALUE);
            complex = complex || !mtx->IsReal();
            if (mtx->Size() != 1)
            {
                if (m == -1)
                {
                    m = mtx->M();
                    n = mtx->N();
                }
                else if (m != mtx->M() || n != mtx->N())
                    throw OML_Error(OML_ERR_ARRAYSIZE);
            }
        }
        else if (cur.IsComplex())
        {
            if (!isint(cur.Complex()))
                throw OML_Error(OML_ERR_INTEGER, (int)i+1, OML_VAR_VALUE);
            complex = true;
        }
        else if (cur.IsScalar())
        {
            if (!isint(cur.Scalar()))
                throw OML_Error(OML_ERR_INTEGER, (int)i+1, OML_VAR_VALUE);
        }
        else
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    size = m * n;
    if (m < 0)
        m = 1;
    if (n < 0)
        n = 1;

    hwMatrix *gcds = EvaluatorInterface::allocateMatrix(m, n, complex ? hwMatrix::COMPLEX : hwMatrix::REAL);
    if (!size)
    {
        outputs.push_back(gcds);
        return true;
    }

    // calculate gcds and coefs
    std::vector<hwMatrix*> coefs;
    for (size_t i = 0; i < inputs.size(); i++)
        coefs.push_back(EvaluatorInterface::allocateMatrix(m, n, complex ? hwMatrix::COMPLEX : hwMatrix::REAL));

    // i is input index
    // j is matrix index
    // k is coefs index
    for (size_t i = 0; i < inputs.size(); i++)
    {
        const Currency &input = inputs[i];
        if (i)
        {
            if (complex)
            {
                hwComplex c;
                const hwMatrix *mtx;
                bool usemtx = false;
                if (input.IsScalar())
                    c = hwComplex((int) input.Scalar(), 0.0);
                else if (input.IsComplex())
                    c = input.Complex();
                else // input.IsMatrix()
                {
                    mtx = input.Matrix();
                    usemtx = true;
                }

                for (int j = 0; j < gcds->Size(); j++)
                {
                    hwComplex s, t;
                    gcds->z(j) = gcd(gcds->z(j), usemtx ? mtx->IsReal() ? hwComplex((*mtx)(j), 0.0) : mtx->z(j) : c, &s, &t);
                    for (size_t k = 0; k < i; k++)
                    {
                        coefs[k]->z(j) *= s;
                    }
                    coefs[i]->z(j) = t;
                }
            }
            else
            {
                double d;
                const hwMatrix *mtx;
                bool usemtx = false;
                if (input.IsScalar())
                    d = input.Scalar();
                else
                {
                    usemtx = true;
                    mtx = input.Matrix();
                }
                for (int j = 0; j < gcds->Size(); j++)
                {
                    int s, t;
                    (*gcds)(j) = gcd((int) (*gcds)(j), (int) (usemtx ? (*mtx)(j) : d), &s, &t);
                    for (size_t k = 0; k < i; k++)
                    {
                        (*coefs[k])(j) *= s;
                    }
                    (*coefs[i])(j) = t;
                }
            }
        }
        else // i == 0
        {
            if (input.IsScalar())
                gcds->SetElements(input.Scalar());
            else if (input.IsComplex())
                gcds->SetElements(input.Complex());
            else // input.IsMatrix()
            {
                *gcds = *input.Matrix();
                if (complex && gcds->IsReal())
                    BuiltInFuncsUtils::CheckMathStatus(eval, gcds->MakeComplex());
            }
            coefs[0]->SetElements(1.0);
        }
    }

    outputs.push_back(gcds);
    if (size == 1 && getNumOutputs(eval) == 2)
    {
        // convert output coefs into matrix
        hwMatrix *out = EvaluatorInterface::allocateMatrix(1, (int)coefs.size(), complex ? hwMatrix::COMPLEX : hwMatrix::REAL);
        for (size_t i = 0; i < coefs.size(); i++)
        {
            hwMatrix *c = coefs[i];
            if (complex)
                out->z((const int)i) = c->z(0);
            else
                (*out)((const int)i) = (*c)(0);
            delete c;
        }
        outputs.push_back(out);
    }
    else
    {
        for (size_t i = 0; i < coefs.size(); i++)
            outputs.push_back(coefs[i]);
    }
    return true;
}
//------------------------------------------------------------------------------
// Displays or sets the separator in the given path [pathsep]
//------------------------------------------------------------------------------
bool oml_pathsep(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
	if (nargin == 0)
	{
		outputs.push_back(addStringMask((double) pathsep));
	}
    else if (nargin == 1)
    {
        std::string temp = readString(toCurrencyStr(eval, inputs[0], false, false));
        if (temp.length() > 1)
            throw OML_Error(HW_ERROR_PATHSEP1CHAR);
        pathsep = temp.length() ? temp[0] : 0;
		outputs.push_back(addStringMask((double) pathsep));
    }
    else if (nargin)
        throw OML_Error(OML_ERR_NUMARGIN);

    return true;
}
//------------------------------------------------------------------------------
// Returns the name of the given function [func2str]
//------------------------------------------------------------------------------
bool oml_func2str(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsFunctionHandle())
    {
        outputs.push_back(inputs[0].FunctionHandle()->FunctionName());
        return true;
    }
    else
        throw OML_Error(OML_ERR_HANDLE, -1, OML_VAR_TYPE);
}
//------------------------------------------------------------------------------
// Returns the function handle for the given function name [str2func]
//------------------------------------------------------------------------------
bool oml_str2func(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency    input    = unnest(inputs[0]);
    std::string funcName = readString(toCurrencyStr(eval, input, false, false));

    FunctionInfo *fi   = nullptr;
    FUNCPTR       fptr = nullptr;

    if (eval.FindFunctionByName(funcName, &fi, &fptr))
    {
		if (fptr)
			fi = new FunctionInfo(funcName, fptr);
		else
			fi = new FunctionInfo(*fi);

        outputs.push_back(fi);
        return true;
    }
	else if (funcName[0] == '@')
	{
		FunctionInfo* fi = eval.FunctionInfoFromString(funcName);
		outputs.push_back(fi);
		return true;
	}
	else
	{
		throw OML_Error("Error: cannot find function '" + funcName + "'");
	}
}
//------------------------------------------------------------------------------
// Performs case-insensitive comparison of first n chars of inputs [strncmpi]
//------------------------------------------------------------------------------
bool oml_strncmpi(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    std::vector<Currency> lower = convertToLower(inputs.cbegin(), inputs.cend());
    return oml_strncmp(eval, lower, outputs);
}
//------------------------------------------------------------------------------
// Performs case-insensitive comparison of inputs [strcmpi]
//------------------------------------------------------------------------------
bool oml_strcmpi(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    std::vector<Currency> lower = convertToLower(inputs.cbegin(), inputs.cend());
    return oml_strcmp(eval, lower, outputs);
}
//------------------------------------------------------------------------------
// Evaluates the given function func on elements of given cell array [cellfun]
//------------------------------------------------------------------------------
bool oml_cellfun(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    bool uniformOutput = true;
    int rawnargout = getNumOutputs(eval);
    int nargout = max(rawnargout, 1);
    std::string func;

    if (input1.IsString())
        func = readString(input1);
    else if (!input1.IsFunctionHandle())
        throw OML_Error(HW_ERROR_FUNCNAMESTR);

    int cellm, celln;
    cellm = celln = -1;
    for (size_t i = 1; i < nargin; ++i)
    {
        const Currency& cur = inputs[i];
        if (cur.IsCellArray())
        {
            HML_CELLARRAY *incell = cur.CellArray();
            if (incell->Size() != 1)
            {
                cellm = incell->M();
                celln = incell->N();
                break;
            }
        }
    }

    if (cellm == -1)
        cellm = celln = 1;

    std::vector<Currency> args;
    for (size_t i = 1; i < nargin; ++i)
    {
        const Currency &in = inputs[i];
        if (in.IsCellArray())
        {
            HML_CELLARRAY *cell = in.CellArray();
            if (cell->Size() == 1)
            {
                Currency &elem = (*cell)(0);
                HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(cellm, celln);
                for (int j = 0; j < cellm * celln; ++j)
                    (*newcell)(j) = elem;
                args.push_back(newcell);
            }
            else if (cell->M() == cellm && cell->N() == celln)
                args.push_back(in);
            else
                throw OML_Error(HW_ERROR_CELLINPSAMESASIZE);
        }
        else if (in.IsString())
        {
            std::string opt = readOption(eval, in);
            if (opt == "uniformoutput")
            {
                if (++i < nargin)
                {
                    if (!inputs[i].IsScalar())
                        throw OML_Error(OML_ERR_REAL, (int)i+1, OML_VAR_VALUE);

                    double tobool = inputs[i].Scalar();
                    uniformOutput = isint(tobool) && (int) tobool == 0 ? false : true;
                }
                else
                    throw OML_Error(HW_ERROR_MISSVALUNIFOUTOPT);
            }
            else
                throw OML_Error(HW_ERROR_INVALIDOPTION(opt));
        }
        else if (nargin == 3 && i == 2 && func == "size")
        {
            HML_CELLARRAY* cell;
            cell = EvaluatorInterface::allocateCellArray(cellm, celln);

            for (int j = 0; j < cell->Size(); j++)
                (*cell)(j) = in;
            args.push_back(cell);
        }
        else
            throw OML_Error(HW_ERROR_INPUTALLCELL);
    }

    size_t size = cellm * celln;
    for (int i = 0; i < size; i++)
    {
        // set up arguments
        std::vector<Currency> newinputs, newoutputs;

        newinputs.push_back(input1);
        for (size_t j = 0; j < args.size(); j++)
            newinputs.push_back((*args[j].CellArray())(i));

        // call the function
        try
        {
            oml_feval(eval, newinputs, newoutputs);
        }
        catch (OML_Error& err)
        {
            if (err.Arg1() != -1)
                err.Arg1(err.Arg1() + 1);

            if (err.Arg2() != -1)
                err.Arg2(err.Arg2() + 1);

            throw err;
        }

        // add newoutputs to outputs
        if (uniformOutput)
        {
            if (newoutputs.size() && !newoutputs[0].IsNothing())
            {
                while (outputs.size() < nargout)
                    outputs.push_back(EvaluatorInterface::allocateMatrix(cellm, celln, 0.0));

                for (size_t j = 0; j < nargout && j < newoutputs.size(); j++)
                {
                    const Currency &cur = newoutputs[j];
                    if (cur.IsScalar())
                        outputs[j].GetWritableMatrix()->SetElement(i, cur.Scalar());
                    else if (cur.IsComplex())
                        outputs[j].GetWritableMatrix()->SetElement(i, cur.Complex());
                    else
                        throw OML_Error(HW_ERROR_OUTNOTUNI);
                }
            }
        }
        else
        {
            if (newoutputs.size() && !newoutputs[0].IsNothing())
            {
                while (outputs.size() < nargout)
                {
                    HML_CELLARRAY *cell = EvaluatorInterface::allocateCellArray(cellm, celln);
                    // in case nargout is different than in other calls
                    for (int i = 0; i < cell->Size(); i++)
                        (*cell)(i) = Currency();
                    outputs.push_back(cell);
                }

                for (size_t j = 0; j < outputs.size() && j < newoutputs.size(); j++)
                {
                    (*outputs[j].CellArray())(i) = newoutputs[j];
                }
            }
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the day/time input as a number of days since January 1, 0000 [datenum]
//------------------------------------------------------------------------------
bool oml_datenum(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	int nargout = getNumOutputs(eval);
	int nargin = (int)inputs.size();
    int m, n;

    std::vector<std::array<double, 6> > dates = datesFromInput(inputs, &m, &n);

    hwMatrix *result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
	hwMatrix *resultSecs = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);

	double daynum;

    for (int i = 0; i < dates.size(); i++)
    {
        std::array<double, 6> date = dates[i];
        double years = date[0];
        double months = date[1] - 1;

        years += (int) (months / 12.0);
        months = doubleMod(months, 12.0);

        // deal with leap years
        daynum = years-- * 365;
        daynum += floor(years / 4.0) - floor(years / 100.0) + floor( years / 400.0) + 1;

        //months
        daynum += daysInMonths((int) months, (int) (years + 1.0));

        // days
        daynum += date[2];

        // hours
        daynum += date[3] / 24.0;

        // minutes
        daynum += date[4] / 1440.0;

        // seconds
        daynum += date[5] / 86400.0;

        (*result)(i) = daynum;
		(*resultSecs)(i) = daynum*60*60*24;
    }
	outputs.push_back(result);
	
	if(nargout>1)
	{
		outputs.push_back(resultSecs);
	}
    return true;
}
//------------------------------------------------------------------------------
// Returns true if the array x sorted by mode [issorted]
//------------------------------------------------------------------------------
bool oml_issorted(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    bool rows = false;

    // these determine which to check for -- they can be both active at once
    bool sortedAscending = false;
    bool sortedDescending = false;

    if (nargin > 1)
    {
        std::string option = readString(toCurrencyStr(eval, inputs[1], false, false));

        if (option == "rows")
        {
            rows = true;
            if (nargin == 2)
                sortedAscending = true;
        }
        else if (option == "ascending")
            sortedAscending = true;
        else if (option == "descending")
            sortedDescending = true;
        else if (option == "either")
            sortedAscending = sortedDescending = true;
        else
            throw OML_Error(HW_ERROR_INVALIDOPTION(option));

        if (nargin > 2)
        {
            bool rowsSet = rows;
            bool modeSet = sortedAscending || sortedDescending;
            const Currency &input3 = inputs[2];
            option = readOption(eval, inputs[2]);
            if (option == "rows")
            {
                if (rowsSet)
                    throw OML_Error(HW_ERROR_NOTSETOPTROWTWICE);
                rows = true;
            }
            else
            {
                if (option == "ascending")
                {
                    sortedAscending = true;
                }
                else if (option == "descending")
                {
                    sortedDescending = true;
                }
                else if (option == "either")
                {
                    sortedAscending = sortedDescending = true;
                }
                else
                    throw OML_Error("Error: invalid option '" + option + "'");

                if (modeSet)
                    throw OML_Error(HW_ERROR_NOTSETDIRECTTWICE);
            }
        }
    }
    else
    {
        sortedAscending = true;
    }

    const Currency &input1 = inputs[0];
    if (input1.IsScalar() || input1.IsComplex())
    {
        outputs.push_back(getTrue());
        return true;
    }
    else if (input1.IsMatrix() || input1.IsString())
    {
        const hwMatrix *mtx = input1.Matrix();
        if (mtx->Size())
        {
            if (rows)
            {
                Currency last = readRow(eval, input1);
                for (int i = 1; i < mtx->M(); i++)
                {
                    Currency current = readRow(eval, input1, i);
                    if (sortedAscending && rowVecGreaterThan(last.Matrix(), current.Matrix()))
                    {
                        sortedAscending = false;
                        if (!sortedDescending)
                            break;
                    }
                    else if (sortedDescending && rowVecLessThan(last.Matrix(), current.Matrix()))
                    {
                        sortedDescending = false;
                        if (!sortedAscending)
                            break;
                    }
                    last = current;
                }
            }
            else // !rows
            {
                if (!mtx->IsVector())
                    throw OML_Error(HW_ERROR_INPVECSORTROW);
                if (mtx->IsReal())
                {
                    double last = (*mtx)(0);
                    for (int i = 1; i < mtx->Size(); i++)
                    {
                        double current = (*mtx)(i);
                        if (sortedAscending && last > current)
                        {
                            sortedAscending = false;
                            if (!sortedDescending)
                                break;
                        }
                        else if (sortedDescending && last < current)
                        {
                            sortedDescending = false;
                            if (!sortedAscending)
                                break;
                        }
                        last = current;
                    }
                }
                else // mtx is complex
                {
                    hwComplex last = mtx->z(0);
                    for (int i = 1; i < mtx->Size(); i++)
                    {
                        hwComplex current = mtx->z(i);
                        if (sortedAscending && complexGreaterThan(last, current))
                        {
                            sortedAscending = false;
                            if (!sortedDescending)
                                break;
                        }
                        else if (sortedDescending && complexLessThan(last, current))
                        {
                            sortedDescending = false;
                            if (!sortedAscending)
                                break;
                        }
                        last = current;
                    }
                }
            }
        }
        else
            throw OML_Error(HW_MATH_MSG_EMPTYMATRIX);
    }
    else if (input1.IsCellArray())
    {
        HML_CELLARRAY *cell = input1.CellArray();
        if (cell->Size())
        {
            if (rows)
            {
                std::string last = concatRowToString(eval, 0, cell);

                for (int i = 1; i < cell->M(); i++)
                {
                    std::string current = concatRowToString(eval, i, cell);

                    // don't break because we still need to verify each element is a string
                    if (sortedAscending && last.compare(current) > 0)
                        sortedAscending = false;
                    else if (sortedDescending && last.compare(current) < 0)
                        sortedDescending = false;

                    last = current;
                }
            }
            else // !rows
            {
                if (!cell->IsVector())
                    throw OML_Error(HW_ERROR_INPVECSORTROW);

                std::string last = readString((*cell)(0));
                for (int i = 1; i < cell->Size(); i++)
                {
                    const Currency &currentCur = (*cell)(i);
                    if (!currentCur.IsString())
                        throw OML_Error(HW_ERROR_CELLELEMSTR);
                    std::string current = readString(currentCur);

                    // don't break because we still need to verify each element is a string
                    if (sortedAscending && last.compare(current) > 0)
                        sortedAscending = false;
                    else if (sortedDescending && last.compare(current) < 0)
                        sortedDescending = false;

                    last = current;
                }
            }
        }
        else
            throw OML_Error(HW_MATH_MSG_EMPTYMATRIX);
    }
    else
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    outputs.push_back(HW_MAKE_BOOL_CURRENCY(sortedAscending || sortedDescending));
    return true;
}
//------------------------------------------------------------------------------
// Extract a field value from a structure [getfield]
//------------------------------------------------------------------------------
bool oml_getfield(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() < 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency nextOutput = inputs[0];
    nextOutput.SetOutputName(std::string());
    std::pair<int, int> nextIndex(0, 0);
    bool indexDefined = false;

    std::vector<Currency>::const_iterator iter = inputs.cbegin() + 1;
    while (iter != inputs.cend())
    {
        const Currency &tocheck = *iter++;
        if (tocheck.IsCellArray())
        {
            std::pair<int, int> tempIndex;

            HML_CELLARRAY *cell = tocheck.CellArray();
            bool indexWasDefined = indexDefined;
            indexDefined = nextOutput.IsStruct();

            switch (cell->Size())
            {
                case 0:
                    tempIndex.first = tempIndex.second = -1;
                    break;
                case 1:
                    if (!(*cell)(0).IsPositiveInteger())
                        throw OML_Error(HW_ERROR_INDEXPOSINT);
                    tempIndex.first = (int)((*cell)(0).Scalar()) - 1;
                    tempIndex.second = -1;
                    if (!indexDefined)
                        nextOutput = getAtIndex(nextOutput, tempIndex.first);
                    break;
                case 2:
                    if (!((*cell)(0).IsPositiveInteger() && (*cell)(1).IsPositiveInteger()))
                        throw OML_Error(HW_ERROR_INDEXPOSINT);
                    tempIndex.first = (int) (*cell)(0).Scalar() - 1;
                    tempIndex.second = (int) (*cell)(1).Scalar() - 1;
                    if (!indexDefined)
                        nextOutput = getAtIndex(nextOutput, tempIndex.first, tempIndex.second);
                    break;
                default:
                    throw OML_Error(OML_ERR_UNSUPPORTDIM, 1);
            }

            if (indexWasDefined)
            {
                StructData *temp = nextOutput.Struct();
                if (tempIndex.first == -1)
                {
                    nextIndex = tempIndex;
                }
                else if (temp->M() * temp->N())
                {
                    if (!(tempIndex.first == 0 && tempIndex.second < 1))
                    {
                        throw OML_Error(HW_MATH_MSG_INVALIDINDEX);
                    }
                }
                else
                {
                    throw OML_Error(HW_MATH_MSG_INVALIDINDEX);
                }
            }
            else
            {
                nextIndex = tempIndex;
            }
        }
        else if (tocheck.IsString())
        {
            if (!indexDefined)
            {
                if (nextOutput.IsStruct() && nextOutput.Struct()->M() * nextOutput.Struct()->N() != 1)
                    throw OML_Error(HW_ERROR_INDEXSTRUCTFIELDS);
                nextIndex.first = nextIndex.second = 0;
            }

            if (!nextOutput.IsStruct())
                throw OML_Error(HW_ERROR_VALNOTSTRUCT);

            StructData *sd = nextOutput.Struct();

			if (nextIndex.second != -1)
			{
				if (nextIndex.first >= sd->M() || nextIndex.second >= sd->N())
					throw OML_Error(HW_ERROR_INDGRDIM);
			}
			else
			{
				if (nextIndex.first >= sd->Size())
					throw OML_Error(HW_ERROR_INDGRDIM);
			}

            std::string field = readString(tocheck);
            if (isField(eval, sd->GetFieldNames(), field))
                nextOutput = sd->GetValue(nextIndex.first, nextIndex.second, field);
            else
                throw OML_Error("Error: struct does not have a member '" + field + "'");

            indexDefined = false;
        }
        else
            throw OML_Error(HW_ERROR_INPSMUSTCELLORSTR);
    }

    if (indexDefined)
    {
        if (nextIndex.first == -1)
        {
            if (nextOutput.IsStruct())
            {
                // convert to struct with same fields but zero elements
                StructData *sd = new StructData(*nextOutput.Struct());
                sd->DimensionNew(0, 0);
                nextOutput = sd;
            }
        }
        else
        {
            Currency temp;
            if (nextIndex.second == -1)
            {
                temp = getAtIndex(nextOutput, nextIndex.first);
            }
            else
            {
                temp = getAtIndex(nextOutput, nextIndex.first, nextIndex.second);
            }

            if (nextOutput.IsCellArray())
            {
                HML_CELLARRAY *cell = EvaluatorInterface::allocateCellArray(1, 1);
                (*cell)(0) = temp;
                nextOutput = cell;
            }
            else
            {
                nextOutput = temp;
            }
        }

    }

    outputs.push_back(nextOutput);
    return true;
}
//------------------------------------------------------------------------------
// Returns a logical matrix first input is a member of second input [ismember]
//------------------------------------------------------------------------------
bool oml_ismember(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> searchfor, searchin;

    const Currency &input1 = inputs[0];
    const Currency &input2 = inputs[1];

    // push all elements of input1 into searchfor, and all elements of input2 into searchin
    // the strategy generates a single search loop, but is inefficient
    int m = 0;
    int n = 0;

    if (nargin > 2)
    {
        Currency input3 = unnest(inputs[2], "Error: option must be a string");

        if (input3.IsString())
        {
            if (readString(input3) != "rows")
                throw OML_Error(HW_ERROR_INVOPTMUSTROW);
        }
        else
        {
            throw OML_Error(HW_ERROR_OPTIONSTRING);
        }
		if ((input1.IsString() && input2.IsString()) || (input1.IsCellArray() && input2.IsCellArray()) || (input1.IsString() && input2.IsCellArray()) || (input1.IsCellArray() && input2.IsString()))
		{
			stringVecFromCurrencyRows(eval, input1, input2, searchfor, searchin, &m, &n);
		}
		else if (input1.IsMatrix() && input2.IsMatrix())
		{
			const hwMatrix *in1 = input1.Matrix();
			const hwMatrix *in2 = input2.Matrix();

			m = in1->M();
			n = in1->N();

			hwMatrix *boolR = EvaluatorInterface::allocateMatrix(1, n, 0.0);
			hwMatrix *idxR  = EvaluatorInterface::allocateMatrix(1, n, 0.0);

			if (!(in1->N() == in2->N()))
				throw OML_Error(HW_ERROR_SAMENUMOFCOLWCOMPROW);
			
			for (int i = 0; i < m; i++)
			{
				for (int j = 0; j < in2->M(); j++)
				{
					bool matches = true;
					for (int k = 0; k<n; k++)
					{
						if((*in1)(i,k) != (*in2)(j,k))
						{
							matches = false;
							break;
						}
					}
					if (matches)
					{
						(*boolR)((const int)i) = 1.0;
						(*idxR)((const int)i) = (double) j;
						break;
					}
				}
			}

			Currency boolcur(boolR), idxcur(idxR);
            boolcur.SetMask(Currency::MASK_LOGICAL);
			outputs.push_back(boolcur);
			outputs.push_back(idxcur);
			return true;
		}
    }
	else if ((input1.IsString() && input2.IsString()) || (input1.IsCellArray() && input2.IsCellArray()) || (input1.IsString() && input2.IsCellArray()) || (input1.IsCellArray() && input2.IsString()))
	{
	    stringVecFromCurrency(eval, input1, input2, searchfor, searchin, &m, &n);
	}
    else
    {
        if (input1.IsScalar())
        {
            searchfor.push_back(input1.Scalar());
        }
        else if (input1.IsComplex())
        {
            searchfor.push_back(input1.Complex());
        }
        else if (input1.IsMatrix())
        {
            const hwMatrix* in1 = input1.Matrix();

            for (int i=0; i < in1->Size(); i++)
            {
                if (in1->IsReal())
                    searchfor.push_back((*in1)(i));
                else
                    searchfor.push_back(in1->z(i));
            }
        }
        else if (input1.IsNDMatrix())
        {
            const hwMatrixN* in1 = input1.MatrixN();

            for (int i=0; i < in1->Size(); i++)
            {
                if (in1->IsReal())
                    searchfor.push_back((*in1)(i));
                else
                    searchfor.push_back(in1->z(i));
            }
        }
        else
        {
            throw OML_Error(HW_ERROR_INPUTSTRCELLMTX);
        }

        if (input2.IsMatrix())
        {
            const hwMatrix* in2 = input2.Matrix();

            for (int j=0; j < in2->Size(); j++)
            {
                if (in2->IsReal())
                    searchin.push_back((*in2)(j));
                else
                    searchin.push_back(in2->z(j));
            }
        }
        else if (input2.IsScalar())
        {
            searchin.push_back(input2.Scalar());
        }
        else if (input2.IsComplex())
        {
            searchin.push_back(input2.Complex());
        }
        else if (input2.IsNDMatrix())
        {
            const hwMatrixN* in2 = input2.MatrixN();

            for (int j=0; j < in2->Size(); j++)
            {
                if (in2->IsReal())
                    searchin.push_back((*in2)(j));
                else
                    searchin.push_back(in2->z(j));
            }
        }
        else
        {
            throw OML_Error(HW_ERROR_INPUTSTRCELLMTX);
        }
	}

    // create ouput objects
    if (!input1.IsNDMatrix())       // 2D output
    {
        if (!input1.IsString() && !input1.IsCellArray())
        {
            if (input1.IsMatrix())
            {
                m = input1.Matrix()->M();
                n = input1.Matrix()->N();
            }
            else
            {
                m = 1;
                n = 1;
            }
        }

        hwMatrix *boolResult = EvaluatorInterface::allocateMatrix(m, n, 0.0);
        hwMatrix *idxResult  = EvaluatorInterface::allocateMatrix(m, n, 0.0);

        // do the search
        for (size_t i = 0; i < searchfor.size(); i++)
        {
            const Currency &tofind = searchfor[i];
            for (size_t j = searchin.size(); j; j--) // find last instance of
            {
                if (isequal(tofind, searchin[j-1]))
                {
                    (*boolResult)((const int)i) = 1.0;
                    (*idxResult)((const int)i) = (double) j;
                    break;
                }
            }
        }

        Currency boolcur(boolResult), idxcur(idxResult);
        boolcur.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(boolcur);
        outputs.push_back(idxcur);
    }
    else if (input1.IsNDMatrix())   // ND output
    {
        const hwMatrixN* in1 = input1.MatrixN();
        const std::vector<int>& dims = in1->Dimensions();
        hwMatrixN *boolResult = EvaluatorInterface::allocateMatrixN();
        hwMatrixN *idxResult  = EvaluatorInterface::allocateMatrixN();
        boolResult->Dimension(dims, hwMatrixN::REAL);
        idxResult->Dimension(dims, hwMatrixN::REAL);
        boolResult->SetElements(0.0);
        idxResult->SetElements(0.0);

        // do the search
        for (size_t i = 0; i < searchfor.size(); i++)
        {
            const Currency &tofind = searchfor[i];
            for (size_t j = searchin.size(); j; j--) // find last instance of
            {
                if (isequal(tofind, searchin[j-1]))
                {
                    (*boolResult)((const int)i) = 1.0;
                    (*idxResult)((const int)i) = (double) j;
                    break;
                }
            }
        }

        Currency boolcur(boolResult), idxcur(idxResult);
        boolcur.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(boolcur);
        outputs.push_back(idxcur);
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSTRCELLMTX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the phase angle of input [angle]
//------------------------------------------------------------------------------
bool oml_angle(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];
    if (input.IsScalar())
    {
        outputs.push_back(hwComplex(input.Scalar(), 0.0).Arg());
    }
    else if (input.IsComplex())
    {
        outputs.push_back(input.Complex().Arg());
    }
    else if (input.IsMatrix() || input.IsString())
    {
        const hwMatrix *mtx = input.Matrix();
        hwMatrix *out = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

        if (mtx->IsReal())
        {
            for (int i = 0; i < out->Size(); i++)
            {
                (*out)(i) = hwComplex((*mtx)(i), 0.0).Arg();
            }
        }
        else
        {
            for (int i = 0; i < out->Size(); i++)
            {
                (*out)(i) = mtx->z(i).Arg();
            }
        }
        outputs.push_back(out);
    }
    else if (input.IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_angle);
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true if given structure has given field(s) [isfield]
//------------------------------------------------------------------------------
bool oml_isfield(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &in1 = inputs[0];
    const Currency &in2 = inputs[1];

    if (!in1.IsStruct())
    {
        Currency val(0.0);
        val.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(val);
        return true;
    }

    StructData *sd = in1.Struct();
    std::map<std::string, int> fieldNames = sd->GetFieldNames();

    if (in2.IsCellArray())
    {
        HML_CELLARRAY *cell = in2.CellArray();
        hwMatrix *result = EvaluatorInterface::allocateMatrix(cell->M(), cell->N(), hwMatrix::REAL);

        for (int i = 0; i < cell->Size(); i++)
        {
            (*result)(i) = isField(eval, fieldNames, (*cell)(i));
        }

        Currency val(result);
        val.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(val);
    }
    else
    {
        Currency val(isField(eval, fieldNames, in2));
        val.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(val);
    }
    return true;
}
//------------------------------------------------------------------------------
// Gives a brief description on the origin of a function [which]
//------------------------------------------------------------------------------
bool oml_which(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    int nargout = getNumOutputs(eval);

    for (int i = 0; i < inputs.size(); i++)
    {
        if (!inputs[i].IsString())
            throw OML_Error(OML_ERR_STRING, i+1, OML_VAR_TYPE);
    }

    std::vector<std::string> todisp;

    for (int i = 0; i < inputs.size(); i++)
    {
        std::string inpstr = readString(inputs[i]);
        std::string topush = "'" + inpstr + "'";
        Currency val = eval.GetValue(inpstr);

        if (eval.IsUserFunction(inpstr))
        {
            todisp.push_back(topush + " is a user function");
            if (nargout)
                outputs.push_back(std::string());
        }
        else if (eval.IsStdFunction(inpstr))
        {
            todisp.push_back(topush + " is a built-in function");
            if (nargout)
                outputs.push_back(std::string());
        }
        else
        {
            // search for file
            std::string filename;
            bool addext = inpstr.rfind(".") == std::string::npos;
            if (addext)
            {
                if (eval.FindFileInPath(inpstr + ".oml", filename) || eval.FindFileInPath(inpstr + ".hml", filename) || eval.FindFileInPath(inpstr + ".m", filename))
                {
                    todisp.push_back(topush + " is a function defined from " + filename);
                    if (nargout)
                        outputs.push_back(filename);
                }
                else if (nargout)
                {
                    outputs.push_back(std::string());
                }
            }
            else
            {
                if (eval.FindFileInPath(inpstr, filename))
                {
                    todisp.push_back(topush + " is a file: " + filename);
                    if (nargout)
                        outputs.push_back(filename);
                }
                else if (nargout)
                {
                    outputs.push_back(std::string());
                }
            }

            topush += " is a function defined from " + filename + "\n";
        }
    }

    if (!nargout)
    {
        std::vector<std::string>::iterator iter;
        for (iter = todisp.begin(); iter != todisp.end(); iter++)
            eval.PrintResult(*iter + '\n');
    }

    return true;
}
//------------------------------------------------------------------------------
// Changes current working directory [cd]
//------------------------------------------------------------------------------
bool oml_cd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    if (nargin == 1)
    {
        std::string dir = getAbsolutePath(eval, inputs[0]);
        cd(dir, eval);
        eval.OnChangeDir(dir); //broadcast change in current dir
    }

    if (getNumOutputs(eval) || !nargin)
    {
        outputs.push_back(BuiltInFuncsUtils::GetCurrentWorkingDir());
    }
    
    if (nargin > 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template <typename T>
static hwMatrix* dolinspace(T start, T end, int n)
{
    T step = (end - start) / (n - 1);
    T current = start;

    hwMatrix *result = EvaluatorInterface::allocateMatrix (1, n, hwMatrix::REAL);

    int i;
    result->SetElement(0, start);
    for (i = 1; i < result->Size() - 1; i++)
        result->SetElement(i, (start += step));
    result->SetElement(i, end);

    return result;
}
//------------------------------------------------------------------------------
// Create a row vector with equally spaced elements [linspace]
//------------------------------------------------------------------------------
bool oml_linspace(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &in1 = inputs[0];
    const Currency &in2 = inputs[1];
    int numElems;

    if (nargin > 2)
    {
        if (!IsInteger(inputs[2].Scalar()).IsOk())
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_PARAMETER);

        numElems = (int) inputs[2].Scalar();
    }
    else
    {
        numElems = 100;
    }

    if (numElems < 2)
    {
        if ((in1.IsScalar() || in1.IsComplex() || in1.IsMatrix()) &&
            (in2.IsScalar() || in2.IsComplex() || in2.IsMatrix()))
        {
            outputs.push_back(in2);
            return true;
        }
    }

    bool complex = false;
    double s1 = 0;
	double s2 = 0;
    hwComplex c1;
	hwComplex c2;
    hwMatrix *m1 = nullptr;

    if (in1.IsScalar())
    {
        s1 = in1.Scalar();
    }
    else if (in1.IsComplex())
    {
        complex = true;
        c1 = in1.Complex();
    }
    else if (in1.IsMatrix())
    {
        m1 = EvaluatorInterface::allocateMatrix(in1.Matrix());
        complex = !m1->IsReal();
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

    if (in2.IsScalar())
    {
        if (complex)
            c2 = hwComplex(in2.Scalar(), 0);
        else
            s2 = in2.Scalar();
    }
    else if (in2.IsComplex())
    {
        c2 = in2.Complex();

        if (!complex)
        {
            complex = true;
            c1 = hwComplex(s1, 0);
            if (m1 != nullptr)
                BuiltInFuncsUtils::CheckMathStatus(eval, m1->MakeComplex());
        }
    }
    else if (in2.IsMatrix())
    {
        if (m1 == nullptr)
        {
            m1 = EvaluatorInterface::allocateMatrix(in2.Matrix());
            if (complex)
            {
                if (m1->IsReal())
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, m1->MakeComplex());
                }
            }
            else
            {
                if (!m1->IsReal())
                {
                    complex = true;
                    c1 = hwComplex(s1, 0);
                }
            }
        }
        else
        {
            hwMatrix *m2 = EvaluatorInterface::allocateMatrix(in2.Matrix());

            if (m1->Size() != m2->Size())
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            complex = checkMakeComplex(eval, m1, m2);
            hwMatrix *result = EvaluatorInterface::allocateMatrix(m1->Size(), numElems, m1->Type());
            Currency out(result);
            if (complex)
            {
                for (int i = 0; i < result->M(); i++)
                {
                    hwMatrix *row = dolinspace(m1->z(i), m2->z(i), numElems);
                    Currency temp(row);
                    BuiltInFuncsUtils::CheckMathStatus(eval, result->WriteRow(i, *row));
                }
            }
            else
            {
                for (int i = 0; i < result->M(); i++)
                {
                    hwMatrix *row = dolinspace((*m1)(i), (*m2)(i), numElems);
                    Currency temp(row);
                    BuiltInFuncsUtils::CheckMathStatus(eval, result->WriteRow(i, *row));
                }
            }
            outputs.push_back(out);
            return true;
        }
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

    if (m1 == nullptr)
    {
        if (complex)
            outputs.push_back(dolinspace(c1, c2, numElems));
        else
            outputs.push_back(dolinspace(s1, s2, numElems));
        return true;
    }
    else
    {
        if (in1.IsMatrix())
        {
            c1 = c2;
            s1 = s2;
        }

        hwMatrix *result = EvaluatorInterface::allocateMatrix(m1->Size(), numElems, m1->Type());
        Currency out(result);
        hwMatrix *row = nullptr;
        Currency temp(row);
        if (complex)
        {
            for (int i = 0; i < result->M(); i++)
            {
                if (in1.IsMatrix())
                    row = dolinspace(m1->z(i), c1, numElems);
                else
                    row = dolinspace(c1, m1->z(i), numElems);
                BuiltInFuncsUtils::CheckMathStatus(eval, result->WriteRow(i, *row));
            }
        }
        else
        {
            for (int i = 0; i < result->M(); i++)
            {
                if (in1.IsMatrix())
                    row = dolinspace((*m1)(i), s1, numElems);
                else
                    row = dolinspace(s1, (*m1)(i), numElems);
                BuiltInFuncsUtils::CheckMathStatus(eval, result->WriteRow(i, *row));
            }
        }
        outputs.push_back(out);
        return true;
    }
}
//------------------------------------------------------------------------------
// Returns true if input is a global variable [isglobal]
//------------------------------------------------------------------------------
bool oml_isglobal(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency input = inputs[0];

    while (input.IsCellArray())
    {
        HML_CELLARRAY *cell = input.CellArray();
        if (cell->Size())
        {
            input = (*cell)(0);
        }
        else
        {
            outputs.push_back(getFalse());
            return true;
        }
    }

    std::string varname = readString(toCurrencyStr(eval, input, false, false));
    Currency out(eval.IsGlobal(varname));
    out.SetMask(Currency::MASK_LOGICAL);
    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Gets file stream type
//------------------------------------------------------------------------------
static FileStreamType determineFileStreamType(FILE*& file, bool rerouteBuiltin = true)
{
    if (file == stdout)
    {
        if (rerouteBuiltin)
            file = makeTempFile();
        return Stdout;
    }
    else if (file == stderr)
    {
        if (rerouteBuiltin)
            file = makeTempFile();
        return Stderr;
    }
    else
    {
        return Other;
    }
}
//------------------------------------------------------------------------------
// Returns true if successul in printing to a file stream [fprintf]
//------------------------------------------------------------------------------
bool oml_fprintf(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)
{
    int numinputs = inputs.empty() ? 0 : static_cast<int>(inputs.size());

    if (numinputs < 1) throw OML_Error(OML_ERR_NUMARGIN);

    std::FILE* file;
    int fileid = getFileFromInput(eval, inputs[0]);
    if (fileid == -1)
        file = stdout;
    
    else
    {
        if (numinputs == 1) throw OML_Error(HW_ERROR_MISOUTPUTTEMPLATE);

        checkFileIndex(eval, fileid, true);
        file = eval.GetFile(fileid);
    }

    FileStreamType ftype = determineFileStreamType(file, false);

    std::string result (sprintf(eval, inputs.cbegin() + (fileid == -1 ? 0 : 1), inputs.cend()));   

    if (ftype == Stdout)
    {
        Currency out(result);
        out.PrintfOutput(); // Set flag that this is coming from printf
        eval.PrintResult(out);
    }
    else if (ftype == Stderr) throw OML_Error(result.c_str());
    else if (ftype == Other)
        fputs(result.c_str(), file);

    if (getNumOutputs(eval))
        outputs.push_back(result.length());

    return true;
}
//------------------------------------------------------------------------------
// Returns true if successul in printing [printf]
//------------------------------------------------------------------------------
bool oml_printf(EvaluatorInterface           eval, 
                const std::vector<Currency>& inputs,
                std::vector<Currency>&       outputs)
{
    std::string result (sprintf(eval, inputs.cbegin(), inputs.cend()));

    Currency out(result);
    out.PrintfOutput(); // Set flag that this is coming from printf
    eval.PrintResult(out);

    if (eval.GetNargoutValue() > 0)
    {
        outputs.push_back(result.length());
    }
    return true;
}
//------------------------------------------------------------------------------
//! Returns true after returning string to double [str2double]
//------------------------------------------------------------------------------
bool oml_str2double(EvaluatorInterface           eval, 
                    const std::vector<Currency>& inputs, 
                    std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency input    = inputs[0];
    bool     isString = input.IsString();
    if (!isString && !input.IsCellArray())
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 1, OML_VAR_TYPE);

    std::vector<Currency> loc_inputs;
    int m = 0;
    int n = 0;
    if (isString)
        loc_inputs = inputs;
    else
    {
        HML_CELLARRAY* cells = input.CellArray();
        assert(cells);

        m = cells->M();
        n = cells->N();

	    for (int i = 0; i < m; ++i)
	    {	
		    for (int j = 0; j < n; ++j)
		    {
                Currency cur = (*cells)(i,j);
                if (!cur.IsString())
                    throw OML_Error(OML_ERR_STRING_STRINGCELL, 1, OML_VAR_TYPE);
                loc_inputs.push_back(cur);
            }
        }
    }

    bool status  = true;
	std::vector<Currency> loc_outputs;

    for (std::vector<Currency>::const_iterator itr = loc_inputs.begin();
         itr != loc_inputs.end(); ++itr)
    {
        std::vector<Currency> in;
        std::vector<Currency> out;
        in.push_back(*itr);
	    oml_str2num(eval, in, out);
        if (out.empty()) 
            continue;
        
        loc_outputs.push_back(out[0]);
        if (out.size() > 1)
        {
            Currency tmp = out[1];
            if (tmp.IsLogical())
            {
                int val = static_cast<int>(tmp.Scalar());
                if (val <= 0)
                    status = false;
            }
        }
    }
    if (isString)
    {
	    for (std::vector<Currency>::const_iterator itr = loc_outputs.begin();
             itr != loc_outputs.end(); ++itr)
	    {
		    Currency cur = *itr;

		    if (cur.IsEmpty())
		    {
			    outputs.push_back(Currency(std::numeric_limits<double>::quiet_NaN()));
                continue;
		    }
		    else if (cur.IsCellArray())
		    {
			    HML_CELLARRAY* cells = cur.CellArray();
			    if (!cells || cells->Size() == 0)
				    outputs.push_back(Currency(std::numeric_limits<double>::quiet_NaN()));
			    else 
				    outputs.push_back(cur);
                continue;
		    }

            outputs.push_back(cur);
	    }
    }
    else
    {
        HML_CELLARRAY* out = EvaluatorInterface::allocateCellArray(m, n);
        int numout = loc_outputs.empty() ? 0 : static_cast<int>(loc_outputs.size());
        for (int j = 0, idx = 0; j < n && idx < numout; ++j)
        {
            for (int i = 0; i < m && idx < numout; ++i)
	        {
		        Currency cur = loc_outputs[idx++];

		        if (cur.IsEmpty())
                {
                    (*out)(i, j) = std::numeric_limits<double>::quiet_NaN();
                    continue;
                }
		        else if (cur.IsCellArray())
		        {
			        HML_CELLARRAY* cells = cur.CellArray();
			        if (!cells || cells->Size() == 0)
                    {
				        (*out)(i, j) = std::numeric_limits<double>::quiet_NaN();
                        continue;
                    }
                }
		        (*out)(i, j) = cur;
            }
	    }

        outputs.push_back(out);
    }

    if (eval.GetNargoutValue() > 1)
    {
        Currency tmp(status);
        tmp.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(tmp);
    }
	return true;
}
//------------------------------------------------------------------------------
// Converts string to number, by evaluating it if needed [str2num]
//------------------------------------------------------------------------------
bool oml_str2num(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];
    if (!input.IsString())
        throw OML_Error(HW_ERROR_INPUTSTRING);

    Interpreter interp(eval);
    std::string todo = readString(input);

    // trim left spaces
    size_t start = todo.find_first_not_of(' ');
    if (start != std::string::npos)
        todo.erase(todo.begin(), todo.begin() + start);

    if (todo.length() && todo[0] != '{')
        todo = '[' + todo + ']';
    Currency first = interp.DoString(todo);
    if (first.IsCellArray())
    {
		// Make sure this is a copy.  Otherwise 'first' and 'outcur' will be holding the
		// same pointer and crash.
        HML_CELLARRAY *out = EvaluatorInterface::allocateCellArray(first.CellArray());
        Currency outcur(out);
        int m = out->M(), n = out->N();
        for (int i = 1; i < input.Matrix()->M(); i++)
        {
            Currency result = interp.DoString(readString(input, i));
            if (result.IsCellArray())
            {
                HML_CELLARRAY *cell = result.CellArray();
                if (n != cell->N())
                {
                    outputs.push_back(EvaluatorInterface::allocateMatrix());
                    return true;
                }
                BuiltInFuncsUtils::CheckMathStatus(eval, out->Resize(m + 1 + cell->M(), n));
                for (int j = 0; j < cell->M(); j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        (*out)(j + m++, k) = (*cell)(j, k);
                    }
                }
            }
            else
            {
                outputs.push_back(EvaluatorInterface::allocateMatrix());
                Currency status(0.0);
                status.SetMask(Currency::MASK_LOGICAL);
                outputs.push_back(status);
                return true;
            }
        }

        outputs.push_back(outcur);
        Currency status(1.0);
        status.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(status);
        return true;
    }
    else
    {
        hwMatrix *out;
        if (first.IsScalar())
            out = EvaluatorInterface::allocateMatrix(1, 1, first.Scalar());
        else if (first.IsComplex())
            out = EvaluatorInterface::allocateMatrix(1, 1, first.Complex());
        else if (first.IsMatrix())
            out = EvaluatorInterface::allocateMatrix(first.Matrix());
        else
        {
            outputs.push_back(EvaluatorInterface::allocateMatrix());
            Currency status(0.0);
            status.SetMask(Currency::MASK_LOGICAL);
            outputs.push_back(status);
            return true;
        }

        int currentrow = out->M() - 1;
        int outn = out->N();
        Currency outcur(out);

        for (int i = 1; i < input.Matrix()->M(); i++)
        {
            std::string todo = '[' + readString(input, i) + ']';
            Currency result = interp.DoString(todo);
            if (result.IsScalar()|| result.IsComplex())
            {
                if (outn != 1)
                {
                    outputs.push_back(EvaluatorInterface::allocateMatrix());
                    return true;
                }
                BuiltInFuncsUtils::CheckMathStatus(eval, out->Resize(out->M() + 1, 1));
                if (result.IsScalar())
                    (*out)(++currentrow) = result.Scalar();
                else
                    out->z(++currentrow) = result.Complex();
            }
            else if (result.IsMatrix())
            {
                const hwMatrix *mtx = result.Matrix();
                if (outn != mtx->N())
                {
                    outputs.push_back(EvaluatorInterface::allocateMatrix());
                    return true;
                }
                BuiltInFuncsUtils::CheckMathStatus(eval, out->Resize(currentrow + 1 + mtx->M(), outn));
                for (int j = 0; j < mtx->M(); j++)
                {
                    hwMatrix *row = readRow(eval, mtx, j);
                    BuiltInFuncsUtils::CheckMathStatus(eval, out->WriteRow(++currentrow, *row));
                    delete row;
                }
            }
            else
            {
                outputs.push_back(EvaluatorInterface::allocateMatrix());
                Currency status(0.0);
                status.SetMask(Currency::MASK_LOGICAL);
                outputs.push_back(status);
                return true;
            }
        }

        outputs.push_back(outcur);
        Currency status(1.0);
        status.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(status);
        return true;
    }
}
//------------------------------------------------------------------------------
// Gets fieldnames of input structure [fieldnames]
//------------------------------------------------------------------------------
bool oml_fieldnames(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];
    if (!input.IsStruct())
        throw OML_Error(HW_ERROR_INPUTSTRUCT);

    StructData *sd = input.Struct();
    const std::map<std::string, int> fieldnames = sd->GetFieldNames();
    HML_CELLARRAY *out = EvaluatorInterface::allocateCellArray((int)fieldnames.size(), 1);

    std::map<std::string, int>::const_iterator iter;
    int cellidx = 0;
    for (iter = fieldnames.cbegin(); iter != fieldnames.cend(); iter++)
    {
        (*out)(cellidx++) = Currency(iter->first);
    }

    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Utility which returns a double during a printf/sprintf command
//------------------------------------------------------------------------------
static double getDoubleForSprintf(std::vector<Currency>::const_iterator &inputIter, int *indexInInput)
{
    if (inputIter->IsScalar())
    {
        *indexInInput = 0;
        return (*inputIter++).Scalar();
    }
    else if (inputIter->IsComplex())
    {
        *indexInInput = 0;
        return (*inputIter++).Complex().Real();
    }
    else if (inputIter->IsMatrix() || inputIter->IsString())
    {
        const hwMatrix *m = inputIter->Matrix();
        double val;
        if (m->IsReal())
        {
            val = (*m)((*indexInInput)++);
        }
        else
        {
            val = m->z((*indexInInput)++).Real();
        }

        if (*indexInInput >= m->Size())
        {
            *indexInInput = 0;
            inputIter++;
        }
        return val;
    }
    else if (inputIter->IsNDMatrix())
    {
        double           val   = 0.0;
        const hwMatrixN* m     = inputIter->MatrixN();
        int              msize = m ? m->Size() : 0;      
        if (m && m->IsReal() && *indexInInput < msize)
            val = (*m)((*indexInInput)++);
        
        else if (m && !m->IsReal() && *indexInInput < msize)
            val = m->z((*indexInInput)++).Real();
        
        if (msize == 0 || *indexInInput >= msize)
        {
            *indexInInput = 0;
            inputIter++;
        }
        return val;
    }
    throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING);
}
//------------------------------------------------------------------------------
// Helper method for sprintf
//------------------------------------------------------------------------------
static std::string getWriteStrForSprintf(EvaluatorInterface& eval, std::vector<Currency>::const_iterator &inputIter, int *indexInInput)
{
    char *towrite = NULL;
    if (inputIter->IsScalar() || inputIter->IsComplex())
    {
        double val = getDoubleForSprintf(inputIter, indexInInput);
        if (BuiltInFuncsUtils::NeedsSpecialPrintfFormat(val))
            return BuiltInFuncsUtils::GetSpecialPrintfFormat("", val);

        towrite = new char[2];
        towrite[0] = BuiltInFuncsUtils::GetValidChar(eval, val, false);
        towrite[1] = 0;
    }
    else if (inputIter->IsMatrix() || inputIter->IsString())
    {
        const hwMatrix *m = (*inputIter++).Matrix();
        towrite = new char[m->Size() - *indexInInput + 1];
        if (m->IsReal())
        {
            for (int i = *indexInInput; i < m->Size(); i++)
            {
                char temp = BuiltInFuncsUtils::GetValidChar(eval, (*m)(i), false);
                towrite[i - *indexInInput] = temp;
                if (!temp)
                    break;
            }
        }
        else
        {
            for (int i = *indexInInput; i < m->Size(); i++)
            {
                char temp = BuiltInFuncsUtils::GetValidChar(eval, m->z(i).Real(), false);
                towrite[i - *indexInInput] = temp;
                if (!temp)
                    break;
            }
        }
        towrite[m->Size() - *indexInInput] = 0;
        *indexInInput = 0;
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING);

    std::string writestr(towrite);
    delete [] towrite;
    return writestr;
}
//------------------------------------------------------------------------------
// Computes length of the string if the double is printed to a string
//------------------------------------------------------------------------------
static size_t computeLength(bool truncate, double d)
{
    if (IsNaN_T(d) || IsInf_T(d)) return 3;  // "NaN"/"Inf"
    if (IsNegInf_T(d))            return 4;  // "-Inf"

	size_t len = 1;
	
	if (abs(d) > 0)
		len = (size_t)(log(abs(d))/log(10.0)+1);

    if (!truncate)
        len += 10;

    // for octal and hexadecimal cases, increase the buffer in case of small values
    return max(len + 1, sizeof(double) * 2);
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
static int tryGetNum(const std::string &str, size_t start, size_t end)
{
	if (str[start] == '-')
		start = start+1;

    std::string sub = str.substr(start, end - start);
    if (sub.length() == std::count_if(sub.cbegin(), sub.cend(), static_cast<int (*)(int)>(isdigit)))
    {
        return atoi(sub.c_str());
    }
    return 0;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
static size_t tryGetBufferIncrement(const std::string &tmplt)
{
    if (tmplt.length() && tmplt[0] == '%')
    {
        size_t end = tmplt.find_first_of(HW_VALID_FORMATCHARS, 1);
        size_t dotLoc = tmplt.find('.', 1);
        size_t firstNumLoc = tmplt.find_first_not_of(" *", 1);
        size_t incr = 0;

        if (firstNumLoc > end || firstNumLoc == std::string::npos)
            return 0;

        if (dotLoc > end || dotLoc == std::string::npos)
        {
            // try to grab number from firstNumLoc to end
            incr += tryGetNum(tmplt, firstNumLoc, end);
        }
        else
        {
            if (firstNumLoc != dotLoc)
            {
                // try to grab number from firstNumLoc to dotLoc
                incr += tryGetNum(tmplt, firstNumLoc, dotLoc);
            }

            // try to grab number from after dotLoc to end
            incr += tryGetNum(tmplt, dotLoc + 1, end);
        }
        return incr;
    }
    return 0;
}
//------------------------------------------------------------------------------
// Returns true if given currency has elements
//------------------------------------------------------------------------------
static bool hasElements(const Currency &cur)
{
    if (cur.IsScalar() || cur.IsComplex() || cur.IsString()) return true;
        
    if (cur.IsMatrix()) return (cur.Matrix()->Size() > 0);
    
    if (cur.IsCellArray()) return (cur.CellArray()->Size() > 0);
    
    if (cur.IsStruct()) return (cur.Struct()->M() && cur.Struct()->N());

    if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtx = cur.MatrixN();
        int   matsize        = mtx ? mtx->Size() : 0;
        return (matsize > 0);
    }
    return false;
}
//------------------------------------------------------------------------------
// Helper method for sprintf
//------------------------------------------------------------------------------
static std::string dosprintf(size_t size, const std::string &tmplt, ...)
{
    if (size)
    {
        va_list vl;
        va_start(vl, tmplt);

        char *buffer = new char[size + 1];
        int numwritten = vsprintf(buffer, tmplt.c_str(), vl);
        std::string rv(buffer);
        delete [] buffer;

        va_end(vl);

        if (numwritten < 0)
            throw OML_Error(HW_ERROR_INVALIDFORMAT);

        return rv;
    }
    return std::string();
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
static std::string dosprintf(EvaluatorInterface& eval, const std::string &tmplt, int *indexInInput, FormatType ft,
    std::vector<Currency>::const_iterator &inputIter, long long dim1, long long dim2)
{
    size_t bufferSize = tmplt.length() + dim1 + dim2 + tryGetBufferIncrement(tmplt);
    std::string writestr;
    double value;
    switch (ft)
    {
        case String:
            writestr = getWriteStrForSprintf(eval, inputIter, indexInInput);
            bufferSize += writestr.length();
            return dosprintf(bufferSize, tmplt, dim1, dim2, writestr.c_str());
        case Percent:
            return dosprintf(bufferSize, tmplt, dim1, dim2);
        case Truncate:
            value = getDoubleForSprintf(inputIter, indexInInput);   

            bufferSize += computeLength(true, value);
            if (BuiltInFuncsUtils::NeedsSpecialPrintfFormat(value))
                return BuiltInFuncsUtils::GetSpecialPrintfFormat(tmplt, value);

            return dosprintf(bufferSize, tmplt, dim1, dim2, (long long) value);
        case NoTruncate:
            value = getDoubleForSprintf(inputIter, indexInInput);
            bufferSize += computeLength(false, value);

            if (BuiltInFuncsUtils::NeedsSpecialPrintfFormat(value))
                return BuiltInFuncsUtils::GetSpecialPrintfFormat(tmplt, value);
            
            return dosprintf(bufferSize, tmplt, dim1, dim2, value);
        default:
            throw OML_Error(HW_ERROR_UNSUPFORMTYPE);
    }
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
static std::string dosprintf(EvaluatorInterface& eval, const std::string &tmplt, int *indexInInput, FormatType ft,
    std::vector<Currency>::const_iterator &inputIter, long long dim1)
{
    size_t bufferSize = tmplt.length() + dim1 + tryGetBufferIncrement(tmplt);
    std::string writestr;
    double value;
    switch (ft)
    {
        case String:
            writestr = getWriteStrForSprintf(eval, inputIter, indexInInput);
            bufferSize += writestr.length();
            return dosprintf(bufferSize, tmplt, dim1, writestr.c_str());
        case Percent:
            return dosprintf(bufferSize, tmplt, dim1);
        case Truncate:
            value = getDoubleForSprintf(inputIter, indexInInput);
            bufferSize += computeLength(true, value);
            if (BuiltInFuncsUtils::NeedsSpecialPrintfFormat(value))
                return BuiltInFuncsUtils::GetSpecialPrintfFormat(tmplt, value);
           
            return dosprintf(bufferSize, tmplt, dim1, (long long) value);
        case NoTruncate:
            value = getDoubleForSprintf(inputIter, indexInInput);
            bufferSize += computeLength(false, value);
            if (BuiltInFuncsUtils::NeedsSpecialPrintfFormat(value))
                return BuiltInFuncsUtils::GetSpecialPrintfFormat(tmplt, value);
            return dosprintf(bufferSize, tmplt, dim1, value);
        default:
            throw OML_Error(HW_ERROR_UNSUPFORMTYPE);
    }
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
static std::string dosprintf(EvaluatorInterface& eval, const std::string &tmplt, int *indexInInput, FormatType ft,
    std::vector<Currency>::const_iterator &inputIter)
{
    size_t bufferSize = tmplt.length() + tryGetBufferIncrement(tmplt);
    std::string writestr;
    double value;
    switch (ft)
    {
        case String:
            writestr = getWriteStrForSprintf(eval, inputIter, indexInInput);
            bufferSize += writestr.length();
            return dosprintf(bufferSize, tmplt, writestr.c_str());
        case Percent:
            return dosprintf(bufferSize, tmplt);
        case Truncate:
            value = getDoubleForSprintf(inputIter, indexInInput);
            bufferSize += computeLength(true, value);
            if (BuiltInFuncsUtils::NeedsSpecialPrintfFormat(value))
                return BuiltInFuncsUtils::GetSpecialPrintfFormat(tmplt, value);
            return dosprintf(bufferSize, tmplt, (long long) value);
        case NoTruncate:
            value = getDoubleForSprintf(inputIter, indexInInput);
            bufferSize += computeLength(false, value);
            if (BuiltInFuncsUtils::NeedsSpecialPrintfFormat(value))
                return BuiltInFuncsUtils::GetSpecialPrintfFormat(tmplt, value);
            
            return dosprintf(bufferSize, tmplt, value);
        default:
            throw OML_Error(HW_ERROR_UNSUPFORMTYPE);
    }
}
//------------------------------------------------------------------------------
// Returns a formatted string [sprintf]
//------------------------------------------------------------------------------
bool oml_sprintf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    std::string result = sprintf(eval, inputs.cbegin(), inputs.cend());
    outputs.push_back(result);
	return true;
}
//------------------------------------------------------------------------------
// Evaluates the specified function with the given input(s) [feval]
//------------------------------------------------------------------------------
bool oml_feval(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    int nargout = getNumOutputs(eval);
    bool multi = nargout > 1;

    FunctionInfo *fi = nullptr;
    FUNCPTR fptr = nullptr;

    const Currency &input1 = inputs[0];
    std::vector<Currency> funcInputs(inputs.begin() + 1, inputs.end());

    // get either fi or fptr via input1
    if (input1.IsString())
    {
        std::string funcName = readString(input1);
        if (funcName.find(' ') != std::string::npos)
            throw OML_Error(HW_ERROR_FUNCNAMENOTSPACE);

        if (multi)
        {
            if (!eval.FindFunctionByName(funcName, &fi, &fptr))
                throw OML_Error("Error: no such function '" + funcName + "'");
        }
        else
        {
            outputs.push_back(eval.CallFunction(funcName, funcInputs));
            return true;
        }
    }
    else if (input1.IsFunctionHandle())
    {
        fi = input1.FunctionHandle();
        if (!multi)
        {
            outputs.push_back(eval.CallInternalFunction(fi, funcInputs));
            return true;
        }
    }
    else
        throw OML_Error(HW_ERROR_INPUTSTRINGFUNC);

    // call the function
    if (fi)
    {
        if (!fi->IsBuiltIn())
        {
			std::vector<const std::string*> parameters = fi->Parameters();

            if (parameters.size() && (*parameters.back() == "varargin"))
            {
                int index = (int)(parameters.size()-1);
                funcInputs[index] = eval.CreateVararginCell(funcInputs, index);
                funcInputs.erase(funcInputs.begin()+index+1, funcInputs.end());
            }
        }
        outputs = eval.DoMultiReturnFunctionCall(fi, funcInputs, (int)funcInputs.size(), nargout, true);
    }
    else if (fptr)
    {
        outputs = eval.DoMultiReturnFunctionCall(fptr, funcInputs, (int)funcInputs.size(), nargout, true);
    }
    else
        throw OML_Error(HW_ERROR_BADFUNCHANDLE);

    return true;
}
//------------------------------------------------------------------------------
// Removes given directories from search path [rmpath]
//------------------------------------------------------------------------------
bool oml_rmpath(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (!nargin)
        throw OML_Error(OML_ERR_NUMARGIN);

    for (int i = 0; i < nargin ; i++)
    {
        std::vector<std::string> toremove = separatePathNames(eval, inputs[i]);

        std::vector<std::string>::reverse_iterator iter;
        for (iter = toremove.rbegin(); iter != toremove.rend(); iter++)
        {
            std::string rmv = BuiltInFuncsUtils::GetAbsolutePath(*iter);
            BuiltInFuncsUtils::StripTrailingSlash(rmv);

            if (!eval.RemovePath(rmv))
            {
                BuiltInFuncsUtils::SetWarning(eval, "Warning: '" + rmv + "' is not in path");
            }
        }
    }

    if (getNumOutputs(eval))
        outputs.push_back(getPathString(eval, pathsep));

    return true;
}
//------------------------------------------------------------------------------
// Sets or displays current search path [path]
//------------------------------------------------------------------------------
bool oml_path(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    int nargout = getNumOutputs(eval);
    char pathSep;
    if (nargout)
        pathSep = pathsep;
    else
        pathSep = '\n';

    if (nargin)
    {
        eval.ClearPath();
        for (int i = 0; i < nargin; i++)
        {
            std::vector<std::string> toadd = separatePathNames(eval, inputs[i]);
            std::vector<std::string>::iterator iter;
            for (iter = toadd.begin(); iter != toadd.end(); iter++)
			{
				std::string abs_path = *iter;
				
				if (abs_path != ".")
					abs_path = BuiltInFuncsUtils::GetAbsolutePath(*iter);

                checkAddPath(eval, abs_path, true);
			}
        }
        eval.ResetFuncSearchCache();
    }
    
    if (nargout)
    {
        outputs.push_back(getPathString(eval, pathSep));
    }
    else if (!nargin)
    {
        const std::vector<std::string> &paths = eval.GetPaths();
        eval.PrintResult(std::string("Current search path:\n"));
        std::vector<std::string>::const_iterator iter = paths.cbegin();
        while(iter != paths.cend())
        {
            eval.PrintResult(*iter++);
        }
    }
    return true;
}
//------------------------------------------------------------------------------
// Adds given path to search path [addpath]
//------------------------------------------------------------------------------
bool oml_addpath(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (!nargin)
        throw OML_Error(OML_ERR_NUMARGIN);

    int nargout = getNumOutputs(eval);
    bool appendToEnd = false;

    const Currency &last = inputs[nargin - 1];
    if (last.IsScalar() || (last.IsMatrix() && last.Matrix()->IsRealData()))
    {
        int appendloc = (int) last.Scalar();

        if (!IsInteger(last.Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, (int)nargin, OML_VAR_VALUE);

        if (appendloc == 1)
            appendToEnd = true;
        else if (appendloc)
            throw OML_Error(HW_ERROR_ADDPATHLOC);
        nargin--;
    }
    else if (last.IsString())
    {
        std::string appendloc = readString(last);
        if (appendloc == "-end")
        {
            appendToEnd = true;
            --nargin;
        }
        else if (appendloc == "-begin")
        {
            --nargin;
        }
    }
    else
        throw OML_Error(HW_ERROR_ADDPATHLOC);

    for (int i = 0; i < nargin; i++)
    {
        std::vector<std::string> toadd = separatePathNames(eval, inputs[i]);
        if (appendToEnd)
        {
            std::vector<std::string>::iterator iter, end = toadd.end();
            for (iter = toadd.begin(); iter != end; ++iter)
			{
				std::string abs_path = BuiltInFuncsUtils::GetAbsolutePath(*iter);
                checkAddPath(eval, abs_path, appendToEnd);
			}
        }
        else
        {
            std::vector<std::string>::reverse_iterator iter, end = toadd.rend();
            for (iter = toadd.rbegin(); iter != end; ++iter)
			{
				std::string abs_path = BuiltInFuncsUtils::GetAbsolutePath(*iter);
                checkAddPath(eval, abs_path, appendToEnd);
			}
        }
    }
    eval.ResetFuncSearchCache();

    if (nargout)
        outputs.push_back(getPathString(eval, pathsep));

    return true;
}
//------------------------------------------------------------------------------
// Modification to add path, which specifies path and function [addpath2]
//------------------------------------------------------------------------------
bool oml_addpath2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (!nargin)
        throw OML_Error(OML_ERR_NUMARGIN);

	Currency new_path   = inputs[0];
	Currency func_names = inputs[1];

	if (!new_path.IsString())
		throw OML_Error(OML_ERR_STRING, 1);

	if (!func_names.IsString() && !func_names.IsCellArray())
		throw OML_Error(OML_ERR_CELLSTRING, 2);

	std::vector<std::string> funcs;

	if (func_names.IsString())
	{
		funcs.push_back(func_names.StringVal());
	}
	else if (func_names.IsCellArray())
	{
		HML_CELLARRAY* cells = func_names.CellArray();
		int size = cells->Size();

		for (int j=0; j<size; j++)
		{
			Currency temp = (*cells)(j);

			if (temp.IsString())
				funcs.push_back(temp.StringVal());
			else
				throw OML_Error(HW_ERROR_CELLINPMUSTSTR);
		}
	}

	eval.AddPath2(new_path.StringVal(), funcs);
	
	return true;
}
//------------------------------------------------------------------------------
// Replicates an input to create a block matrix [repmat]
//------------------------------------------------------------------------------
bool oml_repmat(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];
    const Currency& input2 = inputs[1];

    // determine object type
    bool createND = false;

    if (nargin > 3)
    {
        createND = true;
    }
    else if (input1.IsNDMatrix())
    {
        createND = true;
    }
    else if (input2.IsMatrix())
    {
        if (!input2.Matrix()->IsEmptyOrVector())
            throw OML_Error(OML_ERR_VECTOR, 2, OML_VAR_DATA);

        if (input2.Matrix()->Size() > 2)
            createND = true;
    }

    // construct matrix
    if (createND)
    {
        // count repetitions
        std::vector<int> reps;

        if (nargin > 2)
        {
            for (int i = 1; i < nargin; ++i)
            {
                if (!inputs[i].IsPositiveInteger())
                    throw OML_Error(OML_ERR_NATURALNUM, i+1, OML_VAR_DATA);

                reps.push_back(static_cast<int> (inputs[i].Scalar()));
            }
        }
        else if (input2.IsScalar())
        {
            if (!input2.IsPositiveInteger())
                throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DATA);

            reps.push_back(static_cast<int> (input2.Scalar()));
            reps.push_back(static_cast<int> (input2.Scalar()));
        }
        else if (input2.IsMatrix())
        {
            const hwMatrix* mat = input2.Matrix();
            int numDim = mat->Size();

            if (!mat->IsReal())
                throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DATA);

            for (int i = 0; i < numDim; ++i)
            {
                if (!IsInteger((*mat)(i)).IsOk())
                    throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DATA);

                reps.push_back(static_cast<int> ((*mat)(i)));
            }
        }
        else
        {
            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DATA);
        }

        // get output matrix dimensions and type
        std::vector<int> dims;
        hwMatrixN::DataType type;
        hwMatrixN base;

        if (input1.IsMatrix() || input1.IsScalar() || input1.IsComplex())
        {
            const hwMatrix* mat = input1.ConvertToMatrix();
            dims.push_back(mat->M() * reps[0]);
            dims.push_back(mat->N() * reps[1]);

            for (int i = 2; i < reps.size(); ++i)
                dims.push_back(reps[i]);

            if (mat->IsReal())
                type = hwMatrixN::REAL;
            else
                type = hwMatrixN::COMPLEX;

            base.Convert2DtoND(*mat);
        }
        else if (input1.IsNDMatrix())
        {
            const std::vector<int>& dim1 = input1.MatrixN()->Dimensions();
            int min = _min(static_cast<int> (reps.size()), static_cast<int> (dim1.size()));

            for (int i = 0; i < min; ++i)
                dims.push_back(dim1[i] * reps[i]);

            for (int i = min; i < dim1.size(); ++i)
                dims.push_back(dim1[i]);

            for (int i = min; i < reps.size(); ++i)
                dims.push_back(reps[i]);

            type = input1.MatrixN()->Type();

            base = (*input1.MatrixN());
        }
        else
        {
            throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
        }

        hwMatrixN* result = EvaluatorInterface::allocateMatrixN();
        result->Dimension(dims, type);

        // populate result
        int numreps = 1;

        for (int i = 0; i < reps.size(); ++i)
            numreps *= reps[i];

        const std::vector<int>& dim1 = base.Dimensions();
        std::vector<int> block_index(dims.size());
        std::vector<int> result_index(dims.size());

        for (int i = 0; i < numreps; ++i)
        {
            std::vector<int> base_index(dim1.size());

            // copy base matrix
            for (int ii = 0; ii < base.Size(); ++ii)
            {
                if (result->IsReal())
                    (*result)(result_index) = base(ii);
                else
                    (*result).z(result_index) = base.z(ii);

                // advance base and result indices
                for (int jj = 0; jj < dim1.size(); ++jj)
                {
                    // increment index j if possible
                    if (base_index[jj] < dim1[jj]-1)
                    {
                        ++base_index[jj];
                        ++result_index[jj];
                        break;
                    }

                    // index j is maxed out, so reset and continue to j+1
                    base_index[jj] = 0;
                    result_index[jj] = block_index[jj];
                }
            }

            // advance block and result indices
            for (int j = 0; j < block_index.size(); ++j)
            {
                // increment index j if possible
                if (j < dim1.size())
                {
                    if (block_index[j] + dim1[j] < dims[j])
                    {
                        block_index[j] += dim1[j];
                        result_index[j] = block_index[j];
                        break;
                    }
                }
                else
                {
                    if (block_index[j] < dims[j]-1)
                    {
                        ++block_index[j];
                        result_index[j] = block_index[j];
                        break;
                    }
                }

                // index j is maxed out, so reset and continue to j+1
                block_index[j] = 0;
                result_index[j] = 0;
            }
        }

        outputs.push_back(result);

        return true;
    }

    // 2D matrix
    std::unique_ptr<const HML_CELLARRAY> cell;
    const hwMatrix* mtx;
    int m, n;

    if (input1.IsCellArray())
    {
        cell.reset(EvaluatorInterface::allocateCellArray(input1.CellArray()));
        m = cell->M();
        n = cell->N();
    }
	else if (input1.IsStruct())
	{
		StructData* new_sd = new StructData(*input1.Struct());
		int         new_m = 1;
		int         new_n = 1;

		if (input2.IsScalar())
			new_m = (int)input2.Scalar();

		if (inputs.size() == 3)
		{
			if (inputs[2].IsScalar())
				new_n = (int)inputs[2].Scalar();
		}

		new_sd->DimensionNew(new_m, new_n);
		outputs.push_back(new_sd);
		return true;
	}
    else if (input1.IsMatrix() || input1.IsScalar() || input1.IsComplex() || input1.IsString())
    {
        mtx = input1.ConvertToMatrix();
        if (!mtx)
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING);
        m = mtx->M();
        n = mtx->N();
    }
    else
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
    }

    if (nargin == 2)
    {
        if (input2.IsScalar() || input2.IsComplex())
        {
            int mult = (int) input2.Scalar();

            if (!IsInteger(input2.Scalar()).IsOk())
                throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_DIM);

            m *= mult;
            n *= mult;
        }
        else if (input2.IsMatrix())
        {
            const hwMatrix *dimMtx = input2.Matrix();
            int size = dimMtx->Size();
            if (dimMtx->IsReal())
            {
                if (size == 2)
                {
                    m *= (int) (*dimMtx)(0);
                    n *= (int) (*dimMtx)(1);
                }
                else if (size > 2)
                    throw OML_Error(OML_ERR_UNSUPPORTDIM, 2);
            }
            else
            {
                if (size == 2)
                {
                    m *= (int) dimMtx->z(0).Real();
                    n *= (int) dimMtx->z(1).Real();
                }
                else if (size > 2)
                    throw OML_Error(OML_ERR_UNSUPPORTDIM, 2);
            }
        }
        else
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }
    else
    {
        const Currency& input3 = inputs[2];
        if (input2.IsScalar())
            m *= (int)input2.Scalar();
        else if (input2.IsComplex())
            m *= (int)input2.Complex().Real();
        else
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DIM);

        if (input3.IsScalar())
            n *= (int)input3.Scalar();
        else if (input3.IsComplex())
            n *= (int)input3.Complex().Real();
        else
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);
    }

    if (m < 0)
        m = 0;
    if (n < 0)
        n = 0;

    if (cell == nullptr)
    {
        hwMatrix *result = EvaluatorInterface::allocateMatrix(m, n, mtx->Type());
        int M = mtx->M();
        int N = mtx->N();

        if (result->IsReal())
        {
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    (*result)(i, j) = (*mtx)(i % M, j % N);
                }
            }
        }
        else
        {
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    result->z(i, j) = mtx->z(i % M, j % N);
                }
            }
        }
        Currency out(result);
        if (input1.IsString())
            out.SetMask(Currency::MASK_STRING);
        outputs.push_back(out);
    }
    else
    {
        HML_CELLARRAY *result = EvaluatorInterface::allocateCellArray(m, n);
        int M = cell->M();
        int N = cell->N();
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                (*result)(i, j) = (*cell)(i % M, j % N);
            }
        }
        Currency out(result);
        if (input1.IsString())
            out.SetMask(Currency::MASK_STRING);
        outputs.push_back(out);
    }
    return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template <double (*func)(double)>
static double degreesToFunc(double val)
{
    double calc = (*func)(val * PI / 180);
    return iszero(calc) ? 0.0 : calc;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template <hwComplex (*func)(const hwComplex&)>
static hwComplex degreesToFuncCplx(const hwComplex &val)
{
    hwComplex calc =  val * (PI / 180);
    calc = (*func)(calc);
    if (iszero(calc.Real()))
        calc.Real() = 0.0;
    if (iszero(calc.Imag()))
        calc.Imag() = 0.0;
    return calc;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template <double (*func)(double)>
static double funcToDegrees(double val)
{
    return (*func)(val) * 180 / PI;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template <hwComplex (*func)(const hwComplex&)>
static hwComplex funcToDegreesCplx(const hwComplex &val)
{
    return (*func)(val) * 180 / PI;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template <hwComplex (*func)(double)>
static hwComplex funcToDegreesCplx_c(double val)
{
    return (*func)(val) * 180 / PI;
}
//------------------------------------------------------------------------------
// Converts degrees to radian values [deg2rad]
//------------------------------------------------------------------------------
bool oml_deg2rad(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsScalar())
    {
        outputs.push_back(input.Scalar() * (PI/180.0));
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);
        double scale = PI / 180.0;

        if (mtx->IsReal())
        {
            for (int k = 0; k < mtx->Size(); k++)
                (*result)(k) = (*mtx)(k) * scale;
        }
        else if (mtx->IsRealData())
        {
            for (int k = 0; k < mtx->Size(); k++)
                (*result)(k) = mtx->z(k).Real() * scale;
        }
        else
        {
            delete result;
            throw OML_Error(OML_ERR_SCALARMATRIX, -1, OML_VAR_TYPE);
        }

        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_deg2rad);
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, -1, OML_VAR_TYPE);
    }

    return true;
}
//------------------------------------------------------------------------------
// Converts radian values to degrees [rad2deg]
//------------------------------------------------------------------------------
bool oml_rad2deg(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsScalar())
    {
        outputs.push_back(input.Scalar() * (180.0/PI));
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);
        double scale = 180.0 / PI;

        if (mtx->IsReal())
        {
            for (int k = 0; k < mtx->Size(); k++)
                (*result)(k) = (*mtx)(k) * scale;
        }
        else if (mtx->IsRealData())
        {
            for (int k = 0; k < mtx->Size(); k++)
                (*result)(k) = mtx->z(k).Real() * scale;
        }
        else
        {
            delete result;
            throw OML_Error(OML_ERR_SCALARMATRIX, -1, OML_VAR_TYPE);
        }

        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_rad2deg);
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, -1, OML_VAR_TYPE);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the inverse cosine (in degrees) for each element of x [acosd]
//------------------------------------------------------------------------------
bool oml_acosd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!conditionFunc(inputs, outputs, funcToDegrees<acos>, funcToDegreesCplx<hwComplex::acos>, funcToDegreesCplx_c<hwComplex::acos_c>, absAtMostOne))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_acosd);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the inverse sine in degrees of the input [asind]
//------------------------------------------------------------------------------
bool oml_asind(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!conditionFunc(inputs, outputs, funcToDegrees<asin>, funcToDegreesCplx<hwComplex::asin>, funcToDegreesCplx_c<hwComplex::asin_c>, absAtMostOne))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_asind);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the cosine of the input in degrees [cosd]
//------------------------------------------------------------------------------
bool oml_cosd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, degreesToFunc<cos>, degreesToFuncCplx<hwComplex::cos>))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_cosd);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the sine of the input in degrees [sind]
//------------------------------------------------------------------------------
bool oml_sind(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, degreesToFunc<sin>, degreesToFuncCplx<hwComplex::sin>))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_sind);
    }

    return true;
}
//------------------------------------------------------------------------------
// Vector cross product [cross]
//------------------------------------------------------------------------------
bool oml_cross(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];
    const Currency& input2 = inputs[1];
    int dim = 0;

    if (nargin > 2)
    {
        dim = (int) inputs[2].Scalar();

        if (!IsInteger(inputs[2].Scalar()).IsOk())
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);
        if (dim <= 0)
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);
    }

    if ((input1.IsMatrix() || input1.IsString()) &&  // leave oddball string case here for now
        (input2.IsMatrix() || input2.IsString()))
    {
        if (dim > 2)
            throw OML_Error(OML_ERR_VECLENDIM);

        const hwMatrix* x = input1.Matrix();
        const hwMatrix* y = input2.Matrix();

        if (x->IsVector() && y->IsVector())
        {
            if (x->Size() == 3 && y->Size() == 3)
            {
                if (x->M() != y->M())
                {
                    if (dim)
                        throw OML_Error(HW_ERROR_DIM3ELEM);

                    BuiltInFuncsUtils::SetWarning(eval, "Warning: performing cross product on vectors with different dimensions");
                }
                else if (dim)
                {
                    if ((dim == 1 ? x->M() : x->N()) != 3)
                        throw OML_Error(HW_ERROR_DIM3ELEM);
                }

                hwMatrix *result = EvaluatorInterface::allocateMatrix();
                Currency out(result);
                BuiltInFuncsUtils::CheckMathStatus(eval, result->Cross(*x, *y)); // always returns column vector
                if (x->M() == 1 && y->M() == 1)
                    BuiltInFuncsUtils::CheckMathStatus(eval, result->Transpose());
                outputs.push_back(out);
                return true;
            }
            else
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);
        }

        if (!sameSize(x, y))
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

        if (dim)
        {
            if (3 != (dim == 1 ? x->M() : x->N()))
                throw OML_Error(HW_ERROR_DIM3ELEM);
        }
        else
        {
            if (x->M() == 3)
                dim = 1;
            else if (x->N() == 3)
                dim = 2;
            else
                throw OML_Error(HW_ERROR_INP1DINPMAT3ELE);
        }

        hwMatrix *result = EvaluatorInterface::allocateMatrix(x->M(), x->N(), x->Type());

        if (dim == 1)
        {
            hwMatrix* rcol = EvaluatorInterface::allocateMatrix();
            for (int i = 0; i < x->N(); i++)
            {
                std::unique_ptr<const hwMatrix> xcol(EvaluatorInterface::allocateColumn(x, i));
                std::unique_ptr<const hwMatrix> ycol(EvaluatorInterface::allocateColumn(y, i));
                BuiltInFuncsUtils::CheckMathStatus(eval, rcol->Cross(*xcol, *ycol));
                BuiltInFuncsUtils::CheckMathStatus(eval, result->WriteColumn(i, *rcol));
            }
        }
        else // dim == 2
        {
            std::unique_ptr<hwMatrix> xrow(EvaluatorInterface::allocateMatrix());
            std::unique_ptr<hwMatrix> yrow(EvaluatorInterface::allocateMatrix());
            std::unique_ptr<hwMatrix> rrow(EvaluatorInterface::allocateMatrix());
            for (int i = 0; i < x->M(); i++)
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, x->ReadRow(i, *xrow));
                BuiltInFuncsUtils::CheckMathStatus(eval, y->ReadRow(i, *yrow));
                BuiltInFuncsUtils::CheckMathStatus(eval, rrow->Cross(*xrow, *yrow));
                BuiltInFuncsUtils::CheckMathStatus(eval, rrow->Transpose());
                BuiltInFuncsUtils::CheckMathStatus(eval, result->WriteRow(i, *rrow));
            }
        }

        outputs.push_back(result);
    }
    else if (input1.IsNDMatrix() && input2.IsNDMatrix())
    {
        return oml_MatrixN_VecProd(eval, inputs, outputs, oml_cross, 3);
    }
    else if (!input1.IsMatrix() && !input1.IsNDMatrix())
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
    }
    else if (!input2.IsMatrix() && !input2.IsNDMatrix())
    {
        throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);
    }
    else
    {
        throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is a struct [isstruct]
//------------------------------------------------------------------------------
bool oml_isstruct(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(inputs[0].IsStruct() ? getTrue() : getFalse());
    return true;
}
//------------------------------------------------------------------------------
// Returns true if successful in reading a line from a file [fgetl]
//------------------------------------------------------------------------------
bool oml_fgetl(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    std::string result;
    bool success = dofgets(eval, inputs, result);
    if (success)
    {
        if (!result.empty()) // Strip trailing \r\n
            result.erase(result.find_last_not_of("\r\n")+1);
       
        double len = result.empty() ? 0.0 : static_cast<double>(result.size());

        outputs.push_back(result);
        outputs.push_back(len);
    }
    else
    {
        outputs.push_back(-1.0);
        outputs.push_back(0.0);
    }
    return true;
}
//------------------------------------------------------------------------------
// Converts cell to matrix [cell2mat]
//------------------------------------------------------------------------------
bool oml_cell2mat(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];
    if (!input.IsCellArray())
        throw OML_Error(HW_ERROR_INPUTCELLARRAY);

    enum
    {
        Matrices,
        Cells,
        //Structs,
        Strings
    } inputType;

    HML_CELLARRAY *cell = input.CellArray();

    // get input type
    if (cell->Size())
    {
        const Currency &elem = (*cell)(0, 0);
        if (elem.IsScalar() || elem.IsComplex() || elem.IsMatrix())
        {
            inputType = Matrices;
        }
        else if (elem.IsString())
        {
            inputType = Strings;
        }
        else if (elem.IsCellArray())
        {
            inputType = Cells;
        }
        else if (elem.IsStruct())
        {
            // create struct arrays
            throw OML_Error(HW_ERROR_INPNOTSTRUCT);
        }
        else
            throw OML_Error(HW_MATH_MSG_INVALIDINPUT);
    }
    else
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        return true;
    }

    // dimensions of output
    int m = 0;
    int n = 0;
    int tempn = 0;
    int tempm = 0;
    int *mlist = NULL;
    bool setm = false;
    int i = 0;
    int outrow = -1;

    // verify dimensions
    for (i = 0; i < cell->M(); i++)
    {
        if (setm)
        {
            do
            {
                outrow++;
                tempn = 0;
                increasetemp(mlist, tempn, outrow, n);
            }
            while (tempn == n);
        }
        for (int j = 0; j < cell->N(); j++)
        {
            const Currency &elem = (*cell)(i, j);

            // assume input is valid -- this will be validated later
            if (elem.IsScalar() || elem.IsComplex() || elem.IsCellArray() || elem.IsStruct())
            {
                if (setm)
                {
                    if (tempn >= n)
                    {
                        tempn++;
                        break;
                    }
                    mlist[tempn]++;
                    increasetemp(mlist, ++tempn, outrow, n);
                }
                else
                    n++;
            }
            else if (elem.IsMatrix() || elem.IsString())
            {
                const hwMatrix *mtx = elem.Matrix();

				if (mtx)
					tempm = mtx->M();
				else
					tempm = 0;
                if (mtx && (mtx->M() > 0 || mtx->N() > 0))
                {
                    if (setm)
                    {
                        for (int k = 0; k < mtx->N(); k++)
                        {
                            if (tempn >= n)
                            {
                                tempn++;
                                break;
                            }
                            mlist[tempn] += mtx->M();
                            increasetemp(mlist, ++tempn, outrow, n);
                        }
                    }
                    else
                        n += mtx->N();
                }
				else
				{
					if (mtx && (n == 0))
						tempn = mtx->N(); // get the right size for the empty output
				}
            }
            else
            {
                if (setm)
                    delete [] mlist;

                throw OML_Error(HW_MATH_MSG_INVALIDINPUT);
            }
        }

        if (setm)
        {
            if (tempn != n)
            {
                delete [] mlist;
                throw OML_Error(HW_ERROR_UNEVENDIMENSIONS);
            }
        }
        else
        {
            if (n)
            {
                i = outrow = -1;
                setm = true;
                mlist = new int[n];
                for (int j = 0; j < n; j++)
                    mlist[j] = 0;
            }
            else
            {
                outputs.push_back(EvaluatorInterface::allocateMatrix(tempm, tempn, hwMatrix::REAL));
                return true;
            }
        }
    }

    m = mlist[0];
    mlist[0] = 0;

    // verify number of rows
    for (i = 1; i < n; i++)
    {
        if (mlist[i] != m)
            break;
        mlist[i] = 0; // reset for later use
    }

    if (i != n) // did not finish loop
    {
        delete [] mlist;
        throw OML_Error(HW_ERROR_UNEVENDIMENSIONS);
    }

    if (inputType == Cells)
    {
        delete [] mlist;

        if (cell->M() != m || cell->N() != n)
            throw OML_Error(HW_ERROR_UNEVENDIMENSIONS);

        HML_CELLARRAY *outcell = EvaluatorInterface::allocateCellArray();

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                const Currency &elem = (*cell)(i, j);
                if (elem.IsCellArray())
                {
                    HML_CELLARRAY *elemCell = elem.CellArray();
                    if (elemCell->Size())
                    {
                        HML_CELLARRAY *temp = outcell;
                        outcell = hconcat(temp, elemCell);
                        delete temp;
                    }
                }
                else
                {
                    delete outcell;
                    throw OML_Error(HW_ERROR_MIXEDCELLELEMS);
                }
            }
        }

        outputs.push_back(outcell);
        return true;
    }

    hwMatrix *outMatrix = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
    Currency outCur(outMatrix);
    if (inputType == Strings)
        outCur.SetMask(Currency::MASK_STRING);

    outrow = -1;
    for (i = 0; i < cell->M(); i++)
    {
        do
        {
            outrow++;
            tempn = 0;
            increasetemp(mlist, tempn, outrow, n);
        }
        while (tempn == n);

        for (int j = 0; j < cell->N(); j++)
        {
            const Currency &elem = (*cell)(i, j);
            switch (inputType)
            {
            case Matrices:
                if (elem.IsScalar())
                {
                    outMatrix->SetElement(mlist[tempn]++, tempn, elem.Scalar());
                    increasetemp(mlist, ++tempn, outrow, n);
                }
                else if (elem.IsComplex())
                {
                    outMatrix->SetElement(mlist[tempn]++, tempn, elem.Complex());
                    increasetemp(mlist, ++tempn, outrow, n);
                }
                else if (elem.IsMatrix())
                {
                    const hwMatrix *mtx = elem.Matrix();
                    if (mtx->IsReal())
                    {
                        for (int l = 0; l < mtx->N(); l++)
                        {
                            for (int k = 0; k < mtx->M(); k++)
                            {
                                outMatrix->SetElement(mlist[tempn]++, tempn, (*mtx)(k, l));
                            }
                            increasetemp(mlist, ++tempn, outrow, n);
                        }
                    }
                    else
                    {
                        for (int l = 0; l < mtx->N(); l++)
                        {
                            for (int k = 0; k < mtx->M(); k++)
                            {
                                outMatrix->SetElement(mlist[tempn]++, tempn, mtx->z(k, l));
                            }
                            increasetemp(mlist, ++tempn, outrow, n);
                        }
                    }
                }
                else
                {
                    delete [] mlist;
                    throw OML_Error(HW_ERROR_MIXEDCELLELEMS);
                }
                break;
            case Strings:
                if (elem.IsScalar())
                {
                    outMatrix->SetElement(mlist[tempn]++, tempn,
                        BuiltInFuncsUtils::GetValidChar(eval, elem.Scalar(), false));
                    increasetemp(mlist, ++tempn, outrow, n);
                }
                else if (elem.IsComplex())
                {
                    throw OML_Error(HW_ERROR_MIXEDCELLELEMS);
                }
                else if (elem.IsMatrix() || elem.IsString())
                {
                    const hwMatrix *mtx = elem.Matrix();
                    if (mtx->IsRealData())
                    {
                        for (int l = 0; l < mtx->N(); l++)
                        {
                            for (int k = 0; k < mtx->M(); k++)
                            {
                                outMatrix->SetElement(mlist[tempn]++, tempn,
                                    BuiltInFuncsUtils::GetValidChar(eval, realval(mtx, k, l), false));
                            }
                            increasetemp(mlist, ++tempn, outrow, n);
                        }
                    }
                    else
                    {
                        delete [] mlist;
                        throw OML_Error(HW_ERROR_MIXEDCELLELEMS);
                    }
                }
                else
                {
                    delete [] mlist;
                    throw OML_Error(HW_ERROR_MIXEDCELLELEMS);
                }
                break;
            default:
                delete [] mlist;
                throw OML_Error(HW_MATH_MSG_INTERNALERROR);
            }
        }
    }

    outputs.push_back(outCur);
    return true;
}
//------------------------------------------------------------------------------
// Returns largest floating point that is representable [realmax]
//------------------------------------------------------------------------------
bool oml_realmax(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return limitFunc(eval, inputs, outputs, std::numeric_limits<float>::max, std::numeric_limits<double>::max);
}
//------------------------------------------------------------------------------
// Returns largest floating point that is representable [realmin]
//------------------------------------------------------------------------------
bool oml_realmin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return limitFunc(eval, inputs, outputs, std::numeric_limits<float>::min, std::numeric_limits<double>::min);
}
//------------------------------------------------------------------------------
// Returns true if file exists in path
//------------------------------------------------------------------------------
static bool FileExistsInPath(const EvaluatorInterface& eval, std::string& filename)
{
    
#if OS_WIN
    size_t lastslash = filename.find_last_of("/\\");
#else
    size_t lastslash = filename.find_last_of('/');
#endif
    if (lastslash == std::string::npos)
        lastslash = 0;

    size_t lastdot = filename.find_last_of('.');

    bool hasext = lastdot != std::string::npos && lastdot > lastslash;
    if (!hasext)
    {
        if (eval.FindFileInPath(filename + ".oml", filename) || eval.FindFileInPath(filename + ".hml", filename) || eval.FindFileInPath(filename + ".m", filename))
            return true;
    }
    return eval.FindFileInPath(filename, filename);
}
//------------------------------------------------------------------------------
// Executes given file/command [run]
//------------------------------------------------------------------------------
bool oml_run(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency curstr = toCurrencyStr(eval, inputs[0], false, false);
    std::string filename = readString(curstr);
            
    if (!fileExists(filename))
    {
    if (!FileExistsInPath(eval, filename))
        throw OML_Error(HW_ERROR_NOSUCHFILE(filename));
    }

    std::string full_file_name = BuiltInFuncsUtils::GetAbsolutePath(filename);
    std::string cwd = BuiltInFuncsUtils::GetCurrentWorkingDir();

#ifdef OS_WIN
    size_t slash_index = full_file_name.find_last_of("/\\");
#else
    size_t slash_index = full_file_name.rfind('/');
#endif

    // extract path so we can cd to it
    if (slash_index != std::string::npos)
    {
        std::string pathname(full_file_name.begin(), full_file_name.begin() + slash_index);
        cd(pathname, eval);
    }

    Interpreter interp(eval);
    Currency ret;
    std::string errmsg;
    bool        formatmsg = true;

	size_t dot_index = full_file_name.rfind('.');

	if (dot_index != std::string::npos)
	{
		std::string extension(full_file_name.begin()+dot_index+1, full_file_name.end());

		if (eval.IsExtensionEncrypted(extension))
		{
			eval.RunEncryptedFile(full_file_name, extension);
			return true;
		}
		else if (extension == "omlp")
		{
			eval.RunPrecompiledFile(full_file_name);
			return true;
		}
	}

    try
    {
        ret = interp.DoFile(full_file_name);
        if (ret.IsError())
            errmsg = ret.Message();
    }
    catch (OML_Error& e)
    {
        errmsg = e.GetErrorMessage();
        formatmsg = e.GetFormatMessage();
    }

    std::vector<Currency> results = interp.GetOutputCurrencyList();

	if (slash_index != std::string::npos)
		cd(cwd, eval);

     // to output errors from interp.DoFile()
    if (errmsg.length())
    {
        throw OML_Error(errmsg, formatmsg);
    }
    else
    {
        for (std::vector<Currency>::iterator it = results.begin(); it != results.end(); it++)
        {
            if (it->IsError())
                throw OML_Error(it->Message());

            eval.PushResult(*it);
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Runs the given file [run2]
//------------------------------------------------------------------------------
bool oml_run2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

	std::string filename = inputs[0].StringVal();

	eval.RunFile(filename);
	return true;
}
//------------------------------------------------------------------------------
// Flushes any content written to the file stream [fflush]
//------------------------------------------------------------------------------
bool oml_fflush(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    int fileID = getFileFromInput(eval, inputs[0]);
    checkFileIndex(eval, fileID, true);

    std::string mode = eval.GetFileMode(fileID);
    if (mode == "rb" || mode == "r")
    {
        outputs.push_back(-1.0);
    }
    else
    {
        std::FILE *f = eval.GetFile(fileID);
        int returnval = fflush(f);
        outputs.push_back(returnval);
    }
    return true;
}
//------------------------------------------------------------------------------
// Rewinds to the beginning of the given file [frewind]
//------------------------------------------------------------------------------
bool oml_frewind(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    int fileID = getFileFromInput(eval, inputs[0]);
    checkFileIndex(eval, fileID, false);

    std::FILE *f = eval.GetFile(fileID);
    rewind(f);

    if (getNumOutputs(eval))
        outputs.push_back(ferror(f) ? -1.0 : 0.0);

    return true;
}
//------------------------------------------------------------------------------
// Gets precision from string
//------------------------------------------------------------------------------
static Precision getPrecision(std::string str)
{
    // no support yet for conversion types
    int blockSize = 1;
    size_t index = str.find('*');
    if (index != std::string::npos)
    {
        std::string left = str.substr(0, index);
        const char* cleft = left.c_str();
        char* end;
        double temp = strtod(cleft, &end);
        if (end && end > cleft)
        {
            if (!isint(temp) || temp < 1.0)
                throw OML_Error(HW_ERROR_PRECBLOCKPOSINT);

            blockSize = (int) temp;
            str.erase(0, index + 1);
        }
    }
    str = convertToLower(str);

    if (str == "float")
        return Precision(true, sizeof(float), blockSize, Float);
    else if (str == "ulong" || str == "unsigned long")
        return Precision(false, sizeof(unsigned long), blockSize, Long);
    else if (str == "uint" || str == "unsigned int")
        return Precision(false, sizeof(unsigned int), blockSize, Int);
    else if (str == "ushort" || str == "unsigned short")
        return Precision(false, sizeof(unsigned short), blockSize, Short);
    else if (str == "long")
        return Precision(true, sizeof(long), blockSize, Long);
    else if (str == "int")
        return Precision(true, sizeof(int), blockSize, Int);
    else if (str == "short")
        return Precision(true, sizeof(short), blockSize, Short);
    else if (str == "char" || str == "char*1" || str == "schar" || str == "signed char")
        return Precision(true, sizeof(char), blockSize, Char);
    else if (str == "double" || str == "float64" || str == "real*8")
        return Precision(true, sizeof(double), blockSize, Double);
    else if (str == "single" || str == "float32" || str == "real*4")
        return Precision(true, sizeof(float), blockSize, Float);
    else if (str == "uint64")
        return Precision(false, sizeof(uint64_t), blockSize, Int64);
    else if (str == "uint32")
        return Precision(false, sizeof(uint32_t), blockSize, Int32);
    else if (str == "uint16")
        return Precision(false, sizeof(uint16_t), blockSize, Int16);
    else if (str == "uint8")
        return Precision(false, sizeof(uint8_t), blockSize, Int8);
    else if (str == "int8" || str == "integer*1")
        return Precision(true, sizeof(int8_t), blockSize, Int8);
    else if (str == "int16" || str == "integer*2")
        return Precision(true, sizeof(int16_t), blockSize, Int16);
    else if (str == "int32" || str == "integer*4")
        return Precision(true, sizeof(int32_t), blockSize, Int32);
    else if (str == "int64" || str == "integer*8")
        return Precision(true, sizeof(int64_t), blockSize, Int64);
    else if (str == "uchar" || str == "unsigned char")
        return Precision(false, sizeof(unsigned char), blockSize, Char);
    else
        throw OML_Error(HW_ERROR_INVPRECTYPE);
}
//------------------------------------------------------------------------------
// Gets precision from string
//------------------------------------------------------------------------------
static Precision getPrecision(EvaluatorInterface& eval, Currency input)
{
    return getPrecision(readOption(eval, unnest(input, HW_ERROR_INVPRECTYPE)));
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template <typename T>
static void checkBounds(T *writeTo, double value)
{
    if (std::numeric_limits<T>::is_integer)
    {
        value = round(value);
        T max = (std::numeric_limits<T>::max)();
        T min = (std::numeric_limits<T>::min)();
        if (value < (double) min)
            *writeTo = min;
        else if (value > (double) max)
            *writeTo = max;
        else
            *writeTo = (T) value;
    }
    else
    {
        *writeTo = (T) value;
    }
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template<typename T>
static int doWrite(EvaluatorInterface &eval, const std::vector<Currency> &inputs, int blockSize, size_t dataSize)
{
    T *towrite;
    size_t nargin = inputs.size();
    int skip = 0;
    int numItems;
    int fileID = getFileFromInput(eval, inputs[0]);
    checkFileIndex(eval, fileID, true);
    std::FILE* file = eval.GetFile(fileID);

    if (file == stdin)
        return -1;

    if (nargin > 3)
    {
        skip = (int) inputs[3].Scalar();

        if (!IsInteger(inputs[3].Scalar()).IsOk())
            throw OML_Error(OML_ERR_NATURALNUM, 4, OML_VAR_SKIPVAL);
        if (skip < 0)
            throw OML_Error(OML_ERR_NATURALNUM, 4, OML_VAR_SKIPVAL);
        if (skip && fileID < FIRST_USER_FILE)
            throw OML_Error(HW_ERROR_NOTSKIPBUILTINFUNC);
    }
    
    FileStreamType ftype = determineFileStreamType(file);

    const Currency &input2 = inputs[1];
    if (input2.IsScalar())
    {
        numItems = 1;
        towrite = new T[numItems];
        checkBounds(towrite, input2.Scalar());
    }
    else if (input2.IsComplex())
    {
        numItems = 1;
        towrite = new T[numItems];
        checkBounds(towrite, input2.Complex().Real());
    }
    else if (input2.IsMatrix())
    {
        const hwMatrix *m = input2.Matrix();
        numItems = m->Size();
        if (numItems)
        {
            towrite = new T[numItems];
            if (m->IsReal())
            {
                for (int i = 0; i < numItems; i++)
                {
                    checkBounds(towrite + i, (*m)(i));
                }
            }
            else
            {
                for (int i = 0; i < m->Size(); i++)
                {
                    checkBounds(towrite + i, m->z(i).Real());
                }
            }
        }
        else
        {
            return 0;
        }
    }
    else if (input2.IsString())
    {
        std::string str = orderedStringVal(input2);
        numItems = (int)str.length();
        towrite = new T[numItems];
        for (int i = 0; i < numItems; i++)
        {
            checkBounds(towrite + i, str[i]);
        }
    }
    else
    {
        if (ftype != Other)
            fclose(file);

        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING);
    }

    int numWritten = 0;

    if (skip)
    {
        int numLoops = (int)ceil(numItems / (double) blockSize);
        for (int i = 0; i < numLoops; i++)
        {
            fseek(file, skip, SEEK_CUR);
            // don't write extra
            numWritten += (int)fwrite(towrite + i * blockSize, dataSize,
                i == numLoops - 1 ? ((numItems - 1) % blockSize) + 1 : blockSize, file);
        }
    }
    else
    {
        numWritten = (int)fwrite(towrite, dataSize, numItems, file);
    }

    delete [] towrite;
    if (ferror(file))
    {
        if (ftype != Other)
            fclose(file);
        return -1;
    }

    if (ftype != Other)
    {
        // readTmpFile closes the stream
        std::string result = readTmpFile(file);
        if (ftype == Stdout)
            eval.PrintResult(result);
        else if (ftype == Stderr)
            throw OML_Error(result);
    }
    return numWritten;
}
//------------------------------------------------------------------------------
// Writes to a file [fwrite]
//------------------------------------------------------------------------------
bool oml_fwrite(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    int blockSize = 1;
    bool signedInput = false;
    size_t dataSize = 1;
    DataType dtype = Char;

    // get precision data
    if (nargin > 2)
    {
        const Currency &input3 = inputs[2];
        Precision p = getPrecision(eval, input3);
        signedInput = p.sign;
        dtype = p.dtype;
        blockSize = p.blockSize;
        dataSize = p.numBytes;
    }

    int numWritten;

    if (signedInput)
    {
        switch (dtype)
        {
        case Double:
            numWritten = doWrite<double>(eval, inputs, blockSize, dataSize);
            break;
        case Int:
            numWritten = doWrite<signed int>(eval, inputs, blockSize, dataSize);
            break;
        case Short:
            numWritten = doWrite<signed short>(eval, inputs, blockSize, dataSize);
            break;
        case Long:
            numWritten = doWrite<signed long>(eval, inputs, blockSize, dataSize);
            break;
        case Char:
            numWritten = doWrite<signed char>(eval, inputs, blockSize, dataSize);
            break;
        case Float:
            numWritten = doWrite<float>(eval, inputs, blockSize, dataSize);
            break;
        case LongLong:
            numWritten = doWrite<signed long long>(eval, inputs, blockSize, dataSize);
            break;
        case Int8:
            numWritten = doWrite<int8_t>(eval, inputs, blockSize, dataSize);
            break;
        case Int16:
            numWritten = doWrite<int16_t>(eval, inputs, blockSize, dataSize);
            break;
        case Int32:
            numWritten = doWrite<int32_t>(eval, inputs, blockSize, dataSize);
            break;
        case Int64:
            numWritten = doWrite<int64_t>(eval, inputs, blockSize, dataSize);
            break;
        default:
            throw OML_Error(HW_MATH_MSG_INTERNALERROR);
        }
    }
    else
    {
        switch (dtype)
        {
        case Int:
            numWritten = doWrite<unsigned int>(eval, inputs, blockSize, dataSize);
            break;
        case Short:
            numWritten = doWrite<unsigned short>(eval, inputs, blockSize, dataSize);
            break;
        case Long:
            numWritten = doWrite<unsigned long>(eval, inputs, blockSize, dataSize);
            break;
        case Char:
            numWritten = doWrite<unsigned char>(eval, inputs, blockSize, dataSize);
            break;
        case LongLong:
            numWritten = doWrite<unsigned long long>(eval, inputs, blockSize, dataSize);
            break;
        case Int8:
            numWritten = doWrite<uint8_t>(eval, inputs, blockSize, dataSize);
            break;
        case Int16:
            numWritten = doWrite<uint16_t>(eval, inputs, blockSize, dataSize);
            break;
        case Int32:
            numWritten = doWrite<uint32_t>(eval, inputs, blockSize, dataSize);
            break;
        case Int64:
            numWritten = doWrite<uint64_t>(eval, inputs, blockSize, dataSize);
            break;
        default:
            throw OML_Error(HW_MATH_MSG_INTERNALERROR);
        }
    }
    outputs.push_back(numWritten);
    return true;
}
//------------------------------------------------------------------------------
// Returns infinite value(s) [inf]
//------------------------------------------------------------------------------
bool oml_inf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return createCommonMatrix(eval, inputs, outputs, Currency(std::numeric_limits<double>::infinity()));
}
//------------------------------------------------------------------------------
// Returns values set to NaN [nan]
//------------------------------------------------------------------------------
bool oml_nan(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return createCommonMatrix(eval, inputs, outputs, Currency(std::numeric_limits<double>::quiet_NaN()));
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template<typename T>
static void doRead(int maxLoops, 
                   size_t size, 
                   int blockSize, 
                   int skip, 
                   std::FILE* file, 
                   int nrows, 
                   int ncols, 
                   std::vector<Currency>& outputs)
{
    if (!(nrows && ncols && blockSize))
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        outputs.push_back(0);
        return;
    }

    int count = 0;

    std::vector<T> nums;
    for (int i = 0; i != maxLoops; i++)
    {
        if (feof(file)) break;
        T *output = new T[blockSize];
        assert(output);
        if (!output) break;

        memset(output, 0, blockSize * size);
        size_t numread = fread(output, size, blockSize, file);
        if (numread <= 0 || numread > blockSize)
        {
            delete [] output;
            output = NULL;
            break;
        }

        count += static_cast<int>(numread);
        // Case when numread is less than blocksize
        for (int j = 0; j < numread; ++j) 
            nums.push_back(output[j]);

        if (skip)
            fseek(file, skip, SEEK_CUR);

        delete [] output;
        output = NULL;
    }

    if (ferror(file))
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        outputs.push_back(-1);
        return;
    }
    
    // don't create unnecessary columns -- the size won't necessarily be exactly what the user defined in this case
    if (nrows > -1)
      ncols = (int)min((double) ncols, ceil((double)(nums.size() / (double) nrows)));

    outputs.push_back(containerToMatrix(nums, nrows, ncols));
    outputs.push_back(static_cast<int>(count));
}
//------------------------------------------------------------------------------
// Reads from file [fread]
//------------------------------------------------------------------------------
bool oml_fread(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    int fileID = getFileFromInput(eval, inputs[0]);
    int skip = 0;
    int maxLoops = -1; // unlimited
    int nrows = -1;
    int ncols = -1;
    int blockSize = 1;
    bool signedOutput = false;
    DataType dtype = Char;
    size_t size = sizeof(unsigned char);

    checkFileIndex(eval, fileID, true);

    if (nargin > 1)
    {
        const Currency &input2 = inputs[1];
        if (nargin == 2 && input2.IsString())
        {
            Precision p = getPrecision(eval, input2);
            signedOutput = p.sign;
            blockSize = p.blockSize;
            size = p.numBytes;
            dtype = p.dtype;
        }
        else
        {
            if (input2.IsScalar())
            {
                double temp = input2.Scalar();
                if (temp == std::numeric_limits<double>::infinity())
                    maxLoops = -1;
                else if (temp < 0)
                    throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);
                else
                    maxLoops = (int) round(temp);
            }
            else if (input2.IsMatrix())
            {
                const hwMatrix *m = input2.Matrix();
                if (m->IsReal())
                {
                    if (m->Size() == 2)
                    {
                        nrows = (int) round((*m)(0));
                        ncols = (int) round((*m)(1));
                        maxLoops = nrows * ncols;

                        if ((*m)(0) == std::numeric_limits<double>::infinity())
                            maxLoops = nrows = -1;
                        else if (nrows < 0)
                            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);

                        if ((*m)(1) == std::numeric_limits<double>::infinity())
                            maxLoops = ncols = -1;
                        else if (ncols < 0)
                            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);
                    }
                    else
                        throw OML_Error(OML_ERR_VECTOR2, 2);
                }
                else
                    throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VECTOR);
            }
            else
                throw OML_Error(HW_ERROR_INVINPTYPESIZEINP);
        }
        if (nargin > 2)
        {
            const Currency &input3 = inputs[2];
            Precision p = getPrecision(eval, input3);
            signedOutput = p.sign;
            blockSize = p.blockSize;
            size = p.numBytes;
            dtype = p.dtype;

            if (nargin > 3)
            {
                skip = (int) inputs[3].Scalar();

                if (!IsInteger(inputs[3].Scalar()).IsOk())
                    throw OML_Error(OML_ERR_NATURALNUM, 4, OML_VAR_SKIPVAL);

                if (skip && fileID < FIRST_USER_FILE)
                    throw OML_Error(HW_ERROR_NOTSKIPBUILTINFUNC);
            }
        }
    }

    std::FILE* file = eval.GetFile(fileID);

    if (maxLoops != -1)
    {
        assert(blockSize != 0);
        if (blockSize == 0)
            blockSize = 1;
        maxLoops = (int) ceil(maxLoops / (double) blockSize);
    }

    if (signedOutput)
    {
        switch (dtype)
        {
        case Double: 
            doRead<double>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Int:
            doRead<signed int>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Short:
            doRead<signed short>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Long:
            doRead<signed long>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Char:
            doRead<signed char>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Float:
            doRead<float>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case LongLong:
            doRead<signed long long>(maxLoops, size, blockSize, skip, file, nrows,
                ncols, outputs);
            break;
        case Int8:
            doRead<int8_t>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Int16:
            doRead<int16_t>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Int32:
            doRead<int32_t>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Int64:
            doRead<int64_t>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        default: throw OML_Error(HW_MATH_MSG_INTERNALERROR); break;
        }
    }
    else
    {
        switch (dtype)
        {
        case Int:
            doRead<unsigned int>(maxLoops, size, blockSize, skip, file, nrows,
                ncols, outputs);
            break;
        case Short:
            doRead<unsigned short>(maxLoops, size, blockSize, skip, file, nrows,
                ncols, outputs);
            break;
        case Long:
            doRead<unsigned long>(maxLoops, size, blockSize, skip, file, nrows,
                ncols, outputs);
            break;
        case Char:
            doRead<unsigned char>(maxLoops, size, blockSize, skip, file, nrows,
                ncols, outputs);
            break;
        case LongLong:
            doRead<unsigned long long>(maxLoops, size, blockSize, skip, file, nrows,
                ncols, outputs);
            break;
        case Int8:
            doRead<uint8_t>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Int16:
            doRead<uint16_t>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Int32:
            doRead<uint32_t>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        case Int64:
            doRead<uint64_t>(maxLoops, size, blockSize, skip, file, nrows, ncols, outputs);
            break;
        default:
            throw OML_Error(HW_MATH_MSG_INTERNALERROR); break;
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns file ID of standard input stream [stdin]
//------------------------------------------------------------------------------
bool oml_stdin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return keywordFunc(inputs, outputs, Currency(0.0));
}
//------------------------------------------------------------------------------
// Returns file ID of standard output stream [stdout]
//------------------------------------------------------------------------------
bool oml_stdout(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return keywordFunc(inputs, outputs, Currency(1.0));
}
//------------------------------------------------------------------------------
// Returns file ID of stderr [stderr]
//------------------------------------------------------------------------------
bool oml_stderr(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return keywordFunc(inputs, outputs, Currency(2.0));
}
//------------------------------------------------------------------------------
// Used in fseek to rewind file position to beginning of file [SEEK_SET]
//------------------------------------------------------------------------------
bool oml_seek_set(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return keywordFunc(inputs, outputs, Currency((double) SEEK_SET + SEEK_OFFSET));
}
//------------------------------------------------------------------------------
// Used in fseek to change file position to current position [SEEK_CUR]
//------------------------------------------------------------------------------
bool oml_seek_cur(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return keywordFunc(inputs, outputs, Currency((double) SEEK_CUR + SEEK_OFFSET));
}
//------------------------------------------------------------------------------
// Used in fseek to change file position to end of file [SEEK_END]
//------------------------------------------------------------------------------
bool oml_seek_end(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return keywordFunc(inputs, outputs, Currency((double) SEEK_END + SEEK_OFFSET));
}
//------------------------------------------------------------------------------
// Sets the file pointer to a position in given fileID [fseek]
//------------------------------------------------------------------------------
bool oml_fseek(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    const Currency &input2 = inputs[1];

    int fileID = getFileFromInput(eval, input1);

    checkFileIndex(eval, fileID, true);
    if (fileID < FIRST_USER_FILE)
        throw OML_Error(HW_ERROR_NOTFSEEKONSTDINOUTERR);

    int offset = (int) input2.Scalar();

    if (!IsInteger(input2.Scalar()).IsOk())
        throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_OFFSET);

    int origin = SEEK_SET;
    if (nargin > 2)
    {
        const Currency &input3 = inputs[2];
        if (input3.IsString())
        {
            std::string orig = readString(input3);
            if (orig == "cof")
                origin = SEEK_CUR;
            else if (orig == "bof")
                origin = SEEK_SET;
            else if (orig == "eof")
                origin = SEEK_END;
            else
            {
                outputs.push_back(-1.0);
                return true;
            }
        }
        else
        {
            origin = (int) input3.Scalar() - SEEK_OFFSET;

            if (!IsInteger(input3.Scalar()).IsOk())
                throw OML_Error(OML_ERR_INTEGER, 4, OML_VAR_ORIGIN);

            if (!(origin == SEEK_CUR || origin == SEEK_SET || origin == SEEK_END))
            {
                outputs.push_back(-1.0);
                return true;
            }
        }
    }

    std::FILE *file = eval.GetFile(fileID);

    // to prevent from extending past end of file
    // there may be a better way of doing this
    if (offset > 0)
    {
        long currentPos = ftell(file);
        if (fseek(file, offset, origin))
            outputs.push_back(-1.0);
        else
        {
            long desiredPos = ftell(file);

            if (fseek(file, 0, SEEK_END))
                outputs.push_back(-1.0);
            else
            {
                long endPos = ftell(file);

                if (desiredPos > endPos)
                {
                    outputs.push_back(-1.0);
                    fseek(file, currentPos, SEEK_SET); // go back to where we started
                }
                else
                {
                    outputs.push_back(fseek(file, desiredPos, SEEK_SET)); // go to desired position
                }
            }
        }
    }
    else
        outputs.push_back(fseek(eval.GetFile(fileID), offset, origin) ? -1.0 : 0.0);

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is logical [islogical]
//------------------------------------------------------------------------------
bool oml_islogical(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(HW_MAKE_BOOL_CURRENCY(inputs[0].IsLogical()));
    return true;
}
//------------------------------------------------------------------------------
// Returns true if input file is at the end of the stream [feof]
//------------------------------------------------------------------------------
bool oml_feof(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    int fileID = getFileFromInput(eval, inputs[0]);
    checkFileIndex(eval, fileID, true);
    FILE* f = eval.GetFile(fileID);
    if (feof(f))
    {
        outputs.push_back(getTrue());
    }
    else
    {
        char c = fgetc(f);
        if (c == EOF)
        {
            outputs.push_back(getTrue());
        }
        else
        {
            ungetc(c,f);
            outputs.push_back(getFalse());
        }
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns current position of input file pointer [ftell]
//------------------------------------------------------------------------------
bool oml_ftell(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    int fileID = getFileFromInput(eval, inputs[0]);
    checkFileIndex(eval, fileID, true);

    if (fileID < FIRST_USER_FILE)
        throw OML_Error(HW_ERROR_NOTFTELLONSTDINOUTERR);

    long tell = ftell(eval.GetFile(fileID));
    if (tell == -1L) // error
    {
        std::string error = std::string("Error: error during tell: ") + strerror(errno);
        throw OML_Error(error);
    }

    outputs.push_back((double) tell);
    return true;
}
//------------------------------------------------------------------------------
// Reads characters in a file [fgets]
//------------------------------------------------------------------------------
bool oml_fgets(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    std::string result;
    if (dofgets(eval, inputs, result))
    {
        outputs.push_back(result);
        outputs.push_back((double) result.length());
    }
    else
    {
        outputs.push_back(-1.0);
        outputs.push_back(0.0);
    }
    return true;
}
//------------------------------------------------------------------------------
// Closes specified input file [fclose]
//------------------------------------------------------------------------------
bool oml_fclose(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsString() && readString(input) == "all")
    {        
        double result = eval.CloseAllFiles() ? 0.0 : 1.0; // Return double
        outputs.push_back(result);
        return true;
    }

    int fileID = getFileFromInput(eval, inputs[0]);
    checkFileIndex(eval, fileID, false);
    outputs.push_back(eval.CloseFile(fileID) ? 0.0 : 1.0);
    
    return true;
}
//------------------------------------------------------------------------------
// Opens specified input file [fopen]
//------------------------------------------------------------------------------
bool oml_fopen(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency input = inputs[0];

    if (nargin > 1 || input.IsString() || input.IsCellArray())
    {
        input = unnest(input, "Error: invalid input for mode; should be a string");
        std::string fname = readString(toCurrencyStr(eval, input, false, false));

        if (nargin == 1 && fname == "all")
        {
            std::vector<int> indices = eval.GetFileIndices(FIRST_USER_FILE);
            hwMatrix *m = containerToMatrix(indices);
            outputs.push_back(m);
        }
        else
        {
            std::string mode("rb");

            if (nargin > 1)
            {
                Currency modeCur = unnest(inputs[1], "Error: invalid input for mode; should be a string");
                mode = orderedStringVal(toCurrencyStr(eval, modeCur, false, false));

                size_t len = mode.length();

                // verify mode input
                do
                {
                    if (!len)
                    {
                        throw OML_Error(HW_ERROR_INVMODESTR);
                    }
                    if (len < 4)
                    {
                        char first = mode[0];
                        if (first == 'r' || first == 'w' || first == 'a')
                        {
                            if (len == 2)
                            {
                                char second = mode[1];
                                if (second == '+')
                                {
                                    mode += 'b';
                                    break;
                                }
                                else if (second == 'b')
                                {
                                    break;
                                }
                                else if (second == 't')
                                {
                                    mode.erase(1);
                                    break;
                                }
                            }
                            else if (len == 3)
                            {
                                char second = mode[1];
                                char third = mode[2];
                                if (second == '+')
                                {
                                    if (third == 'b')
                                        break;
                                    else if (third == 't')
                                    {
                                        mode.erase(2);
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                mode += 'b'; // default mode is binary
                                break;
                            }
                        }
                    }

                    throw OML_Error(HW_ERROR_INVMODESTR);

                }
                while (0);
            }

            std::FILE* f = fopen(fname.c_str(), mode.c_str());

            if (f == nullptr)
            {
                outputs.push_back(-1.0);
                outputs.push_back(std::string(strerror(errno)));
            }
            else
            {
                outputs.push_back(eval.AddFile(f, fname, mode));
                outputs.push_back(std::string());
            }
        }
    }
    else
    {
        int fileId = (int) input.Scalar();

        if (!IsInteger(input.Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_FILEID);

        checkFileIndex(eval, (int) fileId, true);

        outputs.push_back(eval.GetFileName(fileId));
        outputs.push_back(eval.GetFileMode(fileId));
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns a string matrix with input arguments, vertically concatenated [char]
//------------------------------------------------------------------------------
bool oml_char(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    std::vector<Currency> strCurrencies;
    int m = 0, n = 0;

    for (int i = 0; i < inputs.size(); i++)
    {
        Currency temp = toCurrencyStr(eval, inputs[i], false, false);

		if (inputs[i].GetMask() == Currency::MASK_STRING)
			m += max(1, temp.Matrix()->M());
		else if (temp.Matrix()->M() == 0)
			;
		else
			m += max(1, temp.Matrix()->M());

        n = max(n, temp.Matrix()->N());
        strCurrencies.push_back(temp);
    }
	
    hwMatrix *out = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);

	if (!out->IsEmpty())
		buildString(strCurrencies, out);

    outputs.push_back(addStringMask(out));
    return true;
}
//------------------------------------------------------------------------------
// Tokenizes input string with the given delimiter [strsplit]
//------------------------------------------------------------------------------
bool oml_strsplit(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    if (size < 1 || size > 6 || size == 5)
        throw OML_Error(OML_ERR_NUMARGIN);

    bool collapseDelimiters = true;
    bool useRegex = false;
    bool collapseSet = false;
    bool regexSet = false;

    std::string tosplit;
    std::vector<std::string> delims;
    const Currency &input1 = inputs[0];

    if (input1.IsString())
        tosplit = readString(input1);
    else
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

    if (size == 3)
    {
        const Currency &input3 = inputs[2];
        if (input3.IsString())
            throw OML_Error("Error: option '" + orderedStringVal(inputs[2]) + "' needs a value");
        else if (input3.IsCellArray())
            throw OML_Error(HW_MATH_MSG_INVALIDINPUT);
    }

    if (size > 1)
    {
        const Currency &input2 = inputs[1];

        if (input2.IsCellArray())
        {
            if (size == 3)
            {
                const Currency &input3 = inputs[2];
                collapseSet = true;
                collapseDelimiters = (bool) inputs[2].Scalar();

                if (!IsInteger(inputs[2].Scalar()).IsOk())
                    throw OML_Error(OML_ERR_FLAG_01, 3, OML_VAR_VALUE);
            }
            HML_CELLARRAY *cell = input2.CellArray();
            for (int i = 0; i < cell->Size(); i++)
            {
                Currency c = (*cell)(i);
                if (c.IsString())
                {
                    std::string str = readString(c);
                    if (str.length())
                        delims.push_back(str);
                }
                else
                    throw OML_Error(HW_ERROR_DELONLYSTR);
            }
        }
        else if (input2.IsString())
        {
            if (size == 3)
            {
                collapseSet = true;
                collapseDelimiters = (bool) inputs[2].Scalar();

                if (!IsInteger(inputs[2].Scalar()).IsOk())
                    throw OML_Error(OML_ERR_FLAG_01, 3, OML_VAR_VALUE);
            }
            std::string str = readString(input2);
            if (str.length())
                delims.push_back(str);
        }
        else if (size == 2 || size == 3)
        {
            // special cases for octave's backwards compatibility
            collapseSet = true;
            collapseDelimiters = (bool) input2.Scalar();

            if (!IsInteger(input2.Scalar()).IsOk())
                throw OML_Error(OML_ERR_FLAG_01, 2, OML_VAR_VALUE);

            delims.push_back(" "); // default delimiter
        }
        else
            throw OML_Error(HW_ERROR_DELONLYSTR);
    }
    else
        delims.push_back(" ");

    for (int i = 2; i < size - 1; i += 2)
    {
        const Currency &option = inputs[i];
        const Currency &value = inputs[i + 1];
        if (option.IsString())
        {
            if (option.Matrix()->M() > 1)
                throw OML_Error(HW_ERROR_INPSTRONLY1D);

            std::string str = convertToLower(orderedStringVal(option));
            if (str == "collapsedelimiters")
            {
                if (collapseSet)
                    throw OML_Error(HW_ERROR_COLAPSDELALLSET);
                collapseSet = true;
                if (value.IsScalar())
                    collapseDelimiters = (bool) value.Scalar();
                else
                    throw OML_Error(HW_ERROR_COLAPSDELALLSET);
            }
            else if (str == "delimitertype")
            {
                if (regexSet)
                    throw OML_Error(HW_ERROR_DELTYPEALLSET);
                regexSet = true;
                if (value.IsString())
                {
                    if (value.Matrix()->M() > 1)
                        throw OML_Error(HW_ERROR_INPSTRONLY1D);

                    std::string str = convertToLower(orderedStringVal(value));
                    if (str == "regularexpression")
                        useRegex = true;
                    else if (str != "simple")
                        throw OML_Error(HW_ERROR_DELTYPESIMPORREGEX);
                }
                else
                    throw OML_Error(HW_ERROR_DELTYPESIMPORREGEX);
            }
            else
                throw OML_Error("Invalid option '" + str + "'");
        }
        else
            throw OML_Error(HW_ERROR_OPTIONSTRING);
    }

    std::vector<std::string> stringPieces;
    size_t index = 0, lastIndex = 0;
    while (1)
    {
        index = std::string::npos;
        size_t delimLen = 1;

        for (int i = 0; i < delims.size(); i++)
        {
            std::string d = delims[i];
            size_t temp = tosplit.find(d, lastIndex);
            if (index == std::string::npos || (temp < index && temp != std::string::npos))
            {
                index = temp;
                delimLen = d.length();
            }
        }

        if (index == std::string::npos)
        {
            std::string sub = tosplit.substr(lastIndex, index);
            stringPieces.push_back(sub);
            break;
        }

        std::string sub = tosplit.substr(lastIndex, index - lastIndex);
        if (!lastIndex || !collapseDelimiters || sub.length())
            stringPieces.push_back(sub);
        lastIndex = index + delimLen;
    }

    HML_CELLARRAY *out = containerToCellArray(stringPieces, true);

    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Concatenates string inputs [strjoin]
//------------------------------------------------------------------------------
bool oml_strjoin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    if (size < 1 || size > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    hwMatrix delimstr(1, 1, hwMatrix::REAL);
    delimstr(0, 0) = ' ';

    Currency input1 = toCurrencyStr(eval, inputs[0], false, true);
    if (input1.IsCellArray())
    {
        HML_CELLARRAY *cell = input1.CellArray();
        int cellsize = cell->Size();
        if (!isstr(cell))
            throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

        if (!cell->Size())
        {
            outputs.push_back(Currency(""));
            return true;
        }
        else if (cell->Size() == 1)
        {
            outputs.push_back((*cell)(0));
            return true;
        }

        int numToConcat = 0;
        bool errorIfNoCell = false, singleDimStrings = true;
        int i;

        for (i = 0; i < cellsize; i++)
        {
            Currency elem = (*cell)(0);
            if (elem.IsString())
            {
                const hwMatrix *mtx = elem.Matrix();
                int m = mtx->M();

                if (!m)
                    continue;

                if (numToConcat)
                {
                    if (numToConcat != m)
                        errorIfNoCell = true;
                }
                else
                    numToConcat = m;

                if (m > 1)
                    singleDimStrings = false;
            }
            else if (elem.IsCellArray())
            {
                break;
            }
            else
                throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);
        }

        // if embedded cells
        if (i != cellsize)
        {
            Currency out((*cell)(0));
            for (i = 1; i < cellsize; i++)
            {
                dostrcat(eval, out, &delimstr); // add delimiter
                Currency elem = (*cell)(i);
                if (elem.IsString())
                    dostrcat(eval, out, elem.GetWritableMatrix());
                else if (elem.IsCellArray())
                    dostrcat(eval, out, elem.CellArray());
                else
                    throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);
            }

            outputs.push_back(out);
            return true;
        }

        if (errorIfNoCell)
            throw OML_Error(HW_ERROR_STRDIM);

        if (!numToConcat)
            numToConcat = cellsize;

        if (size > 1)
        {
            Currency input2 = inputs[1];
            if (input2.IsCellArray())
            {
                HML_CELLARRAY *delimcell = input2.CellArray();
                if (delimcell->Size() == 1)
                    input2 = (*delimcell)(0);
                else if (cellsize != delimcell->Size() + 1)
                    throw OML_Error(HW_ERROR_DELCELLINVAMTOFELE);
                else
                {
                    Currency out((*cell)(0));
                    for (int i = 1; i < cellsize; i++)
                    {
                        out = hconcat(out.Matrix(), (*delimcell)(i - 1).Matrix());
                        out = hconcat(out.Matrix(), (*cell)(i).Matrix());
                    }
                    out.SetMask(Currency::MASK_STRING);
                    outputs.push_back(out);
                    return true;
                }
            }
            if (input2.IsString())
            {
                if (singleDimStrings)
                    delimstr = *readRow(eval, input2).Matrix();
                else
                    delimstr = *input2.Matrix();
            }
            else
                throw OML_Error(HW_ERROR_DELONLYSTR);
        }

        if (!delimstr.IsEmpty() && (delimstr.M() != numToConcat))
            throw OML_Error(HW_ERROR_DELCELLINVAMTOFELE);

        Currency out((*cell)(0));
        for (int i = 1; i < cellsize; i++)
        {
            out = hconcat(out.Matrix(), &delimstr);
            out = hconcat(out.Matrix(), (*cell)(i).Matrix());
        }
        out.SetMask(Currency::MASK_STRING);
        outputs.push_back(out);
        return true;
    }
    else
        throw OML_Error(OML_ERR_CELL, 1, OML_VAR_TYPE);
}
//------------------------------------------------------------------------------
// Searches and replaces all instances of a pattern in input string [strrep]
//------------------------------------------------------------------------------
bool oml_strrep(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    if (size != 3 && size != 5)
        throw OML_Error(OML_ERR_NUMARGIN);

    bool overlap = true;
    const Currency &input1 = inputs[0];
    const Currency &input2 = inputs[1];
    const Currency &input3 = inputs[2];

    if (size > 3)
    {
        const Currency &input4 = inputs[3];
        const Currency &input5 = inputs[4];
        if (input4.IsString())
        {
            std::string str = readString(input4);
            if (str == "overlaps")
            {
                if (input5.IsScalar())
                {
                    overlap = (bool) input5.Scalar();
                }
                else if (input5.IsComplex())
                {
                    hwComplex cplx = input5.Complex();
                    overlap = cplx.Real() ? false : cplx.Imag() ? false : true;
                }
                else
                    throw OML_Error(HW_ERROR_INVINP5THARG);
            }
            else
                throw OML_Error(HW_ERROR_INVALIDOPTION(str));
        }
        else
            throw OML_Error(HW_ERROR_OPTIONSTRING);
    }

    outputs.push_back(tripleCurrencyFunc(eval, input1, input2, input3, &overlap, &doStrRep));
    return true;
}
//------------------------------------------------------------------------------
// Does a case-sensitive comparison of first n characters in inputs [strncmp]
//------------------------------------------------------------------------------
bool oml_strncmp(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    int n = (int) inputs[2].Scalar();

    if (!IsInteger(inputs[2].Scalar()).IsOk())
        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_VALUE);

    if (n <= 0)
        throw OML_Error(OML_ERR_POSINTEGER, 3);

    return dostrcmp(inputs, outputs, n);
}
//------------------------------------------------------------------------------
// Does a case-sensitive comparison of inputs [strcmp]
//------------------------------------------------------------------------------
bool oml_strcmp(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    return dostrcmp(inputs, outputs);
}
//------------------------------------------------------------------------------
// Concatenates input strings [strcat]
//------------------------------------------------------------------------------
bool oml_strcat(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency out(addStringMask(EvaluatorInterface::allocateMatrix()));
    for (int i = 0; i < inputs.size(); i++)
    {
        Currency input = inputs[i];
        if (input.IsCellArray())
        {
            HML_CELLARRAY *cell = input.CellArray();
            if (cell->Size() == 1)
            {
                Currency elem = (*cell)(0);
                if (elem.IsString())
                {
                    dostrcat(eval, out, elem.GetWritableMatrix());
                }
                else
                {
                    HML_CELLARRAY *wrapper = EvaluatorInterface::allocateCellArray(1, 1);
                    (*wrapper)(0) = elem;
                    dostrcat(eval, out, wrapper);
					delete wrapper;
                }
            }
            else
            {
                dostrcat(eval, out, cell);
            }
        }
        else
        {
            input = toCurrencyStr(eval, input, true, true);
            if (input.IsString())
            {
				hwMatrix* temp = trimright(input.Matrix());
                dostrcat(eval, out, temp);
				delete temp;
            }
            else
                throw OML_Error(HW_ERROR_NOTCONVINPTOSTR);
        }
    }

    outputs.push_back(out);

    return true;
}
//------------------------------------------------------------------------------
// Converts input string to lower case characters [lower]
//------------------------------------------------------------------------------
bool oml_lower(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    Currency currstr = toCurrencyStr(eval, input, true, true);
    outputs.push_back(convertToLower(currstr));
    return true;
}
//------------------------------------------------------------------------------
// Converts input string to upper case characters [upper]
//------------------------------------------------------------------------------
bool oml_upper(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    Currency currstr = toCurrencyStr(eval, input, true, true);
    outputs.push_back(convertToUpper(currstr));
    return true;
}
//------------------------------------------------------------------------------
// Returns true if successul in finding a pattern in given input [strfind]
//------------------------------------------------------------------------------
bool oml_strfind(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin != 2 && nargin != 4) throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    if (!input1.IsString() && !input1.IsCellArray()) 
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 1, OML_VAR_TYPE);

    const Currency &input2 = inputs[1];
    if (!input2.IsString() && !input2.IsCellArray()) 
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 2, OML_VAR_TYPE);

    bool allowOverlap = true;
    bool usecell = false;
    std::vector<double> idxs;
    std::vector<std::string> tosearch;
    std::vector<std::string> patterns;
    HML_CELLARRAY *cellidxs = EvaluatorInterface::allocateCellArray();

    if (nargin == 4)
    {
        const Currency &input3 = inputs[2];
        const Currency &input4 = inputs[3];
        if (!input3.IsString())
            throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);
        std::string i3 = readOption(eval, input3);
        if (i3 != "overlaps")
            throw OML_Error(HW_ERROR_INVALIDOPTION(i3));

        if (input4.IsScalar())
            allowOverlap = (bool) input4.Scalar();
        else if (!input4.IsComplex())
            throw OML_Error(HW_ERROR_4THINPBOOL);
    }

    if (input1.IsString())
        tosearch.push_back(orderedStringVal(input1));
    else if (input1.IsCellArray())
    {
        usecell = true;
        HML_CELLARRAY *cell = input1.CellArray();
        if (!isstr(cell))
            throw OML_Error(HW_ERROR_CELLELEMSTR);

        for (int i = 0; i < cell->Size(); i++)
            tosearch.push_back(orderedStringVal((*cell)(i)));

        delete cellidxs;
        cellidxs = EvaluatorInterface::allocateCellArray(cell->M(), cell->N());
    }

    if (input2.IsString())
        patterns.push_back(orderedStringVal(input2));
    else if (input2.IsCellArray())
    {
        HML_CELLARRAY *cell = input2.CellArray();
        if (usecell)
        {
            if (!sameSize(cell, cellidxs))
                throw OML_Error(HW_ERROR_CELLARRAYSSAMESIZE);
        }
        else
        {
            usecell = true;
            delete cellidxs;
            cellidxs = EvaluatorInterface::allocateCellArray(cell->M(), cell->N());
        }
        if (!isstr(cell))
            throw OML_Error(HW_ERROR_CELLELEMSTR);

        for (int i = 0; i < cell->Size(); i++)
            patterns.push_back(orderedStringVal((*cell)(i)));
    }

    bool foundPattern = false;
    if (!tosearch.empty() && !patterns.empty())
    {
        size_t maxval = max(tosearch.size(), patterns.size());
        for (int i = 0; i < maxval; i++)
        {
            idxs = std::vector<double>();

            std::string searchin = tosearch[tosearch.size() == 1 ? 0 : i];
            std::string pattern = patterns[patterns.size() == 1 ? 0 : i];
            size_t start = 0;

            if (!pattern.size())
            {
                for (int j = 0; j <= searchin.size(); j++)
                    idxs.push_back(j + 1);
            }
            else
            {
                while (1)
                {
                    size_t index = searchin.find(pattern, start);
                    if (index == std::string::npos)
                        break;

                    foundPattern = true;
                    if (allowOverlap)
                        start = index + 1;
                    else
                        start = index + pattern.size();

                    idxs.push_back((double) index + 1);
                }
            }

            hwMatrix *mtxidxs = NULL;
			
			if (!idxs.empty() && foundPattern) 
				mtxidxs = containerToMatrix(idxs);
			else // Returns an empty matrix if there is no match
				mtxidxs = EvaluatorInterface::allocateMatrix();

            if (usecell)
            {
                (*cellidxs)(i) = mtxidxs;
            }
            else
            {
                outputs.push_back(mtxidxs);
                return true;
            }
        }
    }

    outputs.push_back(cellidxs);
    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is a valid variable name [isvarname]
//------------------------------------------------------------------------------
bool oml_isvarname(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (!input.IsString())
        throw OML_Error(HW_ERROR_INPUTSTRING);

    std::string str = orderedStringVal(input);

    outputs.push_back(HW_MAKE_BOOL_CURRENCY(isvarname(str)));

    return true;
}
//------------------------------------------------------------------------------
// Returns the unique elements in the given input [unique]
//------------------------------------------------------------------------------
bool oml_unique(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    bool compareRows = false;
    bool forward = false;

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    if (nargin > 1)
    {
        const Currency &input2 = inputs[1];
        bool direcSet = false;
        if (!input2.IsString())
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);

        std::string i2 = readOption(eval, input2);
        if (i2 == "first")
            forward = direcSet = true;
        else if (i2 == "rows")
            compareRows = true;
        else if (i2 == "last")
            direcSet = true;
        else
            throw OML_Error(HW_ERROR_INPONLY1STLASTORROWS);

        if (nargin > 2)
        {
            const Currency &input3 = inputs[2];
            if (!input3.IsString())
                throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);

            std::string i3 = readOption(eval, input3);
            if (i3 == "first")
            {
                if (direcSet)
                    throw OML_Error(HW_ERROR_NOTSETDIRECTTWICE);
                forward = direcSet = true;
            }
            else if (i3 == "rows")
            {
                compareRows = true;
            }
            else if (i3 == "last")
            {
                if (direcSet)
                    throw OML_Error(HW_ERROR_NOTSETDIRECTTWICE);
            }
            else
                throw OML_Error(HW_ERROR_INPONLY1STLASTORROWS);
        }
    }

    int  nargout   = getNumOutputs(eval);
    bool outputIdx = (nargout > 1);
    bool inputIdx  = (nargout > 2);

    UniqueHelperFunc(eval, input1, compareRows, forward, outputIdx, inputIdx, outputs);
    return true;
}
//------------------------------------------------------------------------------
// Creates an empty cell array [cell]
//------------------------------------------------------------------------------
bool oml_cell(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();
    int m = 0, n = 0;

    if (size > 2)
	{
        throw OML_Error(OML_ERR_UNSUPPORTDIM);
	}
    else if (size == 1)
    {
        const Currency &input = inputs[0];
        if (input.IsMatrix())
        {
            const hwMatrix *mtx = input.Matrix();
            if (!mtx->IsRealData())
                throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_VARIABLE);

            if (!mtx->Size())
			{
                throw OML_Error(HW_MATH_MSG_EMPTYMATRIX);
			}
            else if (mtx->Size() == 1)
			{
                m = n = (int) realval(mtx, 0);
			}
            else if (mtx->Size() == 2)
            {
                m = (int) realval(mtx, 0);
                n = (int) realval(mtx, 1);
            }
			else
			{
				throw OML_Error(OML_ERR_UNSUPPORTDIM, 1);
			}
        }
        else if (input.IsScalar())
        {
            m = n = (int)input.Scalar();

            if (!IsInteger(input.Scalar()).IsOk())
                throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_VALUE);
        }
		else
		{
			throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_VALUE);
		}
    }
    else if (size == 2)
    {
        m = (int) inputs[0].Scalar();
        n = (int) inputs[1].Scalar();

        if (!IsInteger(inputs[0].Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_VALUE);

        if (!IsInteger(inputs[1].Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_VALUE);
    }

    HML_CELLARRAY* result = EvaluatorInterface::allocateCellArray(max(0, m), max(0, n));
    outputs.push_back(result);
    return true;
}
//------------------------------------------------------------------------------
// Sets the tic time [tic]
//------------------------------------------------------------------------------
bool oml_tic(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size())
        BuiltInFuncsUtils::SetWarning(eval, "Warning: ignoring arguments");

#if defined (OS_WIN)
    clock_t ct = clock();
#else
    #if defined(_SC_CLK_TCK)
        clock_t ct = times(nullptr);
    #else
        clock_t ct = clock();
    #endif
#endif

    if (getNumOutputs(eval) > 0)
        outputs.push_back((double)ct);
    else
    tictime = ct;

    return true;
}
//------------------------------------------------------------------------------
// Returns true and gets the time in seconds from last tic call [toc]
//------------------------------------------------------------------------------
bool oml_toc(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.empty() ? 0 : inputs.size();

    if (nargin > 1) throw OML_Error(OML_ERR_NUMARGIN);

    // Check if first input is valid
    double clocktime = 0;
    if (nargin == 1)
    {
        const Currency& in = inputs[0];
        if (!in.IsScalar() && !in.IsComplex())
            throw OML_Error(OML_ERR_SCALAR_COMPLEX, 1, OML_VAR_TYPE);

        clocktime = (in.IsScalar()) ? in.Scalar() : (in.Complex().Real());
    }

    clock_t tic = (nargin == 1) ? static_cast<clock_t>(clocktime) : tictime;

    if (tic < 0)
    {
        BuiltInFuncsUtils::SetWarning(eval, "Warning: tic must be called first");
        return true;
    }

#if defined (OS_WIN)
    double diff = (clock() - tic) / static_cast<double>(CLOCKS_PER_SEC);
#else
    #if defined(_SC_CLK_TCK)
        double diff = (times(nullptr) - tic) / static_cast<double>(sysconf(_SC_CLK_TCK));
    #else
        double diff = (clock() - tic) / static_cast<double>(CLOCKS_PER_SEC);
    #endif
#endif

    int numoutputs = eval.GetNargoutValue();
    if (numoutputs != 0)
            outputs.push_back(diff);
        else
    {
        std::ostringstream os;
        os << "Elapsed time is " << diff << " seconds.";
        eval.PrintResult(os.str());
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns the prime factorization of the given input [factor]
//------------------------------------------------------------------------------
bool oml_factor(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];
    bool useTwoVecs = getNumOutputs(eval) > 1;

    if (input.IsScalar())
    {
        std::vector<long long> primes;
        std::vector<int> multiplics;
        double inp = input.Scalar();
        if (!islonglong(inp))
            throw OML_Error(OML_ERR_INTEGER, -1, OML_VAR_VALUE);
        long long x = (long long) inp;

        if (x >= 2)
        {
            long long p = 2;
            while (p <= sqrt((double) x))
            {
                if (islonglong(x / (double) p)) // if p is factor of x
                {
                    x /= p;
                    addValMultiplicity(p, useTwoVecs, primes, multiplics);
                }
                else
                    p++;
            }
            if (x != 1)
                addValMultiplicity(x, useTwoVecs, primes, multiplics);
        }

        hwMatrix *primeMtx, *multMtx;
        if (primes.size())
        {
            primeMtx = EvaluatorInterface::allocateMatrix(1, (int)primes.size(), hwMatrix::REAL);
            for (int i = 0; i < primes.size(); i++)
                (*primeMtx)(i) = (double) primes[i];

            if (useTwoVecs)
            {
                multMtx = EvaluatorInterface::allocateMatrix(1, (int)multiplics.size(), hwMatrix::REAL);
                for (int i = 0; i < multiplics.size(); i++)
                    (*multMtx)(i) = multiplics[i];
            }
        }
        else
        {
            primeMtx = EvaluatorInterface::allocateMatrix(1, 1, inp);
            if (useTwoVecs)
                multMtx = EvaluatorInterface::allocateMatrix(1, 1, 1.0);
        }

        outputs.push_back(primeMtx);
        if (useTwoVecs)
            outputs.push_back(multMtx);
    }
    else if (input.IsComplex())
    {
        hwComplex x = input.Complex();
        if (!islonglong(x.Real()) || !islonglong(x.Imag()))
            throw OML_Error(OML_ERR_INTEGER, -1, OML_VAR_VALUE);

        throw OML_Error(OML_ERR_REAL, -1, OML_VAR_VALUE);
    }
    else
        throw OML_Error(OML_ERR_INTEGER, -1, OML_VAR_VALUE);

    return true;
}
//------------------------------------------------------------------------------
// Evaluates a polynomial [polyval]
//------------------------------------------------------------------------------
bool oml_polyval(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();
    if (size != 2 && size != 4) // the 3rd input is a sparse matrix -- no support yet
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    const Currency &input2 = inputs[1];

    if (!input2.IsMatrix() && !input2.IsScalar() && !input2.IsComplex() && !input2.IsString())
        throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);

    hwMatrix *m1, *mu;
    hwMatrix m2(*input2.ConvertToMatrix());
    const hwMatrix* m3;

    if (size == 4)
    {
        const Currency &input3 = inputs[2];

        if (!input3.IsMatrix() && !input3.IsScalar() && !input3.IsComplex() && !input3.IsString())
            throw OML_Error(OML_ERR_MATRIX, 3, OML_VAR_DATA);

        m3 = input3.ConvertToMatrix();

        if (inputs[3].IsMatrix())
        {
            mu = EvaluatorInterface::allocateMatrix(inputs[3].Matrix());
            if (mu->Size() < 2)
                throw OML_Error(OML_ERR_VECTOR2, 4);
        }
        else
            throw OML_Error(OML_ERR_VECTOR2, 4);

        if (!m3->IsEmpty())
            throw OML_Error(HW_ERROR_3RDINPMUSTEMPMATSPMATUNSUPP);
    }

    if (input1.IsScalar())
    {
        double x = input1.Scalar();
        hwMatrix *result = EvaluatorInterface::allocateMatrix(m2.M(), m2.N(), x);
        outputs.push_back(result);
        return true;
    }
    else if (input1.IsComplex())
    {
        hwComplex cplx = input1.Complex();
        hwMatrix *result = EvaluatorInterface::allocateMatrix(m2.M(), m2.N(), cplx);
        outputs.push_back(result);
        return true;
    }
    else if (input1.IsMatrix())
    {
        m1 = EvaluatorInterface::allocateMatrix(input1.Matrix());
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

    bool complex;
    if (size == 4)
        complex = checkMakeComplex(eval, m1, &m2, mu);
    else
        complex = checkMakeComplex(eval, m1, &m2);

    int numCoefs = m1->Size();

    if (input1.IsMatrix() && !m1->IsVector())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    if (m2.IsEmpty())
        throw OML_Error(HW_MATH_MSG_EMPTYMATRIX);

    if (size == 4)
    {
        if (complex)
        {
            hwComplex c1 = mu->z(0);
            hwComplex c2 = mu->z(1);
            for (int i = 0; i < m2.Size(); i++)
            {
                hwComplex &xval = m2.z(i);
                xval = (xval - c1) / c2;
            }
        }
        else
        {
            double d1 = (*mu)(0);
            double d2 = (*mu)(1);
            for (int i = 0; i < m2.Size(); i++)
            {
                double &xval = m2(i);
                xval = (xval - d1) / d2;
            }
        }
    }

    hwMatrix *y = EvaluatorInterface::allocateMatrix(m2.M(), m2.N(), complex ? hwMatrix::COMPLEX : hwMatrix::REAL);
    if (complex)
    {
        hwComplex xx, yy;
        for (int i = 0; i < m2.Size(); i++)
        {
            xx = m2.z(i);
            int j = 0;
            yy = m1->z(j);
            while (j < numCoefs - 1)
            {
                j++;
                yy = yy * xx + m1->z(j);
            }
            y->z(i) = yy;
        }
    }
    else
    {
        double xx, yy;
        for (int i = 0; i < m2.Size(); i++)
        {
            xx = m2(i);
            int j = 0;
            yy = (*m1)(j);
            while (j < numCoefs - 1)
            {
                j++;
                yy = yy * xx + (*m1)(j);
            }
            (*y)(i) = yy;
        }
    }
    outputs.push_back(y);

    return true;
}
//------------------------------------------------------------------------------
// Returns an upper triangular matrix [triu]
//------------------------------------------------------------------------------
bool oml_triu(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();
    if (size < 1 || size > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency input1, input2;
    input1 = inputs[0];
    int offset = 0;

    if (!input1.IsMatrix() && !input1.IsScalar() && !input1.IsComplex() && !input1.IsString())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);

    hwMatrix* mtx = new hwMatrix(*input1.ConvertToMatrix());
    Currency out(mtx);

    if (size > 1)
    {
        // verify input
        offset = (int) inputs[1].Scalar();

        if (!IsInteger(inputs[1].Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_OFFSET);

        if (offset > 0 && offset > mtx->N())
            throw OML_Error(OML_ERR_INVALID_RANGE, 2, OML_VAR_OFFSET);
        else if (-1 * offset > mtx->M())
            throw OML_Error(OML_ERR_INVALID_RANGE, 2, OML_VAR_OFFSET);
    }

    if (mtx->IsReal())
    {
        for (int i = 0; i < mtx->M(); i++)
        {
            for (int j = 0; j < i + offset && j < mtx->N(); j++)
            {
                if (j >= 0)
                    (*mtx)(i, j) = 0.0;
            }
        }
    }
    else
    {
        for (int i = 0; i < mtx->M(); i++)
        {
            for (int j = 0; j < i + offset && j < mtx->N(); j++)
            {
                if (j >= 0)
                    mtx->z(i, j) = 0.0;
            }
        }
    }

    out.SetMask(inputs[0].GetMask());
    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Returns a lower triangular matrix [tril]
//------------------------------------------------------------------------------
bool oml_tril(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    if (size < 1 || size > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency input1, input2;
    input1 = inputs[0];
    int offset = 0;

    if (!input1.IsMatrix() && !input1.IsScalar() && !input1.IsComplex() && !input1.IsString())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);

    hwMatrix* mtx = new hwMatrix(*input1.ConvertToMatrix());
    Currency out(mtx);

    if (size > 1)
    {
        // verify input
        offset = (int) inputs[1].Scalar();

        if (!IsInteger(inputs[1].Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_OFFSET);

        if (offset > 0 && offset > mtx->N())
            throw OML_Error(OML_ERR_INVALID_RANGE, 2, OML_VAR_OFFSET);
        else if (-1 * offset > mtx->M())
            throw OML_Error(OML_ERR_INVALID_RANGE, 2, OML_VAR_OFFSET);
    }

    if (mtx->IsReal())
    {
        for (int i = 0; i < mtx->M(); i++)
        {
            for (int j = i + offset + 1; j < mtx->N(); j++)
            {
                if (j >= 0)
                    (*mtx)(i, j) = 0.0;
            }
        }
    }
    else
    {
        for (int i = 0; i < mtx->M(); i++)
        {
            for (int j = i + offset + 1; j < mtx->N(); j++)
            {
                if (j >= 0)
                    mtx->z(i, j) = 0.0;
            }
        }
    }

    out.SetMask(inputs[0].GetMask());
    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Returns true in input is a cell array of strings [iscellstr]
//------------------------------------------------------------------------------
bool oml_iscellstr(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<
    Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsCellArray())
    {
        HML_CELLARRAY *cell = input.CellArray();
        outputs.push_back(HW_MAKE_BOOL_CURRENCY(isstr(cell)));
    }
    else
        outputs.push_back(getFalse());

    return true;
}
//------------------------------------------------------------------------------
// Polynomial construction, either as a characteristic polynomial, or from its roots [poly]
//------------------------------------------------------------------------------
bool oml_poly(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input = inputs[0];

    if (!input.IsMatrix() && !input.IsScalar() && !input.IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);

    const hwMatrix* x = input.ConvertToMatrix();
    std::unique_ptr<hwMatrix> p, e;

    if (!x->Size())
    {
        outputs.push_back(1.0);
        return true;
    }

    int m = x->M(), n = x->N();

    if (m == n)
    {
        std::unique_ptr<hwMatrix> V(EvaluatorInterface::allocateMatrix());
        e.reset(EvaluatorInterface::allocateMatrix());

        BuiltInFuncsUtils::CheckMathStatus(eval, x->Eigen(false, V.get(), *e));
        // workaround for a bug
        if (!e->IsReal() && e->IsRealData())
        {
            // convert to real
            hwMatrix *temp = EvaluatorInterface::allocateMatrix(e->M(), e->N(), hwMatrix::REAL);
            for (int i = 0; i < e->Size(); i++)
            {
                (*temp)(i) = e->z(i).Real();
            }
            e.reset(temp);
        }
    }
    else if (input.IsVector())
    {
        e.reset(EvaluatorInterface::allocateMatrix(x));
    }
    else
        throw OML_Error(HW_ERROR_INPSQUAREMATVEC);

    hwMatrix *idx = checkMatrixFinite(e.get());
    Currency idxcur(idx);
    if (e->IsReal())
    {
        std::vector<double> nonzeros;
        for (int i = 0; i < idx->Size(); i++)
        {
            if (!iszero((*idx)(i)))
                nonzeros.push_back((*e)(i));
        }

        e->Dimension((const int)nonzeros.size(), 1, hwMatrix::REAL);

        for (int i = 0; i < nonzeros.size(); i++)
        {
            e->SetElement(i, nonzeros[i]);
        }
    }
    else
    {
        std::vector<hwComplex> nonzeros;
        for (int i = 0; i < idx->Size(); i++)
        {
            if (!iszero((*idx)(i)))
                nonzeros.push_back(e->z(i));
        }

        BuiltInFuncsUtils::CheckMathStatus(eval, e->Dimension((const int)nonzeros.size(), 1, hwMatrix::COMPLEX));

        for (int i = 0; i < nonzeros.size(); i++)
        {
            e->SetElement(i, nonzeros[i]);
        }
    }

    n = max(e->M(), e->N());
    p.reset(EvaluatorInterface::allocateMatrix(1, n + 1, hwMatrix::REAL));
    (*p)(0) = 1.0;
    for (int i = 0; i < n; i++)
    {
        p->SetElement(i + 1, 0.0);
    }

    if (e->IsReal())
    {
        for (int k = 0; k < n; k++)
        {
            double *pvals = new double[k + 1];

            for (int j = 0; j <= k; j++)
            {
                pvals[j] = (*p)(j);
            }

            for (int j = 1; j <= k + 1; j++)
            {
                (*p)(j) -= (*e)(k) * pvals[j - 1];
            }

            delete [] pvals;
        }
    }
    else
    {
        if (p->IsReal())
            BuiltInFuncsUtils::CheckMathStatus(eval, p->MakeComplex());
        for (int k = 0; k < n; k++)
        {
            hwComplex *pvals = new hwComplex[k + 1];

            for (int j = 0; j <= k; j++)
            {
                pvals[j] = p->z(j);
            }

            for (int j = 1; j <= k + 1; j++)
            {
                p->z(j) -= e->z(k) * pvals[j - 1];
            }

            delete [] pvals;
        }
    }

    if (!e->IsReal())
    {
        std::deque<hwComplex> posimag;
        std::deque<hwComplex> negimag;

        for (int i = 0; i < e->Size(); i++)
        {
            hwComplex c = e->z(i);
            hwComplex c2 = c.Conjugate();
            if (c.Imag() > 0)
                posimag.push_back(c);
            else if (c.Imag() < 0)
                negimag.push_back(c2);
        }

        std::sort(posimag.begin(), posimag.end(), &complexLessThan);
        std::sort(negimag.begin(), negimag.end(), &complexLessThan);

        if (posimag.size() == negimag.size())
        {
            int i;
            for (i = 0; i < posimag.size(); i++)
            {
                if (posimag[i] != negimag[i])
                    break;
            }

            // roots are complex conjugates of each other
            if (!p->IsReal() && i == posimag.size())
            {
                hwMatrix *newp = EvaluatorInterface::allocateMatrix(p->M(), p->N(), hwMatrix::REAL);
                for (int j = 0; j < newp->Size(); j++)
                {
                    (*newp)(j) = p->z(j).Real();
                }
                outputs.push_back(newp);
                return true;
            }
        }
    }

    outputs.push_back(p.release());
    return true;
}
//------------------------------------------------------------------------------
// Rounds input to 0 [fix]
//------------------------------------------------------------------------------
bool oml_fix(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];
    if (input.IsScalar())
    {
        outputs.push_back((int) input.Scalar());
    }
    else if (input.IsComplex())
    {
        hwComplex cplx = input.Complex();
        cplx.Set((int) cplx.Real(), (int) cplx.Imag());
        outputs.push_back(cplx);
    }
    else if (input.IsMatrix())
    {
        const hwMatrix *mtx = input.Matrix();
        hwMatrix *result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
        if (mtx->IsReal())
        {
            for (int i = 0; i < result->Size(); i++)
            {
                (*result)(i) = (int) (*mtx)(i);
            }
        }
        else
        {
            for (int i = 0; i < result->Size(); i++)
            {
                hwComplex cplx = mtx->z(i);
                cplx.Set((int) cplx.Real(), (int) cplx.Imag());
                result->z(i) = cplx;
            }
        }
        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_fix);
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

    return true;
}
//------------------------------------------------------------------------------
// The resulting remainder of inputs x / y [rem]
//------------------------------------------------------------------------------
bool oml_rem(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    const Currency &input2 = inputs[1];

    if (input1.IsLogical() || input2.IsLogical())
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    if (input1.IsScalar() || (input1.IsComplex() && iszero(input1.Complex().Imag())))
    {
        double val1 = input1.IsScalar() ? input1.Scalar() : input1.Complex().Real();
        if (input2.IsScalar())
        {
            outputs.push_back(doubleMod(input1.Scalar(), input2.Scalar()));
            return true;
        }
        else if (input2.IsMatrix())
        {
            const hwMatrix *mtx = input2.Matrix();
            hwMatrix *result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);
            for (int i = 0; i < result->Size(); i++)
            {
                (*result)(i) = doubleMod(val1, (*mtx)(i));
            }
            outputs.push_back(result);
            return true;
        }
        else if (inputs[1].IsNDMatrix())
        {
            return oml_MatrixNUtil2(eval, inputs, outputs, oml_rem);
        }
        else
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

    }
    else if (input1.IsMatrix())
    {
        const hwMatrix *m1 = input1.Matrix();
        if (m1->IsRealData())
        {
            if (input2.IsScalar())
            {
                double val = input2.Scalar();
                hwMatrix *result = EvaluatorInterface::allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);
                for (int i = 0; i < result->Size(); i++)
                {
                    (*result)(i) = doubleMod(realval(m1, i), val);
                }
                outputs.push_back(result);
                return true;
            }
            else if (input2.IsMatrix())
            {
                const hwMatrix *m2 = input2.Matrix();
                if (m2->IsRealData())
                {
                    if (!sameSize(m1, m2))
                        throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

                    hwMatrix *result = EvaluatorInterface::allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);
                    for (int i = 0; i < result->Size(); i++)
                    {
                        (*result)(i) = doubleMod(realval(m1, i), realval(m2, i));
                    }
                    outputs.push_back(result);
                    return true;
                }
                else
                    throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
            }
            else
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);
        }
    }
    else if (input1.IsNDMatrix() || input2.IsNDMatrix())
    {
        return oml_MatrixNUtil2(eval, inputs, outputs, oml_rem);
    }

    throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);
}
//------------------------------------------------------------------------------
// Reassigns dimensions of input matrix, retaining data [reshape]
//------------------------------------------------------------------------------
bool oml_reshape(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];
    const Currency& input2 = inputs[1];

    // determine output object type: 2D or ND
    bool NDout = false;

    if (nargin == 2)
    {
        if (input2.IsVector())
        {
            if (input2.Matrix()->Size() > 2)
                NDout = true;
        }
        else
        {
            throw OML_Error(OML_ERR_NNINTVECTOR, 2, OML_VAR_DIMS);
        }
    }
    else if (nargin > 3)
    {
        NDout = true;
    }

    // handle output object type cases
    if (!NDout)     // 2D case
    {
        int m;
        int n;
    
        if (nargin == 2)
        {
            if (input2.IsVector())
            {
                const hwMatrix* dims = input2.Matrix();

                if (!dims->IsRealData())
                    throw OML_Error(OML_ERR_REALVECTOR, 2, OML_VAR_DIMS);

                if (!IsInteger(realval(dims, 0)).IsOk())
                    throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIMS);

                if (!IsInteger(realval(dims, 1)).IsOk())
                    throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIMS);

                m = static_cast<int> (realval(dims, 0));
                n = static_cast<int> (realval(dims, 1));
            }
        }
        else
        {
            const Currency& input3 = inputs[2];

            if (input2.IsPositiveInteger())
                m = static_cast<int>(input2.Scalar());
            else if (input2.IsMatrix() && !input2.Matrix()->M() && !input2.Matrix()->N())
                m = -1;
            else
                throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);

            if (input3.IsPositiveInteger())
                n = static_cast<int>(input3.Scalar());
            else if (input3.IsMatrix() && !input3.Matrix()->M() && !input3.Matrix()->N())
                n = -1;
            else
                throw OML_Error(OML_ERR_NATURALNUM, 3, OML_VAR_DIM);
        }

        if (input1.IsCellArray())
        {
            std::unique_ptr<HML_CELLARRAY> cell(EvaluatorInterface::allocateCellArray(input1.CellArray()));
            BuiltInFuncsUtils::CheckMathStatus(eval, cell->Reshape(m, n));
            outputs.push_back(cell.release());
        }
        else if (input1.IsStruct())
        {
            std::unique_ptr<StructData> strct(EvaluatorInterface::allocateStruct(input1.Struct()));
            BuiltInFuncsUtils::CheckMathStatus(eval, strct->Reshape(m, n));
            outputs.push_back(strct.release());
        }
        else if (input1.IsNDMatrix())
        {
            std::unique_ptr<hwMatrixN> mtxN(EvaluatorInterface::allocateMatrixN(input1.MatrixN()));
            std::vector<int> dims(2);
            dims[0] = m;
            dims[1] = n;

            try
            {
                mtxN->Reshape(dims);
            }
            catch (hwMathException& e)
            {
                e.Status().ResetArgs();
                throw e;
            }

            std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix());
            mtxN->ConvertNDto2D(*mtx);
            Currency out(mtx.release());
            out.SetMask(input1.GetMask());
            outputs.push_back(out);
        }
        else    // hwMatrix, including string
        {
            std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix(input1.ConvertToMatrix()));
            hwMathStatus status = mtx->Reshape(m, n);

            if (!status.IsOk())
            {
                status.ResetArgs();
            }

            BuiltInFuncsUtils::CheckMathStatus(eval, status);
            Currency out(mtx.release());
            out.SetMask(input1.GetMask());
            outputs.push_back(out);
        }
    }
    else    // ND case
    {
        std::vector<int> dims;

        if (nargin == 2)
        {
            const hwMatrix* indx = input2.Matrix();

            for (int i = 0; i < indx->Size(); ++i)
            {
                if (IsInteger((*indx)(i)).IsOk())
                    dims.push_back(static_cast<int>((*indx)(i)));
                else
                    throw OML_Error(OML_ERR_NNINTVECTOR, 2, OML_VAR_DIMS);
            }
        }
        else
        {
            for (int i = 1; i < nargin; ++i)
            {
                const Currency& cur = inputs[i];

                if (cur.IsPositiveInteger())
                    dims.push_back(static_cast<int>(cur.Scalar()));
                else if (inputs[i].IsMatrix() && !inputs[i].Matrix()->M() && !inputs[i].Matrix()->N())
                    dims.push_back(-1);
                else
                    throw OML_Error(OML_ERR_NATURALNUM, i+1, OML_VAR_DIM);
            }
        }

        if (input1.IsNDMatrix())
        {
            std::unique_ptr<hwMatrixN> mtxN(EvaluatorInterface::allocateMatrixN(input1.MatrixN()));

            try
            {
                mtxN->Reshape(dims);
            }
            catch (hwMathException& e)
            {
                e.Status().ResetArgs();
                throw e;
            }

            Currency out(mtxN.release());
            out.SetMask(input1.GetMask());
            outputs.push_back(out);
        }
        else    // hwMatrix, including string
        {
            std::unique_ptr<hwMatrixN> mtxN(EvaluatorInterface::allocateMatrixN());
            mtxN->Convert2DtoND(*input1.Matrix());

            try
            {
                mtxN->Reshape(dims);
            }
            catch (hwMathException& e)
            {
                e.Status().ResetArgs();
                throw e;
            }

            Currency out(mtxN.release());
            out.SetMask(input1.GetMask());
            outputs.push_back(out);
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the generalized transpose for a matrix A using the vector P [permute]
//------------------------------------------------------------------------------
bool oml_permute(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];
    const Currency& input2 = inputs[1];

    // check permutation vector
    if (!input2.IsMatrix())
    {
        throw OML_Error(OML_ERR_POSINTVECTOR, 2, OML_VAR_DIMS);
    }

    const hwMatrix* pvec = input2.Matrix();

    if (!pvec->IsVector() || !pvec->IsRealData())
    {
        throw OML_Error(OML_ERR_POSINTVECTOR, 2, OML_VAR_DIMS);
    }

    std::vector<int> permvec(pvec->Size());

    for (int i = 0; i < pvec->Size(); ++i)
    {
        double dim = realval(pvec, i);

        // if (!isposint(dim))
        if (dim < 1)
        {
            throw OML_Error(hwMathStatus(HW_MATH_ERR_PERMVEC2, 2));
        }

        permvec[i] = static_cast<int>(dim) - 1;
    }

    // determine output object type: 2D or ND
    bool NDout = false;

    if (input1.IsNDMatrix())
        NDout = true;
    else if (permvec.size() > 2)
        NDout = true;

    // handle output object type cases
    if (!NDout)     // 2D case
    {
        if (!input1.IsMatrix() && !input1.IsScalar() && !input1.IsComplex())
            throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_TYPE);

        if (permvec.size() != 2)
            throw OML_Error(hwMathStatus(HW_MATH_ERR_PERMVEC1, 2));

        const hwMatrix* matrix = input1.ConvertToMatrix();

        int dim1 = permvec[0];
        int dim2 = permvec[1];

        if ((dim1 != 0 && dim1 != 1) || (dim2 != 0 && dim2 != 1) || dim1 == dim2)
            throw OML_Error(hwMathStatus(HW_MATH_ERR_PERMVEC3, 2));

        std::unique_ptr<hwMatrix> result(EvaluatorInterface::allocateMatrix());

        if (dim1 == 1)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, result->Transpose(*input1.Matrix()));
            outputs.push_back(result.release());
        }
        else
        {
            outputs.push_back(input1);
        }
    }
    else    // ND case
    {
        if (input1.IsScalar() || input1.IsComplex())
        {
            outputs.push_back(input1);
            return true;
        }

        std::unique_ptr<hwMatrixN> result(EvaluatorInterface::allocateMatrixN());

        if (input1.IsMatrix())
        {
            hwMatrixN temp;

            temp.Convert2DtoND(*input1.ConvertToMatrix());

            result->Permute(temp, permvec);
        }
        else if (input1.IsNDMatrix())
        {
            result->Permute(*input1.MatrixN(), permvec);
        }
        else
            throw OML_Error(OML_ERR_MATRIX, 1);

        outputs.push_back(result.release());
    }

    return true;
}
//------------------------------------------------------------------------------
// Removes singular dimensions from a matrix [squeeze]
//------------------------------------------------------------------------------
bool oml_squeeze(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    if (size != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];

    if (input1.IsMatrix())
    {
        outputs.push_back(input1);
    }
    else if (input1.IsNDMatrix())
    {
        const hwMatrixN* mtxN = input1.MatrixN();
        const std::vector<int>& dims = mtxN->Dimensions();
        std::vector<int> newDims;

        for (int i = 0; i < dims.size(); ++i)
        {
            if (dims[i] != 1)
            {
                newDims.push_back(dims[i]);
            }
        }

        if (newDims.size() == 0 || newDims.size() == dims.size())
        {
            outputs.push_back(input1);
        }
        else
        {
            if (newDims.size() == 1)
                newDims.push_back(1);

            hwMatrixN* out = new hwMatrixN(*mtxN);

            out->Reshape(newDims);
            outputs.push_back(out);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
    }

    return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
static double normForVector(const hwMatrix* vec, double p)
{
    double sum = 0.0;
    if (vec->IsReal())
    {
        for (int i = 0; i < vec->Size(); ++i)
            sum += pow(abs((*vec)(i)), p); 
    }
    else
    {
        for (int i = 0; i < vec->Size(); ++i)
            sum += pow(vec->z(i).Mag(), p); 
    }
    return pow(sum, 1/p);
}
//------------------------------------------------------------------------------
// Computes matrix and vector norms [norm]
//------------------------------------------------------------------------------
bool oml_norm(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();
    if (size < 1 || size > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    const Currency *input3 = size > 2 ? &inputs[2] : nullptr;

    if (input1.IsLogical())
    {
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);
    }
    else if (input1.IsScalar())
    {
        outputs.push_back(abs(input1.Scalar()));
    }
    else if (input1.IsComplex())
    {
        outputs.push_back(abs(input1.Complex().Mag()));
    }
    else if (input1.IsMatrix())
    {
        const hwMatrix *mtx = input1.Matrix();
        double p = 2.0;
        std::string strp, opt;
        bool usedouble = true;
        bool userows;

        if (size > 1)
        {
            const Currency& input2 = inputs[1];
            if (input2.IsScalar())
            {
                p = input2.Scalar();
            }
            else if (input2.IsComplex())
            {
                p = input2.Complex().Real();
            }
            else if (input2.IsString())
            {
                strp = readOption(eval, input2);
                if (!(strp == "inf" || strp == "fro" || strp == "-inf"))
                {
                    if (size == 3)
                        throw OML_Error(HW_ERROR_INVSTRVALINFNEGINFFRO);
                    // emulate size of 3
                    size = 3;
                    input3 = &inputs[1];
                }
                else
                    usedouble = false;
            }
            else
            {
                throw OML_Error(HW_ERROR_INPSCSTR);
            }

            if (usedouble)
            {
                if (p == std::numeric_limits<double>::infinity())
                {
                    strp = std::string("inf");
                    usedouble = false;
                }
                else if (p == -1 * std::numeric_limits<double>::infinity())
                {
                    strp = std::string("-inf");
                    usedouble = false;
                }
            }

            if (size > 2)
            {
                if (!input3->IsString())
                    throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);

                opt = readOption(eval, *input3);
                if (opt == "rows")
                {
                    userows = true;
                }
                else if (opt == "cols" || opt == "columns")
                {
                    userows = false;
                }
                else
                    throw OML_Error(HW_ERROR_INVOPTMUSTROWORCOL);
            }
        }

        if (p < 0.0 && strp != "-inf")
            throw OML_Error(OML_ERR_PNORM, 2);

        // check 3rd argument first
        double norm;
        if (size < 3)
        {
            if (usedouble)
            {
                if (mtx->IsVector())
                {
                    if (isint(p))
                    {
                        if ((int) p == 0)
                        {
                            int size = mtx->Size();
                            int count = 0;

                            if (mtx->IsReal())
                            {
                                for (int i = 0; i < size; ++i)
                                {
                                    if ((*mtx)(i) != 0.0)
                                        ++count;
                                }
                            }
                            else
                            {
                                for (int i = 0; i < size; ++i)
                                {
                                    if (mtx->z(i) != 0.0)
                                        ++count;
                                }
                            }

                            norm = (double) count;
                        }
                        else
                        {
                            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Norm(norm, (int)p));
                        }
                    }
                    else
                    {
                        norm = normForVector(mtx, p);
                    }
                }
                else
                {
                    if (!isint(p))
                        throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_VALUE); // requires constrained optimization

                    BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Norm(norm, (int) p));
                }
            }
            else
                BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Norm(norm, strp.c_str()));

            outputs.push_back(norm);
            return true;
        }

        hwMatrix *result;
        Currency resultCur;
        if (userows)
        {
            result = EvaluatorInterface::allocateMatrix(mtx->M(), 1, hwMatrix::REAL);
            std::unique_ptr<hwMatrix> row(EvaluatorInterface::allocateMatrix());

            for (int i = 0; i < mtx->M(); i++)
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, mtx->ReadRow(i, *row));
                if (usedouble)
                {
                    if (isint(p))
                    {
                        if ((int) p == 0)
                        {
                            int size = row.get()->Size();
                            int count = 0;

                            if (mtx->IsReal())
                            {
                                for (int i = 0; i < size; ++i)
                                {
                                    if ((*row.get())(i) != 0.0)
                                        ++count;
                                }
                            }
                            else
                            {
                                for (int i = 0; i < size; ++i)
                                {
                                    if (row.get()->z(i) != 0.0)
                                        ++count;
                                }
                            }

                            norm = (double) count;
                        }
                        else
                        {
                            BuiltInFuncsUtils::CheckMathStatus(eval, row.get()->Norm(norm, (int)p));
                        }
                    }
                    else
                    {
                        norm = normForVector(row.get(), p);
                    }
                }
                else
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, row->Norm(norm, strp.c_str()));
                }

                result->SetElement(i, norm);
            }
        }
        else
        {
            result = EvaluatorInterface::allocateMatrix(1, mtx->N(), hwMatrix::REAL);

            for (int i = 0; i < mtx->N(); i++)
            {
                std::unique_ptr<const hwMatrix> col(EvaluatorInterface::allocateColumn(mtx, i));
                if (usedouble)
                {
                    if (isint(p))
                    {
                        if ((int) p == 0)
                        {
                            int size = col.get()->Size();
                            int count = 0;

                            if (mtx->IsReal())
                            {
                                for (int i = 0; i < size; ++i)
                                {
                                    if ((*col.get())(i) != 0.0)
                                        ++count;
                                }
                            }
                            else
                            {
                                for (int i = 0; i < size; ++i)
                                {
                                    if (col.get()->z(i) != 0.0)
                                        ++count;
                                }
                            }

                            norm = (double) count;
                        }
                        else
                        {
                            BuiltInFuncsUtils::CheckMathStatus(eval, col.get()->Norm(norm, (int)p));
                        }
                    }
                    else
                    {
                        norm = normForVector(col.get(), p);
                    }
                }
                else
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, col->Norm(norm, strp.c_str()));
                }

                result->SetElement(i, norm);
            }
        }
        outputs.push_back(result);
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }
    return true;
}
//------------------------------------------------------------------------------
// Imaginary number representation [i]
//------------------------------------------------------------------------------
bool oml_i(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return createCommonMatrix(eval, inputs, outputs, Currency(hwComplex(0.0, 1.0)));
}
//------------------------------------------------------------------------------
// Returns the position of non-zero elements in a matrix [find]
//------------------------------------------------------------------------------
bool oml_find(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs)
{
    size_t size = inputs.empty() ? 0: inputs.size();
    if (size == 0 || size > 3) throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    const hwMatrix* findin = 0;   

    if (input1.IsMatrix() || input1.IsScalar() || input1.IsComplex() || input1.IsString())
        findin = input1.ConvertToMatrix();
    else if (input1.IsNDMatrix())
    {
        std::vector<int> newidx(2);

        newidx[0] = inputs[0].MatrixN()->Dimensions()[0];
        newidx[1] = -1;

        hwMatrixN matrix(*inputs[0].MatrixN());
        matrix.Reshape(newidx);

        hwMatrix* matrix2D = new hwMatrix;
        matrix.ConvertNDto2D(*matrix2D);

        std::vector<Currency> inputs2;
        inputs2.push_back(matrix2D);

        for (int i = 1; i < size; ++i)
            inputs2.push_back(inputs[i]);

        return oml_find(eval, inputs2, outputs);
    }
    else
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX);

    int nargout = getNumOutputs(eval);
    int incr = 1;
    double stopAt = 0;

    if (size > 1)
    {
        if (!inputs[1].IsScalar())
            throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

        stopAt = inputs[1].Scalar();
        if (stopAt < 0 || !isint(stopAt)) throw OML_Error(HW_ERROR_NUMVALFINDPOSINT);

        if (size > 2)
        {
            const Currency &input3 = inputs[2];
            if (!input3.IsString()) throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);

            std::string direc = readOption(eval, input3);
            if (direc == "last")
                incr = -1;
            else if (direc != "first")
                throw OML_Error(HW_ERROR_DIRECT1STLAST);
        }
    }

    std::vector<int> ivec; // Vector of indices

    if (findin->IsReal())
    {
        for (int i = (incr == 1 ? 0 : findin->Size() - 1); incr == 1 ? i < findin->Size() : i >= 0;
            i += incr)
        {
            if ((*findin)(i) != 0.0)
            {
                if (size > 1 && ivec.size() >= stopAt)
                    break;
                ivec.push_back(i);
            }
        }
    }
    else
    {
        for (int i = (incr == 1 ? 0 : findin->Size() - 1); incr == 1 ? i < findin->Size() : i >= 0;
            i += incr)
        {
            if (findin->z(i) != hwComplex(0.0, 0.0))
            {
                if (size > 1 && ivec.size() >= stopAt)
                    break;
                ivec.push_back(i);
            }
        }
    }

    if (ivec.empty())  // No pattern match found
    {
		int output_rows = 1;

		if ((findin->M() == 0) && (findin->N() == 0))
			output_rows = 0;

		if (nargout)
		{
			for (int j=0; j<nargout; j++)
				outputs.push_back(EvaluatorInterface::allocateMatrix(output_rows, 0, hwMatrix::REAL));
		}
		else
		{
			outputs.push_back(EvaluatorInterface::allocateMatrix(output_rows, 0, hwMatrix::REAL));
		}

        return true;
    }

    size_t vecsize = ivec.size();

    if (nargout == 1 || !nargout)
    {
        hwMatrix *indices;
        if (findin->M() == 1)
            indices = EvaluatorInterface::allocateMatrix(1, (int)vecsize, hwMatrix::REAL);
        else
            indices = EvaluatorInterface::allocateMatrix((int)vecsize, 1, hwMatrix::REAL);

        for (int i = 0; i < vecsize; i++)
        {
            int idx = ivec[incr == 1 ? i : vecsize-i-1];
            (*indices)(i) = static_cast<double>(idx + 1);
        }
        outputs.push_back(indices);
    }
    else
    {
        hwMatrix *rows, *cols, *vals;
        if (findin->M() == 1)
        {
            rows = EvaluatorInterface::allocateMatrix(1, (int)vecsize, hwMatrix::REAL);
            cols = EvaluatorInterface::allocateMatrix(1, (int)vecsize, hwMatrix::REAL);
            vals = EvaluatorInterface::allocateMatrix(1, (int)vecsize, findin->Type());
        }
        else
        {
            rows = EvaluatorInterface::allocateMatrix((int)vecsize, 1, hwMatrix::REAL);
            cols = EvaluatorInterface::allocateMatrix((int)vecsize, 1, hwMatrix::REAL);
            vals = EvaluatorInterface::allocateMatrix((int)vecsize, 1, findin->Type());
        }

        int m = findin->M();
        if (findin->IsReal())
        {
            for (int i = 0; i < vecsize; i++)
            {
                int idx;
                if (incr == -1)
                    idx = ivec[vecsize - i - 1];
                else
                    idx = ivec[i];

                (*rows)(i) = static_cast<double>(idx % m + 1);
                (*cols)(i) = static_cast<double>(idx / m + 1);
                (*vals)(i) = (*findin)(idx);
            }
        }
        else
        {
            for (int i = 0; i < vecsize; i++)
            {
                int idx;
                if (incr == -1)
                    idx = ivec[vecsize - i - 1];
                else
                    idx = ivec[i];

                (*rows)(i) = static_cast<double>(idx % m + 1);
                (*cols)(i) = static_cast<double>(idx / m + 1);
                vals->z(i) = findin->z(idx);
            }
        }

        outputs.push_back(rows);
        outputs.push_back(cols);
        outputs.push_back(vals);
    }

    return true;
}
//------------------------------------------------------------------------------
// Evaluates given input command [eval]
//------------------------------------------------------------------------------
bool oml_eval(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();
    if (size != 1 && size != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &trycur = inputs[0];
    if (!trycur.IsString())
        throw OML_Error(HW_ERROR_TRYSTR);

    std::string trystr = readString(trycur);
    std::string errmsg;
    Interpreter interp(eval);
    Interpreter interp2(eval);
    Interpreter* pInterp = &interp;
    Currency ret;

    if (trystr.length())
	{
		eval.CacheLineInfomation();
        ret = interp.DoString(trystr);
	}

    if (ret.IsError())
    {
        // if there are any results from this operation, we need to carry them over to the parent
        // and clean out everything in the child -- or do we?
        if (size == 2)
        {
            const Currency &catchcur = inputs[1];
            pInterp = &interp2; // use second interpreter so we don't mix outputs
            if (catchcur.IsString())
            {
                std::string catchstr = readString(catchcur);
                if (catchstr.length())
                {
                    ret = pInterp->DoString(catchstr);
                    if (ret.IsError())
                        errmsg = ret.Message();
                }
                else
                    ret = Currency();
            }
            else
            {
                ret = Currency();
                errmsg = std::string("Error: catch expression must be a string");
        }
        }
        else
            errmsg = ret.Message();
    }

    std::vector<Currency> results = pInterp->GetOutputCurrencyList();
    std::vector<Currency>::iterator it = results.begin(), end = results.end();

    bool outputSomething = getNumOutputs(eval) && results.size();

    if (outputSomething)
        end--;

    for (; it != end; it++)
    {
        if (it->IsError())
            throw OML_Error(it->Message());
        eval.PushResult(*it);
    }

    // to output errors from interp.DoString()
    if (errmsg.length())
        throw OML_Error(errmsg);

    if (outputSomething && !it->IsError())
        outputs.push_back(*it);

	eval.UncacheLineInfomation();

    return true;
}
//------------------------------------------------------------------------------
// Outputs the largest of the dimensions of the input
//------------------------------------------------------------------------------
size_t utf8_get_char_size(unsigned char* ptr)
{
    if ((*ptr & 0x80) == 0)
        return 1;

    if ((*ptr & 0x40) == 0)
        return 1;

    if ((*ptr & 0x20) == 0)
        return 2;

    if ((*ptr & 0x10) == 0)
        return 3;

    if ((*ptr & 0x08) == 0)
        return 4;

    if ((*ptr & 0x04) == 0)
        return 5;

    if ((*ptr & 0x02) == 0)
        return 6;

    return 1;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
unsigned char* utf8_increment_pointer(unsigned char* ptr)
{
    return ptr + utf8_get_char_size(ptr);
}
//------------------------------------------------------------------------------
// Gets length of utf string
//------------------------------------------------------------------------------
size_t utf8_strlen(unsigned char* ptr)
{
    size_t len = 0;

    for ( ; *ptr; ++len)
        ptr = utf8_increment_pointer(ptr);

    return len;
}
//------------------------------------------------------------------------------
// Returns length of input [length]
//------------------------------------------------------------------------------
bool oml_length(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];
    if (input.IsScalar() || input.IsComplex() || input.IsFunctionHandle())
        outputs.push_back(1.0);
    else if (input.IsCellArray())
    {
        HML_CELLARRAY *cells = input.CellArray();

		if (cells->IsEmpty())
			outputs.push_back(0.0);
		else
		   outputs.push_back(max(cells->M(), cells->N()));
    }
    else if (input.IsMatrix())
    {
        const hwMatrix *mtx = input.Matrix();

        if (mtx->IsEmpty())
            outputs.push_back(0.0);
        else
            outputs.push_back(max(mtx->M(), mtx->N()));
    }
    else if (input.IsNDMatrix())
    {
        const hwMatrixN *mtx = input.MatrixN();

        if (mtx->IsEmpty())
        {
            outputs.push_back(0.0);
        }
        else
        {
            const std::vector<int>& dims = mtx->Dimensions();
            int length = dims[0];

            for (int i = 1; i < dims.size(); ++i)
            {
                if (dims[i] > length)
                    length = dims[i];
            }

            outputs.push_back(length);
        }
    }
	else if (input.IsString())
	{
		if (input.IsMultilineString())
		{
		    const hwMatrix *mtx = input.Matrix();
			int max = mtx->M();
			
			if (mtx->N() > max)
				max = mtx->N();

			outputs.push_back(max);
		}
		else
		{
			std::string st = input.StringVal();
			unsigned char* my_ptr = (unsigned char*)st.c_str();
			outputs.push_back(utf8_strlen(my_ptr));
		}
	}
    else if (input.IsStruct())
    {
        StructData *sd = input.Struct();
        outputs.push_back((double) max(sd->M(), sd->N()));
    }
    else
        throw OML_Error(HW_ERROR_INVINPTYPE);

    return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
template <typename T>
static void dosize(const T *inp, int dim, int nargin, int nargout, std::vector<Currency> &outputs)
{
    if (nargin == 1)
    {
        if (nargout < 2)
        {
            hwMatrix *result = EvaluatorInterface::allocateMatrix(1, 2, hwMatrix::REAL);

            if (inp)
            {
                result->SetElement(0, inp->M());
                result->SetElement(1, inp->N());
            }
            else
            {
                result->SetElement(0, 0.0);
                result->SetElement(1, 0.0);
            }
            outputs.push_back(result);
        }
        else
        {
			if (inp)
			{
				outputs.push_back(inp->M());
				outputs.push_back(inp->N());
			}
			else
			{
				outputs.push_back(0.0);
				outputs.push_back(0.0);
			}

            for (int i = 2; i < nargout; ++i)
                outputs.push_back(1.0);
        }
    }
    else // size == 2
    {	
        if ((int) dim == 1)
		{
			if (inp)
				outputs.push_back(inp->M());
			else
				outputs.push_back(0.0);
		}
        else if ((int) dim == 2)
		{
			if (inp)
				outputs.push_back(inp->N());
			else
				outputs.push_back(0.0);
		}
        else // dim > 2
		{
            outputs.push_back(1.0);
		}
    }
}
//------------------------------------------------------------------------------
// Returns size of input [size]
//------------------------------------------------------------------------------
bool oml_size(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();
    int dim = 0;

    if (size < 1 || size > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (size == 2)
    {
        if (!inputs[1].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DIM);

        dim = (int) inputs[1].Scalar();
    }

    const Currency &input1 = inputs[0];
    int nargout = getNumOutputs(eval);

    if (input1.IsScalar() || input1.IsComplex() || input1.IsFunctionHandle())
    {
        if (size == 1)
        {
            if (nargout < 2)
            {
                hwMatrix *result = EvaluatorInterface::allocateMatrix(1, 2, 1.0);
                outputs.push_back(result);
            }
            else
            {
                for (int i = 0; i < nargout; ++i)
                    outputs.push_back(1.0);
            }
        }
        else // size == 2
        {
            outputs.push_back(1.0);
        }
    }
	else if (input1.IsNDMatrix())
	{
		const hwMatrixN* mtxn = input1.MatrixN();
		std::vector<int> dims = mtxn->Dimensions();

        if (size == 1)
        {
            if (nargout < 2)
            {
                hwMatrix *result = EvaluatorInterface::allocateMatrix(1, (int)dims.size(), 1.0);

                for (int i = 0; i < dims.size(); ++i)
                   (*result)(i) = dims[i];

                outputs.push_back(result);
            }
            else
            {
                if (nargout <= dims.size())
                {
                    for (int i = 0; i < nargout - 1; ++i)
                        outputs.push_back(dims[i]);

                    int prod = dims[nargout-1];

                    for (int i = nargout; i < dims.size(); ++i)
                        prod *= dims[i];

                    outputs.push_back(prod);
                }
                else
                {
                    for (int i = 0; i < static_cast<int>(dims.size()); ++i)
                        outputs.push_back(dims[i]);

                    for (int i = static_cast<int>(dims.size()); i < nargout; ++i)
                        outputs.push_back(1.0);
                }
            }
        }
        else // size == 2
        {
            if (dim <= dims.size())
                outputs.push_back(dims[dim-1]);
            else
                outputs.push_back(1.0);
        }
	}
    else if (input1.IsMatrix() || input1.IsString())
        dosize(input1.Matrix(), dim, (int) size, nargout, outputs);
    else if (input1.IsCellArray())
        dosize(input1.CellArray(), dim, (int) size, nargout, outputs);
    else if (input1.IsStruct())
        dosize(input1.Struct(), dim, (int) size, nargout, outputs);
    else
        throw OML_Error(HW_ERROR_INPUTSTRCELLMTX);

    return true;
}
//------------------------------------------------------------------------------
// Returns the number of dimensions of the input matrix [ndims]
//------------------------------------------------------------------------------
bool oml_ndims(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

	if (inputs[0].IsNDMatrix())
	{
		const hwMatrixN* mtxn = inputs[0].MatrixN();
		const std::vector<int> dims = mtxn->Dimensions();
		outputs.push_back(dims.size());
	}
	else
	{
		outputs.push_back(2.0);
	}

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is complex [iscomplex]
//------------------------------------------------------------------------------
bool oml_iscomplex(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(HW_MAKE_BOOL_CURRENCY(inputs[0].GetMask() == Currency::MASK_EXPLICIT_COMPLEX || isActuallyComplex(inputs[0])));
    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is real [isreal]
//------------------------------------------------------------------------------
bool oml_isreal(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(HW_MAKE_BOOL_CURRENCY(inputs[0].GetMask() != Currency::MASK_EXPLICIT_COMPLEX && isActuallyReal(inputs[0])));
    return true;
}
//------------------------------------------------------------------------------
// Returns true if inputs are equal [isequal]
//------------------------------------------------------------------------------
bool oml_isequal(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() < 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &curInput = inputs[0];
    Currency out(1.0);
    for (int i = 1; i < inputs.size(); i++)
    {
        if (!isequal(curInput, inputs[i]))
        {
            out = 0.0;
            break;
        }
    }

    out.SetMask(Currency::MASK_LOGICAL);
    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is empty [isempty]
//------------------------------------------------------------------------------
bool oml_isempty(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(HW_MAKE_BOOL_CURRENCY(isempty(inputs[0])));
    return true;
}
//------------------------------------------------------------------------------
// Compute the inverse matrix [inv]
//------------------------------------------------------------------------------
bool oml_inv(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];
    if (input.IsScalar())
    {
        outputs.push_back(1 / input.Scalar());
        outputs.push_back(1.0);
    }
    else if (input.IsComplex())
    {
        outputs.push_back(1.0 / input.Complex());
        outputs.push_back(1.0);
    }
    else if (input.IsMatrix())
    {
        const hwMatrix *mtx = input.Matrix();
        hwMatrix *result = EvaluatorInterface::allocateMatrix();
        hwMathStatus stat = result->Inverse(*mtx);
        BuiltInFuncsUtils::CheckMathStatus(eval, stat);
        outputs.push_back(result);

        double rcond;
        stat = mtx->RCond(rcond);
        outputs.push_back(rcond);
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

    return true;
}
//------------------------------------------------------------------------------
// Returns the valie of the specified environment variable [getenv]
//------------------------------------------------------------------------------
bool oml_getenv(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    std::string varname = readString(toCurrencyStr(eval, inputs[0], false, false));
    char *env = getenv(varname.c_str());
    if (env == nullptr)
        outputs.push_back(std::string());
    else
        outputs.push_back(std::string(env));
    return true;
}
//------------------------------------------------------------------------------
// Returns the current working directory [pwd]
//------------------------------------------------------------------------------
bool oml_pwd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    outputs.push_back(BuiltInFuncsUtils::GetCurrentWorkingDir());
    return true;
}
//------------------------------------------------------------------------------
// Returns the current date and time [clock]
//------------------------------------------------------------------------------
bool oml_clock(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    hwMatrix *result = EvaluatorInterface::allocateMatrix(1, 6, hwMatrix::REAL);

    time_t t = time(0);
    struct tm *currentTime = localtime(&t);

    result->SetElement(0, currentTime->tm_year + 1900);
    result->SetElement(1, currentTime->tm_mon + 1);
    result->SetElement(2, currentTime->tm_mday);
    result->SetElement(3, currentTime->tm_hour);
    result->SetElement(4, currentTime->tm_min);
    result->SetElement(5, currentTime->tm_sec);

    outputs.push_back(result);
    return true;
}
//------------------------------------------------------------------------------
// Returns a matrix with rows and columns of roughly equal magnitude [balance]
//------------------------------------------------------------------------------
bool oml_balance(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size(), nargout = getNumOutputs(eval);
    if (size < 1 || size > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    bool noperm = true;
    const Currency &input1 = inputs[0];
    Currency input2, input3;

    if (size > 1)
    {
        input2 = inputs[1];
        if (size > 2)
            input3 = inputs[2];
    }

    if (input1.IsScalar() || input1.IsComplex())
    {
        if (size > 1)
        {
            if (input2.IsScalar() || input2.IsComplex())
            {
                if (!nargout)
                    throw OML_Error(OML_ERR_NUMARGOUT);
                else if (nargout == 1)
                    BuiltInFuncsUtils::SetWarning(eval, "Warning: there should be at least 2 outputs");

                outputs.push_back(1.0);
                outputs.push_back(1.0);
                outputs.push_back(input1);
                outputs.push_back(input2);
            }
            else if (input2.IsMatrix())
                throw OML_Error(HW_MATH_MSG_ARRAYSIZE);
            else if ((!input2.IsString() && size == 2) || (input2.IsString() && size == 3))
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

            if (size == 3 && !inputs[2].IsString())
                throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);
        }

        if (size == 1 || input2.IsString())
        {
            if (nargout > 1)
            {
                outputs.push_back(1.0);
                if (nargout > 2)
                    outputs.push_back(1.0);
            }
            outputs.push_back(input1);
        }
        return true;
    }
    else if (input1.IsMatrix())
    {
        if (size > 1)
        {
            input2 = inputs[1];
            if (input2.IsScalar() || input2.IsComplex())
                throw OML_Error(HW_MATH_MSG_ARRAYSIZE);
            else if (!input2.IsMatrix())
            {
                std::string str;
                if (input2.IsString())
                {
                    if (size == 3)
                        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
                    else
                        str = readString(input2);
                }
                else if (input3.IsString())
                    str = readString(input3);
                else
                    throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);

                if (str == "p" || str == "noscal")
                    noperm = false;
                else if (!(str == "s" || str == "noperm"))
                    throw OML_Error(HW_ERROR_INVSTRVAL);
            }
        }

        const hwMatrix *A = input1.Matrix();
        hwMatrix *AA = EvaluatorInterface::allocateMatrix();
        hwMatrix *P = EvaluatorInterface::allocateMatrix();
        hwMatrix *D = EvaluatorInterface::allocateMatrix();
        Currency aacur(AA), pcur(P), dcur(D);

        hwMathStatus stat = A->Balance(noperm, *D, *P, *AA);
        BuiltInFuncsUtils::CheckMathStatus(eval, stat);

        hwMatrix *DD = vectorToDiag(D);
        Currency ddcur(DD);

        if (size == 1 || input2.IsString())
        {
            if (nargout < 2)
            {
                outputs.push_back(aacur);
            }
            else if (nargout == 2)
            {
                outputs.push_back(ddcur);
                outputs.push_back(aacur);
            }
            else
            {
                outputs.push_back(dcur);
                outputs.push_back(pcur);
                outputs.push_back(aacur);
            }

            return true;
        }

        throw OML_Error(HW_MATH_MSG_NOTIMPLEMENT);

        /*
        hwMatrix *DDI = EvaluatorInterface::allocateMatrix();
        hwMatrix *AI = EvaluatorInterface::allocateMatrix();
        hwMatrix *B = input2.Matrix();
        Currency ddicur(DDI), aicur(AI);

        if (!sameSize(A, B))
            throw OML_Error(HW_MATH_MSG_ARRAYSIZE);

        if (!nargout)
            throw OML_Error(OML_ERR_NUMARGOUT);
        else if (nargout == 1)
            warning(eval, "there should be at least 2 outputs");

        stat = DDI->Inverse(*DD);
        checkMathStatus(eval, stat);
        stat = AI->Inverse(*A);
        checkMathStatus(eval, stat);

        hwMatrix *CC = EvaluatorInterface::allocateMatrix();
        Currency cccur(CC);
        CC->Mult(*DDI, *AI);
        CC->Mult(*AA, hwMatrix(*CC));
        if (CC->IsReal())
        {
            for (int i=0; i<CC->M(); i++)
            {
                for (int j=0; j<CC->N(); j++)
                {
                    if (i != j)
                        (*CC)(i,j) = 0.0;
                }
            }
        }
        else
        {
            for (int i=0; i<CC->M(); i++)
            {
                for (int j=0; j<CC->N(); j++)
                {
                    if (i != j)
                        CC->z(i,j) = hwComplex(0.0, 0.0);
                }
            }
        }

        hwMatrix *BB = EvaluatorInterface::allocateMatrix();
        Currency bbcur(BB);
        stat = A->Balance(noperm, *D, *P, *BB);
        checkMathStatus(eval, stat);

        outputs.push_back(ccur);
        outputs.push_back(ddcur);
        outputs.push_back(aacur);
        outputs.push_back(bbcur);

        return true;
        */
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
}
//------------------------------------------------------------------------------
// Returns the elements that are in the input matrices [union]
//------------------------------------------------------------------------------
bool oml_union(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return sortBasedOperation(eval, inputs, outputs, false, &dounion, &dounion, &dounion, &dounion);
}
//------------------------------------------------------------------------------
// Returns all primes up to given input [primes]
//------------------------------------------------------------------------------
bool oml_primes(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    int max;
    const Currency &input = inputs[0];

    if (input.IsScalar())
        max = (int)input.Scalar();
    else if (input.IsComplex())
        max = (int)input.Complex().Real();
    else
        throw OML_Error(OML_ERR_NATURALNUM);

    std::deque<int> p = getprimes(max);

    hwMatrix *result = EvaluatorInterface::allocateMatrix(1, (int)p.size(), hwMatrix::REAL);
    for (int i = 0; i < p.size(); i++)
    {
        result->SetElement(i, p[i]);
    }

    outputs.push_back(result);
    return true;
}
//------------------------------------------------------------------------------
// Returns the element-wise square root of the sum of the squares of inputs [hypot]
//------------------------------------------------------------------------------
bool oml_hypot(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    const Currency &input2 = inputs[1];

    if (input1.IsLogical() || input2.IsLogical())
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    double    singleton     = 0.0;
    bool      singlePresent = false;
    hwMatrix* m1            = nullptr;
    hwMatrix* m2            = nullptr;
    hwMatrix* result        = nullptr;

    if (input1.IsScalar())
    {
        singleton = input1.Scalar();
        singlePresent = true;
    }
    else if (input1.IsComplex())
    {
        singleton = input1.Complex().Mag();
        singlePresent = true;
    }
    else if (input1.IsMatrix())
    {
        m1 = EvaluatorInterface::allocateMatrix(input1.Matrix());
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

    if (input2.IsScalar())
    {
        if (singlePresent)
        {
            outputs.push_back(hypot(singleton, input2.Scalar()));
            return true;
        }
        singlePresent = true;
        singleton = input2.Scalar();
    }
    else if (input2.IsComplex())
    {
        if (singlePresent)
        {
            outputs.push_back(hypot(singleton, input2.Complex().Mag()));
            return true;
        }
        singlePresent = true;
        singleton = input2.Complex().Mag();
    }
    else if (input2.IsMatrix())
    {
        if (singlePresent)
            m1 = EvaluatorInterface::allocateMatrix(input2.Matrix());
        else
            m2 = EvaluatorInterface::allocateMatrix(input2.Matrix());
    }
    else
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

    if (singlePresent)
    {
        result = EvaluatorInterface::allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);
        if (m1->IsReal())
        {
            for (int i = 0; i < m1->Size(); i++)
            {
                result->SetElement(i, hypot(singleton, (*m1)(i)));
            }
        }
        else
        {
            for (int i = 0; i < m1->Size(); i++)
            {
                result->SetElement(i, hypot(singleton, m1->z(i).Mag()));
            }
        }
    }
    else
    {
        if (sameSize(m1, m2))
        {
            result = EvaluatorInterface::allocateMatrix(m1->M(), m1->N(), hwMatrix::REAL);
            if (m1->IsReal())
            {
                if (m2->IsReal())
                {
                    for (int i = 0; i < m1->Size(); i++)
                    {
                        result->SetElement(i, hypot((*m1)(i), (*m2)(i)));
                    }
                }
                else
                {
                    for (int i = 0; i < m1->Size(); i++)
                    {
                        result->SetElement(i, hypot((*m1)(i), m2->z(i).Mag()));
                    }
                }
            }
            else
            {
                if (m2->IsReal())
                {
                    for (int i = 0; i < m1->Size(); i++)
                    {
                        result->SetElement(i, hypot(m1->z(i).Mag(), (*m2)(i)));
                    }
                }
                else
                {
                    for (int i = 0; i < m1->Size(); i++)
                    {
                        result->SetElement(i, hypot(m1->z(i).Mag(), m2->z(i).Mag()));
                    }
                }
            }
        }
        else
            throw OML_Error(HW_MATH_MSG_ARRAYSIZE);
    }

    outputs.push_back(result);
    return true;
}
//------------------------------------------------------------------------------
// Returns epsilon [eps]
//------------------------------------------------------------------------------
bool oml_eps(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    bool usedouble = true;
    bool scalarOut = false;

    // determine output value
    if (nargin < 2)
        scalarOut = true;

    if (nargin)
    {
        const Currency &last = inputs[nargin - 1];
        if (last.IsString())
        {
            --nargin;
            std::string className = readOption(eval, last);
            if (className == "float" || className == "single")
                usedouble = false;
            else if (className != "double")
                throw OML_Error(HW_ERROR_INVCLASSNAME);
        }
    }

    double x;
    double epsilonVal;

    if (scalarOut && nargin == 1)
    {
        if (!inputs[0].IsScalar())
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

        x = inputs[0].Scalar();
    }
    else
    {
        x = 1.0;
    }

    if (usedouble)
    {
        epsilonVal = MachPrecision((double) x);
    }
    else
    {
        epsilonVal = (double) MachPrecision((float) x);
    }

    // determine object type
    bool createND = false;

    if (nargin > 2)
        createND = true;

    // construct matrix
    if (!createND)
    {
        if (scalarOut)
        {
            outputs.push_back(epsilonVal);
        }
        else 
        {
            int m, n;

            if (!inputs[0].IsInteger())
                throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);

            m = static_cast<int> (inputs[0].Scalar());
            m = max(m, 0);

            if (nargin == 1)
            {
                n = m;
            }
            else    // nargin == 2
            {
                if (!inputs[1].IsInteger())
                    throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);

                n = static_cast<int> (inputs[1].Scalar());
                n = max(n, 0);
            }

            hwMatrix *out = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            out->SetElements(epsilonVal);
            outputs.push_back(out);
        }
    }
    else
    {
		std::vector<int> dimensions;
        double dim;

        for (int j = 0; j < nargin; ++j)
        {
            if (inputs[j].IsScalar())
            {
                dim = inputs[j].Scalar();
            }
            else if (inputs[j].IsComplex())
            {
                dim = inputs[j].Complex().Real();
            }
            else
            {
                throw OML_Error(OML_ERR_NATURALNUM, j+1, OML_VAR_DIM);
            }

            if (!(::isfinite(dim)))
                throw OML_Error(HW_ERROR_DIMFINITE);

            if (dim < 0.0)
                dim = 0.0;

            dimensions.push_back(static_cast<int>(dim));
        }

        hwMatrixN* out = EvaluatorInterface::allocateMatrixN();
        out->Dimension(dimensions, hwMatrixN::REAL);
        out->SetElements(epsilonVal);
        outputs.push_back(out);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns a complex number [complex]
//------------------------------------------------------------------------------
bool oml_complex(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    Currency result;

    if (nargin == 1)
    {
        if (input1.IsScalar())
        {
            result = hwComplex(input1.Scalar(), 0.0);
        }
        else if (input1.IsComplex())
        {
            result = input1.Complex();
        }
        else if (input1.IsMatrix())
        {
            hwMatrix *resmtx = EvaluatorInterface::allocateMatrix(input1.Matrix());
            result = resmtx;
            if (resmtx->IsReal())
                BuiltInFuncsUtils::CheckMathStatus(eval, resmtx->MakeComplex());
        }
        else
        {
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
        }
    }
    else
    {
        result = mtxFun(eval, inputs, 1, &makeAComplex)[0];
    }

    result.SetMask(Currency::MASK_EXPLICIT_COMPLEX);
    outputs.push_back(result);
    return true;
}
//------------------------------------------------------------------------------
// Returns a matrix whose elements are 1 [ones]
//------------------------------------------------------------------------------
bool oml_ones(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	return createCommonMatrix(eval, inputs, outputs, Currency(1.0));
}
//------------------------------------------------------------------------------
// Returns pi [pi]
//------------------------------------------------------------------------------
bool oml_pi(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return createCommonMatrix(eval, inputs, outputs, Currency(PI));
}
//------------------------------------------------------------------------------
// Symbol for e [e]
//------------------------------------------------------------------------------
bool oml_e(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return createCommonMatrix(eval, inputs, outputs, Currency(exp(1.0)));
}
//------------------------------------------------------------------------------
// Returns a matrix whose elements are all zeros [zeros]
//------------------------------------------------------------------------------
bool oml_zeros(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	return createCommonMatrix(eval, inputs, outputs, Currency(0.0));
}
//------------------------------------------------------------------------------
// Returns the min value in the given input [min]
//------------------------------------------------------------------------------
bool oml_min(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();
    int dim = 1;
    Currency input1, input2, input3;

    if (size < 1 || size > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    input1 = inputs[0];
    if (size > 1)
    {
        input2 = inputs[1];
        if (size > 2)
            input3 = inputs[2];
    }

    if (size == 3)
    {
        if (!input2.IsMatrix())
            throw OML_Error(OML_ERR_EMPTYMATRIX, 2, OML_VAR_INPUT);

        if (input2.Matrix()->M() || input2.Matrix()->N())
            throw OML_Error(OML_ERR_EMPTYMATRIX, 2, OML_VAR_INPUT);

        if (!input3.IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);

        dim = static_cast<int>(input3.Scalar());
    }

    if (size == 1)
    {
        if (input1.IsScalar() || input1.IsComplex())
        {
            outputs.push_back(input1);
            outputs.push_back(1.0);
        }
        else if (input1.IsMatrix() || input1.IsString())
        {
            const hwMatrix *mtx = input1.Matrix();
            int index = 0;
            if (mtx->IsVector())
            {
                if (mtx->IsReal())
                {
                    double loc_min = (*mtx)(0);
                    for (int k = 1; k < mtx->Size(); k++)
                    {
                        double temp = (*mtx)(k);
                        if (temp < loc_min)
                        {
                            loc_min = temp;
                            index = k;
                        }
                    }
                    outputs.push_back(loc_min);
                    outputs.push_back(index + 1);
                }
                else
                {
                    hwComplex loc_min = mtx->z(0);
                    for (int k = 1; k < mtx->Size(); k++)
                    {
                        hwComplex temp = mtx->z(k);
                        if (complexLessThan(temp, loc_min))
                        {
                            loc_min = temp;
                            index = k;
                        }
                    }
                    outputs.push_back(loc_min);
                    outputs.push_back(index + 1);
                }
            }
            else
            {
                if (mtx->M() == 0)
                {
                    outputs.push_back(EvaluatorInterface::allocateMatrix(0, mtx->N(), hwMatrix::REAL)); // maximums
                    outputs.push_back(EvaluatorInterface::allocateMatrix(0, mtx->N(), hwMatrix::REAL)); // indices
                }
                else
                {
                    hwMatrix *result = EvaluatorInterface::allocateMatrix(1, mtx->N(), mtx->Type());
                    hwMatrix *indices = EvaluatorInterface::allocateMatrix(1, mtx->N(), mtx->Type());

                    if (mtx->IsReal())
                    {
                        double loc_min;
                        for (int i = 0; i < mtx->N(); ++i)
                        {
                            loc_min = (*mtx)(0, i);
                            index = 0;
                            for (int j = 1; j < mtx->M(); ++j)
                            {
                                double temp = (*mtx)(j, i);
                                if (temp < loc_min)
                                {
                                    loc_min = temp;
                                    index = j;
                                }
                            }
                            result->SetElement(i, loc_min);
                            indices->SetElement(i, index + 1);
                        }
                    }
                    else
                    {
                        hwComplex min;
                        for (int i = 0; i < mtx->N(); ++i)
                        {
                            min = mtx->z(0, i);
                            index = 0;
                            for (int j = 1; j < mtx->M(); ++j)
                            {
                                hwComplex temp = mtx->z(j, i);
                                if (complexLessThan(temp, min))
                                {
                                    min = temp;
                                    index = j;
                                }
                            }
                            result->SetElement(i, min);
                            indices->SetElement(i, index + 1);
                        }
                    }
                    outputs.push_back(result);
                    outputs.push_back(indices);
                }
            }
        }
        else if (input1.IsNDMatrix())
        {
            return oml_MatrixNUtil3(eval, inputs, outputs, oml_min);
        }
        else
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }
    else if (input1.IsNDMatrix() && size == 3)
    {
        return oml_MatrixNUtil3(eval, inputs, outputs, oml_min, 3);
    }
    else if (input1.IsNDMatrix() || input2.IsNDMatrix())
    {
        return oml_MatrixNUtil2(eval, inputs, outputs, oml_min);
    }
    else if (!((input1.IsMatrix() || input1.IsString()) || (input2.IsMatrix() || input2.IsString())))
    {
        if (input1.IsScalar())
        {
            if (input2.IsScalar())
            {
                outputs.push_back(min(input1.Scalar(), input2.Scalar()));
            }
            else if (input2.IsComplex())
            {
                double di1 = input1.Scalar();
                hwComplex i1 = hwComplex(di1, 0.0);
                hwComplex i2 = input2.Complex();
                if (complexLessThan(i1, i2))
                    outputs.push_back(di1);
                else
                    outputs.push_back(i2);
            }
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else if (input1.IsComplex())
        {
            if (input2.IsScalar())
            {
                hwComplex i1 = input1.Complex();
                double di2 = input2.Scalar();
                hwComplex i2 = hwComplex(di2, 0.0);
                if (complexLessThan(i2, i1))
                    outputs.push_back(di2);
                else
                    outputs.push_back(i1);
            }
            else if (input2.IsComplex())
            {
                hwComplex i1 = input1.Complex(), i2 = input2.Complex();
                if (complexLessThan(i1, i2))
                    outputs.push_back(i1);
                else
                    outputs.push_back(i2);
            }
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else
        {
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
        }
    }
    else if ((input1.IsMatrix() || input1.IsString()) && (input2.IsMatrix() || input2.IsString()))
    {
        const hwMatrix *i1 = input1.Matrix(), *i2 = input2.Matrix();

        if (i1->M() != i2->M() || i1->N() != i2->N()) // matrices have different dimensions
        {
            if (!i2->M() && !i2->N() && size == 3 && dim < 3)
            {
                hwMatrix *result;
                int m = i1->M(), n = i1->N();
                if (dim == 1)
                    result = EvaluatorInterface::allocateMatrix(1, n, i1->Type());
                else
                    result = EvaluatorInterface::allocateMatrix(m, 1, i1->Type());
                hwMatrix* indices = EvaluatorInterface::allocateMatrix(result->M(), result->N(), hwMatrix::REAL);

                if (i1->IsReal())
                {
                    double loc_min;
                    int min_index;
                    for (int i = 0; i < (dim == 1 ? n : m); ++i)
                    {
                        min_index = 0;
                        if (dim == 1)
                            loc_min = (*i1)(0, i);
                        else
                            loc_min = (*i1)(i, 0);

                        for (int j = 1; j < (dim == 1 ? m : n); ++j)
                        {
                            double val = (dim == 1) ? (*i1)(j,i) : (*i1)(i,j);
                            if (val < loc_min)
                            {
                                loc_min = val;
                                min_index = j;
                            }
                        }
                        result->SetElement(i, loc_min);
                        (*indices)(i) = min_index + 1.0;
                    }
                }
                else
                {
                    hwComplex min;
                    int min_index;
                    for (int i = 0; i < (dim == 1 ? n : m); ++i)
                    {
                        min_index = 0;
                        if (dim == 1)
                            min = i1->z(0, i);
                        else
                            min = i1->z(i, 0);

                        for (int j = 1; j < (dim == 1 ? m : n); ++j)
                        {
                            if (dim == 1)
                            {
                                hwComplex val = i1->z(j, i);
                                if (complexLessThan(val, min))
                                {
                                    min = val;
                                    min_index = j;
                                }
                            }
                            else
                            {
                                hwComplex val = i1->z(i, j);
                                if (complexLessThan(val, min))
                                {
                                    min = val;
                                    min_index = j;
                                }
                            }
                        }
                        result->SetElement(i, min);
                        (*indices)(i) = min_index + 1.0;
                    }
                }
                outputs.push_back(result);
                outputs.push_back(indices);
            }
            else if (!i2->M() && !i2->N() && size == 3 && dim > 2)
            {
                outputs.push_back(input1);
            }
            else
            {
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);
            }
        }
        else // matrices have same dimensions
        {
            if (!i1->Size())
            {
                outputs.push_back(EvaluatorInterface::allocateMatrix(0, 0, hwMatrix::REAL));
            }
            else
            {
                hwMatrix *result = EvaluatorInterface::allocateMatrix(i1->M(), i1->N(), i1->Type());

                if (i1->IsReal())
                {
                    if (i2->IsReal())
                    {
                        for (int k = 0; k < i1->Size(); k++)
                        {
                            result->SetElement(k, min((*i1)(k), (*i2)(k)));
                        }
                    }
                    else
                    {
                        for (int k = 0; k < i1->Size(); k++)
                        {
                            double v1 = (*i1)(k);
                            hwComplex val1 = hwComplex(v1, 0.0);
                            hwComplex val2 = i2->z(k);
                            if (complexLessThan(val1, val2))
                                result->SetElement(k, v1);
                            else
                                result->SetElement(k, val2);
                        }
                    }
                }
                else
                {
                    if (i2->IsReal())
                    {
                        for (int k = 0; k < i1->Size(); k++)
                        {
                            hwComplex val1 = i1->z(k);
                            double v2 = (*i2)(k);
                            hwComplex val2 = hwComplex(v2, 0.0);
                            if (complexLessThan(val2, val1))
                                result->SetElement(k, v2);
                            else
                                result->SetElement(k, val1);
                        }
                    }
                    else
                    {
                        for (int k = 0; k < i1->Size(); k++)
                        {
                            hwComplex val1 = i1->z(k);
                            hwComplex val2 = i2->z(k);
                            if (complexLessThan(val1, val2))
                                result->SetElement(k, val1);
                            else
                                result->SetElement(k, val2);
                        }
                    }
                }
                outputs.push_back(result);
            }
        }
    }
    else // one input is a matrix, another is either scalar or complex
    {
        const hwMatrix* mtx;
        hwMatrix* result;
        double dbl;
        hwComplex cplx;
        bool scalar = false;

        if (input1.IsMatrix() || input1.IsString())
        {
            mtx = input1.Matrix();
            if (input2.IsScalar())
            {
                dbl = input2.Scalar();
                scalar = true;
            }
            else if (input2.IsComplex())
                cplx = input2.Complex();
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else
        {
            mtx = input2.Matrix();
            if (input1.IsScalar())
            {
                dbl = input1.Scalar();
                scalar = true;
            }
            else if (input1.IsComplex())
                cplx = input1.Complex();
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }

        result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
        if (mtx->IsReal())
        {
            double val;
            if (scalar)
            {
                for (int k = 0; k < mtx->Size(); ++k)
                {
                    val = (*mtx)(k);
                    result->SetElement(k, min(dbl, val));
                }
            }
            else
            {
                for (int k = 0; k < mtx->Size(); ++k)
                {
                    val = (*mtx)(k);
                    if (complexLessThan(hwComplex(val, 0.0), cplx))
                        result->SetElement(k, val);
                    else
                        result->SetElement(k, cplx);
                }
            }
        }
        else
        {
            hwComplex val;
            if (scalar)
            {
                for (int k = 0; k < mtx->Size(); ++k)
                {
                    val = mtx->z(k);
                    if (complexLessThan(hwComplex(dbl, 0.0), val))
                        result->SetElement(k, dbl);
                    else
                        result->SetElement(k, val);
                }
            }
            else
            {
                for (int k = 0; k < mtx->Size(); ++k)
                {
                    val = mtx->z(k);
                    if (complexLessThan(cplx, val))
                        result->SetElement(k, cplx);
                    else
                        result->SetElement(k, val);
                }
            }
        }
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the max value in the given input [max]
//------------------------------------------------------------------------------
bool oml_max(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();
    int dim = 1;
    Currency input1, input2, input3;

    if (size < 1 || size > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    input1 = inputs[0];
    if (size > 1)
    {
        input2 = inputs[1];
        if (size > 2)
            input3 = inputs[2];
    }

    if (size == 3)
    {
        if (!input2.IsMatrix())
            throw OML_Error(OML_ERR_EMPTYMATRIX, 2, OML_VAR_INPUT);

        if (input2.Matrix()->M() || input2.Matrix()->N())
            throw OML_Error(OML_ERR_EMPTYMATRIX, 2, OML_VAR_INPUT);

        if (!input3.IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);

        dim = static_cast<int>(input3.Scalar());
    }

    if (size == 1)
    {
        if (input1.IsScalar() || input1.IsComplex())
        {
            outputs.push_back(input1);
            outputs.push_back(1.0);
        }
        else if (input1.IsMatrix() || input1.IsString())
        {
            const hwMatrix *mtx = input1.Matrix();
            int index = 0;
            if (mtx->IsVector())
            {
                if (mtx->IsReal())
                {
                    double max = (*mtx)(0);
                    for (int k = 1; k < mtx->Size(); k++)
                    {
                        double temp = (*mtx)(k);
                        if (temp > max)
                        {
                            max = temp;
                            index = k;
                        }
                    }
                    outputs.push_back(max);
                    outputs.push_back(index + 1);
                }
                else
                {
                    hwComplex max = mtx->z(0);
                    for (int k = 1; k < mtx->Size(); k++)
                    {
                        hwComplex temp = mtx->z(k);
                        if (complexGreaterThan(temp, max))
                        {
                            max = temp;
                            index = k;
                        }
                    }
                    outputs.push_back(max);
                    outputs.push_back(index + 1);
                }
            }
            else
            {
                if (mtx->M() == 0)
                {
                    outputs.push_back(EvaluatorInterface::allocateMatrix(0, mtx->N(), hwMatrix::REAL)); // maximums
                    outputs.push_back(EvaluatorInterface::allocateMatrix(0, mtx->N(), hwMatrix::REAL)); // indices
                }
                else
                {
                    hwMatrix *result = EvaluatorInterface::allocateMatrix(1, mtx->N(), mtx->Type());
                    hwMatrix *indices = EvaluatorInterface::allocateMatrix(1, mtx->N(), mtx->Type());

                    if (mtx->IsReal())
                    {
                        double loc_max;
                        for (int i = 0; i < mtx->N(); ++i)
                        {
                            loc_max = (*mtx)(0, i);
                            index = 0;
                            for (int j = 1; j < mtx->M(); ++j)
                            {
                                double temp = (*mtx)(j, i);
                                if (temp > loc_max)
                                {
                                    loc_max = temp;
                                    index = j;
                                }
                            }
                            result->SetElement(i, loc_max);
                            indices->SetElement(i, index + 1);
                        }
                    }
                    else
                    {
                        hwComplex max;
                        for (int i = 0; i < mtx->N(); ++i)
                        {
                            max = mtx->z(0, i);
                            index = 0;
                            for (int j = 1; j < mtx->M(); ++j)
                            {
                                hwComplex temp = mtx->z(j, i);
                                if (complexGreaterThan(temp, max))
                                {
                                    max = temp;
                                    index = j;
                                }
                            }
                            result->SetElement(i, max);
                            indices->SetElement(i, index + 1);
                        }
                    }
                    outputs.push_back(result);
                    outputs.push_back(indices);
                }
            }
        }
        else if (input1.IsNDMatrix())
        {
            return oml_MatrixNUtil3(eval, inputs, outputs, oml_max);
        }
        else
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }
    else if (input1.IsNDMatrix() && size == 3)
    {
        return oml_MatrixNUtil3(eval, inputs, outputs, oml_max, 3);
    }
    else if (input1.IsNDMatrix() || input2.IsNDMatrix())
    {
        return oml_MatrixNUtil2(eval, inputs, outputs, oml_max);
    }
    else if (!((input1.IsMatrix() || input1.IsString()) || (input2.IsMatrix() || input2.IsString())))
    {
        if (input1.IsScalar())
        {
            if (input2.IsScalar())
            {
                outputs.push_back(max(input1.Scalar(), input2.Scalar()));
            }
            else if (input2.IsComplex())
            {
                double di1 = input1.Scalar();
                hwComplex i1 = hwComplex(di1, 0.0);
                hwComplex i2 = input2.Complex();
                if (complexGreaterThan(i1, i2))
                    outputs.push_back(di1);
                else
                    outputs.push_back(i2);
            }
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else if (input1.IsComplex())
        {
            if (input2.IsScalar())
            {
                hwComplex i1 = input1.Complex();
                double di2 = input2.Scalar();
                hwComplex i2 = hwComplex(di2, 0.0);
                if (complexGreaterThan(i2, i1))
                    outputs.push_back(di2);
                else
                    outputs.push_back(i1);
            }
            else if (input2.IsComplex())
            {
                hwComplex i1 = input1.Complex(), i2 = input2.Complex();
                if (complexGreaterThan(i1, i2))
                    outputs.push_back(i1);
                else
                    outputs.push_back(i2);
            }
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else
        {
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
        }
    }
    else if ((input1.IsMatrix() || input1.IsString()) && (input2.IsMatrix() || input2.IsString()))
    {
        const hwMatrix *i1 = input1.Matrix(), *i2 = input2.Matrix();

        if (i1->M() != i2->M() || i1->N() != i2->N()) // matrices have different dimensions
        {
            if (!i2->M() && !i2->N() && size == 3 && dim < 3)
            {
                hwMatrix *result;
                int m = i1->M(), n = i1->N();
                if (dim == 1)
                    result = EvaluatorInterface::allocateMatrix(1, n, i1->Type());
                else
                    result = EvaluatorInterface::allocateMatrix(m, 1, i1->Type());
                hwMatrix* indices = EvaluatorInterface::allocateMatrix(result->M(), result->N(), hwMatrix::REAL);

                if (i1->IsReal())
                {
                    double loc_max;
                    int max_index;
                    for (int i = 0; i < (dim == 1 ? n : m); ++i)
                    {
                        max_index = 0;
                        if (dim == 1)
                            loc_max = (*i1)(0, i);
                        else
                            loc_max = (*i1)(i, 0);

                        for (int j = 1; j < (dim == 1 ? m : n); ++j)
                        {
                            double val = (dim == 1) ? (*i1)(j,i) : (*i1)(i,j);
                            if (val > loc_max)
                            {
                                loc_max = val;
                                max_index = j;
                            }
                        }
                        result->SetElement(i, loc_max);
                        (*indices)(i) = max_index + 1.0;
                    }
                }
                else
                {
                    hwComplex max;
                    int max_index;
                    for (int i = 0; i < (dim == 1 ? n : m); ++i)
                    {
                        max_index = 0;
                        if (dim == 1)
                            max = i1->z(0, i);
                        else
                            max = i1->z(i, 0);

                        for (int j = 1; j < (dim == 1 ? m : n); ++j)
                        {
                            if (dim == 1)
                            {
                                hwComplex val = i1->z(j, i);
                                if (complexGreaterThan(val, max))
                                {
                                    max = val;
                                    max_index = j;
                                }
                            }
                            else
                            {
                                hwComplex val = i1->z(i, j);
                                if (complexGreaterThan(val, max))
                                {
                                    max = val;
                                    max_index = j;
                                }
                            }
                        }
                        result->SetElement(i, max);
                        (*indices)(i) = max_index + 1.0;
                    }
                }
                outputs.push_back(result);
                outputs.push_back(indices);
            }
            else if (!i2->M() && !i2->N() && size == 3 && dim > 2)
            {
                outputs.push_back(input1);
            }
            else
            {
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);
            }
        }
        else // matrices have same dimensions
        {
            if (!i1->Size())
            {
                outputs.push_back(EvaluatorInterface::allocateMatrix(0, 0, hwMatrix::REAL));
            }
            else
            {
                hwMatrix *result = EvaluatorInterface::allocateMatrix(i1->M(), i1->N(), i1->Type());

                if (i1->IsReal())
                {
                    if (i2->IsReal())
                    {
                        for (int k = 0; k < i1->Size(); k++)
                        {
                            result->SetElement(k, max((*i1)(k), (*i2)(k)));
                        }
                    }
                    else
                    {
                        for (int k = 0; k < i1->Size(); k++)
                        {
                            double v1 = (*i1)(k);
                            hwComplex val1 = hwComplex(v1, 0.0);
                            hwComplex val2 = i2->z(k);
                            if (complexGreaterThan(val1, val2))
                                result->SetElement(k, val1);
                            else
                                result->SetElement(k, val2);
                        }
                    }
                }
                else
                {
                    if (i2->IsReal())
                    {
                        for (int k = 0; k < i1->Size(); k++)
                        {
                            hwComplex val1 = i1->z(k);
                            double v2 = (*i2)(k);
                            hwComplex val2 = hwComplex(v2, 0.0);
                            if (complexGreaterThan(val2, val1))
                                result->SetElement(k, v2);
                            else
                                result->SetElement(k, val1);
                        }
                    }
                    else
                    {
                        for (int k = 0; k < i1->Size(); k++)
                        {
                            hwComplex val1 = i1->z(k);
                            hwComplex val2 = i2->z(k);
                            if (complexGreaterThan(val1, val2))
                                result->SetElement(k, val1);
                            else
                                result->SetElement(k, val2);
                        }
                    }
                }
                outputs.push_back(result);
            }
        }
    }
    else // one input is a matrix, another is either scalar or complex
    {
        const hwMatrix* mtx;
        hwMatrix* result;
        double dbl;
        hwComplex cplx;
        bool scalar = false;

        if (input1.IsMatrix() || input1.IsString())
        {
            mtx = input1.Matrix();
            if (input2.IsScalar())
            {
                dbl = input2.Scalar();
                scalar = true;
            }
            else if (input2.IsComplex())
                cplx = input2.Complex();
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else
        {
            mtx = input2.Matrix();
            if (input1.IsScalar())
            {
                dbl = input1.Scalar();
                scalar = true;
            }
            else if (input1.IsComplex())
                cplx = input1.Complex();
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }

        result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
        if (mtx->IsReal())
        {
            double val;
            if (scalar)
            {
                for (int k = 0; k < mtx->Size(); ++k)
                {
                    val = (*mtx)(k);
                    result->SetElement(k, max(dbl, val));
                }
            }
            else
            {
                for (int k = 0; k < mtx->Size(); ++k)
                {
                    val = (*mtx)(k);
                    if (complexGreaterThan(hwComplex(val, 0.0), cplx))
                        result->SetElement(k, val);
                    else
                        result->SetElement(k, cplx);
                }
            }
        }
        else
        {
            hwComplex val;
            if (scalar)
            {
                for (int k = 0; k < mtx->Size(); ++k)
                {
                    val = mtx->z(k);
                    if (complexGreaterThan(hwComplex(dbl, 0.0), val))
                        result->SetElement(k, dbl);
                    else
                        result->SetElement(k, val);
                }
            }
            else
            {
                for (int k = 0; k < mtx->Size(); ++k)
                {
                    val = mtx->z(k);
                    if (complexGreaterThan(cplx, val))
                        result->SetElement(k, cplx);
                    else
                        result->SetElement(k, val);
                }
            }
        }
        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Vector dot product [dot]
//------------------------------------------------------------------------------
bool oml_dot(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    if (size < 2 || size > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];
    const Currency& input2 = inputs[1];

    if (input1.IsScalar())
    {
        if (input2.IsScalar())
            outputs.push_back(input1.Scalar() * input2.Scalar());
        else if (input2.IsComplex())
            outputs.push_back(input2.Complex() * input1.Scalar());
        else if (!input2.IsMatrix())
            throw OML_Error(HW_ERROR_INPUTSCALARMATRIX);
    }
    else if (input1.IsComplex())
    {
        if (input2.IsScalar())
            outputs.push_back(input1.Complex().Conjugate() * input2.Scalar());
        else if (input2.IsComplex())
            outputs.push_back(input1.Complex().Conjugate() * input2.Complex());
        else if (!input2.IsMatrix())
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }
    else if (input1.IsMatrix())
    {
        if (!input2.IsMatrix())
            throw OML_Error(HW_ERROR_MATRIXDIM);

        const hwMatrix *m1 = input1.Matrix(), *m2 = input2.Matrix();
        hwMatrix* result;

        int dim = 0;
        if (size == 3)
        {
            const Currency& input3 = inputs[2];
            if (input3.IsScalar())
                dim = (int) input3.Scalar();
            else if (input3.IsComplex())
                dim = (int) input3.Complex().Real();
            else
            {
                BuiltInFuncsUtils::SetWarning(eval, "Warning: matrix used as matrix dimension; using first element's real component");
                const hwMatrix *mtx = input3.Matrix();

                if (!mtx->Size())
                    throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIMS);

                if (mtx->IsReal())
                    dim = (int) (*mtx)(0);
                else
                    dim = (int) mtx->z(0).Real();
            }
        }
        hwMathStatus stat;

        // if 3 inputs, do it as a matrix in case the dimension would split the vector
        if (size != 3 && m1->IsVector() && m2->IsVector())
        {
            if (m1->IsReal() && m2->IsReal())
            {
                double dot;
                BuiltInFuncsUtils::CheckMathStatus(eval, hwMatrix::Dot(*m1, *m2, dot));
                outputs.push_back(dot);
            }
            else if (!m1->IsReal() && !m2->IsReal())
            {
                hwComplex dot;
                BuiltInFuncsUtils::CheckMathStatus(eval, hwMatrix::Dot(*m1, *m2, dot));
                outputs.push_back(dot);
            }
            else if (m1->IsReal())
            {
                hwMatrix* temp = new hwMatrix;
                temp->PackComplex(*m1, NULL);
                std::vector<Currency> newInputs;
                newInputs.push_back(temp);
                newInputs.push_back(input2);
                return oml_dot(eval, newInputs, outputs);
            }
            else // m2->IsReal()
            {
                hwMatrix* temp = new hwMatrix;
                temp->PackComplex(*m2, NULL);
                std::vector<Currency> newInputs;
                newInputs.push_back(input1);
                newInputs.push_back(temp);
                return oml_dot(eval, newInputs, outputs);
            }

            return true;
        }
        else if (m1->M() == m2->M() && m1->N() == m2->N())
        {
            if (!dim)
                dim = 1;
            else if (dim < 1 || dim > 2)
                throw OML_Error(OML_ERR_UNSUPPORTDIM, 3);

            if (m1->IsReal() && m2->IsReal())
            {
                double dot;
                if (dim == 1)
                {
                    result = EvaluatorInterface::allocateMatrix(1, m1->N(), hwMatrix::REAL);

                    for (int i = 0; i < m1->N(); ++i)
                    {
                        std::unique_ptr<const hwMatrix> vec1(EvaluatorInterface::allocateColumn(m1, i));
                        std::unique_ptr<const hwMatrix> vec2(EvaluatorInterface::allocateColumn(m2, i));
                        BuiltInFuncsUtils::CheckMathStatus(eval, hwMatrix::Dot(*vec1, *vec2, dot));
                        result->SetElement(i, dot);
                    }
                }
                else // dim == 2
                {
                    std::unique_ptr<hwMatrix> vec1(EvaluatorInterface::allocateMatrix());
                    std::unique_ptr<hwMatrix> vec2(EvaluatorInterface::allocateMatrix());
                    result = EvaluatorInterface::allocateMatrix(m1->M(), 1, hwMatrix::REAL);

                    for (int i = 0; i < m1->M(); ++i)
                    {
                        BuiltInFuncsUtils::CheckMathStatus(eval, m1->ReadRow(i, *vec1));
                        BuiltInFuncsUtils::CheckMathStatus(eval, m2->ReadRow(i, *vec2));
                        BuiltInFuncsUtils::CheckMathStatus(eval, hwMatrix::Dot(*vec1, *vec2, dot));
                        result->SetElement(i, dot);
                    }
                }
            }
            else if (!m1->IsReal() && !m2->IsReal())
            {
                hwComplex dot;
                if (dim == 1)
                {
                    result = EvaluatorInterface::allocateMatrix(1, m1->N(), hwMatrix::COMPLEX);

                    for (int i = 0; i < m1->N(); ++i)
                    {
                        std::unique_ptr<const hwMatrix> vec1(EvaluatorInterface::allocateColumn(m1, i));
                        std::unique_ptr<const hwMatrix> vec2(EvaluatorInterface::allocateColumn(m2, i));
                        BuiltInFuncsUtils::CheckMathStatus(eval, hwMatrix::Dot(*vec1, *vec2, dot));
                        result->SetElement(i, dot);
                    }
                }
                else // dim == 2
                {
                    std::unique_ptr<hwMatrix> vec1(EvaluatorInterface::allocateMatrix());
                    std::unique_ptr<hwMatrix> vec2(EvaluatorInterface::allocateMatrix());
                    result = EvaluatorInterface::allocateMatrix(m1->M(), 1, hwMatrix::COMPLEX);

                    for (int i = 0; i < m1->M(); ++i)
                    {
                        BuiltInFuncsUtils::CheckMathStatus(eval, m1->ReadRow(i, *vec1));
                        BuiltInFuncsUtils::CheckMathStatus(eval, m2->ReadRow(i, *vec2));
                        BuiltInFuncsUtils::CheckMathStatus(eval, hwMatrix::Dot(*vec1, *vec2, dot));
                        result->SetElement(i, dot);
                    }
                }
            }
            else if (m1->IsReal())
            {
                hwMatrix* temp = new hwMatrix;
                temp->PackComplex(*m1, NULL);
                std::vector<Currency> newInputs;
                newInputs.push_back(temp);
                newInputs.push_back(input2);
                newInputs.push_back(dim);
                return oml_dot(eval, newInputs, outputs);
            }
            else // m2->IsReal()
            {
                hwMatrix* temp = new hwMatrix;
                temp->PackComplex(*m2, NULL);
                std::vector<Currency> newInputs;
                newInputs.push_back(input1);
                newInputs.push_back(temp);
                newInputs.push_back(dim);
                return oml_dot(eval, newInputs, outputs);
            }
        }
        else
        {
            throw OML_Error(HW_ERROR_MATRIXDIM);
        }

        outputs.push_back(result);
    }
    else if (input1.IsNDMatrix() && input2.IsNDMatrix())
    {
        return oml_MatrixN_VecProd(eval, inputs, outputs, oml_dot, 1);
    }
	else
	{
		throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
	}

    return true;
}
//------------------------------------------------------------------------------
// Kronecker matrix product [kron]
//------------------------------------------------------------------------------
bool oml_kron(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    if (size < 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (size > 2)
    {
        // std::vector<Currency> pair_inputs(inputs.begin(), inputs.begin()+1);
        std::vector<Currency> pair_inputs;
        pair_inputs.push_back(inputs[0]);

        for (size_t i = 1; i < size; ++i)
        {
            pair_inputs.push_back(inputs[i]);

            try
            {
                oml_kron(eval, pair_inputs, outputs);
            }
            catch (OML_Error& err)
            {
                if (err.Arg1() == 2)
                    err.Arg1((int)i+1);

                throw err;
            }

            if (i != size - 1)
            {
                pair_inputs.clear();
                pair_inputs.push_back(outputs[0]);
                outputs.clear();
            }
        }
    }

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar() && !inputs[1].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);

    const hwMatrix* mat1 = inputs[0].ConvertToMatrix();
    const hwMatrix* mat2 = inputs[1].ConvertToMatrix();
    std::unique_ptr<hwMatrix> result(EvaluatorInterface::allocateMatrix());
    
    hwMathStatus status = result->Kronecker(*mat1, *mat2);
    BuiltInFuncsUtils::CheckMathStatus(eval, status);

    outputs.push_back(result.release());

    return true;
}
//------------------------------------------------------------------------------
// Computes the rank of the input matrix [rank]
//------------------------------------------------------------------------------
bool oml_rank(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    if (size != 1 && size != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsLogical() || (size > 1 && inputs[1].IsLogical()))
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    const Currency &input1 = inputs[0];
    Currency input2;
    hwMathStatus stat;

    if (size == 1)
    {
        if (input1.IsScalar() || input1.IsComplex())
        {
            outputs.push_back(Currency(1.0));
            return true;
        }
        else if (input1.IsMatrix())
        {
            const hwMatrix *mtx = input1.Matrix();
            int rank;
            stat = mtx->Rank(rank);
            if (stat.IsOk())
            {
                outputs.push_back(rank);
                return true;
            }
        }
        else
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }
    else // size == 2
    {
        input2 = inputs[1];
        hwMatrix *mtx;

        double tol;
        int rank;
        bool ok = true;

        if (input1.IsScalar() || input1.IsComplex())
        {
            mtx = EvaluatorInterface::allocateMatrix(1, 1, hwMatrix::REAL);

            if (input1.IsScalar())
                mtx->SetElement(0, input1.Scalar());
            else if (input1.IsComplex())
                mtx->SetElement(0, input1.Complex().Mag());

        }
        else if (!input1.IsMatrix())
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
        else
            mtx = EvaluatorInterface::allocateMatrix(input1.Matrix());

        if (input2.IsScalar())
        {
            tol = input2.Scalar();
            stat = mtx->Rank(tol, rank);
            if (!stat.IsOk())
                ok = false;
        }
        else if (input2.IsComplex())
        {
            tol = input2.Complex().Real();
            stat = mtx->Rank(tol, rank);
            if (!stat.IsOk())
                ok = false;
        }
        else if (input2.IsMatrix())
        {
            throw OML_Error(HW_MATH_MSG_NOTIMPLEMENT);
        }
        //	hwMatrix *tolerances = input2.Matrix();
        //
        //	if ((!mtx->Size()) != (!tolerances->Size()))
        //	{
        //		outputs.push_back(Currency(0.0));
        //		return true;
        //	}

        //	if (mtx->IsVector())
        //	{
        //		if (tolerances->IsVector())
        //		{
        //			rank = 0;
        //			for (int i=0; i<tolerances->Size(); ++i)
        //			{
        //				if (tolerances->IsReal())
        //					tol = (*tolerances)(i);
        //				else
        //					tol = tolerances->z(i).Real();

        //				int tempRank = 0;
        //				if (tol < 0 && mtx->Size() == 1)
        //				{
        //					double val;
        //					if (mtx->IsReal())
        //						val = (*mtx)(0);
        //					else
        //						val = mtx->z(0).Real();
        //					if (val < 0)
        //						tempRank = (int) (tol >= val);
        //					else
        //						tempRank = (int) (tol > val);
        //				}
        //				else
        //				{
        //					stat = mtx->Rank(tol, tempRank);
        //					if (!stat.IsOk())
        //					{
        //						ok=false;
        //						break;
        //					}
        //				}
        //				rank += tempRank;
        //			}
        //		}
        //		else
        //		{
        //			hwMatrix *result = EvaluatorInterface::allocateMatrix(1, tolerances->N(), hwMatrix::REAL);
        //          Currency resultCur(result);
        //			for (int c=0; c<tolerances->N(); ++c)
        //			{
        //				rank = 0;
        //				for (int r=0; r<tolerances->M();++r)
        //				{
        //					if (tolerances->IsReal())
        //						tol = (*tolerances)(r);
        //					else
        //						tol = tolerances->z(r).Real();

        //					int tempRank;
        //					if (tol < 0 && mtx->Size() == 1)
        //					{
        //						double val;
        //						if (mtx->IsReal())
        //							val = (*mtx)(0);
        //						else
        //							val = mtx->z(0).Real();
        //						if (val < 0)
        //							tempRank = (int) (tol >= val);
        //						else
        //							tempRank = (int) (tol > val);
        //					}
        //					else
        //					{
        //						stat = mtx->Rank(tol, tempRank);
        //						if (!stat.IsOk())
        //						{
        //							ok=false;
        //							break;
        //						}
        //					}
        //					rank += tempRank;
        //				}
        //				if (!ok)
        //					break;
        //				result->SetElement(c, rank);
        //			}
        //			if (ok)
        //			{
        //				outputs.push_back(resultCur);
        //				return true;
        //			}
        //		}
        //	}
        //	else // !mtx->IsVector()
        //	{
        //		int m=mtx->M(), n=mtx->N();
        //		if (!m || !n)
        //			throw OML_Error("Error: cannot rank an empty matrix");
        //		if (m != tolerances->M() || tolerances->N() != 1)
        //			throw OML_Error("Error: bad matrix tolerance size; should be a vector with the same number of rows are the matrix to rank");

        //		// sum ranks of each row of mtx with the tolerance value at the same row
        //		hwMatrix *row = EvaluatorInterface::allocateMatrix();
        //      Currency rowCur(row);
        //		rank = 0;
        //		for (int i=0; i<m; ++i)
        //		{
        //			stat = mtx->ReadRow(i,*row);
        //			if (tolerances->IsReal())
        //				tol = (*tolerances)(i);
        //			else
        //				tol = tolerances->z(i).Real();
        //			int temp;
        //			stat = row->Rank(tol, temp);
        //			if (!stat.IsOk())
        //			{
        //				ok = false;
        //				break;
        //			}
        //			rank += temp;
        //		}
        //	}

        //}
        //else
        //	throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);

        if (ok)
        {
            outputs.push_back(Currency(rank));
            return true;
        }
    }

    BuiltInFuncsUtils::CheckMathStatus(eval, stat);
    return true;
}
//------------------------------------------------------------------------------
// Singular value decomposition [svd]
//------------------------------------------------------------------------------
bool oml_svd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int flag;
    int size = (int)inputs.size();

    if (size != 1 && size != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    int nargout = getNumOutputs(eval);

    if (nargout < 2)
        flag = 2;
    else
        flag = (int) (size == 2);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_VARIABLE);

    const hwMatrix* mtx = inputs[0].ConvertToMatrix();
    std::unique_ptr<hwMatrix> u, v, s(EvaluatorInterface::allocateMatrix());

    if (nargout > 1)
    {
        u.reset(EvaluatorInterface::allocateMatrix());
        if (nargout > 2)
            v.reset(EvaluatorInterface::allocateMatrix());
    }

    hwMathStatus stat = mtx->SVD(flag, u.get(), *s, v.get());

    BuiltInFuncsUtils::CheckMathStatus(eval, stat);

    if (nargout < 2)
    {
        outputs.push_back(s.release());
        return true;
    }

    int m = mtx->M(), n = mtx->N(), local_min = min(m, n);
    hwMatrix *w;

    if (flag)
        w = EvaluatorInterface::allocateMatrix(local_min, local_min, hwMatrix::REAL);
    else
        w = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);

    w->SetElements(0.0);

    for (int i = 0; i < local_min; i++)
        (*w)(i, i) = (*s)(i);

    outputs.push_back(u.release());
    outputs.push_back(w);

    if (nargout > 2)
        outputs.push_back(v.release());

    return true;
}
//------------------------------------------------------------------------------
// Gets schur decomposition of input [schur]
//------------------------------------------------------------------------------
bool oml_schur(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];
    int nargout = getNumOutputs(eval);

    if (input.IsScalar() || input.IsComplex())
    {
        if (nargout > 1)
            outputs.push_back(1.0);
        outputs.push_back(input);
        return true;
    }
    else if (input.IsMatrix())
    {
        const hwMatrix *mtx;
        hwMatrix *U, *result;
        mtx = input.Matrix();

        U = EvaluatorInterface::allocateMatrix();
        result = EvaluatorInterface::allocateMatrix();

        bool real = mtx->IsRealData();

        if (nargin > 1)
        {
            std::string opt = readOption(eval, inputs[1]);
            if (opt == "complex")
                real = false;
            else if (opt == "real")
                real = true;
            else
                throw OML_Error(HW_ERROR_INVALIDOPTION(opt));
        }

        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Schur(real, *U, *result));

        if (nargout > 1)
            outputs.push_back(U);
        outputs.push_back(result);
        return true;
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }
}
//------------------------------------------------------------------------------
// Returns qr decomposition of input matrix [qr]
//------------------------------------------------------------------------------
bool oml_qr(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() < 1 || inputs.size() > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_VARIABLE);

    const hwMatrix* A = inputs[0].ConvertToMatrix();
    std::unique_ptr<hwMatrix> Q(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> R(EvaluatorInterface::allocateMatrix());
    hwMathStatus status;

    if (inputs.size() == 1)
    {
        int nCol = A->M() - A->N();

        if (nCol > 0)
        {
            hwMatrix AA;

            status = AA.InsertColumns(*A, A->N(), nCol);
            if (!status.IsOk())
            {
                status.ResetArgs();
                BuiltInFuncsUtils::CheckMathStatus(eval, status);
            }

            status = AA.QR(*Q, *R);
            if (!status.IsOk())
            {
                status.ResetArgs();
                BuiltInFuncsUtils::CheckMathStatus(eval, status);
            }

            status = R->DeleteColumns(A->N(), nCol);
            if (!status.IsOk())
            {
                status.ResetArgs();
                BuiltInFuncsUtils::CheckMathStatus(eval, status);
            }
        }
        else
        {
            status = A->QR(*Q, *R);
        }
    }
    else // (inputs.size() == 2)
    {
        if (!inputs[1].IsScalar())
            throw OML_Error(OML_ERR_SCALAR, 2);

        if (inputs[1].Scalar() != 0)
            throw OML_Error(HW_ERROR_INVINPVAL);

        status = A->QR(*Q, *R);
    }

    BuiltInFuncsUtils::CheckMathStatus(eval, status);
    outputs.push_back(Q.release());
    outputs.push_back(R.release());

    return true;
}
//------------------------------------------------------------------------------
// Gets lu decomposition of input matrix [lu]
//------------------------------------------------------------------------------
bool oml_lu(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();
    if (size != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_VARIABLE);

    const hwMatrix* mtx = inputs[0].ConvertToMatrix();
    std::unique_ptr<hwMatrix> L(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> U(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrixI> Piv(new hwMatrixI);

    BuiltInFuncsUtils::CheckMathStatus(eval, mtx->LU(*Piv, *L, *U));

    int m = mtx->M();
    hwMatrix* P = EvaluatorInterface::allocateMatrix(m, m, hwMatrix::REAL);
    Currency pcur(P);
    int nargout = getNumOutputs(eval);

    if (nargout < 2)
    {
        if (L->M() < U->M() || L->N() < U->N())
        {
            // copy L to U
            if (checkMakeComplex(eval, L.get(), U.get()))
            {
                for (int j = 0; j < L->N()-1; j++)
                {
                    for (int i = j+1; i < L->M(); i++)
                    {
                        U->z(i, j) = L->z(i, j);
                    }
                }
            }
            else
            {
                for (int j = 0; j < L->N()-1; j++)
                {
                    for (int i = j+1; i < L->M(); i++)
                    {
                        (*U)(i, j) = (*L)(i, j);
                    }
                }
            }

            outputs.push_back(U.release());
        }
        else
        {
            // copy U to L
            if (checkMakeComplex(eval, L.get(), U.get()))
            {
                for (int j = 0; j < U->N(); j++)
                {
                    for (int i = 0; i <= j; i++)
                    {
                        L->z(i, j) = U->z(i, j);
                    }
                }
            }
            else
            {
                for (int j = 0; j < U->N(); j++)
                {
                    for (int i = 0; i <= j; i++)
                    {
                        (*L)(i, j) = (*U)(i, j);
                    }
                }
            }

            outputs.push_back(L.release());        }

        return true;
    }

    if (P->Size() != 0)
    {
        P->SetElements(0.0);
        for (int i = 0; i < m; i++)
        {
            (*P)((*Piv)(i), i) = 1.0;
        }
    }

    if (nargout == 2)
        *L = *P * *L;
    else
        BuiltInFuncsUtils::CheckMathStatus(eval, P->Transpose());

    outputs.push_back(L.release());
    outputs.push_back(U.release());
    outputs.push_back(pcur);

    return true;
}
//------------------------------------------------------------------------------
// Returns a vector or matrix of differences [diff]
//------------------------------------------------------------------------------
bool oml_diff(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int cycles = 1;
    int dim = 1;
    size_t size = inputs.size();

    if (size < 1 || size > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];
    if (size > 1)
    {
        const Currency& input2 = inputs[1];

        // set up number of cycles -- default is 1
        if (input2.IsScalar())
            cycles = (int) input2.Scalar();
        else if (input2.IsComplex())
            cycles = (int) input2.Complex().Real();
        else
            throw OML_Error(OML_ERR_NATURALNUM, 2);

        if (cycles < 0)
            throw OML_Error(OML_ERR_NATURALNUM, 2);

        // set up dimension -- default is first dimension above 0
        if (size > 2)
        {
            const Currency& input3 = inputs[2];
            if (input3.IsScalar())
            {
                dim = (int) input3.Scalar();
                if (dim < 0)
                {
                    throw OML_Error(OML_ERR_UNSUPPORTDIM, 3, OML_VAR_DIM);
                }
            }
            else
            {
                throw OML_Error(OML_ERR_NATURALNUM, 3, OML_VAR_DIM);
            }
        }
    }

    if (input1.IsScalar() || input1.IsComplex())
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
    }
    else
    {
        if (input1.IsMatrix())
        {
            int m = input1.Matrix()->M();
            int n = input1.Matrix()->N();

            if (dim > 2)
            {
                std::vector<int> dims;
                dims.push_back(m);
                dims.push_back(n);

                for (int i = 2; i < dim-1; ++i)
                    dims.push_back(1);

                dims.push_back(0);
                hwMatrixN* result = EvaluatorInterface::allocateMatrixN();
                result->Dimension(dims, hwMatrixN::REAL);
                outputs.push_back(result);
                return true;
            }

            hwMatrix* result = EvaluatorInterface::allocateMatrix(input1.Matrix());
            Currency resultCur(result);

            bool real = result->IsReal();
            bool recalcDim = (size != 3);
            // loop through, diff across specified dimensions, then delete last row/column

            for (int c = 0; c < cycles; c++)
            {
                //if (!m || !n)
                    //break;
    
                if (recalcDim)
                {
                    if (m == 1 || n == 1)
                        dim = (m != 1 ? 1 : 2);
                    else
                        dim = (n != 1 ? 1 : 2);
                }

                if (m && n)
                {
                    for (int i = 0; i < (dim == 1 ? n : m); ++i)
                    {
                        for (int j = 1; j < (dim == 1 ? m : n); ++j)
                        {
                            // find difference between corresponding elements
                            if (real)
                            {
                                if (dim == 1)
                                {
                                    (*result)(j - 1, i) = (*result)(j, i) - (*result)(j - 1, i);
                                }
                                else
                                { // dim == 2
                                    (*result)(i, j - 1) = (*result)(i, j) - (*result)(i, j - 1);

                                }
                            }
                            else
                            {
                                if (dim == 1)
                                {
                                    result->z(j - 1, i) = result->z(j, i) - result->z(j - 1, i);
                                }
                                else
                                { // dim == 2
                                    result->z(i, j - 1) = result->z(i, j) - result->z(i, j - 1);
                                }
                            }
                        }
                    }
                }

                if (dim == 2)
                    --n;
                else
                    --m;
            }

            safeResizeMatrix(eval, max(0,m), max(0,n), *result, false);
            outputs.push_back(resultCur);
        }
        else if (input1.IsNDMatrix())
        {
            return oml_MatrixN_diff(eval, inputs, outputs, oml_diff);
        }
        else
        {
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns an identity matrix [eye]
//------------------------------------------------------------------------------
bool oml_eye(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    if (!size)
    {
        outputs.push_back(1.0);
        return true;
    }

    if (size != 1 && size != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    int m, n;
    bool valSet = false; // whether the second dimension was also set

    getDimensionsFromInput(inputs, &m, &n);

    hwMatrix *result = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);

    if (result->Size() != m * n)
        throw OML_Error(HW_ERROR_NOTCREATEMATRIX);

    result->Identity();
    outputs.push_back(result);

    return true;
}
//------------------------------------------------------------------------------
// Returns memory usage [memoryuse]
//------------------------------------------------------------------------------
#ifdef OS_WIN
#include "psapi.h"    // MS SDK header, PPROCESS_MEMORY_COUNTERS
typedef BOOL (WINAPI *PFNGetPSMemInfo)(HANDLE, PPROCESS_MEMORY_COUNTERS, DWORD);
bool oml_memoryuse(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    PFNGetPSMemInfo s_pfnGetProcessMemoryInfo;
    HMODULE s_hModPSAPI = ::LoadLibraryA("psapi.dll");

    if (s_hModPSAPI)
    {
        s_pfnGetProcessMemoryInfo = (PFNGetPSMemInfo)GetProcAddress(
            s_hModPSAPI, TEXT("GetProcessMemoryInfo"));
    }

    std::string memstat;
    DWORD processID = GetCurrentProcessId();
    PROCESS_MEMORY_COUNTERS pmc;
    HANDLE hProcess = OpenProcess(PROCESS_QUERY_INFORMATION |
        PROCESS_VM_READ,
        FALSE, processID);

    if (s_pfnGetProcessMemoryInfo &&
        (*s_pfnGetProcessMemoryInfo)(hProcess, &pmc, sizeof(pmc)) )
    {
        char temp[64];
        sprintf(temp, "\tProcessPageFaultCount: %d\n", pmc.PageFaultCount);
        memstat += temp;
        sprintf(temp, "\tProcessPeakWorkingSetSize: %d\n", pmc.PeakWorkingSetSize);
        memstat += temp;
        sprintf(temp, "\tProcessWorkingSetSize: %d\n", pmc.WorkingSetSize);
        memstat += temp;
        sprintf(temp, "\tProcessQuotaPeakPagedPoolUsage: %d\n", pmc.QuotaPeakPagedPoolUsage);
        memstat += temp;
        sprintf(temp, "\tProcessQuotaPagedPoolUsage: %d\n", pmc.QuotaPagedPoolUsage);
        memstat += temp;
        sprintf(temp, "\tProcessQuotaPeakNonPagedPoolUsage: %d\n", pmc.QuotaPeakNonPagedPoolUsage);
        memstat += temp;
        sprintf(temp, "\tProcessQuotaNonPagedPoolUsage: %d\n", pmc.QuotaNonPagedPoolUsage);
        memstat += temp;
        sprintf(temp, "\tProcessPagefileUsage: %d\n", pmc.PagefileUsage);
        memstat += temp;
        sprintf(temp, "\tProcessPeakPagefileUsage: %d\n", pmc.PeakPagefileUsage);
        memstat += temp;
    }
    CloseHandle(hProcess);

    outputs.push_back(memstat);
    return true;
}
#else
bool oml_memoryuse(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return true;
}
#endif
//------------------------------------------------------------------------------
// Cholesky decomposition [chol]
//------------------------------------------------------------------------------
bool oml_chol(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    Currency input1, input2, input3;
    bool upper = true;
    size_t size = inputs.size();

    if (size < 1 || size > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    input1 = inputs[0];
    if (size > 1)
    {
        input2 = inputs[1];
        if (input2.IsString())
        {
            std::string str = readOption(eval, input2);
            if (str == "lower")
            {
                upper = false;
            }
            else if (str != "upper")
            {
                throw OML_Error(OML_ERR_TRIANGMATTYPE, 2, OML_VAR_STRING);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
        }
    }

    hwMatrix *mtx;

    if (input1.IsScalar())
    {
        mtx = EvaluatorInterface::allocateMatrix(1, 1, hwMatrix::REAL);
        mtx->SetElement(0, input1.Scalar());
    }
    else if (input1.IsMatrix())
    {
        mtx = EvaluatorInterface::allocateMatrix((input1.Matrix()));
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, 1, OML_VAR_TYPE);
    }

    hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx);
    hwMathStatus stat = mtx->Csky(*result, upper);

    int index = 0;
    int nargout = getNumOutputs(eval);

    if (!stat.IsOk())
    {
        hwMathMsgCode c = stat.GetMsgCode();

        if (c == HW_MATH_ERR_MTXNOTSPD && nargout > 1)
            index = result->M() + 1;
        else
            BuiltInFuncsUtils::CheckMathStatus(eval, stat);
    }

    outputs.push_back(result);

    if (nargout > 1)
        outputs.push_back(index);

    return true;
}
//------------------------------------------------------------------------------
// Returns the eigen values/vectors of the input(s) [eig]
//------------------------------------------------------------------------------
bool oml_eig(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();
    if (size < 1 || size > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsLogical() || (size > 1 && inputs[1].IsLogical()))
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_VARIABLE);

    const hwMatrix* i1 = inputs[0].ConvertToMatrix();
    const hwMatrix* i2;

    if (size > 1)
    {
        if (!inputs[1].IsMatrix() && !inputs[1].IsScalar() && !inputs[1].IsComplex() && !inputs[1].IsString())
            throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);

        i2 = inputs[1].ConvertToMatrix();
    }

    std::unique_ptr<hwMatrix> D(EvaluatorInterface::allocateMatrix());
    hwMatrix *V = EvaluatorInterface::allocateMatrix();
    Currency vcur(V);

    hwMathStatus stat;

    if (size == 2)
    {
        if (inputs[1].IsString())
        {
            if (inputs[1].StringVal() == "nobalance")
            {
                hwMathStatus status;

                if (i1->IsHermitian())
                    status = i1->EigenSH(V, *D);
                else
                    status = i1->Eigen(false, V, *D);

                if (!status.IsOk())
                    throw OML_Error(status.GetMessage());
            }
            else if (inputs[1].StringVal() == "balance")
            {
                hwMathStatus status;

                if (i1->IsHermitian())
                    status = i1->EigenSH(V, *D);
                else
                    status = i1->Eigen(true, V, *D);

                if (!status.IsOk())
                    throw OML_Error(status.GetMessage());
            }
            else
            {
                throw OML_Error(HW_ERROR_INVBALFLAG);
            }
        }
        else if (sameSize(i1, i2))
            stat = i1->Eigen(*i1, *i2, *V, *D);
        else
            throw OML_Error(HW_ERROR_INPMUSTSAMESIZE);
    }
    else
    {
        hwMathStatus status;

        if (i1->IsHermitian())
            status = i1->EigenSH(V, *D);
        else
            status = i1->Eigen(true, V, *D);

        if (!status.IsOk())
            throw OML_Error(status.GetMessage());
    }

    BuiltInFuncsUtils::CheckMathStatus(eval, stat);
    if (getNumOutputs(eval) <= 1)
    {
        outputs.push_back(D.release());
    }
    else
    {
        outputs.push_back(vcur);
        outputs.push_back(vectorToDiag(D.get()));
    }

    return true;
}
//------------------------------------------------------------------------------
// Creates diagonal matrix or returns the diagonal elements of input [diag]
//------------------------------------------------------------------------------
bool oml_diag(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    Currency input1, input2, input3;

    size_t size = inputs.size();

    if (size < 1 || size > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    input1 = inputs[0];
    if (size > 1)
    {
        input2 = inputs[1];
        if (size > 2)
            input3 = inputs[2];
    }

    // First get/convert data into proper types
    // Then perform diag

    if (!input1.IsMatrix() && !input1.IsScalar() && !input1.IsComplex() && !input1.IsString())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);

    const hwMatrix* i1 = input1.ConvertToMatrix();
    hwMatrix *result = EvaluatorInterface::allocateMatrix();
    Currency resultCur(result);
    int i2, i3;

    if (size == 1)
        i2 = 0;
    else
    {
        if (!input2.IsScalar())
            throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_VALUE);

        if (!IsInteger(input2.Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_VALUE);

        i2 = (int) input2.Scalar();

        if (size == 3)
        {
            if (!i1->IsVector())
                throw OML_Error(HW_ERROR_3INP1STVEC);

            if (!input3.IsScalar())
                throw OML_Error(OML_ERR_NATURALNUM, 3, OML_VAR_DIM);

            if (!IsInteger(input3.Scalar()).IsOk())
                throw OML_Error(OML_ERR_NATURALNUM, 3, OML_VAR_DIM);

            i3 = (int) input3.Scalar();
        }
    }

    // now calculate/create diagonal

    bool ok = true;
    hwMathStatus stat;
    if (size < 3)
    {
        if (i1->Size())
        {
            stat = result->Diag(*i1, i2);
            ok = stat.IsOk();
        }
    }
    else // size == 3
    {
        if (i2 < 0)
            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);

        if (i3 < 0)
            throw OML_Error(OML_ERR_NATURALNUM, 3, OML_VAR_DIM);

        stat = result->Diag(*i1, 0);
        ok = stat.IsOk();

        if (ok)
        {
            // adjust size to be i2 x i3 with new elements as 0
            safeResizeMatrix(eval, i2, i3, *result, true);
        }
    }

    if (!ok)
        BuiltInFuncsUtils::CheckMathStatus(eval, stat);

    if (input1.IsString())
        resultCur.SetMask(Currency::MASK_STRING);
    outputs.push_back(resultCur);
    return true;
}
//------------------------------------------------------------------------------
// Returns transpose of input [transpose]
//------------------------------------------------------------------------------
bool oml_transpose(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    if (input1.IsScalar() || input1.IsComplex() || input1.IsMatrix() || input1.IsCellArray() || input1.IsStruct())
        outputs.push_back(transpose(eval, input1));
    else
        throw OML_Error(HW_ERROR_INVINPTYPE);

    return true;
}
//------------------------------------------------------------------------------
// Returns convolution of input vectors [conv]
//------------------------------------------------------------------------------
bool oml_conv(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex() && !inputs[0].IsString())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar() && !inputs[1].IsComplex() && !inputs[1].IsString())
        throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);

    const hwMatrix* m1 = inputs[0].ConvertToMatrix();
    const hwMatrix* m2 = inputs[1].ConvertToMatrix();
    hwMatrix *result = EvaluatorInterface::allocateMatrix();
    Currency out(result);

    if (!(m1->IsVector() && m2->IsVector()))
        throw OML_Error(HW_ERROR_MATINPMUSTVEC);

    BuiltInFuncsUtils::CheckMathStatus(eval, result->ConvLin(*m2, *m1));
    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// Computes the trace of the input matrix [trace]
//------------------------------------------------------------------------------
bool oml_trace(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsScalar())
    {
        outputs.push_back(input.Scalar());
    }
    else if (input.IsComplex())
    {
        outputs.push_back(input.Complex());
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();

        if (!mtx->IsSquare())
        {
            throw OML_Error(HW_MATH_MSG_MTXNOTSQUARE);
        }

        if (mtx->IsReal())
        {
			double trace = 0.0;

			for (int j=0; j<mtx->M(); j++)
				trace += (*mtx)(j, j);

			outputs.push_back(trace);
        }
        else
        {
			hwComplex trace(0.0, 0.0);

			for (int j=0; j<mtx->M(); j++)
				trace += mtx->z(j, j);

			outputs.push_back(trace);
        }
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes matix determinant [det]
//------------------------------------------------------------------------------
bool oml_det(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsScalar())
    {
        outputs.push_back(Currency(input.Scalar()));
        outputs.push_back(Currency(0.0));
    }
    else if (input.IsComplex())
    {
        outputs.push_back(Currency(input.Complex()));
        outputs.push_back(Currency(0.0));
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMathStatus stat;
        bool ok;

		if (!mtx)
		{
			outputs.push_back(1.0);
			return true;
		}

        if (!mtx->IsSquare())
            throw OML_Error(HW_MATH_MSG_MTXNOTSQUARE);

        if (mtx->IsReal())
        {
            double det;
            stat = mtx->Determinant(det);
            ok = stat.IsOk();
            if (ok)
            {
                outputs.push_back(Currency(det));
                return oml_rcond(eval, inputs, outputs);
            }
        }
        else
        {
            hwComplex det;
            stat = mtx->Determinant(det);
            ok = stat.IsOk();
            if (ok)
            {
                outputs.push_back(det);
                return oml_rcond(eval, inputs, outputs);
            }
        }

        if (!ok)
            BuiltInFuncsUtils::CheckMathStatus(eval, stat);
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Compute the 1-norm estimate of the reciprocal condition number [rcond]
//------------------------------------------------------------------------------
bool oml_rcond(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsScalar())
        outputs.push_back(Currency(input.Scalar()));
    else if (input.IsComplex())
        outputs.push_back(Currency(input.Complex()));
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();

        if (!mtx->IsSquare())
        {
            throw OML_Error(HW_MATH_MSG_MTXNOTSQUARE);
        }

        double rcond = 1.0;
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->RCond(rcond));
        outputs.push_back(Currency(rcond));
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns x modulo y where x and y are the inputs [mod]
//------------------------------------------------------------------------------
bool oml_mod(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    const Currency &input2 = inputs[1];

    if (input1.IsLogical() || input2.IsLogical())
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    if (input1.IsComplex())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

    if (input2.IsComplex())
        throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

    if (input1.IsScalar())
    {
        double dividend = input1.Scalar();

        if (input2.IsScalar())
        {
            double modulo = mod(dividend, input2.Scalar());
            outputs.push_back(Currency(modulo));
        }
        else if (input2.IsMatrix() && !input2.IsString())
        {
            const hwMatrix* mtx = input2.Matrix();

            if (!mtx->IsRealData())
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

            //hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

            for (int k = 0; k < mtx->Size(); k++)
            {
                double divisor = realval(mtx, k);
                result->SetElement(k, mod(dividend, divisor));
            }

            outputs.push_back(result);
        }
        else if (inputs[1].IsNDMatrix())
        {
            return oml_MatrixNUtil2(eval, inputs, outputs, oml_mod);
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }
    }
    else if (input1.IsMatrix() && !input1.IsString())
    {
        const hwMatrix* mtx1 = input1.Matrix();

        if (!mtx1->IsRealData())
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

        //hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx1->M(), mtx1->N(), mtx1->Type());
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx1->M(), mtx1->N(), hwMatrix::REAL);

        if (input2.IsScalar())
        {
            double divisor = input2.Scalar();

            for (int k = 0; k < mtx1->Size(); k++)
            {
                double dividend = realval(mtx1, k);
                result->SetElement(k, mod(dividend, divisor));
            }
        }
        else if (input2.IsMatrix() && !input2.IsString())
        {
            const hwMatrix* mtx2 = input2.Matrix();

            if (!mtx2->IsRealData())
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

            if (!sameSize(mtx1, mtx2))
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

            for (int k = 0; k < mtx1->Size(); k++)
            {
                double divisor  = realval(mtx1, k);
                double dividend = realval(mtx2, k);
                result->SetElement(k, mod(divisor, dividend));
            }
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARMATRIX, 2, OML_VAR_TYPE);
        }

        outputs.push_back(result);
    }
    else if (input1.IsNDMatrix() || input2.IsNDMatrix())
    {
        return oml_MatrixNUtil2(eval, inputs, outputs, oml_mod);
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARMATRIX, 1, OML_VAR_TYPE);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns real part of input [real]
//------------------------------------------------------------------------------
bool oml_real(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsScalar())
    {
        outputs.push_back(Currency(input.Scalar()));
    }
    else if (input.IsComplex())
    {
        hwComplex cplx = input.Complex();
        outputs.push_back(Currency(cplx.Real()));
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result;
        if (mtx->IsReal())
        {
            result = EvaluatorInterface::allocateMatrix(mtx);
        }
        else
        {
            result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

            for (int k = 0; k < mtx->Size(); k++)
            {
                hwComplex cplx = mtx->z(k);
                result->SetElement(k, cplx.Real());
            }
        }

        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_real);
    }
    else if (input.IsString())
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix(input.Matrix()));
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns imaginary part of input [imag]
//------------------------------------------------------------------------------
bool oml_imag(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsScalar())
    {
        outputs.push_back(Currency(0.0));
    }
    else if (input.IsComplex())
    {
        hwComplex cplx = input.Complex();
        outputs.push_back(cplx.Imag());
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

        if (mtx->IsReal())
        {
            for (int k = 0; k < mtx->Size(); k++)
                result->SetElement(k, 0.0);
        }
        else
        {
            for (int k = 0; k < mtx->Size(); k++)
            {
                hwComplex cplx = mtx->z(k);
                result->SetElement(k, cplx.Imag());
            }
        }

        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_imag);
    }
    else if (input.IsString())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), 0.0);
        outputs.push_back(result);
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING);
    }

    return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
void rotHelper(hwMatrix*& mat, bool ccw)
{
    // rewrite to perform out of place
	hwMatrix* temp = mat;
	mat = new hwMatrix(temp->N(), temp->M(), temp->Type());

	if (!ccw)
	{
		if (mat->IsReal())
		{
			for (int j=0; j<mat->M(); j++)
			{
				for (int k=0; k<mat->N(); k++)
					(*mat)(j,k) = (*temp)(mat->N()-k-1, j);
			}
		}
		else
		{
			for (int j=0; j<mat->M(); j++)
			{
				for (int k=0; k<mat->N(); k++)
					mat->z(j,k) = temp->z(mat->N()-k-1, j);
			}
		}
	}
	else
	{
		if (mat->IsReal())
		{
			for (int j=0; j<mat->M(); j++)
			{
				for (int k=0; k<mat->N(); k++)
					(*mat)(mat->M()-j-1,mat->N()-k-1) = (*temp)(mat->N()-k-1, j);
			}
		}
		else
		{
			for (int j=0; j<mat->M(); j++)
			{
				for (int k=0; k<mat->N(); k++)
					mat->z(mat->M()-j-1,mat->N()-k-1) = temp->z(mat->N()-k-1, j);
			}
		}
	}

    delete temp;
}
//------------------------------------------------------------------------------
// Rotates the input matrix in 90 degree succession(s) [rot90]
//------------------------------------------------------------------------------
bool oml_rot90(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (!inputs.size())
		throw OML_Error(OML_ERR_NUMARGIN);

	if (inputs[0].IsScalar())
	{
		outputs.push_back(inputs[0].Scalar());
		return true;
	}
	else if (inputs[0].IsComplex())
	{
		outputs.push_back(inputs[0].Complex());
		return true;
	}
	else if (inputs[0].IsMatrix())
	{
		int  num_rots = 1;
		bool do_ccw   = true;

		if (inputs.size() > 1)
		{
			int input = (int)inputs[1].Scalar();

			if (input < 0)
				do_ccw = false;
		
			if (fabs((double)input))
				num_rots = (int)fabs((double)input); 

			num_rots = num_rots % 4;
		}

		// have to make a copy because of const-ness
		hwMatrix*    result   = EvaluatorInterface::allocateMatrix(inputs[0].Matrix());

		for (int j=0; j<num_rots; j++)
			rotHelper(result, do_ccw);

		outputs.push_back(result);
		return true;
	}
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, 1, OML_VAR_DATA);
    }

	return false;
}
//------------------------------------------------------------------------------
// Computes powers of 2 [pow2]
//------------------------------------------------------------------------------
bool oml_pow2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    Currency coef, exp;

    if (size == 1)
        exp = inputs[0];
    else if (size == 2)
    {
        coef = inputs[0];
        exp = inputs[1];
    }
    else
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    if (size == 2)
    {
        if (!(coef.IsScalar() || coef.IsComplex()))
        {
            if (coef.IsMatrix())
            {
                if (exp.IsMatrix())
                {
                    const hwMatrix *m1 = coef.Matrix(), *m2 = exp.Matrix();
                    if (!sameSize(m1, m2))
                    {
                        throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);
                    }
                }
            }
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
    }

    Currency calcExp;
    if (exp.IsScalar())
        calcExp = Currency(pow(2, exp.Scalar()));
    else if (exp.IsComplex())
        calcExp = Currency(hwComplex::pow(2.0, exp.Complex()));
    else if (exp.IsMatrix())
    {
        hwMatrix *mtx = EvaluatorInterface::allocateMatrix();
        hwMathStatus stat = mtx->PowerByElems(2.0, *(exp.Matrix()));
        BuiltInFuncsUtils::CheckMathStatus(eval, stat);
        calcExp = Currency(mtx);
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    if (size == 2)
        calcExp = curMultElems(coef, calcExp);

    outputs.push_back(calcExp);
    return true;
}
//------------------------------------------------------------------------------
// Returns sign of the value passed in [sign]
//------------------------------------------------------------------------------
bool oml_sign(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, signum, signum))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_sign);
    }

    return true;
}
//------------------------------------------------------------------------------
// For each element of the input x, returns e ^ x [exp]
//------------------------------------------------------------------------------
bool oml_exp(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, exp, hwComplex::exp))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_exp);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the complex conjugate [conj]
//------------------------------------------------------------------------------
bool oml_conj(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsScalar())
    {
        outputs.push_back(Currency(input.Scalar()));
    }
    else if (input.IsComplex())
    {
        hwComplex cplx = input.Complex();
        outputs.push_back(cplx.Conjugate());
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

        if (mtx->IsReal())
        {
            for (int k = 0; k < mtx->Size(); k++)
            {
                result->SetElement(k, (*mtx)(k));
            }
        }
        else
        {
            for (int k = 0; k < mtx->Size(); k++)
            {
                hwComplex cplx = mtx->z(k);
                result->SetElement(k, cplx.Conjugate());
            }
        }

        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_conj);
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the product of elements along given dimension if specified [prod]
//------------------------------------------------------------------------------
bool oml_prod(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int dim = 0;
    size_t size = inputs.size();
    if (size != 1 && size != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    if (size == 2)
    {
		const Currency &temp = inputs[1];
		if (temp.IsScalar() || temp.IsComplex() || (temp.IsMatrix() && !temp.IsString()))
		{
            dim = (int) temp.Scalar();

            if (!IsInteger(temp.Scalar()).IsOk())
                throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_DIM);

            if (dim < 0)
                throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);
		}
		else
		{
			throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
		}
    }

    const Currency &input1 = inputs[0];

    if (input1.IsScalar())
        outputs.push_back(input1.Scalar());
    else if (input1.IsComplex())
        outputs.push_back(input1.Complex());
    else if (input1.IsMatrix() && !input1.IsString())
    {
        const hwMatrix* mtx = input1.Matrix();
        hwMatrix* result;
        if (!mtx->Size())
            outputs.push_back(Currency(1.0));
        else if (dim > 2)
        {
            result = EvaluatorInterface::allocateMatrix(mtx);
            outputs.push_back(result);
        }
        else
        {
            int m = mtx->M(), n = mtx->N();
            if (!dim)
                dim = (m > 1 ? 1 : 2);

            if (dim == 1)
                result = EvaluatorInterface::allocateMatrix(1, n, mtx->Type());
            else if (dim == 2)
                result = EvaluatorInterface::allocateMatrix(m, 1, mtx->Type());

            if (mtx->IsReal())
            {
                double prod;
                for (int i = 0; i < (dim == 1 ? n : m); ++i)
                {
                    prod = 1.0;
                    for (int j = 0; j < (dim == 1 ? m : n); ++j)
                    {
                        if (dim == 1)
                            prod *= (*mtx)(j, i);
                        else
                            prod *= (*mtx)(i, j);
                    }
                    result->SetElement(i, prod);
                }
            }
            else
            {
                hwComplex prod;
                for (int i = 0; i < (dim == 1 ? n : m); ++i)
                {
                    prod = hwComplex(1, 0);
                    for (int j = 0; j < (dim == 1 ? m : n); ++j)
                    {
                        if (dim == 1)
                            prod *= mtx->z(j, i);
                        else
                            prod *= mtx->z(i, j);
                    }
                    result->SetElement(i, prod);
                }
            }
            outputs.push_back(result);
        }
    }
    else if (input1.IsNDMatrix())
    {
        if (size == 1)
        {
            oml_MatrixNUtil3(eval, inputs, outputs, oml_prod);
        }
        else
        {
            oml_MatrixNUtil3(eval, inputs, outputs, oml_prod, 2);
        }
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Sums elements of input, across a dimension if specified [sum]
//------------------------------------------------------------------------------
bool oml_sum(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int dim = 0;
    size_t size = inputs.size();
    if (size != 1 && size != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    if (size == 2)
    {
        if (!inputs[1].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DIM);

        dim = static_cast<int>(inputs[1].Scalar());
    }

    const Currency &input1 = inputs[0];

    if (input1.IsScalar())
        outputs.push_back(input1.Scalar());
    else if (input1.IsComplex())
        outputs.push_back(input1.Complex());
    else if (input1.IsMatrix() || input1.IsString())
    {
        const hwMatrix* mtx = input1.Matrix();
        hwMatrix* result;
        if (!mtx->Size())
            outputs.push_back(Currency(0.0));
        else if (dim > 2)
        {
            result = EvaluatorInterface::allocateMatrix(mtx);
            outputs.push_back(result);
        }
        else
        {
            int m = mtx->M(), n = mtx->N();

            if (!dim)
                dim = (m > 1 ? 1 : 2);

            if (dim == 1)
                result = EvaluatorInterface::allocateMatrix(1, n, mtx->Type());
            else if (dim == 2)
                result = EvaluatorInterface::allocateMatrix(m, 1, mtx->Type());

            if (mtx->IsReal())
            {
                double sum;
                for (int i = 0; i < (dim == 1 ? n : m); ++i)
                {
                    sum = 0.0;
                    for (int j = 0; j < (dim == 1 ? m : n); ++j)
                    {
                        if (dim == 1)
                            sum += (*mtx)(j, i);
                        else
                            sum += (*mtx)(i, j);
                    }
                    result->SetElement(i, sum);
                }
            }
            else
            {
                hwComplex sum;
                for (int i = 0; i < (dim == 1 ? n : m); ++i)
                {
                    sum = hwComplex(0, 0);
                    for (int j = 0; j < (dim == 1 ? m : n); ++j)
                    {
                        if (dim == 1)
                            sum += mtx->z(j, i);
                        else
                            sum += mtx->z(i, j);
                    }
                    result->SetElement(i, sum);
                }
            }
            outputs.push_back(result);
        }
    }
    else if (input1.IsNDMatrix())
    {
        if (size == 1)
        {
            oml_MatrixNUtil3(eval, inputs, outputs, oml_sum);
        }
        else
        {
            oml_MatrixNUtil3(eval, inputs, outputs, oml_sum, 2);
        }
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Gets cumulative sum of elements of input [cumsum]
//------------------------------------------------------------------------------
bool oml_cumsum(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int dim = 0;
    size_t nargin = inputs.size();

    if (nargin != 1 && nargin != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (nargin == 2)
    {
        dim = (int) inputs[1].Scalar();

        if (!IsInteger(inputs[1].Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_DIM);

        if (dim < 0)
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DIM);
    }

    const Currency &input1 = inputs[0];

    if (input1.IsScalar())
    {
        outputs.push_back(input1.Scalar());
    }
    else if (input1.IsComplex())
    {
        outputs.push_back(input1.Complex());
    }
    else if (input1.IsMatrix() || input1.IsString())
    {
        const hwMatrix* mtx = input1.Matrix();
        hwMatrix* result;
        if (!mtx->Size())
            outputs.push_back(EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL));
        else if (dim > 2)
        {
            result = EvaluatorInterface::allocateMatrix(mtx);
            outputs.push_back(result);
        }
        else
        {
            int m = mtx->M(), n = mtx->N();

            if (!dim)
                dim = (m > 1 ? 1 : 2);

            result = EvaluatorInterface::allocateMatrix(m, n, mtx->Type());

            if (mtx->IsReal())
            {
                double sum;
                for (int i = 0; i < (dim == 1 ? n : m); ++i)
                {
                    sum = 0.0;
                    for (int j = 0; j < (dim == 1 ? m : n); ++j)
                    {
                        if (dim == 1)
                        {
                            sum += (*mtx)(j, i);
                            result->SetElement(j, i, sum);
                        }
                        else
                        {
                            sum += (*mtx)(i, j);
                            result->SetElement(i, j, sum);
                        }
                    }
                }
            }
            else
            {
                hwComplex sum;
                for (int i = 0; i < (dim == 1 ? n : m); ++i)
                {
                    sum = hwComplex(0, 0);
                    for (int j = 0; j < (dim == 1 ? m : n); ++j)
                    {
                        if (dim == 1)
                        {
                            sum += mtx->z(j, i);
                            result->SetElement(j, i, sum);
                        }
                        else
                        {
                            sum += mtx->z(i, j);
                            result->SetElement(i, j, sum);
                        }
                    }
                }
            }
            outputs.push_back(result);
        }
    }
    else if (input1.IsNDMatrix())
    {
        if (nargin == 1)
        {
            return oml_MatrixNUtil4(eval, inputs, outputs, oml_cumsum);
        }
        else
        {
            return oml_MatrixNUtil4(eval, inputs, outputs, oml_cumsum, 2);
        }
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Gets cumulative product of input elements along a specified dimension [cumprod]
//------------------------------------------------------------------------------
bool oml_cumprod(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int dim = 0;
    size_t nargin = inputs.size();

    if (nargin != 1 && nargin != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (nargin == 2)
    {
        dim = (int) inputs[1].Scalar();

        if (!IsInteger(inputs[1].Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_DIM);

        if (dim < 0)
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_DIM);
    }

    const Currency &input1 = inputs[0];

    if (input1.IsScalar())
    {
        outputs.push_back(input1.Scalar());
    }
    else if (input1.IsComplex())
    {
        outputs.push_back(input1.Complex());
    }
    else if (input1.IsMatrix() || input1.IsString())
    {
        const hwMatrix* mtx = input1.Matrix();
        hwMatrix* result;
        if (!mtx->Size())
            outputs.push_back(EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL));
        else if (dim > 2)
        {
            result = EvaluatorInterface::allocateMatrix(mtx);
            outputs.push_back(result);
        }
        else
        {
            int m = mtx->M(), n = mtx->N();

            if (!dim)
                dim = (m > 1 ? 1 : 2);

            result = EvaluatorInterface::allocateMatrix(m, n, mtx->Type());

            if (mtx->IsReal())
            {
                double prod;
                for (int i = 0; i < (dim == 1 ? n : m); ++i)
                {
                    prod = 1.0;
                    for (int j = 0; j < (dim == 1 ? m : n); ++j)
                    {
                        if (dim == 1)
                        {
                            prod *= (*mtx)(j, i);
                            result->SetElement(j, i, prod);
                        }
                        else
                        {
                            prod *= (*mtx)(i, j);
                            result->SetElement(i, j, prod);
                        }
                    }
                }
            }
            else
            {
                hwComplex prod;
                for (int i = 0; i < (dim == 1 ? n : m); ++i)
                {
                    prod = hwComplex(1, 0);
                    for (int j = 0; j < (dim == 1 ? m : n); ++j)
                    {
                        if (dim == 1)
                        {
                            prod *= mtx->z(j, i);
                            result->SetElement(j, i, prod);
                        }
                        else
                        {
                            prod *= mtx->z(i, j);
                            result->SetElement(i, j, prod);
                        }
                    }
                }
            }
            outputs.push_back(result);
        }
    }
    else if (input1.IsNDMatrix())
    {
        if (nargin == 1)
        {
            return oml_MatrixNUtil4(eval, inputs, outputs, oml_cumprod);
        }
        else
        {
            return oml_MatrixNUtil4(eval, inputs, outputs, oml_cumprod, 2);
        }
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns smallest integer greater than or equal to input [ceil]
//------------------------------------------------------------------------------
bool oml_ceil(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!roundingFunc(inputs, outputs, ceil))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_ceil);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns largest integer that is not greater than given input [floor]
//------------------------------------------------------------------------------
bool oml_floor(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!roundingFunc(inputs, outputs, floor))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_floor);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns nearest integer to X [round]
//------------------------------------------------------------------------------
bool oml_round(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!roundingFunc(inputs, outputs, round))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_round);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns sin of input in radians [sin]
//------------------------------------------------------------------------------
bool oml_sin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, sin, hwComplex::sin))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_sin);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns cosine of input in radians [cos]
//------------------------------------------------------------------------------
bool oml_cos(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, cos, hwComplex::cos))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_cos);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns tangent of input in radians [tan]
//------------------------------------------------------------------------------
bool oml_tan(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, tan, hwComplex::tan))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_tan);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns inverse sine of input in radians [asin]
//------------------------------------------------------------------------------
bool oml_asin(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!conditionFunc(inputs, outputs, asin, hwComplex::asin, hwComplex::asin_c, absAtMostOne))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_asin);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns inverse cosine of input in radians [acos]
//------------------------------------------------------------------------------
bool oml_acos(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!conditionFunc(inputs, outputs, acos, hwComplex::acos, hwComplex::acos_c, absAtMostOne))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_acos);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns inverse tangent of input in radians [atan]
//------------------------------------------------------------------------------
bool oml_atan(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, atan, hwComplex::atan))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_atan);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns hyperbolic sine of input in radians [sinh]
//------------------------------------------------------------------------------
bool oml_sinh(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, sinh, hwComplex::sinh))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_sinh);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns hyperbolic cosine of input in radians [cosh]
//------------------------------------------------------------------------------
bool oml_cosh(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, cosh, hwComplex::cosh))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_cosh);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns hyperbolic tangent of input in radians [tanh]
//------------------------------------------------------------------------------
bool oml_tanh(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, tanh, hwComplex::tanh))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_tanh);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns inverse hyperbolic sine of input in radians [asinh]
//------------------------------------------------------------------------------
bool oml_asinh(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!noConditionFunc(inputs, outputs, asinh, hwComplex::asinh))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_asinh);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns inverse hyperbolic cosine of input in radians [acosh]
//------------------------------------------------------------------------------
bool oml_acosh(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!conditionFunc(inputs, outputs, acosh, hwComplex::acosh,
                       hwComplex::acosh_c, atLeastOne))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_acosh);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns inverse hyperbolic tangent of input in radians [atanh]
//------------------------------------------------------------------------------
bool oml_atanh(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!conditionFunc(inputs, outputs, atanh, hwComplex::atanh,
                       hwComplex::atanh_c, absAtMostOne))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_atanh);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes the base-2 logarithm of input [log2]
//------------------------------------------------------------------------------
bool oml_log2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (eval.GetNargoutValue() > 1)
    {
        if (inputs[0].IsScalar())
        {
            int exp;
            double f = std::frexp(inputs[0].Scalar(), &exp);
            outputs.push_back(f);
            outputs.push_back(static_cast<double>(exp));
        }
        else if (inputs[0].IsComplex())
        {
            hwComplex cplx = inputs[0].Complex();
            int exp;
            std::frexp(cplx.Mag(), &exp);
            double p2 = std::pow(2.0, exp);

            outputs.push_back(hwComplex(cplx.Real() / p2, cplx.Imag() / p2));
            outputs.push_back(static_cast<double>(exp));
        }
        else if (inputs[0].IsMatrix() || inputs[0].IsString())
        {
            const hwMatrix* mtx = inputs[0].Matrix();
            if (mtx->IsRealData())
            {
                hwMatrix* fmtx = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);
                hwMatrix* emtx = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);
                int exp;
                for (int i = 0; i < mtx->Size(); ++i)
                {
                    (*fmtx)(i) = std::frexp(realval(mtx,i), &exp);
                    (*emtx)(i) = static_cast<double>(exp);
                }
                outputs.push_back(fmtx);
                outputs.push_back(emtx);
            }
            else
            {
                hwMatrix* fmtx = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::COMPLEX);
                hwMatrix* emtx = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);
                int exp;
                for (int i = 0; i < mtx->Size(); ++i)
                {
                    const hwComplex& cplx = mtx->z(i);
                    int exp;
                    std::frexp(cplx.Mag(), &exp);
                    double p2 = std::pow(2.0, exp);
                    fmtx->z(i) = hwComplex(cplx.Real() / p2, cplx.Imag() / p2);
                    (*emtx)(i) = static_cast<double>(exp);
                }
                outputs.push_back(fmtx);
                outputs.push_back(emtx);
            }
        }
        else if (inputs[0].IsNDMatrix())
        {
            return oml_MatrixNUtil1(eval, inputs, outputs, oml_log2);
        }
        else
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }
    else
    {
        if (!conditionFunc(inputs, outputs, log2, hwComplex::log2, hwComplex::log2_c, nonnegative))
        {
            return oml_MatrixNUtil1(eval, inputs, outputs, oml_log2);
        }
    }
    return true;
}
//------------------------------------------------------------------------------
// Computes the base-10 logarithm of input/elements of input [log10]
//------------------------------------------------------------------------------
bool oml_log10(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!conditionFunc(inputs, outputs, log10, hwComplex::log10, hwComplex::log10_c, nonnegative))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_log10);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes natural log of input/each element of input [log]
//------------------------------------------------------------------------------
bool oml_log(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!conditionFunc(inputs, outputs, log, hwComplex::log, hwComplex::log_c, nonnegative))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_log);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns square root of input [sqrt]
//------------------------------------------------------------------------------
bool oml_sqrt(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (!conditionFunc(inputs, outputs, sqrt, hwComplex::sqrt, hwComplex::sqrt_c, nonnegative))
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_sqrt);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns absolute value of input [abs]
//------------------------------------------------------------------------------
bool oml_abs(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    const Currency &input = inputs[0];

    if (input.IsScalar())
    {
        double operand = input.Scalar();
        outputs.push_back(Currency(abs(operand)));

    }
    else if (input.IsComplex())
    {
        hwComplex cplx = input.Complex();
        outputs.push_back(cplx.Mag());
    }
    else if (input.IsMatrix() || input.IsString())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

        result->Abs(*mtx);
        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_abs);
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true if successful in displaying the result [disp]
//------------------------------------------------------------------------------
bool oml_print(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1) throw OML_Error(OML_ERR_NUMARGIN);

    // Need to mark the currency as coming from 'disp' command so we know not to
    // prepend 'ans = ' at output time. This is mainly for issues like eval('disp(3)')
    Currency cur (inputs[0]);
    cur.DispOutput();

    eval.PrintResult(cur);
        
    return true;
}
//------------------------------------------------------------------------------
// Returns input to logical type [logical]
//------------------------------------------------------------------------------
bool oml_logical(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (!input.IsMatrix() && !input.IsScalar() && !input.IsComplex())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_DATA);

    const hwMatrix* mtx = input.ConvertToMatrix();
    hwMatrix* res = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

    if (mtx->IsRealData())
    {
        for (size_t j = 0; j < mtx->Size(); j++)
        {
            double d = (*mtx)((const int)j);
            if (IsNaN_T(d))
                throw OML_Error(HW_ERROR_NOTCONVNANTOLOG);

            if (d == 0.0)
                (*res)((const int)j) = (0.0);
			else
                (*res)((const int)j) = (1.0);
        }
    }
    else
        throw OML_Error(HW_ERROR_NOTCONVCOMPTOLOG);

    Currency result(res);
    result.SetMask(Currency::MASK_LOGICAL);
    outputs.push_back(result);

    return true;
}
//------------------------------------------------------------------------------
// Called when 'quit' or 'exit' is called [exit/quit]
//------------------------------------------------------------------------------
bool oml_exit(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    eval.SetQuit(true);
	eval.OnSaveOnExit();
	return true;
}
//------------------------------------------------------------------------------
// Returns true if given string is a variable name with wild cards
//------------------------------------------------------------------------------
static bool IsVarnameWithWildcard(const std::string& str)
{
    std::string safecopy = str;
    safecopy.erase(std::remove_if(safecopy.begin(), safecopy.end(), [](char c)->bool { return c=='?'|| c=='*'; }),
                   safecopy.end());
    return isvarname(safecopy);
}
//------------------------------------------------------------------------------
// Helper method for regex
//------------------------------------------------------------------------------
static std::string WildcardToRegex(const std::string& str)
{
    std::stringstream ss;
    std::string::const_iterator iter;
    for (iter = str.cbegin(); iter != str.cend(); ++iter)
    {
        if (*iter == '*')
            ss << ".*";
        else if (*iter == '?')
            ss << '.';
        else
            ss << *iter;
    }
    return ss.str();
}
//------------------------------------------------------------------------------
// Obliterates variable(s) from memory in the current session [clear]
//------------------------------------------------------------------------------
bool oml_clear(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    enum
    {
        ALL,
        FUNCS,
        GLOBALS,
        VARIABLES
    } toclear = ALL;

    if (nargin)
    {
        // read through options

        size_t i = 1;
        std::string str = readString(inputs[0]);
        if (str == "-a" || str == "-all" || (str == "all" && !eval.Contains(str)))
            toclear = ALL;
        else if (str == "-f" || str == "-functions" || str == "functions")
            toclear = FUNCS;
        else if (str == "-g" || str == "-global" || str == "global")
            toclear = GLOBALS;
        else if (str == "-v" || str == "-variables" || str == "variables")
            toclear = VARIABLES;
        else
            i = 0;

        if (i == inputs.size())
        {
            switch (toclear)
            {
            case ALL:
                eval.ClearFunctions();
                eval.ClearGlobals();
                eval.ClearVariables();
                break;
            case FUNCS:
                eval.ClearFunctions();
                break;
            case GLOBALS:
                eval.ClearGlobals();
                break;
            case VARIABLES:
                eval.ClearVariables();
                break;
            }
        }
        else
        {
            for (; i < inputs.size(); ++i)
            {
                str = readString(inputs[i]);
                // if str has wildcards just clear the string -- it's more efficient to search through
                // a map for one than check a regex against every element of one
                if (std::count_if(str.cbegin(), str.cend(), [](char c){ return c == '?' || c == '*'; }))
                {
                    if (IsVarnameWithWildcard(str))
                    {
                        std::regex regexp(WildcardToRegex(str));
                        switch (toclear)
                        {
                        case ALL:
                            eval.Clear(regexp);
                            break;
                        case FUNCS:
                            eval.ClearFromFunctions(regexp);
                            break;
                        case GLOBALS:
                            eval.ClearFromGlobals(regexp);
                            break;
                        case VARIABLES:
                            eval.ClearFromVariables(regexp);
                            break;
                        }
                    }
                }
                else
                {
                    switch (toclear)
                    {
                    case ALL:
                        eval.Clear(str);
                        break;
                    case FUNCS:
                        eval.ClearFromFunctions(str);
                        break;
                    case GLOBALS:
                        eval.ClearFromGlobals(str);
                        break;
                    case VARIABLES:
                        eval.ClearFromVariables(str);
                        break;
                    }
                }
            }
        }
    }
    else
    {
        eval.ClearVariables();
        eval.ClearGlobals();
    }

    return true;
}
bool oml_end(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    outputs.push_back(eval.GetContextEndValue());
    return true;
}
//------------------------------------------------------------------------------
// List variables in the current session matching given pattern [who]
//------------------------------------------------------------------------------
bool oml_who(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
	size_t nargin = inputs.size();

	std::vector<std::string> names    = eval.GetVariableNames();
    size_t                   numnames = (names.empty()) ? 0 : names.size();

    HML_CELLARRAY* ret_val = EvaluatorInterface::allocateCellArray();
    assert(ret_val);
	ret_val->Dimension(static_cast<const int>(numnames), 1, HML_CELLARRAY::REAL);

	bool isGlob = false;
	std::string str;

	if (nargin==0) // Return all variables if no pattern is specified
    {
        if (numnames == 0)
            BuiltInFuncsUtils::SetWarning(eval, "No variables in current scope");

		if (getNumOutputs(eval))
        {
			for (int j = 0; j < static_cast<int>(numnames); ++j)
				(*ret_val)(j, 0) = Currency(names[j]);				
			outputs.push_back(ret_val);
		} 
        else if (numnames > 0)
        {
			eval.PrintResult(std::string("Variables in current scope:\n"));
			for (size_t i = 0; i < numnames; ++i)
				eval.PrintResult(names[i] + '\n');
		}
		return true;
	} 
    
    if(nargin == 1)
    {
		if(!inputs[0].IsString())
			throw OML_Error(OML_ERR_STRING);

		str = readString(inputs[0]);
		if (str == "-g" || str == "-global" || str == "global"){
			isGlob = true;
			str = "";
		}
	} else if(nargin == 2){
		if(!inputs[0].IsString() || !inputs[1].IsString())
			throw OML_Error(OML_ERR_STRING);

		str = readString(inputs[1]);

		if (str == "-g" || str == "-global" || str == "global")
			isGlob = true;
		else
			throw OML_Error("Error: 2nd string input must be 'global'");

		str = readString(inputs[0]); // set str to the searched string
	}
	
	bool searchStr = false;

	if(str != ""){
		searchStr = true;
	}

	if(isGlob){
		if (getNumOutputs(eval))
		{
			int k = 0;
			for (int j = 0; j < names.size(); j++){ // how big does the cell have to be
				if(eval.IsGlobal(names[j])){
					if(searchStr){
						std::regex regexp(WildcardToRegex(str));
						if(std::regex_match(names[j].cbegin(), names[j].cend(), regexp)){
							k++;
						}
					} else {
						k++;
					}
				}
			}
			ret_val->Dimension(k, 1, HML_CELLARRAY::REAL);
			k=0;
			for(int j = 0; j < names.size(); j++){ // put values in cell
				if(eval.IsGlobal(names[j])){
					if(searchStr){
						std::regex regexp(WildcardToRegex(str));
						if(std::regex_match(names[j].cbegin(), names[j].cend(), regexp)){
							(*ret_val)(k, 0) = Currency(names[j]);
							k++;
						}
					} else {
						(*ret_val)(k, 0) = Currency(names[j]);
						k++;
					}
				}
			}
			outputs.push_back(ret_val);
		}
		else
		{
			if (names.size())
			{
				eval.PrintResult(std::string("Variables in current scope:\n"));
				for (size_t i = 0; i < names.size(); ++i){
					if(eval.IsGlobal(names[i])){
						if(searchStr){
							std::regex regexp(WildcardToRegex(str));
							if(std::regex_match(names[i].cbegin(), names[i].cend(), regexp)){
								eval.PrintResult(names[i] + '\n');
							}
						} else {
							eval.PrintResult(names[i] + '\n');
						}
					}
				}
			} 
            else 
				BuiltInFuncsUtils::SetWarning(eval, "No variables in current scope");
		}
	} else {

		std::vector<std::string> names = eval.GetVariableNames();

		HML_CELLARRAY* ret_val = EvaluatorInterface::allocateCellArray();

		if (getNumOutputs(eval))
		{
			int k = 0; 
			for (int j = 0; j < names.size(); j++){ // how big does the cell have to be
				if(searchStr){
					std::regex regexp(WildcardToRegex(str));
					if(std::regex_match(names[j].cbegin(), names[j].cend(), regexp)){
						k++;
					}
				} else {
					k++;
				}
			}
			ret_val->Dimension(k, 1, HML_CELLARRAY::REAL);
			k = 0;
			for (int j = 0; j < names.size(); j++){ // put values in cell
				if(searchStr){
					std::regex regexp(WildcardToRegex(str));
					if(std::regex_match(names[j].cbegin(), names[j].cend(), regexp)){
						(*ret_val)(k, 0) = Currency(names[j]);
						k++;
					}
				} else {
					(*ret_val)(j, 0) = Currency(names[j]);
					k++;
				}
			}
			outputs.push_back(ret_val);
		}
		else
		{
			if (names.size())
			{
				eval.PrintResult(std::string("Variables in current scope:\n"));
				for (size_t i = 0; i < names.size(); ++i){
					if(searchStr){
						std::regex regexp(WildcardToRegex(str));
						if(std::regex_match(names[i].cbegin(), names[i].cend(), regexp)){
							eval.PrintResult(names[i] + '\n');
						}
					} else {
						eval.PrintResult(names[i] + '\n');
					}
				}
			} 
            else 
				BuiltInFuncsUtils::SetWarning(eval, "No variables in current scope.");
		}
	}
    return true;
}
//------------------------------------------------------------------------------
// Returns true values [true]
//------------------------------------------------------------------------------
bool oml_true(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return createCommonMatrix(eval, inputs, outputs, getTrue());
}
//------------------------------------------------------------------------------
// Returns false values [false]
//------------------------------------------------------------------------------
bool oml_false(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    return createCommonMatrix(eval, inputs, outputs, getFalse());
}
//------------------------------------------------------------------------------
// Returns ref count - for test purposes only [refcnt]
//------------------------------------------------------------------------------
bool oml_refcnt(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int count = 0;

    if (inputs.size() == 1)
    {
        if (inputs[0].IsMatrixOrString())
        {
            const hwMatrix* mtx = inputs[0].Matrix();
            count = mtx->GetRefCount();

            // There's a strange problem here.  Because the matrix is in a currency that's in the inputs vector, the ref count
            // is artificially high. -- JDS

            count--;
        }
        else if (inputs[0].IsNDMatrix())
        {
            const hwMatrixN* mtx = inputs[0].MatrixN();
            count = mtx->GetRefCount();

            // There's a strange problem here.  Because the matrix is in a currency that's in the inputs vector, the ref count
            // is artificially high. -- JDS

            count--;
        }
        else if (inputs[0].IsStruct())
        {
            StructData* sd = inputs[0].Struct();
            count = sd->GetRefCount();

            count--;
        }
    }

    outputs.push_back(count);
    return true;
}
//------------------------------------------------------------------------------
// Get/set last error message [lasterr]
//------------------------------------------------------------------------------
bool oml_lasterr(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() > 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs.size() == 0)
    {
        int nargout = eval.GetNargoutValue();

        if (nargout > 1)
            throw OML_Error(OML_ERR_NUMARGOUT);

        outputs.push_back(EvaluatorInterface::GetLastErrorMessage());
    }
    else /* (inputs.size() == 1) */
    {
        int nargout = eval.GetNargoutValue();

        if (nargout > 0)
            throw OML_Error(OML_ERR_NUMARGOUT);

        if (!inputs[0].IsString())
            throw OML_Error(HW_ERROR_INPUTSTRING);

        EvaluatorInterface::SetLastErrorMessage(inputs[0].StringVal());
    }
    return true;
}
//------------------------------------------------------------------------------
// Get/set last warning [lastwarn]
//------------------------------------------------------------------------------
bool oml_lastwarn(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() > 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs.size() == 0)
    {
        int nargout = eval.GetNargoutValue();

        if (nargout > 1)
            throw OML_Error(OML_ERR_NUMARGOUT);

        outputs.push_back(EvaluatorInterface::GetLastWarning());
    }
    else /* (inputs.size() == 1) */
    {
        int nargout = eval.GetNargoutValue();

        if (nargout > 0)
            throw OML_Error(OML_ERR_NUMARGOUT);

        if (!inputs[0].IsString())
            throw OML_Error(HW_ERROR_INPUTSTRING);

        EvaluatorInterface::SetLastWarning(inputs[0].StringVal());
    }
    return true;
}
//------------------------------------------------------------------------------
// Continue command in a loop [continue]
//------------------------------------------------------------------------------
bool oml_continue(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    outputs.push_back(Currency(-1.0, Currency::TYPE_CONTINUE));
    return true;
}
//------------------------------------------------------------------------------
// Break command in a loop [break]
//------------------------------------------------------------------------------
bool oml_break(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    outputs.push_back(Currency(-1.0, Currency::TYPE_BREAK));
    return true;
}
//------------------------------------------------------------------------------
// Given vectors of x and y coordinates, return matrices xx and yy corresponding 
// to a full 2D grid [meshgrid]
//------------------------------------------------------------------------------
bool oml_meshgrid(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int nargin = eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin != 1 && nargout != nargin)
        throw OML_Error(OML_ERR_NUMARGINOUT);

    for (int i = 0; i < nargin; ++i)
    {
        if (!inputs[i].IsMatrix() && !inputs[i].IsScalar())
            throw OML_Error(OML_ERR_REALVECTOR, i+1, OML_VAR_TYPE);

        if (inputs[i].IsMatrix() && (!inputs[i].Matrix()->IsVector() || !inputs[i].Matrix()->IsReal()))
            throw OML_Error(OML_ERR_REALVECTOR, i+1, OML_VAR_TYPE);
    }

    if (nargin < 3 && nargout < 3)  // create 2D object
    {
        const hwMatrix* x = NULL;
        const hwMatrix* y = NULL;

        if (nargin == 2)
        {
            Currency x_cur = inputs[0];
            Currency y_cur = inputs[1];

            x = x_cur.ConvertToMatrix();
            y = y_cur.ConvertToMatrix();
        }
        else if (nargin == 1)
        {
            Currency x_cur = inputs[0];

            x = x_cur.ConvertToMatrix();
            y = x_cur.ConvertToMatrix();
        }

        int size1 = y->Size();
        int size2 = x->Size();

        hwMatrix* result1 = EvaluatorInterface::allocateMatrix(size1, size2, hwMatrix::REAL);
        hwMatrix* result2 = EvaluatorInterface::allocateMatrix(size1, size2, hwMatrix::REAL);

        for (int j = 0; j < size1; j++)
        {
            for (int k = 0; k < size2; k++)
            {
                (*result1)(j,k) = (*x)(k);
                (*result2)(j,k) = (*y)(j);
            }
        }

        outputs.push_back(result1);
        outputs.push_back(result2);
    }
    else  // create ND object
    {
        // determine object dimensions
        std::vector<int> dims(nargout);

        const hwMatrix* x = NULL;
        const hwMatrix* y = NULL;
        const hwMatrix* z = NULL;

        if (nargin == 1)
        {
            x = inputs[0].Matrix();
            y = x;

            if (nargout == 3)
                z = x;
        }
        else if (nargin == 2)
        {
            x = inputs[0].Matrix();
            y = inputs[1].Matrix();
        }
        else
        {
            x = inputs[0].Matrix();
            y = inputs[1].Matrix();
            z = inputs[2].Matrix();
        }

        dims[0] = y->Size();
        dims[1] = x->Size();

        if (nargout == 3)
            dims[2] = z->Size();

        // process first output
        hwMatrixN* vec1 = EvaluatorInterface::allocateMatrixN();
        vec1->Convert2DtoND(*x);

        std::vector<int> reshapedim(2);
        reshapedim[0] = 1;
        reshapedim[1] = dims[1];
        vec1->Reshape(reshapedim);

        std::vector<double> repdim(nargout);
        repdim[0] = static_cast<double> (dims[0]);
        repdim[1] = 1.0;

        if (nargout == 3)
            repdim[2] = static_cast<double> (dims[2]);

        std::vector<Currency> newinputs;
        std::vector<Currency> newoutputs;
        newinputs.push_back(vec1);
        newinputs.push_back(repdim);
        oml_repmat(eval, newinputs, newoutputs);
        outputs.push_back(newoutputs[0]);

        // process second output
        hwMatrixN* vec2 = EvaluatorInterface::allocateMatrixN();
        vec2->Convert2DtoND(*y);
        reshapedim[0] = dims[0];
        reshapedim[1] = 1;
        vec2->Reshape(reshapedim);
        repdim[0] = 1.0;
        repdim[1] = static_cast<double> (dims[1]);
        newinputs.clear();
        newoutputs.clear();
        newinputs.push_back(vec2);
        newinputs.push_back(repdim);
        oml_repmat(eval, newinputs, newoutputs);
        outputs.push_back(newoutputs[0]);

        // process third output
        if (nargout == 3)
        {
            hwMatrixN* vec = EvaluatorInterface::allocateMatrixN();
            vec->Convert2DtoND(*z);

            reshapedim[0] = 1;
            reshapedim[1] = 1;
            reshapedim.push_back(dims[2]);
            vec->Reshape(reshapedim);

            repdim[0] = static_cast<double> (dims[0]);
            repdim[1] = static_cast<double> (dims[1]);
            repdim[2] = 1.0;

            newinputs.clear();
            newoutputs.clear();
            newinputs.push_back(vec);
            newinputs.push_back(repdim);
            oml_repmat(eval, newinputs, newoutputs);
            outputs.push_back(newoutputs[0]);
        }
    }

	return true;
}
//------------------------------------------------------------------------------
// Returns nd arrays for given input vectors [ndgrid]
//------------------------------------------------------------------------------
bool oml_ndgrid(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int nargin = eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin < 1)
        throw OML_Error(OML_ERR_NUMARGIN);
    
    if (nargin != 1 && nargout > nargin)
        throw OML_Error(OML_ERR_NUMARGINOUT);
    
    for (int i = 0; i < nargin; ++i)
    {
        if (!inputs[i].IsMatrix() && !inputs[i].IsScalar())
            throw OML_Error(OML_ERR_REALVECTOR, i+1, OML_VAR_TYPE);

        if (inputs[i].IsMatrix() && (!inputs[i].Matrix()->IsVector() || !inputs[i].Matrix()->IsReal()))
            throw OML_Error(OML_ERR_REALVECTOR, i+1, OML_VAR_TYPE);
    }

    if (nargin < 3 && nargout < 3)  // create 2D object
    {
        const hwMatrix* x = NULL;
        const hwMatrix* y = NULL;

        if (nargin == 2)
        {
            Currency x_cur = inputs[0];
            Currency y_cur = inputs[1];

            x = x_cur.ConvertToMatrix();
            y = y_cur.ConvertToMatrix();
        }
        else if (nargin == 1)
        {
            Currency x_cur = inputs[0];

            x = x_cur.ConvertToMatrix();
            y = x_cur.ConvertToMatrix();
        }

        int size1 = x->Size();
        int size2 = y->Size();

        hwMatrix* result1 = EvaluatorInterface::allocateMatrix(size1, size2, hwMatrix::REAL);
        hwMatrix* result2 = EvaluatorInterface::allocateMatrix(size1, size2, hwMatrix::REAL);

        for (int j = 0; j < size2; j++)
        {
            for (int k = 0; k < size1; k++)
            {
                (*result1)(k,j) = (*x)(k);
                (*result2)(k,j) = (*y)(j);
            }
        }

        outputs.push_back(result1);
        outputs.push_back(result2);
    }
    else  // create ND object
    {
        // determine object dimensions
        std::vector<int> dims((nargin == 1) ? nargout : nargin);

        if (nargin == 1)
        {
            for (int i = 0; i < nargout; ++i)
                dims[i] = inputs[0].Matrix()->Size();
        }
        else
        {
            for (int i = 0; i < nargin; ++i)
                dims[i] = inputs[i].Matrix()->Size();
        }

        // process first argument
        hwMatrixN* vec = EvaluatorInterface::allocateMatrixN();
        vec->Convert2DtoND(*inputs[0].Matrix());

        std::vector<int> reshapedim(2);
        reshapedim[0] = dims[0];
        reshapedim[1] = 1;
        vec->Reshape(reshapedim);

        std::vector<double> repdim(dims.size());

        repdim[0] = 1.0;

        for (int i = 1; i < dims.size(); ++i)
            repdim[i] = static_cast<double> (dims[i]);

        std::vector<Currency> newinputs;
        std::vector<Currency> newoutputs;
        newinputs.push_back(vec);
        newinputs.push_back(repdim);
        oml_repmat(eval, newinputs, newoutputs);
        outputs.push_back(newoutputs[0]);

        // process remaining arguments
        for (int i = 1; i < nargout; ++i)
        {
            hwMatrixN* vec2 = EvaluatorInterface::allocateMatrixN();

            if (nargin == 1)
                vec2->Convert2DtoND(*inputs[0].Matrix());
            else
                vec2->Convert2DtoND(*inputs[i].Matrix());

            std::vector<int> reshapedim2(i+1);

            for (int j = 0; j < i; ++j)
                reshapedim2[j] = 1;

            reshapedim2[i] = dims[i];
            vec2->Reshape(reshapedim2);

            for (int j = 0; j < dims.size(); ++j)
            {
                if (j == i)
                    repdim[j] = 1.0;
                else
                    repdim[j] = static_cast<double> (dims[j]);
            }

            newinputs.clear();
            newoutputs.clear();
            newinputs.push_back(vec2);
            newinputs.push_back(repdim);
            oml_repmat(eval, newinputs, newoutputs);
            outputs.push_back(newoutputs[0]);
        }
    }

	return true;
}
//------------------------------------------------------------------------------
// Returns current cpu time used by the application [cputime]
//------------------------------------------------------------------------------
bool oml_cputime(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	double cputime=getCPUTime();
	outputs.push_back(cputime);
	return true;
}
//------------------------------------------------------------------------------
// Records input/output commands [diary]
//------------------------------------------------------------------------------
bool oml_diary(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if(inputs.size()<1)
		throw OML_Error(OML_ERR_NUMARGIN);

	if(inputs.size()==0)
		eval.SetDiary(!eval.IsDiaryOpen());
	else
	{
		if(inputs[0].IsString())
		{
			std::string diary = readString(inputs[0]);
			if(diary=="on")
				eval.SetDiary(true);
			else if(diary=="off")
				eval.SetDiary(false);
			else
			{
				if (diary.find(".txt") != std::string::npos) {
					eval.SetDiary(diary);
				} else if (diary.find(".") != std::string::npos)
				{
					throw OML_Error(HW_ERROR_INPUTONOROFF);
				} else if (diary.find(":") != std::string::npos || diary.find("/") != std::string::npos || diary.find("\\") != std::string::npos) {
						throw OML_Error(HW_ERROR_ILLEGALCHARFORFILE);
				} else {
					eval.SetDiary(diary + ".txt");
				}
			}
		}
		else 
			throw OML_Error(HW_ERROR_INPUTONOROFF);
	}
	
	return true;
}
//------------------------------------------------------------------------------
// Returns true and gives the number of elements in the given input [numel]
//------------------------------------------------------------------------------
bool oml_numel(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)  //\todo: Multiple input arguments not implemented
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    Currency input = inputs[0];

    if (input.IsScalar() || input.IsComplex() || input.IsFunctionHandle())
    {
		outputs.push_back(1);
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();
        outputs.push_back(mtx->Size());
    }
    else if (input.IsNDMatrix())
    {
        const hwMatrixN* mtx = input.MatrixN();
        outputs.push_back(mtx->Size());
    }
	else if (input.IsString())
	{
		std::string st = input.StringVal();
        if (st.empty())
        {
            outputs.push_back(0);
        }
        else 
        {
            if (input.IsMultilineString())
            {
                // Strip new lines from the string
                for (std::string::iterator itr = st.begin(); itr != st.end();)
                {
                    if (*itr == '\n')
                    {
                        itr = st.erase(itr);
                    }
                    else
                    {
                        ++itr;
                    }
                }
            }
            unsigned char* my_ptr = (unsigned char*)st.c_str();
		    outputs.push_back(utf8_strlen(my_ptr));
        }
	}
	else if (input.IsCellArray())
    {
        HML_CELLARRAY* cells = input.CellArray();
        outputs.push_back(cells->Size());
    }
	else if (input.IsStruct())
    {
        StructData* sd = input.Struct();
        outputs.push_back(sd->M()*sd->N());
    }
    else
    {
        throw OML_Error(HW_ERROR_UNKOWNTYPE);
    }

    return true;
}
//------------------------------------------------------------------------------
// Circularly shifts given input [circshift]
//------------------------------------------------------------------------------
bool oml_circshift(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	int nargin = static_cast<int> (inputs.size());

	if (nargin != 2 && nargin != 3)
		throw OML_Error(OML_ERR_NUMARGIN);

    Currency target = inputs[0];

    if (target.IsMatrixOrString())
    {
        const Currency& shift = inputs[1];
        Currency dim;

        if (!shift.IsInteger() && !shift.IsVector())
            throw OML_Error(OML_ERR_SCALARVECTOR, 2);

        int dim_val = 1;

        if (nargin == 3)
        {
            dim = inputs[2];

            if (!shift.IsInteger())
                throw OML_Error(OML_ERR_SCALAR, 2);

            if (!dim.IsPositiveInteger())
                throw OML_Error(OML_ERR_INVALID_INDEX, 3);

            dim_val = (int)dim.Scalar();
        }

        std::vector<double> shift_vals;

        if (shift.IsRealVector())
        {
            shift_vals = shift.Vector();
        }
        else if (shift.IsScalar())
        {
            for (int j=1; j<dim_val; j++)
                shift_vals.push_back(0);

            shift_vals.push_back(shift.Scalar());
        }
        else
        {
            throw OML_Error("Invalid shift value");
        }

        if (shift_vals.size() > 2)
            throw OML_Error(OML_ERR_UNSUPPORTDIM, 2);

        int shift_val_1 = (int)shift_vals[0];
        int shift_val_2 = 0;

        if (shift_vals.size() > 1)
            shift_val_2 = (int)shift_vals[1];

        const hwMatrix* target_mtx = target.Matrix();
        hwMatrix*       result     = EvaluatorInterface::allocateMatrix(target_mtx);

        if (target_mtx->IsReal())
        {
            for (int j=0; j<target_mtx->M(); j++)
            {
                int shift_index_1 = (j+shift_val_1) % target_mtx->M();

                if (shift_index_1 < 0)
                    shift_index_1 += target_mtx->M();

                for (int k=0; k<target_mtx->N(); k++)
                {
                    int shift_index_2 = (k+shift_val_2) % target_mtx->N();

                    if (shift_index_2 < 0)
                        shift_index_2 += target_mtx->N();

                    (*result)(shift_index_1, shift_index_2) = (*target_mtx)(j,k); 
                }
            }
        }
        else
        {
            for (int j=0; j<target_mtx->M(); j++)
            {
                int shift_index_1 = (j+shift_val_1) % target_mtx->M();

                if (shift_index_1 < 0)
                    shift_index_1 += target_mtx->M();

                for (int k=0; k<target_mtx->N(); k++)
                {
                    int shift_index_2 = (k+shift_val_2) % target_mtx->N();

                    if (shift_index_2 < 0)
                        shift_index_2 += target_mtx->N();

                    result->z(shift_index_1, shift_index_2) = target_mtx->z(j,k); 
                }
            }
        }

        Currency temp(result);

        if (target.IsString())
            temp.SetMask(Currency::MASK_STRING);

        outputs.push_back(temp);
    }
    else if (target.IsNDMatrix())
    {
        if (nargin == 2)
        {
            if (inputs[1].IsVector())
                return oml_MatrixN_circshift(eval, inputs, outputs, oml_circshift);
            else if (inputs[1].IsInteger())
                return oml_MatrixNUtil4(eval, inputs, outputs, oml_circshift);
            else
                throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_TYPE);
        }
        else
        {
            return oml_MatrixNUtil4(eval, inputs, outputs, oml_circshift, 3);
        }
    }
    else
    {
        throw OML_Error(OML_ERR_MATRIX, 1);
    }

	return true;
}
//------------------------------------------------------------------------------
// Returns true if input is a valid string [checksyntax]
//------------------------------------------------------------------------------
bool oml_checksyntax(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error(OML_ERR_NUMARGIN);

	Currency input = inputs[0];

	if (!input.IsString())
		throw OML_Error(OML_ERR_STRING);

	std::string temp = eval.IsValidString(input.StringVal());

	if (temp.empty())
		outputs.push_back(true);
	else
		outputs.push_back(false);

	outputs.push_back(temp);

	return true;
}
//------------------------------------------------------------------------------
// ast command [ast]
//------------------------------------------------------------------------------
bool oml_ast(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error(OML_ERR_NUMARGIN);

	Currency input = inputs[0];

	if (!input.IsFunctionHandle() && !input.IsString())
	{
		throw OML_Error("Error: invalid input; must be a string or function handle");
	}

	FunctionInfo* fi   = NULL;
	FUNCPTR       fptr = NULL;

	if (input.IsString())
		eval.FindFunctionByName(input.StringVal(), &fi, &fptr);
	else if (input.IsFunctionHandle())
		fi = input.FunctionHandle();

	if (fi)
	{
		if (fi->IsBuiltIn())
			outputs.push_back("builtin function");
		else
			outputs.push_back(fi->GetAST());
	}
	else if (fptr)
	{
		outputs.push_back("builtin function");
	}
	else
	{
		throw OML_Error(HW_ERROR_UNKOWNFUN);
	}

	return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
bool oml_writepfile(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs[0].IsString() && inputs[1].IsString())
		eval.WritePFile(inputs[0].StringVal(), inputs[1].StringVal());
	return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
void analyze_helper(EvaluatorInterface eval, const std::vector<std::string>& targets, std::vector<std::string>& files_needed)
{
	for (int j=0; j<targets.size(); j++)
	{
		std::string target = targets[j];

		if (std::find(files_needed.begin(), files_needed.end(), target) != files_needed.end())
			continue;
		else
			files_needed.push_back(target);

		Currency result = eval.Analyze(target);

		HML_CELLARRAY* cells = result.CellArray();
		std::vector<std::string> new_targets;

		for (int k=0; k<cells->Size(); k++)
		{
			std::string target = (*cells)(k).StringVal();

			if (std::find(files_needed.begin(), files_needed.end(), target) != files_needed.end())
				continue; // we already have this one

			FunctionInfo* fi = NULL;

			eval.FindFunctionByName(target, &fi, NULL);

			if (fi)
				new_targets.push_back(fi->FileName());
		}

		analyze_helper(eval, new_targets, files_needed);
	}
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
bool oml_analyze(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error(OML_ERR_NUMARGIN);

	if (inputs[0].IsString())
	{
		std::vector<std::string> needed_files;
		std::vector<std::string> targets;

		targets.push_back(inputs[0].StringVal());

		analyze_helper(eval, targets, needed_files);

		HML_CELLARRAY* cells = eval.allocateCellArray(1, (int)needed_files.size());

		for (int j=0; j<needed_files.size(); j++)
			(*cells)(j) = Currency(needed_files[j]);

		outputs.push_back(cells);
	}

	return true;
}
//------------------------------------------------------------------------------
// Gets help module name [helpmodule]
//------------------------------------------------------------------------------
bool oml_helpmodule(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	Interpreter intBase(eval);
	std::string func_name = inputs[0].StringVal();

	std::string module = intBase.GetHelpModule(func_name);
	outputs.push_back(module);

	return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
bool oml_helptest(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	Interpreter intBase(eval);

	const char* ROOT_DIRECTORY = "OML_HELP";
	char *c = getenv(ROOT_DIRECTORY);
	if(c == nullptr)
		throw OML_Error(HW_ERROR_ENVVARFAIL);

	std::string base = c;

	eval.PrintResult("Help files missing or misplaced for the functions:");

	std::vector<std::string> funclist = intBase.GetFunctionNames();
	std::string funcstring;
	for (std::vector<std::string>::iterator it = funclist.begin() ; it != funclist.end(); ++it)
	{
	std::string userDocLocation = base + intBase.GetHelpModule((*it)) + "/" + (*it) + ".htm";

		std::ifstream myLocation;
		myLocation.open(userDocLocation);
		std::string shortHelp;

		if(!myLocation.is_open())
		{
			shortHelp.append(userDocLocation);
		    eval.PrintResult(shortHelp);

		}
		myLocation.close();
	}
	return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
std::string help_format_output(const std::string &word)
{
	std::string result;
	result = word;

	std::string find_list[] = { "&lt;", "&lt", "&amp;", "&amp", "&gt;", "&gt" };
	vector<std::string> find_strings(find_list, find_list +6);

	std::string replace_list[] = { "<", "<", "&", "&", ">", ">" };
	vector<std::string> replace_strings(replace_list, replace_list +6);

	std::size_t word_position = 0;
	std::size_t subst_length = 0;

	for (unsigned ix = 0; ix < find_strings.size(); ++ix)
    {
		std::string findstr;
		findstr.assign(find_strings[ix]);
		std::string replstr;
		replstr.assign(replace_strings[ix]);

		word_position = 0;
		while((word_position = result.find(findstr, word_position)) != std::string::npos)
		{
			result.replace (word_position,findstr.length(),replstr);
		}
    }
	return result;
}
//------------------------------------------------------------------------------
// updated help function - attempting to reduce the complex number of specific strings 
// being searched for in the original (now) help1() function.
//------------------------------------------------------------------------------
bool oml_help(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	Interpreter intBase(eval);
	std::stringstream ss;
	std::string hwRootDir;
	if(inputs.size() == 0)
		throw OML_Error("Please enter a function name as a string.");

	if (!inputs[0].IsString())
		throw OML_Error("Please enter a function name as a string.");

	FunctionInfo* fi = NULL;
	FUNCPTR       fp = NULL;

	if (eval.GetExperimental())
	{
		eval.FindFunctionByName(inputs[0].StringVal(), &fi, &fp);

		if (fi)
		{
			std::string help_str = fi->HelpString();

			if (help_str.size())
			{
				eval.PrintResult(help_str);
				return true;
			}
		}
	}

	const char* ROOT_DIRECTORY = "OML_HELP";
	char *c = getenv(ROOT_DIRECTORY);

	std::string userDocLocation;

	if (c)
	{
		std::string base = c;
		userDocLocation = base + intBase.GetHelpModule(inputs[0].StringVal()) + "/" + inputs[0].StringVal() + ".htm";
	}

	std::ifstream myLocation;
	std::string word, wordout;
	std::string shortHelp;
	bool endShortdesc = false;
	bool begin_skip = false;
	bool word_done = false;
	bool file_done = false;
	std::size_t found_begin, found_end;
	bool wordout_print = false;
	bool header_found = false;				// ignore all content before the "topictitle1"
	bool header_completed = false;			// ignore all content before the end of "topictitle1" tag
	bool title_CRLF = false;				// used to add blank line before and after titles
#if INPUTVAR_COLON
	bool InputVar_Colon = false;			// Disable \n 
#endif 

	myLocation.open(userDocLocation);
	if(!myLocation.is_open())
	{
        std::string msg = HW_ERROR_NOFMATCH;
        msg += "\n" + BuiltInFuncsUtils::Normpath(userDocLocation);
		throw OML_Error(msg);
	}

	while(file_done != true){
		myLocation >> word;
		word_done = false;
		std::size_t word_position = 0;

		if(myLocation.eof())
		{
			file_done = true;
			continue;
		}
		if(!header_found)
		{
			if(word.find("topictitle1") == string::npos)
				continue;  // ignore everything until the topictitle1 tag is detected.
			else
			{
				header_found = true;
				myLocation >> word;	// Handle topictitle1 tag for help()
				word_position = word.find("title1") + 8;
				title_CRLF = true;  // try this to get crlf after function name
			}
		}

		while(!word_done)
		{
			if(begin_skip)
			{   
				found_end = word.find(">", word_position);
				if(found_end == string::npos) // std::string::npos)  // not found
				{
					word_done = true;
					continue;
				}
				else
				{
					begin_skip = false;
					word_position = found_end+1;  // need to make sure we don't go past end of word.
				}
				if(word.find("sectiontitle\">") != string::npos)
				{
					title_CRLF = true;
					// wordout.append("\n");
				}
#if INPUTVAR_COLON
				if(word.find("varname\">") != string::npos)
				{
					InputVar_Colon = true;
					// wordout.append("\n");
				}
#endif
			}
			
			// look for next tag
			found_begin = word.find("<", word_position);

			// done when we hit this tag
			if(word.find("</html>", found_begin) != std::string::npos) // this will need attention later
			{
				file_done = true;
				break;
			}
			// done when we hit this tag
			if(word.find("<?STOPOUT?>", found_begin) != std::string::npos) // this will need attention later
			{
				file_done = true;
				break;
			}
	
			if(found_begin == 0)
			{
				begin_skip = true;
				word_position = found_begin + 1;
				continue;
			}
			else
			{
				if(found_begin == string::npos) // std::string::npos)  // not found
				{
#if 0
					// Test code to detect the end of the help text we want to display and terminate the help display 
					std::string temp = word.substr(word_position, word.length() - word_position);
					// End help display at the Example/Examples or History section if not previously terminated.
					if(temp.find("Example") != string::npos || 
					   temp.find("Examples") != string::npos ||
					   temp.find("History") != string::npos) 
					{
						file_done = true;
						break;
					}
#endif					
					wordout.append(word.substr(word_position, word.length() - word_position));
					if (myLocation.peek() != '\n')
						wordout.append(" ");
					else
					{
						wordout_print = true;
						if(title_CRLF)
						{ 
#if INPUTVAR_COLON
							if(InputVar_Colon)
								InputVar_Colon = false;
							else
#endif 
							{
								wordout.append(":\n");
								title_CRLF = false;
							}
						}
#if INPUTVAR_COLON
						if(InputVar_Colon)  // try to put colon after Input variable name.
						{
							wordout.append(": ");
							// InputVar_Colon = false;
						}
#endif 
					}
					word_done = true;
					continue;
				}
				else if(found_begin > word_position)
				{
					if(title_CRLF)
					{   
#if INPUTVAR_COLON
						if(InputVar_Colon)
							InputVar_Colon = false;
						else
#endif
							wordout.append("\n");  // append \n before sectiontitle
					}
					wordout.append(word.substr(word_position, found_begin - word_position));
				}
				word_position = found_begin+1; // need to make sure we don't go past end of word.
				begin_skip = true;
				// wordout_print = true;
			}
				
		}
		if(wordout_print && wordout.length())
		{
			if(wordout.find("\nExample:\n") != std::string::npos)
			{
				file_done = true;
				continue;
			}
			wordout = help_format_output(wordout);
			shortHelp = shortHelp + wordout + "\n";
			wordout.clear();
			wordout_print = false;
		}
	}
	int nargout = eval.GetNargoutValue();
	if (nargout == 0)
		eval.PrintResult(shortHelp);
	else
		outputs.push_back(shortHelp);
	return true;
}
//------------------------------------------------------------------------------
// Create a row vector with logarithmically spaced elements [logspace]
//------------------------------------------------------------------------------
bool oml_logspace(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	int nargin = (int)inputs.size();

	if(nargin > 3 || nargin < 1){
		throw OML_Error(OML_ERR_NUMARGIN);
	}

	if(!inputs[0].IsScalar() && !inputs[0].IsComplex())
		throw OML_Error(OML_ERR_SCALAR,1,OML_VAR_INDEX);
	
	if(!inputs[1].IsScalar() && !inputs[1].IsComplex())
		throw OML_Error(OML_ERR_SCALAR,2,OML_VAR_INDEX);

	double numpts;

	if(nargin==3){
		if(!inputs[2].IsScalar())
			throw OML_Error(OML_ERR_SCALAR,3,OML_VAR_INDEX);
		numpts = inputs[2].Scalar();
	} else {
		numpts = 50;
	}

	hwMatrix *out;

	if(inputs[0].IsScalar() && inputs[1].IsScalar()){
		double startpt; double endpt; 
		startpt = inputs[0].Scalar();
		endpt   = inputs[1].Scalar();


		out = EvaluatorInterface::allocateMatrix(1, (int)numpts, hwMatrix::REAL);

		for(int i = 0; i < numpts; i++){
			(*out)(0,i) = pow(10,(startpt+((endpt-startpt)*(i/((floor(numpts)-1))))));
		}
	} else {
		hwComplex startpt; hwComplex endpt; 
		startpt = inputs[0].Complex();
		endpt   = inputs[1].Complex();

		
		out = EvaluatorInterface::allocateMatrix(1, (int)numpts, hwMatrix::COMPLEX);
		
		
		for(int i = 0; i < numpts; i++){
			out->z(i) = hwComplex::pow(10.0,(startpt+((endpt-startpt)*(i/((floor(numpts)-1))))));
		}
	}
	

	

	outputs.push_back(out);

	return true;
}
// helper methods

//------------------------------------------------------------------------------
// function replacement method
//------------------------------------------------------------------------------
bool keywordFunc(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, Currency out)
{
    if (inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
bool noConditionFunc(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, double (*realFunc)(double),
    hwComplex (*cplxFunc)(const hwComplex&))
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsScalar())
    {
        double operand = input.Scalar();
        outputs.push_back(Currency((*realFunc)(operand)));
    }
    else if (input.IsComplex())
    {
        hwComplex cplx = input.Complex();
        hwComplex result = (*cplxFunc)(cplx);
        outputs.push_back(result);
    }
    else if (input.IsMatrix() || input.IsString())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

        if (mtx->IsReal())
        {
            for (int k = 0; k < mtx->Size(); k++)
            {
                result->SetElement(k, (*realFunc)((*mtx)(k)));
            }
        }
        else
        {
            for (int k = 0; k < mtx->Size(); k++)
            {
                hwComplex cplx = mtx->z(k);
                result->SetElement(k, (*cplxFunc)(cplx));
            }
        }

        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return false;   // ND case is handled upon return
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
bool conditionFunc(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, double (*realFunc)(double),
    hwTComplex<double> (*cplxFunc)(const hwComplex&), hwComplex (*cnvrtFunc)(double), bool (*conditionFunc)(double))
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsScalar())
    {
        double operand = input.Scalar();
        if (IsNaN_T(operand))
        {
            outputs.push_back(std::numeric_limits<double>::quiet_NaN());
        }
        else
        {
            if ((*conditionFunc)(operand))
                outputs.push_back(Currency((*realFunc)(operand)));
            else
                outputs.push_back(Currency((*cnvrtFunc)(operand)));
        }
    }
    else if (input.IsComplex())
    {
        hwComplex cplx = input.Complex();
        if (IsNaN_T(cplx))
            outputs.push_back(hwComplex(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
        else
            outputs.push_back((*cplxFunc)(cplx));
    }
    else if (input.IsMatrix() || input.IsString())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

        if (mtx->IsReal())
        {
            for (int k = 0; k < mtx->Size(); k++)
            {
                double dbl = (*mtx)(k);
                if (IsNaN_T(dbl))
                {
                    result->SetElement(k, std::numeric_limits<double>::quiet_NaN());
                }
                else
                {
                    if ((*conditionFunc)(dbl))
                        result->SetElement(k, (*realFunc)(dbl));
                    else
                        result->SetElement(k, (*cnvrtFunc)(dbl));
                }
            }
        }
        else
        {
            for (int k = 0; k < mtx->Size(); k++)
            {
                hwComplex cplx = mtx->z(k);
                if (IsNaN_T(cplx))
                    result->SetElement(k, hwComplex(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()));
                else
                    result->SetElement(k, (*cplxFunc)(cplx));
            }
        }

        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return false;   // ND case is handled upon return
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool roundingFunc(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, double (*roundFunc)(double))
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input = inputs[0];

    if (input.IsLogical())
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    if (input.IsScalar())
    {
        double operand = input.Scalar();
        outputs.push_back((*roundFunc)(operand));
    }
    else if (input.IsComplex())
    {
        hwComplex cplx = input.Complex();
        double real = (*roundFunc)(cplx.Real());
        double imag = (*roundFunc)(cplx.Imag());
        outputs.push_back(hwComplex(real, imag));
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

        if (mtx->IsReal())
        {
            for (int k = 0; k < mtx->Size(); k++)
            {
                result->SetElement(k, (*roundFunc)((*mtx)(k)));
            }
        }
        else
        {
            for (int k = 0; k < mtx->Size(); k++)
            {
                hwComplex cplx = mtx->z(k);
                double real = (*roundFunc)(cplx.Real());
                double imag = (*roundFunc)(cplx.Imag());
                result->SetElement(k, hwComplex(real, imag));
            }
        }

        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return false;   // ND case is handled upon return
    }
    else if (input.IsString())
    {
        // nothing to round
        Currency out = input;
        out.SetMask(Currency::MASK_DOUBLE);
        outputs.push_back(out);
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING);
    }

    return true;
}
//------------------------------------------------------------------------------
// if no inputs are provided, a singular newval is returned
// if a matrix input is provided, it's values are used as dimensions of the output matrix where every element is newval
// if singular values are given as input, they are treated individually as matrix dimensions
// first dimension is length of columns, second is length of rows, etc.
// if complex values are given, their real components are used
// if decimal values are given, they are truncated
// negative dimensions are treated as 0
//------------------------------------------------------------------------------
bool createCommonMatrix(EvaluatorInterface& eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, Currency newval)
{
    size_t size = inputs.size();

    if (!size)
    {
        outputs.push_back(newval);
        return true;
    }

    const Currency &input1 = inputs[0];

    // determine matrix class: hwMatrix or hwMatrixN
    bool createND = false;

    if (size == 1 && input1.IsVector() && input1.Matrix()->Size() > 2)
    {
        createND = true;
    }
    else if (size > 2)
    {
        createND = true;
    }

    // construct matrix
    if (!createND)
    {
        // first determine matrix dimensions
        double m, n;
        if (size == 1)
        {
            const Currency &input1 = inputs[0];
            if (input1.IsScalar())
                m = n = input1.Scalar();
            else if (input1.IsComplex())
                m = n = input1.Complex().Real();
            else if (input1.IsMatrix())
            {
                const hwMatrix *mtx = input1.Matrix();
                int matrixSize = mtx->Size();
                if (!matrixSize)
                {
                    throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIMS);
                }
                else if (!mtx->IsRealData())
                {
                    throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);
                }
                else if (matrixSize == 2)
                {
                    m = realval(mtx, 0);
                    n = realval(mtx, 1);
                }
                else
                {
                    throw OML_Error(OML_ERR_UNSUPPORTDIM, 1);
                }
            }
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else if (size == 2)
        {
            const Currency &input1 = inputs[0];
            const Currency &input2 = inputs[1];

            if (input1.IsScalar())
                m = input1.Scalar();
            else if (input1.IsComplex())
                m = input1.Complex().Real();
            else
                throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);

            if (input2.IsScalar())
                n = input2.Scalar();
            else if (input2.IsComplex())
                n = input2.Complex().Real();
            else
                throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);
        }
        else
        {
            throw OML_Error(OML_ERR_UNSUPPORTDIM);
        }

        if (!(::isfinite(m) && ::isfinite(n)))
            throw OML_Error(HW_ERROR_DIMFINITE);

        // take care of negative inputs
        if (m < 0)
            m = 0;
        if (n < 0)
            n = 0;

        // now set all elements in new matrix
        hwMatrix *result;
        if (newval.IsScalar())
        {
            result = EvaluatorInterface::allocateMatrix((int) m, (int) n, newval.Scalar());
        }
        else if (newval.IsComplex())
        {
            result = EvaluatorInterface::allocateMatrix((int) m, (int) n, newval.Complex());
        }
        else
        {
            throw OML_Error(HW_ERROR_INVTYPESETELE);
        }

        Currency out(result);
        out.SetMask(newval.GetMask()); // preserve logical/string mask
        outputs.push_back(out);
    }
    else    // createND
    {
		std::vector<int> dimensions;
        double dim;

        if (size == 1)
        {
            const hwMatrix* dimens = input1.Matrix();
            int numDim = dimens->Size();

            for (int j = 0; j < numDim; ++j)
            {
                dim = realval(dimens, j);

                if (!(::isfinite(dim)))
                    throw OML_Error(HW_ERROR_DIMFINITE);

                if (dim < 0.0)
                    dim = 0.0;

                dimensions.push_back(static_cast<int>(dim));
            }
        }
        else    // size > 2
        {
            for (int j = 0; j < size; ++j)
            {
                if (inputs[j].IsScalar())
                {
                    dim = inputs[j].Scalar();
                }
                else if (inputs[j].IsComplex())
                {
                    dim = inputs[j].Complex().Real();
                }
                else
                {
                    throw OML_Error(OML_ERR_NATURALNUM, j+1, OML_VAR_DIM);
                }

                if (!(::isfinite(dim)))
                    throw OML_Error(HW_ERROR_DIMFINITE);

                if (dim < 0.0)
                    dim = 0.0;

                dimensions.push_back(static_cast<int>(dim));
            }
        }

        hwMatrixN* result = EvaluatorInterface::allocateMatrixN();

        if (newval.IsScalar())
        {
            result->Dimension(dimensions, hwMatrixN::REAL);
            result->SetElements(newval.Scalar());
        }
        else if (newval.IsComplex())
        {
            result->Dimension(dimensions, hwMatrixN::COMPLEX);
            result->SetElements(newval.Complex());
        }
        else
        {
            throw OML_Error(HW_ERROR_INVTYPESETELE);
        }

        Currency out(result);
        out.SetMask(newval.GetMask()); // preserve logical/string mask
        outputs.push_back(out);
    }

    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool limitFunc(EvaluatorInterface& eval, const std::vector<Currency> &inputs, std::vector<Currency> &outputs, float (*floatFunc)(), double (*doubleFunc)())
{
    size_t nargin = inputs.size();
    bool usedouble = true;

    if (nargin)
    {
        const Currency &last = inputs[nargin - 1];
        if (last.IsString())
        {
            --nargin;
            std::string className = readOption(eval, last);
            if (className == "float" || className == "single")
                usedouble = false;
            else if (className != "double")
                throw OML_Error(HW_ERROR_INVALIDOPTION(className));
        }
    }

    double val;
    if (usedouble)
        val = (*doubleFunc)();
    else
        val = (double) (*floatFunc)();

    // determine object type
    bool createND = false;

    if (nargin > 2)
        createND = true;

    // construct matrix
    if (!createND)
    {
        if (nargin == 0)
        {
            outputs.push_back(val);
        }
        else 
        {
            int m, n;

            if (!inputs[0].IsInteger())
                throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);

            m = static_cast<int> (inputs[0].Scalar());
            m = max(m, 0);

            if (nargin == 1)
            {
                n = m;
            }
            else    // nargin == 2
            {
                if (!inputs[1].IsInteger())
                    throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);

                n = static_cast<int> (inputs[1].Scalar());
                n = max(n, 0);
            }

            hwMatrix *out = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
            out->SetElements(val);
            outputs.push_back(out);
        }
    }
    else
    {
        std::vector<int> dimensions;
        double dim;

        for (int j = 0; j < nargin; ++j)
        {
            if (inputs[j].IsScalar())
            {
                dim = inputs[j].Scalar();
            }
            else if (inputs[j].IsComplex())
            {
                dim = inputs[j].Complex().Real();
            }
            else
            {
                throw OML_Error(OML_ERR_NATURALNUM, j+1, OML_VAR_DIM);
            }

            if (!(::isfinite(dim)))
                throw OML_Error(HW_ERROR_DIMFINITE);

            if (dim < 0.0)
                dim = 0.0;

            dimensions.push_back(static_cast<int>(dim));
        }

        hwMatrixN* out = EvaluatorInterface::allocateMatrixN();
        out->Dimension(dimensions, hwMatrixN::REAL);
        out->SetElements(val);
        outputs.push_back(out);
    }

    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool scalarToLogicalFunc(const std::vector<Currency> &inputs, std::vector<Currency> &outputs, double (*checker)(double))
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(doScalarToLogical(inputs[0], checker));
    return true;
}
//------------------------------------------------------------------------------
// for functions like union, intersect, setxor, etc.
//------------------------------------------------------------------------------
bool sortBasedOperation(EvaluatorInterface &eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, bool allowRepeatIndices,
    std::deque<std::string> (*workerCell)(std::deque<std::string>&,std::deque<std::string>&),
    std::deque<hwMatrix> (*workerRow) (std::deque<hwMatrix>& ,std::deque<hwMatrix>&),
    std::deque<double> (*workerReal)(std::deque<double>&, std::deque<double>&),
    std::deque<hwComplex> (*workerCplx)(std::deque<hwComplex>& ,std::deque<hwComplex>&))
{
    size_t size = inputs.size();

    if (size < 2 || size > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];
    const Currency &input2 = inputs[1];
    bool rows = false;
    if (size == 3)
    {
        if (inputs[2].IsString())
        {
            std::string opt = readOption(eval, inputs[2]);
            if (opt == "rows")
                rows = true;
            else
                throw OML_Error(HW_ERROR_INVALIDOPTION(opt));
        }
    }

    int nargout = getNumOutputs(eval);
    Currency u;
    HML_CELLARRAY *acell = nullptr;
    HML_CELLARRAY *bcell = nullptr;

    std::unique_ptr<hwMatrix>a, b;

    if (input1.IsCellArray())
    {
        acell = input1.CellArray();
    }
    else
    {
        a.reset(safeMatrixCopyFromInput(input1, true));
    }

    if (input2.IsCellArray())
    {
        bcell = input2.CellArray();
    }
    else
    {
        b.reset(safeMatrixCopyFromInput(input2, true));
    }

    if (acell || bcell)
    {
        if (rows)
            throw OML_Error(HW_ERROR_NOUNIONBYROWCELLINP);

        if ((acell && !isstr(acell)) || (bcell && !isstr(bcell)))
            throw OML_Error(HW_ERROR_CELLINPMUSTSTR);

        std::deque<std::string> avals, bvals, vals;
        if (acell)
        {
            for (int i = 0; i < acell->Size(); ++i)
            {
                std::string str = readString((*acell)(i));
                avals.push_back(str);
            }
        }
        else
        {
            for (int i = 0; i < a->M(); ++i)
            {
                hwMatrix *row = readRow(eval, a.get(), i);
                std::string str = Currency(row).StringVal();
                avals.push_back(str);
            }
        }

        if (bcell)
        {
            for (int i = 0; i < bcell->Size(); ++i)
            {
                std::string str = readString((*bcell)(i));
                bvals.push_back(str);
            }
        }
        else
        {
            for (int i = 0; i < b->M(); ++i)
            {
                hwMatrix *row = readRow(eval, b.get(), i);
                std::string str = Currency(row).StringVal();
                bvals.push_back(str);
            }
        }

        // index vectors -- these will be disjoint
        bool returnrow = false;
        if (acell && acell->M() == 1)
            returnrow = true;
        else if (bcell && bcell->M() == 1)
            returnrow = true;

        vals = (*workerCell)(avals, bvals);

        u = containerToCellArray(vals, returnrow);
        outputs.push_back(u);

        if (nargout > 1)
        {
            std::vector<int> ai, bi;

            for (size_t i = 0; i < vals.size(); ++i)
            {
                std::string tofind = vals[i];
                int j;
                for (j = (int)bvals.size() - 1; j >= 0; --j)
                {
                    std::string str = bvals[j];
                    if (str == tofind)
                    {
                        bi.push_back(j + 1);
                        break;
                    }
                }

                if (allowRepeatIndices || j == -1)
                {
                    for (j = (int)avals.size() - 1; j >= 0; --j)
                    {
                        std::string str = avals[j];
                        if (str == tofind)
                        {
                            ai.push_back(j + 1);
                            break;
                        }
                    }
                }
            }

            outputs.push_back(containerToMatrix(ai, returnrow));
            outputs.push_back(containerToMatrix(bi, returnrow));
        }
        return true;
    }

    if (rows)
    {
        int n = a->N();
        if (n != b->N())
            throw OML_Error(HW_ERROR_WHENUNIONINDROWVECSAMECOL);

        std::deque<hwMatrix> avals, bvals, vals;

        hwMatrix *val = EvaluatorInterface::allocateMatrix();
        Currency valcur(val);
        for (int i = 0; i < a->M(); i++)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, a->ReadRow(i, *val));
            avals.push_back(*val);
        }

        for (int i = 0; i < b->M(); i++)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, b->ReadRow(i, *val));
            bvals.push_back(*val);
        }

        vals = (*workerRow)(avals, bvals);
        u = dequeToMatrix(vals);

        if (input1.IsString() || input2.IsString())
            u.SetMask(Currency::MASK_STRING);

        outputs.push_back(u);

        if (nargout > 1)
        {
            std::vector<int> ai, bi;
            bool returnrow = a->M() == 1 && b->M() == 1;

            for (size_t i = 0; i < vals.size(); ++i)
            {
                const hwMatrix &tofind = vals[i];
                int j;
                for (j = b->M() - 1; j >= 0; --j)
                {
                    hwMatrix *row = readRow(eval, b.get(), j);
                    if (*row == tofind)
                    {
                        bi.push_back(j + 1);
                        delete row;
                        break;
                    }
                    delete row;
                }

                if (allowRepeatIndices || j == -1)
                {
                    for (j = a->M() - 1; j >= 0; --j)
                    {
                        hwMatrix *row = readRow(eval, a.get(), j);
                        if (*row == tofind)
                        {
                            ai.push_back(j + 1);
                            delete row;
                            break;
                        }
                        delete row;
                    }
                }
            }

            outputs.push_back(containerToMatrix(ai, returnrow));
            outputs.push_back(containerToMatrix(bi, returnrow));
        }
    }
    else // !rows
    {
        if (checkMakeComplex(eval, a.get(), b.get()))
        {
            std::deque<hwComplex> avals, bvals, vals;

            for (int i = 0; i < a->Size(); i++)
            {
                avals.push_back(a->z(i));
            }

            for (int i = 0; i < b->Size(); i++)
            {
                bvals.push_back(b->z(i));
            }

            vals = (*workerCplx)(avals, bvals);

            u = containerToMatrix(vals, a->M() == 1 || b->M() == 1);
            outputs.push_back(u);

            if (nargout > 1)
            {
                std::vector<int> ai, bi;
                bool returnrow = a->M() == 1 && b->M() == 1;

                for (size_t i = 0; i < vals.size(); ++i)
                {
                    hwComplex tofind = vals[i];
                    int j;
                    for (j = b->Size() - 1; j >= 0; --j)
                    {
                        if (b->z(j) == tofind)
                        {
                            bi.push_back(j + 1);
                            break;
                        }
                    }

                    if (allowRepeatIndices || j == -1)
                    {
                        for (j = a->Size() - 1; j >= 0; --j)
                        {
                            if (a->z(j) == tofind)
                            {
                                ai.push_back(j + 1);
                                break;
                            }
                        }
                    }
                }

                outputs.push_back(containerToMatrix(ai, returnrow));
                outputs.push_back(containerToMatrix(bi, returnrow));
            }
        }
        else
        {
            std::deque<double> avals, bvals, vals;

            for (int i = 0; i < a->Size(); i++)
            {
                avals.push_back((*a)(i));
            }

            for (int i = 0; i < b->Size(); i++)
            {
                bvals.push_back((*b)(i));
            }

            vals = (*workerReal)(avals, bvals);

            u = containerToMatrix(vals, a->M() == 1 || b->M() == 1);
            if (input1.IsString() || input2.IsString())
                u.SetMask(Currency::MASK_STRING);
            outputs.push_back(u);

            if (nargout > 1)
            {
                std::vector<int> ai, bi;
                bool returnrow = a->M() == 1 && b->M() == 1;

                for (size_t i = 0; i < vals.size(); ++i)
                {
                    double tofind = vals[i];
                    int j;
                    for (j = b->Size() - 1; j >= 0; --j)
                    {
                        if ((*b)(j) == tofind)
                        {
                            bi.push_back(j + 1);
                            break;
                        }
                    }

                    if (allowRepeatIndices || j == -1)
                    {
                        for (j = a->Size() - 1; j >= 0; --j)
                        {
                            if ((*a)(j) == tofind)
                            {
                                ai.push_back(j + 1);
                                break;
                            }
                        }
                    }
                }

                outputs.push_back(containerToMatrix(ai, returnrow));
                outputs.push_back(containerToMatrix(bi, returnrow));
            }
        }
    }

    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool allowComplexToLogical(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, bool (*checker)(double), bool complexBoth)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];

    if (cur.IsScalar())
    {
        outputs.push_back((*checker)(cur.Scalar()) ? getTrue() : getFalse());
    }
    else if (cur.IsComplex())
    {
        const hwComplex& cplx = cur.Complex();
        if (complexBoth)
            outputs.push_back((*checker)(cplx.Real()) && (*checker)(cplx.Imag()) ? getTrue() : getFalse());
        else
            outputs.push_back((*checker)(cplx.Real()) || (*checker)(cplx.Imag()) ? getTrue() : getFalse());
    }
    else if (cur.IsMatrix() || cur.IsString())
    {
        const hwMatrix* mtx = cur.Matrix();
        hwMatrix* out = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

        if (mtx->IsRealData())
        {
            for (int i = 0; i < mtx->Size(); ++i)
            {
                (*out)(i) = (*checker)(realval(mtx, i));
            }
        }
        else
        {
            for (int i = 0; i < mtx->Size(); ++i)
            {
                const hwComplex& cplx = mtx->z(i);
                if (complexBoth)
                    (*out)(i) = (*checker)(cplx.Real()) && (*checker)(cplx.Imag());
                else
                    (*out)(i) = (*checker)(cplx.Real()) || (*checker)(cplx.Imag());
            }
        }

        Currency outcur(out);
        outcur.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(outcur);
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mtx = cur.MatrixN();
        hwMatrixN* out = EvaluatorInterface::allocateMatrixN();

        out->Dimension(mtx->Dimensions(), hwMatrixN::REAL);

        if (mtx->IsRealData())
        {
            for (int i = 0; i < mtx->Size(); ++i)
            {
                (*out)(i) = (*checker)((*mtx)(i));
            }
        }
        else
        {
            for (int i = 0; i < mtx->Size(); ++i)
            {
                const hwComplex& cplx = mtx->z(i);
                if (complexBoth)
                    (*out)(i) = (*checker)(cplx.Real()) && (*checker)(cplx.Imag());
                else
                    (*out)(i) = (*checker)(cplx.Real()) || (*checker)(cplx.Imag());
            }
        }

        Currency outcur(out);
        outcur.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(outcur);
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING);
    }

    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool ineqOperatorFunc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(eval.InequalityOperator(inputs[0], inputs[1], op));
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool eqOperatorFunc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(eval.EqualityOperator(inputs[0], inputs[1], op));
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool unOperatorFunc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    outputs.push_back(eval.UnaryOperator(inputs[0], op));
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool binOperatorFunc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency in1 = inputs[0];
	Currency in2 = inputs[1];

	if (in1.GetMask() == Currency::MASK_STRING)
		in1.SetMask(Currency::MASK_NONE);

	if (in2.GetMask() == Currency::MASK_STRING)
		in2.SetMask(Currency::MASK_NONE);

    outputs.push_back(eval.BinaryOperator(in1, in2, op));
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool binOperatorFuncVararg(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op)
{
    if (inputs.size() < 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency result = inputs[0];

    for (size_t i = 1; i < inputs.size(); ++i)
        result = eval.BinaryOperator(result, inputs[i], op);

    outputs.push_back(result);
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool logOperatorFunc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int op)
{
    if (!inputs.size())
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency result = inputs[0];

    for (size_t i = 1; i < inputs.size(); ++i)
        result = eval.LogicalOperator(result, inputs[i], op);

    outputs.push_back(result);
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool anyall(EvaluatorInterface& eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, double (*func)(EvaluatorInterface&, const hwMatrix*), bool defaultval)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& in1 = inputs[0];
    if (in1.IsCellArray())
    {
        Currency out;
        out.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(out);
        return true;
    }
    else if (in1.IsEmpty())
    {
        outputs.push_back(HW_MAKE_BOOL_CURRENCY(defaultval));
        return true;
    }
    else if (in1.IsNDMatrix())
    {
        return false;   // ND case is handled upon return
    }
    else if (!(in1.IsScalar() || in1.IsComplex() || in1.IsMatrix() || in1.IsString()))
    {
        outputs.push_back(getFalse());
        return true;
    }

    const hwMatrix* data = in1.ConvertToMatrix();
    int dim;

    if (nargin > 1)
    {
        if (!inputs[1].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

        dim = static_cast<int> (inputs[1].Scalar());
    }
    else
    {
        dim = data->M() == 1 ? 2 : 1;
    }

    outputs.push_back(oml_MatrixUtil(eval, data, dim, func));
    outputs[0].SetMask(Currency::MASK_LOGICAL);
    return true;
}
//------------------------------------------------------------------------------
// condition methods
//------------------------------------------------------------------------------
bool nonnegative(double inpt)
{
    return (inpt >= 0);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool absAtMostOne(double inpt)
{
    return (abs(inpt) <= 1);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool atLeastOne(double inpt)
{
    return inpt >= 1;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isint(const hwComplex &c)
{
    return isint(c.Real()) && isint(c.Imag());
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isnonnegint(double d)
{ 
    return isint(d) && (int) d >= 0; 
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool iszero(double d)
{
    return isint(d) && ((int) d) == 0;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool iszero(const hwComplex &c)
{
    return isint(c) && ((int) c.Real()) == 0 && ((int) c.Imag()) == 0;
}
//------------------------------------------------------------------------------
// input-reading methods
//------------------------------------------------------------------------------
bool boolFromCur(const Currency& cur)
{
    if (cur.IsScalar())
        return !iszero(cur.Scalar());
    throw OML_Error(HW_ERROR_OPTNOTSING);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwMatrix* matrixCopyFromInput(const Currency &input, bool allowString)
{
    if (input.IsScalar())
        return EvaluatorInterface::allocateMatrix(1, 1, input.Scalar());
    else if (input.IsComplex())
        return EvaluatorInterface::allocateMatrix(1, 1, input.Complex());
    else if (input.IsMatrix() || (allowString && input.IsString()))
        return EvaluatorInterface::allocateMatrix(input.Matrix());
    return nullptr;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwMatrix* safeMatrixCopyFromInput(const Currency& input, bool allowString)
{
    hwMatrix* mtx = matrixCopyFromInput(input, allowString);
    if (!mtx)
        throw OML_Error(allowString ? HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING : HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    return mtx;
}
//------------------------------------------------------------------------------
// Returns file id from input
//------------------------------------------------------------------------------
int getFileFromInput(EvaluatorInterface &eval, const Currency &input)
{
    if (input.IsString())
    {
        std::string nameInput = readString(input);
        int numfiles = eval.GetNumFiles();
        for (int i = 0; i < numfiles; i++)
        {
            std::string fname = eval.GetFileName(i);
            if (fname.length() && fname == nameInput)
                return i;
        }
        return -1;
    }
    else
    {
        if (!IsInteger(input.Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_FILEID);

        return (int) input.Scalar();
    }
}
//------------------------------------------------------------------------------
// Creates temporary file and returns file pointer
//------------------------------------------------------------------------------
FILE* makeTempFile()
{
#ifdef OS_WIN
    DWORD len = GetTempPath(0, nullptr);

    if (len != 0)
    {
        LPSTR dirpath = new CHAR[len + 1];
        len = GetTempPath(len + 1, dirpath);
        if (len != 0)
        {
            LPSTR path = new CHAR[MAX_PATH + 1];
            if (GetTempFileName(dirpath, "HW_TEMP_FILE_BUFFER", 0, path))
            {
                FILE* f = fopen(path, "w+b");
                delete [] path;
                if (f)
                {
                    delete [] dirpath;
                    return f;
                }
            }
        }
        delete [] dirpath;
    }
#else
    FILE* f = tmpfile();
    if (f)
        return f;
#endif
    throw OML_Error(HW_ERROR_PROBCREATTEMPF);
}
//------------------------------------------------------------------------------
// Reads tmp file
//------------------------------------------------------------------------------
std::string readTmpFile(FILE* file)
{
    std::stringstream ss;
    char buffer[256];
    fflush(file);
    rewind(file);
    while(!feof(file))
    {
        size_t numread = fread(buffer, sizeof(char), 256, file);
        if (!numread)
            break;
        
        ss.write(buffer, numread);
    }
    
    fclose(file);
    // TODO: delete file if on windows (flcose handles this on *nix because we used tmpfile() to create it)
    return ss.str();
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::string readOption(EvaluatorInterface& eval, const Currency &input)
{
    return convertToLower(readString(toCurrencyStr(eval, input, false, false)));
}
//------------------------------------------------------------------------------
// if matrices are used as inputs, an array is returned instead of just one
// assumes order of year, month, day, hour, minute, second
// the latter three are optional
// no support yet for date vector
//------------------------------------------------------------------------------
std::vector<std::array<double, 6> > datesFromInput(const std::vector<Currency> &inputs, int *m, int *n)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 6 || nargin == 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    // get dimensions of date
    *m = *n = 1;

    for (size_t i = 0; i < nargin; i++)
    {
        const Currency &cur = inputs[i];
        if (cur.IsMatrix() || cur.IsString())
        {
            const hwMatrix *mtx = cur.Matrix();

            if (!mtx->IsRealData())
                throw OML_Error(HW_ERROR_DATENOTCOMP);
            if (mtx->Size() != 1)
            {
                if (*m == 1 && *n == 1)
                {
                    *m = mtx->M();
                    *n = mtx->N();
                }
                else if (!(*m == mtx->M() && *n == mtx->N()))
                    throw OML_Error(HW_ERROR_MATRIXDIM);
            }
        }
    }

    size_t size = (size_t) (*m * *n);
    std::vector<std::array<double, 6> > dates;

    for (size_t i =  0; i < size; i++)
    {
        std::array<double, 6> d = {0, 0, 0, 0, 0, 0};
        dates.push_back(d);
    }

    for (size_t i = 0; i < nargin; i++)
    {
        const Currency &cur = inputs[i];
        if (cur.IsScalar() || (cur.IsComplex() && iszero(cur.Complex().Imag())))
        {
            for (size_t j = 0; j < size; j++)
            {
                dates[j][i] = cur.IsScalar() ? cur.Scalar() : cur.Complex().Real();
            }
        }
        else if (cur.IsComplex())
        {
            throw OML_Error(HW_ERROR_DATENOTCOMP);
        }
        else if (cur.IsMatrix() || cur.IsString())
        {
            const hwMatrix *mtx = cur.Matrix();
                for (size_t j = 0; j < size; j++)
                    dates[j][i] = (*mtx)((int) j);
            }
            else
                throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);
        }

    // month-related corrections
    for (size_t i = 0; i < size; i++)
    {
        // don't allow months less than 1
        if (dates[i][1] < 1.0)
            dates[i][1] = 1.0;

        if (!isint(dates[i][1]))
            throw OML_Error(HW_ERROR_MONTHMUSTINT);
    }

    return dates;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void get1DStringsFromInput(const Currency& cur, std::vector<std::string>& vec)
{
    if (cur.IsString())
    {
        if (cur.Matrix()->M() != 1)
            throw OML_Error(HW_ERROR_STRINPMUST1DIM);
        vec.push_back(cur.StringVal());
    }
    else if (cur.IsCellArray())
    {
        HML_CELLARRAY* cell = cur.CellArray();
        if (!isstr(cell))
            throw OML_Error(HW_ERROR_CELLELEMSTR);
        for (int i = 0; i < cell->Size(); ++i)
        {
            const Currency &elem = (*cell)(i);
            if (elem.Matrix()->M() != 1)
                throw OML_Error(HW_ERROR_TYPE1DSTR);
            vec.push_back((*cell)(i).StringVal());
        }
    }
    else
        throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);
}

// output-related methods
//------------------------------------------------------------------------------
// Gets number of outputs
//------------------------------------------------------------------------------
int getNumOutputs(const EvaluatorInterface &eval)
{
    return eval.GetNargoutValue();
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwMatrix* vectorToDiag(hwMatrix *mtx)
{
    int size = mtx->Size();
    hwMatrix *result;

    if (mtx->IsReal())
    {
        result = EvaluatorInterface::allocateMatrix(size, size, 0.0);
        for (int k = 0; k < size; k++)
        {
            (*result)(k, k) = (*mtx)(k);
        }
    }
    else
    {
        result = EvaluatorInterface::allocateMatrix(size, size, hwComplex(0.0, 0.0));
        for (int k = 0; k < size; k++)
        {
            result->z(k, k) = mtx->z(k);
        }
    }

    return result;
}
//------------------------------------------------------------------------------
// assumes each matrix in deq is a row vector of the same size
//------------------------------------------------------------------------------
hwMatrix* dequeToMatrix(const std::deque<hwMatrix> &deq)
{
    if (deq.empty())
        return EvaluatorInterface::allocateMatrix();

    hwMatrix *ret = EvaluatorInterface::allocateMatrix((int)deq.size(), deq[0].N(), hwMatrix::REAL);

    for (int i = 0; i < ret->M(); i++)
        ret->WriteRow(i, deq[i]);

    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwMatrix* dequeToMatrix(const std::deque<hwComplex> &deq, bool row)
{
    hwMatrix *ret;
    if (row)
        ret = EvaluatorInterface::allocateMatrix(1, (int)deq.size(), hwMatrix::COMPLEX);
    else
        ret = EvaluatorInterface::allocateMatrix((int)deq.size(), 1, hwMatrix::COMPLEX);

    int size = (int)deq.size();
    for (int i = 0; i < size; i++)
        ret->z(i) = deq[i];

    return ret;
}
//------------------------------------------------------------------------------
// Adds string mask and returns currency
//------------------------------------------------------------------------------
Currency addStringMask(const hwMatrix *m)
{
    Currency c((hwMatrix*)m);
    c.SetMask(Currency::MASK_STRING);
    return c;
}

// math methods
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double log2(double x)
{
    return (log(x) / log(2.0));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double round(double x)
{
    if (x >= 0.0)
        return floor(x + 0.5);
    else
        return ceil(x - 0.5);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double asinh(double x)
{
    return log(x + sqrt(x * x + 1));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double acosh(double x)
{
    return log(x + sqrt(x - 1) * sqrt(x + 1));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double atanh(double x)
{
    return 0.5 * (log(1 + x) - log(1 - x));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double mod(double x, double y)
{
    if (y == 0.0)
        return x;
    return x - y * floor(x / y);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double signum(double x)
{
    if (x < 0)
        return -1.0;
    if (x > 0)
        return 1.0;
    return 0;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double doubleMod(double a, double b)
{
    return a - b * (int) (a / b);
}
//------------------------------------------------------------------------------
// calculated via Euclidean algorithm
//------------------------------------------------------------------------------
int gcd(int q, int r, int *s, int *t)
{
    bool swapped = false;
    int temp, intdiv, s1 = 1, s2 = 0, t1 = 0, t2 = 1;
    if (q > r)
    {
        temp = q;
        q = r;
        r = temp;
        swapped = true;
    }
    while (r)
    {
        intdiv = q / r;

        temp = r;
        r = q % r;
        q = temp;

        temp = s1;
        s1 = s2;
        s2 = temp - intdiv * s2;

        temp = t1;
        t1 = t2;
        t2 = temp - intdiv * t2;
    }

    if (q < 0)
    {
        q *= -1;
        s1 *= -1;
        t1 *= -1;
    }

    if (swapped)
    {
        *s = t1;
        *t = s1;
    }
    else
    {
        *s = s1;
        *t = t1;
    }
    return q;
}
//------------------------------------------------------------------------------
// calculated via Euclidean algorithm
//------------------------------------------------------------------------------
hwComplex gcd(hwComplex q, hwComplex r, hwComplex *s, hwComplex *t)
{
    if (q.IsReal(HW_TOLERANCE) && r.IsReal(HW_TOLERANCE))
    {
        int ss, tt;
        int g = gcd((int) q.Real(), (int) r.Real(), &ss, &tt);
        s->Real() = ss;
        t->Real() = tt;
        s->Imag() = 0.0;
        t->Imag() = 0.0;
        return hwComplex(g, 0.0);
    }

    bool swapped = false;
    hwComplex temp, intdiv, s1, s2, t1, t2;
    s1 = t2 = hwComplex(1.0, 0.0);
    s2 = t1 = hwComplex(0.0, 0.0);

    if (q.Mag() > r.Mag())
    {
        temp = q;
        q = r;
        r = temp;
        swapped = true;
    }
    while (!isint(r) || ((int) r.Real() || (int) r.Imag()))
    {
        intdiv = q / r;
        makeInt(intdiv);

        temp = r;
        r = q - intdiv * r;
        q = temp;

        temp = s1;
        s1 = s2;
        s2 = temp - intdiv * s2;

        temp = t1;
        t1 = t2;
        t2 = temp - intdiv * t2;
    }

    if (swapped)
    {
        *s = t1;
        *t = s1;
    }
    else
    {
        *s = s1;
        *t = t1;
    }
    return q;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwComplex signum(const hwComplex &cplx)
{
    return cplx / cplx.Mag();
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency curMultElems(Currency &x, Currency &y)
{
    // Assumes both inputs are scalar, complex, or matrices
    // if x and y are matrices, assumes they are the same size
    if (x.IsScalar())
    {
        if (y.IsScalar())
            return Currency(x.Scalar() * y.Scalar());
        else if (y.IsComplex())
            return Currency(y.Complex() * x.Scalar());
        else if (y.IsMatrix())
        {
            hwMatrix mtx = *y.Matrix() * x.Scalar();
            return EvaluatorInterface::allocateMatrix(&mtx);
        }
    }
    else if (x.IsComplex())
    {
        if (y.IsScalar())
            return Currency(x.Complex() * y.Scalar());
        else if (y.IsComplex())
            return Currency(x.Complex() * y.Complex());
        else if (y.IsMatrix())
        {
            hwMatrix mtx = *y.Matrix() * x.Complex();
            return EvaluatorInterface::allocateMatrix(&mtx);
        }
    }
    else if (x.IsMatrix())
    {
        hwMatrix *result;
        if (y.IsScalar())
        {
            hwMatrix mtx = *x.Matrix() * y.Scalar();
            result = EvaluatorInterface::allocateMatrix(&mtx);
            return Currency(result);
        }
        else if (y.IsComplex())
        {
            hwMatrix mtx = *x.Matrix() * y.Complex();
            result = EvaluatorInterface::allocateMatrix(&mtx);
            return Currency(result);
        }
        else if (y.IsMatrix())
        {
            result = EvaluatorInterface::allocateMatrix();
            Currency resultCur(result);
            result->MultByElems(*(x.Matrix()), *(y.Matrix()));
            return resultCur;
        }
    }
    return Currency();
}

// sort-related methods
//------------------------------------------------------------------------------
// a sorted complex vector is sorted primarily based on magnitude,
// secondarily on the sign of the imaginary component (negative first),
// tertiarily on the real component (increasing order if the complex
// component is negative, decreasing otherwise)
//
// this is to format the union output appropriately
//------------------------------------------------------------------------------
int indexOf(const std::deque<hwComplex> &d, const hwComplex &val, bool reverseReturn, bool allowDuplicate)
{
    int size = (int)d.size(), minIndex = 0, maxIndex = size, midIndex = 0;
    while (minIndex <= maxIndex && minIndex < size)
    {
        midIndex = (minIndex + maxIndex) / 2;
        hwComplex curVal = d[midIndex];

        if (curVal == val)
            return (!allowDuplicate && reverseReturn ? -1 : midIndex);
        else if (complexLessThan(curVal, val))
            minIndex = midIndex + 1;
        else
            maxIndex = midIndex - 1;
    }

    return (allowDuplicate || reverseReturn ? minIndex : -1);
}
//------------------------------------------------------------------------------
// Assumes both matrices are row vectors of the same size
//------------------------------------------------------------------------------
int indexOf(const std::deque<hwMatrix> &d, hwMatrix &val, bool reverseReturn, bool allowDuplicate)
{
    int size = (int)d.size(), minIndex = 0, maxIndex = size, midIndex = 0;
    while (minIndex <= maxIndex && minIndex < size)
    {
        midIndex = (minIndex + maxIndex) / 2;
        hwMatrix curVal = d[midIndex];

        if (curVal == val)
            return (!allowDuplicate && reverseReturn ? -1 : midIndex);
        else if (rowVecLessThan(&curVal, &val))
            minIndex = midIndex + 1;
        else
            maxIndex = midIndex - 1;
    }
    return (allowDuplicate || reverseReturn ? minIndex : -1);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
int indexOf(const std::deque<std::string> &d, const std::string &val, bool reverseReturn, bool allowDuplicate)
{
    int size = (int)d.size(), minIndex = 0, maxIndex = size, midIndex = 0;
    while (minIndex <= maxIndex && minIndex < size)
    {
        midIndex = (minIndex + maxIndex) / 2;
        int comparison = d[midIndex].compare(val);

        if (!comparison)
            return (!allowDuplicate && reverseReturn ? -1 : midIndex);
        else if (comparison < 0)
            minIndex = midIndex + 1;
        else
            // comparison > 0
            maxIndex = midIndex - 1;
    }
    return (allowDuplicate || reverseReturn ? minIndex : -1);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool complexLessThan(const hwComplex &cplx1, const hwComplex &cplx2)
{
    return complexLessThanTol(cplx1, cplx2, 0.0);
}
//------------------------------------------------------------------------------
// Assumes both vectors are of the same size
//------------------------------------------------------------------------------
bool rowVecLessThan(const hwMatrix *mtx1, const hwMatrix *mtx2)
{
    if (mtx1->M() != 1 || mtx2->M() != 1)
        throw OML_Error(HW_ERROR_INPROWVECT);

    bool cplx1 = !mtx1->IsReal();
    bool cplx2 = !mtx2->IsReal();

    if (cplx1 || cplx2)
    {
        for (int i = 0; i < mtx1->Size(); i++)
        {
            const hwComplex c1 = cplx1 ? mtx1->z(i) : hwComplex((*mtx1)(i), 0.0);
            const hwComplex c2 = cplx2 ? mtx2->z(i) : hwComplex((*mtx2)(i), 0.0);
            if (c1 == c2)
                continue;
            if (complexLessThan(c1, c2))
                return true;
            return false;
        }
        return false;
    }
    else
    {
        for (int i = 0; i < mtx1->Size(); i++)
        {
            double d1 = (*mtx1)(i), d2 = (*mtx2)(i);
            if (d1 == d2)
                continue;
            if ((*mtx1)(i) < (*mtx2)(i))
                return true;
            return false;
        }
        return false;
    }
}

// string-related methods
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::string doEscapeSequences(const std::string& str)
{
    std::stringstream ss;
    std::string::const_iterator iter = str.cbegin();
    bool escape = false;
    while (iter != str.cend())
    {
        if (escape)
        {
            ss << escapeChar(*iter++);
            escape = false;
        }
        else 
        {
            if (*iter == '\\')
                escape = true;
            else
                ss << *iter;
            ++iter;
        }

    }
    if (escape)
        ss << '\\';

    return ss.str();
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
char escapeChar(char c)
{
    switch(c)
    {
        case '\\':
            return '\\';
        case '"':
            return '"';
        case '\'':
            return '\'';
        case '0':
            return '\0';
        case 'a':
            return '\a';
        case 'b':
            return '\b';
        case 'f':
            return '\f';
        case 'n':
            return '\n';
        case 'r':
            return '\r';
        case 't':
            return '\t';
        case 'v':
            return '\v';
        case '?':
            return '\?';
        default:
            throw OML_Error("Error: Invalid escape sequence: '\\"  + std::string(1, c) + "'");
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isstr(HML_CELLARRAY *cell)
{
    for (int i = 0; i < cell->Size(); i++)
        if (!(*cell)(i).IsString())
            return false;
    return true;
}
//------------------------------------------------------------------------------
// still throws error for complex matrix
//------------------------------------------------------------------------------
hwMatrix* tostr(EvaluatorInterface& eval, const hwMatrix *m, bool throwError)
{
    if (!m)
        return EvaluatorInterface::allocateMatrix();

    hwMatrix *newm = EvaluatorInterface::allocateMatrix(m->M(), m->N(), hwMatrix::REAL);
    try
    {
        if (m->IsRealData())
        {
            for (int i = 0; i < m->Size(); i++)
                (*newm)(i) = (double) BuiltInFuncsUtils::GetValidChar(eval, realval(m, i), throwError);
        }
        else
            throw OML_Error(HW_ERROR_NOTCOMPTOSTR);
    }
    catch (OML_Error& e)
    {
        delete newm;
        throw;
    }
    return newm;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwMatrix* tostr(EvaluatorInterface& eval, HML_CELLARRAY *cell, bool throwError)
{
    int m = 0, n = 0;
    std::vector<Currency> curVec;
    for (int i = 0; i < cell->Size(); i++)
    {
        Currency &elem = (*cell)(i);
        elem = toCurrencyStr(eval, elem, throwError, false);
        const hwMatrix *currentStr = elem.Matrix();
        n = max(currentStr->N(), n);
        m += max(1, currentStr->M());
        curVec.push_back(elem);
    }

    hwMatrix *str = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
    buildString(curVec, str);
    return str;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
HML_CELLARRAY* tocellstr(EvaluatorInterface& eval, HML_CELLARRAY *c, bool throwError)
{
    HML_CELLARRAY *strcell = EvaluatorInterface::allocateCellArray(c->M(), c->N());
    try
    {
        for (int i = 0; i < c->Size(); i++)
            (*strcell)(i) = toCurrencyStr(eval, (*c)(i), throwError, true);
    }
    catch (OML_Error& e)
    {
        delete strcell;
        throw;
    }
    return strcell;
}
//------------------------------------------------------------------------------
// still throws error for complex values
// returned currency and will be a string or cell array with all elements (and their elements) as strings
//------------------------------------------------------------------------------
Currency toCurrencyStr(EvaluatorInterface& eval, const Currency &c, bool throwError, bool keepcell)
{
    hwMatrix *toreturn;
    if (c.IsScalar())
        toreturn = EvaluatorInterface::allocateMatrix(1, 1, (double) BuiltInFuncsUtils::GetValidChar(eval, c.Scalar(), throwError));
    else if (c.IsMatrix())
        toreturn = tostr(eval, c.Matrix(), throwError);
    else if (c.IsCellArray())
        if (keepcell)
            return tocellstr(eval, c.CellArray(), throwError);
        else
            toreturn = tostr(eval, c.CellArray(), throwError);
    else if (c.IsString())
        toreturn = EvaluatorInterface::allocateMatrix((c.Matrix()));
    else
        throw OML_Error(HW_ERROR_INVINPSTRUCTANDCOMPNOTCONVSTR);

    return addStringMask(toreturn);
}
//------------------------------------------------------------------------------
// Converts to lower
//------------------------------------------------------------------------------
std::string convertToLower(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(), static_cast<int (*)(int)>(tolower));
    return str;
}
//------------------------------------------------------------------------------
// Conversts to lower
//------------------------------------------------------------------------------
hwMatrix* convertToLower(const hwMatrix *str)
{
    hwMatrix *newstr = EvaluatorInterface::allocateMatrix(str->M(), str->N(), hwMatrix::REAL);
    for (int i = 0; i < str->Size(); i++)
        (*newstr)(i) = tolower((int) (*str)(i));
    return newstr;
}
//------------------------------------------------------------------------------
// Converts to upper
//------------------------------------------------------------------------------
hwMatrix* convertToUpper(const hwMatrix *str)
{
    hwMatrix *newstr = EvaluatorInterface::allocateMatrix(str->M(), str->N(), hwMatrix::REAL);
    for (int i = 0; i < str->Size(); i++)
        (*newstr)(i) = toupper((int) (*str)(i));
    return newstr;
}
//------------------------------------------------------------------------------
// assumes all elements are strings or cell arrays of strings
//------------------------------------------------------------------------------
HML_CELLARRAY* convertToUpper(const HML_CELLARRAY *cell)
{
    HML_CELLARRAY *returnval = EvaluatorInterface::allocateCellArray(cell->M(), cell->N());
    for (int i = 0; i < cell->Size(); i++)
    {
        (*returnval)(i) = convertToUpper((*cell)(i));
    }
    return returnval;
}
//------------------------------------------------------------------------------
// assumes all elements are strings or cell arrays of strings
//------------------------------------------------------------------------------
HML_CELLARRAY* convertToLower(const HML_CELLARRAY *cell)
{
    HML_CELLARRAY *returnval = EvaluatorInterface::allocateCellArray(cell->M(), cell->N());
    for (int i = 0; i < cell->Size(); i++)
    {
        (*returnval)(i) = convertToLower((*cell)(i));
    }
    return returnval;
}
//-----------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency convertToUpper(const Currency &cur)
{
    if (cur.IsString())
        return addStringMask(convertToUpper(cur.Matrix()));
    else if (cur.IsCellArray())
        return convertToUpper(cur.CellArray());
    else
        return cur;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency convertToLower(const Currency &cur)
{
    if (cur.IsString())
        return addStringMask(convertToLower(cur.Matrix()));
    else if (cur.IsCellArray())
        return convertToLower(cur.CellArray());
    else
        return cur;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::vector<Currency> convertToLower(std::vector<Currency>::const_iterator start, std::vector<Currency>::const_iterator end)
{
    std::vector<Currency> ret;
    while (start != end)
        ret.push_back(convertToLower(*start++));
    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void cellAppendStr(EvaluatorInterface& eval, HML_CELLARRAY *outcell, const hwMatrix *newstr, bool throwError)
{
    Currency holder = addStringMask(newstr);

    if (newstr->Size())
    {
        for (int j = 0; j < outcell->Size(); j++)
        {
            Currency &toappend = (*outcell)(j);
            toappend = hconcat(eval, toappend, holder);
        }
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwMatrix* trimright(const hwMatrix *totrim)
{
    if (!totrim)
        return EvaluatorInterface::allocateMatrix();

    int n = totrim->N();
    int m = totrim->M();
    for (int i = n - 1; i > -1; i--)
    {
        int j;
        for (j = m - 1; j > -1; j--)
        {
            int val = (int) (*totrim)(j, i);
            if (val != ' ')
                break;
        }

        if (j != -1)
            break;

        n--;
    }

    hwMatrix *trimmed = EvaluatorInterface::allocateMatrix(m, n, hwMatrix::REAL);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            (*trimmed)(i, j) = (*totrim)(i, j);
        }
    }

    return trimmed;
}
//------------------------------------------------------------------------------
// returns whether the input was valid
//------------------------------------------------------------------------------
bool getSingularStringOrCell(Currency input, hwMatrix *&outstr, HML_CELLARRAY *&outcell, bool &usecell)
{
    if (input.IsString())
    {
        outstr = input.GetWritableMatrix();
        usecell = false;
    }
    else if (input.IsCellArray())
    {
        HML_CELLARRAY *temp = input.CellArray();
        if (temp->Size() == 1)
        {
            Currency c = (*temp)(0);
            if (c.IsString())
            {
                outstr = c.GetWritableMatrix();
                usecell = false;
            }
            else
            {
                return false;
            }
        }
        else
        {
            outcell = EvaluatorInterface::allocateCellArray(temp);
            usecell = true;
        }
    }
    else
    {
        return false;
    }
    return true;
}
//------------------------------------------------------------------------------
// every currency should be a string
// this just concatenates them
//------------------------------------------------------------------------------
void buildString(std::vector<Currency> vec, hwMatrix *str)
{
    int m = 0;
    for (int i = 0; i < vec.size(); i++)
    {
        const hwMatrix *elem = vec[i].Matrix();
        int j;
        for (j = 0; j < elem->M(); j++)
        {
            int k;
            for (k = 0; k < elem->N(); k++)
                (*str)(m, k) = (double) (*elem)(j, k);

            for (int l = k; l < str->N(); l++)
                (*str)(m, l) = ' ';

            m++;
        }
        if (!j) //&& str->M() != 0)
        {
            for (int l = 0; l < str->N(); l++)
                (*str)(m, l) = ' ';

            m++;
        }
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::string orderedStringVal(const Currency &cur)
{
    return orderedStringVal(cur.Matrix());
}
//------------------------------------------------------------------------------
// assumes all elements of mtx are within unsigned char range
//------------------------------------------------------------------------------
std::string orderedStringVal(const hwMatrix *strmtx)
{
    std::string returnval;
    for (int i = 0; i < strmtx->Size(); i++)
        returnval += (unsigned char) (*strmtx)(i);
    return returnval;
}

// date-related methods

//------------------------------------------------------------------------------
// assumes no leap year
//------------------------------------------------------------------------------
int getDaysInMonth(int month)
{

    switch (month)
    {
        case 1:
        case 3:
        case 5:
        case 7:
        case 8:
        case 10:
        case 12:
            return 31;
        case 2:
            return 28;
        case 4:
        case 6:
        case 9:
        case 11:
            return 30;
        default:
            throw OML_Error(HW_MATH_MSG_INTERNALERROR);
    }
}
//------------------------------------------------------------------------------
// assumes months is between 1 and 12 inclusively
//------------------------------------------------------------------------------
double daysInMonths(int monthsPassed, int year)
{
    double rv = 0.0;
    //int monthsPassed = (int) months;
    for (int i = 1; i <= monthsPassed; i++)
    {
        rv += getDaysInMonth(i);
    }
    // rv += getDaysInMonth(monthsPassed + 1) * (months - monthsPassed);
    if (monthsPassed >= 2 && isLeapYear(year))
        rv++;
    return rv;
}
//------------------------------------------------------------------------------
// Returns true if leap year
//------------------------------------------------------------------------------
bool isLeapYear(int year)
{
    return !year || (isint(year / 4.0) && (isint(year / 400.0) || !isint(year / 100.0)));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency deblankHelper(EvaluatorInterface& eval, const Currency& cur)
{
    if (cur.IsCellArray())
    {
        HML_CELLARRAY* cell = cur.CellArray();
        HML_CELLARRAY* result = EvaluatorInterface::allocateCellArray(cell->M(), cell->N());
        Currency out(result);
        for (int i = 0; i < cell->Size(); ++i)
        {
            (*result)(i) = deblankHelper(eval, (*cell)(i));
        }
        return out;
    }
    else if (cur.IsString())
    {
        const hwMatrix* str = cur.Matrix();
        int maxlen = 0;
        for (int i = 0; i < str->M(); ++i)
        {
            for (int j = maxlen; j < str->N(); ++j)
            {
                int ch = (int) (*str)(i,j);
                if (ch && !isspace(ch))
                    maxlen = j + 1;
            }
        }

        if (maxlen)
        {
            hwMatrix* result = EvaluatorInterface::allocateMatrix(str->M(), maxlen, hwMatrix::REAL);
            for (int i = 0; i < str->M(); ++i)
            {
                for (int j = 0; j < maxlen; ++j)
                {
                    (*result)(i,j) = (*str)(i,j);
                }
            }
            return addStringMask(result);
        }
        else
        {
            return addStringMask(EvaluatorInterface::allocateMatrix());
        }
    }
    else
        throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void removeFields(EvaluatorInterface& eval, StructData* strct, const hwMatrix* fields)
{
    for (int i = 0; i < fields->M(); ++i)
    {
        std::string field = orderedStringVal(readRow(eval, fields, i));
        if (strct->Contains(field))
            strct->removeField(field);
        else
            throw OML_Error("Error: struct does not contain field: " + field);
    }
}
//------------------------------------------------------------------------------
// if n is 0, a full comparison of the strings is used
//------------------------------------------------------------------------------
bool _strncmpHelper(const hwMatrix *s1, const hwMatrix *s2, int n)
{
    if (n)
    {
        if (s1->Size() < n || s2->Size() < n)
            return false;

        for (int i = 0; i < n; i++)
        {
            if ((*s1)(i) != (*s2)(i))
                return false;
        }
        return true;
    }
    return *s1 == *s2;
}
//------------------------------------------------------------------------------
// if n is 0, does a full comparison
//------------------------------------------------------------------------------
bool dostrcmp(const std::vector<Currency>& inputs, std::vector<Currency>& outputs, int n)
{
    const Currency &input1 = inputs[0];
    const Currency &input2 = inputs[1];

    bool cellInput1 = false;
    bool cellInput2 = false;
    HML_CELLARRAY *cell1, *cell2;
    hwMatrix *str1, *str2;
    if (!getSingularStringOrCell(input1, str1, cell1, cellInput1) || !getSingularStringOrCell(input2, str2, cell2, cellInput2))
    {
        outputs.push_back(getFalse());
        return true;
    }

    if (cellInput1)
    {
        if (cellInput2)
        {
            if (sameSize(cell1, cell2))
            {
                hwMatrix *outmtx = EvaluatorInterface::allocateMatrix(cell1->M(), cell1->N(), hwMatrix::REAL);
                for (int i = 0; i < outmtx->Size(); i++)
                {
                    Currency cur1 = (*cell1)(i);
                    Currency cur2 = (*cell2)(i);
                    if (cur1.IsString() && cur2.IsString())
                    {
                        (*outmtx)(i) = _strncmpHelper(cur1.Matrix(), cur2.Matrix(), n);
                    }
                    else
                    {
                        (*outmtx)(i) = 0.0;
                    }
                }

                Currency out(outmtx);
                out.SetMask(Currency::MASK_LOGICAL);
                outputs.push_back(out);
            }
            else
                throw OML_Error(HW_ERROR_INPMUSTSAMESIZE);
        }
        else
        {
            hwMatrix *outmtx = EvaluatorInterface::allocateMatrix(cell1->M(), cell1->N(), hwMatrix::REAL);
            for (int i = 0; i < outmtx->Size(); i++)
            {
                Currency cur1 = (*cell1)(i);
                if (cur1.IsString())
                {
                    (*outmtx)(i) = _strncmpHelper(cur1.Matrix(), str2, n);
                }
                else
                {
                    (*outmtx)(i) = 0.0;
                }
            }

            Currency out(outmtx);
            out.SetMask(Currency::MASK_LOGICAL);
            outputs.push_back(out);
        }
    }
    else
    {
        if (cellInput2)
        {
            hwMatrix *outmtx = EvaluatorInterface::allocateMatrix(cell2->M(), cell2->N(), hwMatrix::REAL);
            for (int i = 0; i < outmtx->Size(); i++)
            {
                Currency cur2 = (*cell2)(i);
                if (cur2.IsString())
                {
                    (*outmtx)(i) = _strncmpHelper(str1, cur2.Matrix(), n);
                }
                else
                {
                    (*outmtx)(i) = 0.0;
                }
            }

            Currency out(outmtx);
            out.SetMask(Currency::MASK_LOGICAL);
            outputs.push_back(out);
        }
        else
        {
            Currency out(_strncmpHelper(str1, str2, n));
            out.SetMask(Currency::MASK_LOGICAL);
            outputs.push_back(out);
        }
    }
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency tripleCurrencyFunc(EvaluatorInterface& eval, const Currency &cur1, const Currency &cur2, const Currency &cur3, void* overlap,
    Currency (*func)(EvaluatorInterface&, const Currency*, const Currency*, const Currency*, void*))
{
    HML_CELLARRAY *searchCell, *patternCell, *replaceCell;
    bool usecell1 = cur1.IsCellArray();
    bool usecell2 = cur2.IsCellArray();
    bool usecell3 = cur3.IsCellArray();
    bool singular1, singular2, singular3;
    const Currency *firstSearch, *firstPattern, *firstReplace; // only needed if singular
    firstSearch = firstPattern = firstReplace = nullptr;

    if (usecell1)
    {
        searchCell = cur1.CellArray();
        if (singular1 = (searchCell->Size() == 1))
            firstSearch = &(*searchCell)(0);
    }
    else
    {
        singular1 = true;
        firstSearch = &cur1;
    }

    if (usecell2)
    {
        patternCell = cur2.CellArray();
        if (singular2 = (patternCell->Size() == 1))
            firstPattern = &(*patternCell)(0);
    }
    else
    {
        singular2 = true;
        firstPattern = &cur2;
    }

    if (usecell3)
    {
        replaceCell = cur3.CellArray();
        if (singular3 = (replaceCell->Size() == 1))
            firstReplace = &(*replaceCell)(0);
    }
    else
    {
        singular3 = true;
        firstReplace = &cur3;
    }

    // check dimensions
    if (singular1)
    {
        if (!(singular2 || singular3 || sameSize(patternCell, replaceCell)))
            throw OML_Error(HW_ERROR_CELLARRAYSSAMESIZE);
    }
    else if (singular2)
    {
        if (!(singular1 || singular3 || sameSize(searchCell, replaceCell)))
            throw OML_Error(HW_ERROR_CELLARRAYSSAMESIZE);
    }
    else if (singular3)
    {
        if (!(singular1 || singular2 || sameSize(searchCell, patternCell)))
            throw OML_Error(HW_ERROR_CELLARRAYSSAMESIZE);
    }
    else
    {
        if (!(sameSize(searchCell, patternCell) && sameSize(patternCell, replaceCell)))
            throw OML_Error(HW_ERROR_CELLARRAYSSAMESIZE);
    }

    HML_CELLARRAY *output;

    int singularCode = singular1 + singular2 * 2 + singular3 * 4; // possible results are 0-7

    // for each possibility, loop through as appropriate
    switch (singularCode)
    {
        case 0:
            // no singular inputs
            output = EvaluatorInterface::allocateCellArray(searchCell->M(), searchCell->N());
            for (int i = 0; i < output->Size(); i++)
                (*output)(i) = tripleCurrencyFunc(eval, (*searchCell)(i), (*patternCell)(i), (*replaceCell)(i), overlap, func);
            return output;
        case 1:
            // only first input is singular
            output = EvaluatorInterface::allocateCellArray(patternCell->M(), patternCell->N());
            for (int i = 0; i < output->Size(); i++)
                (*output)(i) = tripleCurrencyFunc(eval, *firstSearch, (*patternCell)(i), (*replaceCell)(i), overlap, func);
            return output;
        case 2:
            // only second input is singular
            output = EvaluatorInterface::allocateCellArray(searchCell->M(), searchCell->N());
            for (int i = 0; i < output->Size(); i++)
                (*output)(i) = tripleCurrencyFunc(eval, (*searchCell)(i), *firstPattern, (*replaceCell)(i), overlap, func);
            return output;
        case 3:
            // first and second inputs are singular
            output = EvaluatorInterface::allocateCellArray(replaceCell->M(), replaceCell->N());
            for (int i = 0; i < output->Size(); i++)
                (*output)(i) = tripleCurrencyFunc(eval, *firstSearch, *firstPattern, (*replaceCell)(i), overlap, func);
            return output;
        case 4:
            // only third input is singular
            output = EvaluatorInterface::allocateCellArray(searchCell->M(), searchCell->N());
            for (int i = 0; i < output->Size(); i++)
                (*output)(i) = tripleCurrencyFunc(eval, (*searchCell)(i), (*patternCell)(i), *firstReplace, overlap, func);
            return output;
        case 5:
            // first and third inputs are singular
            output = EvaluatorInterface::allocateCellArray(patternCell->M(), patternCell->N());
            for (int i = 0; i < output->Size(); i++)
                (*output)(i) = tripleCurrencyFunc(eval, *firstSearch, (*patternCell)(i), *firstReplace, overlap, func);
            return output;
        case 6:
            // second and third inputs are singular
            output = EvaluatorInterface::allocateCellArray(searchCell->M(), searchCell->N());
            for (int i = 0; i < output->Size(); i++)
                (*output)(i) = tripleCurrencyFunc(eval, (*searchCell)(i), *firstPattern, *firstReplace, overlap, func);
            return output;
        case 7:
            // all inputs are singular
            if (!(usecell1 || usecell2 || usecell3))
                return (*func)(eval, &cur1, &cur2, &cur3, overlap);

            output = EvaluatorInterface::allocateCellArray(1, 1);
            (*output)(0) = tripleCurrencyFunc(eval, *firstSearch, *firstPattern, *firstReplace, overlap, func);
            return output;
        default:
            // shouldn't ever be reached
            throw OML_Error(HW_MATH_MSG_INTERNALERROR);
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency doStrRep(EvaluatorInterface& eval, const Currency *searchcur, const Currency *patcur, const Currency *repcur, void *dooverlap)
{
    bool overlap = *static_cast<bool*>(dooverlap);
    if (!(searchcur->IsString() && patcur->IsString() && repcur->IsString()))
        throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

    std::string searchStr = orderedStringVal(*searchcur);
    std::string pattern = orderedStringVal(*patcur);
    std::string replacement = orderedStringVal(*repcur);

    if (!pattern.length())
        return Currency(std::string());

    if (overlap)
    {
        std::vector<int> indices, lengths;
        size_t start = -1;
        int offset = 0;
        int offsetincr = (int)replacement.length() - (int)pattern.length();
        while ((start = searchStr.find(pattern, start + 1)) != std::string::npos)
        {
            if (indices.size())
            {
                int prev = indices.back();
                int endlast = prev + (int)replacement.length();
                int startIndex = (int)start + offset;
                if (startIndex < endlast)
                {
                    lengths.back() -= endlast - startIndex;
                    offset += endlast - startIndex;
                }

                indices.push_back((int)start + offset);
                lengths.push_back((int)pattern.length());
                offset += offsetincr;
            }
            else
            {
                indices.push_back((int)start);
                lengths.push_back((int)pattern.length());
                offset += offsetincr;
            }
        }

        if (indices.size())
        {
            for (int i = 0; i < indices.size(); i++)
            {
                searchStr.replace(indices[i], lengths[i], replacement);
            }
            return Currency(searchStr);
        }
        else
        {
            return *searchcur;
        }
    }
    else
    {
        size_t start = 0;
        bool madeReplacement = false;
        while ((start = searchStr.find(pattern, start)) != std::string::npos)
        {
            madeReplacement = true;
            searchStr.replace(start, pattern.length(), replacement);
            start += 1 + replacement.length() - pattern.length();
        }

        if (madeReplacement)
        {
            return Currency(searchStr);
        }
        else
        {
            return *searchcur;
        }
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void dostrcat(EvaluatorInterface& eval, Currency &out, HML_CELLARRAY *cell)
{
    if (out.IsCellArray())
    {
        HML_CELLARRAY *outcell = out.CellArray();
        if (sameSize(cell, outcell))
        {
            for (int j = 0; j < cell->Size(); j++)
            {
                Currency &newc = (*cell)(j);
                Currency &oldc = (*outcell)(j);
                oldc = hconcat(eval, oldc, newc);
            }
        }
        else
            throw OML_Error(HW_ERROR_CELLARRAYSSAMESIZE);
    }
    else
    {
        const hwMatrix *strout = out.Matrix();
        HML_CELLARRAY *outcell = EvaluatorInterface::allocateCellArray(cell->M(), cell->N());
        Currency temp(outcell);
        for (int j = 0; j < cell->Size(); j++)
        {
            Currency &toappend = (*cell)(j);
            if (toappend.IsCellArray())
            {
                if (strout->Size())
                    (*outcell)(j) = hconcat(eval, out, toappend.CellArray());
                else
                    (*outcell)(j) = EvaluatorInterface::allocateCellArray(toappend.CellArray());
            }
            else
                (*outcell)(j) = dostrcat(eval, strout, toappend.Matrix());
        }
        out = temp;
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency dostrcat(EvaluatorInterface& eval, const hwMatrix *first, const hwMatrix *second)
{
    int firstm = first->M(), secondm = second->M();
    if (firstm == secondm || !first->Size() || !second->Size())
        return addStringMask(hconcat(first, second));

    if (firstm == 1)
    {
        hwMatrix *toreturn = EvaluatorInterface::allocateMatrix(secondm, first->N() + second->N(), hwMatrix::REAL);
        Currency returnCur(toreturn);
        hwMatrix row;
        for (int i = 0; i < secondm; i++)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, second->ReadRow(i, row));
            hwMatrix *newrow = hconcat(first, &row);
            Currency newrowCur(newrow);
            BuiltInFuncsUtils::CheckMathStatus(eval, toreturn->WriteRow(i, *newrow));
        }
        returnCur.SetMask(Currency::MASK_STRING);
        return returnCur;
    }
    else if (secondm == 1)
    {
        hwMatrix *toreturn = EvaluatorInterface::allocateMatrix(firstm, first->N() + second->N(), hwMatrix::REAL);
        Currency returnCur(toreturn);
        hwMatrix row;
        for (int i = 0; i < firstm; i++)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, first->ReadRow(i, row));
            hwMatrix *newrow = hconcat(&row, second);
            Currency newrowCur(newrow);
            BuiltInFuncsUtils::CheckMathStatus(eval, toreturn->WriteRow(i, *newrow));
        }
        returnCur.SetMask(Currency::MASK_STRING);
        return returnCur;
    }
    else
        throw OML_Error(HW_ERROR_STRSAMENUMROWCONCAT);
}

void dostrcat(EvaluatorInterface& eval, Currency &out, const hwMatrix *new_str)
{
    if (out.IsCellArray())
	{		
		hwMatrix* newstr = EvaluatorInterface::allocateMatrix(new_str);
        cellAppendStr(eval, out.CellArray(), newstr, true);
		//delete newstr; // uncommenting this causes some tests to fail
	}
    else
	{
        out = dostrcat(eval, out.Matrix(), new_str);
	}
}
//------------------------------------------------------------------------------
// Returns true if successful in reading till end of newline, EOF or a 
// specified number of characters
//------------------------------------------------------------------------------
bool dofgets(EvaluatorInterface&          eval, 
             const std::vector<Currency>& inputs, 
             std::string&                 readline)
{
    size_t nargin = inputs.empty() ? 0 : inputs.size();
    if (nargin < 1 || nargin > 2) throw OML_Error(OML_ERR_NUMARGIN);

    readline = "";

    // Get the file name or file id
    const Currency &input1 = inputs[0];  
    int fid = getFileFromInput(eval, input1);

    checkFileIndex(eval, fid, true);
    std::FILE *file = eval.GetFile(fid);

    if (nargin > 1)  // Reads specified number of characters
    {
        int len = (int) inputs[1].Scalar();

        if (!IsInteger(inputs[1].Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_FILEID);

        if (len < 0) throw OML_Error(HW_ERROR_LENNOTNEG);

        if (!len)
        {
            if (feof(file))
                return false;
            return true;
        }

        char *data = new char[len + 1];
        memset(data, 0, sizeof(char)*(len+1));
        if (fgets(data, len + 1, file) == nullptr)
        {
            delete [] data;
            return false;
        }
        readline = std::string(data);

        delete [] data;
        return true;
    }
   
    // Read file till the first '\r' or '\n' or '\r\n' is encountered
    while (1)
    {
        int c = fgetc (file);
        if (c == EOF) break;  

        if (c == '\n')                     // First newline character
        {
            readline += c;
            break;   
        }

        else if (c == '\r')               // First carriage return
        {
            readline += c;
            int nextch = fgetc (file);
            if (nextch == EOF) break;
            else if (nextch == '\n')     // '\r\n'
            {
                readline += nextch;
                break;
            }
            else                        // Rewind here as we want only '\r'
            {
                fseek(file, static_cast<long>(-1), SEEK_CUR);
                break;
            }
        }
        else
            readline += c;
    }
    return (!readline.empty());
}
//------------------------------------------------------------------------------
// Helper function for oml_celldisp
//------------------------------------------------------------------------------
void celldisp(EvaluatorInterface& eval, 
              HML_CELLARRAY*      cell, 
              const std::string&  cellname)
{
    if (cell->IsVector())
    {
        int cellsize = cell->Size();
        for (int i = 0; i < cellsize; ++i)
        {
            std::string index("{");
            index += std::to_string(static_cast<long long>(i + 1)) + "}";

            Currency &elem = (*cell)(i);
            elem.DispOutput();
            if (elem.IsCellArray())
                celldisp(eval, elem.CellArray(), cellname + index);
            else
            {
                Currency label (cellname + index + " =");
                label.DispOutput();
                eval.PrintResult(label);
                eval.PrintResult(elem);
            }
        }
    }
    else
    {
        int numcols = cell->N();
        int numrows = cell->M();
        for (int j = 0; j < numcols; ++j)
        {
            for (int i = 0; i < numrows; ++i)
            {
                std::string index("{");
                index += std::to_string(static_cast<long long>(i + 1)) + ",";
                index += std::to_string(static_cast<long long>(j + 1)) + "}";

                Currency &elem = (*cell)(i, j);
                elem.DispOutput();
                if (elem.IsCellArray())
                    celldisp(eval, elem.CellArray(), cellname + index);
                else
                {
                    Currency label (cellname + index + " =");
                    label.DispOutput();
                    eval.PrintResult(label);
                    eval.PrintResult(elem);
                }
            }
        }
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency doScalarToLogical(const Currency &cur, double (*checker)(double))
{
    if (cur.IsScalar())
    {
        Currency out = (*checker)(cur.Scalar());
        out.SetMask(Currency::MASK_LOGICAL);
        return out;
    }
    else if (cur.IsMatrix() || cur.IsString())
    {
        const hwMatrix *mtx = cur.Matrix();
        if (mtx->IsRealData())
        {
            hwMatrix *m = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);

            for (int i = 0; i < mtx->Size(); i++)
            {
                (*m)(i) = (*checker)(realval(mtx, i));
            }

            Currency out(m);
            out.SetMask(Currency::MASK_LOGICAL);
            return out;
        }
        else
            throw OML_Error(HW_ERROR_INPUTREALSTR);
    }
    else if (cur.IsCellArray())
    {
        HML_CELLARRAY *cell = cur.CellArray();
        HML_CELLARRAY *outcell = EvaluatorInterface::allocateCellArray(cell->M(), cell->N());
        Currency out(outcell);

        for (int i = 0; i < outcell->Size(); i++)
        {
            (*outcell)(i) = doScalarToLogical((*cell)(i), checker);
        }

        out.SetMask(Currency::MASK_LOGICAL);
        return out;
    }
    else
        throw OML_Error(HW_ERROR_INPUTREALSTR);
}
//------------------------------------------------------------------------------
// Changes current working directory
//------------------------------------------------------------------------------
void cd(std::string &dir, EvaluatorInterface &eval)
{
    checkIsDirectory(eval, dir, true);
#ifdef OS_WIN
    if (!SetCurrentDirectory((LPCSTR) dir.c_str()))
#else
    if (chdir(dir.c_str()))
#endif
    {
        throw OML_Error(HW_ERROR_PROBCHANGCURDIR);
    }
    eval.ResetFuncSearchCache();
}
//------------------------------------------------------------------------------
// Returns true if empty
//------------------------------------------------------------------------------
bool isempty(const Currency &input)
{
    if (input.IsScalar() || input.IsComplex())
    {
        return false;
    }
    else if (input.IsMatrix() || input.IsString())
    {
        if (input.Matrix())
        return input.Matrix()->IsEmpty();
        else
            return true;
    }
    else if (input.IsCellArray())
    {
        return input.CellArray()->IsEmpty();
    }
    else if (input.IsStruct())
    {
        StructData *sd = input.Struct();
        return !(sd->M() * sd->N());
    }
    else
        throw OML_Error(HW_MATH_MSG_INVALIDINPUT);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::string strrep(EvaluatorInterface& eval, const std::string &search, const std::string &pattern, const std::string &replace, bool overlap)
{
    Currency searchcur(search), patterncur(pattern), replacecur(replace);
    return doStrRep(eval, &searchcur, &patterncur, &replacecur, &overlap).StringVal();
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::vector<Currency> mtxFun(EvaluatorInterface& eval, const std::vector<Currency>& inputs, int numOutputs,
    std::vector<Currency> (*method)(EvaluatorInterface&, const std::vector<Currency>&), const std::vector<Currency>& extras)
{
    size_t nargin = inputs.size();

    // get dimensions of inputs (and therefore outputs)
    const Currency& input1 = inputs[0];
    bool complex = false;
    int m, n;
    m = n = -1;
    for (size_t i = 0; i < nargin; ++i)
    {
        const Currency& cur = inputs[i];
        if (cur.IsComplex())
            complex = true;
        else if (cur.IsMatrix())
        {
            const hwMatrix* inmtx = cur.Matrix();

            if (!inmtx->IsReal())
                complex = true;

            if (m == -1)
            {
                if (inmtx->Size() != 1)
                {
                    m = inmtx->M();
                    n = inmtx->N();
                }
            }
            else if (complex)
                break;
        }
    }

    if (m == -1)
        m = n = 1;

    // initialize output matrices
    std::vector<Currency> outvec;
    for (int i = 0; i < numOutputs; ++i)
        outvec.push_back(EvaluatorInterface::allocateMatrix(m, n, complex ? hwMatrix::COMPLEX : hwMatrix::REAL));

    // set up args
    std::vector<Currency> args;
    for (size_t i = 0; i < nargin; ++i)
    {
        const Currency& in = inputs[i];
        if (in.IsScalar())
        {
            if (complex)
                args.push_back(EvaluatorInterface::allocateMatrix(m, n, hwComplex(in.Scalar(), 0.0)));
            else
                args.push_back(EvaluatorInterface::allocateMatrix(m, n, in.Scalar()));
        }
        else if (in.IsComplex())
        {
            args.push_back(EvaluatorInterface::allocateMatrix(m, n, in.Complex()));
        }
        else if (in.IsMatrix())
        {
            hwMatrix* mtx = EvaluatorInterface::allocateMatrix(in.Matrix());
            if (complex && mtx->IsReal())
                BuiltInFuncsUtils::CheckMathStatus(eval, mtx->MakeComplex());

            if (mtx->Size() == 1)
            {
                hwMatrix* newmtx = EvaluatorInterface::allocateMatrix(m, n, (*mtx)(0));
                args.push_back(newmtx);
            }
            else if (mtx->M() == m && mtx->N() == n)
                args.push_back(mtx);
            else
                throw OML_Error(HW_ERROR_MATRIXDIM);
        }
        else
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    // call method
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            // set up individual arguments
            std::vector<Currency> newinputs(extras);

            if (complex)
            {
                for (size_t k = 0; k < nargin; ++k)
                    newinputs.push_back(args[k].Matrix()->z(i, j));
            }
            else
            {
                for (size_t k = 0; k < nargin; ++k)
                    newinputs.push_back((*args[k].Matrix())(i, j));
            }

            // store results
            std::vector<Currency> results = (*method)(eval, newinputs);
            for (int k = 0; k < results.size(); ++k)
            {
                if (results[k].IsScalar())
                    outvec[k].GetWritableMatrix()->SetElement(i, j, results[k].Scalar());
                else
                    outvec[k].GetWritableMatrix()->SetElement(i, j, results[k].Complex());
            }
        }
    }

    return outvec;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool twoMatrixCaller(EvaluatorInterface& eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, hwMathStatus(hwMatrix::*func)(const hwMatrix&, const hwMatrix&))
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_VARIABLE);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_VARIABLE);

    const hwMatrix* mtx1 = inputs[0].ConvertToMatrix();
    const hwMatrix* mtx2 = inputs[1].ConvertToMatrix();
    hwMatrix* result = EvaluatorInterface::allocateMatrix();
    Currency out(result);
    BuiltInFuncsUtils::CheckMathStatus(eval, (result->*func)(*mtx1, *mtx2));
    outputs.push_back(out);
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool oneMatrixCaller(EvaluatorInterface& eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs, hwMathStatus(hwMatrix::*func)(const hwMatrix&))
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_VARIABLE);

    const hwMatrix* mtx = inputs[0].ConvertToMatrix();
    hwMatrix* result = EvaluatorInterface::allocateMatrix();
    Currency out(result);
    BuiltInFuncsUtils::CheckMathStatus(eval, (result->*func)(*mtx));
    outputs.push_back(out);
    return true;
}

// general common helpers
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void getMatchingSizes(const hwMatrix* m1, const hwMatrix* m2, int* m, int* n)
{
    *m = *n = -1;
    if (m1->Size() == 1)
    {
        if (m2->Size() != 1)
        {
            *m = m2->M();
            *n = m2->N();
        }
    }
    else
    {
        if (m2->Size() != 1 && !sameSize(m1, m2))
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2, OML_VAR_DIMS);

        *m = m1->M();
        *n = m1->N();
    }
}
//------------------------------------------------------------------------------
// Returns true if input is a variable name
//------------------------------------------------------------------------------
bool isvarname(const std::string& str)
{
    if (!str.length() || isdigit(str[0]))
        return false;
    return str.find_first_of(HW_BAD_VARNAME_CHARS) == std::string::npos;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void sleep(double sec)
{
#ifdef OS_WIN
        Sleep((DWORD)sec * 1000); // sleep specified amount of time
#else
        usleep(sec * 1000000);
#endif
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwMatrix* makeRealCopy(const hwMatrix* cplxMtx)
{
    hwMatrix *copy = EvaluatorInterface::allocateMatrix(cplxMtx->M(), cplxMtx->N(), hwMatrix::REAL);
    for (int i = 0; i < cplxMtx->Size(); i++)
        (*copy)(i) = cplxMtx->z(i).Real();
    return copy;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isActuallyComplex(const Currency& cur)
{
    if (cur.IsMatrix())
        return !cur.Matrix()->IsRealData();
    else if (cur.IsNDMatrix())
        return !cur.MatrixN()->IsRealData();
    else if (cur.IsComplex())
        return cur.Complex().Imag() != 0.0;
    return false;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isActuallyReal(const Currency& cur)
{
    if (cur.IsScalar() || cur.IsString())
        return true;
    else if (cur.IsMatrix() && cur.Matrix())
        return cur.Matrix()->IsRealData();
    else if (cur.IsNDMatrix())
        return cur.MatrixN()->IsRealData();
    else if (cur.IsComplex())
        return cur.Complex().Imag() == 0.0;
    return false;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::string readString(const Currency& str, int index)
{
    if (str.IsString())
    {
        std::stringstream ss;
        const hwMatrix* mat = str.Matrix();

		if (index < mat->M())
		{
			for (int i = 0; i < mat->N(); ++i)
				ss << (unsigned char) (*mat)(index, i);

			return ss.str();
		}
    }

    return std::string();
}
//------------------------------------------------------------------------------
// assumes input is a matrix or string
//------------------------------------------------------------------------------
Currency readRow(EvaluatorInterface& eval, const Currency &input, int index)
{
    Currency out = readRow(eval, input.Matrix(), index);
    out.SetMask(input.GetMask());
    return out;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwMatrix* readRow(EvaluatorInterface& eval, const hwMatrix *mtx, int index)
{
    const hwMatrix* loc   = mtx;
    bool            alloc = false;

    if (!loc)
    {
        alloc = true;
        loc = EvaluatorInterface::allocateMatrix();
    }

    std::unique_ptr<hwMatrix>row(EvaluatorInterface::allocateMatrix(1, loc->N(), loc->Type()));

    if (index < loc->M())
    {
        hwMathStatus stat = loc->ReadRow(index, *row);
        if (alloc && !stat.IsOk())
        {
            delete loc;
            loc = nullptr;
        }

        BuiltInFuncsUtils::CheckMathStatus(eval, stat);
    }

    if (alloc && loc)
    {
        delete loc;
        loc = nullptr;
    }

    return row.release();
}
//------------------------------------------------------------------------------
// Returns true if input is a directory
//------------------------------------------------------------------------------
bool checkIsDirectory(EvaluatorInterface& eval, std::string &str, bool throwError)
{
    std::string errmsg;
    if (isDirectory(str, &errmsg))
        return true;

    if (throwError)
        throw OML_Error("Error: " + errmsg);
    else
        BuiltInFuncsUtils::SetWarning(eval, "Warning: " + errmsg);
    return false;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isDirectory(std::string& str, std::string* errmsg)
{
    BuiltInFuncsUtils::StripTrailingSlash(str);
    struct stat st;
#if OS_WIN
    // handle root drives
    if (str.length() == 2 && isalpha(str[0]) && str[1] == ':')
        str += '\\';
#endif
    if (stat(str.c_str(), &st))
    {
        if (errmsg)
        {
            if (errno == ENOENT || errno == ENOTDIR)
                *errmsg = "no such directory exists";
            else
                *errmsg = "problem getting file info";
        }
    }
    else
    {
        if (S_ISDIR(st.st_mode))
            return true;
        else if (errmsg)
            *errmsg = str + " is not a directory";
    }
    return false;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::string getAbsolutePath(EvaluatorInterface& eval, const Currency &cur)
{
    if (cur.IsString())
        return BuiltInFuncsUtils::GetAbsolutePath(readString(cur));
    throw OML_Error(HW_ERROR_INPSTRMUSTFILEDIR);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::string sprintf(EvaluatorInterface& eval, std::vector<Currency>::const_iterator iter, const std::vector<Currency>::const_iterator enditer)
{
    if (iter == enditer)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency input1 = *iter;

    if (!input1.IsString())
        throw OML_Error(HW_ERROR_INPUTSTRING);

    std::string tmplt = readString(input1);
    return sprintf(eval, tmplt, iter + 1, enditer);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::string sprintf(EvaluatorInterface& eval, const std::string &origTemplate, std::vector<Currency>::const_iterator rawiter, const std::vector<Currency>::const_iterator rawenditer)
{
    std::string tmplt = doEscapeSequences(origTemplate);

    struct FormatSegment
    {
        std::string str;
        FormatType ft;
        int numReplacements;
        FormatSegment(std::string segstr, FormatType formatType, int replacements)
            : str(segstr), ft(formatType), numReplacements(replacements) {}
    };

    // separate template into segments with a maximum of 1 replacement string each
    std::vector<FormatSegment> segments;
    std::string stopatChars(HW_VALID_FORMATCHARS);

    size_t currentPerc = -1;
    size_t prevPerc = 0;
    do
    {
        std::string seg;
        currentPerc = tmplt.find('%', currentPerc + 1);
        size_t stopat = tmplt.find_first_of(stopatChars, prevPerc + 1);
        if (segments.size() && currentPerc == stopat)
            currentPerc = tmplt.find('%', stopat + 1);

        if (currentPerc == std::string::npos)
            seg = tmplt.substr(prevPerc);
        else
            seg = tmplt.substr(prevPerc, currentPerc - prevPerc);

        if (segments.size())
        {
            if (stopat == std::string::npos)
                throw OML_Error(HW_ERROR_INVALIDFORMAT);

	// temporary situation cause by lack of std::regex support on gcc 4.7.x
#ifdef OS_WIN
            // verify syntax for sprintf
            std::string temp = "% *[-+0# ]? *(\\*|([0-9]*))?(\\.(\\*|([0-9]*)))?(([lh]?[diouxXcs])|(L?[feEgG])|%)(.|\r|\n)*";
            std::regex pattern(temp);
            if (!std::regex_match(seg.cbegin(), seg.cend(), pattern))
                throw OML_Error(HW_ERROR_INVALIDFORMAT);
#endif

            // find number of extra replacements to be made in this segment (number of *'s)
            int count = (int) (std::count(tmplt.begin() + prevPerc, tmplt.begin() + stopat, '*'));

            switch (tmplt[stopat])
            {
                case '%':
                    segments.push_back(FormatSegment(seg, Percent, count));
                    break;
                case 'd':
                case 'i':
                case 'u':
                case 'c':
                case 'x':
                case 'X':
                case 'o':
                case 'O':
                    segments.push_back(FormatSegment(seg, Truncate, count));
                    break;
                case 's':
                    segments.push_back(FormatSegment(seg, String, count));
                    break;
                default:
                    segments.push_back(FormatSegment(seg, NoTruncate, count));
                    break;
            }
        }
        else
        {
            //segments.push_back(FormatSegment(seg, false, false, 0));
            // treat first segment as percent so no values are passed
            segments.push_back(FormatSegment(seg, Percent, 0));
        }
        prevPerc = currentPerc;
    }
    while (currentPerc != std::string::npos);

    if (segments.size() == 0)
    {
        return std::string();
    }
    else if (segments.size() == 1)
    {
        std::string seg = segments[0].str;
        return dosprintf(seg.length(), seg.c_str());
    }

    // Get a vector of valid inputs - don't print empty cells, matrices, structs
    std::vector<Currency> validInputs;
    for (; rawiter != rawenditer; rawiter++)
    {
        Currency inp = *rawiter;
        if (hasElements(inp))
            validInputs.push_back(inp);
    }

    std::string result;
    int inputindex = 0;
    FormatSegment fs = segments[0];
    size_t currentSegment = 0;
    std::vector<Currency>::const_iterator iter = validInputs.cbegin(), enditer = validInputs.cend();

    do
    {
        if (!currentSegment)
        {
            result += fs.str;
        }
        else
        {
            int numrep = fs.numReplacements;
            if (numrep > 0)
            {
                long long firstReplacement = (long long) getDoubleForSprintf(iter, &inputindex);
                if (iter == enditer && fs.ft != Percent)
                    break;

                if (numrep > 1)
                {
                    long long secondReplacement = (long long) getDoubleForSprintf(iter, &inputindex);
                    if (iter == enditer && fs.ft != Percent)
                        break;

                    // call dosprintf 2 extra replacement args
                    result += dosprintf(eval, fs.str, &inputindex, fs.ft, iter, firstReplacement, secondReplacement);
                }
                else
                {
                    // call dosprintf 1 extra replacement arg
                    result += dosprintf(eval, fs.str, &inputindex, fs.ft, iter, firstReplacement);
                }
            }
            else
            {
                // call dosprintf no extra replacement args
                result += dosprintf(eval, fs.str, &inputindex, fs.ft, iter);
            }
        }

        if (currentSegment > segments.size() - 2)
        {
            if (iter == enditer || (iter == validInputs.cbegin() && !inputindex))
                break;
            currentSegment = 0;
        }
        else
            currentSegment++;

        fs = segments[currentSegment];
    }
    while (iter != enditer || (fs.ft == Percent && !fs.numReplacements));

    return result;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void getDimensionsFromInput(const std::vector<Currency> &inputs, int *m, int *n)
{
    // temporary double storage to check if finite
    double dbm = 0.0;
    double dbn = 0.0;
    size_t size = inputs.size();
    if (size == 0)
    {
        *m = *n = -1;
        return;
    }
    else if (size == 1)
    {
        if (size == 1)
        {
            const Currency &input1 = inputs[0];
            if (input1.IsInteger())
            {
                dbm = dbn = input1.Scalar();
            }
            else if (input1.IsMatrix())
            {
                const hwMatrix *mtx = input1.Matrix();
                int matrixSize = mtx->Size();
                if (!matrixSize)
                {
                    throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIMS);
                }
                else if (!mtx->IsRealData())
                {
                    throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);
                }
                else if (matrixSize == 2)
                {
                    if (!IsInteger(realval(mtx, 0)).IsOk())
                        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);
                    
                    if (!IsInteger(realval(mtx, 1)).IsOk())
                        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);
                    
                    dbm = realval(mtx, 0);
                    dbn = realval(mtx, 1);
                }
                else
                {
                    throw OML_Error(OML_ERR_UNSUPPORTDIM);
                }
            }
            else
            {
                throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);
            }
        }
    }
    else if (size == 2)
    {
        if (!inputs[0].IsInteger())
            throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);

        if (!inputs[1].IsInteger())
            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);

        dbm = (int) inputs[0].Scalar();
        dbn = (int) inputs[1].Scalar();
    }
    else
    {
        throw OML_Error(OML_ERR_UNSUPPORTDIM, 1);
    }

    if (!(::isfinite(dbm) && ::isfinite(dbn)))
        throw OML_Error(HW_ERROR_DIMFINITE);

    *m = (int) dbm;
    *n = (int) dbn;

    if (*m < 0)
        throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);

    if (*n < 0)
        throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
static bool isequal(double exp, double obs, double tol)
{
    bool use_abs = iszero(exp) || tol > 0;
    double t = std::abs(use_abs ? tol : tol*exp);
    return std::abs(obs - exp) <= t; 
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
static bool isequal(double exp, double obs, const Currency& tol)
{
    if (tol.IsScalar())
        return isequal(exp, obs, tol.Scalar());
    else if (tol.IsComplex())
        return isequal(exp, obs, tol.Complex().Mag());
    else if (tol.IsMatrix() || tol.IsString())
        return false;
    else
        throw OML_Error(HW_ERROR_INVTOL);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
static bool isequal(const hwComplex& exp, const hwComplex& obs, double tol)
{
    bool use_abs = iszero(exp) || tol > 0;
    double t = use_abs ?  std::abs(tol) : (tol*exp).Mag();
    return (obs - exp).Mag() <= t; 
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
static bool isequal(const hwComplex& exp, const hwComplex& obs, const Currency& tol)
{
    if (tol.IsScalar())
        return isequal(exp, obs, tol.Scalar());
    else if (tol.IsComplex())
        return isequal(exp, obs, tol.Complex().Mag());
    else if (tol.IsMatrix() || tol.IsString())
        return false;
    else
        throw OML_Error(HW_ERROR_INVTOL);
}
//------------------------------------------------------------------------------
// if expected is 0 or tol is positive, it is used as an absolute tolerance.
//------------------------------------------------------------------------------
bool isequal(const Currency &observed, const Currency &expected, const Currency& tol)
{
    if (expected.IsScalar())
    {
        if (observed.IsScalar())
            return isequal(expected.Scalar(), observed.Scalar(), tol);
        else if (observed.IsComplex())
            return isequal(hwComplex(expected.Scalar(), 0.0), observed.Complex(), tol);
        return false;
    }
    else if (expected.IsComplex())
    {
        if (observed.IsScalar())
            return isequal(expected.Complex(), hwComplex(observed.Scalar(), 0.0), tol);
        else if (observed.IsComplex())
            return isequal(expected.Complex(), observed.Complex(), tol);
        return false;
    }
    else if (expected.IsMatrix() || expected.IsString())
    {
        if ((observed.IsMatrix() || observed.IsString()) && observed.GetMask() == expected.GetMask())
        {
            const hwMatrix* e = expected.Matrix();
            const hwMatrix* o = observed.Matrix();

            if (sameSize(e,o))
            {
                bool tol_singular = tol.IsScalar() || tol.IsComplex();
                if (!tol_singular && !tol.IsMatrix())
                    throw OML_Error(HW_ERROR_INVTOL);

                if (tol_singular || sameSize(e, tol.Matrix()))
                {
                    if (e->IsReal())
                    {
                        if (o->IsReal())
                        {
                            for (int i = 0; i < e->Size(); ++i)
                            {
                                if (tol_singular)
                                {
                                    if (!isequal((*e)(i), (*o)(i), tol))
                                        return false;
                                }
                                else
                                {
                                    double t = tol.Matrix()->IsReal() ? (*tol.Matrix())(i) : tol.Matrix()->z(i).Mag();
                                    if (!isequal((*e)(i), (*o)(i), t))
                                        return false;
                                }
                            }
                        }
                        else
                        {
                            for (int i = 0; i < e->Size(); ++i)
                            {
                                if (tol_singular)
                                {
                                    if (!isequal(hwComplex((*e)(i),0.0), o->z(i), tol))
                                        return false;
                                }
                                else
                                {
                                    double t = tol.Matrix()->IsReal() ? (*tol.Matrix())(i) : tol.Matrix()->z(i).Mag();
                                    if (!isequal(hwComplex((*e)(i),0.0), o->z(i), t))
                                        return false;
                                }
                            }
                        }
                    }
                    else // e is complex
                    {
                        if (o->IsReal())
                        {
                            for (int i = 0; i < e->Size(); ++i)
                            {
                                if (tol_singular)
                                {
                                    if (!isequal(e->z(i), hwComplex((*o)(i),0.0), tol))
                                        return false;
                                }
                                else
                                {
                                    double t = tol.Matrix()->IsReal() ? (*tol.Matrix())(i) : tol.Matrix()->z(i).Mag();
                                    if (!isequal(e->z(i), hwComplex((*o)(i),0.0), t))
                                        return false;
                                }
                            }
                        }
                        else
                        {
                            for (int i = 0; i < e->Size(); ++i)
                            {
                                if (tol_singular)
                                {
                                    if (!isequal(e->z(i), o->z(i), tol))
                                        return false;
                                }
                                else
                                {
                                    double t = tol.Matrix()->IsReal() ? (*tol.Matrix())(i) : tol.Matrix()->z(i).Mag();
                                    if (!isequal(e->z(i), o->z(i), t))
                                        return false;
                                }
                            }
                        }
                    }
                }
                else
                    return false;
                return true;
            }
        }
        return false;
    }
    else if (expected.IsNDMatrix() || observed.IsNDMatrix())
    {
        return oml_MatrixNUtil10(observed, expected, tol, isequal);
    }
    else if (expected.IsCellArray())
    {
        return observed.IsCellArray() && cellArraysEqual(observed.CellArray(), expected.CellArray(), &tol);
    }
    else if (expected.IsStruct())
    {
        return observed.IsStruct() && structsEqual(observed.Struct(), expected.Struct(), &tol);
    }
    else if (expected.IsFunctionHandle())
    {
        return observed.IsFunctionHandle() && funcsEqual(observed.FunctionHandle(), expected.FunctionHandle());
    }
    else
        throw OML_Error(HW_MATH_MSG_INVALIDINPUT);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isequal(const Currency &c1, const Currency &c2)
{
    if (c1.IsScalar())
    {
        return c2.IsScalar() && c1.Scalar() == c2.Scalar();
    }
    else if (c1.IsComplex())
    {
        return c2.IsComplex() && c1.Complex() == c2.Complex();
    }
    else if (c1.IsMatrix())
    {
		if (c1.Matrix())
		{
			return c2.IsMatrix() && *c1.Matrix() == *c2.Matrix();
		}
		else if (c2.IsMatrix())
		{
			const hwMatrix* m2 = c2.Matrix();
			return ((m2->M() == 0) && (m2->N() == 0));
		}
		else
		{
			return false;
		}
    }
    else if (c1.IsNDMatrix() || c2.IsNDMatrix())
    {
        return oml_MatrixNUtil9(c1, c2, isequal);
    }
    else if (c1.IsString())
    {
        return c2.IsString() && *c1.Matrix() == *c2.Matrix();
    }
    else if (c1.IsCellArray())
    {
        return c2.IsCellArray() && cellArraysEqual(c1.CellArray(), c2.CellArray());
    }
    else if (c1.IsStruct())
    {
        return c2.IsStruct() && structsEqual(c1.Struct(), c2.Struct());
    }
	else if (c1.IsObject())
    {
		if (c2.IsObject())
		{
			if (c1.GetClassname() == c2.GetClassname())
				return structsEqual(c1.Struct(), c2.Struct());
		}

		return false;
    }
    else if (c1.IsFunctionHandle())
    {
        return c2.IsFunctionHandle() && funcsEqual(c1.FunctionHandle(), c2.FunctionHandle());
    }
    else
	{
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING);
	}
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency getAtIndex(const Currency &cur, int i, int j)
{
    if (i < 0 || j < 0)
        throw OML_Error(HW_ERROR_INDEXPOSINT);

    if (cur.IsScalar() || cur.IsComplex())
    {
        if (!(i == 0 && j == 0))
            throw OML_Error(HW_MATH_MSG_INVALIDINDEX);
        return cur;
    }
    else if (cur.IsMatrix() || cur.IsString())
    {
        const hwMatrix *mtx = cur.Matrix();
        if (i >= mtx->M() || j >= mtx->N())
            throw OML_Error(HW_MATH_MSG_INVALIDINDEX);

        Currency rv;
        if (cur.IsMatrix())
        {
            if (mtx->IsReal())
                rv = (*mtx)(i, j);
            else
                rv = mtx->z(i, j);
        }
        else
        {
            rv = EvaluatorInterface::allocateMatrix(1, 1, (*mtx)(i, j));
            rv.SetMask(Currency::MASK_STRING);
        }

        return rv;
    }
    else if (cur.IsCellArray())
    {
        HML_CELLARRAY *cell = cur.CellArray();

        if (i >= cell->M() || j >= cell->N())
            throw OML_Error(HW_MATH_MSG_INVALIDINDEX);
        HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(1, 1);
        (*newcell)(0) = (*cell)(i, j);
        return newcell;
    }
    else if (cur.IsStruct())
    {
        StructData *sd = cur.Struct();
        if (i >= sd->M() || j >= sd->N())
            throw OML_Error(HW_MATH_MSG_INVALIDINDEX);
        return sd->GetElement(i + 1, j + 1);
    }
    else
        throw OML_Error(HW_ERROR_NOTINDTYPEOBJ);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency getAtIndex(const Currency &cur, int index)
{
    if (index < 0)
        throw OML_Error(HW_ERROR_INDEXPOSINT);

    if (cur.IsScalar() || cur.IsComplex())
    {
        if (index != 0)
            throw OML_Error(HW_MATH_MSG_INVALIDINDEX);
        return cur;
    }
    else if (cur.IsMatrix() || cur.IsString())
    {
        const hwMatrix *mtx = cur.Matrix();
        if (index >= mtx->Size())
            throw OML_Error(HW_MATH_MSG_INVALIDINDEX);

        Currency rv;
        if (cur.IsMatrix())
        {
            if (mtx->IsReal())
                rv = (*mtx)(index);
            else
                rv = mtx->z(index);
        }
        else
        {
            rv = EvaluatorInterface::allocateMatrix(1, 1, (*mtx)(index));
            rv.SetMask(Currency::MASK_STRING);
        }

        return rv;
    }
    else if (cur.IsCellArray())
    {
        HML_CELLARRAY *cell = cur.CellArray();
        if (index >= cell->Size())
            throw OML_Error(HW_MATH_MSG_INVALIDINDEX);
        HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(1, 1);
        (*newcell)(0) = (*cell)(index);
        return newcell;
    }
    else if (cur.IsStruct())
    {
        StructData *sd = cur.Struct();
        if (index >= sd->M() * sd->N())
            throw OML_Error(HW_MATH_MSG_INVALIDINDEX);
        return sd->GetElement(index + 1, -1);
    }
    else
        throw OML_Error(HW_ERROR_NOTINDTYPEOBJ);
}
//------------------------------------------------------------------------------
// warning: reads top row of strings -- not full matrices or StringVal()
//------------------------------------------------------------------------------
std::string concatRowToString(EvaluatorInterface& eval, int row, const HML_CELLARRAY *cell)
{
    std::stringstream rv;
    for (int j = 0; j < cell->N(); j++)
    {
        Currency temp = (*cell)(row, j);
        if (!temp.IsString())
            throw OML_Error(HW_ERROR_CELLELEMSTR);
        rv << readString(temp);
    }
    return rv.str();
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency unnest(Currency nested, const std::string &errmsg)
{
    while (nested.IsCellArray())
    {
        HML_CELLARRAY *cell = nested.CellArray();
        if (cell->IsEmpty())
            throw OML_Error(errmsg);
        nested = (*cell)(0);
    }
    return nested;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void makeInt(hwComplex &c)
{
    c.Real() = round(c.Real());
    c.Imag() = round(c.Imag());
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency transpose(EvaluatorInterface& eval, const Currency &cur)
{
    if (cur.IsMatrix() || cur.IsString())
    {
        Currency ret = _transpose(eval, cur.Matrix());
        ret.SetMask(cur.GetMask());
        return ret;
    }
    else if(cur.IsCellArray())
    {
        return _transpose(eval, cur.CellArray());
    }
    else if(cur.IsStruct())
    {
        StructData *s = cur.Struct();
        StructData *result = new StructData();
        result->DimensionNew(s->N(), s->M());
        Currency ret(result);

        std::map<std::string, int> fields = s->GetFieldNames();
        std::map<std::string, int>::const_iterator iter;
        for (iter = fields.cbegin(); iter != fields.cend(); iter++)
        {
            for (int i = 0; i < s->M(); i++)
            {
                for (int j = 0; j < s->N(); j++)
                {
                    result->SetValue(j, i, iter->first, s->GetValue(i, j, iter->first));
                }
            }
        }
        return ret;
    }
    return cur;
}

// removes zeros padding the outside of the matrix
hwMatrix* removePadding(const hwMatrix *mtx)
{
    // these stand for the number of dimensions to remove

    int m = mtx->M(), n = mtx->N();
    int top = 0, bottom = m, left = 0, right = n;

    void (*countMethod)(const hwMatrix*, int*, int, int, int, bool) = mtx->IsReal() ? &_countZeros : &_countZerosComplex;

    // determine how many rows/cols from each side to remove
    (*countMethod)(mtx, &top, m, n, 1, true);
    if (top == m)
        return EvaluatorInterface::allocateMatrix();
    (*countMethod)(mtx, &bottom, m, n, -1, true);
    (*countMethod)(mtx, &left, n, m, 1, false);
    (*countMethod)(mtx, &right, n, m, -1, false);

    return getInnerMatrix(mtx, top, bottom, left, right);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwMatrix *getInnerMatrix(const hwMatrix *mtx, int top, int bottom, int left, int right)
{
    if (top > bottom || left > right)
        return EvaluatorInterface::allocateMatrix (max(bottom - top, 0), max(right - left, 0), hwMatrix::REAL);

    hwMatrix *ret;
    if (mtx->IsReal())
        ret = EvaluatorInterface::allocateMatrix(bottom - top, right - left, 0.0);
    else
        ret = EvaluatorInterface::allocateMatrix(bottom - top, right - left, hwComplex());

    // correct invalid dimensions -- values outside of the matrix will be padded with zeros
    if (top < 0)
        top = 0;
    if (left < 0)
        left = 0;
    if (bottom > mtx->M())
        bottom = mtx->M();
    if (right > mtx->N())
        right = mtx->N();

    if (mtx->IsReal())
    {
        for (int i = top; i < bottom; i++)
            for (int j = left; j < right; j++)
                (*ret)(i - top, j - left) = (*mtx)(i, j);
    }
    else
    {
        for (int i = top; i < bottom; i++)
            for (int j = left; j < right; j++)
                ret->z(i - top, j - left) = mtx->z(i, j);
    }

    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::deque<int> getprimes(int max)
{
    bool *isprime;
    std::deque<int> p;

    if (max > 0)
        isprime = new bool[max - 1]; // offset to range from 2 to max inclusively

    for (int i = 0; i < max - 1; ++i)
        isprime[i] = true;

    int sqrtmax = (int) sqrt((double) max);
    for (int i = 2; i <= sqrtmax; i++)
    {
        if (isprime[i - 2])
        {
            for (int j = i * i; j <= max; j += i)
            {
                isprime[j - 2] = false;
            }
        }
    }

    for (int i = 0; i < max - 1; i++)
    {
        if (isprime[i])
            p.push_back(i + 2);
    }

    if (max > 0)
        delete [] isprime;
    return p;
}
//------------------------------------------------------------------------------
// This algorithm takes advantage of the fact that all primes except 2 and 3 are of the form 6k+1 or 6k-1.
// So mod 2 eliminates half of numbers, leaving only odd numbers to check.
// mod 3 eliminates 1/3 of all odd numbers-- since every third odd is a multiple of 3. 
// Therefore, the increment switches between 2 and 4, to skip all odd numbers that are multiples of 3. 
// The result is therefore a bit more effecient than checking all odd numbers. 
//------------------------------------------------------------------------------
bool isprime(int p)
{
	if(p==1 || p==0) // 1 and 0 are not prime
		return false;

	if(p==2 || p==3) // checks 2 and 3, which are prime
		return true;

	if(p%2 == 0 || p%3 == 0) // eliminates multiples of 2 and 3
		return false;

	int odd = 5;
	int inc = 2;

	int max = (int)sqrt((double)p);

	while(odd <= max) // this loops checks all odd numbers except for multiples of 3, already checked above
	{
		if(p%odd == 0) //if mod 0, not prime
			return false;

		odd += inc;
		inc = 6 - inc;
	}
	return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool fileExists(const std::string &file_name)
{
#ifdef OS_WIN
    DWORD dwAttrib = GetFileAttributes((LPCSTR) file_name.c_str());
    return (dwAttrib != INVALID_FILE_ATTRIBUTES && !(dwAttrib & FILE_ATTRIBUTE_DIRECTORY));
#else
    struct stat st;
    return (!stat(file_name.c_str(), &st));
#endif
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isEscaped(const std::string &str, size_t index)
{
    int num_back_slahes = 0;
    while (index && str[--index] == '\\')
        ++num_back_slahes;
    return num_back_slahes % 2;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency findUnsorted(EvaluatorInterface& eval, std::deque<hwMatrix> vals, const hwMatrix *searchin, bool forward)
{
    int numrows = (int) vals.size();
    hwMatrix *row = EvaluatorInterface::allocateMatrix();
    hwMatrix *inpIndices = EvaluatorInterface::allocateMatrix(numrows, 1, hwMatrix::REAL);
    Currency ret(inpIndices), rowcur(row);
    if (forward)
    {
        for (int i = 0; i < numrows; i++)
        {
            hwMatrix *tofind = &vals[i];
            for (int j = 0; j < searchin->M(); j++)
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, searchin->ReadRow(j, *row));
                if (*row == *tofind)
                {
                    (*inpIndices)(i) = j + 1;
                    break;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < numrows; i++)
        {
            hwMatrix *tofind = &vals[i];
            for (int j = searchin->M() - 1; j >= 0; j--)
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, searchin->ReadRow(j, *row));
                if (*row == *tofind)
                {
                    (*inpIndices)(i) = j + 1;
                    break;
                }
            }
        }
    }
    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency findUnsorted(std::deque<double> vals, const hwMatrix *searchin, bool forward)
{
    int valsize = (int) vals.size();
    hwMatrix *inpIndices;
    if (searchin->M() == 1)
        inpIndices = EvaluatorInterface::allocateMatrix(1, valsize, hwMatrix::REAL);
    else
        inpIndices = EvaluatorInterface::allocateMatrix(valsize, 1, hwMatrix::REAL);

    Currency ret(inpIndices);

    if (forward)
    {
        for (int i = 0; i < valsize; i++)
        {
            double tofind = vals[i];
            for (int j = 0; j < searchin->Size(); j++)
            {
                double temp = (*searchin)(j);
                if (temp == tofind)
                {
                    (*inpIndices)(i) = j + 1;
                    break;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < valsize; i++)
        {
            double tofind = vals[i];
            for (int j = searchin->Size() - 1; j >= 0; j--)
            {
                double temp = (*searchin)(j);
                if (temp == tofind)
                {
                    (*inpIndices)(i) = j + 1;
                    break;
                }
            }
        }
    }
    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency findUnsorted(std::deque<hwComplex> vals, const hwMatrix *searchin, bool forward)
{
    int valsize = (int) vals.size();
    hwMatrix *inpIndices;
    if (searchin->M() == 1)
        inpIndices = EvaluatorInterface::allocateMatrix(1, valsize, hwMatrix::REAL);
    else
        inpIndices = EvaluatorInterface::allocateMatrix(valsize, 1, hwMatrix::REAL);

    Currency ret(inpIndices);

    if (forward)
    {
        for (int i = 0; i < valsize; i++)
        {
            hwComplex tofind = vals[i];
            for (int j = 0; j < searchin->Size(); j++)
            {
                hwComplex temp = searchin->z(j);
                if (temp == tofind)
                {
                    (*inpIndices)(i) = j + 1;
                    break;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < valsize; i++)
        {
            hwComplex tofind = vals[i];
            for (int j = searchin->Size() - 1; j >= 0; j--)
            {
                hwComplex temp = searchin->z(j);
                if (temp == tofind)
                {
                    (*inpIndices)(i) = j + 1;
                    break;
                }
            }
        }
    }
    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency findUnsorted(EvaluatorInterface& eval, std::deque<std::string> vals, const HML_CELLARRAY *searchin, bool forward)
{
    hwMatrix *inpIndices;
    int valsize = (int)vals.size();
    if (searchin->M() == 1)
        inpIndices = EvaluatorInterface::allocateMatrix(1, valsize, hwMatrix::REAL);
    else
        inpIndices = EvaluatorInterface::allocateMatrix(valsize, 1, hwMatrix::REAL);

    Currency ret(inpIndices);

    if (forward)
    {
        for (int i = 0; i < valsize; i++)
        {
            std::string tofind = vals[i];
            for (int j = 0; j < searchin->Size(); j++)
            {
                std::string temp = readString((*searchin)(j));
                if (temp == tofind)
                {
                    (*inpIndices)(i) = j + 1;
                    break;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < valsize; i++)
        {
            std::string tofind = vals[i];
            for (int j = searchin->Size() - 1; j >= 0; j--)
            {
                std::string temp = readString((*searchin)(j));
                if (temp == tofind)
                {
                    (*inpIndices)(i) = j + 1;
                    break;
                }
            }
        }
    }
    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency findSorted(EvaluatorInterface& eval, std::deque<hwMatrix> vals, const hwMatrix *searchin)
{
    hwMatrix *outIndices = EvaluatorInterface::allocateMatrix(searchin->M(), 1, hwMatrix::REAL);
    hwMatrix row;
    Currency ret(outIndices);

    for (int i = 0; i < searchin->M(); i++)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, searchin->ReadRow(i, row));
        (*outIndices)(i) = indexOf(vals, row, false, true) + 1;
    }
    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency findSorted(std::deque<double> vals, const hwMatrix *searchin)
{
    hwMatrix *outIndices;
    if (searchin->M() == 1)
        outIndices = EvaluatorInterface::allocateMatrix(1, searchin->Size(), hwMatrix::REAL);
    else
        outIndices = EvaluatorInterface::allocateMatrix(searchin->Size(), 1, hwMatrix::REAL);

    Currency ret(outIndices);

    for (int i = 0; i < searchin->Size(); i++)
        (*outIndices)(i) = indexOf(vals, (*searchin)(i), false, true) + 1;
    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency findSorted(std::deque<hwComplex> vals, const hwMatrix *searchin)
{
    hwMatrix *outIndices;
    if (searchin->M() == 1)
        outIndices = EvaluatorInterface::allocateMatrix(1, searchin->Size(), hwMatrix::REAL);
    else
        outIndices = EvaluatorInterface::allocateMatrix(searchin->Size(), 1, hwMatrix::REAL);

    Currency ret(outIndices);

    for (int i = 0; i < searchin->Size(); i++)
        (*outIndices)(i) = indexOf(vals, searchin->z(i), false, true) + 1;
    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency findSorted(EvaluatorInterface& eval, std::deque<std::string> vals, const HML_CELLARRAY *searchin)
{
    hwMatrix *outIndices;
    if (searchin->M() == 1)
        outIndices = EvaluatorInterface::allocateMatrix(1, searchin->Size(), hwMatrix::REAL);
    else
        outIndices = EvaluatorInterface::allocateMatrix(searchin->Size(), 1, hwMatrix::REAL);

    Currency ret(outIndices);

    for (int i = 0; i < searchin->Size(); i++)
        (*outIndices)(i) = indexOf(vals, readString((*searchin)(i)), false, true) + 1;
    return ret;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
int posIntFromDouble(double d)
{
    return isposint(d) ? (int) d : 0;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void writeCol(EvaluatorInterface& eval, hwMatrix* mtx, hwMatrix* col, int index)
{
    BuiltInFuncsUtils::CheckMathStatus(eval, mtx->WriteColumn(index, *col));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void writeRow(EvaluatorInterface& eval, hwMatrix* mtx, hwMatrix* row, int index)
{
    BuiltInFuncsUtils::CheckMathStatus(eval, mtx->WriteRow(index, *row));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void writeCol(EvaluatorInterface& eval, HML_CELLARRAY* cell, HML_CELLARRAY* col, int index)
{
    for (int i = 0; i < cell->M(); ++i)
        (*cell)(i, index) = (*col)(i);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void writeRow(EvaluatorInterface& eval, HML_CELLARRAY* cell, HML_CELLARRAY* row, int index)
{
    for (int i = 0; i < cell->N(); ++i)
        (*cell)(index, i) = (*row)(i);
}

// matrix operations
//------------------------------------------------------------------------------
// checks if either input matrix is complex
// if one is, converts them both to complex and returns true
//------------------------------------------------------------------------------
bool checkMakeComplex(EvaluatorInterface& eval, hwMatrix *mtx1, hwMatrix *mtx2)
{
    if (!mtx1->IsReal())
    {
        if (mtx2->IsReal())
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx2->MakeComplex());
    }
    else if (!mtx2->IsReal())
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx1->MakeComplex());
    }
    else
        return false;
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool checkMakeComplex(EvaluatorInterface& eval, hwMatrix *mtx1, hwMatrix *mtx2, hwMatrix *mtx3)
{
    if (!mtx1->IsReal())
    {
        if (mtx2->IsReal())
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx2->MakeComplex());
        if (mtx3->IsReal())
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx3->MakeComplex());
    }
    else if (!mtx2->IsReal())
    {
        if (mtx1->IsReal())
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx1->MakeComplex());
        if (mtx3->IsReal())
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx3->MakeComplex());
    }
    else if (!mtx3->IsReal())
    {
        if (mtx1->IsReal())
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx1->MakeComplex());
        if (mtx2->IsReal())
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx2->MakeComplex());
    }
    else
        return false;
    return true;
}

// data checkers
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool checkisfinite(double d)
{
    double inf = std::numeric_limits<double>::infinity();
    if (IsNaN_T(d) || d == inf || d == -inf)
        return false;
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isfinite(hwComplex c)
{
    if (!(checkisfinite(c.Real()) && checkisfinite(c.Imag())))
        return false;
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isfinite(hwMatrix *m)
{
    if (m->IsReal())
    {
        for (int i = 0; i < m->Size(); i++)
        {
            if (!checkisfinite((*m)(i)))
                return false;
        }
    }
    else
    {
        for (int i = 0; i < m->Size(); i++)
        {
            if (!isfinite(m->z(i)))
                return false;
        }
    }
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isposint(double d)
{ 
    return isint(d) && d > 0.0;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool islonglong(double d)
{ 
    return abs(((long long) d) - d) < HW_TOLERANCE; 
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isint(double d)
{ 
    return abs(((int) d) - d) < HW_TOLERANCE;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isint(const hwMatrix *mtx)
{
    if (mtx->IsReal())
    {
        for (int i = 0; i < mtx->Size(); i++)
        {
            if (!isint((*mtx)(i)))
                return false;
        }
        return true;
    }
    else
    {
        for (int i = 0; i < mtx->Size(); i++)
        {
            if (!isint(mtx->z(i)))
                return false;
        }
        return true;
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isinfinity(double d)
{
    return std::numeric_limits<double>::has_infinity && (abs(d) == std::numeric_limits<double>::infinity());
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool isfinite(double d)
{
    return !(IsNaN_T(d) || isinfinity(d));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double realval(const hwMatrix* mtx, int i, int j)
{
    return mtx->IsReal() ? (*mtx)(i, j) : mtx->z(i, j).Real();
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double realval(const hwMatrix* mtx, int index)
{
    return mtx->IsReal() ? (*mtx)(index) : mtx->z(index).Real();
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double realvalorscalar(const hwMatrix* mtx, int i, int j)
{
    if (mtx->Size() == 1)
        return mtx->IsReal() ? (*mtx)(0) : mtx->z(0).Real();
    return realval(mtx, i, j);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double realvalorscalar(const hwMatrix* mtx, int index)
{
    if (mtx->Size() == 1)
        return mtx->IsReal() ? (*mtx)(0) : mtx->z(0).Real();
    return realval(mtx, index);
}
//------------------------------------------------------------------------------
// Checks if given string is a directory, strips trailing slash and adds to path
//------------------------------------------------------------------------------
void checkAddPath(EvaluatorInterface &eval, std::string &str, bool appendToEnd)
{
    if (checkIsDirectory(eval, str, false))
    {
        BuiltInFuncsUtils::StripTrailingSlash(str);
        eval.AddPath(str, appendToEnd);
    }
}
//------------------------------------------------------------------------------
// Returns path as string 
//------------------------------------------------------------------------------
std::string getPathString(EvaluatorInterface &eval, char pathSep)
{
    std::string out;
    const std::vector<std::string> &paths = eval.GetPaths();
    std::vector<std::string>::const_iterator iter;
    for (iter = paths.begin(); iter != paths.end(); iter++)
    {
        if (iter != paths.begin())
            out += pathSep;
        out += *iter;
    }
    return out;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::vector<std::string> separatePathNames(EvaluatorInterface& eval, Currency cur)
{
    std::vector<std::string> pathlist;
    if (cur.IsString() || cur.IsMatrix())
        cur = readRow(eval, cur);

    std::string inputstr = orderedStringVal(toCurrencyStr(eval, cur, false, false));

    while (inputstr.length())
    {
        size_t sepIndex = inputstr.find(pathsep);
        pathlist.push_back(inputstr.substr(0, sepIndex));
        if (sepIndex != std::string::npos)
            sepIndex++;
        inputstr.erase(0, sepIndex);
    }
    return pathlist;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void checkFileIndex(EvaluatorInterface &eval, int i, bool checkStdStreams)
{
    if (i < (checkStdStreams ? 0 : FIRST_USER_FILE) || i >= eval.GetNumFiles()
        || nullptr == eval.GetFile(i))
    {
        std::stringstream str;
        str << "Error: invalid file stream id: " << i;
        throw OML_Error(str.str().c_str());
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void safeResizeMatrix(EvaluatorInterface& eval, int m, int n, hwMatrix &A, bool initZero)
{
    hwMathStatus stat = hwMathStatus();
    int oldm = A.M(), oldn = A.N();
    if (m > oldm)
    {
        if (n > oldn)
            stat = A.Resize(m, n, initZero);
        else
            stat = A.Resize(m, oldn, initZero);
    }
    else if (n > oldn)
        stat = A.Resize(oldm, n, initZero);

    if (stat.IsOk())
    {
        if (m < oldm)
            stat = A.DeleteRows(m, oldm - m);
        if (n < oldn)
            stat = A.DeleteColumns(n, oldn - n);
    }
    BuiltInFuncsUtils::CheckMathStatus(eval, stat);

    if (A.Size() != m*n)
        throw OML_Error(HW_ERROR_OUTMEM);
}
//------------------------------------------------------------------------------
// helper for oml_factor
// assumes primes and multiplics have same size and indices for each value are lined up
//------------------------------------------------------------------------------
void addValMultiplicity(long long value, bool useTwoVecs, std::vector<long long> &primes, std::vector<int> &multiplics)
{
    if (useTwoVecs)
    {
        for (int i = 0; i < primes.size(); i++)
        {
            if (value == primes[i])
            {
                multiplics[i]++;
                return;
            }
        }
        primes.push_back(value);
        multiplics.push_back(1);
    }
    else
        primes.push_back(value);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
hwMatrix* checkMatrixFinite(hwMatrix *m)
{
    hwMatrix *returnVal = EvaluatorInterface::allocateMatrix(m->M(), m->N(), hwMatrix::REAL);
    if (m->IsReal())
    {
        for (int i = 0; i < m->Size(); i++)
        {
            if (checkisfinite((*m)(i)))
                (*returnVal)(i) = 1.0;
            else
                (*returnVal)(i) = 0.0;
        }
    }
    else
    {
        for (int i = 0; i < m->Size(); i++)
        {
            if (isfinite(m->z(i)))
                (*returnVal)(i) = 1.0;
            else
                (*returnVal)(i) = 0.0;
        }
    }
    return returnVal;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool cellArraysEqual(const HML_CELLARRAY *c1, const HML_CELLARRAY *c2, const Currency* tol)
{
    if (!sameSize(c1, c2))
        return false;

    if (tol)
    {
        for (int i = 0; i < c1->Size(); i++)
        {
            if (!isequal((*c1)(i), (*c2)(i), *tol))
                return false;
        }
    }
    else
    {
        for (int i = 0; i < c1->Size(); i++)
        {
            if (!isequal((*c1)(i), (*c2)(i)))
                return false;
        }
    }

    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool structsEqual(const StructData* s1, const StructData* s2, const Currency* tol)
{
    if (!sameSize(s1, s2))
        return false;

    const std::map<std::string, int> fnames1 = s1->GetFieldNames();
    const std::map<std::string, int> fnames2 = s2->GetFieldNames();

    if (fnames1.size() != fnames2.size())
        return false;

    std::map<std::string, int>::const_iterator iter;
    
    for (iter = fnames1.begin(); iter != fnames1.end(); ++iter)
    {
        const std::string key = iter->first;        
        if (!fnames2.count(key))
            return false;

        for (int i = 0; i < s1->M(); ++i)
        {
            for (int j = 0; j < s1->N(); ++j)
            {
                if (tol)
                {
                    if (!isequal(s1->GetValue(i, j, key), s2->GetValue(i, j, key), *tol))
                        return false;
                }
                else
                {
                    if (!isequal(s1->GetValue(i, j, key), s2->GetValue(i, j, key)))
                        return false;
                }
            }
        }
    }
    return true;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool funcsEqual(FunctionInfo* f1, FunctionInfo* f2)
{
    if (f1->IsBuiltIn() != f2->IsBuiltIn())
        return false;

    if (f1->IsBuiltIn() || f2->IsBuiltIn())
        return f1->FunctionName() == f2->FunctionName();

    if (f1->IsAnonymous() != f2->IsAnonymous())
        return false;

    return f1->Statements() == f2->Statements();
}
//------------------------------------------------------------------------------
// helper for cell2mat
//------------------------------------------------------------------------------
void increasetemp(int *list, int &temp, int rowindex, int max)
{
    while (temp < max && list[temp] > rowindex)
        temp++;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::string makeString(EvaluatorInterface& eval, double d, int precision)
{		
	if (IsNaN_T(d))
		return "NaN";
	else if (IsInf_T(d))
		return "Inf";
	else if (IsNegInf_T(d))
		return "-Inf";

    bool usell = islonglong(d);
    if (precision > -1)
    {
        int numdigits = 0;
        long long temp = (long long) d;

        while (temp > 1)
        {
            numdigits++;
            temp /= 10;
        }

        std::string tmplt("%");
        if (!(usell && precision >= numdigits))
            tmplt += '.';
        tmplt += std::to_string((long long) precision) + 'g';
        return makeString(eval, d, tmplt);
    }

    if (usell)
        return std::to_string((long long) d);
    return std::to_string((long double) d);

}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::string makeString(EvaluatorInterface& eval, double d, const std::string &format)
{
	if (IsNaN_T(d))
		return "NaN";
	else if (IsInf_T(d))
		return "Inf";
	else if (IsNegInf_T(d))
		return "-Inf";

    std::vector<Currency> tosub;
    tosub.push_back(Currency(d));
    std::string result = sprintf(eval, format, tosub.cbegin(), tosub.cend());

    // trim left
    size_t start = result.find_first_not_of(' ');
    if (start == std::string::npos)
        return std::string();

    result.erase(result.begin(), result.begin() + start);
    return result;
}
//------------------------------------------------------------------------------
// Returns string for a complex number
//------------------------------------------------------------------------------
template <typename T>
std::string makeString(EvaluatorInterface& eval, const hwComplex &c, T precForm)
{
    std::string out (makeString(eval, c.Real(), precForm));

    double imag = c.Imag();
    if (IsInf_T(imag) || IsNaN_T(imag) || imag >= 0) 
        out += '+';   // Add '+' when imag is Inf/Nan/non-negative

    out += makeString(eval, imag, precForm) + 'i';
    return out;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void stringVecFromCurrencyRows(EvaluatorInterface& eval, const Currency &input1, const Currency &input2, std::vector<Currency> &searchfor,
    std::vector<Currency> &searchin, int *m, int *n)
{
    if (input1.IsCellArray() || input2.IsCellArray())
        throw OML_Error(HW_ERROR_NOTUSECELLSEARCHFORROW);

    if (!input1.IsString())
        OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

    if (!input2.IsString())
        OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);

    const hwMatrix *mtx1 = input1.Matrix();
    const hwMatrix *mtx2 = input2.Matrix();

    if (mtx1->N() != mtx2->N())
        throw OML_Error(HW_ERROR_STRSAMENUMOFCOLWCOMPROW);

    tryPushBackStringByRows(eval, mtx1, searchfor);
    tryPushBackStringByRows(eval, mtx2, searchin);

    *n = 1;
    *m = mtx1->M();
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void stringVecFromCurrency(EvaluatorInterface& eval, const Currency &input1, const Currency &input2, std::vector<Currency> &searchfor,
    std::vector<Currency> &searchin, int *m, int *n)
{
    if (input1.IsString())
    {
        const hwMatrix *mtx = input1.Matrix();

        if (input2.IsCellArray())
        {
            *m = mtx->M();
            *n = 1;

            tryPushBackStringByRows(eval, mtx, searchfor);

            HML_CELLARRAY *cell = input2.CellArray();
            for (int i = 0; i < cell->Size(); i++)
            {
                const Currency &elem = (*cell)(i);
                if (elem.IsString())
                {
                    const hwMatrix *mtx2 = elem.Matrix();
                    for (int j = 0; j < mtx2->M(); j++)
                        searchin.push_back(readRow(eval, elem, j));
                }
                else
                    throw OML_Error(HW_ERROR_CELLELEMSTR);
            }
        }
        else if (input2.IsString())
        {
            const hwMatrix *mtx2 = input2.Matrix();
            *m = mtx->M();
            *n = mtx->N();

            tryPushBackString(mtx, searchfor);
            tryPushBackString(mtx2, searchin);
        }
        else
            throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);
    }
    else if (input1.IsCellArray())
    {
        HML_CELLARRAY *cell = input1.CellArray();
        *m = cell->M();
        *n = cell->N();

        tryPushBackString(cell, searchfor);

        if (input2.IsCellArray())
        {
            tryPushBackString(input2.CellArray(), searchin);
        }
        else if (input2.IsString())
        {
            tryPushBackStringByRows(eval, input2.Matrix(), searchin);
        }
        else
            throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);
    }
    else
        throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void tryPushBackString(const HML_CELLARRAY *cell, std::vector<Currency> &topush)
{
    for (int i = 0; i < cell->Size(); i++)
    {
        const Currency &cur = (*cell)(i);
        if (cur.IsString())
            topush.push_back(addStringMask(EvaluatorInterface::allocateMatrix(cur.Matrix())));
        else
            throw OML_Error(HW_ERROR_CELLELEMSTR);
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void tryPushBackStringByRows(EvaluatorInterface& eval, const hwMatrix *mtx, std::vector<Currency> &topush)
{
    for (int i = 0; i < mtx->M(); i++)
        topush.push_back(addStringMask(readRow(eval, mtx, i)));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void tryPushBackString(const hwMatrix *mtx, std::vector<Currency> &topush)
{
    for (int i = 0; i < mtx->Size(); i++)
        topush.push_back((*mtx)(i));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template <typename T>
Currency _transpose(EvaluatorInterface& eval, const T *source)
{
    T *trans = new T;
    Currency outcur(trans);
    BuiltInFuncsUtils::CheckMathStatus(eval, trans->Transpose(*source));
    return outcur;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
Currency getFirstInputOrCell(const std::vector<Currency> &inputs)
{
    size_t nargin = inputs.size();
    for (size_t i = 0 ; i < nargin; i++)
    {
        if (inputs[i].IsCellArray())
        {
            return EvaluatorInterface::allocateCellArray();
        }
    }
    return inputs[0];
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
bool complexLessThanTol(const hwComplex &cplx1, const hwComplex &cplx2, double tol)
{
    double comp = cplx1.Mag() - cplx2.Mag();
    if (comp < -tol)
        return true;
    else if (comp > tol)
        return false;

    comp = cplx1.Arg() - cplx2.Arg();
    return comp < -tol;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
int indexOfTol(const std::deque<hwComplex> &d, const hwComplex &val, double tol)
{
    int size = (int)d.size(), minIndex = 0, maxIndex = size, midIndex = 0;
    while (minIndex <= maxIndex && minIndex < size)
    {
        midIndex = (minIndex + maxIndex) / 2;
        hwComplex curVal = d[midIndex];

        if (curVal.IsEqual(val, tol))
            return midIndex;
        else if (complexLessThanTol(curVal, val, tol))
            minIndex = midIndex + 1;
        else
            maxIndex = midIndex - 1;
    }

    return -1;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::vector<Currency> makeAComplex(EvaluatorInterface& eval, const std::vector<Currency>& inputs)
{
    if (!inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

    if (!inputs[1].IsScalar())
        throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

    double real = inputs[0].Scalar();
    double imag = inputs[1].Scalar();
    std::vector<Currency> results;
    results.push_back(hwComplex(real, imag));
    return results;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double all(EvaluatorInterface& eval, const hwMatrix* mtx)
{
    if (mtx->IsReal())
    {
        for (int i = 0; i < mtx->Size(); ++i)
        {
            if (iszero((*mtx)(i)))
                return 0.0;
        }
        return 1.0;
    }
    else
    {
        for (int i = 0; i < mtx->Size(); ++i)
        {
            if (iszero(mtx->z(i)))
                return 0.0;
        }
        return 1.0;
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
double any(EvaluatorInterface& eval, const hwMatrix* mtx)
{
    if (mtx->IsReal())
    {
        for (int i = 0; i < mtx->Size(); ++i)
        {
            if (!iszero((*mtx)(i)))
                return 1.0;
        }
        return 0.0;
    }
    else
    {
        for (int i = 0; i < mtx->Size(); ++i)
        {
            if (!iszero(mtx->z(i)))
                return 1.0;
        }
        return 0.0;
    }
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
template <bool ASCEND>
hwMatrix* sort(EvaluatorInterface& eval, const hwMatrix* vec, std::pair<int&, hwMatrix*>* index_data)
{
    // find NaN locations and set values to Inf
    std::vector<int> NaNidx;
    hwMatrix* copy = nullptr;

    if (vec->IsReal())
    {
        for (int i = 0; i < vec->Size(); ++i)
        {
            if (IsNaN_T((*vec)(i)))
            {
                NaNidx.push_back(i);
            }
        }
    }
    else
    {
        for (int i = 0; i < vec->Size(); ++i)
        {
            if (IsNaN_T(vec->z(i)))
            {
                NaNidx.push_back(i);
            }
        }
    }

    if (NaNidx.size())
    {
        copy = new hwMatrix(*vec);

        // set the NaN values to Inf when sorting
        if (copy->IsReal())
        {
            for (int i = 0; i < NaNidx.size(); ++i)
                (*copy)(NaNidx[i]) = std::numeric_limits<double>::infinity();
        }
        else
        {
            for (int i = 0; i < NaNidx.size(); ++i)
                copy->z(NaNidx[i]) = std::numeric_limits<double>::infinity();
        }
    }

    // create index vector
    std::vector<int> indexvec(vec->Size());
    for (int i = 0 ; i < vec->Size(); ++i)
        indexvec[i] = i;

    bool row = vec->M() == 1;
    hwMatrix* result = EvaluatorInterface::allocateMatrix(vec->M(), vec->N(), vec->Type());

    if (!copy)
    {
        if (vec->IsReal())
        {
            if (ASCEND)
                std::sort(indexvec.begin(), indexvec.end(), [&indexvec, vec] (int i1, int i2) { return (*vec)(i1) < (*vec)(i2); });
            else
                std::sort(indexvec.begin(), indexvec.end(), [&indexvec, vec] (int i1, int i2) { return (*vec)(i1) > (*vec)(i2); });

            int indexIndex = index_data->first;
            hwMatrix* indices = index_data->second;
            for (int i = 0; i < result->Size(); ++i)
            {
                if (row)
                    (*indices)(indexIndex, i) = indexvec[i] + 1;
                else
                    (*indices)(i, indexIndex) = indexvec[i] + 1;
                (*result)(i) = (*vec)(indexvec[i]);
            }
            ++index_data->first;
            return result;
        }
        else
        {
            if (ASCEND)
                std::sort(indexvec.begin(), indexvec.end(), [&indexvec, vec] (int i1, int i2) { return complexLessThan(vec->z(i1), vec->z(i2)); });
            else
                std::sort(indexvec.begin(), indexvec.end(), [&indexvec, vec] (int i1, int i2) { return complexGreaterThan(vec->z(i1), vec->z(i2)); });

            int indexIndex = index_data->first;
            hwMatrix* indices = index_data->second;
            for (int i = 0; i < result->Size(); ++i)
            {
                if (row)
                    (*indices)(indexIndex, i) = indexvec[i] + 1;
                else
                    (*indices)(i, indexIndex) = indexvec[i] + 1;
                result->z(i) = vec->z(indexvec[i]);
            }
            ++index_data->first;
        }
    }
    else    // use copy to manage NaN values
    {
        if (copy->IsReal())
        {
            if (ASCEND)
                std::sort(indexvec.begin(), indexvec.end(), [&indexvec, copy] (int i1, int i2) { return (*copy)(i1) < (*copy)(i2); });
            else
                std::sort(indexvec.begin(), indexvec.end(), [&indexvec, copy] (int i1, int i2) { return (*copy)(i1) > (*copy)(i2); });
        }
        else
        {
            if (ASCEND)
                std::sort(indexvec.begin(), indexvec.end(), [&indexvec, copy] (int i1, int i2) { return complexLessThan(copy->z(i1), copy->z(i2)); });
            else
                std::sort(indexvec.begin(), indexvec.end(), [&indexvec, copy] (int i1, int i2) { return complexGreaterThan(copy->z(i1), copy->z(i2)); });
        }

        // manage NaN indices
        size_t nanCount = 1;
        int nanIndex;
        int infIndex;

        if (ASCEND)
        {
            for (int i = 0; i < indexvec.size(); ++i)
            {
                // locate first NaN
                if (indexvec[i] == NaNidx[0])
                {
                    nanIndex = i;
                    break;
                }
            }

            for (int i = nanIndex; i < indexvec.size()-NaNidx.size(); ++i)
            {
                // locate next Inf
                infIndex = static_cast<int> (indexvec.size());

                for (int j = nanIndex+1; j < indexvec.size(); ++j)
                {
                    if (nanCount < NaNidx.size() && indexvec[j] == NaNidx[nanCount])
                    {
                        ++nanCount;
                        continue;
                    }

                    infIndex = j;
                    break;
                }

                if (infIndex == indexvec.size())
                    break;

                // locate next NaN
                nanIndex = static_cast<int> (indexvec.size());

                for (int j = infIndex+1; j < indexvec.size(); ++j)
                {
                    if (indexvec[j] == NaNidx[nanCount])
                    {
                        ++nanCount;
                        nanIndex = j;
                        break;
                    }
                }

                if (nanIndex == indexvec.size())
                    ++nanCount;     // pretend that another NaN was just beyond at the end

                // shift indices of NaN values
                for (int j = 0; j < (nanIndex-infIndex); ++j)
                    indexvec[i+j] = indexvec[i+j+nanCount-1];

                // jump ahead
                i += nanIndex-infIndex-1;
            }

            nanCount = 0;

            // restore NaN values and populate outputs
            for (size_t i = indexvec.size()-NaNidx.size(); i < indexvec.size(); ++i)
                indexvec[i] = NaNidx[nanCount++];
        }
        else    // DESCEND
        {
            for (int i = static_cast<int> (indexvec.size()-1); i > -1; --i)
            {
                // locate last NaN
                if (indexvec[i] == NaNidx[NaNidx.size()-1])
                {
                    nanIndex = i;
                    break;
                }
            }

            for (int i = nanIndex; i >= NaNidx.size(); --i)
            {
                // locate next Inf
                infIndex = -1;

                for (int j = nanIndex-1; j >= 0; --j)
                {
                    if (nanCount < NaNidx.size() && indexvec[j] == NaNidx[NaNidx.size()-1-nanCount])
                    {
                        ++nanCount;
                        continue;
                    }

                    infIndex = j;
                    break;
                }

                if (infIndex == -1)
                    break;

                // locate next NaN
                nanIndex = -1;

                for (int j = infIndex-1; j >= 0; --j)
                {
                    if (indexvec[j] == NaNidx[NaNidx.size()-1-nanCount])
                    {
                        ++nanCount;
                        nanIndex = j;
                        break;
                    }
                }

                if (nanIndex == -1)
                    ++nanCount;     // pretend that another NaN was just prior to the beginning

                // shift indices of NaN values
                for (int j = 0; j < (infIndex-nanIndex); ++j)
                    indexvec[i-j] = indexvec[i-j-nanCount+1];

                // jump back
                i -= infIndex-nanIndex-1;
            }

            nanCount = 0;

            // restore NaN values and populate outputs
            for (int i = 0; i < NaNidx.size(); ++i)
                indexvec[i] = NaNidx[nanCount++];
        }

        if (copy->IsReal())
        {
            for (size_t i = 0; i < NaNidx.size(); ++i)
                (*copy)(NaNidx[i]) = std::numeric_limits<double>::quiet_NaN();

            int indexIndex = index_data->first;
            hwMatrix* indices = index_data->second;
            for (int i = 0; i < result->Size(); ++i)
            {
                if (row)
                    (*indices)(indexIndex, i) = indexvec[i] + 1;
                else
                    (*indices)(i, indexIndex) = indexvec[i] + 1;
                (*result)(i) = (*copy)(indexvec[i]);
            }
            ++index_data->first;
        }
        else
        {
            for (size_t i = 0; i < NaNidx.size(); ++i)
                copy->z(NaNidx[i]) = std::numeric_limits<double>::quiet_NaN();

            int indexIndex = index_data->first;
            hwMatrix* indices = index_data->second;
            for (int i = 0; i < result->Size(); ++i)
            {
                if (row)
                    (*indices)(indexIndex, i) = indexvec[i] + 1;
                else
                    (*indices)(i, indexIndex) = indexvec[i] + 1;
                result->z(i) = copy->z(indexvec[i]);
            }
            ++index_data->first;
        }

        delete copy;
    }

    return result;
}
//------------------------------------------------------------------------------
// assumes all elements of cell are strings and that cell is a vector
//------------------------------------------------------------------------------
template <bool ASCEND>
HML_CELLARRAY* sort(EvaluatorInterface& eval, const HML_CELLARRAY* cell, std::pair<int&, hwMatrix*>* index_data)
{
    std::vector<hwMatrix*> rowvec;
    for (int i = 0; i < cell->Size(); ++i)
        rowvec.push_back(readRow(eval, (*cell)(i).Matrix()));

    std::vector<int> indexvec(rowvec.size());
    for (int i = 0 ; i < rowvec.size(); ++i)
        indexvec[i] = i;

    if (ASCEND)
        std::sort(indexvec.begin(), indexvec.end(), [&indexvec, &rowvec] (int i1, int i2) { return rowVecLessThan(rowvec[i1], rowvec[i2]); });
    else                                                      
        std::sort(indexvec.begin(), indexvec.end(), [&indexvec, &rowvec] (int i1, int i2) { return rowVecGreaterThan(rowvec[i1], rowvec[i2]); });
    
    bool row = cell->M() == 1;
    HML_CELLARRAY* result = EvaluatorInterface::allocateCellArray(row ? 1 : (int)rowvec.size(), row ? (int)rowvec.size() : 1);

    int indexIndex = index_data->first;
    hwMatrix* indices = index_data->second;
    for (int i = 0; i < result->Size(); ++i)
    {
        if (row)
            (*indices)(indexIndex, i) = indexvec[i] + 1;
        else
            (*indices)(i, indexIndex) = indexvec[i] + 1;
        (*result)(i) = rowvec[indexvec[i]];
    }
    ++index_data->first;
    return result;
}

// unique helpers for helpers

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void sph2cartHelper(double az, double elev, double r, double& x, double& y, double& z)
{
    x = r * cos(az) * sin(elev);
    y = r * sin(az) * sin(elev);
    z = r * cos(elev);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void cart2sphHelper(double x, double y, double z, double& az, double& elev, double& r)
{
    r    = sqrt(x*x + y*y + z*z);
    az   = atan2(y,x);
    elev = atan2(sqrt(x*x + y*y),z);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::vector<Currency> cart2sphHelper(EvaluatorInterface& eval, const std::vector<Currency>& inputs)
{
    if (!inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

    if (!inputs[1].IsScalar())
        throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

    if (!inputs[2].IsScalar())
        throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);

    std::vector<Currency> results;
    double x = inputs[0].Scalar();
    double y = inputs[1].Scalar();
    double z = inputs[2].Scalar();
    double az, elev, r;
    cart2sphHelper(x,y,z,az,elev,r);
    results.push_back(az);
    results.push_back(elev);
    results.push_back(r);

    return results;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::vector<Currency> sph2cartHelper(EvaluatorInterface& eval, const std::vector<Currency>& inputs)
{
    if (!inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

    if (!inputs[1].IsScalar())
        throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

    if (!inputs[2].IsScalar())
        throw OML_Error(OML_ERR_REAL, 3, OML_VAR_VALUE);

    std::vector<Currency> results;
    double x = inputs[0].Scalar();
    double y = inputs[1].Scalar();
    double z = inputs[2].Scalar();
    double az, elev, r;
    sph2cartHelper(x,y,z,az,elev,r);
    results.push_back(az);
    results.push_back(elev);
    results.push_back(r);

    return results;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::vector<Currency> pol2cartHelper(EvaluatorInterface& eval, const std::vector<Currency>& inputs)
{
    if (!inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

    if (!inputs[1].IsScalar())
        throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

    std::vector<Currency> results;
    std::pair<double,double> p = pol2cartHelper(inputs[0].Scalar(), inputs[1].Scalar());
    results.push_back(p.first);
    results.push_back(p.second);
    
    if (inputs.size() > 2)
        results.push_back(inputs[2]);

    return results;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::pair<double, double> pol2cartHelper(double theta, double r)
{
    return std::pair<double, double>(r*cos(theta), r*sin(theta));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::vector<Currency> cart2polHelper(EvaluatorInterface& eval, const std::vector<Currency>& inputs)
{
    if (!inputs[0].IsScalar())
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);

    if (!inputs[1].IsScalar())
        throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

    std::vector<Currency> results;
    std::pair<double,double> p = cart2polHelper(inputs[0].Scalar(), inputs[1].Scalar());
    results.push_back(p.first);
    results.push_back(p.second);

    if (inputs.size() > 2)
        results.push_back(inputs[2]);

    return results;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
std::pair<double, double> cart2polHelper(double x, double y)
{
    return std::pair<double, double>(atan2(y,x), sqrt(x*x + y*y));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void _countZeros(const hwMatrix *mtx, int *count, int stop_count, int stop_inner, int incr, bool horiz)
{
    for (; *count < stop_count; count += incr)
        for (int i = 0; i < stop_inner; i++)
            if (!iszero(horiz ? (*mtx)(*count, i) : (*mtx)(i, *count)))
                return;
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void _countZerosComplex(const hwMatrix *mtx, int *count, int stop_count, int stop_inner, int incr, bool horiz)
{
    for (; *count < stop_count; count += incr)
        for (int i = 0; i < stop_inner; i++)
            if (!iszero(horiz ? mtx->z(*count, i) : mtx->z(i, *count)))
                return;
}

// set-related helpers

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void dounique(std::deque<hwMatrix> &vals)
{
    std::sort(vals.begin(), vals.end(), static_cast<bool (*)(const hwMatrix&,const hwMatrix&)>(&rowVecLessThan));
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
}

void dounique(std::deque<double> &vals)
{
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
}

void dounique(std::deque<hwComplex> &vals)
{
    std::sort(vals.begin(), vals.end(), &complexLessThan);
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
}

void dounique(std::deque<std::string> &vals)
{
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
}

// all implementations of doSetFunc also rely on dounique to remove any duplicates

std::deque<hwMatrix> doSetFunc(std::deque<hwMatrix> &a, std::deque<hwMatrix> &b,
    std::back_insert_iterator<std::deque<hwMatrix> > (*stdmethod)(std::deque<hwMatrix>::iterator,
    std::deque<hwMatrix>::iterator, std::deque<hwMatrix>::iterator, std::deque<hwMatrix>::iterator,
    std::back_insert_iterator<std::deque<hwMatrix> >, bool (*)(const hwMatrix&, const hwMatrix&)))
{
    bool (*lessThan)(const hwMatrix&, const hwMatrix&) = static_cast<bool (*)(const hwMatrix&, const hwMatrix&)>(&rowVecLessThan);
    std::sort(a.begin(), a.end(), lessThan);
    std::sort(b.begin(), b.end(), lessThan);
    std::deque<hwMatrix> ret;
    (*stdmethod)(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(ret), lessThan);
    dounique(ret);
    return ret;
}

std::deque<hwComplex> doSetFunc(std::deque<hwComplex> &a, std::deque<hwComplex> &b,
    std::back_insert_iterator<std::deque<hwComplex> > (*stdmethod)(std::deque<hwComplex>::iterator,
    std::deque<hwComplex>::iterator, std::deque<hwComplex>::iterator, std::deque<hwComplex>::iterator,
    std::back_insert_iterator<std::deque<hwComplex> >, bool (*)(const hwComplex&, const hwComplex&)))
{
    std::sort(a.begin(), a.end(), &complexLessThan);
    std::sort(b.begin(), b.end(), &complexLessThan);
    std::deque<hwComplex> ret;
    (*stdmethod)(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(ret), &complexLessThan);
    dounique(ret);
    return ret;
}

std::deque<double> doSetFunc(std::deque<double> &a, std::deque<double> &b,
    std::back_insert_iterator<std::deque<double> > (*stdmethod)(std::deque<double>::iterator,
    std::deque<double>::iterator, std::deque<double>::iterator, std::deque<double>::iterator,
    std::back_insert_iterator<std::deque<double> >))
{
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    std::deque<double> ret;
    (*stdmethod)(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(ret));
    dounique(ret);
    return ret;
}

std::deque<std::string> doSetFunc(std::deque<std::string> &a, std::deque<std::string> &b,
    std::back_insert_iterator<std::deque<std::string> > (*stdmethod)(std::deque<std::string>::iterator,
    std::deque<std::string>::iterator, std::deque<std::string>::iterator, std::deque<std::string>::iterator,
    std::back_insert_iterator<std::deque<std::string> >))
{
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    std::deque<std::string> ret;
    (*stdmethod)(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(ret));
    dounique(ret);
    return ret;
}
/**
 * Returns the amount of CPU time used by the current process,
 * in seconds, or -1.0 if an error occurred.
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */
double getCPUTime()
{
#if defined(_WIN32)
	/* Windows -------------------------------------------------- */
	FILETIME createTime;
	FILETIME exitTime;
	FILETIME kernelTime;
	FILETIME userTime;
	if ( GetProcessTimes( GetCurrentProcess( ),
		&createTime, &exitTime, &kernelTime, &userTime ) != -1 )
	{
		SYSTEMTIME userSystemTime;
		if ( FileTimeToSystemTime( &userTime, &userSystemTime ) != -1 )
			return (double)userSystemTime.wHour * 3600.0 +
			(double)userSystemTime.wMinute * 60.0 +
			(double)userSystemTime.wSecond +
			(double)userSystemTime.wMilliseconds / 1000.0;
	}

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	/* AIX, BSD, Cygwin, HP-UX, Linux, OSX, and Solaris --------- */

#if defined(_POSIX_TIMERS) && (_POSIX_TIMERS > 0)
	/* Prefer high-res POSIX timers, when available. */
	{
		clockid_t id;
		struct timespec ts;
#if _POSIX_CPUTIME > 0
		/* Clock ids vary by OS.  Query the id, if possible. */
		if ( clock_getcpuclockid( 0, &id ) == -1 )
#endif
#if defined(CLOCK_PROCESS_CPUTIME_ID)
			/* Use known clock id for AIX, Linux, or Solaris. */
			id = CLOCK_PROCESS_CPUTIME_ID;
#elif defined(CLOCK_VIRTUAL)
			/* Use known clock id for BSD or HP-UX. */
			id = CLOCK_VIRTUAL;
#else
			id = (clockid_t)-1;
#endif
		if ( id != (clockid_t)-1 && clock_gettime( id, &ts ) != -1 )
			return (double)ts.tv_sec +
			(double)ts.tv_nsec / 1000000000.0;
	}
#endif

#if defined(RUSAGE_SELF)
	{
		struct rusage rusage;
		if ( getrusage( RUSAGE_SELF, &rusage ) != -1 )
			return (double)rusage.ru_utime.tv_sec +
			(double)rusage.ru_utime.tv_usec / 1000000.0;
	}
#endif

#if defined(_SC_CLK_TCK)
	{
		const double ticks = (double)sysconf( _SC_CLK_TCK );
		struct tms tms;
		if ( times( &tms ) != (clock_t)-1 )
			return (double)tms.tms_utime / ticks;
	}
#endif

#if defined(CLOCKS_PER_SEC)
	{
		clock_t cl = clock( );
		if ( cl != (clock_t)-1 )
			return (double)cl / (double)CLOCKS_PER_SEC;
	}
#endif

#endif

	return -1;		/* Failed. */
}
//------------------------------------------------------------------------------
//! Returns ordered vector of output options to display for regexp
//! \param[in] eval   Evaluator interface
//! \param[in] inputs Inputs
//------------------------------------------------------------------------------
std::vector<std::string> GetRegExpOptions(EvaluatorInterface           eval,
                                          const std::vector<Currency>& inputs)
{
    assert(!inputs.empty());

    std::vector<std::string> options;
    options.reserve(7);

    if (inputs.size() > 2)
    {
        // First currency in input is search string
        // Second currency in input is pattern string
        // Rest are options for outputs specified by user
        for (std::vector<Currency>::const_iterator itr = inputs.begin() + 2;
             itr != inputs.end(); ++itr)
        {
            std::string option (readString(toCurrencyStr(eval, *itr, false, false)));
            if (option.empty()) continue;

            std::transform(option.begin(), option.end(), option.begin(), ::tolower);

            if (option != "start"        &&
                option != "end"          &&
                option != "tokenextents" &&
                option != "match"        &&
                option != "tokens"       &&
                option != "names"        &&
                option != "split") OML_Error(HW_ERROR_INVALIDOPTION(option));

            options.push_back(option);
        }
    }

    // Now push in other options

    // Start indices of each matching substring
    if (options.empty() || 
        std::find(options.begin(), options.end(), "start") == options.end())
        options.push_back("start"); 

    // End indices of each matching substring
    if (std::find(options.begin(), options.end(), "end") == options.end())
        options.push_back("end");   

    // Extents of each matching token surrounded by (..) in pattern
    if (std::find(options.begin(), options.end(), "tokenextents") == options.end())
        options.push_back("tokenextents"); //

    // Cell array of the text of each matching string
    if (std::find(options.begin(), options.end(), "match") == options.end())
        options.push_back("match");

    // Cell array of each token matched
    if (std::find(options.begin(), options.end(), "tokens") == options.end())
        options.push_back("tokens");

    // Text of each matched named token
    if (std::find(options.begin(), options.end(), "names") == options.end())
        options.push_back("names");

    // Cell array of what is not returned by match
    if (std::find(options.begin(), options.end(), "split") == options.end())
        options.push_back("split");

    return options;
}
//------------------------------------------------------------------------------
//! Returns regular expression string matching results
//! \param[in] search  String to search in
//! \param[in] pattern Pattern to search for
//! \param[in] options Output options, in user specified order if applicable
//------------------------------------------------------------------------------
std::vector<Currency> DoRegExp(const std::string&              search, 
                               const std::string&              pattern, 
                               const std::vector<std::string>& options)
{
    std::string str = search;
    std::string pat = pattern;
    std::vector<Currency> outputs;
    std::vector<std::string> sub_names;

    // std::regex_search does not work with multiline strings in msvc10. The '.' 
    // in the pattern ignores carriage returns or new lines. So, modify the 
    // pattern if the search string is multilined
    if (!str.empty() && str.find("\n") != std::string::npos && !pattern.empty())
    {
        size_t pos = pattern.find(".*");
        if (pos != std::string::npos)
        {
            std::regex  r1(pattern);
            std::smatch m1;

            bool result = std::regex_search(str, m1, r1);
            if (!result)
            {
                pat = "";
                std::string::const_iterator itr = pattern.begin();
                for ( size_t i = 0; itr != pattern.end(); ++itr, ++i)
                {
                    char c = *itr;
                    if (c != '.' ||    
                       (c == '.' && i > 0 && pattern[i-1] == '\\')) // Literal
                    {
                        pat += c;
                        continue; 
                    }
                    pat += "(.|\r|\n)";
                }
            }
        }
    }

    size_t offset = 0;
    size_t sub_start = pat.find('(');
    while ((sub_start = pat.find('(', offset)) != std::string::npos)
    {
        offset = sub_start + 1;
        if (!isEscaped(pat, sub_start))
        {
            // doesn't allow empty names (e.g. ?<>)
            if (sub_start + 2 < pat.length())
            {
                if (pat[offset] == '?' && pat[++offset] == '<')
                {
                    offset = pat.find('>', ++offset);
                    // if we do this we should also escape \\
                    //while (isEscaped(pat, offset))
                    //{
                    //    pat.erase(pat.begin() + offset - 1);
                    //    offset = pat.find('>', offset);
                    //}

                    if (offset == std::string::npos)
                        throw OML_Error(HW_ERROR_SUBEXPNAMNOCLOSINGGR);

                    if (offset == sub_start + 3)
                        throw OML_Error(HW_ERROR_NOTEMPSUBEXPNAME);

                    std::string name = pat.substr(sub_start + 3, offset - (sub_start + 3));
                    sub_names.push_back(name);
                    pat.erase(pat.begin() + sub_start + 1, pat.begin() + offset + 1);
                    offset = sub_start + 1;
                    continue;
                }
            }
            sub_names.push_back(std::string());
        }
    }

    std::vector<std::string> match_text, unmatched_text;

    std::vector<double> startindicies; 
    std::vector<double> endindices;    

    // each outer container represents a matrix in the cell array output
    // each inner vector represents a row in said matrix
    std::vector<std::vector<std::pair<int, int> > > token_ranges;
    std::vector<std::vector<std::string> > token_text;
    std::map<std::string, std::vector<std::string> > token_names;

    try
    {
        int index = 0;

        std::regex reg(pat);
        std::smatch mtch;

        if (pat.length())
        {
            // find all matches
            while (std::regex_search(str, mtch, reg, std::regex_constants::match_default))
            {
                std::string matchstr = mtch.str();
                size_t matchpos = mtch.position();

                double start = static_cast<double>(index + matchpos + 1);
                double end   = static_cast<double>(index + matchpos + matchstr.length());

                startindicies.push_back(start);
                endindices.push_back(end);

                match_text.push_back(matchstr);
                unmatched_text.push_back(mtch.prefix().str());

                std::vector<std::pair<int, int> > matched_token_ranges;
                std::vector<std::string> matched_token_text;

                if (mtch.size() != sub_names.size() + 1)
                    throw OML_Error(HW_ERROR_NOTLOCALLSUBEXP);

                for (size_t i = 1; i < mtch.size(); i++)
                {
                    if (mtch[i].matched)
                    {
                        int token_start = index + (int)mtch.position(i);
                        matched_token_ranges.push_back(std::pair<int, int>(token_start + 1, token_start + (int) mtch[i].length()));
                        matched_token_text.push_back(mtch[i]);
                        if (sub_names[i-1].length())
                            token_names[sub_names[i-1]].push_back(mtch.str(i));
                    }
                    else if (sub_names[i-1].length())
                        token_names[sub_names[i-1]].push_back(std::string());
                }
                token_ranges.push_back(matched_token_ranges);
                token_text.push_back(matched_token_text);

                index += (int)(matchpos + matchstr.length());

                std::string updatedSearchStr (mtch.suffix().str());

                // Fix for memory leak if input contains '\n'
                if (updatedSearchStr == str) break; 

                str = updatedSearchStr;
            }
        }
        unmatched_text.push_back(str);
    }
    catch (std::regex_error& err)
    {
        BuiltInFuncsString::ThrowRegexError(err.code());
    }

    // create outputs
    HML_CELLARRAY *te = EvaluatorInterface::allocateCellArray(1, (int)token_ranges.size());
    for (size_t i = 0; i < token_ranges.size(); i++)
    {
        std::vector<std::pair<int, int> > sub_ranges = token_ranges[i];
        hwMatrix *range_matrix = EvaluatorInterface::allocateMatrix((int)sub_ranges.size(), 2, hwMatrix::REAL);
        for (size_t j = 0; j < sub_ranges.size(); j++)
        {
            std::pair<int, int> p = sub_ranges[j];
            (*range_matrix)((int)j, 0) = (double) p.first;
            (*range_matrix)((int)j, 1) = (double) p.second;
        }
        (*te)((int)i) = range_matrix;
    }
    HML_CELLARRAY *m = containerToCellArray(match_text, true);
    HML_CELLARRAY *t = EvaluatorInterface::allocateCellArray(1, (int)token_text.size());
    for (size_t i = 0; i < token_text.size(); i++)
    {
        std::vector<std::string> &sub_text = token_text[i];
        if (!sub_text.empty())
        {
            HML_CELLARRAY* textCell = containerToCellArray(sub_text, true);

            for (size_t j = 0; j < sub_text.size(); ++j)
                (*textCell)((int)j) = sub_text[j];
            
            (*t)((int)i) = textCell;
        }
        else
            (*t)((int)i) = EvaluatorInterface::allocateCellArray(1, 0);
    }

    // Option for printing names
    StructData *nm = new StructData();
    std::map<std::string, std::vector<std::string> >::const_iterator itr1 = token_names.begin();
    // if we care about the order of the field names of nm, we should iterate through
    // sub_names instead and check to make sure the size of sub_names[i] > 0
    for (; itr1 != token_names.end(); ++itr1)
    {
        std::string              name (itr1->first);
        std::vector<std::string> vec(itr1->second);
        if (vec.size() == 1)
        {
            nm->SetValue(0, 0, name, vec[0]);
            continue;
        }
        for (std::vector<std::string>::iterator itr2 = vec.begin();
             itr2 != vec.end();)
        {
            std::string val = *itr2;
            if (val.empty())
                itr2 = vec.erase(itr2); // Don't return empty elements
            else
                ++itr2;
        }

        HML_CELLARRAY *cell = containerToCellArray(vec, true);
        nm->SetValue(0, 0, name, cell);
    }
    HML_CELLARRAY *sp = containerToCellArray(unmatched_text, true);

    Currency tecur(te), mcur(m), tcur(t), spcur(sp);

    // deal with user-specified output order
    for (std::vector<std::string>::const_iterator itr = options.begin();
         itr != options.end(); ++itr)
    {
        std::string option (*itr);
       
        if (option == "start")
            outputs.push_back(Currency(startindicies));
        
        else if (option == "end")
            outputs.push_back(Currency(endindices));
        
        else if (option == "tokenextents")
            outputs.push_back(Currency(tecur));
        
        else if (option == "match")
            outputs.push_back(Currency(mcur));
        
        else if (option == "tokens")
            outputs.push_back(Currency(tcur));
        
        else if (option == "names")
            outputs.push_back(Currency(nm));
        
        else if (option == "split")
            outputs.push_back(Currency(spcur));
        
        else
            throw OML_Error(HW_ERROR_INVALIDOPTION(option)); // should be verified elsewhere
    }

    return outputs;
}
//------------------------------------------------------------------------------
//! Returns true if successul in regular expression string matching
//! \param[in] eval       Evaluator interface
//! \param[in] searchCur  Search currency
//! \param[in] patternCur Pattern currency
//! \param[in] outputs    Outputs
//------------------------------------------------------------------------------
bool GetRegExpOutput( EvaluatorInterface              eval,
                      const Currency&                 searchCur,
                      const Currency&                 patternCur,
                      const std::vector<std::string>& options,
                      std::vector<Currency>&          outputs)
{
    HML_CELLARRAY *cell1 = 0;
    HML_CELLARRAY *cell2 = 0;
    std::string str1, str2;

    Currency cur1 = toCurrencyStr(eval, searchCur, false, true);
    Currency cur2 = toCurrencyStr(eval, patternCur, false, true);

    if (cur1.IsString())
        str1 = readString(cur1);
    else if (cur1.IsCellArray())
        cell1 = cur1.CellArray();
    else
        throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

    if (cur2.IsString())
        str2 = readString(cur2);
    else if (cur2.IsCellArray())
        cell2 = cur2.CellArray();
    else
        throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

    if (cell1 && cell1->Size() == 1)
    {
        Currency temp = unnest(cur1, HW_ERROR_INPUTSTRINGCELLARRAY);
        if (temp.IsString())
        {
            str1 = readString(temp);
            cell1 = nullptr;
        }
        else
        {
            cur1 = temp;
            cell1 = cur1.CellArray();
        }
    }

    if (cell2 && cell2->Size() == 1)
    {
        Currency temp = unnest(cur2, HW_ERROR_INPUTSTRINGCELLARRAY);
        if (temp.IsString())
        {
            str2 = readString(temp);
            cell2 = nullptr;
        }
        else
        {
            cur2 = temp;
            cell2 = cur2.CellArray();
        }
    }

    if (!cell1 && !cell2) // There are only strings
    {
        outputs = DoRegExp(str1, str2, options);
        return true;
    }

    if (cell1)
    {
        if (cell2)
        {
            if (!sameSize(cell1, cell2))
                throw OML_Error(HW_ERROR_CELLARRAYSSAMESIZE);

            for (int i = 0; i < cell1->Size(); i++)
            {
                Currency elem1 = unnest((*cell1)(i), HW_ERROR_INPUTSTRINGCELLARRAY);
                Currency elem2 = unnest((*cell2)(i), HW_ERROR_INPUTSTRINGCELLARRAY);

                if (!(elem1.IsString() && elem2.IsString()))
                    throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

                std::vector<Currency> outvec = DoRegExp(readString(elem1), readString(elem2), options);

                if (!i)
                {
                    for (int j = 0; j < outvec.size(); j++)
                        outputs.push_back(EvaluatorInterface::allocateCellArray(cell1->M(), cell1->N()));
                }

                for (int j = 0; j < outvec.size(); j++)
                {
                    (*outputs[j].CellArray())(i) = outvec[j];
                }
            }
        }
        else
        {
            for (int i = 0; i < cell1->Size(); i++)
            {
                Currency elem1 = unnest((*cell1)(i), HW_ERROR_INPUTSTRINGCELLARRAY);
                if (!elem1.IsString())
                    throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

                std::vector<Currency> outvec = DoRegExp(readString(elem1), str2, options);

                if (!i)
                {
                    for (int j = 0; j < outvec.size(); j++)
                        outputs.push_back(EvaluatorInterface::allocateCellArray(cell1->M(), cell1->N()));
                }

                for (int j = 0; j < outvec.size(); j++)
                {
                    (*outputs[j].CellArray())(i) = outvec[j];
                }
            }
        }
    }
    else if (cell2)
    {
        for (int i = 0; i < cell2->Size(); i++)
        {
            Currency elem2 = unnest((*cell2)(i), HW_ERROR_INPUTSTRINGCELLARRAY);
            if (!elem2.IsString())
                throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

            std::vector<Currency> outvec = DoRegExp(str1, readString(elem2), options);

            if (!i)
            {
                for (int j = 0; j < outvec.size(); j++)
                    outputs.push_back(EvaluatorInterface::allocateCellArray(cell2->M(), cell2->N()));
            }

            for (int j = 0; j < outvec.size(); j++)
            {
                (*outputs[j].CellArray())(i) = outvec[j];
            }
        }
    }

    return true;
}
//------------------------------------------------------------------------------
//! Concatenates n-dimensional array objects along given dimension.
//! \param[in] eval     Evaluator interface
//! \param[in] inputs   Vector of inputs - string to be displayed to user
//! \param[out] outputs Vector of outputs - string returned by user
//------------------------------------------------------------------------------
bool oml_cat(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    size_t nargin = (inputs.empty()) ? 0 : inputs.size();

    if (nargin == 0)
        throw OML_Error(OML_ERR_NUMARGIN);

    // Get the concatenation dimension
    if (!inputs[0].IsPositiveInteger())
        throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_DIM);

    int dim = static_cast<int> (inputs[0].Scalar());

    if (nargin == 1)
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix(0, 0, hwMatrix::REAL));
        return true;
    }

    // determine dimensions and output object type: 2D or ND
    std::vector<int> dims;
	bool all_inputs_string = true;

    for (int i = 1; i < nargin; ++i)
    {
		const Currency& cur = inputs[i];

        if (cur.IsMatrix() || cur.IsScalar() || cur.IsComplex() || cur.IsString())
        {
			if (!cur.IsString())
				all_inputs_string = false;

            const hwMatrix* matrix = cur.ConvertToMatrix();
            int m = matrix->M();
            int n = matrix->N();

            if (!m && !n)
                continue;

            if (!dims.size())
            {
                dims.push_back(m);
                dims.push_back(n);

                if (nargin > 2)
                {
                    for (int i = 2; i < dim; ++i)
                        dims.push_back(1);
                }
            }
            else if (dim == 1)
            {
                dims[0] += m;

                if (dims[1] != n)
                    throw OML_Error(OML_ERR_ARRAYCATDIM, 2, i+1);
            }
            else if (dim == 2)
            {
                dims[1] += n;

                if (dims[0] != m)
                    throw OML_Error(OML_ERR_ARRAYCATDIM, 2, i+1);
            }
            else
            {
                if (dims[0] != m || dims[1] != n)
                    throw OML_Error(OML_ERR_ARRAYCATDIM, 2, i+1);

                for (int i = 2; i < dim-1; ++i)
                {
                    if (dims[i] != 1)
                        throw OML_Error(OML_ERR_ARRAYCATDIM, 2, i+1);
                }

                dims[dim-1] += 1;
            }
        }
        else if (cur.IsNDMatrix())
        {
            const hwMatrixN* matrix = cur.MatrixN();
            const std::vector<int>& curdims = matrix->Dimensions();
            int numDims = static_cast<int> (curdims.size());

            if (!dims.size())
            {
                dims = curdims;

                if (nargin > 2)
                {
                    for (int i = static_cast<int> (dims.size()); i < dim; ++i)
                        dims.push_back(1);
                }
            }
            else
            {
                for (int i = 0; i < numDims; ++i)
                {
                    if (i == dim-1)
                    {
                        dims[dim-1] += curdims[i];
                    }
                    else
                    {
                        if (dims[i] != curdims[i])
                            throw OML_Error(OML_ERR_ARRAYCATDIM, 2, i+1);
                    }
                }

                if (dim-1 >= numDims)
                    dims[dim-1] += 1;
            }
        }
        else
        {
            throw OML_Error(OML_ERR_MATRIX, i+1, OML_VAR_DIMS);
        }
    }

    bool NDout = false;

    if (!dims.size())
        dims.resize(2);
    else if (dims.size() > 2)
        NDout = true;

    // handle output object type cases
    if (!NDout)     // 2D case
    {
        hwMatrix* result = EvaluatorInterface::allocateMatrix(dims[0], dims[1], hwMatrix::REAL);

        if (dim == 1)
        {
            int m = 0;

            for (int i = 1; i < nargin; ++i)
            {
                const hwMatrix* matrix = inputs[i].ConvertToMatrix();
                hwMathStatus status = result->WriteSubmatrix(m, 0, *matrix);

                if (!status.IsOk())
                    throw OML_Error(status);

                m += matrix->M();
            }
        }
        else    // dim == 2
        {
            int n = 0;

            for (int i = 1; i < nargin; ++i)
            {
                const hwMatrix* matrix = inputs[i].ConvertToMatrix();
                hwMathStatus status = result->WriteSubmatrix(0, n, *matrix);

                if (!status.IsOk())
                    throw OML_Error(status);

                n += matrix->N();
            }
        }

        Currency ret(result);

        if (all_inputs_string)
            ret.SetMask(Currency::MASK_STRING);

        outputs.push_back(ret);
    }
    else    // ND case
    {
        hwMatrixN* result = EvaluatorInterface::allocateMatrixN();
        result->Dimension(dims, hwMatrixN::REAL);
        int dimCount = 0;

        for (int i = 1; i < nargin; ++i)
        {
            const Currency& cur = inputs[i];
            std::vector<hwSliceArg> sliceArgs;

            if (!cur.IsNDMatrix())
            {
                const hwMatrix* matrix = inputs[i].ConvertToMatrix();
                hwMatrixN temp;
                temp.Convert2DtoND(*matrix);
                const std::vector<int>& curDims = temp.Dimensions();
                const std::vector<int>& resultDims = result->Dimensions();

                for (int j = 0; j < resultDims.size(); ++j)
                {
                    if (j == dim-1)
                    {
                        if (dim-1 < curDims.size())
                        {
                            int length = curDims[dim-1];
                            std::vector<int> dimVec(length);

                            for (int k = 0; k < length; ++k)
                                dimVec[k] = dimCount++;

                            sliceArgs.push_back(dimVec);
                        }
                        else
                        {
                            sliceArgs.push_back(dimCount++);
                        }
                    }
                    else
                    {
                        sliceArgs.push_back(hwSliceArg());
                    }
                }

                try
                {
                    result->SliceLHS(sliceArgs, temp);
                }
                catch (hwMathException& )
                {
                    throw OML_Error(OML_ERR_ARRAYCATDIM, 2, i+1);
                }
            }
            else    // cur.IsNDMatrix()
            {
                const hwMatrixN* matrixN = inputs[i].MatrixN();
                const std::vector<int>& curDims = matrixN->Dimensions();
                const std::vector<int>& resultDims = result->Dimensions();

                for (int j = 0; j < resultDims.size(); ++j)
                {
                    if (j == dim-1)
                    {
                        if (dim-1 < curDims.size())
                        {
                            int length = curDims[dim-1];
                            std::vector<int> dimVec(length);

                            for (int k = 0; k < length; ++k)
                                dimVec[k] = dimCount++;

                            sliceArgs.push_back(dimVec);
                        }
                        else
                        {
                            sliceArgs.push_back(dimCount++);
                        }
                    }
                    else
                    {
                        sliceArgs.push_back(hwSliceArg());
                    }
                }

                try
                {
                    result->SliceLHS(sliceArgs, *matrixN);
                }
                catch (hwMathException& )
                {
                    throw OML_Error(OML_ERR_ARRAYCATDIM, 2, i+1);
                }
            }
        }

        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
//! Helper function for unique command
//! \param[in]  eval      Evaluator interface
//! \param[in]  x         Input currency
//! \param[in]  cmprows   True if whole rows are compared
//! \param[in]  forward   True if forward search
//! \param[in]  outputIdx True if idx(i) of each element of output(y) in given 
//!                       curr(x) is returned such that y = x(i)
//! \param[in]  inputIdx  True if idx(j) of each element of input(x) in 
//!                       output(y) such that x = y(j)
//! \param[out] outputs   Vector of outputs
//------------------------------------------------------------------------------
void UniqueHelperFunc(EvaluatorInterface&    eval,
                      const Currency&        x,
                      bool                   cmprows,
                      bool                   forward,
                      bool                   outputIdx,
                      bool                   inputIdx,
                      std::vector<Currency>& outputs)
{
    if (x.IsScalar() || x.IsComplex())
    {
        outputs.push_back(x);
        if (outputIdx)
            outputs.push_back(1);
        if (inputIdx)
            outputs.push_back(1);
        return;
    }
    
    if (x.IsMatrix() || x.IsString())      // Matrix
    {
        UniqueHelperFuncMtx(eval, x, cmprows, forward, outputIdx, inputIdx, outputs);
        return;
    }
    
    if (!x.IsCellArray() && !x.IsNDMatrix()) 
        throw OML_Error(HW_ERROR_INPUTSTRCELLMTX);

    std::string rwarn ("Warning: Rows can be compared only for 2-D matrices;");
    rwarn += " ignoring option 'rows'";

    if (cmprows)
        BuiltInFuncsUtils::SetWarning(eval, rwarn);

    if (x.IsCellArray())                   // Cell
        UniqueHelperFuncCell(eval, x, forward, outputIdx, inputIdx, outputs);
    
    else if (x.IsNDMatrix())
        BuiltInFuncsElemMath::UniqueHelperFuncMtxN(eval, x, forward, 
                              outputIdx, inputIdx, outputs);        
}
//------------------------------------------------------------------------------
//! Helper function for unique command for 2D matrices or strings
//! \param[in]  eval      Evaluator interface
//! \param[in]  x         Input currency
//! \param[in]  cmprows   True if whole rows are compared
//! \param[in]  forward   True if forward search
//! \param[in]  outputIdx True if idx(i) of each element of output(y) in given 
//!                       curr(x) is returned such that y = x(i)
//! \param[in]  inputIdx  True if idx(j) of each element of input(x) in 
//!                       output(y) such that x = y(j)
//! \param[out] outputs   Vector of outputs
//------------------------------------------------------------------------------
void UniqueHelperFuncMtx(EvaluatorInterface&    eval,
                         const Currency&        x,
                         bool                   cmprows,
                         bool                   forward,
                         bool                   outputIdx,
                         bool                   inputIdx,
                         std::vector<Currency>& outputs) 
{
    const hwMatrix* mtx = x.Matrix();
    if (!mtx || mtx->Size() <= 1)  // Handle empty matrices
    {
        BuiltInFuncsElemMath::UniqueHelperFuncMtx(eval, x, forward, outputIdx, 
                                                  inputIdx, outputs);        
        return;
    }

    int rows = mtx->M();
    int cols = mtx->N();
    
    if (cmprows && cols != 1)
    {
        std::deque<hwMatrix> vals;
        for (int i = 0; i < rows; ++i)
        {
            Currency rowcur = BuiltInFuncsUtils::ReadRow(eval, x, i);
            const hwMatrix* rowmtx = rowcur.Matrix();
            if (rowmtx)
                vals.push_back(*rowmtx);
        }
        dounique(vals);
        Currency out = dequeToMatrix(vals);
        if (x.IsString())
            out.SetMask(Currency::MASK_STRING);
        outputs.push_back(out);

        if (outputIdx)
            outputs.push_back(findUnsorted(eval, vals, mtx, forward));

        if (inputIdx)
            outputs.push_back(findSorted(eval, vals, mtx));
        return;
    }
    else
        BuiltInFuncsElemMath::UniqueHelperFuncMtx(eval, x, forward, outputIdx, 
                                                  inputIdx, outputs);
}
//------------------------------------------------------------------------------
//! Helper function for unique command for cell arrays
//! \param[in]  eval      Evaluator interface
//! \param[in]  x         Input currency
//! \param[in]  cmprows   True if whole rows are compared
//! \param[in]  forward   True if forward search
//! \param[in]  outputIdx True if idx(i) of each element of output(y) in given 
//!                       curr(x) is returned such that y = x(i)
//! \param[in]  inputIdx  True if idx(j) of each element of input(x) in 
//!                       output(y) such that x = y(j)
//! \param[out] outputs   Vector of outputs
//------------------------------------------------------------------------------
void UniqueHelperFuncCell(EvaluatorInterface&    eval,
                          const Currency&        x,
                          bool                   forward,
                          bool                   outputIdx,
                          bool                   inputIdx,
                          std::vector<Currency>& outputs)
{
    HML_CELLARRAY* cell = x.CellArray();
    if (!cell || cell->IsEmpty())
    {
        outputs.push_back(EvaluatorInterface::allocateCellArray());
        if (outputIdx)
            outputs.push_back(EvaluatorInterface::allocateMatrix());
        if (inputIdx)
            outputs.push_back(EvaluatorInterface::allocateMatrix());
        return;
    }

    int cellsize = cell->Size();
    if (!isstr(cell)) throw OML_Error(HW_ERROR_FIELDNAMECELLSTR);

    std::deque<std::string> vals;
    for (int i = 0; i < cellsize; ++i)
        vals.push_back(readString((*cell)(i)));
    
    dounique(vals);
    outputs.push_back(containerToCellArray(vals, cell->M() == 1));

    if (outputIdx)
        outputs.push_back(findUnsorted(eval, vals, cell, forward));

    if (inputIdx)
        outputs.push_back(findSorted(eval, vals, cell));           
}
//------------------------------------------------------------------------------
//! Dummy method, returning 1 if needed - implemented in GUI only
//------------------------------------------------------------------------------
bool dummyint(EvaluatorInterface           eval,
              const std::vector<Currency>& inputs,
              std::vector<Currency>&       outputs)
{ 
    if (eval.GetNargoutValue() > 0)
        outputs.push_back("1");

    return true; 
}
//------------------------------------------------------------------------------
//! Dummy method returning empty string if needed - implemented in GUI only
//------------------------------------------------------------------------------
bool dummystring(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs)
{ 
    int numoutputs = eval.GetNargoutValue();
    if (numoutputs > 0)
        outputs.push_back("");

    if (numoutputs > 1)
        outputs.push_back("");

    return true; 
}
//------------------------------------------------------------------------------
//! Dummy method returning empty cell if needed - implemented in GUI only
//------------------------------------------------------------------------------
bool dummycell(EvaluatorInterface           eval,
               const std::vector<Currency>& inputs,
               std::vector<Currency>&       outputs)
{ 
    if (eval.GetNargoutValue())
        outputs.push_back(EvaluatorInterface::allocateCellArray());

    return true; 
}
//------------------------------------------------------------------------------
// Throws an error without stack information [erromsgonly]
//------------------------------------------------------------------------------
bool oml_erromsgonly(EvaluatorInterface           eval, 
                     const std::vector<Currency>& inputs, 
                     std::vector<Currency>&       outputs)
{
    int nargin = (inputs.empty()) ? 0 : static_cast<int>(inputs.size());
    if (nargin < 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    Currency cur = inputs[0];
    if (!cur.IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }
    std::string msg = (nargin == 1) ? cur.StringVal() : 
                      sprintf(eval, inputs.cbegin(), inputs.cend());
    if (msg.empty())
    {
        return true;
    }

    std::string lower (msg);
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower.find("error: ", 0, 7) == std::string::npos)
    {
        msg.insert(0, "Error: ");
    }

    throw OML_Error(msg, false); // Don't add stack information to message
    return true;
}
