﻿/**
* @file BuiltInFuncs.cpp
* @date October 2013
* Copyright (C) 2013-2021 Altair Engineering, Inc.  
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
#if OS_WIN
#	if _MSC_VER >= 1900
#define _CRT_NO_VA_START_VALIDATION
#   endif
#endif

#ifdef OS_WIN
#    define NOMINMAX
#endif

#include "BuiltInFuncs.h"

#include <algorithm>
#include <cassert>
#include <chrono>
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
#include "BuiltInFuncsMKL.h"
#include "BuiltInFuncsString.h"
#include "BuiltInFuncsSystem.h"
#include "BuiltInFuncsTime.h"
#include "BuiltInFuncsUtils.h"
#include "CurrencyDisplay.h"
#include "MemoryScope.h"
#include "Interpreter.h"
#include "FunctionInfo.h"
#include "OMLTree.h"
#include "ErrorInfo.h"
#include "OML_Error.h"
#include "hwComplex.h"
#include "StructData.h"
#include "MatrixNUtils.h"
#include "SignalHandlerBase.h"
#include "utf8utils.h"

#include <math.h>
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

// Regular expression matching helpers

// Enum for regular expression matching outputs order
enum OMLREGEXP
{
    OMLREGEXP_START,       // Start indicies for matches
    OMLREGEXP_END,         // End indicies for matches
    OMLREGEXP_EXT,         // Extent of each matching token in (..)
    OMLREGEXP_MATCH,       // Matches
    OMLREGEXP_TOK,         // Text of matched token
    OMLREGEXP_NAMES,       // Named tokens
    OMLREGEXP_SPLIT        // Remainder of unmatched token
};
// Returns regular expression string matching results
std::vector<Currency> DoRegExp(EvaluatorInterface&             eval,
                               const std::string&              search,
                               const std::string&              pattern,
                               const std::vector<OMLREGEXP>&   options,
	                           const std::vector<std::string>& flags);
// Returns regular expression string matching results
void DoRegExp(EvaluatorInterface&, 
    HML_CELLARRAY*, 
    const std::string&,
    const std::vector<OMLREGEXP>&, 
    const std::vector<std::string>&, 
    int, 
    int,
    int,
    std::vector<Currency>&);

// Gets ordered vector of output options to display for regexp
std::vector<OMLREGEXP> GetRegexOptions(const std::vector<Currency>& inputs,
                                       int                          nargout,
	                                   std::vector<std::string>&    flags);
// Returns true if successful in regular expression string matching
void GetRegExpOutput(EvaluatorInterface              eval,
                     const Currency&                 searchCur,
                     const Currency&                 patternCur,
                     const std::vector<OMLREGEXP>&   options,
	                 const std::vector<std::string>& flags,
                     std::vector<Currency>&          outputs);
// Returns the first element in a nested cell array
Currency Unnest(Currency, omlMathErrCode, int idx = -1);
// Helper method to read from file
std::vector<Currency> FreadData(EvaluatorInterface, DataType, bool, int, 
    size_t, int, int, std::FILE*, int, int);
// Helper method to read signed/unsigned character block from file
bool ReadCharBlock(EvaluatorInterface, std::FILE*, bool, bool, size_t, int,
    hwMatrix*, int&, int&);
// Helper method to read double block from file
bool ReadDoubleBlock(EvaluatorInterface, std::FILE*, bool, size_t, int,
    hwMatrix*, int&, int&);
// Helper method to read float block from file
bool ReadFloatBlock(EvaluatorInterface, std::FILE*, bool, size_t, int,
    hwMatrix*, int&, int&);
// Helper method to read signed/unsigned int block from file
bool ReadIntBlock(EvaluatorInterface, std::FILE*, bool, bool, size_t, int,
    hwMatrix*, int&, int&);
// Helper method to read signed/unsigned int8_t block from file
bool ReadInt8Block(EvaluatorInterface, std::FILE*, bool, bool, size_t, int,
    hwMatrix*, int&, int&);
// Helper method to read signed/unsigned int16_t block from file
bool ReadInt16Block(EvaluatorInterface, std::FILE*, bool, bool, size_t, int,
    hwMatrix*, int&, int&);
// Helper method to read signed/unsigned int32_t block from file
bool ReadInt32Block(EvaluatorInterface, std::FILE*, bool, bool, size_t, int,
    hwMatrix*, int&, int&);
// Helper method to read signed/unsigned int64_t block from file
bool ReadInt64Block(EvaluatorInterface, std::FILE*, bool, bool, size_t, int,
    hwMatrix*, int&, int&);

// Helper method to read signed/unsigned long block from file
bool ReadLongBlock(EvaluatorInterface, std::FILE*, bool, bool, size_t, int,
    hwMatrix*, int&, int&);
// Helper method to read signed/unsigned long long block from file
bool ReadLongLongBlock(EvaluatorInterface, std::FILE*, bool, bool, size_t, int,
    hwMatrix*, int&, int&);
// Helper method to read signed/unsigned short block from file
bool ReadShortBlock(EvaluatorInterface, std::FILE*, bool, bool, size_t, int,
    hwMatrix*, int&, int&);

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
    (*std_functions)["isdiag"]             = BuiltinFunc(oml_isdiag, FunctionMetaData(1, 1, LINA));
    (*std_functions)["istril"]             = BuiltinFunc(oml_istril, FunctionMetaData(1, 1, LINA));
    (*std_functions)["istriu"]             = BuiltinFunc(oml_istriu, FunctionMetaData(1, 1, LINA));
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
    (*std_functions)["isfield"]            = BuiltinFunc(oml_isfield, FunctionMetaData(2, 1, DATA));
    (*std_functions)["which"]              = BuiltinFunc(oml_which, FunctionMetaData(-1, -1, CORE));
    (*std_functions)["linspace"]           = BuiltinFunc(oml_linspace, FunctionMetaData(3, 1, CORE));
    (*std_functions)["isglobal"]           = BuiltinFunc(oml_isglobal, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["printf"]             = BuiltinFunc(oml_printf, FunctionMetaData(-2, 1, FILEIO));
    (*std_functions)["fprintf"]            = BuiltinFunc(oml_fprintf, FunctionMetaData(-3, 1, FILEIO));
    (*std_functions)["str2num"]            = BuiltinFunc(oml_str2num, FunctionMetaData(1, 2, STNG));
    (*std_functions)["fieldnames"]         = BuiltinFunc(oml_fieldnames, FunctionMetaData(1, 1, DATA));
    (*std_functions)["sprintf"]            = BuiltinFunc(oml_sprintf, FunctionMetaData(-2, 2, STNG));
    (*std_functions)["feval"]              = BuiltinFunc(oml_feval, FunctionMetaData(-2, -1, CORE));
    (*std_functions)["rmpath"]             = BuiltinFunc(oml_rmpath, FunctionMetaData(-1, 1, CORE));
	(*std_functions)["rmhiddenpath"]       = BuiltinFunc(oml_rmhiddenpath, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["addpath"]            = BuiltinFunc(oml_addpath, FunctionMetaData(-1, 1, CORE));
	(*std_functions)["addhiddenpath"]      = BuiltinFunc(oml_addhiddenpath, FunctionMetaData(-1, 1, CORE));
	(*std_functions)["registerpath"]       = BuiltinFunc(oml_addpath2, FunctionMetaData(2, 1, CORE));
    (*std_functions)["path"]               = BuiltinFunc(oml_path, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["repmat"]             = BuiltinFunc(oml_repmat, FunctionMetaData(-2, 1, ELEM));
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
	(*std_functions)["lcm"]                = BuiltinFunc(oml_lcm, FunctionMetaData(1, 2, ELEM));
    (*std_functions)["polyval"]            = BuiltinFunc(oml_polyval, FunctionMetaData(4, 1, POLY));
    (*std_functions)["triu"]               = BuiltinFunc(oml_triu, FunctionMetaData(2, 1, LINA));
    (*std_functions)["tril"]               = BuiltinFunc(oml_tril, FunctionMetaData(2, 1, LINA));
    (*std_functions)["iscellstr"]          = BuiltinFunc(oml_iscellstr, FunctionMetaData(1, 1, DATA));   
    (*std_functions)["poly"]               = BuiltinFunc(oml_poly, FunctionMetaData(1, 1, POLY));
    (*std_functions)["reshape"]            = BuiltinFunc(oml_reshape, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["permute"]            = BuiltinFunc(oml_permute, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["ipermute"]           = BuiltinFunc(oml_ipermute, FunctionMetaData(2, 1, ELEM));
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
	(*std_functions)["ctranspose"]         = BuiltinFunc(oml_ctranspose, FunctionMetaData(1, 1, LINA));
    (*std_functions)["eye"]                = BuiltinFunc(oml_eye, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["sign"]               = BuiltinFunc(oml_sign, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["pow2"]               = BuiltinFunc(oml_pow2, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["dot"]                = BuiltinFunc(oml_dot, FunctionMetaData(3, 1, LINA));
    (*std_functions)["kron"]               = BuiltinFunc(oml_kron, FunctionMetaData(2, 1, LINA));
    (*std_functions)["diff"]               = BuiltinFunc(oml_diff, FunctionMetaData(3, 1, ELEM));
    (*std_functions)["diag"]               = BuiltinFunc(oml_diag, FunctionMetaData(-2, 1, ELEM));
    (*std_functions)["sparse"]             = BuiltinFunc(oml_sparse, FunctionMetaData(-2, 1, LINA));
    (*std_functions)["issparse"]           = BuiltinFunc(oml_issparse, FunctionMetaData(1, 1, LINA));
    (*std_functions)["nnz"]                = BuiltinFunc(oml_nnz, FunctionMetaData(1, 1, LINA));
    (*std_functions)["speye"]              = BuiltinFunc(oml_speye, FunctionMetaData(-2, 1, LINA));
    (*std_functions)["spones"]             = BuiltinFunc(oml_spones, FunctionMetaData(2, 1, LINA));
    (*std_functions)["full"]               = BuiltinFunc(oml_full, FunctionMetaData(1, 1, LINA));
    (*std_functions)["conv"]               = BuiltinFunc(oml_conv, FunctionMetaData(3, 1, LINA));
    (*std_functions)["e"]                  = BuiltinFunc(oml_e, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["pi"]                 = BuiltinFunc(oml_pi, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["ones"]               = BuiltinFunc(oml_ones, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["zeros"]              = BuiltinFunc(oml_zeros, FunctionMetaData(-1, 1, ELEM));
    (*std_functions)["trace"]              = BuiltinFunc(oml_trace, FunctionMetaData(1, 2, LINA));
    (*std_functions)["det"]                = BuiltinFunc(oml_det, FunctionMetaData(1, 2, LINA));
    (*std_functions)["rcond"]              = BuiltinFunc(oml_rcond, FunctionMetaData(1, 1, LINA));
    (*std_functions)["real"]               = BuiltinFunc(oml_real, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["imag"]               = BuiltinFunc(oml_imag, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["conj"]               = BuiltinFunc(oml_conj, FunctionMetaData(1, 1, LINA));
    (*std_functions)["sum"]                = BuiltinFunc(oml_sum, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["cumsum"]             = BuiltinFunc(oml_cumsum, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["prod"]               = BuiltinFunc(oml_prod, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["cumprod"]            = BuiltinFunc(oml_cumprod, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["accumarray"]         = BuiltinFunc(oml_accumarray, FunctionMetaData(-3, 1, ELEM));
    (*std_functions)["deg2rad"]            = BuiltinFunc(oml_deg2rad, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["rad2deg"]            = BuiltinFunc(oml_rad2deg, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["cos"]                = BuiltinFunc(BuiltInFuncsMKL::Cos, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["cosd"]               = BuiltinFunc(BuiltInFuncsMKL::CosD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["sin"]                = BuiltinFunc(BuiltInFuncsMKL::Sin, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["sind"]               = BuiltinFunc(BuiltInFuncsMKL::SinD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["tan"]                = BuiltinFunc(BuiltInFuncsMKL::Tan, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["tand"]               = BuiltinFunc(BuiltInFuncsMKL::TanD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["sec"]                = BuiltinFunc(BuiltInFuncsMKL::Sec, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["secd"]               = BuiltinFunc(BuiltInFuncsMKL::SecD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["csc"]                = BuiltinFunc(BuiltInFuncsMKL::Csc, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["cscd"]               = BuiltinFunc(BuiltInFuncsMKL::CscD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["cot"]                = BuiltinFunc(BuiltInFuncsMKL::Cot, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["cotd"]               = BuiltinFunc(BuiltInFuncsMKL::CotD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["acos"]               = BuiltinFunc(BuiltInFuncsMKL::aCos, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["acosd"]              = BuiltinFunc(BuiltInFuncsMKL::aCosD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["asin"]               = BuiltinFunc(BuiltInFuncsMKL::aSin, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["asind"]              = BuiltinFunc(BuiltInFuncsMKL::aSinD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["atan"]               = BuiltinFunc(BuiltInFuncsMKL::aTan, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["atand"]              = BuiltinFunc(BuiltInFuncsMKL::aTanD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["asec"]               = BuiltinFunc(BuiltInFuncsMKL::aSec, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["asecd"]              = BuiltinFunc(BuiltInFuncsMKL::aSecD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["acsc"]               = BuiltinFunc(BuiltInFuncsMKL::aCsc, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["acscd"]              = BuiltinFunc(BuiltInFuncsMKL::aCscD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["acot"]               = BuiltinFunc(BuiltInFuncsMKL::aCot, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["acotd"]              = BuiltinFunc(BuiltInFuncsMKL::aCotD, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["atan2"]              = BuiltinFunc(BuiltInFuncsMKL::aTan2, FunctionMetaData(2, 1, TRIG));
    (*std_functions)["atan2d"]             = BuiltinFunc(BuiltInFuncsMKL::aTan2D, FunctionMetaData(2, 1, TRIG));
    (*std_functions)["cosh"]               = BuiltinFunc(BuiltInFuncsMKL::Cosh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["sinh"]               = BuiltinFunc(BuiltInFuncsMKL::Sinh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["tanh"]               = BuiltinFunc(BuiltInFuncsMKL::Tanh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["acosh"]              = BuiltinFunc(BuiltInFuncsMKL::aCosh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["asinh"]              = BuiltinFunc(BuiltInFuncsMKL::aSinh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["atanh"]              = BuiltinFunc(BuiltInFuncsMKL::aTanh, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["log10"]              = BuiltinFunc(BuiltInFuncsMKL::Log10, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["log2"]               = BuiltinFunc(BuiltInFuncsMKL::Log2, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["log"]                = BuiltinFunc(BuiltInFuncsMKL::Log, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["exp"]                = BuiltinFunc(BuiltInFuncsMKL::Exp, FunctionMetaData(1, 1, TRIG));
    (*std_functions)["nextpow2"]           = BuiltinFunc(oml_nextpow2, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["sqrt"]               = BuiltinFunc(BuiltInFuncsMKL::Sqrt, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["abs"]                = BuiltinFunc(BuiltInFuncsMKL::Abs, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["arg"]                = BuiltinFunc(BuiltInFuncsMKL::Arg, FunctionMetaData(1, 1, CORE));
    (*std_functions)["angle"]              = BuiltinFunc(BuiltInFuncsMKL::Arg, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["round"]              = BuiltinFunc(oml_round, FunctionMetaData(-2, 1, ELEM));
    (*std_functions)["ceil"]               = BuiltinFunc(BuiltInFuncsMKL::Ceil, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["floor"]              = BuiltinFunc(BuiltInFuncsMKL::Floor, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["fix"]                = BuiltinFunc(BuiltInFuncsMKL::Fix, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["rem"]                = BuiltinFunc(BuiltInFuncsMKL::Rem, FunctionMetaData(2, 1, ELEM));
    (*std_functions)["mod"]                = BuiltinFunc(BuiltInFuncsMKL::Mod, FunctionMetaData(2, 2, ELEM));
    (*std_functions)["max"]                = BuiltinFunc(BuiltInFuncsMKL::Max, FunctionMetaData(3, 2, ELEM));
    (*std_functions)["min"]                = BuiltinFunc(BuiltInFuncsMKL::Min, FunctionMetaData(3, 2, ELEM));
    (*std_functions)["disp"]               = BuiltinFunc(oml_print, FunctionMetaData(1, 0, CORE));
    (*std_functions)["inv"]                = BuiltinFunc(oml_inv, FunctionMetaData(1, 2, LINA));
    (*std_functions)["logical"]            = BuiltinFunc(oml_logical, FunctionMetaData(1, 1, ELEM));
    (*std_functions)["clear"]              = BuiltinFunc(oml_clear, FunctionMetaData(-1, 0, CORE));
    (*std_functions)["end"]                = BuiltinFunc(oml_end, FunctionMetaData(OML_NO_NARG, OML_NO_NARG, CORE));
    (*std_functions)["nargout"]            = BuiltinFunc(oml_nargout, FunctionMetaData(1, 1, CORE));
    (*std_functions)["nargin"]             = BuiltinFunc(oml_nargin, FunctionMetaData(1, 1, CORE));
    (*std_functions)["who"]                = BuiltinFunc(oml_who, FunctionMetaData(-1, -1, CORE));
    (*std_functions)["true"]               = BuiltinFunc(oml_true, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["false"]              = BuiltinFunc(oml_false, FunctionMetaData(-1, 1, CORE));
    (*std_functions)["refcnt"]             = BuiltinFunc(oml_refcnt,FunctionMetaData( 1, 1, CORE));
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
	(*std_functions)["cputime"]            = BuiltinFunc(oml_cputime, FunctionMetaData(0, 1, SYSTEM));
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
    (*std_functions)["de2bi"]     = BuiltinFunc(BuiltInFuncsConvert::De2Bi, FunctionMetaData(1, -4, ELEM));
    (*std_functions)["bi2de"]     = BuiltinFunc(BuiltInFuncsConvert::Bi2De, FunctionMetaData(1, -3, ELEM));
    (*std_functions)["vec2mat"]     = BuiltinFunc(BuiltInFuncsConvert::Vec2Mat, FunctionMetaData(-2, -3, ELEM));
    (*std_functions)["single"] = BuiltinFunc(BuiltInFuncsElemMath::Single, FunctionMetaData(1, 0, ELEM));

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
	(*std_functions)["addlibrary"] = BuiltinFunc(BuiltInFuncsCore::AddToolbox, FunctionMetaData(1, 0, CORE));
	(*std_functions)["registerlibrary"] = BuiltinFunc(BuiltInFuncsCore::RegisterLibrary, FunctionMetaData(2, 0, CORE));
	(*std_functions)["removelibrary"] = BuiltinFunc(BuiltInFuncsCore::RemoveToolbox, FunctionMetaData(1, 0, CORE));
    (*std_functions)["funccount"]  = BuiltinFunc(BuiltInFuncsCore::Funccount,   FunctionMetaData(0, 1, CORE));
    (*std_functions)["funclist"]   = BuiltinFunc(BuiltInFuncsCore::Funclist,  FunctionMetaData(-1, 0, CORE));
    (*std_functions)["varlist"]    = BuiltinFunc(BuiltInFuncsCore::Varlist,   FunctionMetaData(0, 1, CORE));
	(*std_functions)["keywordlist"] = BuiltinFunc(BuiltInFuncsCore::Keywordlist, FunctionMetaData(0, 1, CORE));
	(*std_functions)["omlfilename"] = BuiltinFunc(BuiltInFuncsCore::Omlfilename, FunctionMetaData(1, 1, FILEIO));
	(*std_functions)["omllinenumber"] = BuiltinFunc(BuiltInFuncsCore::Omllinenumber, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["exist"]       = BuiltinFunc(BuiltInFuncsCore::Exist, FunctionMetaData(2, 1, CORE));
    (*std_functions)["type"]        = BuiltinFunc(BuiltInFuncsCore::Type,  FunctionMetaData(-1, 1, CORE));
    (*std_functions)["omlpaginate"] = BuiltinFunc(BuiltInFuncsCore::OmlPaginate, FunctionMetaData(1, 0, CORE));
    (*std_functions)["sleep"]       = BuiltinFunc(BuiltInFuncsCore::OmlSleep, FunctionMetaData(1, 0, TIME));
    (*std_functions)["errormsgonly"] = BuiltinFunc(oml_erromsgonly, FunctionMetaData(-2, 0, CORE));
    (*std_functions)["warningmsgonly"] = BuiltinFunc(oml_warningmsgonly, FunctionMetaData(-1, 0, CORE));
    (*std_functions)["clc"]            = BuiltinFunc(BuiltInFuncsCore::Clc, FunctionMetaData(0, 0, CORE));
    (*std_functions)["display"] = BuiltinFunc(BuiltInFuncsCore::Display, FunctionMetaData(1, 0, CORE));
    (*std_functions)["skipformat"] = BuiltinFunc(BuiltInFuncsCore::SkipFormat, FunctionMetaData(1, 0, CORE));
    (*std_functions)["exit"]           = BuiltinFunc(BuiltInFuncsCore::Exit, FunctionMetaData(-1, 0, CORE));
    (*std_functions)["quit"]           = BuiltinFunc(BuiltInFuncsCore::Exit, FunctionMetaData(-1, 0, CORE));
	(*std_functions)["properties"] = BuiltinFunc(oml_properties, FunctionMetaData(-1, 0, CORE));
	(*std_functions)["methods"] = BuiltinFunc(oml_methods, FunctionMetaData(-1, 0, CORE));

    // Data structures
    (*std_functions)["cat"]       = BuiltinFunc(oml_cat, FunctionMetaData(-2, 1, DATA));
    (*std_functions)["setfield"]  = BuiltinFunc(BuiltInFuncsData::Setfield, FunctionMetaData(-3, 1, DATA));
	(*std_functions)["mat2cell"]  = BuiltinFunc(BuiltInFuncsData::Mat2Cell, FunctionMetaData(-2, 1, DATA));
	(*std_functions)["num2cell"]  = BuiltinFunc(BuiltInFuncsData::Num2Cell, FunctionMetaData(-2, 1, DATA));
	(*std_functions)["isrow"]     = BuiltinFunc(BuiltInFuncsData::IsRow, FunctionMetaData(1, 1, DATA));
    (*std_functions)["iscolumn"]  = BuiltinFunc(BuiltInFuncsData::IsColumn, FunctionMetaData(1, 1, DATA));
    (*std_functions)["sortrows"] = BuiltinFunc(BuiltInFuncsData::Sortrows, FunctionMetaData(-2, -2, DATA));

    // File I/O
	(*std_functions)["textread"]  = BuiltinFunc(BuiltInFuncsFile::Textread, FunctionMetaData(1, 1, FILEIO));
    (*std_functions)["dlmwrite"]  = BuiltinFunc(BuiltInFuncsFile::Dlmwrite, FunctionMetaData(-5, 1, FILEIO));
    (*std_functions)["copyfile"]  = BuiltinFunc(BuiltInFuncsFile::Copyfile, FunctionMetaData(3, 3, FILEIO));
    (*std_functions)["importdata"] = BuiltinFunc(BuiltInFuncsFile::Importdata, FunctionMetaData(-1, -1, FILEIO));
	(*std_functions)["textscan"]   = BuiltinFunc(BuiltInFuncsFile::Textscan, FunctionMetaData(-1, 1, FILEIO));
    (*std_functions)["fscanf"]     = BuiltinFunc(BuiltInFuncsFile::Fscanf, FunctionMetaData(-3, 1, FILEIO));
    (*std_functions)["rename"]     = BuiltinFunc(BuiltInFuncsFile::Rename, FunctionMetaData(2, 2, FILEIO));
    (*std_functions)["dlmread"]    = BuiltinFunc(BuiltInFuncsFile::Dlmread, FunctionMetaData(-6, 1, FILEIO));
    (*std_functions)["pwd"]        = BuiltinFunc(oml_pwd, FunctionMetaData(0, 1, FILEIO));
    (*std_functions)["fopen"]      = BuiltinFunc(BuiltInFuncsFile::Fopen, FunctionMetaData(-3, -2, FILEIO));
    (*std_functions)["csvread"]    = BuiltinFunc(BuiltInFuncsFile::Csvread, FunctionMetaData(-6, 1, FILEIO));
    (*std_functions)["csvwrite"]   = BuiltinFunc(BuiltInFuncsFile::Csvwrite, FunctionMetaData(-5, 1, FILEIO));

    // String functions
    (*std_functions)["blanks"]   = BuiltinFunc(BuiltInFuncsString::hml_blanks, 
                                       FunctionMetaData(1, 1, STNG));
    (*std_functions)["strtok"]   = BuiltinFunc(BuiltInFuncsString::hml_strtok, 
                                       FunctionMetaData(2, 2, STNG));
    (*std_functions)["strvcat"]  = BuiltinFunc(BuiltInFuncsString::hml_strvcat,
                                       FunctionMetaData(1, 1, STNG));
    (*std_functions)["sscanf"]   = BuiltinFunc(BuiltInFuncsString::Sscanf,
                                       FunctionMetaData(-3, 1, STNG));
    (*std_functions)["toascii"]     = BuiltinFunc(BuiltInFuncsString::ToAscii, FunctionMetaData(1, 1, STNG));
    (*std_functions)["strtrim"]     = BuiltinFunc(BuiltInFuncsString::Strtrim, FunctionMetaData(1, 1, STNG));
    (*std_functions)["mat2str"]     = BuiltinFunc(BuiltInFuncsString::Mat2Str, FunctionMetaData(-2, 1, STNG));
    (*std_functions)["regexprep"]   = BuiltinFunc(BuiltInFuncsString::Regexprep, FunctionMetaData(1, -4, STNG));
    (*std_functions)["circshift"]   = BuiltinFunc(oml_circshift, FunctionMetaData(-3, 1, ELEM));
    (*std_functions)["shiftdim"]    = BuiltinFunc(oml_shiftdim, FunctionMetaData(-2, -2, ELEM));
    (*std_functions)["checksyntax"] = BuiltinFunc(oml_checksyntax, FunctionMetaData(1, 1, CORE));
    (*std_functions)["ast"]         = BuiltinFunc(oml_ast, FunctionMetaData(1, 1, CORE));
    (*std_functions)["num2str"]     = BuiltinFunc(BuiltInFuncsString::Num2Str, FunctionMetaData(2, 1, STNG));
    (*std_functions)["str2mat"]     = BuiltinFunc(BuiltInFuncsString::Str2mat, FunctionMetaData(1, -1, STNG));
	(*std_functions)["cellstr"] = BuiltinFunc(BuiltInFuncsString::CellStr, FunctionMetaData(1, -1, STNG));
    (*std_functions)["str2double"]  = BuiltinFunc(BuiltInFuncsString::Str2Double, FunctionMetaData(1, 1, STNG));
    (*std_functions)["contains"]    = BuiltinFunc(BuiltInFuncsString::Contains, FunctionMetaData(-2, 1, STNG));
    (*std_functions)["strip"] = BuiltinFunc(BuiltInFuncsString::Strip, FunctionMetaData(2, 1, STNG));

    // System functions
    (*std_functions)["dir"]      = BuiltinFunc(BuiltInFuncsSystem::Dir, FunctionMetaData(1, 1, SYSTEM));
    (*std_functions)["ls"]       = BuiltinFunc(BuiltInFuncsSystem::Ls, FunctionMetaData(-2, -1, SYSTEM));
    (*std_functions)["system"]   = BuiltinFunc(BuiltInFuncsSystem::System, FunctionMetaData(-2, -2, SYSTEM));
    (*std_functions)["dos"]      = BuiltinFunc(BuiltInFuncsSystem::System, FunctionMetaData(-2, -2, SYSTEM));
    (*std_functions)["unix"]     = BuiltinFunc(BuiltInFuncsSystem::Unix,   FunctionMetaData(2, 2, SYSTEM));
    (*std_functions)["delete"]   = BuiltinFunc(BuiltInFuncsSystem::Delete, FunctionMetaData(-1, 0, SYSTEM));
    (*std_functions)["chdir"]   = BuiltinFunc(BuiltInFuncsSystem::Cd,     FunctionMetaData(-1, 1, SYSTEM));
    (*std_functions)["cd"]      = BuiltinFunc(BuiltInFuncsSystem::Cd,     FunctionMetaData(-1, 1, SYSTEM));
    (*std_functions)["rmdir"]   = BuiltinFunc(BuiltInFuncsSystem::Rmdir,  FunctionMetaData(-2, 3, SYSTEM));
    (*std_functions)["genpath"] = BuiltinFunc(BuiltInFuncsSystem::Genpath, FunctionMetaData(-2, 1, SYSTEM));
    (*std_functions)["getpid"]  = BuiltinFunc(BuiltInFuncsSystem::GetPid, FunctionMetaData(0, 1, SYSTEM));
    (*std_functions)["make_absolute_filename"] = BuiltinFunc(BuiltInFuncsSystem::MakeAbsFilename, FunctionMetaData(1, 1, SYSTEM));
    (*std_functions)["filesep"] = BuiltinFunc(BuiltInFuncsSystem::Filesep, FunctionMetaData(-1, 1, SYSTEM));
    (*std_functions)["diary"]   = BuiltinFunc(BuiltInFuncsSystem::Diary, FunctionMetaData(1, 0, SYSTEM));
    (*std_functions)["outputlog"] = BuiltinFunc(BuiltInFuncsSystem::OutputLog, FunctionMetaData(1, 0, SYSTEM));

    // Time functions
    (*std_functions)["time"]    = BuiltinFunc(BuiltInFuncsSystem::Time, FunctionMetaData(0, 1, TIME));
    (*std_functions)["ctime"]   = BuiltinFunc(BuiltInFuncsSystem::CTime, FunctionMetaData(1, 1, TIME));
    (*std_functions)["year"]    = BuiltinFunc(BuiltInFuncsTime::Year, FunctionMetaData(-2, 1, TIME));
    (*std_functions)["day"]     = BuiltinFunc(BuiltInFuncsTime::Day, FunctionMetaData(-2, 1, TIME));
    (*std_functions)["month"]   = BuiltinFunc(BuiltInFuncsTime::Month, FunctionMetaData(-2, -2, TIME));
    (*std_functions)["minute"]  = BuiltinFunc(BuiltInFuncsTime::Minute, FunctionMetaData(-2, 1, TIME));
    (*std_functions)["hour"]    = BuiltinFunc(BuiltInFuncsTime::Hour, FunctionMetaData(-2, 1, TIME));
    (*std_functions)["second"]  = BuiltinFunc(BuiltInFuncsTime::Second, FunctionMetaData(-2, 1, TIME));
    (*std_functions)["now"]     = BuiltinFunc(oml_now, FunctionMetaData(0, 1, TIME));
    (*std_functions)["quarter"] = BuiltinFunc(BuiltInFuncsTime::Quarter, FunctionMetaData(-2, 1, TIME));

    // Client specific environment related functions
    (*std_functions)["getbaseenv"]    = BuiltinFunc(BuiltInFuncsSystem::GetBaseEnv, 
                                        FunctionMetaData(0, 1));
    (*std_functions)["getcurrentenv"] = BuiltinFunc(BuiltInFuncsSystem::GetCurrentEnv, 
                                        FunctionMetaData(0, 1));
    (*std_functions)["getnewenv"]     = BuiltinFunc(BuiltInFuncsSystem::GetNewEnv, 
                                        FunctionMetaData(0, 1));
    (*std_functions)["getenvvalue"]   = BuiltinFunc(BuiltInFuncsSystem::GetEnvVal, 
                                        FunctionMetaData(2, 1));
    (*std_functions)["setenvvalue"]   = BuiltinFunc(BuiltInFuncsSystem::SetEnvVal, 
                                        FunctionMetaData(3, 0));
    (*std_functions)["clearenvvalue"] = BuiltinFunc(BuiltInFuncsSystem::ClearEnvVal, 
                                        FunctionMetaData(2, 0));
    (*std_functions)["importenv"]     = BuiltinFunc(BuiltInFuncsSystem::ImportEnv, 
                                        FunctionMetaData(1, 0));
    (*std_functions)["importenvin"]   = BuiltinFunc(BuiltInFuncsSystem::ImportEnvIn, 
                                        FunctionMetaData(2, 0));
    (*std_functions)["cloneenv"]      = BuiltinFunc(BuiltInFuncsSystem::CloneEnv, 
                                        FunctionMetaData(1, 1));
    // End of client specific environment related functions
	 
	(*std_functions)["memoryuse"]      = BuiltinFunc(oml_memoryuse, FunctionMetaData(0, 1, CORE));
	(*std_functions)["analyze"]        = BuiltinFunc(oml_analyze, FunctionMetaData(1, 0, CORE));
	(*std_functions)["getmetadata"]    = BuiltinFunc(oml_getmetadata, FunctionMetaData(1, 0, CORE));
	(*std_functions)["rehash"]         = BuiltinFunc(oml_rehash, FunctionMetaData(0, 0, CORE));
	(*std_functions)["verbose"]        = BuiltinFunc(oml_verbose, FunctionMetaData(0, 1, CORE));
	(*std_functions)["parcluster"]     = BuiltinFunc(oml_parcluster, FunctionMetaData(0, 1, CORE));

    // Function place holders have been added for the following. 
    // The implementation is on the client side
    (*std_functions)["getargc"] = BuiltinFunc(oml_getargc, FunctionMetaData(0, 1, CORE));
    (*std_functions)["getargv"] = BuiltinFunc(oml_getargv, FunctionMetaData(1, 1, CORE));

    (*std_functions)["inputname"] = BuiltinFunc(oml_inputname, FunctionMetaData(1, 1, CORE));

    // MKL control parameter settings (hidden)
    BuiltinFunc mklsdpt = BuiltinFunc(BuiltInFuncsMKL::oml_mkl_sdpt, FunctionMetaData(1, 0, CORE));
    mklsdpt.Lock();
    (*std_functions)["mklsdpt"] = mklsdpt;

	bool experimental = true;

#if !defined(_DEBUG)
	experimental = Currency::GetExperimental();
#endif

	(*std_functions)["getprofiledata"] = BuiltinFunc(oml_getprofiledata, FunctionMetaData(0, 1, CORE));
	(*std_functions)["clearprofiledata"] = BuiltinFunc(oml_clearprofiledata, FunctionMetaData(0, 1, CORE));
	(*std_functions)["profile"] = BuiltinFunc(oml_profile, FunctionMetaData(0, 1, CORE));

	if (experimental)
	{
		(*std_functions)["writepfile"] = BuiltinFunc(oml_writepfile, FunctionMetaData(2, 0, CORE));

		// helptest function is restricted to experimental mode
	    (*std_functions)["helptest"]      = BuiltinFunc(oml_helptest, FunctionMetaData(1,1, CORE));
	    (*std_functions)["helpmodule"]      = BuiltinFunc(oml_helpmodule, FunctionMetaData(1,1, CORE));

		(*std_functions)["p_definefunction"] = BuiltinFunc(oml_p_definefunction, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_registerfunction"] = BuiltinFunc(oml_p_registerfunction, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_assign"] = BuiltinFunc(oml_p_assign, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_multassign"] = BuiltinFunc(oml_p_multassign, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_struct"] = BuiltinFunc(oml_p_struct, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_func"] = BuiltinFunc(oml_p_func, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_add"] = BuiltinFunc(oml_p_add, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_statement"] = BuiltinFunc(oml_p_statement, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_handle"] = BuiltinFunc(oml_p_handle, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_string"] = BuiltinFunc(oml_p_string, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_equal"] = BuiltinFunc(oml_p_equal, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_if"] = BuiltinFunc(oml_p_if, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_elseif"] = BuiltinFunc(oml_p_elseif, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_else"] = BuiltinFunc(oml_p_else, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_conditional"] = BuiltinFunc(oml_p_conditional, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_divide"] = BuiltinFunc(oml_p_divide, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_bitand"] = BuiltinFunc(oml_p_bitand, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_entrypow"] = BuiltinFunc(oml_p_entrypow, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_entrydivide"] = BuiltinFunc(oml_p_entrydivide, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_entryleftdivide"] = BuiltinFunc(oml_p_entryleftdivide, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_greaterequal"] = BuiltinFunc(oml_p_greaterequal, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_greaterthan"] = BuiltinFunc(oml_p_greaterthan, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_logand"] = BuiltinFunc(oml_p_logand, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_leftdivide"] = BuiltinFunc(oml_p_leftdivide, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_lessthan"] = BuiltinFunc(oml_p_lessthan, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_lessequal"] = BuiltinFunc(oml_p_lessequal, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_logor"] = BuiltinFunc(oml_p_logor, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_subtract"] = BuiltinFunc(oml_p_subtract, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_notequal"] = BuiltinFunc(oml_p_notequal, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_or"] = BuiltinFunc(oml_p_or, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_pow"] = BuiltinFunc(oml_p_pow, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_multiply"] = BuiltinFunc(oml_p_multiply, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_negate"] = BuiltinFunc(oml_p_negate, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_uminus"] = BuiltinFunc(oml_p_uminus, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_transpose"] = BuiltinFunc(oml_p_transpose, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_ctranspose"] = BuiltinFunc(oml_p_ctranspose, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_global"] = BuiltinFunc(oml_p_global, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_persistent"] = BuiltinFunc(oml_p_persistent, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_range"] = BuiltinFunc(oml_p_range, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_forloop"] = BuiltinFunc(oml_p_forloop, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_trycatch"] = BuiltinFunc(oml_p_trycatch, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_try"] = BuiltinFunc(oml_p_try, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_catch"] = BuiltinFunc(oml_p_catch, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_cellvalue"] = BuiltinFunc(oml_p_cellvalue, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_whileloop"] = BuiltinFunc(oml_p_whileloop, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_switch"] = BuiltinFunc(oml_p_switch, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_case"] = BuiltinFunc(oml_p_case, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_otherwise"] = BuiltinFunc(oml_p_otherwise, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_matrix"] = BuiltinFunc(oml_p_matrix, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_vector"] = BuiltinFunc(oml_p_vector, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_cellarray"] = BuiltinFunc(oml_p_cellarray, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_colon"] = BuiltinFunc(oml_p_colon, FunctionMetaData(0, 0, CORE));
		(*std_functions)["p_return"] = BuiltinFunc(oml_p_return, FunctionMetaData(0, 0, CORE));
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
    BuiltInFuncsUtils::CheckFileIndex(eval, fileID, 1, true);
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
            throw OML_Error(OML_ERR_MULTILINE_STRING, i + 1);

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
            throw OML_Error(OML_ERR_MULTILINE_STRING, i + 1);

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
bool oml_fileparts(EvaluatorInterface eval, 
                   const std::vector<Currency>& inputs, 
                   std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    if (inputs[0].IsMultilineString())
    {
        throw OML_Error(OML_ERR_MULTILINE_STRING, 1);
    }

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
        if (lastslash == 0)
        {
            outputs.push_back(filename.substr(0, 1));
        }
        else
        {
            outputs.push_back(filename.substr(0, lastslash));
        }
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

	int min = (int)inputs[0].Scalar();
	int max;

	if (inputs[1].IsScalar())
	{
		if (!inputs[1].IsInteger())
		{
			double val = inputs[1].Scalar();

			if (IsInf_T(val))
				max = INT_MAX;
			else
				throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_VALUE);
		}
		else
		{
			max = (int)inputs[1].Scalar();
		}

		if (max == 0.0)
			throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_VALUE);
	}

	if (IsInf_T(inputs[1].Scalar()))
		max = INT_MAX;

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
// Returns true and joins path names to build complete filename [fullfile]
//------------------------------------------------------------------------------
bool oml_fullfile(EvaluatorInterface           eval, 
                  const std::vector<Currency>& inputs, 
                  std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    int nargin = static_cast<int>(inputs.size());

    // Check the last element as this determines output type
    if (!inputs[nargin - 1].IsString() && !inputs[nargin - 1].IsCellArray())
    {
        throw OML_Error(OML_ERR_STRING_STRINGCELL, nargin);
    }

    HML_CELLARRAY* cell = nullptr;
    std::unique_ptr<HML_CELLARRAY> out(nullptr);
    if (inputs[nargin - 1].IsCellArray())
    {
        cell = inputs[nargin - 1].CellArray();
        if (!cell || cell->IsEmpty())
        {
            out.reset(EvaluatorInterface::allocateCellArray());
            outputs.push_back(out.release());
            return true;
        }
        out.reset(EvaluatorInterface::allocateCellArray(1, cell->Size()));
    }

    std::string       path;
    BuiltInFuncsUtils utils;
    for (int i = 0; i < nargin; ++i)
    {
        const Currency& cur = inputs[i];

        // Everything preceding last input must be a string. Checks have already
        // been made for last input
        if (i < nargin - 1)
        {
            if (!cur.IsString())
            {
                throw OML_Error(OML_ERR_STRING, i + 1, OML_VAR_TYPE);
            }
        }

        std::string val;
        if (cur.IsString())
        {
            if (!cur.IsMultilineString())
            {
                val = cur.StringVal();
            }
            else
            {
                Currency tmp = utils.ReadRow(eval, cur, 0);
                if (!tmp.IsString())
                {
                    throw OML_Error(OML_ERR_STRING, i + 1, OML_VAR_TYPE);
                }
                val = tmp.StringVal();
            }
            val = utils.RTrim(val);
            val = utils.LTrim(val);

            if (val.empty())
            {
                continue;
            }
            if (!path.empty())
            {
                path += "/";
            }
            path += val;
        }
        else if (cell)
        {
            int cellsize = cell->Size();
            for (int j = 0; j < cellsize; ++j)
            {
                const Currency& element = (*cell)(j);
                if (!element.IsString())
                {
                    throw OML_Error(OML_ERR_STRING_STRINGCELL, nargin);
                }

                std::string val;
                if (!element.IsMultilineString())
                {
                    val = element.StringVal();
                }
                else
                {
                    Currency tmp = utils.ReadRow(eval, element, 0);
                    if (!tmp.IsString())
                    {
                        throw OML_Error(OML_ERR_STRING, i + 1, OML_VAR_TYPE);
                    }
                    val = tmp.StringVal();
                }
                val = utils.RTrim(val);
                val = utils.LTrim(val);
                if (val.empty())
                {
                    continue;
                }
                std::string newpath = path;
                if (!path.empty())
                {
                    newpath += "/";
                }
                newpath += val;
                newpath = utils.StripMultipleSlashesAndNormalize(newpath);
                (*out)(j) = newpath;
            }
        }
    }

    if (out)
    {
        outputs.push_back(out.release());
    }
    else
    {
        outputs.push_back(utils.StripMultipleSlashesAndNormalize(path));
    }

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
    else if (in1.IsSparse())
    {
        outputs.push_back(HW_MAKE_BOOL_CURRENCY(in1.MatrixS()->IsVector()));
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
        if (!eval.FindFunctionByName(fname, &fi, &fptr, NULL))
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
        if (!eval.FindFunctionByName(fname, &fi, &fptr, NULL))
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
bool oml_isdir(EvaluatorInterface eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsString())
	{
		throw OML_Error(OML_ERR_STRING, 1);
	}
    std::string path (readString(inputs[0]));
    outputs.push_back(BuiltInFuncsUtils::IsDir(path));
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
    return allowComplexToLogical(inputs, outputs, &checkisfinite, true);
}
//------------------------------------------------------------------------------
// Unary plus sign [uplus]
//------------------------------------------------------------------------------
bool oml_uplus(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

	Currency temp = inputs[0];
	temp.SetMask(Currency::MASK_DOUBLE);

    outputs.push_back(temp);

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
    size_t nargin = inputs.size();

    if (nargin != 1 && nargin != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_VARIABLE);

    double tol = -999.0;

    if (nargin == 2)
    {
        if (!inputs[1].IsScalar())
            throw OML_Error(OML_ERR_NONNEGATIVE_SCALAR, 2, OML_VAR_VARIABLE);

        tol = inputs[1].Scalar();

        if (tol < 0.0)
            throw OML_Error(OML_ERR_NONNEGATIVE_SCALAR, 2, OML_VAR_VARIABLE);
    }

    const hwMatrix* mtx = inputs[0].ConvertToMatrix();
    hwMatrix* result = EvaluatorInterface::allocateMatrix();
    BuiltInFuncsUtils::CheckMathStatus(eval, (result->Pinv)(*mtx, tol));
    outputs.push_back(result);
    return true;
}
//------------------------------------------------------------------------------
// Normalizes vectors [normalize]
//------------------------------------------------------------------------------
bool oml_normalize(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsString())
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    const hwMatrix* mtx = inputs[0].ConvertToMatrix();

    if (!mtx)
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    if (mtx->IsVector())
    {
        if (mtx->IsReal())
        {
            std::unique_ptr<hwMatrix> mtxcopy(new hwMatrix(*mtx));
            BuiltInFuncsUtils::CheckMathStatus(eval, mtxcopy->Normalize());
            outputs.push_back(mtxcopy.release());
        }
        else if (mtx->IsRealData())
        {
            std::unique_ptr<hwMatrix> mtxcopy(makeRealCopy(mtx));
            BuiltInFuncsUtils::CheckMathStatus(eval, mtxcopy->Normalize());
            outputs.push_back(mtxcopy.release());
        }
        else
            throw OML_Error(OML_ERR_REAL);
    }
    else
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_DATA);

    return true;
}
//------------------------------------------------------------------------------
// Returns true after making directories [mkdir]
//------------------------------------------------------------------------------
bool oml_mkdir(EvaluatorInterface           eval, 
               const std::vector<Currency>& inputs, 
               std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }

    size_t nargin = inputs.size();
    if (nargin > 1 && !inputs[1].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 2);
    }

    BuiltInFuncsUtils utils;

#ifdef OS_WIN
    std::wstring parentdir;
    std::wstring newdir;
    if (nargin == 1)
    {
        newdir = utils.StdString2WString(inputs[0].StringVal());
    }
    else
    {
        parentdir = utils.StdString2WString(inputs[0].StringVal());
        newdir = utils.StdString2WString(inputs[1].StringVal());
    }

    std::wstring fullpath(newdir);

	if (!parentdir.empty())
	{
		parentdir = utils.StripTrailingSlashW(parentdir);
        fullpath = parentdir + L"/" + newdir;
	}
    else
    {
        bool isabsolutepath = utils.IsAbsolutePathW(newdir);
        if (!isabsolutepath)
        {
            fullpath = utils.GetCurrentWorkingDirW();
            utils.StripTrailingSlashW(fullpath);
            fullpath += L"/" + newdir;
        }
    }
    fullpath = utils.GetNormpathW(fullpath);
    if (fullpath.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, static_cast<int>(nargin));
    }

    // Check to see all parent directories exist
    size_t pos = fullpath.find_last_of(L'\\');
    if (pos != std::wstring::npos)
    {
        parentdir = fullpath.substr(0, pos);
        if (!utils.DoesPathExistW(parentdir))
        {
            utils.StripTrailingSlashW(parentdir);
            pos = parentdir.find_first_not_of(L":\\", 1);
            if (pos != std::wstring::npos)
            {
                std::wstring dir = parentdir.substr(0, pos);
                std::wstring path = parentdir.substr(pos);
                wchar_t* tok = wcstok((wchar_t*)path.c_str(), L"\\");
                while (tok)
                {
                    utils.StripTrailingSlashW(dir);
                    dir += L"\\" + std::wstring(tok);

                    if (!utils.DoesPathExistW(dir))
                    {
                        BOOL result = CreateDirectoryW(dir.c_str(), nullptr);
                        if (!result)
                        {
                            break;
                        }
                    }
                    tok = wcstok(nullptr, L"\\");
                }
            }
        }
    }
    BOOL result = CreateDirectoryW(fullpath.c_str(), nullptr);
    if (!result)
    {
        DWORD err = GetLastError();
        outputs.push_back(Currency(false));
        if (err == ERROR_ALREADY_EXISTS)
            outputs.push_back("Directory already exists");
        else if (err == ERROR_PATH_NOT_FOUND)
            outputs.push_back("No such file or directory");
        else
            outputs.push_back("Could not create directory");
        outputs.push_back("mkdir");
    }
#else
    std::string parentdir;
    std::string newdir;
    if (nargin == 1)
    {
        newdir = inputs[0].StringVal();
    }
    else
    {
        parentdir = inputs[0].StringVal();
        newdir = inputs[1].StringVal();
    }

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
    std::string strcmd ("mkdir -p \"" + fullpath + "\" > /dev/null");
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
    {
        throw OML_Error(HW_ERROR_ENVVARFAIL);
    }
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
                                            in1.IsMatrix() || in1.IsNDMatrix() || in1.IsSparse())));
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
		if (inputs[0].IsNDMatrix())
		{
			if (nargin == 1)
			{
				oml_MatrixNUtil3(eval, inputs, outputs, oml_all);
			}
			else
			{
				oml_MatrixNUtil3(eval, inputs, outputs, oml_all, 2);
			}
		}
		else if (inputs[0].IsSparse())
		{
			const hwMatrixS* spm = inputs[0].MatrixS();
			int dim;

			if (nargin > 1)
			{
				if (!inputs[1].IsPositiveInteger())
					throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

				dim = static_cast<int> (inputs[1].Scalar());
			}
			else
			{
				dim = spm->M() == 1 ? 2 : 1;
			}

			if (dim == 1)
			{
				MKL_INT* pb = const_cast<MKL_INT*> (spm->pointerB());
				MKL_INT* pe = const_cast<MKL_INT*> (spm->pointerE());
				std::vector<int> ivec;
				std::vector<int> jvec;
				hwMatrix V;
				hwMatrixS* ret = new hwMatrixS(ivec, jvec, V, 1, spm->N());

				for (int i = 0; i < spm->N(); ++i)
				{
					if (pe[i] - pb[i] == spm->M())
						(*ret)(0, i) = 1.0;
				}

				outputs.push_back(ret);
			}
			else if (dim == 2)
			{
				MKL_INT* pr = const_cast<MKL_INT*> (spm->rows());
				std::vector<int> ivec;
				std::vector<int> jvec;
				hwMatrix V;
				hwMatrixS* ret = new hwMatrixS(ivec, jvec, V, spm->M(), 1);
				std::vector<int> count(spm->M());

				for (int i = 0; i < spm->NNZ(); ++i)
				{
					count[pr[i]]++;
				}

				for (int i = 0; i < spm->M(); ++i)
				{
					if (count[i] == spm->N())
						(*ret)(i, 0) = 1.0;
				}

				outputs.push_back(ret);
			}
			else
			{
				std::vector<Currency> inputs2;
				inputs2.push_back(inputs[0]);
				oml_spones(eval, inputs2, outputs);
			}
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
		if (inputs[0].IsNDMatrix())
		{
			if (nargin == 1)
			{
				oml_MatrixNUtil3(eval, inputs, outputs, oml_any);
			}
			else
			{
				oml_MatrixNUtil3(eval, inputs, outputs, oml_any, 2);
			}
		}
		else if (inputs[0].IsSparse())
		{
			const hwMatrixS* spm = inputs[0].MatrixS();
			int dim;

			if (nargin > 1)
			{
				if (!inputs[1].IsPositiveInteger())
					throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VALUE);

				dim = static_cast<int> (inputs[1].Scalar());
			}
			else
			{
				dim = spm->M() == 1 ? 2 : 1;
			}

			if (dim == 1)
			{
				MKL_INT* pb = const_cast<MKL_INT*> (spm->pointerB());
				MKL_INT* pe = const_cast<MKL_INT*> (spm->pointerE());
				std::vector<int> ivec;
				std::vector<int> jvec;
				hwMatrix V;
				hwMatrixS* ret = new hwMatrixS(ivec, jvec, V, 1, spm->N());

				for (int i = 0; i < spm->N(); ++i)
				{
					if (pe[i] != pb[i])
						(*ret)(0, i) = 1.0;
				}

				outputs.push_back(ret);
			}
			else if (dim == 2)
			{
				MKL_INT* pr = const_cast<MKL_INT*> (spm->rows());
				std::vector<int> ivec;
				std::vector<int> jvec;
				hwMatrix V;
				hwMatrixS* ret = new hwMatrixS(ivec, jvec, V, spm->M(), 1);

				for (int i = 0; i < spm->NNZ(); ++i)
				{
					(*ret)(pr[i], 0) = 1.0;

					if (ret->NNZ() == ret->M())
						break;
				}

				outputs.push_back(ret);
			}
			else
			{
				std::vector<Currency> inputs2;
				inputs2.push_back(inputs[0]);
				oml_spones(eval, inputs2, outputs);
			}
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
    else if (cur.IsSparse())
        return cur.MatrixS()->Size() == 1;
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

	if (nargin == 1)
	{
		if (!inputs[0].IsString())
			throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

		std::string msg = readString(inputs[0]);

		if (msg.length())
			BuiltInFuncsUtils::SetWarning(eval, "Warning: " + msg);
	}
	else if (nargin == 2)
	{
		if (!inputs[0].IsString())
			throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);

		std::string instruction = inputs[0].StringVal();
		std::string value;
		
		if (inputs[1].IsString())
			value = inputs[1].StringVal();
		
		if (instruction == "on")
		{
			eval.EnableWarning(value);
		}
		else if (instruction == "off")
		{
			eval.DisableWarning(value);
		}
		else if (instruction == "query")
		{
			bool output = !eval.IsWarningDisabled(value);
			outputs.push_back(output);
		}
		else
		{
			std::string msg = sprintf(eval, inputs.cbegin(), inputs.cend());

			if (msg.length())
				BuiltInFuncsUtils::SetWarning(eval, "Warning: " + msg);
		}
	}
	else
	{
		std::string msg = sprintf(eval, inputs.cbegin(), inputs.cend());

		if (msg.length())
			BuiltInFuncsUtils::SetWarning(eval, "Warning: " + msg);
	}

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
    outputs.push_back(HW_MAKE_BOOL_CURRENCY(cur.IsScalar() || cur.IsComplex() || cur.IsMatrix() ||
                      cur.IsString() || cur.IsNDMatrix() || cur.IsSparse()));
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
    }
    else if (input1.IsMatrixOrString())
    {
        const hwMatrix* matrix = input1.Matrix();

        if (matrix->IsSquare())
            outputs.push_back(true);
        else
            outputs.push_back(false);
    }
    else if (input1.IsSparse())
    {
        const hwMatrixS* matrix = input1.MatrixS();

        if (matrix->IsSquare())
            outputs.push_back(true);
        else
            outputs.push_back(false);
    }
    else
    {
        outputs.push_back(false);
    }

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
// Returns true if input is a diagonal matrix [isdiag]
//------------------------------------------------------------------------------
bool oml_isdiag(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];

    if (input1.IsScalar() || input1.IsComplex())
    {
        outputs.push_back(true);
    }
    else if (input1.IsMatrix())
    {
        outputs.push_back(input1.Matrix()->IsDiag());
    }
    else if (input1.IsSparse())
    {
        const hwMatrixS* sparse = input1.MatrixS();
        int nnz = sparse->NNZ();
        int m = sparse->M();
        std::vector<int> index;
        sparse->NZinfo(0, nnz - 1, index);

        // check for non-zeros with row != col
        for (int k = 0; k < nnz; ++k)
        {
            if (index[k] % m != index[k] / m)
            {
                outputs.push_back(false);
                return true;
            }
        }

        outputs.push_back(true);
    }
    else
    {
        outputs.push_back(false);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is a lower traingular matrix [istril]
//------------------------------------------------------------------------------
bool oml_istril(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];

    if (input1.IsScalar() || input1.IsComplex())
    {
        outputs.push_back(true);
    }
    else if (input1.IsMatrix())
    {
        outputs.push_back(input1.Matrix()->IsLowerTriag());
    }
    else if (input1.IsSparse())
    {
        const hwMatrixS* sparse = input1.MatrixS();
        int nnz = sparse->NNZ();
        int m = sparse->M();
        std::vector<int> index;
        sparse->NZinfo(0, nnz - 1, index);

        // check for non-zeros with row < col
        for (int k = 0; k < nnz; ++k)
        {
            if (index[k] % m < index[k] / m)
            {
                outputs.push_back(false);
                return true;
            }
        }
    
        outputs.push_back(true);
    }
    else
    {
        outputs.push_back(false);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is an upper traingular matrix [istriu]
//------------------------------------------------------------------------------
bool oml_istriu(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input1 = inputs[0];

    if (input1.IsScalar() || input1.IsComplex())
    {
        outputs.push_back(true);
    }
    else if (input1.IsMatrix())
    {
        outputs.push_back(input1.Matrix()->IsUpperTriag());
    }
    else if (input1.IsSparse())
    {
        const hwMatrixS* sparse = input1.MatrixS();
        int nnz = sparse->NNZ();
        int m = sparse->M();
        std::vector<int> index;
        sparse->NZinfo(0, nnz - 1, index);

        // check for non-zeros with row > col
        for (int k = 0; k < nnz; ++k)
        {
            if (index[k] % m > index[k] / m)
            {
                outputs.push_back(false);
                return true;
            }
        }

        outputs.push_back(true);
    }
    else
    {
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
    if (inputs.empty() || inputs.size() < 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    // Get output options, in user specified order if applicable
	std::vector<std::string> flags;
    std::vector<OMLREGEXP> options (GetRegexOptions(
        inputs, eval.GetNargoutValue(), flags));

    GetRegExpOutput(eval, inputs[0], inputs[1], options, flags, outputs);
    return true;
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
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);

    if (!inputs[1].IsMatrix() && !inputs[1].IsScalar() && !inputs[1].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);

    const hwMatrix* a = inputs[0].ConvertToMatrix();
    const hwMatrix* b = inputs[1].ConvertToMatrix();
    const hwMatrix* mtx = nullptr;
    hwMatrix *result = EvaluatorInterface::allocateMatrix();

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
            mtx = in3.ConvertToMatrix();
        }
    }
    else if (nargin == 4)
    {
        const Currency &in3 = inputs[2];

        if (in3.IsString())
            throw OML_Error(OML_ERR_MATRIX, 3, OML_VAR_DATA);

        mtx = in3.ConvertToMatrix();
        shape = readOption(eval, inputs[3]);
    }

    if (mtx && !(a->IsVector() && b->IsVector()))
        throw OML_Error(HW_ERROR_IFCONVMATINPMUSTVEC);

    if (shape == "full")
        shapeCode = Full;
    else if (shape == "same")
        shapeCode = Same;
    else if (shape == "valid")
        shapeCode = Valid;
    else
        throw OML_Error(HW_ERROR_INVSHAPEFULLSAMEVALID);

    if (!(b->Size() && a->Size()))
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        return true;
    }

    int top, left, m, n;

    if (mtx)
    {
        int am = mtx->M();
        int an = mtx->N();
        int bm = a->Size();
        int bn = b->Size();

        if (bm * bn == 0)
        {
            outputs.push_back(EvaluatorInterface::allocateMatrix());
            return true;
        }

        if (shapeCode == Same)
        {
            top = bm / 2;
            left = bn / 2;
            m = am;
            n = an;
            BuiltInFuncsUtils::CheckMathStatus(eval, result->Conv2D(*a, *b, *mtx, top, left, top + m - 1, left + n - 1));
        }
        else if (shapeCode == Valid)
        {
            top = bm - 1;
            left = bn - 1;
            m = max(am - bm + 1, 0);
            n = max(an - bn + 1, 0);
            BuiltInFuncsUtils::CheckMathStatus(eval, result->Conv2D(*a, *b, *mtx, top, left, top + m - 1, left + n - 1));
        }
        else
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, result->Conv2D(*a, *b, *mtx));
        }
    }
    else
    {
        int am = a->M();
        int an = a->N();
        int bm = b->M();
        int bn = b->N();

        if (shapeCode == Same)
        {
            top = bm / 2;
            left = bn / 2;
            m = am;
            n = an;
            BuiltInFuncsUtils::CheckMathStatus(eval, result->Conv2D(*a, *b, top, left, top + m - 1, left + n - 1));
        }
        else if (shapeCode == Valid)
        {
            top = bm - 1;
            left = bn - 1;
            m = max(am - bm + 1, 0);
            n = max(an - bn + 1, 0);
            BuiltInFuncsUtils::CheckMathStatus(eval, result->Conv2D(*a, *b, top, left, top + m - 1, left + n - 1));
        }
        else    // shapeCode == Full
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, result->Conv2D(*a, *b));
        }
    }

    outputs.push_back(result);
    return true;
}

//------------------------------------------------------------------------------
// hwTMatrix partial template specialization for Currency, with T2 = void*
// Set every matrix element to a specified real value
//------------------------------------------------------------------------------
template<>
inline void hwTMatrix<Currency, void*>::SetElements(Currency value)
{
    int size = Size();

    if (m_real)
    {
        for (int i = 0; i < size; ++i)
            m_real[i] = value;
    }
}

//------------------------------------------------------------------------------
// hwTMatrix partial template specialization for Currency, with T2 = void*
// Convert the matrix to complex
//------------------------------------------------------------------------------
template<>
inline hwMathStatus hwTMatrix<Currency, void*>::MakeComplex()
{
    return hwMathStatus();
}

//------------------------------------------------------------------------------
// Helper method for matrix concatenation - assumes both inputs of the same type
// real or complex
//------------------------------------------------------------------------------
template<bool HORIZ, typename T1, typename T2>
static hwTMatrix<T1, T2>* concat(const hwTMatrix<T1, T2> *m1, const hwTMatrix<T1, T2> *m2)
{
    if (m1->Is0x0())
    {
        return new hwTMatrix<T1, T2>(*m2);
    }
    else if (m2->Is0x0())
    {
        return new hwTMatrix<T1, T2>(*m1);
    }

    if ((m1->IsReal() == m2->IsReal()) || m1->IsEmpty() || m2->IsEmpty())
    {
        hwTMatrix<T1, T2>* newmtx = nullptr;

        if (HORIZ)
        {
            if (m1->M() != m2->M())
            {
                throw OML_Error(HW_ERROR_INPUTHCATDIM);
            }

            newmtx = new hwTMatrix<T1, T2>(m1->M(), m1->N() + m2->N(),
                m1->Type());

            hwMathStatus status = newmtx->WriteSubmatrix(0, 0, *m1);

            if (!status.IsOk())
                throw OML_Error(status);

            status = newmtx->WriteSubmatrix(0, m1->N(), *m2);

            if (!status.IsOk())
                throw OML_Error(status);
        }
        else
        {
            if (m1->N() != m2->N())
            {
                throw OML_Error(HW_ERROR_INPUTVCATDIM);
            }

            newmtx = new hwTMatrix<T1, T2>(m1->M() + m2->M(), m1->N(),
                m1->Type());

            hwMathStatus status = newmtx->WriteSubmatrix(0, 0, *m1);

            if (!status.IsOk())
                throw OML_Error(status);

            status = newmtx->WriteSubmatrix(m1->M(), 0, *m2);

            if (!status.IsOk())
                throw OML_Error(status);
        }

        return newmtx;
    }
    else if (m1->IsReal())
    {
        hwTMatrix<T1, T2> m1C;
        m1C.PackComplex(*m1);
        return concat<HORIZ, T1, T2>(&m1C, m2);
    }
    else // m2->IsReal()
    {
        hwTMatrix<T1, T2> m2C;
        m2C.PackComplex(*m2);
        return concat<HORIZ, T1, T2>(m1, &m2C);
    }
}
//------------------------------------------------------------------------------
// Helper method for cell array concatenation
//------------------------------------------------------------------------------
template <bool HORIZ>
static HML_CELLARRAY* concat(EvaluatorInterface& eval, const HML_CELLARRAY *cell, Currency right)
{
    if (right.IsCellArray())
        return concat<HORIZ>(cell, right.CellArray());

    if ((HORIZ && cell->M() != 1) || (!HORIZ && cell->N() != 1))
    {
        if (!cell->Size())
        {
            HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(1, 1);
            (*newcell)(0) = right;
            (*newcell)(0).ClearOutputName();
            return newcell;
        }
        else
            throw OML_Error(HORIZ ? HW_ERROR_INPUTHCATDIM : HW_ERROR_INPUTVCATDIM);
    }

    HML_CELLARRAY* newcell = nullptr;

    if (HORIZ)
    {
        newcell = EvaluatorInterface::allocateCellArray(cell->M(), cell->N() + 1);
    }
    else
    {
        newcell = EvaluatorInterface::allocateCellArray(cell->M() + 1, cell->N());
    }

    int i;
    for (i = 0; i < cell->Size(); i++)
    {
        (*newcell)(i) = (*cell)(i);
    }

    (*newcell)(i) = right;
    (*newcell)(i).ClearOutputName();

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

    if (cell->Size())
    {
        if (HORIZ && cell->M() != 1)
            throw OML_Error(HW_ERROR_INPUTHCATDIM);
        else if (!HORIZ && cell->N() != 1)
            throw OML_Error(HW_ERROR_INPUTVCATDIM);

        HML_CELLARRAY* newcell = nullptr;

        if (HORIZ)
            newcell = EvaluatorInterface::allocateCellArray(cell->M(), cell->N() + 1);
        else
            newcell = EvaluatorInterface::allocateCellArray(cell->M() + 1, cell->N());

        (*newcell)(0) = left;
        (*newcell)(0).ClearOutputName();

        for (int i = 0; i < cell->Size(); i++)
        {
            (*newcell)(i + 1) = (*cell)(i);
            (*newcell)(i + 1).ClearOutputName();
        }

        return newcell;
    }
    else
    {
        HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(1, 1);
        (*newcell)(0) = left;
        (*newcell)(0).ClearOutputName();
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
        const hwMatrix* lm = nullptr;
        const hwMatrix* rm = nullptr;

        if (left.IsString())
            lm = left.Matrix();
        else
            lm = left.ConvertToMatrix();

        if (right.IsString())
            rm = right.Matrix();
        else
            rm = right.ConvertToMatrix();

        if (!lm || !rm)
            throw OML_Error(HW_ERROR_CONCATMATWSCALCOMPMATORCELL);

        if (left.IsString() || right.IsString())
        {
            if (!lm->IsReal() || !rm->IsReal())
                throw OML_Error(HW_ERROR_NOTCOMPTOSTR);

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
            int lm = l->M();
            int ln = l->N();
            int rm = r->M();
            int rn = r->N();

            if (lm * ln == 0)
                return new StructData(*r);
            else if (rm * rn == 0)
                return new StructData(*l);

            if ((HORIZ && lm != rm) || (!HORIZ && ln != rn))
            {
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

            if (HORIZ)
            {
                result->DimensionNew(lm, ln + rn);

                int i, size = lm * ln;
                for (i = 0; i < size; i++)
                {
                    for (iter = lfields.cbegin(); iter != lfields.cend(); iter++)
                    {
                        result->SetValue(i, -1, iter->first, l->GetValue(i, -1, iter->first));
                    }
                }

                size = lm * (ln + rn);
                for (int j = i; j < size; j++)
                {
                    for (iter = lfields.cbegin(); iter != lfields.cend(); iter++)
                    {
                        result->SetValue(j, -1, iter->first, r->GetValue(j - i, -1, iter->first));
                    }
                }
            }
            else
            {
                result->DimensionNew(lm + rm, ln);

                for (int j = 0; j < ln; j++)
                {
                    for (int i = 0; i < lm; i++)
                    {
                        for (iter = lfields.cbegin(); iter != lfields.cend(); iter++)
                        {
                            result->SetValue(i, j, iter->first, l->GetValue(i, j, iter->first));
                        }
                    }

                    for (int i = 0; i < rm; i++)
                    {
                        for (iter = lfields.cbegin(); iter != lfields.cend(); iter++)
                        {
                            result->SetValue(i + lm, j, iter->first, r->GetValue(i, j, iter->first));
                        }
                    }
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
        current.ClearOutputName();

        for (size_t i = current.IsCellArray() ? 0 : 1; i < nargin; i++)
        {
            current = vconcat(eval, current, inputs[i]);
        }
        outputs.push_back(current);
    }
    else
        outputs.push_back(EvaluatorInterface::allocateMatrix());

    return true;
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
// 
//------------------------------------------------------------------------------
template <int (*F)(int)>
static double wrapLogical(double d)
{
    static const int modby = (const int)pow(2.0, (int) sizeof(char) * 8);
    int i = ((int)d) % modby;

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

    std::string cellname (*input.GetOutputNamePtr());

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
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
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
		FunctionInfo* fi   = inputs[0].FunctionHandle();
		OMLTree*      tree = fi->Statements();

		if (fi->IsAnonymous())
			outputs.push_back(tree->GetStringRepresentation());
		else
			outputs.push_back(fi->FunctionName());

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
	ALT_FUNCPTR   aptr = nullptr;

    if (eval.FindFunctionByName(funcName, &fi, &fptr, &aptr))
    {
		if (fptr)
			fi = new FunctionInfo(funcName, fptr);
		else if (aptr)
			fi = new FunctionInfo(funcName, aptr);
		// else
		//	fi->IncrRefCount(); The ref count will be increased by pushing it into a Currency

        outputs.push_back(fi);
        return true;
    }
	else if (funcName[0] == '@')
	{
		FunctionInfo* fi = eval.FunctionInfoFromString(funcName);

		if (!fi)
			throw OML_Error("Invalid function");
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
// Clears all user-defined functions
//------------------------------------------------------------------------------
bool oml_rehash(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 0)
		throw OML_Error(OML_ERR_NUMARGIN);

	eval.ClearFunctions();
	return true;
}
//------------------------------------------------------------------------------
// Sets the verbosity
//------------------------------------------------------------------------------
bool oml_verbose(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error(OML_ERR_NUMARGIN);

	if (!inputs[0].IsScalar())
		throw OML_Error(OML_ERR_SCALAR);

	eval.SetVerbose((int)inputs[0].Scalar());
	return true;
}
//------------------------------------------------------------------------------
// Evaluates the given function func on elements of given cell array [cellfun]
//------------------------------------------------------------------------------
bool oml_cellfun_nd_helper(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	size_t nargin = inputs.size();

	if (nargin < 2)
		throw OML_Error(OML_ERR_NUMARGIN);

	const Currency& input1 = inputs[0];

	std::string func;

	if (!input1.IsString() && !input1.IsFunctionHandle())
		throw OML_Error(HW_ERROR_FUNCNAMESTR);

	std::vector<Currency> newinputs;
	std::vector<Currency> newoutputs;

	newinputs.push_back(input1);

	int num_inputs         = (int)inputs.size();
	int size_of_each_input = 0;
	
	std::vector<int> dims;

	Currency first_function_input = inputs[1];

	if (first_function_input.IsNDCellArray())
	{
		HML_ND_CELLARRAY* cells = first_function_input.CellArrayND();
		size_of_each_input = cells->Size();
		dims = cells->Dimensions();
	}

	hwMatrixN* result = EvaluatorInterface::allocateMatrixN();
	result->Dimension(dims, hwMatrixN::REAL);

	for (int k = 0; k < size_of_each_input; k++)
	{
		newinputs.clear();
		newoutputs.clear();

		newinputs.push_back(inputs[0]);

		for (size_t j = 1; j < num_inputs; j++)
		{
			Currency function_input = inputs[j];

			if (function_input.IsNDCellArray())
			{
				HML_ND_CELLARRAY* cells = function_input.CellArrayND();

				if (cells->Dimensions() != dims)
					throw OML_Error(HW_ERROR_CELLINPSAMESASIZE);

				newinputs.push_back((*cells)(k));
			}
			else if (function_input.IsCellArray())
			{
				HML_CELLARRAY* cells = function_input.CellArray();

				if (cells->Size() == 1)
					newinputs.push_back((*cells)(0));
				else
					throw OML_Error(HW_ERROR_INPUTALLCELL);
			}
			else
			{
				throw OML_Error(HW_ERROR_INPUTALLCELL);
			}
		}

		// call the function
		oml_feval(eval, newinputs, newoutputs);

		if ((newoutputs.size() == 1) && (newoutputs[0].IsScalar()))
			(*result)(k) = newoutputs[0].Scalar();
		else
			throw OML_Error(HW_ERROR_OUTNOTUNI);
	}

	outputs.push_back(result);

	return true;
}
//------------------------------------------------------------------------------
// Evaluates given function on a cell array [cellfun]
//------------------------------------------------------------------------------
bool oml_cellfun(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)
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
            if (i == 1 && incell->IsEmpty())
            {
                outputs.push_back(EvaluatorInterface::allocateCellArray());
                return true;
            }
            if (incell->Size() != 1)
            {
                cellm = incell->M();
                celln = incell->N();
                break;
            }
        }
		else if (cur.IsNDCellArray())
		{
			// This isn't what I wanted to do, but the logic below is already quite complicated
			// and including this functionality inline would have made it so much worse.
			return oml_cellfun_nd_helper(eval, inputs, outputs);
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
					{
						outputs[j].GetWritableMatrix()->SetElement(i, cur.Scalar());

						if (cur.IsLogical())
							outputs[j].SetMask(Currency::MASK_LOGICAL);
					}
					else if (cur.IsComplex())
					{
						outputs[j].GetWritableMatrix()->SetElement(i, cur.Complex());
					}
					else
					{
						throw OML_Error(HW_ERROR_OUTNOTUNI);
					}
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
        months = rem(months, 12.0);

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
            throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_EMPTYMATRIX));
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
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_EMPTYMATRIX));
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
						throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINDEX));
                    }
                }
                else
                {
					throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINDEX));
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
    if (nargin < 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    const Currency& input1 = inputs[0];
    const Currency& input2 = inputs[1];

    bool isCell1 = input1.IsCellArray();
    bool isCell2 = input2.IsCellArray();

    bool isString1 = input1.IsString();
    bool isString2 = input2.IsString();
    if (isCell1 && !(isString2 || isCell2))
    {
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 2);
    }
    else if (isCell2 && !(isString1 || isCell1))
    {
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
    }

    bool byrows = false;
    if (nargin > 2)
    {
        const Currency& input3 = inputs[2];
        if (!input3.IsString())
        {
            throw OML_Error(OML_ERR_STRING, 3);
        }

        std::string val(input3.StringVal());
        if (!val.empty())
        {
            std::transform(val.begin(), val.end(), val.begin(), ::tolower);
            if (val != "rows")
            {
                throw OML_Error(OML_ERR_BAD_STRING, 3);
            }
            if (isCell1 || isCell2)
            {
                throw OML_Error(OML_Error(OML_ERR_BAD_STRING, 3).GetErrorMessage()
                    + "; 'rows' is not applicable for cell arrays");
            }
            byrows = true;
        }
    }

    int nargout = eval.GetNargoutValue();
    int m       = 0;
    int n       = 0;

    std::unique_ptr<hwMatrix> boolMtx(nullptr);
    std::unique_ptr<hwMatrix> idxMtx(nullptr);

    // Performance improvement for cells inputs
    if (isCell1 && isCell2)
    {
        HML_CELLARRAY* cell1 = input1.CellArray();
        HML_CELLARRAY* cell2 = input2.CellArray();
        int cell1Size = 0;
        if (cell1)
        {
            m = cell1->M();
            n = cell1->N();
            cell1Size = cell1->Size();
        }
        int cell2Size = (cell2) ? cell2->Size() : 0;

        boolMtx.reset(EvaluatorInterface::allocateMatrix(m, n, 0.0));
        if (nargout > 1)
        {
            idxMtx.reset(EvaluatorInterface::allocateMatrix(m, n, 0.0));
        }

        for (int i = 0; i < cell1Size; ++i)
        {
            const Currency& cur1 = (*cell1)(i);
            if (!cur1.IsString())
            {
                throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
            }
            const hwMatrix* mtx1 = cur1.Matrix();
            if (!mtx1)
            {
                continue;
            }

            for (int j = cell2Size; j; --j) // find last instance of
            {
                const Currency& cur2 = (*cell2)(j - 1);
                if (!cur2.IsString())
                {
                    throw OML_Error(OML_ERR_STRING_STRINGCELL, 2);
                }
                if (mtx1->IsEqual(*cur2.Matrix()))
                {
                    (*boolMtx)(i) = 1;
                    if (idxMtx)
                    {
                        (*idxMtx)(i) = j;
                    }
                    break;
                }
            }
        }

        Currency boolcur(boolMtx.release());
        boolcur.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(boolcur);

        if (idxMtx)
        {
            outputs.push_back(idxMtx.release());
        }
        return true;
    }
    else if ((isString1 && isString2) ||
             (isString1 && isCell2)   || 
             (isCell1   && isString2))
    {
        std::vector <Currency> searchfor;
        std::vector <Currency> searchin;

        // push all elements of input1 into searchfor, and all elements of input2 into searchin
        // the strategy generates a single search loop, but is inefficient
        if (byrows && isString1 && isString2)
        {
            n = 1;
            m = stringVecFromCurrencyRows(eval, input1, input2, searchfor, searchin);
        }
        else
        {
            stringVecFromCurrency(eval, input1, input2, searchfor, searchin, &m, &n);
        }

        boolMtx.reset(EvaluatorInterface::allocateMatrix(m, n, 0.0));
        if (nargout > 1)
        {
            idxMtx.reset(EvaluatorInterface::allocateMatrix(m, n, 0.0));
        }

        int searchforsize = static_cast<int>(searchfor.size());
        int searchinsize  = static_cast<int>(searchin.size());

        // do the search
        for (int i = 0; i < searchforsize; ++i)
        {
            const Currency& tofind = searchfor[i];
            for (int j = searchinsize; j; j--) // find last instance of
            {
                if (isequal(tofind, searchin[j-1]))
                {
                    (*boolMtx)(i) = 1;
                    if (idxMtx)
                    {
                        (*idxMtx)(i) = j;
                    }
                    break;
                }
            }
        }

        Currency boolcur(boolMtx.release());
        boolcur.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(boolcur);

        if (idxMtx)
        {
            outputs.push_back(idxMtx.release());
        }
        return true;
    }
    else if (byrows && input1.IsMatrix() && input2.IsMatrix())
    {
        const hwMatrix* in1 = input1.Matrix();
        const hwMatrix* in2 = input2.Matrix();
        assert(in1);
        assert(in2);

        m = in1->M();
        n = in1->N();

        if (n != in2->N())
        {
            std::string msg("Error: invalid inputs in arguments 1 and 2; ");
            msg += "; columns must match when searching by rows";
            throw OML_Error(msg);
        }


        boolMtx.reset(EvaluatorInterface::allocateMatrix(m, 1, 0.0));
        if (nargout > 1)
        {
            idxMtx.reset(EvaluatorInterface::allocateMatrix(m, 1, 0.0));
        }

        int m2 = in2->M();

        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < m2; ++j)
            {
                bool matches = true;
                for (int k = 0; k< n; k++)
                {
                    if ((*in1)(i,k) != (*in2)(j,k))
                    {
                        matches = false;
                        break;
                    }
                }
                if (matches)
                {
                    (*boolMtx)(i) = 1;
                    if (idxMtx)
                    {
                        (*idxMtx)(i) = j + 1;
                    }
                    break;
                }
            }
        }

        Currency boolcur(boolMtx.release());
        boolcur.SetMask(Currency::MASK_LOGICAL);
        outputs.push_back(boolcur);

        if (idxMtx)
        {
            outputs.push_back(idxMtx.release());
        }
        return true;
    }
    else
    {
        static bool unsorted1 = true;

        if (unsorted1 && (input1.IsMatrix() || input1.IsNDMatrix()))
        {
            // input1_col = input1(:);
            std::vector<Currency> colon;
            colon.push_back(Currency(0.0, Currency::TYPE_COLON));
            Currency input1_col = eval.VariableIndex(input1, colon);

            // if (!issorted(input1_col)
            std::vector<Currency> newinputs;
            std::vector<Currency> newoutputs;
            newinputs.push_back(input1_col);
            newinputs.push_back("ascending");

            oml_issorted(eval, newinputs, outputs);

            if (outputs[0].Scalar() == 0)
            {
                // [input1_col_sorted, idx] = sort(input1_col);
                outputs.clear();
                newinputs.pop_back();
                newinputs.push_back("ascend");
                oml_sort(eval, newinputs, outputs);

                // [m_sortedv, m_idxv] = oml_ismember(input1_col_sorted, input2);
                newinputs.clear();
                newinputs.push_back(outputs[0]);
                newinputs.push_back(input2);
                unsorted1 = false;
                // oml_ismember(eval, newinputs, newoutputs);
                newoutputs = eval.DoMultiReturnFunctionCall(oml_ismember, newinputs, 2, 2, true);

                // v_resort(idx) = m_sortedv([1:n]);
                // idx_resort(idx) = m_idxv([1:n]);
                const hwMatrix* idx = outputs[1].Matrix();
                hwMatrix* m_sortedv = newoutputs[0].GetWritableMatrix();
                hwMatrix* m_idxv = newoutputs[1].GetWritableMatrix();

                std::unique_ptr<hwMatrix> m_resort(EvaluatorInterface::allocateMatrix(input1_col.Matrix()->M(), 1, 0.0));
                std::unique_ptr<hwMatrix> m_idx_resort(EvaluatorInterface::allocateMatrix(input1_col.Matrix()->M(), 1, 0.0));
                int size = m_resort->Size();

                for (int i = 0; i < size; ++i)
                {
                    (*m_resort)(static_cast<int> ((*idx)(i))-1) = (*m_sortedv)(i);
                    (*m_idx_resort)(static_cast<int> ((*idx)(i))-1) = (*m_idxv)(i);
                }

                // output[0] = reshape(v_resort, input1_dims)
                // output[1] = reshape(idx_resort, input1_dims)
                std::vector<Currency> dims;
                newinputs.clear();
                newinputs.push_back(input1);
                dims = eval.DoMultiReturnFunctionCall(oml_size, newinputs, 1, 1, true);

                newinputs.clear();
				Currency temp(m_resort.release());
				temp.SetMask(Currency::MASK_LOGICAL);
                newinputs.push_back(temp);
                newinputs.push_back(dims[0]);
                outputs.clear();
                oml_reshape(eval, newinputs, outputs);

                newinputs.clear();
                newinputs.push_back(m_idx_resort.release());
                newinputs.push_back(dims[0]);
                newoutputs.clear();
                oml_reshape(eval, newinputs, newoutputs);
                outputs.push_back(newoutputs[0]);

                return true;
            }

            unsorted1 = false;
            outputs.clear();
        }

        static bool unsorted2 = true;

        if (unsorted2 && (input2.IsMatrix() || input2.IsNDMatrix()))
        {
            // input2_col = input2(:);
            std::vector<Currency> colon;
            colon.push_back(Currency(0.0, Currency::TYPE_COLON));
            Currency input2_col = eval.VariableIndex(input2, colon);

            // if (!issorted(input2_col)
            std::vector<Currency> newinputs;
            std::vector<Currency> newoutputs;
            newinputs.push_back(input2_col);
            newinputs.push_back("ascending");

            oml_issorted(eval, newinputs, newoutputs);

            if (newoutputs[0].Scalar() == 0)
            {
                // [input2_col_sorted, idx] = sort(input2_col);
                newoutputs.clear();
                newinputs.pop_back();
                newinputs.push_back("ascend");
                oml_sort(eval, newinputs, newoutputs);

                // [m_sortedv, m_idxv] = oml_ismember(input1, input2_col_sorted);
                newinputs.clear();
                newinputs.push_back(input1);
                newinputs.push_back(newoutputs[0]);
                unsorted2 = false;
                // oml_ismember(eval, newinputs, outputs);
                outputs = eval.DoMultiReturnFunctionCall(oml_ismember, newinputs, 2, 2, true);

                // v_resort = idx(m_sortedv([1:n]));
                const hwMatrix* idx = newoutputs[1].Matrix();
                const hwMatrix* m_idxv = outputs[1].Matrix();
                hwMatrix* m_idx_resort = EvaluatorInterface::allocateMatrix(m_idxv->M(), m_idxv->N(), 0.0);
                int size = m_idxv->Size();

                for (int i = 0; i < size; ++i)
                {
                    double j = ((*m_idxv)(i));

                    if (static_cast<int> (j) != 0)
                        (*m_idx_resort)(i) = (*idx)(static_cast<int> (j-1));
                }

                outputs.pop_back();
                outputs.push_back(m_idx_resort);

                return true;
            }

            unsorted2 = false;
            outputs.clear();
        }

        // reset flags for next time
        unsorted1 = true;
        unsorted2 = true;

        // get raw pointers for algorithm
        const double* real_in1 = nullptr;
        const double* real_in2 = nullptr;
        const hwComplex* complex_in1 = nullptr;
        const hwComplex* complex_in2 = nullptr;
        int size1;
        int size2;
        hwMatrix* out1 = nullptr;
        hwMatrix* out2 = nullptr;
        hwMatrixN* outN1 = nullptr;
        hwMatrixN* outN2 = nullptr;
        double* real_out1 = nullptr;
        double* real_out2 = nullptr;

        if (input1.IsMatrix() || input1.IsScalar() || input1.IsComplex())
        {
            const hwMatrix* mtx1 = input1.ConvertToMatrix();
            size1 = mtx1->Size();

            if (mtx1->IsReal())
                real_in1 = mtx1->GetRealData();
            else
                complex_in1 = mtx1->GetComplexData();

            out1 = EvaluatorInterface::allocateMatrix(mtx1->M(), mtx1->N(), 0.0);
            real_out1 = out1->GetRealData();

            if (nargout == 2)
            {
                out2 = EvaluatorInterface::allocateMatrix(mtx1->M(), mtx1->N(), 0.0);
                real_out2 = out2->GetRealData();
            }            
        }
        else if (input1.IsNDMatrix())
        {
            const hwMatrixN* mtxN1 = input1.MatrixN();
            size1 = mtxN1->Size();

            if (mtxN1->IsReal())
                real_in1 = mtxN1->GetRealData();
            else
                complex_in1 = mtxN1->GetComplexData();

            outN1 = EvaluatorInterface::allocateMatrixN();
            outN1->Dimension(mtxN1->Dimensions(), hwMatrixN::REAL);
            outN1->SetElements(0.0);
            real_out1 = outN1->GetRealData();

            if (nargout == 2)
            {
                outN2 = EvaluatorInterface::allocateMatrixN();
                outN2->Dimension(mtxN1->Dimensions(), hwMatrixN::REAL);
                outN2->SetElements(0.0);
                real_out2 = outN2->GetRealData();
            }
        }
        else
        {
            throw OML_Error(HW_ERROR_INPUTSTRCELLMTX);
        }

        if (input2.IsMatrix() || input2.IsScalar() || input2.IsComplex())
        {
            const hwMatrix* mtx2 = input2.ConvertToMatrix();
            size2 = mtx2->Size();

            if (mtx2->IsReal())
                real_in2 = mtx2->GetRealData();
            else
                complex_in2 = mtx2->GetComplexData();
        }
        else if (input2.IsNDMatrix())
        {
            const hwMatrixN* mtxN2 = input2.MatrixN();
            size2 = mtxN2->Size();

            if (mtxN2->IsReal())
                real_in2 = mtxN2->GetRealData();
            else
                complex_in2 = mtxN2->GetComplexData();
        }
        else
        {
            throw OML_Error(HW_ERROR_INPUTSTRCELLMTX);
        }

        // search algorithm
        if (real_in1)
        {
            if (real_in2)
            {
                int index = size2-1;

                for (int i = size1-1; i > -1; --i)
                {
                    if (real_in1[i] > real_in2[index])
                        continue;

                    if (real_in1[i] == real_in2[index])
                    {
                        real_out1[i] = 1.0;

                        if (real_out2)
                            real_out2[i] = index + 1.0;
                    }
                    else
                    {
                        if (--index < 0)
                            break;

                        ++i;
                    }
                }
            }
            else if (complex_in2)
            {
                int index = size2-1;
                int sameMagCount = 0;

                for (int i = size1-1; i > -1; --i)
                {
                    if (real_in1[i] > complex_in2[index].Mag())
                        continue;

                    if (complex_in2[index] == real_in1[i])
                    {
                        real_out1[i] = 1.0;

                        if (real_out2)
                            real_out2[i] = index + 1.0;

                        index += sameMagCount;
                        sameMagCount = 0;
                    }
                    else if (complex_in2[index].Mag() == real_in1[i])
                    {
                        if (index == 0)
                            continue;
                            
                        --index;
                        sameMagCount++;
                        ++i;
                    }
                    else
                    {
                        index += sameMagCount;
                        sameMagCount = 0;

                        if (--index < 0)
                            break;

                        ++i;
                    }
                }
            }
        }
        else if (complex_in1)
        {
            if (real_in2)
            {
                int index = size2-1;

                for (int i = size1-1; i > -1; --i)
                {
                    if (complex_in1[i].Mag() > real_in2[index])
                        continue;

                    if (complex_in1[i] == real_in2[index])
                    {
                        real_out1[i] = 1.0;

                        if (real_out2)
                            real_out2[i] = index + 1.0;
                    }
                    else
                    {
                        if (--index < 0)
                            break;

                        ++i;
                    }
                }
            }
            else if (complex_in2)
            {
                int index = size2-1;
                int sameMagCount = 0;

                for (int i = size1-1; i > -1; --i)
                {
                    if (complex_in1[i].Mag() > complex_in2[index].Mag())
                        continue;

                    if (complex_in1[i] == complex_in2[index])
                    {
                        real_out1[i] = 1.0;

                        if (real_out2)
                            real_out2[i] = index + 1.0;

                        index += sameMagCount;
                        sameMagCount = 0;
                    }
                    else if (complex_in1[i].Mag() == complex_in2[index].Mag())
                    {
                        if (index == 0)
                            continue;
                            
                        --index;
                        sameMagCount++;
                        ++i;
                    }
                    else
                    {
                        index += sameMagCount;
                        sameMagCount = 0;

                        if (--index < 0)
                            break;

                        ++i;
                    }
                }
            }
        }

        // pack outputs
        if (out1)
        {
            Currency boolcur(out1);
            boolcur.SetMask(Currency::MASK_LOGICAL);
            outputs.push_back(boolcur);

            if (out2)
                outputs.push_back(out2);
        }
        else if (outN1)
        {
            Currency boolcur(outN1);
            boolcur.SetMask(Currency::MASK_LOGICAL);
            outputs.push_back(boolcur);

            if (outN2)
                outputs.push_back(outN2);
        }
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
		HML_CELLARRAY* cell = in2.CellArray();
		hwMatrix* result = EvaluatorInterface::allocateMatrix(cell->M(), cell->N(), hwMatrix::REAL);

		for (int i = 0; i < cell->Size(); i++)
		{
			if (!(*cell)(i).IsString())
				(*result)(i) = false;
			else
				(*result)(i) = sd->Contains((*cell)(i).StringVal());
		}

		Currency val(result);
		val.SetMask(Currency::MASK_LOGICAL);
		outputs.push_back(val);
	}
	else
	{
		Currency val(false);

		if (in2.IsString())
			val = sd->Contains(in2.StringVal());

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
                if (eval.FindFileInPath(inpstr + ".oml", filename) || eval.FindFileInPath(inpstr + ".m", filename))
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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

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

    if (numinputs < 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    std::FILE* file = nullptr;

    int fileid = getFileFromInput(eval, inputs[0]);
    if (fileid == -1)
    {
        file = stdout;
    }
    else
    {
        if (numinputs == 1)
        {
            throw OML_Error(OML_ERR_NUMARGIN);
        }
        BuiltInFuncsUtils::CheckFileIndex(eval, fileid, 1, true);
        file = eval.GetFile(fileid);
    }

    FileStreamType ftype = determineFileStreamType(file, false);
    std::vector<Currency>::const_iterator itr = inputs.cbegin() + (fileid == -1 ? 0 : 1);

    if (!(*itr).IsString())
    {
        int index = (fileid == -1) ? 1 : 2;
        throw OML_Error(OML_ERR_STRING_FILESTREAM, index, OML_VAR_VALUE);
    }
    std::string result(sprintf(eval, itr, inputs.cend()));

    switch (ftype)
    {
        case Stdout:
        {
            Currency out(result);
            out.PrintfOutput();    // Set flag that this is coming from printf
            eval.PrintResult(out);
            break;
        }

        case Stderr: throw OML_Error(result);     break;
        case Other:  fputs(result.c_str(), file); break;
        default: break;
    }

    if (getNumOutputs(eval))
    {
        outputs.push_back(result.length());
    }
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
// Converts string to number, by evaluating it if needed [str2num]
//------------------------------------------------------------------------------
bool oml_str2num(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    const Currency& input = inputs[0];
    if (!input.IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    if (!input.IsMultilineString()) // Try with strod before evaluating
    {
        std::string str(input.StringVal()); 
        if (str.empty())
        {
            outputs.push_back(EvaluatorInterface::allocateMatrix());
            outputs.push_back(true);
            return true;
        }

        double rval     = 0;
        double ival     = 0;
        bool   isscalar = true;
        BuiltInFuncsString funcs;
        if (funcs.Str2Num(str, rval, ival, isscalar))
        {
            if (isscalar)
            {
                outputs.push_back(rval);
            }
            else
            {
                outputs.push_back(hwComplex(rval, ival));
            }
            outputs.push_back(true);  // Status
            return true;
        }
    }

    // Read line by line, so input.StringVal() cannot be used below
    std::string todo(readString(input)); 

    BuiltInFuncsUtils utils;
    todo = utils.LTrim(todo);  // Trim left space

    if (!todo.empty() && todo[0] != '{')
    {
        todo = '[' + todo + ']';
    }

    Interpreter __interp(eval);
    Currency first = __interp.DoString(todo);
    if (first.IsCellArray())
    {
		// Make sure this is a copy.  Otherwise 'first' and 'outcur' will be holding the
		// same pointer and crash.
        HML_CELLARRAY *out = EvaluatorInterface::allocateCellArray(first.CellArray());
        Currency outcur(out);
        int m = out->M();
        int n = out->N();
        for (int i = 1; i < input.Matrix()->M(); ++i) // Multiline strings
        {
            Currency result = __interp.DoString(readString(input, i));
            if (result.IsCellArray())
            {
                HML_CELLARRAY *cell = result.CellArray();
                if (n != cell->N())
                {
                    outputs.push_back(EvaluatorInterface::allocateMatrix());
                    outputs.push_back(false);
                    return true;
                }
                utils.CheckMathStatus(eval, out->Resize(m + 1 + cell->M(), n));
                for (int j = 0; j < cell->M(); ++j)
                {
                    for (int k = 0; k < n; ++k)
                    {
                        (*out)(j + m++, k) = (*cell)(j, k);
                    }
                }
            }
            else
            {
                outputs.push_back(EvaluatorInterface::allocateMatrix());
                outputs.push_back(false);
                return true;
            }
        }

        outputs.push_back(outcur);
        outputs.push_back(true);
        return true;
    }
    else
    {
        hwMatrix* out = nullptr;
        if (first.IsScalar())
            out = EvaluatorInterface::allocateMatrix(1, 1, first.Scalar());
        else if (first.IsComplex())
            out = EvaluatorInterface::allocateMatrix(1, 1, first.Complex());
        else if (first.IsMatrix())
            out = EvaluatorInterface::allocateMatrix(first.Matrix());
        else
        {
            outputs.push_back(EvaluatorInterface::allocateMatrix());
            outputs.push_back(false);
            return true;
        }

        int currentrow = out->M() - 1;
        int outn = out->N();
        Currency outcur(out);

        for (int i = 1; i < input.Matrix()->M(); ++i)
        {
            std::string todo = '[' + readString(input, i) + ']';
            Currency result = __interp.DoString(todo);
            if (result.IsScalar()|| result.IsComplex())
            {
                if (outn != 1)
                {
                    outputs.push_back(EvaluatorInterface::allocateMatrix());
                    outputs.push_back(false);
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
                    outputs.push_back(false);
                    return true;
                }
                utils.CheckMathStatus(eval, out->Resize(currentrow + 1 + mtx->M(), outn));
                for (int j = 0; j < mtx->M(); ++j)
                {
                    hwMatrix *row = readRow(eval, mtx, j);
                    utils.CheckMathStatus(eval, out->WriteRow(++currentrow, *row));
                    // Memory leak if row is invalid
                    delete row;
                }
            }
            else
            {
                outputs.push_back(EvaluatorInterface::allocateMatrix());
                outputs.push_back(false);
                return true;
            }
        }

        outputs.push_back(outcur);
        outputs.push_back(true);
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

		if (!m->IsEmpty())
		{
			if (m->IsReal())
			{
				val = (*m)((*indexInInput)++);
			}
			else
			{
				val = m->z((*indexInInput)++).Real();
			}
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
    if (IsNaN_T(d) || IsInf_T(d)) 
    {
        return 3;  // "NaN"/"Inf"
    }
    if (IsNegInf_T(d)) 
    {
        return 4;  // "-Inf"
    }

	size_t len = 1;
	
	if (abs(d) > 0)
    {
        len = static_cast<size_t>(abs(log(abs(d)) / log(10.0) + 1));
    }

#ifndef OS_WIN
    if (len >= 18) // Max digits displayed by 
    {
        char dummy[1];
        int  siz = 0;

        siz = snprintf(dummy, sizeof dummy, "%lf", d);
        len = static_cast<size_t>(siz);
    }
#endif

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
// Helper method for sprintf
//------------------------------------------------------------------------------
static std::string dosprintf(size_t size, const std::string &tmplt, ...)
{
    if (size)
    {
        va_list vl;
        va_start(vl, tmplt);

#ifdef OS_WIN
        size_t count    = 0;
        int    new_size = vsnprintf(nullptr, count, tmplt.c_str(), vl);

		if (new_size > size)
			size = new_size;
#endif

        char *buffer = new char[size + 1];
        int numwritten = vsprintf(buffer, tmplt.c_str(), vl);

        std::string rv(buffer);

        delete [] buffer;
        buffer = nullptr;

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
// Helper method for printing
//------------------------------------------------------------------------------
static std::string dosprintf(EvaluatorInterface& eval, const std::string &tmplt, int *indexInInput, FormatType ft,
    std::vector<Currency>::const_iterator &inputIter)
{
    size_t bufferSize = tmplt.length() + tryGetBufferIncrement(tmplt);
    std::string writestr;
    double value = 0;
    switch (ft)
    {
        case String:
            writestr = getWriteStrForSprintf(eval, inputIter, indexInInput);
            bufferSize += writestr.length();
            return dosprintf(bufferSize, tmplt, writestr.c_str());
        case Percent:
            return dosprintf(bufferSize, tmplt);
        case Truncate:
        {
            if (!inputIter->IsEmpty())
            {
                value = getDoubleForSprintf(inputIter, indexInInput);
                if (BuiltInFuncsUtils::NeedsSpecialPrintfFormat(value))
                    return BuiltInFuncsUtils::GetSpecialPrintfFormat(tmplt, value);
                bufferSize += computeLength(true, value);
            }
            else
            {
                bufferSize = 0;
            }
            // Check for very large numbers
            if (fabs(value) > INT_MAX)
            {
                if (fabs(value) > LLONG_MAX)
                {
                    std::string fmt(tmplt);
                    size_t pos = fmt.find("d");
                    if (pos != std::string::npos)
                    {
                        fmt.replace(pos, 1, "e");
                    }
                    return CurrencyDisplay::GetFormattedString(fmt.c_str(), value);
                }
                else
                {
                    // Don't apply formatting
                    return std::to_string(static_cast<long long>(value));
                }
            }
            return dosprintf(bufferSize, tmplt, (long long)value);
        }

        case NoTruncate:
			if (!inputIter->IsEmpty())
			{
				value = getDoubleForSprintf(inputIter, indexInInput);
				if (BuiltInFuncsUtils::NeedsSpecialPrintfFormat(value))
					return BuiltInFuncsUtils::GetSpecialPrintfFormat(tmplt, value);
                bufferSize += computeLength(false, value);
			}
			else
			{
				bufferSize = 0;
			}
			return dosprintf(bufferSize, tmplt, value);

        default:
            throw OML_Error(HW_ERROR_UNSUPFORMTYPE);
    }
}
//------------------------------------------------------------------------------
// Returns a formatted string [sprintf]
//------------------------------------------------------------------------------
bool oml_sprintf(EvaluatorInterface           eval, 
                 const std::vector<Currency>& inputs, 
                 std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }

    std::string tmplt(readString(inputs[0]));

    std::string result (sprintf(
        eval, tmplt, inputs.cbegin() + 1, inputs.cend()));

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
	ALT_FUNCPTR aptr = nullptr;

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
            if (!eval.FindFunctionByName(funcName, &fi, &fptr, &aptr))
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
	else if (aptr)
	{
		throw OML_Error(HW_ERROR_UNSUPOP);
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
// Removes given directories from search path [rmpath]
//------------------------------------------------------------------------------
bool oml_rmhiddenpath(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	size_t nargin = inputs.size();

	if (!nargin)
		throw OML_Error(OML_ERR_NUMARGIN);

	for (int i = 0; i < nargin; i++)
	{
		std::vector<std::string> toremove = separatePathNames(eval, inputs[i]);

		std::vector<std::string>::reverse_iterator iter;
		for (iter = toremove.rbegin(); iter != toremove.rend(); iter++)
		{
			std::string rmv = BuiltInFuncsUtils::GetAbsolutePath(*iter);
			BuiltInFuncsUtils::StripTrailingSlash(rmv);

			if (!eval.RemoveHiddenPath(rmv))
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
        std::vector<std::string> paths = eval.GetPaths();
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
// Adds given path to search path [addpath]
//------------------------------------------------------------------------------
bool oml_addhiddenpath(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	size_t nargin = inputs.size();

	if (!nargin)
		throw OML_Error(OML_ERR_NUMARGIN);

	Currency path = inputs[0];

	if (!path.IsString())
		throw OML_Error(OML_ERR_STRING, 1);

	eval.AddHiddenPath(path.StringVal());

	return true;
}

//------------------------------------------------------------------------------
// Modification to add path, which specifies path and function [addpath2]
//------------------------------------------------------------------------------
bool oml_addpath2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2)
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
				throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
		}
	}

	std::string abs_path = BuiltInFuncsUtils::GetAbsolutePath(new_path.StringVal());
	BuiltInFuncsUtils::StripTrailingSlash(abs_path);

	eval.AddPath2(abs_path, funcs);

	if (nargin == 3)
	{
		if (inputs[2].IsString())
			eval.RegisterLibraryAlias(abs_path, inputs[2].StringVal());
		else
			throw OML_Error(OML_ERR_STRING, 3);
	}
	
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
        std::vector<int> dimR;
        hwMatrixN::DataType type;
        hwMatrixN base;

        if (input1.IsMatrix() || input1.IsScalar() || input1.IsComplex())
        {
            const hwMatrix* mat = input1.ConvertToMatrix();
            dimR.push_back(mat->M() * reps[0]);
            dimR.push_back(mat->N() * reps[1]);

            for (int i = 2; i < reps.size(); ++i)
                dimR.push_back(reps[i]);

            if (mat->IsReal())
                type = hwMatrixN::REAL;
            else
                type = hwMatrixN::COMPLEX;

            base.Convert2DtoND(*mat, false);
        }
        else if (input1.IsNDMatrix())
        {
            const std::vector<int>& dimB = input1.MatrixN()->Dimensions();
            int min = _min(static_cast<int> (reps.size()), static_cast<int> (dimB.size()));

            for (int i = 0; i < min; ++i)
                dimR.push_back(dimB[i] * reps[i]);

            for (int i = min; i < dimB.size(); ++i)
                dimR.push_back(dimB[i]);

            for (int i = min; i < reps.size(); ++i)
                dimR.push_back(reps[i]);

            type = input1.MatrixN()->Type();

            base = (*input1.MatrixN());
        }
        else
        {
            throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
        }

        //// populate result
        hwMatrixN* result = EvaluatorInterface::allocateMatrixN();
        result->Dimension(dimR, type);

        if (base.Size() == 1)
        {
            if (base.IsReal())
                result->SetElements(base(0));
            else
                result->SetElements(base.z(0));
        }
        else
        {
            const std::vector<int>& dimB = base.Dimensions();
            std::vector<int> base_index(dimB.size());
            int numElemsToCopy = dimB[0];
            int numBlocks = base.Size() / numElemsToCopy;

            if (result->IsReal() && !result->IsEmpty())
            {
                // copy base to result
                double* src = base.GetRealData();
                double* dest = result->GetRealData();

                for (int j = 0; j < numBlocks; ++j)
                {
                    if (numElemsToCopy == 1)
                    {
                        (*dest) = (*src);
                    }
                    else
                    {
                        memcpy_s(dest, numElemsToCopy * sizeof(double), src, numElemsToCopy * sizeof(double));
                    }

                    // advance base indices
                    for (int k = 1; k < dimB.size(); ++k)
                    {
                        // increment index k if possible
                        if (base_index[k] != dimB[k] - 1)
                        {
                            ++base_index[k];
                            break;
                        }

                        // index k is maxed out, so reset and continue to k+1
                        base_index[k] = 0;
                    }

                    src = base.GetRealData() + base.Index(base_index);
                    dest = result->GetRealData() + result->Index(base_index);
                }

                // replicate the copied base in all dimensions
                numBlocks = base.Size();
                numElemsToCopy = 1;
                std::vector<int> result_index(dimR.size());

                for (int i = 0; i < reps.size(); ++i)
                {
                    if (i < dimB.size())
                    {
                        numBlocks /= dimB[i];
                        numElemsToCopy *= dimB[i];
                    }

                    double* src = result->GetRealData();

                    for (int j = 0; j < numBlocks; ++j)
                    {
                        double* dest = src;

                        if (numElemsToCopy == 1)
                        {
                            for (int k = 1; k < reps[i]; ++k)
                            {
                                ++dest;
                                (*dest) = (*src);
                            }
                        }
                        else
                        {
                            for (int k = 1; k < reps[i]; ++k)
                            {
                                dest += numElemsToCopy;
                                memcpy_s(dest, numElemsToCopy * sizeof(double), src, numElemsToCopy * sizeof(double));
                                src = dest;
                            }
                        }

                        // advance result indices
                        for (int k = i + 1; k < dimB.size(); ++k)
                        {
                            // increment index k if possible
                            if (result_index[k] < dimB[k] - 1)
                            {
                                ++result_index[k];
                                break;
                            }

                            // index k is maxed out, so reset and continue to k+1
                            result_index[k] = 0;
                        }

                        src = result->GetRealData() + result->Index(result_index);
                    }

                    numElemsToCopy *= reps[i];
                }
            }
            else if (!result->IsEmpty())  // complex
            {
                // copy base to result
                hwComplex* src = base.GetComplexData();
                hwComplex* dest = result->GetComplexData();

                for (int j = 0; j < numBlocks; ++j)
                {
                    if (numElemsToCopy == 1)
                    {
                        (*dest) = (*src);
                    }
                    else
                    {
                        memcpy_s(dest, numElemsToCopy * sizeof(hwComplex), src, numElemsToCopy * sizeof(hwComplex));
                    }

                    // advance base indices
                    for (int k = 1; k < dimB.size(); ++k)
                    {
                        // increment index k if possible
                        if (base_index[k] != dimB[k] - 1)
                        {
                            ++base_index[k];
                            break;
                        }

                        // index k is maxed out, so reset and continue to k+1
                        base_index[k] = 0;
                    }

                    src = base.GetComplexData() + base.Index(base_index);
                    dest = result->GetComplexData() + result->Index(base_index);
                }

                // replicate the copied base in all dimensions
                numBlocks = base.Size();
                numElemsToCopy = 1;
                std::vector<int> result_index(dimR.size());

                for (int i = 0; i < reps.size(); ++i)
                {
                    if (i < dimB.size())
                    {
                        numBlocks /= dimB[i];
                        numElemsToCopy *= dimB[i];
                    }

                    hwComplex* src = result->GetComplexData();

                    for (int j = 0; j < numBlocks; ++j)
                    {
                        hwComplex* dest = src;

                        if (numElemsToCopy == 1)
                        {
                            for (int k = 1; k < reps[i]; ++k)
                            {
                                ++dest;
                                (*dest) = (*src);
                            }
                        }
                        else
                        {
                            for (int k = 1; k < reps[i]; ++k)
                            {
                                dest += numElemsToCopy;
                                memcpy_s(dest, numElemsToCopy * sizeof(hwComplex), src, numElemsToCopy * sizeof(hwComplex));
                                src = dest;
                            }
                        }

                        // advance result indices
                        for (int k = i + 1; k < dimB.size(); ++k)
                        {
                            // increment index k if possible
                            if (result_index[k] < dimB[k] - 1)
                            {
                                ++result_index[k];
                                break;
                            }

                            // index k is maxed out, so reset and continue to k+1
                            result_index[k] = 0;
                        }

                        src = result->GetComplexData() + result->Index(result_index);
                    }

                    numElemsToCopy *= reps[i];
                }
            }
        }

        outputs.push_back(result);

        return true;
    }

    // 2D matrix
    std::unique_ptr<const HML_CELLARRAY> cell;
    const hwMatrix* base;
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
        base = input1.ConvertToMatrix();
        if (!base)
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMTXSTRING);
        m = base->M();
        n = base->N();
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
        hwMatrix *result = EvaluatorInterface::allocateMatrix(m, n, base->Type());
        int M = base->M();
        int N = base->N();

        if (base->Size() == 1)
        {
            if (base->IsReal())
                result->SetElements((*base)(0));
            else
                result->SetElements(base->z(0));
        }
        else if (result->IsReal() && !result->IsEmpty())
        {
            // copy base to result
            int numBytesToCopy = M * sizeof(double);
            double* src = const_cast<double*> (base->GetRealData());
            double* dest = result->GetRealData();

            if (M == 1)
            {
                for (int j = 0; j < N; ++j)
                {
                    (*dest) = src[j];
                    dest += m;
                }
            }
            else
            {
                for (int j = 0; j < N; ++j)
                {
                    memcpy_s(dest, numBytesToCopy, src, numBytesToCopy);
                    src += M;
                    dest += m;
                }
            }

            // replicate the copied base in both dimensions
            for (int j = 0; j < N; ++j)
            {
                src = result->GetRealData() + j * m;
                dest = src;

                if (M == 1)
                {
                    for (int k = 1; k < m / M; ++k)
                    {
                        ++dest;
                        (*dest) = (*src);
                    }
                }
                else
                {
                    for (int k = 1; k < m / M; ++k)
                    {
                        dest += M;
                        memcpy_s(dest, numBytesToCopy, src, numBytesToCopy);
                        src = dest;
                    }
                }
            }

            src = result->GetRealData();
            dest = src;
            numBytesToCopy = m * N * sizeof(double);

            for (int k = 1; k < n / N; ++k)
            {
                dest += m * N;
                memcpy_s(dest, numBytesToCopy, src, numBytesToCopy);
                src = dest;
            }
        }
        else if (!result->IsEmpty())  // complex
        {
            // copy base to result
            int numBytesToCopy = M * sizeof(hwComplex);
            hwComplex* src = const_cast<hwComplex*> (base->GetComplexData());
            hwComplex* dest = result->GetComplexData();

            if (M == 1)
            {
                for (int j = 0; j < N; ++j)
                {
                    (*dest) = src[j];
                    dest += m;
                }
            }
            else
            {
                for (int j = 0; j < N; ++j)
                {
                    memcpy_s(dest, numBytesToCopy, src, numBytesToCopy);
                    src += M;
                    dest += m;
                }
            }

            // replicate the copied base in both dimensions
            for (int j = 0; j < N; ++j)
            {
                src = result->GetComplexData() + j * m;
                dest = src;

                if (M == 1)
                {
                    for (int k = 1; k < m / M; ++k)
                    {
                        ++dest;
                        (*dest) = (*src);
                    }
                }
                else
                {
                    for (int k = 1; k < m / M; ++k)
                    {
                        dest += M;
                        memcpy_s(dest, numBytesToCopy, src, numBytesToCopy);
                        src = dest;
                    }
                }
            }

            src = result->GetComplexData();
            dest = src;
            numBytesToCopy = m * N * sizeof(hwComplex);

            for (int k = 1; k < n / N; ++k)
            {
                dest += m * N;
                memcpy_s(dest, numBytesToCopy, src, numBytesToCopy);
                src = dest;
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
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < m; i++)
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
    else if (input.IsComplex())
    {
        outputs.push_back(input.Complex() * (PI / 180.0));
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
        double scale = PI / 180.0;

        if (mtx->IsReal())
        {
            for (int k = 0; k < mtx->Size(); k++)
                (*result)(k) = (*mtx)(k) * scale;
        }
        else
        {
            for (int k = 0; k < mtx->Size(); k++)
                result->z(k) = mtx->z(k) * scale;
        }

        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_deg2rad);
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, -1, OML_VAR_TYPE);
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
        outputs.push_back(input.Scalar() * (180.0 / PI));
    }
    else if (input.IsComplex())
    {
        outputs.push_back(input.Complex() * (180.0 / PI));
    }
    else if (input.IsMatrix())
    {
        const hwMatrix* mtx = input.Matrix();
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());
        double scale = 180.0 / PI;

        if (mtx->IsReal())
        {
            for (int k = 0; k < mtx->Size(); k++)
                (*result)(k) = (*mtx)(k) * scale;
        }
        else
        {
            for (int k = 0; k < mtx->Size(); k++)
                result->z(k) = mtx->z(k) * scale;
        }

        outputs.push_back(result);
    }
    else if (input.IsNDMatrix())
    {
        return oml_MatrixNUtil1(eval, inputs, outputs, oml_rad2deg);
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX, -1, OML_VAR_TYPE);
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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
        }

        if (!sameSize(x, y))
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

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
        throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
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
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    int fid = BuiltInFuncsUtils::GetFileId(eval, inputs[0], 1);
    BuiltInFuncsUtils::CheckFileIndex(eval, fid, 1, true);

    if (eval.GetFile(fid) == stdin)
    {
       throw OML_Error("Error: invalid file stream; cannot read from stdin; use command [input]");
    }

    // Get the file name or file id

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
bool oml_cell2mat(EvaluatorInterface           eval, 
                  const std::vector<Currency>& inputs, 
                  std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsCellArray())
    {
        throw OML_Error(OML_ERR_CELLARRAY, 1);
    }

    HML_CELLARRAY* cell = inputs[0].CellArray();
    assert(cell);
    if (!cell || cell->Size() == 0)
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        return true;
    }

    enum
    {
        Matrices,
        Cells,
        //Structs,
        Strings
    } inputType;


    // get input type
    const Currency& elem = (*cell)(0, 0);
    if (elem.IsScalar() || elem.IsComplex())
    {
        if (cell->Size() == 1)
        {
            outputs.push_back(elem);
            return true;
        }
        inputType = Matrices;
    }
    else if (elem.IsMatrix())
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
    else
    {
        throw OML_Error(OML_ERR_SCAL_COMP_MTX_STR_CELL, 1);
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
                if (mlist)
                {
                    delete [] mlist;
                    mlist = nullptr;
                }

				throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINPUT));
            }
        }

        if (setm)
        {
            if (tempn != n)
            {
                if (mlist)
                {
                    delete [] mlist;
                    mlist = nullptr;
                }
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
                if (mlist)
                {
                    delete [] mlist;
                    mlist = nullptr;
                }
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
        if (mlist)
        {
            delete [] mlist;
            mlist = nullptr;
        }
        throw OML_Error(HW_ERROR_UNEVENDIMENSIONS);
    }

    if (inputType == Cells)
    {
        if (mlist)
        {
            delete [] mlist;
            mlist = nullptr;
        }

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
                    if (mlist)
                    {
                        delete [] mlist;
                        mlist = nullptr;
                    }
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
                    if (mlist)
                    {
                        delete [] mlist;
                        mlist = nullptr;
                    }
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
                        if (mlist)
                        {
                            delete [] mlist;
                            mlist = nullptr;
                        }
                        throw OML_Error(HW_ERROR_MIXEDCELLELEMS);
                    }
                }
                else
                {
                    if (mlist)
                    {
                        delete [] mlist;
                        mlist = nullptr;
                    }
                    throw OML_Error(HW_ERROR_MIXEDCELLELEMS);
                }
                break;
            default:
                if (mlist)
                {
                    delete [] mlist;
                    mlist = nullptr;
                }
                throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
            }
        }
    }
    outputs.push_back(outCur);
    if (mlist)
    {
        delete [] mlist;
        mlist = nullptr;
    }

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
        if (eval.FindFileInPath(filename + ".oml", filename) || eval.FindFileInPath(filename + ".m", filename))
            return true;
    }
    return eval.FindFileInPath(filename, filename);
}
//------------------------------------------------------------------------------
// Executes given file/command [run]
//------------------------------------------------------------------------------
bool oml_run(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency curstr = toCurrencyStr(eval, inputs[0], false, false);
    std::string filename(readString(curstr));

    BuiltInFuncsUtils utils;
    if (!utils.FileExists(filename))
    {
        if (!FileExistsInPath(eval, filename))  // Check for utf8
        {
            throw OML_Error(OML_ERR_FILE_NOTFOUND, 1);
        }
    }

    std::string full_file_name = BuiltInFuncsUtils::GetAbsolutePath(filename);
    std::string cwd = BuiltInFuncsUtils::GetCurrentWorkingDir();
	std::string pathname;
#ifdef OS_WIN
    size_t slash_index = full_file_name.find_last_of("/\\");
#else
    size_t slash_index = full_file_name.rfind('/');
#endif

    // extract path so we can cd to it
    if (slash_index != std::string::npos)
    {
        pathname = std::string(full_file_name.begin(), full_file_name.begin() + slash_index);
        cd(pathname, eval);
    }

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

    // Create a child intepreter as this is in the middle of an eval
    Interpreter interp(eval);
    SignalHandlerBase* childHandler  = interp.GetSignalHandler();
    SignalHandlerBase* parentHandler = eval.GetSignalHandler();
    if (!childHandler)
    {
        childHandler = parentHandler;
        interp.SetSignalHandler(childHandler);
    }
    if (childHandler)
    {
        // Signals like exit, cd need to be copied but not printing as all
        // printing is handled in parent, when results are transfered over
        childHandler->CopyNonPrintSignals();
    }

    // Evaluate the try string
	eval.RegisterChildEvaluator(interp.GetEvaluator());

#if 0
	// Line below is correct, but it breaks other products in combination with 
    // other changes
	int old_val = eval_ptr->SetNestedFunctionMarker(0);
#endif

    bool wasopen = BuiltInFuncsUtils::IsOutputLogOpen();
    if (parentHandler)
    {
        // Disconnect printing till results are copied over. Otherwise this 
        // results in printing being done twice in the child and parent in batch
        parentHandler->DisconnectPrintSignal();
    }

    try
    {
        Currency ret;
        
        if (!utils.HasWideChars(full_file_name))
        {
            ret = interp.DoFile(full_file_name);
        }
        else // ANTLR does not support running of files with wide characters
        {
            std::string contents(utils.GetFileContents(full_file_name));
            ret = interp.DoString(contents);
        }
        if (ret.IsError())
        {
            errmsg = ret.Message();
        }
    }
    catch (OML_Error& e)
    {
        errmsg    = e.GetErrorMessage();
        formatmsg = e.GetFormatMessage();
    }

    if (parentHandler)
    {
        // Reconnect printing since all execution in the child is complete
        parentHandler->ConnectPrintSignal();
    }

    if (!wasopen)  // If the log was opened in the child evaluator, close it
    {
        BuiltInFuncsUtils::CloseOutputLog();
    }

	eval.RemoveChildEvaluator();

    if (slash_index != std::string::npos)
    {
		std::string cwd2 = BuiltInFuncsUtils::GetCurrentWorkingDir();

		if (cwd2 == pathname) // if the run script changed the cwd, leave it
			cd(cwd, eval);
    }

    // Push results
    bool logerror = false;

    std::vector<Currency> results (interp.GetOutputCurrencyList());
    for (std::vector<Currency>::const_iterator itr = results.begin(); 
         itr != results.end(); ++itr)
    {
        const Currency& cur = *itr;
        if (cur.IsError())
        {
            errmsg = cur.Message(); // To display later
            if (cur.IsLogOutput())
            {
                logerror = true;
            }
        }
        else
        {
            eval.PushResult(cur);
        }
    }

    // Output errors
    if (!errmsg.empty())
    {
        if (logerror)
        {
            BuiltInFuncsUtils::OpenOutputLogForAppend();
            BuiltInFuncsUtils::SaveToOutputLog(errmsg);
            BuiltInFuncsUtils::CloseOutputLog();
        }
        throw OML_Error(errmsg, formatmsg);
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
    BuiltInFuncsUtils::CheckFileIndex(eval, fileID, 1, true);

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
    BuiltInFuncsUtils::CheckFileIndex(eval, fileID, 1, false);

    std::FILE *f = eval.GetFile(fileID);
    rewind(f);

    if (getNumOutputs(eval))
        outputs.push_back(ferror(f) ? -1.0 : 0.0);

    return true;
}
//------------------------------------------------------------------------------
// Gets precision from string
//------------------------------------------------------------------------------
static Precision getPrecision(const std::string& in)
{
    std::string str(in);

    // no support yet for conversion types
    int blockSize = 1;
    size_t index = str.find('*');
    if (index != std::string::npos)
    {
        std::string left (str.substr(0, index));
        const char* cleft = left.c_str();
        char* end = nullptr;
        double temp = strtod(cleft, &end);
        if (end && end > cleft)
        {
            if (!isint(temp) || temp < 1.0)
            {
                throw OML_Error(HW_ERROR_PRECBLOCKPOSINT);
            }
            blockSize = (int) temp;
            str.erase(0, index + 1);
        }
    }
    
    index = str.find("=>");
    if (index != std::string::npos)
    {
        str = str.substr(0, index);
    }

    if (!str.empty())
    {
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    }
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
static Precision getPrecision(EvaluatorInterface& eval, const Currency& input)
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
    BuiltInFuncsUtils::CheckFileIndex(eval, fileID, 1, true);
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
            throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
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
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
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

    BuiltInFuncsUtils::CheckFileIndex(eval, fileID, 1, true);
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
    BuiltInFuncsUtils::CheckFileIndex(eval, fileID, 1, true);
    bool flushcout = BuiltInFuncsUtils::IsFlushCout(eval, fileID);
    FILE* f = eval.GetFile(fileID);
    if (feof(f))
    {
        outputs.push_back(getTrue());
    }
    else
    {
        char c = fgetc(f);
        if (flushcout)
        {
            std::cout << std::flush;
        }
        if (c == EOF)
        {
            outputs.push_back(getTrue());
        }
        else
        {
            ungetc(c,f);
            if (flushcout)
            {
                std::cout << std::flush;
            }
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
    BuiltInFuncsUtils::CheckFileIndex(eval, fileID, 1, true);

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
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    int fid = BuiltInFuncsUtils::GetFileId(eval, inputs[0], 1);
    BuiltInFuncsUtils::CheckFileIndex(eval, fid, 1, true);

    if (eval.GetFile(fid) == stdin)
    {
        throw OML_Error("Error: invalid file stream; cannot read from stdin; use command [input]");
    }


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
    BuiltInFuncsUtils::CheckFileIndex(eval, fileID, 1, false);
    outputs.push_back(eval.CloseFile(fileID) ? 0.0 : 1.0);
    
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
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINPUT));
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
                collapseDelimiters = (inputs[2].Scalar()) ? true : false;

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
                collapseDelimiters = (inputs[2].Scalar()) ? true : false;

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
            collapseDelimiters = (input2.Scalar()) ? true : false;

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
                {
                    collapseDelimiters = (value.Scalar()) ? true : false;
                }
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
            throw OML_Error(OML_ERR_STRING, i + 1);
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
        int total_size = 0;

        for (i = 0; i < cellsize; i++)
        {
            Currency elem = (*cell)(0);
            if (elem.IsString())
            {
                const hwMatrix *mtx = elem.Matrix();
                int m = mtx->M();

                if (!m)
                    continue;

                total_size += mtx->N();

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

        if (!delimstr.IsEmpty() && (delimstr.M() != 1) && (delimstr.M() != numToConcat))
            throw OML_Error(HW_ERROR_DELCELLINVAMTOFELE);

        total_size += delimstr.Size() * (cellsize - 1);

        //allocate one matrix of total_size, then use CopyBlock to fill it in

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
                    overlap = (input5.Scalar()) ? true : false;
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
            throw OML_Error(OML_ERR_STRING, 4);
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
        {
            allowOverlap = (input4.Scalar()) ? true : false;
        }
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
    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    bool compareRows = false;
    bool forward = false;


    const Currency &input1 = inputs[0];
    if (nargin > 1)
    {
        bool direcSet = false;
        if (!inputs[1].IsString())
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);

        std::string i2(readOption(eval, inputs[1]));
        if (i2 == "first")
            forward = direcSet = true;
        else if (i2 == "rows")
            compareRows = true;
        else if (i2 == "last")
            direcSet = true;
        else
            throw OML_Error("Error: invalid input in argument 2; must be 'first', 'last' or 'rows'");

        if (nargin > 2)
        {
            if (!inputs[2].IsString())
                throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);

            std::string i3 (readOption(eval, inputs[2]));
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
                throw OML_Error("Error: invalid input in argument 3; must be 'first', 'last' or 'rows'");
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
		std::vector<int> dims;

		for (int j = 0; j < size; j++)
		{
			Currency input = inputs[j];

			if (input.IsPositiveInteger())
			{
				m = (int)input.Scalar();

				if (m < 0)
					throw OML_Error(OML_ERR_INTEGER, j, OML_VAR_VALUE);

				dims.push_back(m);
			}
			else
			{
				throw OML_Error(OML_ERR_INTEGER, j, OML_VAR_VALUE);
			}
		}

		HML_ND_CELLARRAY* ret = new HML_ND_CELLARRAY(dims, HML_ND_CELLARRAY::REAL);
		outputs.push_back(ret);
		return true;
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
				throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_EMPTYMATRIX));
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
		throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_EMPTYMATRIX));

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

        BuiltInFuncsUtils::CheckMathStatus(eval, x->Eigen(true, V.get(), *e));
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

            if (input2.IsInteger())
			{
                m = static_cast<int>(input2.Scalar());

	            if (m < 0)
		            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);
			}
            else if (input2.IsMatrix() && !input2.Matrix()->M() && !input2.Matrix()->N())
			{
                m = -1;
			}
            else
			{
                throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);
			}

            if (input3.IsInteger())
			{
                n = static_cast<int>(input3.Scalar());

				if (n < 0)
					throw OML_Error(OML_ERR_NATURALNUM, 3, OML_VAR_DIM);
			}
            else if (input3.IsMatrix() && !input3.Matrix()->M() && !input3.Matrix()->N())
			{
                n = -1;
			}
            else
			{
                throw OML_Error(OML_ERR_NATURALNUM, 3, OML_VAR_DIM);
			}
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

                if (cur.IsInteger())
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
        else if (input1.IsMatrix() || input1.IsScalar() || input1.IsComplex() || input1.IsString())
        {
            const hwMatrix* temp = input1.ConvertToMatrix();
            std::unique_ptr<hwMatrixN> mtxN(EvaluatorInterface::allocateMatrixN());
            mtxN->Convert2DtoND(*temp);

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
		else
		{
			throw OML_Error(OML_ERR_MATRIX, 1);
		}
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the generalized transpose of a matrix A using the vector P [permute]
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

    int n = pvec->Size();
    std::vector<int> permvec(n);

    for (int i = 0; i < n; ++i)
    {
        double dim = realval(pvec, i);
        int idim = static_cast<int>(dim) - 1;

        if (idim < 0)
        {
            throw OML_Error(hwMathStatus(HW_MATH_ERR_PERMVEC2, 2));
        }

        permvec[i] = idim;
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
        if (!input1.IsMatrix() && !input1.IsScalar() && !input1.IsComplex() && !input1.IsString())
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

        if (input1.IsMatrix() || input1.IsString())
        {
            hwMatrixN temp;

            temp.Convert2DtoND(*input1.Matrix(), false);

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

    outputs[0].SetMask(input1.GetMask());

    return true;
}
//------------------------------------------------------------------------------
// Returns the reversed permutation of a matrix A using the vector P [ipermute]
//------------------------------------------------------------------------------
bool oml_ipermute(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& input2 = inputs[1];

    // reverse the permutation and call permute
    if (!input2.IsMatrix())
    {
        throw OML_Error(OML_ERR_POSINTVECTOR, 2, OML_VAR_DIMS);
    }

    const hwMatrix* ipvec = input2.Matrix();

    if (!ipvec->IsVector() || !ipvec->IsRealData())
    {
        throw OML_Error(OML_ERR_POSINTVECTOR, 2, OML_VAR_DIMS);
    }

    int n = ipvec->Size();
    std::unique_ptr<hwMatrix> permvec(EvaluatorInterface::allocateMatrix(n, 1, hwMatrix::REAL));

    for (int i = 0; i < n; ++i)
    {
        double dim = realval(ipvec, i);
        int idim = static_cast<int>(dim) - 1;

        if (idim < 0)
        {
            throw OML_Error(hwMathStatus(HW_MATH_ERR_PERMVEC2, 2));
        }

        if (idim > n - 1)
        {
            throw OML_Error(hwMathStatus(HW_MATH_ERR_PERMVEC3, 2));
        }

        (*permvec)(idim) = static_cast<double>(i + 1);
    }

    std::vector<Currency> inputs2;
    inputs2.push_back(inputs[0]);
    inputs2.push_back(permvec.release());
    
    return oml_permute(eval, inputs2, outputs);
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
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];

    if (!input1.IsMatrix() && !input1.IsScalar() && !input1.IsComplex())
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    const hwMatrix *mtx = input1.ConvertToMatrix();

    // unpack arguments 2 and 3
    double p = 2;
    std::string type;   // "inf", "-inf" or "fro"
    std::string opt;

    if (nargin > 1)
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
            std::string strarg2 = readOption(eval, input2);

            if (strarg2 == "inf" || strarg2 == "fro" || strarg2 == "-inf")
            {
                type = strarg2;
            }
            else if (strarg2 == "rows" || strarg2 == "cols" || strarg2 == "columns")
            {
                opt = strarg2;
            }
            else
            {
                throw OML_Error(OML_ERR_BAD_STRING, 2);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_SCALARSTRING);
        }

        if (p == std::numeric_limits<double>::infinity())
        {
            type = "inf";
        }
        else if (p == -1 * std::numeric_limits<double>::infinity())
        {
            type = "-inf";
        }
    }

    if (nargin > 2)
    {
        const Currency& input3 = inputs[2];

        if (input3.IsString())
        {
            std::string strarg3 = readOption(eval, input3);

            if (strarg3 == "rows" || strarg3 == "cols" || strarg3 == "columns")
            {
                if (!opt.empty())
                {
                    throw OML_Error(OML_ERR_BAD_STRING, 3);
                }

                opt = strarg3;
            }
            else
            {
                throw OML_Error(OML_ERR_NORM_STRING3, 3);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_STRING, 3);
        }
    }

    if (p < 0.0 && type != "-inf")
        throw OML_Error(OML_ERR_PNORM, 2);

    // process the options
    double norm;

    if (opt.empty())
    {
        if (!type.empty())  // type = "inf", "-inf" or "fro"
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Norm(norm, type.c_str()));
        }
        else if (isint(p)) // integer p
        {
            if (mtx->IsVector())
            {
                if (static_cast<int> (p) == 0)
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

                    norm = static_cast<double> (count);
                }
                else
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Norm(norm, static_cast<int> (p)));
                }
            }
            else // matrix norm
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Norm(norm, static_cast<int> (p)));
            }
        }
        else // non-integer p
        {
            if (mtx->IsVector())
            {
                norm = normForVector(mtx, p);
            }
            else if (mtx->M() == 0 && mtx->N() == 0)
            {
                norm = 0.0;
            }
            else
            {
                throw OML_Error(OML_ERR_INTSTRING, 2, OML_VAR_VALUE); // requires constrained optimization
            }
        }

        outputs.push_back(norm);
    }
    else if (opt == "rows")
    {
        Currency resultCur;
        hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), 1, hwMatrix::REAL);
        std::unique_ptr<hwMatrix> row(EvaluatorInterface::allocateMatrix());

        for (int i = 0; i < mtx->M(); i++)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->ReadRow(i, *row));

            if (!type.empty())  // type = "inf", "-inf" or "fro"
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, row->Norm(norm, type.c_str()));
            }
            else if (isint(p))
            {
                if (static_cast<int> (p) == 0)
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

                    norm = static_cast<double> (count);
                }
                else
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, row.get()->Norm(norm, static_cast<int> (p)));
                }
            }
            else // non-integer p
            {
                norm = normForVector(row.get(), p);
            }

            result->SetElement(i, norm);
        }

        outputs.push_back(result);
    }
    else if (opt == "cols" || opt == "columns")
    {
        Currency resultCur;
        hwMatrix* result = EvaluatorInterface::allocateMatrix(1, mtx->N(), hwMatrix::REAL);

        for (int i = 0; i < mtx->N(); i++)
        {
            std::unique_ptr<const hwMatrix> col(EvaluatorInterface::allocateColumn(mtx, i));

            if (!type.empty())  // type = "inf", "-inf" or "fro"
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, col->Norm(norm, type.c_str()));
            }
            else if (isint(p)) // integer p
            {
                if (static_cast<int> (p) == 0)
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

                    norm = static_cast<double> (count);
                }
                else
                {
                    BuiltInFuncsUtils::CheckMathStatus(eval, col.get()->Norm(norm, static_cast<int> (p)));
                }
            }
            else // non-integer p
            {
                norm = normForVector(col.get(), p);
            }

            result->SetElement(i, norm);
        }

        outputs.push_back(result);
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
    size_t nargin = inputs.size();
    int    nargout = eval.GetNargoutValue();

    if (nargin == 0 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &input1 = inputs[0];

    if (input1.IsMatrix() || input1.IsScalar() || input1.IsComplex() || input1.IsString())
    {
        const hwMatrix* findin = findin = input1.ConvertToMatrix();
        int incr = 1;
        double stopAt = 0;

        if (nargin > 1)
        {
            if (!inputs[1].IsScalar())
                throw OML_Error(OML_ERR_REAL, 2, OML_VAR_VALUE);

            stopAt = inputs[1].Scalar();
            if (stopAt < 0 || !isint(stopAt))
            {
                throw OML_Error(OML_ERR_POSINTEGER, 2);
            }
            if (nargin > 2)
            {
                const Currency& input3 = inputs[2];
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
                    if (nargin > 1 && ivec.size() >= stopAt)
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
                    if (nargin > 1 && ivec.size() >= stopAt)
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
                for (int j = 0; j < nargout; j++)
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
            hwMatrix* indices;
            if (findin->M() == 1)
                indices = EvaluatorInterface::allocateMatrix(1, (int)vecsize, hwMatrix::REAL);
            else
                indices = EvaluatorInterface::allocateMatrix((int)vecsize, 1, hwMatrix::REAL);

            for (int i = 0; i < vecsize; i++)
            {
                int idx = ivec[incr == 1 ? i : vecsize - i - 1];
                (*indices)(i) = static_cast<double>(idx + 1);
            }
            outputs.push_back(indices);
        }
        else
        {
            hwMatrix* rows, *cols, *vals;
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
    }
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

        for (int i = 1; i < nargin; ++i)
            inputs2.push_back(inputs[i]);

        return oml_find(eval, inputs2, outputs);
    }
    else if (input1.IsSparse())
    {
        const hwMatrixS* sparse = input1.MatrixS();
        int nnz;

        if (nargin > 1)
        {
            if (!inputs[1].IsPositiveInteger())
                throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_VALUE);

            nnz = static_cast<int> (inputs[1].Scalar());
        }
        else
        {
            nnz = sparse->NNZ();
        }

        int first;
        int last;

        if (nargin > 2)
        {
            const Currency& input3 = inputs[2];
            if (!input3.IsString())
                throw OML_Error(OML_ERR_STRING, 3, OML_VAR_TYPE);

            std::string direc = readOption(eval, input3);
            if (direc == "first")
            {
                first = 0;
                last = min(nnz, sparse->NNZ()) - 1;
            }
            else if (direc == "last")
            {
                last = sparse->NNZ() - 1;
                first = max(0, last - nnz + 1);
            }
            else
            {
                throw OML_Error(HW_ERROR_DIRECT1STLAST);
            }
        }
        else
        {
            first = 0;
            last = min(nnz, sparse->NNZ()) - 1;
        }

        if (nargout < 2)
        {
            std::vector<int> indexI;

            sparse->NZinfo(first, last, indexI);

            nnz = static_cast<int> (indexI.size());
            hwMatrix* index = EvaluatorInterface::allocateMatrix(nnz, 1, hwMatrix::REAL);

            for (int i = 0; i < nnz; ++i)
                (*index)(i) = static_cast<double> (indexI[i] + 1);

            if (sparse->M() == 1)
            {
                index->Transpose();
            }

            outputs.push_back(index);
        }
        else
        {
            std::vector<int> rowI;
            std::vector<int> colI;
            std::unique_ptr<hwMatrix> val(EvaluatorInterface::allocateMatrix());

            sparse->NZinfo(first, last, rowI, colI, *val);

            nnz = val->Size();
            hwMatrix* row = EvaluatorInterface::allocateMatrix(nnz, 1, hwMatrix::REAL);
            hwMatrix* col = EvaluatorInterface::allocateMatrix(nnz, 1, hwMatrix::REAL);

            for (int i = 0; i < nnz; ++i)
            {
                (*row)(i) = static_cast<double> (rowI[i] + 1);
                (*col)(i) = static_cast<double> (colI[i] + 1);
            }

            if (sparse->M() == 1)
            {
                row->Transpose();
                col->Transpose();
                val->Transpose();
            }

            outputs.push_back(row);
            outputs.push_back(col);

            if (nargout == 3)
                outputs.push_back(val.release());
        }
    }
    else
        throw OML_Error(OML_ERR_SCALARCOMPLEXMTX);

    return true;
}
//------------------------------------------------------------------------------
// Evaluates given input command [eval]
//------------------------------------------------------------------------------
bool oml_eval(EvaluatorInterface           eval, 
              const std::vector<Currency>& inputs, 
              std::vector<Currency>&       outputs)
{
    size_t size = inputs.size();
    if (size != 1 && size != 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }
    std::string trystr(readString(inputs[0]));
    if (trystr.empty())
    {
        return true;
    }

    bool hasSemiColon    = (trystr.back() == ';');
    bool outputSomething = (getNumOutputs(eval) > 0);
    Currency result;

    // Try if input can be pre-parsed and evaluated without ANTLR
    if (trystr.find(' ') == std::string::npos)
    {
        if (hasSemiColon)
        {
            trystr.pop_back();
        }

		bool evaluated   = false;

		result = eval.GetValue(trystr);
		if (!result.IsError())
		{
			evaluated = !result.IsNothing();
		}

		if (!evaluated)
			evaluated = BuiltInFuncsString::IsNumber(trystr, result);  // strtod

        if (!evaluated)
        {
            if (!evaluated && !trystr.empty() && trystr.length() > 1 &&
                trystr[0] == '\'' && trystr.back() == '\'')   // String with multiple leading/trailing single quotes
            {
                std::string tmp(BuiltInFuncsUtils::LTrim(trystr, "'"));
                tmp = BuiltInFuncsUtils::RTrim(tmp, "'");
                if (tmp.find('\'') == std::string::npos)
                {
                    result = Currency(tmp);
                    evaluated = true;
                }
            }
        }
        
        if (evaluated)
        {
            if (!hasSemiColon || outputSomething)
            {
                outputs.push_back(result);
            }
            return true;
        }
        if (hasSemiColon)
        {
            trystr += ';';
        }
    }
   


    // Create a child intepreter as this is in the middle of an eval
    Interpreter interp1(eval);
    SignalHandlerBase* handler1 = interp1.GetSignalHandler();
    if (handler1)
    {
        handler1->CopyNonPrintSignals(); // Printing will be handled in eval.PushResult
    }

    // Evaluate the try string
    bool storeSupressedResult = true;
    eval.CacheLineInfomation();
    eval.RegisterChildEvaluator(interp1.GetEvaluator());
    result = interp1.DoString(trystr, storeSupressedResult);
    std::vector<Currency> results (interp1.GetOutputCurrencyList());
    eval.RemoveChildEvaluator();

    std::string errmsg;
    if (result.IsError())
    {
        if (size < 2)
        {
            errmsg = result.Message();
        }
        else // Evaluate the catch string
        {
            if (!inputs[1].IsString())
            {
                throw OML_Error(OML_ERR_STRING, 2);
            }

            std::string catchstr (readString(inputs[1]));
            if (!catchstr.empty())
            {
                Interpreter interp2(eval); // Create a second interpreter so outputs don't get mixed up
                SignalHandlerBase* handler2 = interp2.GetSignalHandler();
                if (handler2)
                {
                    handler2->CopyNonPrintSignals(); 
                }
                hasSemiColon = (catchstr.back() == ';');
                result = interp2.DoString(catchstr, storeSupressedResult);
                if (result.IsError())
                {
                    errmsg = result.Message();
                }
                results = interp2.GetOutputCurrencyList();
            }
            else
            {
                results.clear();
            }
        }
    }

    if (!results.empty())
    {
        for (std::vector<Currency>::const_iterator itr = results.begin();
             itr != results.end(); ++itr)
        {
            const Currency& r = *itr;
            if (r.IsError())
            {
                throw OML_Error(r.Message());
            }
            if ((!hasSemiColon || outputSomething) && itr == results.end()-1)
            {
                outputs.push_back(r);
            }
            else
            {
                eval.PushResult(r);
            }
        }        
    }
    else if (!result.IsError()) // If the try/catch strings have semicolon, results are empty
    {
        // Copy the result over only if there is an output and there is an explicit supression of output
        // If this has results and it is not suppressed, we would not be in this section of code
        if (hasSemiColon && outputSomething)
        {
            outputs.push_back(result);
        }
    }
    
    if (!errmsg.empty())
    {
        throw OML_Error(errmsg);
    }

	eval.UncacheLineInfomation();
    return true;
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

		if (!cells || cells->IsEmpty())
			outputs.push_back(0.0);
		else
		   outputs.push_back(max(cells->M(), cells->N()));
    }
    else if (input.IsMatrix())
    {
        const hwMatrix *mtx = input.Matrix();

        if (!mtx || mtx->IsEmpty())
            outputs.push_back(0.0);
        else
            outputs.push_back(max(mtx->M(), mtx->N()));
    }
    else if (input.IsNDMatrix())
    {
        const hwMatrixN *mtx = input.MatrixN();

        if (!mtx || mtx->IsEmpty())
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
    else if (input.IsSparse())
    {
        const hwMatrixS *mtx = input.MatrixS();

        if (!mtx || mtx->IsEmpty())
            outputs.push_back(0.0);
        else
            outputs.push_back(max(mtx->M(), mtx->N()));
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
    else if (input.IsStruct() || input.IsObject())
    {
        StructData *sd = input.Struct();
        if (!sd)
        {
            outputs.push_back(0);
        }
        else
        {
            outputs.push_back((double)max(sd->M(), sd->N()));
        }
    }
	else if (input.IsPointer())
	{
		Currency ref = *input.Pointer();

		if (ref.IsObject())
		{
			StructData* sd = ref.Struct();
            if (!sd)
            {
                outputs.push_back(0);
            }
            else
            {
                outputs.push_back(sd->M() * sd->N());
            }
		}
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
	else if (input1.IsNDMatrix() || input1.IsNDCellArray())
	{
		std::vector<int> dims;

		if (input1.IsNDMatrix())
		{
			const hwMatrixN* mtxn = input1.MatrixN();
			dims = mtxn->Dimensions();
		}
		else
		{
			HML_ND_CELLARRAY* cells = input1.CellArrayND();
			dims = cells->Dimensions();
		}

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
    else if (input1.IsSparse())
        dosize(input1.MatrixS(), dim, (int)size, nargout, outputs);
    else if (input1.IsCellArray())
        dosize(input1.CellArray(), dim, (int) size, nargout, outputs);
    else if (input1.IsStruct() || input1.IsObject())
        dosize(input1.Struct(), dim, (int) size, nargout, outputs);
	else if (input1.IsPointer() && input1.Pointer()->IsObject())
		dosize(input1.Pointer()->Struct(), dim, (int)size, nargout, outputs);
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
    if (inputs[0].IsNDCellArray())
    {
        HML_ND_CELLARRAY* cells_nd = inputs[0].CellArrayND();
        const std::vector<int> dims = cells_nd->Dimensions();
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

    int nargout = getNumOutputs(eval);

    if (nargout > 2)
        throw OML_Error(OML_ERR_NUMARGOUT);

    const Currency &input = inputs[0];

    if (input.IsMatrix() || input.IsScalar() || input.IsComplex())
    {
        const hwMatrix *mtx = input.ConvertToMatrix();
        hwMatrix *result = EvaluatorInterface::allocateMatrix();
        hwMathStatus status = result->Inverse(*mtx);
        double rcond;

        if (!status.IsOk() && !status.IsWarning())  // error
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, status);
        }

        if (nargout < 2)
        {
            if (!status.IsWarning())
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, mtx->RCond(rcond));

                if (rcond < MACHEP2)
                {
                    status(HW_MATH_WARN_SINGMATRIX, 1);
                }
            }

            BuiltInFuncsUtils::CheckMathStatus(eval, status);
            outputs.push_back(result);
        }
        else
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->RCond(rcond));
            outputs.push_back(result);
            outputs.push_back(rcond);
        }
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
bool oml_pwd(EvaluatorInterface           eval, 
             const std::vector<Currency>& inputs, 
             std::vector<Currency>&       outputs)
{
    BuiltInFuncsUtils utils;
    std::string       pwd;

#ifdef OS_WIN
    std::wstring wstr = utils.GetCurrentWorkingDirW();
    pwd = utils.WString2StdString(wstr);
#else
    pwd = utils.GetCurrentWorkingDir();
#endif
    outputs.push_back(pwd);
    return true;
}
//------------------------------------------------------------------------------
// Returns the current date and time a matrix [clock]
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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
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
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
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
                    throw OML_Error(OML_ERR_BAD_STRING, 3);
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

        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_NOTIMPLEMENT));

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
// Returns the least common multiple of the inputs
//------------------------------------------------------------------------------
bool oml_lcm(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	int max_val = 0;

	if (inputs.size() == 1)
	{
		if (inputs[0].IsInteger())
		{
			outputs.push_back(inputs[0].Scalar());
			return true;
		}
		else
		{
			throw OML_Error(OML_ERR_INTEGER, 1);
		}
	}

	for (int j=0; j<inputs.size(); j++)
	{
		const Currency& temp = inputs[j];

		if (temp.IsScalar())
		{
			double val = temp.Scalar();

			if (val < 0)
				throw OML_Error(OML_ERR_NATURALNUM, j + 1);

			if (val == 0)
			{
				outputs.push_back(0);
				return true;
			}

			if (!temp.IsInteger())
				throw OML_Error(OML_ERR_INTEGER, j + 1);

			if (temp.Scalar() > INT_MAX)
				throw OML_Error(OML_ERR_INVALID_RANGE, j + 1);

			if (temp.Scalar() > max_val)
				max_val = (int)temp.Scalar();
		}
		else
		{
			throw OML_Error(OML_ERR_SCALAR, j + 1);
		}
	}

	std::deque<int> primes = getprimes(max_val);

	std::vector<std::map<int, int> > factor_list;

	for (int j=0; j<inputs.size(); j++)
	{
		std::map<int, int> factors;

		Currency temp_cur = inputs[j];

		int test_val = (int)temp_cur.Scalar();

		for (int k=0; k<primes.size(); k++)
		{
			int factor = primes[k];
			int count  = 0;

			if ((test_val % factor) == 0)
			{
				count = 1;

				bool keep_going = true;

				int new_val = test_val;

				while (keep_going)
				{
					new_val /= factor;

					if ((new_val % factor) == 0)
						count++;
					else
						keep_going = false;
				}

				factors[factor] = count;
			}
		}

		factor_list.push_back(factors);
	}

	std::map<int, int> result;
	std::map<int, int>::iterator iter;

	for (int j=0; j<factor_list.size(); j++)
	{
		std::map<int, int> cur_map = factor_list[j];

		for (iter = cur_map.begin(); iter != cur_map.end(); iter++)
		{
			int factor = iter->first;
			int multiplicity = iter->second;

			std::map<int, int>::iterator iter2 = result.find(factor);

			if (iter2 == result.end())
			{
				result[factor] = multiplicity;
			}
			else
			{
				if (multiplicity > iter2->second)
					result[factor] = multiplicity;
			}
		}
	}

	int product = 1;

	for (iter = result.begin(); iter != result.end(); iter++)
			product *= (int)pow((double)iter->first, iter->second);

	outputs.push_back(product);

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
                    BuiltInFuncsMKL::Hypot(*m1, *m2, *result);
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
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
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

            if (!(checkisfinite(dim)))
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
    size_t nargin = inputs.size();
    int nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency input1 = inputs[0];
    int dim = 1;

    if (nargin == 2)
    {
        // handle cases with no dimension argument
        Currency input2 = inputs[1];

        if (input1.IsScalar())
        {
            if (input2.IsScalar())
            {
                if (IsNaN_T(input1.Scalar()))
                    outputs.push_back(input2.Scalar());
                else
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
            else if (input2.IsMatrix() || input2.IsString())
            {
                double dbl = input1.Scalar();
                const hwMatrix* mtx = input2.Matrix();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

                if (mtx->IsReal())
                {
                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        result->SetElement(k, min(dbl, (*mtx)(k)));
                    }
                }
                else
                {
                    hwComplex val;

                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        val = mtx->z(k);
                        if (complexLessThan(hwComplex(dbl, 0.0), val))
                            result->SetElement(k, dbl);
                        else
                            result->SetElement(k, val);
                    }
                }

                outputs.push_back(result);
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
            else if (input2.IsMatrix() || input2.IsString())
            {
                hwComplex cplx = input1.Complex();
                const hwMatrix* mtx = input2.Matrix();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::COMPLEX);

                if (mtx->IsReal())
                {
                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        if (complexLessThan(cplx, hwComplex((*mtx)(k), 0.0)))
                            result->SetElement(k, cplx);
                        else
                            result->SetElement(k, (*mtx)(k));
                    }
                }
                else
                {
                    hwComplex val;

                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        val = mtx->z(k);
                        if (complexLessThan(cplx, val))
                            result->SetElement(k, cplx);
                        else
                            result->SetElement(k, val);
                    }
                }

                outputs.push_back(result);
            }
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else if (input1.IsMatrix() || input1.IsString())
        {
            if (input2.IsScalar())
            {
                const hwMatrix* mtx = input1.Matrix();
                double dbl = input2.Scalar();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

                if (mtx->IsReal())
                {
                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        result->SetElement(k, min((*mtx)(k), dbl));
                    }
                }
                else
                {
                    hwComplex val;

                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        val = mtx->z(k);
                        if (complexLessThan(val, hwComplex(dbl, 0.0)))
                            result->SetElement(k, val);
                        else
                            result->SetElement(k, dbl);
                    }
                }

                outputs.push_back(result);
            }
            else if (input2.IsComplex())
            {
                const hwMatrix* mtx = input1.Matrix();
                hwComplex cplx = input2.Complex();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::COMPLEX);

                if (mtx->IsReal())
                {
                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        if (complexLessThan(hwComplex((*mtx)(k), 0.0), cplx))
                            result->SetElement(k, (*mtx)(k));
                        else
                            result->SetElement(k, cplx);
                    }
                }
                else
                {
                    hwComplex val;

                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        val = mtx->z(k);
                        if (complexLessThan(val, cplx))
                            result->SetElement(k, val);
                        else
                            result->SetElement(k, cplx);
                    }
                }

                outputs.push_back(result);
            }
            else if (input2.IsMatrix() || input2.IsString())
            {
                const hwMatrix* i1 = input1.Matrix();
                const hwMatrix* i2 = input2.Matrix();

                if (i1->M() == i2->M() && i1->N() == i2->N())
                {
                    hwMatrix* result = EvaluatorInterface::allocateMatrix(i1->M(), i1->N(), i1->Type());

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
                else
                {
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                }
            }
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else if (input1.IsNDMatrix() || input2.IsNDMatrix())
        {
            return oml_MatrixNUtil2(eval, inputs, outputs, oml_min);
        }
        else
        {
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
        }

        return true;
    }
    else if (nargin == 1)
    {
        if (input1.IsMatrix())
        {
            const hwMatrix* mtx = input1.Matrix();

            if (mtx->IsVector())
            {
                if (mtx->M() == 1)
                    dim = 2;
                else
                    dim = 1;
            }
        }
        else if (input1.IsString())
        {
            dim = 2;
        }
    }
    else if (nargin == 3)
    {
        Currency input2 = inputs[1];
        Currency input3 = inputs[2];

        if (!input2.IsMatrix())
            throw OML_Error(OML_ERR_EMPTYMATRIX, 2, OML_VAR_INPUT);

        if (!input2.Matrix()->Is0x0())
            throw OML_Error(OML_ERR_EMPTYMATRIX, 2, OML_VAR_INPUT);

        if (!input3.IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);

        dim = static_cast<int>(input3.Scalar());
    }

    // handle cases with a dimension argument
    if (input1.IsScalar() || input1.IsComplex())
    {
        outputs.push_back(input1);
        outputs.push_back(1.0);
    }
    else if (input1.IsMatrix() || input1.IsString())
    {
        const hwMatrix* mtx = input1.Matrix();
        int index = 0;
        if (mtx->IsVector() && (mtx->M() == 1 && dim == 2) || (mtx->N() == 1 && dim == 1))
        {
			if (mtx->IsEmpty())
			{
				outputs.push_back(EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL));
			}
            else if (mtx->IsReal())
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
        else if (dim == 1)
        {
            if (mtx->M() == 0)
            {
                outputs.push_back(EvaluatorInterface::allocateMatrix(0, mtx->N(), hwMatrix::REAL)); // maximums
                outputs.push_back(EvaluatorInterface::allocateMatrix(0, mtx->N(), hwMatrix::REAL)); // indices
            }
            else
            {
                hwMatrix* result = EvaluatorInterface::allocateMatrix(1, mtx->N(), mtx->Type());
                hwMatrix* indices = EvaluatorInterface::allocateMatrix(1, mtx->N(), mtx->Type());

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
        else if (dim == 2)
        {
            const hwMatrix* i1 = input1.Matrix();
            hwMatrix* result;
            int m = i1->M(), n = i1->N();

            if (n > 0)
                result = EvaluatorInterface::allocateMatrix(m, 1, i1->Type());
            else
                result = EvaluatorInterface::allocateMatrix(m, 0, i1->Type());

            hwMatrix* indices = EvaluatorInterface::allocateMatrix(result->M(), result->N(), hwMatrix::REAL);

            if (result->IsEmpty())
            {
                outputs.push_back(result);
                outputs.push_back(indices);
                return true;
            }

            if (i1->IsReal())
            {
                double loc_min;
                int min_index;
                for (int i = 0; i < m; ++i)
                {
                    min_index = 0;
                    loc_min = (*i1)(i, 0);

                    for (int j = 1; j < n; ++j)
                    {
                        double val = (*i1)(i, j);
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
                for (int i = 0; i < m; ++i)
                {
                    min_index = 0;
                    min = i1->z(i, 0);

                    for (int j = 1; j < n; ++j)
                    {
                        hwComplex val = i1->z(i, j);
                        if (complexLessThan(val, min))
                        {
                            min = val;
                            min_index = j;
                        }
                    }
                    result->SetElement(i, min);
                    (*indices)(i) = min_index + 1.0;
                }
            }

            outputs.push_back(result);
            outputs.push_back(indices);
        }
        else    // dim > 2
        {
            outputs.push_back(input1);
        }
    }
    else if (input1.IsNDMatrix())
    {
        if (nargin == 1)
            return oml_MatrixNUtil3(eval, inputs, outputs, oml_min);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, oml_min, 3);
    }
    else if (input1.IsSparse())
    {
        const hwMatrixS* spInput = inputs[0].MatrixS();
        hwMatrixS* spOutput = new hwMatrixS;
        hwMatrixI* indexI = nullptr;

        if (nargout > 1)
            indexI = new hwMatrixI;

        if (dim == 1)
            spInput->Min(*spOutput, indexI, true);
        else
            spInput->Min(*spOutput, indexI, false);

        outputs.push_back(spOutput);

        if (nargout > 1)
        {
            hwMatrix* index = EvaluatorInterface::allocateMatrix(indexI->M(), indexI->N(), hwMatrix::REAL);

            for (int i = 0; i < index->Size(); ++i)
                (*index)(i) = static_cast<double> ((*indexI)(i) + 1);

            delete indexI;
            outputs.push_back(index);
        }
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
    }

    return true;
}
//------------------------------------------------------------------------------
// Returns the max value in the given input [max]
//------------------------------------------------------------------------------
bool oml_max(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    int nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargout > 2)
        throw OML_Error(OML_ERR_NUMARGOUT);

    Currency input1 = inputs[0];
    int dim = 1;

    if (nargin == 2)
    {
        // handle cases with no dimension argument
        Currency input2 = inputs[1];

        if (input1.IsScalar())
        {
            if (input2.IsScalar())
            {
                if (IsNaN_T(input1.Scalar()))
                    outputs.push_back(input2.Scalar());
                else
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
            else if (input2.IsMatrix() || input2.IsString())
            {
                double dbl = input1.Scalar();
                const hwMatrix* mtx = input2.Matrix();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

                if (mtx->IsReal())
                {
                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        result->SetElement(k, max(dbl, (*mtx)(k)));
                    }
                }
                else
                {
                    hwComplex val;

                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        val = mtx->z(k);
                        if (complexGreaterThan(hwComplex(dbl, 0.0), val))
                            result->SetElement(k, dbl);
                        else
                            result->SetElement(k, val);
                    }
                }

                outputs.push_back(result);
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
            else if (input2.IsMatrix() || input2.IsString())
            {
                hwComplex cplx = input1.Complex();
                const hwMatrix* mtx = input2.Matrix();
                hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::COMPLEX);

                if (mtx->IsReal())
                {
                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        if (complexGreaterThan(cplx, hwComplex((*mtx)(k), 0.0)))
                            result->SetElement(k, cplx);
                        else
                            result->SetElement(k, (*mtx)(k));
                    }
                }
                else
                {
                    hwComplex val;

                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        val = mtx->z(k);
                        if (complexGreaterThan(cplx, val))
                            result->SetElement(k, cplx);
                        else
                            result->SetElement(k, val);
                    }
                }

                outputs.push_back(result);
            }
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else if (input1.IsMatrix() || input1.IsString())
        {
            if (input2.IsScalar())
            {
                const hwMatrix* mtx = input1.Matrix();
                double dbl = input2.Scalar();

                hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), mtx->Type());

                if (mtx->IsReal())
                {
                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        result->SetElement(k, max((*mtx)(k), dbl));
                    }
                }
                else
                {
                    hwComplex val;

                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        val = mtx->z(k);
                        if (complexGreaterThan(val, hwComplex(dbl, 0.0)))
                            result->SetElement(k, val);
                        else
                            result->SetElement(k, dbl);
                    }
                }

                outputs.push_back(result);
            }
            else if (input2.IsComplex())
            {
                const hwMatrix* mtx = input1.Matrix();
                hwComplex cplx = input2.Complex();

                hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::COMPLEX);

                if (mtx->IsReal())
                {
                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        if (complexGreaterThan(hwComplex((*mtx)(k), 0.0), cplx))
                            result->SetElement(k, (*mtx)(k));
                        else
                            result->SetElement(k, cplx);
                    }
                }
                else
                {
                    hwComplex val;

                    for (int k = 0; k < mtx->Size(); ++k)
                    {
                        val = mtx->z(k);
                        if (complexGreaterThan(val, cplx))
                            result->SetElement(k, val);
                        else
                            result->SetElement(k, cplx);
                    }
                }

                outputs.push_back(result);
            }
            else if (input2.IsMatrix() || input2.IsString())
            {
                const hwMatrix* i1 = input1.Matrix();
                const hwMatrix* i2 = input2.Matrix();

                if (i1->M() == i2->M() && i1->N() == i2->N())
                {
                    hwMatrix* result = EvaluatorInterface::allocateMatrix(i1->M(), i1->N(), i1->Type());

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
                else
                {
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
                }
            }
            else
            {
                throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
            }
        }
        else if (input1.IsNDMatrix() || input2.IsNDMatrix())
        {
            return oml_MatrixNUtil2(eval, inputs, outputs, oml_max);
        }
        else
        {
            throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
        }

        return true;
    }
    else if (nargin == 1)
    {
        if (input1.IsMatrix())
        {
            const hwMatrix* mtx = input1.Matrix();

            if (mtx->IsVector())
            {
                if (mtx->M() == 1)
                    dim = 2;
                else
                    dim = 1;
            }
        }
        else if (input1.IsString())
        {
            dim = 2;
        }
    }
    else // if (nargin == 3)
    {
        Currency input2 = inputs[1];
        Currency input3 = inputs[2];

        if (!input2.IsMatrix())
            throw OML_Error(OML_ERR_EMPTYMATRIX, 2, OML_VAR_INPUT);

        if (!input2.Matrix()->Is0x0())
            throw OML_Error(OML_ERR_EMPTYMATRIX, 2, OML_VAR_INPUT);

        if (!input3.IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);

        dim = static_cast<int>(input3.Scalar());
    }

    // handle cases with a dimension argument
    if (input1.IsScalar() || input1.IsComplex())
    {
        outputs.push_back(input1);
        outputs.push_back(1.0);
    }
    else if (input1.IsMatrix() || input1.IsString())
    {
        const hwMatrix* mtx = input1.Matrix();
        int index = 0;
        if (mtx->IsVector() && (mtx->M() == 1 && dim == 2) || (mtx->N() == 1 && dim == 1))
        {
			if (mtx->IsEmpty())
			{
				outputs.push_back(EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL));
			}
            else if (mtx->IsReal())
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
        else if (dim == 1)
        {
            if (mtx->M() == 0)
            {
                outputs.push_back(EvaluatorInterface::allocateMatrix(0, mtx->N(), hwMatrix::REAL)); // maximums
                outputs.push_back(EvaluatorInterface::allocateMatrix(0, mtx->N(), hwMatrix::REAL)); // indices
            }
            else
            {
                hwMatrix* result = EvaluatorInterface::allocateMatrix(1, mtx->N(), mtx->Type());
                hwMatrix* indices = EvaluatorInterface::allocateMatrix(1, mtx->N(), mtx->Type());

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
        else if (dim == 2)
        {
            const hwMatrix* i1 = input1.Matrix();
            hwMatrix* result;
            int m = i1->M(), n = i1->N();

            if (n > 0)
                result = EvaluatorInterface::allocateMatrix(m, 1, i1->Type());
            else
                result = EvaluatorInterface::allocateMatrix(m, 0, i1->Type());

            hwMatrix* indices = EvaluatorInterface::allocateMatrix(result->M(), result->N(), hwMatrix::REAL);

            if (result->IsEmpty())
            {
                outputs.push_back(result);
                outputs.push_back(indices);
                return true;
            }

            if (i1->IsReal())
            {
                double loc_max;
                int max_index;
                for (int i = 0; i < m; ++i)
                {
                    max_index = 0;
                    loc_max = (*i1)(i, 0);

                    for (int j = 1; j < n; ++j)
                    {
                        double val = (*i1)(i, j);
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
                for (int i = 0; i < m; ++i)
                {
                    max_index = 0;
                    max = i1->z(i, 0);

                    for (int j = 1; j < n; ++j)
                    {
                        hwComplex val = i1->z(i, j);
                        if (complexGreaterThan(val, max))
                        {
                            max = val;
                            max_index = j;
                        }
                    }
                    result->SetElement(i, max);
                    (*indices)(i) = max_index + 1.0;
                }
            }

            outputs.push_back(result);
            outputs.push_back(indices);
        }
        else    // dim > 2
        {
            outputs.push_back(input1);
        }
    }
    else if (input1.IsNDMatrix())
    {
        if (nargin == 1)
            return oml_MatrixNUtil3(eval, inputs, outputs, oml_max);
        else
            return oml_MatrixNUtil3(eval, inputs, outputs, oml_max, 3);
    }
    else if (input1.IsSparse())
    {
        const hwMatrixS* spInput = inputs[0].MatrixS();
        hwMatrixS* spOutput = new hwMatrixS;
        hwMatrixI* indexI = nullptr;

        if (nargout > 1)
            indexI = new hwMatrixI;

        if (dim == 1)
            spInput->Max(*spOutput, indexI, true);
        else
            spInput->Max(*spOutput, indexI, false);

        outputs.push_back(spOutput);

        if (nargout > 1)
        {
            hwMatrix* index = EvaluatorInterface::allocateMatrix(indexI->M(), indexI->N(), hwMatrix::REAL);

            for (int i = 0; i < index->Size(); ++i)
                (*index)(i) = static_cast<double> ((*indexI)(i) + 1);

            delete indexI;
            outputs.push_back(index);
        }
    }
    else
    {
        throw OML_Error(HW_ERROR_INPUTSCALARCOMPLEXMATRIX);
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
            throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_NOTIMPLEMENT));
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

            BuiltInFuncsUtils::CheckMathStatus(eval, result->Resize(max(0, m), max(0, n), false));
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
bool oml_memoryuse(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    std::string function_def;
    function_def = 
        "function memorystring = memoryused()\n "
            "_pid = getpid();\n "
            "if ispc()\n "
                "cmd = sprintf('tasklist /FI \"PID eq %d\" /NH /V | sort /+150', _pid);\n "
                "[R1, R2] = system(cmd);\n "
                "data = strsplit(R2);\n "
                "memory = cell2mat(data(5));\n "
            "else\n "
                "cmd = sprintf('ps aux | awk ''$2 ~ /%d/ {print $6}''', _pid);\n "
                "[R1, R2] = system(cmd);\n "
                "memory = R2(1:length(R2)-1);\n "   // trim trailing linefeed 
            "end\n "
            "memorystring = sprintf('Memory Used: %s KB', memory);\n "
        "end\n ";    Interpreter __interp(eval);
    Currency first = __interp.DoString(function_def);
    Currency second = __interp.DoString("memoryused");
    if (second.IsString())
    {
        outputs.push_back(second);
    }
    return true;
}
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
    int nargin = eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsLogical() || (nargin > 1 && inputs[1].IsLogical()))
        throw OML_Error(HW_ERROR_INPUTISLOGICAL);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_VARIABLE);

    const hwMatrix* i1 = inputs[0].ConvertToMatrix();
    bool balance = false;

    if (nargin == 2)
    {
        if (inputs[1].IsString())
        {
            if (inputs[1].StringVal() == "balance")
            {
                balance = true;
            }
            else if (inputs[1].StringVal() != "nobalance")
            {
                throw OML_Error(HW_ERROR_INVBALFLAG);
            }
        }
    }

    std::unique_ptr<hwMatrix> D(EvaluatorInterface::allocateMatrix());
    std::unique_ptr<hwMatrix> V = nullptr;

    if (nargout == 2)
    {
        V.reset(EvaluatorInterface::allocateMatrix());
    }
    else if (nargout > 2)
    {
        throw OML_Error(OML_ERR_NUMARGOUT);
    }

    hwMathStatus status;

    if (nargin == 1 || (nargin == 2 && inputs[1].IsString()))
    {
        hwTMatrix<double> T;

        if (i1->IsHermitian() && i1->Csky(T).IsOk())
            status = i1->EigenSH(V.get(), *D);
        else
            status = i1->Eigen(balance, V.get(), *D);
    }
    else
    {
        if (!inputs[1].IsMatrix() && !inputs[1].IsScalar() && !inputs[1].IsComplex() && !inputs[1].IsString())
            throw OML_Error(OML_ERR_MATRIX, 2, OML_VAR_DATA);

        const hwMatrix* i2 = inputs[1].ConvertToMatrix();

        if (sameSize(i1, i2))
            status = i1->Eigen(*i1, *i2, V.get(), *D);
        else
            throw OML_Error(HW_ERROR_INPMUSTSAMESIZE);
    }

    BuiltInFuncsUtils::CheckMathStatus(eval, status);
    if (nargout <= 1)
    {
        outputs.push_back(D.release());
    }
    else
    {
        outputs.push_back(V.release());
        outputs.push_back(vectorToDiag(D.get()));
    }

    return true;
}
//------------------------------------------------------------------------------
// Creates diagonal matrix or returns the diagonal elements of input [diag]
//------------------------------------------------------------------------------
bool oml_diag(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t size = inputs.size();

    if (size < 1 || size > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    // check input matrix
    const Currency& input1 = inputs[0];
    const hwMatrix* i1 = nullptr;

    if (input1.IsMatrix() || input1.IsScalar() || input1.IsComplex() || input1.IsString())
    {
        i1 = input1.ConvertToMatrix();
    }
    else if (!input1.IsSparse())
    {
        throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);
    }

    // check dimensions
    int i2, i3;

    if (size == 1)
    {
        i2 = 0;
    }
    else
    {
        const Currency& input2 = inputs[1];

        if (!input2.IsScalar())
            throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_VALUE);

        if (!IsInteger(input2.Scalar()).IsOk())
            throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_VALUE);

        i2 = static_cast<int> (input2.Scalar());

        if (size == 3)
        {
            const Currency& input3 = inputs[2];

            if (i1)
            {
                if (!i1->IsVector())
                    throw OML_Error(HW_ERROR_3INP1STVEC);
            }
            else
            {
                const hwMatrixS* spm = input1.MatrixS();

                if (!spm->IsVector())
                    throw OML_Error(HW_ERROR_3INP1STVEC);
            }

            if (!input3.IsScalar())
                throw OML_Error(OML_ERR_NATURALNUM, 3, OML_VAR_DIM);

            if (!IsInteger(input3.Scalar()).IsOk())
                throw OML_Error(OML_ERR_NATURALNUM, 3, OML_VAR_DIM);

            i3 = static_cast<int> (input3.Scalar());
        
            if (i2 < 0)
                throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);

            if (i3 < 0)
                throw OML_Error(OML_ERR_NATURALNUM, 3, OML_VAR_DIM);
        }
    }

    // now calculate/create diagonal
    if (i1)
    {
        hwMatrix* result = EvaluatorInterface::allocateMatrix();
        Currency resultCur(result);
        hwMathStatus status;

        if (size < 3)
        {
            status = result->Diag(*i1, i2);
        }
        else // size == 3
        {
            status = result->Diag(*i1, 0);
            BuiltInFuncsUtils::CheckMathStatus(eval, status);

            // adjust size to be i2 x i3 with new elements as 0
            status = result->Resize(i2, i3, true);
        }

        BuiltInFuncsUtils::CheckMathStatus(eval, status);

        if (input1.IsString())
            resultCur.SetMask(Currency::MASK_STRING);

        outputs.push_back(resultCur);
    }
    else    // sparse
    {
        const hwMatrixS* spm = input1.MatrixS();
        hwMatrixS* result = new hwMatrixS;

        if (size < 3)
        {
            result->Diag(*spm, i2);
        }
        else // size == 3
        {
            result->Diag(*spm, i2, i3);
        }

        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Creates a sparse matrix [sparse]
//------------------------------------------------------------------------------
bool oml_sparse(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int nargin = static_cast<int> (inputs.size());

    if (nargin < 1 || nargin == 4 || nargin > 6)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin == 1)
    {
        if (inputs[0].IsSparse())
        {
            outputs.push_back(inputs[0]);
            return true;
        }

        if (!inputs[0].IsMatrix() && !inputs[0].IsScalar() && !inputs[0].IsComplex())
            throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_DATA);

        const hwMatrix* full = inputs[0].ConvertToMatrix();
        hwMatrixS* sparse = new hwMatrixS(*full);
        outputs.push_back(sparse);

        return true;
    }

    const hwMatrix* ivec = nullptr;
    const hwMatrix* jvec = nullptr;
    const hwMatrix* vals = nullptr;
    int             m = -1;
    int             n = -1;
    std::string     option;

    if (nargin == 2)
    {
        if (!inputs[0].IsInteger())
            throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_VALUE);

        if (!inputs[1].IsInteger())
            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_VALUE);

        m = static_cast<int> (inputs[0].Scalar());
        n = static_cast<int> (inputs[1].Scalar());

        if (m < 0)
            throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_VALUE);

        if (n < 0)
            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_VALUE);
    }

    if (nargin > 2)
    {
        if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
            throw OML_Error(OML_ERR_POSINTVECTOR, 1, OML_VAR_INDEX);

        if (!inputs[1].IsMatrix() && !inputs[1].IsScalar())
            throw OML_Error(OML_ERR_POSINTVECTOR, 2, OML_VAR_INDEX);

        if (!inputs[2].IsMatrix() && !inputs[2].IsScalar())
            throw OML_Error(OML_ERR_VECTOR, 3, OML_VAR_DATA);

        ivec = inputs[0].ConvertToMatrix();
        jvec = inputs[1].ConvertToMatrix();
        vals = inputs[2].ConvertToMatrix();

        if (!ivec->IsReal())
            throw OML_Error(OML_ERR_POSINTVECTOR, 1, OML_VAR_INDEX);

        if (!jvec->IsReal())
            throw OML_Error(OML_ERR_POSINTVECTOR, 2, OML_VAR_INDEX);
    }

    if (nargin > 4)
    {
        if (inputs[3].IsInteger())
        {
            m = static_cast<int> (inputs[3].Scalar());

            if (m < 0)
                throw OML_Error(OML_ERR_NATURALNUM, 4, OML_VAR_VALUE);
        }
        else if (inputs[3].IsMatrix())
        {
            if (!inputs[3].Matrix()->Is0x0())
                throw OML_Error(OML_ERR_NATURALNUM, 4, OML_VAR_VALUE);
        }
        else
        {
            throw OML_Error(OML_ERR_NATURALNUM, 4, OML_VAR_VALUE);
        }

        if (inputs[4].IsInteger())
        {
            n = static_cast<int> (inputs[4].Scalar());

            if (n < 0)
                throw OML_Error(OML_ERR_NATURALNUM, 5, OML_VAR_VALUE);
        }
        else if (inputs[4].IsMatrix())
        {
            if (!inputs[4].Matrix()->Is0x0())
                throw OML_Error(OML_ERR_NATURALNUM, 5, OML_VAR_VALUE);
        }
        else
        {
            throw OML_Error(OML_ERR_NATURALNUM, 5, OML_VAR_VALUE);
        }
    }

    if (nargin > 5)
    {
        if (!inputs[5].IsString())
            throw OML_Error(OML_ERR_POSINTEGER, 4, OML_VAR_VALUE);

        option = inputs[5].StringVal();
    }

    std::vector<int> iVec;
    std::vector<int> jVec;

    if (ivec && jvec)
    {
        iVec.resize(ivec->Size());
        jVec.resize(jvec->Size());

        for (int i = 0; i < ivec->Size(); ++i)
        {
            if (!IsInteger( (*ivec)(i)).IsOk() )
                throw OML_Error(OML_ERR_POSINTVECTOR, 1, OML_VAR_INDEX);

            iVec[i] = static_cast<int>((*ivec)(i)) - 1;
        }

        for (int i = 0; i < jvec->Size(); ++i)
        {
            if (!IsInteger((*jvec)(i)).IsOk())
                throw OML_Error(OML_ERR_POSINTVECTOR, 2, OML_VAR_INDEX);

            jVec[i] = static_cast<int>((*jvec)(i)) - 1;
        }
    }

    bool emptyVals = false;

    if (!vals)
    {
        vals = new hwMatrix;
        emptyVals = true;
    }

    hwMatrixS* sparse = new hwMatrixS(iVec, jVec, *vals, m, n, option.c_str());
    outputs.push_back(sparse);

    if (emptyVals)
        delete vals;

    return true;
}
//------------------------------------------------------------------------------
// Returns true if input is a sparse matrix [issparse]
//------------------------------------------------------------------------------
bool oml_issparse(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];
    outputs.push_back(HW_MAKE_BOOL_CURRENCY(cur.IsSparse()));
    return true;
}
//------------------------------------------------------------------------------
// Returns the number of non-zero elements in a matrix [nnz]
//------------------------------------------------------------------------------
bool oml_nnz(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& cur = inputs[0];
    int nnz = 0;

    if (cur.IsMatrix())
    {
        const hwMatrix* mat = cur.Matrix();

        if (mat->IsReal())
        {
            for (int i = 0; i < mat->Size(); ++i)
            {
                if ((*mat)(i) != 0.0)
                    ++nnz;
            }
        }
        else
        {
            for (int i = 0; i < mat->Size(); ++i)
            {
                if (mat->z(i) != 0.0)
                    ++nnz;
            }
        }
    }
    else if (cur.IsSparse())
    {
        const hwMatrixS* mat = cur.MatrixS();
        nnz = mat->NNZ();
    }
    else if (cur.IsNDMatrix())
    {
        const hwMatrixN* mat = cur.MatrixN();

        if (mat->IsReal())
        {
            for (int i = 0; i < mat->Size(); ++i)
            {
                if ((*mat)(i) != 0.0)
                    ++nnz;
            }
        }
        else
        {
            for (int i = 0; i < mat->Size(); ++i)
            {
                if (mat->z(i) != 0.0)
                    ++nnz;
            }
        }
    }
    else if (cur.IsScalar())
    {
        if (cur.Scalar() != 0.0)
            nnz = 1;
    }
    else if (cur.IsComplex())
    {
        if (cur.Complex() != 0.0)
            nnz = 1;
    }
    else
    {
        throw OML_Error(OML_ERR_MATRIX, 1);
    }

    outputs.push_back(nnz);
    return true;
}
//------------------------------------------------------------------------------
// Returns a sparse identity matrix [speye]
//------------------------------------------------------------------------------
bool oml_speye(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();
    int m;
    int n;

    if (nargin == 1)
    {
        if (inputs[0].IsPositiveInteger())
        {
            m = static_cast<int> (inputs[0].Scalar());
            n = m;
        }
        else if (inputs[0].IsPositiveIntegralVector())
        {
            const hwMatrix* dims = inputs[0].Matrix();

            if (dims->Size() != 2)
            {
                throw OML_Error(OML_ERR_VECTOR2, 1, OML_VAR_DIMS);
            }

            if (!dims->IsReal())
            {
                throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_DIMS);
            }

            m = static_cast<int> ((*dims)(0));
            n = static_cast<int> ((*dims)(1));
        }
        else
        {
            throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);
        }
    }
    else if (nargin == 2)
    {
        if (!inputs[0].IsPositiveInteger())
            throw OML_Error(OML_ERR_NATURALNUM, 1, OML_VAR_DIM);

        if (!inputs[1].IsPositiveInteger())
            throw OML_Error(OML_ERR_NATURALNUM, 2, OML_VAR_DIM);

        m = static_cast<int> (inputs[0].Scalar());
        n = static_cast<int> (inputs[1].Scalar());
    }
    else
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    int k = _min(m, n);
    std::vector<int> ivec(k);
    hwMatrix V(k, hwMatrix::REAL);

    for (int i = 0; i < k; ++i)
        ivec[i] = i;

    std::vector<int> jvec(ivec);
    V.SetElements(1.0);
    hwMatrixS* eye = new hwMatrixS(ivec, jvec, V, m, n);

    outputs.push_back(eye);

    return true;
}
//------------------------------------------------------------------------------
// Returns a sparse matrix with all non zero elements set to 1 [spones]
//------------------------------------------------------------------------------
bool oml_spones(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    hwMatrixS* ones = nullptr;

    if (inputs[0].IsSparse())
        ones = new hwMatrixS(*inputs[0].MatrixS());
    else if (inputs[0].IsMatrix() || inputs[0].IsScalar() || inputs[0].IsComplex())
        ones = new hwMatrixS(*inputs[0].ConvertToMatrix());
    else
        throw OML_Error(HW_ERROR_UNSUPOP);

    ones->SetElements(1.0);
    outputs.push_back(ones);

    return true;
}
//------------------------------------------------------------------------------
// Creates a dense matrix from a sparse matrix [full]
//------------------------------------------------------------------------------
bool oml_full(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs[0].IsMatrix())
    {
        outputs.push_back(inputs[0]);
        return true;
    }

    if (!inputs[0].IsSparse())
        throw OML_Error(OML_ERR_MATRIX, 1);

    const hwMatrixS* sparse = inputs[0].MatrixS();
    hwMatrix*        dense  = new hwMatrix;

    sparse->Full(*dense);

    outputs.push_back(dense);

    return true;
}
//------------------------------------------------------------------------------
// Returns transpose of input [transpose]
//------------------------------------------------------------------------------
bool oml_transpose(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency &cur = inputs[0];

    if (cur.IsScalar() || cur.IsComplex() || cur.IsMatrix() ||
        cur.IsCellArray() || cur.IsStruct() || cur.IsSparse())
    {
        if (cur.IsMatrix() || cur.IsString())
        {
            Currency ret = _transpose(eval, cur.Matrix());
            ret.SetMask(cur.GetMask());
            outputs.push_back(ret);
        }
        else if (cur.IsCellArray())
        {
            outputs.push_back(_transpose(eval, cur.CellArray()));
        }
        else if (cur.IsSparse())
        {
            const hwMatrixS* source = cur.MatrixS();
            hwMatrixS* trans = new hwMatrixS;
            trans->Transpose(*source);
            outputs.push_back(trans);
        }
        else if (cur.IsStruct())
        {
            StructData* s = cur.Struct();
            StructData* result = new StructData();
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
            outputs.push_back(ret);
        }
		else
			outputs.push_back(cur);
    }
    else
        throw OML_Error(HW_ERROR_INVINPTYPE);

    return true;
}
//------------------------------------------------------------------------------
// Returns transpose of input [transpose]
//------------------------------------------------------------------------------
bool oml_ctranspose(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error(OML_ERR_NUMARGIN);

	const Currency& cur = inputs[0];

	if (cur.IsScalar() || cur.IsComplex() || cur.IsMatrix() ||
		cur.IsCellArray() || cur.IsStruct() || cur.IsSparse())
	{
		if (cur.IsMatrix() || cur.IsString())
		{
			Currency ret = _transpose(eval, cur.Matrix());

			ret.GetWritableMatrix()->Conjugate();

			ret.SetMask(cur.GetMask());
			outputs.push_back(ret);
		}
		else if (cur.IsCellArray())
		{
			outputs.push_back(_transpose(eval, cur.CellArray()));
		}
		else if (cur.IsComplex())
		{
			outputs.push_back(cur.Complex().Conjugate());
		}
		else if (cur.IsSparse())
		{
			const hwMatrixS* source = cur.MatrixS();
			hwMatrixS* trans = new hwMatrixS;
			trans->Transpose(*source);
			trans->Conjugate();
			outputs.push_back(trans);
		}
		else if (cur.IsStruct())
		{
			StructData* s = cur.Struct();
			StructData* result = new StructData();
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
			outputs.push_back(ret);
		}
		else
			outputs.push_back(cur);
	}
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

    BuiltInFuncsUtils::CheckMathStatus(eval, result->ConvLin(*m2, *m1));
    outputs.push_back(result);
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
            throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_MTXNOTSQUARE));
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
            throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_MTXNOTSQUARE));

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
            throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_MTXNOTSQUARE));
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
                        throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
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
        hwMatrix* result = EvaluatorInterface::allocateMatrix();
        BuiltInFuncsMKL::Conj(*mtx, *result);
        outputs.push_back(result);
    }
    else if (input.IsSparse())
    {
        const hwMatrixS* mtx = input.MatrixS();
        hwMatrixS* result = new hwMatrixS;
        result->Conjugate(*mtx);
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
    {
        outputs.push_back(input1.Scalar());
    }
    else if (input1.IsComplex())
    {
        outputs.push_back(input1.Complex());
    }
    else if (input1.IsMatrix() && !input1.IsString())
    {
        const hwMatrix* mtx = input1.Matrix();

        if (mtx->Is0x0() || (mtx->IsEmpty() && mtx->IsVector()))
        {
            outputs.push_back(1.0);
        }
        else if (dim > 2)
        {
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx);
            outputs.push_back(result);
        }
        else
        {
            int m = mtx->M(), n = mtx->N();
            if (!dim)
                dim = (m > 1 ? 1 : 2);

            int start = 0;

            if ((m == 1 && n != 0 && dim == 2) || (n == 1 && m != 0 && dim == 1))
            {
                int size = mtx->Size();

                if (mtx->IsReal())
                {
                    const double* real = mtx->GetRealData();
                    double prod = 1.0;

                    for (int i = 0; i < size; ++i)
                        prod *= *(real++);

                    outputs.push_back(prod);
                }
                else
                {
                    const hwComplex* cmplx = mtx->GetComplexData();
                    hwComplex prod(1.0, 0.0);

                    for (int i = 0; i < size; ++i)
                        prod *= *(cmplx++);

                    outputs.push_back(prod);
                }
            }
            else if (dim == 1)
            {
                hwMatrix* result = EvaluatorInterface::allocateMatrix(1, n, mtx->Type());

                if (m == 0)
                {
                    result->SetElements(1.0);
                }
                else
                {
                    if (mtx->IsReal())
                    {
                        for (int j = 0; j < n; ++j)
                        {
                            double prod = 1.0;

                            for (int i = 0; i < m; ++i)
                                prod *= (*mtx)(start + i);

                            (*result)(j) = prod;
                            start += m;
                        }
                    }
                    else
                    {
                        for (int j = 0; j < n; ++j)
                        {
                            hwComplex prod(1.0, 0.0);

                            for (int i = 0; i < m; ++i)
                                prod *= mtx->z(start + i);

                            result->z(j) = prod;
                            start += m;
                        }
                    }
                }

                outputs.push_back(result);
            }
            else // dim == 2
            {
                hwMatrix* result = EvaluatorInterface::allocateMatrix(m, 1, mtx->Type());

                if (n == 0)
                {
                    result->SetElements(1.0);
                }
                else
                {
                    if (mtx->IsReal())
                    {
                        for (int i = 0; i < m; ++i)
                            (*result)(i) = (*mtx)(i);

                        for (int j = 1; j < n; ++j)
                        {
                            start += m;

                            for (int i = 0; i < m; ++i)
                                (*result)(i) *= (*mtx)(start + i);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < m; ++i)
                            result->z(i) = mtx->z(i);

                        for (int j = 1; j < n; ++j)
                        {
                            start += m;

                            for (int i = 0; i < m; ++i)
                                result->z(i) *= mtx->z(start + i);
                        }
                    }
                }

                outputs.push_back(result);
            }
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

        if (mtx->Is0x0() || (mtx->IsEmpty() && mtx->IsVector()))
        {
            outputs.push_back(0.0);
        }
        else if (dim > 2)
        {
            hwMatrix* result = EvaluatorInterface::allocateMatrix(mtx);
            outputs.push_back(result);
        }
        else
        {
            // TODO: push the summation code into hwTMatrix
            int m = mtx->M(), n = mtx->N();

            if (!dim)
                dim = (m > 1 ? 1 : 2);

            if ((m == 1 && n != 0 && dim == 2) || (n == 1 && m != 0 && dim == 1))
            {
                // vector case
                int size = mtx->Size();

                if (mtx->IsReal())
                {
                    const double* real = mtx->GetRealData();
                    double sum = 0.0;

                    for (int i = 0; i < size; ++i)
                        sum += *(real++);

                    outputs.push_back(sum);
                }
                else
                {
                    const hwComplex* cmplx = mtx->GetComplexData();
                    hwComplex sum;

                    for (int i = 0; i < size; ++i)
                        sum += *(cmplx++);

                    outputs.push_back(sum);
                }
            }
            else if (dim == 1)
            {
                hwMatrix* result = EvaluatorInterface::allocateMatrix(1, n, mtx->Type());

                if (m == 0)
                {
                    result->SetElements(0.0);
                }
                else
                {
                    if (mtx->IsReal())
                    {
                        const double* real = mtx->GetRealData();

                        for (int j = 0; j < n; ++j)
                        {
                            double sum = 0.0;

                            for (int i = 0; i < m; ++i)
                                sum += *(real++);

                            (*result)(j) = sum;
                        }
                    }
                    else
                    {
                        const hwComplex* cmplx = mtx->GetComplexData();

                        for (int j = 0; j < n; ++j)
                        {
                            hwComplex sum;

                            for (int i = 0; i < m; ++i)
                                sum += *(cmplx++);

                            result->z(j) = sum;
                        }
                    }
                }

                outputs.push_back(result);
            }
            else // dim == 2
            {
                hwMatrix* result = EvaluatorInterface::allocateMatrix(m, 1, mtx->Type());

                if (n == 0)
                {
                    result->SetElements(0.0);
                }
                else if (m > n || m > 16)
                {
                    if (mtx->IsReal())
                    {
                        const double* real = mtx->GetRealData();
                        hwMatrix dataVec(m, 1, (void*)real, hwMatrix::REAL);
                        *result = dataVec;

                        for (int j = 1; j < n; ++j)
                        {
                            real += m;
                            hwMatrix dataVec(m, 1, (void*)real, hwMatrix::REAL);
                            *result += dataVec;
                        }
                    }
                    else
                    {
                        const hwComplex* cmplx = mtx->GetComplexData();
                        hwMatrix dataVec(m, 1, (void*)cmplx, hwMatrix::COMPLEX);
                        *result = dataVec;

                        for (int j = 1; j < n; ++j)
                        {
                            cmplx += m;
                            hwMatrix dataVec(m, 1, (void*)cmplx, hwMatrix::COMPLEX);
                            *result += dataVec;
                        }
                    }
                }
                else
                {
                    if (mtx->IsReal())
                    {
                        for (int i = 0; i < m; ++i)
                        {
                            (*result)(i) = mtx->RowSumReal(i);
                        }
                    }
                    else
                    {
                        for (int i = 0; i < m; ++i)
                        {
                            result->z(i) = mtx->RowSumComplex(i);
                        }
                    }
                }

                outputs.push_back(result);
            }
        }
    }
    else if (input1.IsNDMatrix())
    {
        if (size == 1)
        {
            oml_MatrixN_sum(eval, inputs, outputs);
        }
        else
        {
            oml_MatrixN_sum(eval, inputs, outputs, 2);
        }
    }
    else if (input1.IsSparse())
    {
        const hwMatrixS* mtx = input1.MatrixS();

        if (dim == 1)
        {
            hwMatrixS* result = new hwMatrixS;
            result->Sum(*mtx, true);
            outputs.push_back(result);
        }
        else if (dim == 2)
        {
            hwMatrixS* result = new hwMatrixS;
            result->Sum(*mtx, false);
            outputs.push_back(result);
        }
        else // dim > 2
        {
            hwMatrixS* copy = new hwMatrixS(*mtx);
            outputs.push_back(copy);
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
        {
            outputs.push_back(EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL));
        }
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

            if ((m == 1 && n != 0 && dim == 2) || (n == 1 && m != 0 && dim == 1))
            {
                result = EvaluatorInterface::allocateMatrix(m, n, mtx->Type());
                int size = mtx->Size();

                if (mtx->IsReal())
                {
                    const double* real = mtx->GetRealData();
                    double* cumsum = result->GetRealData();
                    double sum = 0.0;

                    for (int i = 0; i < size; ++i)
                    {
                        sum += *(real++);
                        *(cumsum++) = sum;
                    }
                }
                else
                {
                    const hwComplex* cmplx = mtx->GetComplexData();
                    hwComplex* cumsum = result->GetComplexData();
                    hwComplex sum;

                    for (int i = 0; i < size; ++i)
                    {
                        sum += *(cmplx++);
                        *(cumsum++) = sum;
                    }
                }
            }
            else if (dim == 1)
            {
                result = EvaluatorInterface::allocateMatrix(m, n, mtx->Type());

                if (mtx->IsReal())
                {
                    const double* real = mtx->GetRealData();
                    double* cumsum = result->GetRealData();

                    for (int j = 0; j < n; ++j)
                    {
                        double sum = 0.0;

                        for (int i = 0; i < m; ++i)
                        {
                            sum += *(real++);
                            *(cumsum++) = sum;
                        }
                    }
                }
                else
                {
                    const hwComplex* cmplx = mtx->GetComplexData();
                    hwComplex* cumsum = result->GetComplexData();

                    for (int j = 0; j < n; ++j)
                    {
                        hwComplex sum;

                        for (int i = 0; i < m; ++i)
                        {
                            sum += *(cmplx++);
                            *(cumsum++) = sum;
                        }
                    }
                }
            }
            else    // dim == 2
            {
                result = EvaluatorInterface::allocateMatrix(mtx);

                if (mtx->IsReal())
                {
                    const double* realR = result->GetRealData();

                    for (int j = 1; j < n; ++j)
                    {
                        hwMatrix resultVec1(m, 1, (void*)realR, hwMatrix::REAL);
                        realR += m;
                        hwMatrix resultVec2(m, 1, (void*)realR, hwMatrix::REAL);

                        resultVec2.AddEquals(resultVec1);
                    }
                }
                else
                {
                    const hwComplex* cmplxR = result->GetComplexData();

                    for (int j = 1; j < n; ++j)
                    {
                        hwMatrix resultVec1(m, 1, (void*)cmplxR, hwMatrix::COMPLEX);
                        cmplxR += m;
                        hwMatrix resultVec2(m, 1, (void*)cmplxR, hwMatrix::COMPLEX);

                        resultVec2.AddEquals(resultVec1);
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
            oml_MatrixN_cumsum(eval, inputs, outputs);
        }
        else
        {
            oml_MatrixN_cumsum(eval, inputs, outputs, 2);
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
        {
            outputs.push_back(EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL));
        }
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

            if ((m == 1 && n != 0 && dim == 2) || (n == 1 && m != 0 && dim == 1))
            {
                int size = mtx->Size();

                if (mtx->IsReal())
                {
                    const double* real = mtx->GetRealData();
                    double* cumprod = result->GetRealData();
                    double prod = 1.0;

                    for (int i = 0; i < size; ++i)
                    {
                        prod *= *(real++);
                        *(cumprod++) = prod;
                    }
                }
                else
                {
                    const hwComplex* cmplx = mtx->GetComplexData();
                    hwComplex* cumprod = result->GetComplexData();
                    hwComplex prod(1.0, 0.0);

                    for (int i = 0; i < size; ++i)
                    {
                        prod *= *(cmplx++);
                        *(cumprod++) = prod;
                    }
                }
            }
            else if (dim == 1)
            {
                int start = 0;

                if (mtx->IsReal())
                {
                    for (int j = 0; j < n; ++j)
                    {
                        double prod = 1.0;

                        for (int i = 0; i < m; ++i)
                        {
                            prod *= (*mtx)(start + i);
                            (*result)(start + i) = prod;
                        }

                        start += m;
                    }
                }
                else
                {
                    for (int j = 0; j < n; ++j)
                    {
                        hwComplex prod(1.0, 0.0);

                        for (int i = 0; i < m; ++i)
                        {
                            prod *= mtx->z(start + i);
                            result->z(start + i) = prod;
                        }

                        start += m;
                    }
                }
            }
            else    // dim == 2
            {
                int start = 0;

                if (mtx->IsReal())
                {
                    for (int i = 0; i < m; ++i)
                        (*result)(i) = (*mtx)(i);

                    for (int j = 1; j < n; ++j)
                    {
                        start += m;

                        for (int i = 0; i < m; ++i)
                            (*result)(start + i) = (*result)(start - m + i) * (*mtx)(start + i);
                    }
                }
                else
                {
                    for (int i = 0; i < m; ++i)
                        result->z(i) = mtx->z(i);

                    for (int j = 1; j < n; ++j)
                    {
                        start += m;

                        for (int i = 0; i < m; ++i)
                            result->z(start + i) = result->z(start - m + i) * mtx->z(start + i);
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
// Returns a matrix with accumulated results [accumarray]
//------------------------------------------------------------------------------
bool oml_accumarray(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 2 || nargin > 5)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    const Currency& indexC = inputs[0];
    const Currency& valueC = inputs[1];

    if (!indexC.IsMatrix() && !indexC.IsScalar())
    {
        if (indexC.IsCellArray())
        {
            HML_CELLARRAY* cell = indexC.CellArray();
            hwMatrix* matrix = new hwMatrix;
            int length;

            for (int i = 0; i < cell->Size(); ++i)
            {
                const Currency& idx = (*cell)(i);

                if (!idx.IsMatrix() && !idx.IsScalar())
                {
                    throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_CELL);
                }

                const hwMatrix* vec = idx.ConvertToMatrix();

                if (!vec->IsReal())
                {
                    throw OML_Error(OML_ERR_REALVECTOR, 1, OML_VAR_CELL);
                }

                if (!vec->IsVector())
                {
                    throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_CELL);
                }

                if (i == 0)
                {
                    length = vec->Size();
                    hwMathStatus status = matrix->Dimension(1, length, hwMatrix::REAL);
                }
                else
                {
                    if (vec->Size() != length)
                        throw OML_Error(OML_ERR_CELLSIZE, 1, OML_VAR_DIMS);

                    hwMathStatus status = matrix->Resize(i+1, length);
                }

                for (int j = 0; j < length; ++j)
                    (*matrix)(i, j) = (*vec)(j);
            }

            std::vector<Currency> inputs2;
            inputs2.push_back(matrix);

            for (int i = 1; i < nargin; ++i)
                inputs2.push_back(inputs[i]);

            return oml_accumarray(eval, inputs2, outputs);
        }
        else
        {
            throw OML_Error(OML_ERR_MATRIX, 1, OML_VAR_TYPE);
        }
    }

    const hwMatrix* indexM = indexC.ConvertToMatrix();
    int numVals            = indexM->M();
    int numDims            = indexM->N();

    if (numDims == 0)
    {
        outputs.push_back(new hwMatrix(0, 1, hwMatrix::REAL));
        return true;
    }

    const hwMatrix* valueM = valueC.ConvertToMatrix();

    if (!indexM->IsReal())
    {
        throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_TYPE);
    }

    for (int i = 0; i < indexM->Size(); ++i)
    {
        if (!isint((*indexM)(i)))
            throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_INDEX);
    }

    if (!valueM->IsEmptyOrVector())
    {
        throw OML_Error(OML_ERR_VECTOR, 2, OML_VAR_TYPE);
    }

    if (valueM->Size() != numVals)
    {
        throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
    }

    // get dimension information
    std::vector<Currency> inputs2;
    inputs2.push_back(indexC);
    inputs2.push_back(new hwMatrix);
    inputs2.push_back(1);
    std::vector<Currency> outputs2;
    outputs2 = eval.DoMultiReturnFunctionCall(oml_max, inputs2, 3, 1, true);
    hwMatrix dimVec = (*outputs2[0].ConvertToMatrix());

    if (nargin > 2)
    {
        const Currency& dims3 = inputs[2];

        if (!dims3.IsMatrix())
        {
            throw OML_Error(OML_ERR_VECTOR, 3, OML_VAR_TYPE);
        }

        const hwMatrix* dimVec3 = dims3.Matrix();

        if (!dimVec3->IsEmptyOrVector())
        {
            throw OML_Error(OML_ERR_VECTOR, 3, OML_VAR_TYPE);
        }

        if (numDims == 1 && !dimVec3->IsEmpty())
        {
            double size = (*dimVec3)(0);

            for (int i = 1; i < dimVec3->Size(); ++i)
                size *= (*dimVec3)(i);

            if (size < dimVec(0))
            {
                throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);
            }

            dimVec = (*dimVec3);
        }
        else if (dimVec3->Size() == numDims)
        {
            for (int i = 0; i < numDims; ++i)
            {
                if ((*dimVec3)(i) > dimVec(i))
                {
                    dimVec(i) = (*dimVec3)(i);
                }
                else if ((*dimVec3)(i) < dimVec(i))
                {
                    throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);
                }
            }
        }
        else if (!dimVec3->Is0x0())
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 3);
        }
    }

    FunctionInfo* funcInfo = nullptr;
    FUNCPTR       funcPntr = nullptr;
    bool          anonFunc = false;

    if (nargin > 3)
    {
        const Currency& cur4 = inputs[3];

        if (cur4.IsFunctionHandle())
        {
            std::string funcName = cur4.FunctionHandle()->FunctionName();

            if (funcName == "anonymous")
            {
                funcInfo = cur4.FunctionHandle();
                anonFunc = true;
            }
            else if (!eval.FindFunctionByName(funcName, &funcInfo, &funcPntr, NULL))
            {
                throw OML_Error(OML_ERR_FUNCNAME, 4);
            }

            if (funcInfo && funcInfo->Parameters().size() != 1)
            {
                throw OML_Error(OML_ERR_NUMARGIN, 4);
            }
        }
        else if (cur4.IsMatrix())
        {
            const hwMatrix* empty = cur4.Matrix();

            if (!empty->Is0x0())
                throw OML_Error(OML_ERR_EMPTYMATRIX, 4, OML_VAR_TYPE);

            if (!eval.FindFunctionByName("sum", &funcInfo, &funcPntr, NULL))
            {
                throw OML_Error(OML_ERR_FUNCNAME, 4);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_HANDLE, 4, OML_VAR_TYPE);
        }
    }
    else if (!eval.FindFunctionByName("sum", &funcInfo, &funcPntr, NULL))
    {
        throw OML_Error(OML_ERR_FUNCNAME, 4);
    }

    // convert indices to memory position and sort
    hwMatrix* singleIndx = new hwMatrix(numVals, hwMatrix::REAL);

    for (int i = 0; i < numVals; ++i)
    {
        double pos = (*indexM)(i, numDims-1) - 1.0;

        for (int j = numDims-2; j > -1; --j)
            pos = pos * dimVec(j) + (*indexM)(i, j) - 1.0;

        if (pos < 0)
        {
            throw OML_Error(OML_ERR_INVALID_INDEX, 1, OML_VAR_VALUE);
        }

        (*singleIndx)(i) = pos;
    }

    inputs2.clear();
    outputs2.clear();
    inputs2.push_back(singleIndx);
    inputs2.push_back("ascend");
    oml_sort(eval, inputs2, outputs2);   // sort(singleIndx)

    hwMatrix* indxPntr = outputs2[1].GetWritableMatrix();

    (*indxPntr) -= 1.0;

    // get fill value
    bool      fillReal   = true;
    double    fillValueD = 0.0;
    hwComplex fillValueC;

    if (nargin > 4)
    {
        if (inputs[4].IsScalar())
        {
            fillValueD = inputs[4].Scalar();
        }
        else if (inputs[4].IsMatrix())
        {
            const hwMatrix* empty = inputs[4].Matrix();

            if (!empty->Is0x0())
                throw OML_Error(OML_ERR_EMPTYMATRIX, 5, OML_VAR_TYPE);
        }
        else if (inputs[4].IsComplex())
        {
            fillValueC = inputs[4].Complex();
            fillReal = false;
        }
        else
        {
            throw OML_Error(OML_ERR_SCALAR, 5, OML_VAR_TYPE);
        }
    }

    // read values and accumulate
    if (numDims < 3)
    {
        int m = static_cast<int> (dimVec(0));
        int n = 1;

        if (dimVec.Size() == 2)
            n = static_cast<int> (dimVec(1));

        // 2D case
        hwMatrix* result = new hwMatrix(m, n, hwMatrix::REAL);
        hwMatrixI set(m, n, hwMatrixI::REAL);
        set.SetElements(0);

        if (valueM->IsReal())
        {
            for (int i = 0; i < numVals; ++i)
            {
                int       p1    = static_cast<int> ((*indxPntr)(i));
                int       count = 1;
                hwMatrix* vec   = new hwMatrix(1, hwMatrix::REAL);
                (*vec)(0)       = (*valueM)(p1);

                int j = i + 1;

                while (j < numVals)
                {
                    int p2 = static_cast<int> ((*indxPntr)(j));

                    if ((*singleIndx)(p1) != (*singleIndx)(p2))
                    {
                        break;
                    }

                    vec->Resize(j-i+1);
                    (*vec)(j-i) = (*valueM)(p2);
                    ++j;
                }

                std::vector<Currency> inputs3;
                std::vector<Currency> outputs3;
                inputs3.push_back(vec);

                try
                {
                    if (anonFunc)
                    {
                        Currency result = eval.CallInternalFunction(funcInfo, inputs3);
                        outputs3.push_back(result);
                    }
                    else if (funcInfo)
                    {
                        outputs3 = eval.DoMultiReturnFunctionCall(funcInfo, inputs3, 1, 1, true);
                    }
                    else if (funcPntr)
                    {
                        outputs3 = eval.DoMultiReturnFunctionCall(funcPntr, inputs3, 1, 1, true);
                    }
                }
                catch (OML_Error&)
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                if (outputs3.size() < 1)
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                Currency& out3 = outputs3[0];
                int pos = static_cast<int>((*singleIndx)(static_cast<int>(p1)));

                if (out3.IsScalar())
                {
                    if (result->IsReal())
                        (*result)(pos) = out3.Scalar();
                    else
                        result->z(pos) = out3.Scalar();
                }
                else if (out3.IsComplex())
                {
                    if (result->IsReal())
                        result->MakeComplex();

                    result->z(pos) = out3.Complex();
                }
                else
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                set(pos) = 1;
                i += j - i - 1;
            }

            if (fillReal)
            {
                for (int i = 0; i < set.Size(); ++i)
                {
                    if (set(i) == 0)
                        (*result)(i) = fillValueD;
                }
            }
            else
            {
                result->MakeComplex();

                for (int i = 0; i < set.Size(); ++i)
                {
                    if (set(i) == 0)
                        result->z(i) = fillValueC;
                }
            }
        }
        else
        {
            for (int i = 0; i < numVals; ++i)
            {
                int       p1    = static_cast<int> ((*indxPntr)(i));
                int       count = 1;
                hwMatrix* vec   = new hwMatrix(1, hwMatrix::COMPLEX);
                vec->z(0)       = valueM->z(p1);

                int j = i + 1;

                while (j < numVals)
                {
                    int p2 = static_cast<int> ((*indxPntr)(j));

                    if ((*singleIndx)(p1) != (*singleIndx)(p2))
                    {
                        break;
                    }

                    vec->Resize(j-i+1);
                    vec->z(j-i) = valueM->z(p2);
                    ++j;
                }

                std::vector<Currency> inputs3;
                std::vector<Currency> outputs3;
                inputs3.push_back(vec);

                try
                {
                    if (anonFunc)
                    {
                        Currency result = eval.CallInternalFunction(funcInfo, inputs3);
                        outputs3.push_back(result);
                    }
                    else if (funcInfo)
                    {
                        outputs3 = eval.DoMultiReturnFunctionCall(funcInfo, inputs3, 1, 1, true);
                    }
                    else if (funcPntr)
                    {
                        outputs3 = eval.DoMultiReturnFunctionCall(funcPntr, inputs3, 1, 1, true);
                    }
                }
                catch (OML_Error&)
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                if (outputs3.size() < 1)
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                Currency& out3 = outputs3[0];
                int pos = static_cast<int>((*singleIndx)(static_cast<int>(p1)));

                if (out3.IsScalar())
                {
                    if (result->IsReal())
                        (*result)(pos) = out3.Scalar();
                    else
                        result->z(pos) = out3.Scalar();
                }
                else if (out3.IsComplex())
                {
                    if (result->IsReal())
                        result->MakeComplex();

                    result->z(pos) = out3.Complex();
                }
                else
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                set(pos) = 1;
                i += j - i - 1;
            }

            if (fillReal)
            {
                for (int i = 0; i < set.Size(); ++i)
                {
                    if (set(i) == 0)
                        result->z(i) = fillValueD;
                }
            }
            else
            {
                result->MakeComplex();

                for (int i = 0; i < set.Size(); ++i)
                {
                    if (set(i) == 0)
                        result->z(i) = fillValueC;
                }
            }
        }

        outputs.push_back(result);
    }
    else
    {
        // ND case
        std::vector<int> dimsI;

        for (int i = 0; i < numDims; ++i)
            dimsI.push_back(static_cast<int> (dimVec(i)));

        hwMatrixN* result = new hwMatrixN(dimsI, hwMatrixN::REAL);
        hwMatrixNI set(dimsI, hwMatrixNI::REAL);
        set.SetElements(0);

        if (valueM->IsReal())
        {
            for (int i = 0; i < numVals; ++i)
            {
                int       p1    = static_cast<int> ((*indxPntr)(i));
                int       count = 1;
                hwMatrix* vec   = new hwMatrix(1, hwMatrix::REAL);
                (*vec)(0)       = (*valueM)(p1);

                int j = i + 1;

                while (j < numVals)
                {
                    int p2 = static_cast<int> ((*indxPntr)(j));

                    if ((*singleIndx)(p1) != (*singleIndx)(p2))
                    {
                        break;
                    }

                    vec->Resize(j-i+1);
                    (*vec)(j-i) = (*valueM)(p2);
                    ++j;
                }

                std::vector<Currency> inputs3;
                std::vector<Currency> outputs3;
                inputs3.push_back(vec);

                try
                {
                    if (anonFunc)
                    {
                        Currency result = eval.CallInternalFunction(funcInfo, inputs3);
                        outputs3.push_back(result);
                    }
                    else if (funcInfo)
                    {
                        outputs3 = eval.DoMultiReturnFunctionCall(funcInfo, inputs3, 1, 1, true);
                    }
                    else if (funcPntr)
                    {
                        outputs3 = eval.DoMultiReturnFunctionCall(funcPntr, inputs3, 1, 1, true);
                    }
                }
                catch (OML_Error&)
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                if (outputs3.size() < 1)
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                Currency& out3 = outputs3[0];
                int pos = static_cast<int>((*singleIndx)(static_cast<int>(p1)));

                if (out3.IsScalar())
                {
                    if (result->IsReal())
                        (*result)(pos) = out3.Scalar();
                    else
                        result->z(pos) = out3.Scalar();
                }
                else if (out3.IsComplex())
                {
                    if (result->IsReal())
                        result->MakeComplex();

                    result->z(pos) = out3.Complex();
                }
                else
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                set(pos) = 1;
                i += j - i - 1;
            }

            if (fillReal)
            {
                for (int i = 0; i < set.Size(); ++i)
                {
                    if (set(i) == 0)
                        (*result)(i) = fillValueD;
                }
            }
            else
            {
                result->MakeComplex();

                for (int i = 0; i < set.Size(); ++i)
                {
                    if (set(i) == 0)
                        result->z(i) = fillValueC;
                }
            }
        }
        else
        {
            for (int i = 0; i < numVals; ++i)
            {
                int       p1    = static_cast<int> ((*indxPntr)(i));
                int       count = 1;
                hwMatrix* vec   = new hwMatrix(1, hwMatrix::COMPLEX);
                vec->z(0)       = valueM->z(p1);

                int j = i + 1;

                while (j < numVals)
                {
                    int p2 = static_cast<int> ((*indxPntr)(j));

                    if ((*singleIndx)(p1) != (*singleIndx)(p2))
                    {
                        break;
                    }

                    vec->Resize(j-i+1);
                    vec->z(j-i) = valueM->z(p2);
                    ++j;
                }

                std::vector<Currency> inputs3;
                std::vector<Currency> outputs3;
                inputs3.push_back(vec);

                try
                {
                    if (anonFunc)
                    {
                        Currency result = eval.CallInternalFunction(funcInfo, inputs3);
                        outputs3.push_back(result);
                    }
                    else if (funcInfo)
                    {
                        outputs3 = eval.DoMultiReturnFunctionCall(funcInfo, inputs3, 1, 1, true);
                    }
                    else if (funcPntr)
                    {
                        outputs3 = eval.DoMultiReturnFunctionCall(funcPntr, inputs3, 1, 1, true);
                    }
                }
                catch (OML_Error&)
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                if (outputs3.size() < 1)
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                Currency& out3 = outputs3[0];
                int pos = static_cast<int>((*singleIndx)(static_cast<int>(p1)));

                if (out3.IsScalar())
                {
                    if (result->IsReal())
                        (*result)(pos) = out3.Scalar();
                    else
                        result->z(pos) = out3.Scalar();
                }
                else if (out3.IsComplex())
                {
                    if (result->IsReal())
                        result->MakeComplex();

                    result->z(pos) = out3.Complex();
                }
                else
                {
                    throw OML_Error(OML_ERR_ACCUMFUNC, 4, OML_VAR_FUNC);
                }

                set(pos) = 1;
                i += j - i - 1;
            }

            if (fillReal)
            {
                for (int i = 0; i < set.Size(); ++i)
                {
                    if (set(i) == 0)
                        result->z(i) = fillValueD;
                }
            }
            else
            {
                result->MakeComplex();

                for (int i = 0; i < set.Size(); ++i)
                {
                    if (set(i) == 0)
                        result->z(i) = fillValueC;
                }
            }
        }

        outputs.push_back(result);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes the base-2 logarithm of input [log2]
//------------------------------------------------------------------------------
bool oml_log2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    // This function gets called from BuiltInFuncsMKL::Log2 for the two-argument return case.
    // Prior to MKL it supported all of the cases.
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
// Computes the next higher power of 2 of input [nextpow2]
//------------------------------------------------------------------------------
bool oml_nextpow2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    bool retv;
    std::vector<Currency> outputs2;
    std::vector<Currency> outputs3;

    if (!BuiltInFuncsMKL::Abs(eval, inputs, outputs2))
        return false;

    if (!oml_log2(eval, outputs2, outputs3))
        return false;

    outputs3 = eval.DoMultiReturnFunctionCall(oml_log2, outputs2, 1, 2, true);

    if (outputs3[0].IsScalar())
    {
        double f = outputs3[0].Scalar();
        double n = outputs3[1].Scalar();

        // if (IsNanT(f))
        if (isnan(f))
        {
            n = std::numeric_limits<double>::quiet_NaN();
        }
        else if (IsInf_T(f) || IsNegInf_T(f))
        {
            n = std::numeric_limits<double>::infinity();
        }
        else if (f == 0.5)
        {
            n -= 1.0;
        }

        outputs.push_back(n);
    }
    else if (outputs3[0].IsMatrix())
    {
        hwMatrix* F = outputs3[0].GetWritableMatrix();
        hwMatrix* N = outputs3[1].GetWritableMatrix();
        int size = F->Size();

        for (int i = 0; i < size; ++i)
        {
            double f = (*F)(i);

            if (isnan(f))
            {
                (*N)(i) = std::numeric_limits<double>::quiet_NaN();
            }
            else if (IsInf_T(f) || IsNegInf_T(f))
            {
                (*N)(i) = std::numeric_limits<double>::infinity();
            }
            else if (f == 0.5)
            {
                (*N)(i) -= 1.0;
            }
        }

        N->IncrRefCount();
        outputs.push_back(N);
    }
    else // ND case
    {
        hwMatrixN* F = outputs3[0].GetWritableMatrixN();
        hwMatrixN* N = outputs3[1].GetWritableMatrixN();
        int size = F->Size();

        for (int i = 0; i < size; ++i)
        {
            double f = (*F)(i);

            if (isnan(f))
            {
                (*N)(i) = std::numeric_limits<double>::quiet_NaN();
            }
            else if (IsInf_T(f) || IsNegInf_T(f))
            {
                (*N)(i) = std::numeric_limits<double>::infinity();
            }
            else if (f == 0.5)
            {
                (*N)(i) -= 1.0;
            }
        }

        N->IncrRefCount();
        outputs.push_back(N);
    }

    return true;
}
//------------------------------------------------------------------------------
// Computes the next higher power of 2 of input [round]
//------------------------------------------------------------------------------
bool oml_round(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (nargin == 1)
    {
        BuiltInFuncsMKL::Round(eval, inputs, outputs);
    }
    else
    {
        if (!inputs[1].IsInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2);

        std::string type;

        if (nargin == 3)
        {
            if (!inputs[2].IsString())
                throw OML_Error(OML_ERR_STRING, 3);

            type = inputs[2].StringVal();
        }
        else
        {
            type = "decimal";
        }

        if (inputs[0].IsScalar())
        {
            double value = inputs[0].Scalar();

            if (type == "decimal")
            {
                double shift = pow(10.0, inputs[1].Scalar());
                outputs.push_back(round(value * shift) / shift);
            }
            else if (type == "significant")
            {
                if (!inputs[1].IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 2);

                double ndright = inputs[1].Scalar() - (static_cast<int>(log10(abs(value))) + 1);
                double shift = pow(10.0, ndright);
                outputs.push_back(round(value * shift) / shift);
            }
            else if (type == "binary")
            {
                double shift = pow(2.0, inputs[1].Scalar());
                outputs.push_back(round(value * shift) / shift);
            }
            else if (type == "sigbits")
            {
                if (!inputs[1].IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 2);

                double ndright = inputs[1].Scalar() - (static_cast<int>(log2(abs(value))) + 1);
                double shift = pow(2.0, ndright);
                outputs.push_back(round(value * shift) / shift);
            }
            else
            {
                throw OML_Error(OML_ERR_BAD_STRING, 3);
            }
        }
        else if (inputs[0].IsComplex())
        {
            hwComplex value = inputs[0].Complex();
            double real = value.Real();
            double imag = value.Imag();
            std::vector<Currency> inputs2;
            std::vector<Currency> outputs2;

            inputs2.push_back(real);
            inputs2.push_back(inputs[1]);
            inputs2.push_back(inputs[2]);
            oml_round(eval, inputs2, outputs2);
            real = outputs2[0].Scalar();
            outputs2.clear();

            inputs2[0] = imag;
            oml_round(eval, inputs2, outputs2);
            imag = outputs2[0].Scalar();

            outputs.push_back(hwComplex(real, imag));
        }
        else if (inputs[0].IsMatrix())
        {
            const hwMatrix* mtx = inputs[0].ConvertToMatrix();

            if (type == "decimal")
            {
                double shift = pow(10.0, inputs[1].Scalar());
                hwMatrix* temp = EvaluatorInterface::allocateMatrix();
                (*temp) = (*mtx) * shift;
                std::vector<Currency> inputs2;
                inputs2.push_back(temp);
                BuiltInFuncsMKL::Round(eval, inputs2, outputs);
                (*outputs[0].GetWritableMatrix()) /= shift;
            }
            else if (type == "significant" || type == "sigbits")
            {
                if (!inputs[1].IsPositiveInteger())
                    throw OML_Error(OML_ERR_POSINTEGER, 2);

                if (mtx->IsReal())
                {
                    std::vector<Currency> inputs2;
                    inputs2.push_back(inputs[0]);
                    std::vector<Currency> absc;
                    BuiltInFuncsMKL::Abs(eval, inputs2, absc);
                    inputs2.clear();

                    inputs2.push_back(absc[0]);

                    if (type == "significant")
                    {
                        std::vector<Currency> log10c;
                        BuiltInFuncsMKL::Log10(eval, inputs2, log10c);
                        inputs2.clear();
                        inputs2.push_back(log10c[0]);
                    }
                    else
                    {
                        std::vector<Currency> log2c;
                        BuiltInFuncsMKL::Log2(eval, inputs2, log2c);
                        inputs2.clear();
                        inputs2.push_back(log2c[0]);
                    }

                    std::vector<Currency> ndright;
                    BuiltInFuncsMKL::Fix(eval, inputs2, ndright);
                    (*ndright[0].GetWritableMatrix()) += 1.0;
                    (*ndright[0].GetWritableMatrix()) *= -1.0;
                    (*ndright[0].GetWritableMatrix()) += inputs[1].Scalar();
                    hwMatrix shift;

                    if (type == "significant")
                    {
                        hwMatrix tens(mtx->M(), mtx->N(), hwMatrix::REAL);
                        tens.SetElements(10.0);
                        BuiltInFuncsMKL::PowerByElems(tens, (*ndright[0].GetWritableMatrix()), shift);
                    }
                    else
                    {
                        hwMatrix twos(mtx->M(), mtx->N(), hwMatrix::REAL);
                        twos.SetElements(10.0);
                        BuiltInFuncsMKL::PowerByElems(twos, (*ndright[0].GetWritableMatrix()), shift);
                    }

                    inputs2.clear();
                    hwMatrix* temp = EvaluatorInterface::allocateMatrix();
                    BuiltInFuncsMKL::MultByElems(*mtx, shift, *temp);
                    inputs2.push_back(temp);
                    std::vector<Currency> rounded;
                    BuiltInFuncsMKL::Round(eval, inputs2, rounded);
                    hwMatrix* result = EvaluatorInterface::allocateMatrix();
                    BuiltInFuncsMKL::DivideByElems((*rounded[0].GetWritableMatrix()), shift, *result);

                    outputs.push_back(result);
                }
                else    // complex
                {
                    hwMatrix* Real = EvaluatorInterface::allocateMatrix();
                    hwMatrix* Imag = EvaluatorInterface::allocateMatrix();
                    mtx->UnpackComplex(Real, Imag);

                    std::vector<Currency> inputs2;
                    std::vector<Currency> outputs2;
                    inputs2.push_back(Real);
                    inputs2.push_back(inputs[1]);
                    inputs2.push_back(inputs[2]);
                    oml_round(eval, inputs2, outputs2);
                    hwMatrix* RealR = outputs2[0].GetWritableMatrix();

                    std::vector<Currency> outputs3;
                    inputs2.clear();
                    inputs2.push_back(Imag);
                    inputs2.push_back(inputs[1]);
                    inputs2.push_back(inputs[2]);
                    oml_round(eval, inputs2, outputs3);
                    hwMatrix* ImagR = outputs3[0].GetWritableMatrix();

                    hwMatrix* cmplx = EvaluatorInterface::allocateMatrix();
                    cmplx->PackComplex(*RealR, ImagR);
                    outputs.push_back(cmplx);
                }
            }
            else if (type == "binary")
            {
                double shift = pow(2.0, inputs[1].Scalar());
                hwMatrix* temp = EvaluatorInterface::allocateMatrix();
                (*temp) = (*mtx) * shift;
                std::vector<Currency> inputs2;
                inputs2.push_back(temp);
                BuiltInFuncsMKL::Round(eval, inputs2, outputs);
                (*outputs[0].GetWritableMatrix()) /= shift;
            }
            else
            {
                throw OML_Error(OML_ERR_BAD_STRING, 3);
            }
        }
        else if (inputs[0].IsNDMatrix())
        {
            oml_MatrixNUtil1(eval, inputs, outputs, oml_round);
        }
        else
        {
            throw OML_Error(OML_ERR_MATRIX, 1);
        }
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
        VARIABLES,
        CLASSES
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
        else if (str == "-c" || str == "-classes" || str == "classes")
            toclear = CLASSES;
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
            case CLASSES:
                eval.ClearClassesAndObjects();
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
		else if (inputs[0].IsFunctionHandle())
		{
			FunctionInfo* fi = inputs[0].FunctionHandle();
			count = fi->GetRefCount();

			count--;
		}
    }

    outputs.push_back(count);
    return true;
}
//------------------------------------------------------------------------------
// Get/set last error message [lasterr]
//------------------------------------------------------------------------------
#undef FormatMessage

bool oml_lasterr(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() > 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (inputs.size() == 0)
    {
        int nargout = eval.GetNargoutValue();

        if (nargout > 1)
            throw OML_Error(OML_ERR_NUMARGOUT);

        if (eval.GetVerbose())
            outputs.push_back(eval.FormatMessage(EvaluatorInterface::GetLastErrorMessage()));
        else
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
            const Currency& x_cur = inputs[0];
            const Currency& y_cur = inputs[1];

            x = x_cur.ConvertToMatrix();
            y = y_cur.ConvertToMatrix();
        }
        else if (nargin == 1)
        {
            const Currency& x_cur = inputs[0];

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
            x = inputs[0].ConvertToMatrix();
            y = x;

            if (nargout == 3)
                z = x;
        }
        else if (nargin == 2)
        {
            x = inputs[0].ConvertToMatrix();
            y = inputs[1].ConvertToMatrix();
        }
        else
        {
            x = inputs[0].ConvertToMatrix();
            y = inputs[1].ConvertToMatrix();
            z = inputs[2].ConvertToMatrix();
        }

        dims[0] = y->Size();
        dims[1] = x->Size();

        if (nargout == 3)
            dims[2] = z->Size();

        // process first output
        hwMatrixN* vec1 = EvaluatorInterface::allocateMatrixN();
        vec1->Convert2DtoND(*x, false);

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
        vec2->Convert2DtoND(*y, false);
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
            vec->Convert2DtoND(*z, false);

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
            const Currency& x_cur = inputs[0];
            const Currency& y_cur = inputs[1];

            x = x_cur.ConvertToMatrix();
            y = y_cur.ConvertToMatrix();
        }
        else if (nargin == 1)
        {
            const Currency& x_cur = inputs[0];

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
        vec->Convert2DtoND(*inputs[0].Matrix(), false);

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
                vec2->Convert2DtoND(*inputs[0].Matrix(), false);
            else
                vec2->Convert2DtoND(*inputs[i].Matrix(), false);

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
    else if (input.IsSparse())
    {
        const hwMatrixS* mtx = input.MatrixS();
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
	else if (input.IsStruct() || input.IsObject())
    {
        StructData* sd = input.Struct();
        outputs.push_back(sd->M()*sd->N());
    }
	else if (input.IsPointer())
	{
		Currency ref = *input.Pointer();

		if (ref.IsObject())
		{
			StructData* sd = ref.Struct();
			outputs.push_back(sd->M() * sd->N());
		}
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

    if (target.IsMatrixOrString() || target.IsScalar() || target.IsComplex())
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
        else
        {
            const hwMatrix* mat = target.ConvertToMatrix();

            if (mat->M() == 1 && mat->N() > 1)
                dim_val = 2;
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

        const hwMatrix* target_mtx = target.ConvertToMatrix();
        int m = target_mtx->M();
        int n = target_mtx->N();
        hwMatrix*       result     = EvaluatorInterface::allocateMatrix(m, n, target_mtx->Type());

        if (target_mtx->IsReal())
        {
            for (int j=0; j<m; j++)
            {
                int shift_index_1 = (j+shift_val_1) % target_mtx->M();

                if (shift_index_1 < 0)
                    shift_index_1 += m;

                for (int k=0; k<n; k++)
                {
                    int shift_index_2 = (k+shift_val_2) % target_mtx->N();

                    if (shift_index_2 < 0)
                        shift_index_2 += n;

                    (*result)(shift_index_1, shift_index_2) = (*target_mtx)(j,k); 
                }
            }
        }
        else
        {
            for (int j=0; j<m; j++)
            {
                int shift_index_1 = (j+shift_val_1) % target_mtx->M();

                if (shift_index_1 < 0)
                    shift_index_1 += m;

                for (int k=0; k<n; k++)
                {
                    int shift_index_2 = (k+shift_val_2) % target_mtx->N();

                    if (shift_index_2 < 0)
                        shift_index_2 += n;

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
// Shifts dimensions of the input [shiftdim]
//------------------------------------------------------------------------------
bool oml_shiftdim(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    int nargin = static_cast<int> (inputs.size());

    if (nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (nargin == 0)
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
        return true;
    }

    std::vector<Currency> inputs2;
    inputs2.push_back(inputs[0]);
    std::vector<Currency> outputs2;
    outputs2 = eval.DoMultiReturnFunctionCall(oml_size, inputs2, 1, 1, true);
    const Currency& dims = outputs2[0];
    const hwMatrix* dimVec = dims.ConvertToMatrix();
    int n = dimVec->Size();
    int k = 0;

    if (nargin == 1)
    {
        for (int i = 0; i < n; ++i)
        {
            if ((*dimVec)(i) != 1.0)
            {
                k = i;
                break;
            }
        }
    }
    else
    {
        if (!inputs[1].IsInteger())
            throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_TYPE);

        k = static_cast<int>(inputs[1].Scalar());
    }

    if (k == 0)
    {
        outputs.push_back(inputs[0]);
    }
    else if (k > 0)
    {
        inputs2.clear();

        hwMatrix* permVec = new hwMatrix(n, hwMatrix::REAL);

        for (int i = 0; i < n; ++i)
            (*permVec)(i) = static_cast<double>(i+1);

        inputs2.push_back(permVec);
        inputs2.push_back(-k);
        outputs2.clear();
        oml_circshift(eval, inputs2, outputs2);

        inputs2.clear();
        inputs2.push_back(inputs[0]);
        inputs2.push_back(outputs2[0]);
        oml_permute(eval, inputs2, outputs);
    }
    else
    {
        inputs2.clear();
        inputs2.push_back(1.0);
        inputs2.push_back(-k);
        std::vector<Currency> outputs3;
        oml_ones(eval, inputs2, outputs3);

        inputs2.clear();
        inputs2.push_back(outputs3[0]);
        inputs2.push_back(dims);
        outputs2.clear();
        oml_horzcat(eval, inputs2, outputs2);

        inputs2.clear();
        inputs2.push_back(inputs[0]);
        inputs2.push_back(outputs2[0]);
        oml_reshape(eval, inputs2, outputs);
    }

    if (nargin == 1)
        outputs.push_back(k);

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
	ALT_FUNCPTR   aptr = NULL;

	if (input.IsString())
		eval.FindFunctionByName(input.StringVal(), &fi, &fptr, &aptr);
	else if (input.IsFunctionHandle())
		fi = input.FunctionHandle();

	if (fi)
	{
		if (fi->IsBuiltIn())
			outputs.push_back("builtin function");
		else if (fi->IsEncrypted())
			outputs.push_back("encrypted function");
		else
			outputs.push_back(fi->GetAST());
	}
	else if (fptr || aptr)
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

			eval.FindFunctionByName(target, &fi, NULL, NULL);

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

bool oml_getmetadata(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error(OML_ERR_NUMARGIN);

	if (inputs[0].IsString())
	{
		Currency ret = eval.GetMetadata(inputs[0].StringVal());
		outputs.push_back(ret);
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
	if (inputs.empty())
	{
		throw OML_Error(OML_ERR_NUMARGIN);
	}
	if (!inputs[0].IsString())
	{
		throw OML_Error(OML_ERR_STRING, 1);
	}

	FunctionInfo* fi = nullptr;
	FUNCPTR       fp = nullptr;
	eval.FindFunctionByName(inputs[0].StringVal(), &fi, &fp, NULL);

	const char* ROOT_DIRECTORY = "OML_HELP";
	char *c = getenv(ROOT_DIRECTORY);

	std::string userDocLocation;

	if (c)
	{
		Interpreter intBase(eval);
		std::string base = c;
		std::string helpdir;
		helpdir = intBase.GetHelpDirectory(inputs[0].StringVal());

#ifdef _DEBUG
#    define DEBUG_PRINT(s) { std::cout << s << std::endl; }
#else
#    define DEBUG_PRINT(s)
#endif

		// DIRECTORY_DELIM
		//if(helpdir.find_last_of("/\\") !=  string::npos)
		if (helpdir.find_last_of(DIRECTORY_DELIM) != string::npos)
		{
			userDocLocation = helpdir + DIRECTORY_DELIM + inputs[0].StringVal() + ".htm";
			if (!BuiltInFuncsUtils::FileExists(userDocLocation))
			{
				DEBUG_PRINT("DEBUG1: FILE NOT FOUND")
				DEBUG_PRINT(userDocLocation)
				std::string lower(inputs[0].StringVal());
				std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

				userDocLocation = helpdir + DIRECTORY_DELIM + lower + ".htm";
				DEBUG_PRINT("DEBUG1: CHECK FOR FILE")
				DEBUG_PRINT(userDocLocation)

			}

		}
		else
		{
			userDocLocation = base + intBase.GetHelpModule(inputs[0].StringVal()) + DIRECTORY_DELIM + inputs[0].StringVal() + ".htm";
			if (!BuiltInFuncsUtils::FileExists(userDocLocation))
			{
				DEBUG_PRINT("DEBUG2: FILE NOT FOUND")
				DEBUG_PRINT(userDocLocation)
				std::string lower(inputs[0].StringVal());
				std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

				userDocLocation = base + intBase.GetHelpModule(inputs[0].StringVal()) + DIRECTORY_DELIM + lower + ".htm";
				DEBUG_PRINT("DEBUG2: CHECK FOR FILE")
				DEBUG_PRINT(userDocLocation)
			}
		}

		// userDocLocation = base + intBase.GetHelpModule(inputs[0].StringVal()) + "/" + inputs[0].StringVal() + ".htm";
	}

	std::stringstream ss;
	std::string hwRootDir;

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
		if (fi)
		{
			std::string help_str = fi->HelpString();
			if (!help_str.empty())
			{
				eval.PrintResult(help_str);
				return true;
			}
		}

		std::string msg = "Could not locate help document for '" + inputs[0].StringVal() + "'";
#if defined(_DEBUG)
		msg.append(" ");
		msg.append(userDocLocation);
#endif

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
    else if (input.IsMatrix())
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

        if (!(checkisfinite(m) && checkisfinite(n)))
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

                if (!(checkisfinite(dim)))
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

                if (!(checkisfinite(dim)))
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

            if (!(checkisfinite(dim)))
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
            throw OML_Error(OML_ERR_STRING_STRINGCELL);

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
	else if (in1.IsSparse())
	{
		return false;   // sparse case is handled upon return
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
hwMatrix* safeMatrixCopyFromInput(const Currency& input, bool allowString)
{
    hwMatrix* mtx = nullptr;

    if (input.IsScalar())
        mtx = EvaluatorInterface::allocateMatrix(1, 1, input.Scalar());
    else if (input.IsComplex())
        mtx = EvaluatorInterface::allocateMatrix(1, 1, input.Complex());
    else if (input.IsMatrix() || (allowString && input.IsString()))
        mtx = EvaluatorInterface::allocateMatrix(input.Matrix());
    else
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
            throw OML_Error(OML_ERR_POSINTEGER);
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
double rem(double a, double b)
{
    return a - b * (long long) (a / b);
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
		double d1;
		double d2;

		int size = mtx1->Size();

		if (mtx2->Size() > size)
			size = mtx2->Size();

        for (int i = 0; i < size; i++)
        {
			if (i < mtx1->Size())
				d1 = (*mtx1)(i);
			else
				return true;
			
			if (i < mtx2->Size())
				d2 = (*mtx2)(i);
			else
				return false;

            if (d1 == d2)
                continue;
            if (d1 < d2)
                return true;
            return false;
        }
        return false;
    }
}

// string-related methods
//------------------------------------------------------------------------------
// Returns string with escape characters
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
// Returns appropriate escaped character
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
            return '\v';  // Set a warning later as vertical tabs are not supported
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
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
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
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
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
			else if (toappend.IsMatrixOrString() || toappend.IsScalar())
			{
				(*outcell)(j) = dostrcat(eval, strout, toappend.ConvertToMatrix());
			}
			else
			{
				throw OML_Error("Invalid concatentation");
			}
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
    BuiltInFuncsUtils::CheckFileIndex(eval, fid, 1, true);
    std::FILE *file = eval.GetFile(fid);

    bool flushcout = BuiltInFuncsUtils::IsFlushCout(eval, fid);

    if (nargin > 1)  // Reads specified number of characters
    {
        if (!inputs[1].IsInteger())
        {
            throw OML_Error(OML_ERR_NATURALNUM, 2);
        }
        int len = static_cast<int>(inputs[1].Scalar());
        if (len < 0)
        {
            throw OML_Error(OML_ERR_NATURALNUM, 2);
        }

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
            if (flushcout)
            {
                std::cout << std::flush;
            }
            return false;
        }
        readline = std::string(data);
        if (flushcout)
        {
            std::cout << std::flush;
        }
        delete [] data;
        return true;
    }
   
    // Read file till the first '\r' or '\n' or '\r\n' is encountered
    while (1)
    {
        int c = fgetc (file);
        if (flushcout)
        {
            std::cout << std::flush;
        }

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
            if (flushcout)
            {
                std::cout << std::flush;
            }

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
    else if (input.IsSparse())
    {
        return input.MatrixS()->IsEmpty();
    }
    else if (input.IsCellArray())
    {
        return input.CellArray()->IsEmpty();
    }
	else if (input.IsNDMatrix())
	{
		return input.MatrixN()->IsEmpty();
	}
	else if (input.IsNDCellArray())
	{
		return input.CellArrayND()->IsEmpty();
	}
    else if (input.IsStruct() || input.IsObject())
    {
        StructData *sd = input.Struct();
        return !(sd->M() * sd->N());
    }
	else if (input.IsFunctionHandle())
	{
		return false;
	}
    else
		throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINPUT));
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
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);

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
    else if (cur.IsSparse())
        return !cur.MatrixS()->IsRealData();
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
    else if (cur.IsSparse())
        return cur.MatrixS()->IsRealData();
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
// Returns true if given path is a directory
//------------------------------------------------------------------------------
bool isDirectory(std::string& str, std::string* errmsg)
{
#ifndef OS_WIN
    if (str == "/")  // Handle root drive on Linux
    {
        return true;
    }
#endif

    BuiltInFuncsUtils::StripTrailingSlash(str);

#if OS_WIN
    // handle root drives
    if (str.length() == 2 && isalpha(str[0]) && str[1] == ':')
        str += '\\';

    std::wstring wstr(BuiltInFuncsUtils::StdString2WString(str));
    struct _stat64i32 st;
    if (_wstat(wstr.c_str(), &st) == 0)
    {
        if (st.st_mode & S_IFDIR)
            return true;  // This is a valid directory

        else if (errmsg)
            *errmsg = str + " is not a directory";
        return false;
    }
#else
    struct stat st;
    if (stat(str.c_str(), &st) == 0)
    {
        if (S_ISDIR(st.st_mode))
            return true;
        else if (errmsg)
            *errmsg = str + " is not a directory";
        return false;
    }
#endif
    if (errmsg)
    {
        if (errno == ENOENT || errno == ENOTDIR)
            *errmsg = "no such directory exists";
        else
            *errmsg = "problem getting file info";
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
            try
            {
                if (!std::regex_match(seg.cbegin(), seg.cend(), pattern))
                    throw OML_Error(HW_ERROR_INVALIDFORMAT);
            }
            catch (std::regex_error& err)
            {
                // Ignore memory stack error from regex
                if (err.code() != std::regex_constants::error_stack)
                {
                    BuiltInFuncsUtils utils;
                    utils.ThrowRegexError(err.code());
                }
            }
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

    std::string result;
    int inputindex = 0;
    FormatSegment fs = segments[0];
    size_t currentSegment = 0;
    std::vector<Currency>::const_iterator iter    = rawiter;
    std::vector<Currency>::const_iterator enditer = rawenditer;

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
            if (iter == enditer || (iter == rawiter && !inputindex))
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

    if (!(checkisfinite(dbm) && checkisfinite(dbn)))
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
		throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINPUT));
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
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINDEX));
        return cur;
    }
    else if (cur.IsMatrix() || cur.IsString())
    {
        const hwMatrix *mtx = cur.Matrix();
        if (i >= mtx->M() || j >= mtx->N())
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINDEX));

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
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINDEX));
        HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(1, 1);
        (*newcell)(0) = (*cell)(i, j);
        return newcell;
    }
    else if (cur.IsStruct())
    {
        StructData *sd = cur.Struct();
        if (i >= sd->M() || j >= sd->N())
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINDEX));
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
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINDEX));
        return cur;
    }
    else if (cur.IsMatrix() || cur.IsString())
    {
        const hwMatrix *mtx = cur.Matrix();
        if (index >= mtx->Size())
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINDEX));

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
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINDEX));
        HML_CELLARRAY *newcell = EvaluatorInterface::allocateCellArray(1, 1);
        (*newcell)(0) = (*cell)(index);
        return newcell;
    }
    else if (cur.IsStruct())
    {
        StructData *sd = cur.Struct();
        if (index >= sd->M() * sd->N())
			throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INVALIDINDEX));
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
// Returns the first element in a nested cell array
//------------------------------------------------------------------------------
Currency unnest(Currency nested, const std::string& errmsg)
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
// removes zeros padding the outside of the matrix
//------------------------------------------------------------------------------
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

	try
	{
		if (max > 0)
			isprime = new bool[max - 1]; // offset to range from 2 to max inclusively
		else
			return p; 	// Zero and negatives return empty matrix.
	}
	catch (...)
	{ 
		throw OML_Error(HW_ERROR_OUTMEM);
	}

    for (int i = 0; i < max - 1; ++i)
        isprime[i] = true;

    int sqrtmax = (int) sqrt((double) max);
    for (int i = 2; i <= sqrtmax; i++)
    {
        if (isprime[i - 2])
        {
            for (int j = i * i; j <= max && j  > 0; j += i)  // j will routinely overflow to negative.  Just break this inner loop.
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
#if JDS
	const char* orig = file_name.data();

	wchar_t  wszDest[2048];
	MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, orig, file_name.length(), wszDest, 2048);

    DWORD dwAttrib = GetFileAttributesW(wszDest);
#endif
	DWORD dwAttrib = GetFileAttributes((LPCSTR)file_name.c_str());
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

    if (num_back_slahes % 2)
    {
        return true;
    }
    return false;
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
#if 0
bool isfinite(double d)
{
    return !(IsNaN_T(d) || isinfinity(d));
}
#endif
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
// Helper method for ismember when inputs are strings
//------------------------------------------------------------------------------
int stringVecFromCurrencyRows(EvaluatorInterface& eval, 
                              const Currency&     input1, 
                              const Currency&     input2, 
                              std::vector<Currency>& searchfor,
                               std::vector<Currency>& searchin)
{
    assert(input1.IsString());  // Checked earlier
    assert(input2.IsString());  // Checked earlier

    const hwMatrix* mtx1 = input1.Matrix();
    const hwMatrix* mtx2 = input2.Matrix();

    if (!mtx1 || !mtx2 || mtx1->N() != mtx2->N())
    {
        std::string msg("Error: invalid inputs in arguments 1 and 2; ");
        msg += "; string columns must match when searching by rows";
        throw OML_Error(msg);
    }

    int m = mtx1->M();

    int rows1 = mtx1->M();
    searchfor.reserve(rows1);
    for (int i = 0; i < rows1; ++i)
    {
        searchfor.emplace_back(BuiltInFuncsUtils::ReadRow(eval, input1, i));
    }
            
    int rows2 = mtx2->M();
    searchin.reserve(rows2);
    for (int i = 0; i < rows2; ++i)
    {
        searchin.emplace_back(BuiltInFuncsUtils::ReadRow(eval, input2, i));
    }

    return m;
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
                    throw OML_Error(OML_ERR_STRING_STRINGCELL, 2);
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
            throw OML_Error(OML_ERR_STRING_STRINGCELL, 2);
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
            throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
    }
    else
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void tryPushBackString(const HML_CELLARRAY *cell, std::vector<Currency> &topush)
{
    int csize = cell->Size();
    topush.reserve(csize);

    for (int i = 0; i < csize; ++i)
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
    int m = mtx->M();
    topush.reserve(m);

    for (int i = 0; i <m; ++i)
        topush.emplace_back(addStringMask(readRow(eval, mtx, i)));
}
//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void tryPushBackString(const hwMatrix *mtx, std::vector<Currency> &topush)
{
    int msize = mtx->Size();
    topush.reserve(msize);

    for (int i = 0; i < msize; ++i)
        topush.emplace_back((*mtx)(i));
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
                std::stable_sort(indexvec.begin(), indexvec.end(), [&indexvec, vec] (int i1, int i2) { return (*vec)(i1) < (*vec)(i2); });
            else
                std::stable_sort(indexvec.begin(), indexvec.end(), [&indexvec, vec] (int i1, int i2) { return (*vec)(i1) > (*vec)(i2); });

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
                std::stable_sort(indexvec.begin(), indexvec.end(), [&indexvec, vec] (int i1, int i2) { return complexLessThan(vec->z(i1), vec->z(i2)); });
            else
                std::stable_sort(indexvec.begin(), indexvec.end(), [&indexvec, vec] (int i1, int i2) { return complexGreaterThan(vec->z(i1), vec->z(i2)); });

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
                std::stable_sort(indexvec.begin(), indexvec.end(), [&indexvec, copy] (int i1, int i2) { return (*copy)(i1) < (*copy)(i2); });
            else
                std::stable_sort(indexvec.begin(), indexvec.end(), [&indexvec, copy] (int i1, int i2) { return (*copy)(i1) > (*copy)(i2); });
        }
        else
        {
            if (ASCEND)
                std::stable_sort(indexvec.begin(), indexvec.end(), [&indexvec, copy] (int i1, int i2) { return complexLessThan(copy->z(i1), copy->z(i2)); });
            else
                std::stable_sort(indexvec.begin(), indexvec.end(), [&indexvec, copy] (int i1, int i2) { return complexGreaterThan(copy->z(i1), copy->z(i2)); });
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
        std::stable_sort(indexvec.begin(), indexvec.end(), [&indexvec, &rowvec] (int i1, int i2) { return rowVecLessThan(rowvec[i1], rowvec[i2]); });
    else                                                      
        std::stable_sort(indexvec.begin(), indexvec.end(), [&indexvec, &rowvec] (int i1, int i2) { return rowVecGreaterThan(rowvec[i1], rowvec[i2]); });
    
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
// Returns regular expression string matching results
//------------------------------------------------------------------------------
std::vector<Currency> DoRegExp(EvaluatorInterface&             eval,
                               const std::string&              search, 
                               const std::string&              pattern, 
                               const std::vector<OMLREGEXP>&   options,
	                           const std::vector<std::string>& flags)
{
    std::string str (search);
    std::string pat (pattern);
    std::vector<Currency> outputs;
	// To cache regex, in Linux regex compilation at runtime is very costly
    static std::string previous_pat;
    static std::regex  reg(previous_pat, std::regex_constants::ECMAScript);
    static std::regex_constants::syntax_option_type  previous_syntaxflag = std::regex_constants::ECMAScript;

    // std::regex_search does not seem to work with '<' in pattern. So, for now remove it
    if (pattern.find("<") != std::string::npos)
    {
        size_t len = (pat.empty()) ? 0 : pat.size();
        std::string tmp;

        for (size_t i = 0; i < len; ++i)
        {
            if (pat[i] == '<' && i < len - 1 && pat[i + 1] == '=')
            {
                // Skip this as it leads to syntax errors in regex
            }
            else
            {
                tmp += pat[i];
            }
        }
        if (!tmp.empty())
        {
            pat = tmp;
        }
    }
    bool hasnames = (std::find(
        options.begin(), options.end(), OMLREGEXP_NAMES) != options.end());
    std::vector<std::string>    fieldnames;
    std::map< std::string, HML_CELLARRAY* > namesdata;
    std::unique_ptr<StructData> nm = nullptr;
    if (hasnames)
    {
        nm.reset(EvaluatorInterface::allocateStruct());
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

                        if (offset == std::string::npos)
                            throw OML_Error(HW_ERROR_SUBEXPNAMNOCLOSINGGR);

                        if (offset == sub_start + 3)
                            throw OML_Error(HW_ERROR_NOTEMPSUBEXPNAME);

                        std::string name(pat.substr(sub_start + 3, offset - (sub_start + 3)));
                        if (!name.empty())
                        {
                            fieldnames.push_back(name);
                            nm->addField(name);
                            namesdata[name] = nullptr;
                        }
                        pat.erase(pat.begin() + sub_start + 1, pat.begin() + offset + 1);
                        offset = sub_start + 1;
                        continue;
                    }
                }
            }
        }
    }
    BuiltInFuncsUtils utils;

    bool hasstart = (std::find(
        options.begin(), options.end(), OMLREGEXP_START) != options.end());
    std::unique_ptr<hwMatrix> startmtx = nullptr;

    bool hasend = (std::find(
        options.begin(), options.end(), OMLREGEXP_END) != options.end());
    std::unique_ptr<hwMatrix> endmtx = nullptr;

    bool hasmatch = (std::find(
        options.begin(), options.end(), OMLREGEXP_MATCH) != options.end());
    std::unique_ptr<HML_CELLARRAY> matchcell = nullptr;

    bool hassplit = (std::find(
        options.begin(), options.end(), OMLREGEXP_SPLIT) != options.end());
    std::unique_ptr<HML_CELLARRAY> splitcell = nullptr;

    bool hasext = (std::find(
        options.begin(), options.end(), OMLREGEXP_EXT) != options.end());
    std::unique_ptr<HML_CELLARRAY> extcell = nullptr;

    bool hastok = (std::find(
        options.begin(), options.end(), OMLREGEXP_TOK) != options.end());
    std::unique_ptr<HML_CELLARRAY> tokcell = nullptr;
            
    try
    {
		// Process syntax and match options

		// Case insensitive 
		std::regex_constants::syntax_option_type syntaxflag =
			std::regex_constants::ECMAScript;
		if (!flags.empty() &&
			std::find(flags.begin(), flags.end(), "ignorecase") != flags.end())
		{
			syntaxflag = std::regex_constants::ECMAScript | std::regex_constants::icase;
		}

		// First match only - std::regex_constants::format_first_only is ignored
		// by std::regex_search
		std::regex_constants::match_flag_type matchflag = std::regex_constants::match_default;
		bool once = (!flags.empty() &&
			std::find(flags.begin(), flags.end(), "once") != flags.end());

		int index = 0;
        std::smatch mtch;
		// To cache regex, in Linux regex compilation at runtime is very costly
        if ((0 != previous_pat.compare(pat)) || (previous_syntaxflag != syntaxflag))
        {
            reg.assign(pat, syntaxflag);
            previous_pat = pat;
            previous_syntaxflag = syntaxflag;
        }

        if (pat.length())
        {
            // find all matches
            while (1)
            {
                bool result = std::regex_search(str, mtch, reg, matchflag);
                if (!result)
                {
                    break;
                }
                std::string matchstr (mtch.str());
                size_t      matchpos = mtch.position();

                if (hasstart)
                {
                    double val = static_cast<double>(index + matchpos + 1);
                    if (!startmtx)
                    {
                        startmtx.reset(EvaluatorInterface::allocateMatrix(1, 1, val));
                    }
                    else
                    {
                        int n = startmtx->N();
                        utils.CheckMathStatus(eval, startmtx->Resize(1, n + 1));
                        (*startmtx)(n) = val;
                    }
                }

                if (hasend)
                {
                    double val = static_cast<double>(index + matchpos + matchstr.length());
                    if (!endmtx)
                    {
                        endmtx.reset(EvaluatorInterface::allocateMatrix(1, 1, val));
                    }
                    else
                    {
                        int n = endmtx->N();
                        utils.CheckMathStatus(eval, endmtx->Resize(1, n + 1));
                        (*endmtx)(n) = val;
                    }
                }

                if (hasmatch)
                {
                    int n = 0;
                    if (!matchcell)
                    {
                        matchcell.reset(EvaluatorInterface::allocateCellArray(1, 1));
                    }
                    else
                    {
                        n = matchcell->N();
                        utils.CheckMathStatus(eval, matchcell->Resize(1, n + 1));
                    }
                    (*matchcell)(0, n) = matchstr;
                }

                if (hassplit)
                {
                    int n = 0;
                    if (!splitcell)
                    {
                        splitcell.reset(EvaluatorInterface::allocateCellArray(1, 1));
                    }
                    else
                    {
                        n = splitcell->N();
                        utils.CheckMathStatus(eval, splitcell->Resize(1, n + 1));
                    }
                    (*splitcell)(0, n) = mtch.prefix().str();
                }

                std::unique_ptr<hwMatrix>      mtx  = nullptr;
                std::unique_ptr<HML_CELLARRAY> cell = nullptr;

                for (size_t i = 1; i < mtch.size(); ++i)
                {
                    if (mtch[i].matched)
                    {
                        int val = index + static_cast<int>(mtch.position(i));
                        if (hasext)
                        {
                            int m = 0;
                            double tokstart = static_cast<double>(index + mtch.position(i));
                            if (!mtx)
                            {
                                mtx.reset(EvaluatorInterface::allocateMatrix(1, 2, hwMatrix::REAL));
                            }
                            else
                            {
                                m = mtx->M();
                                utils.CheckMathStatus(eval, mtx->Resize(m + 1, 2));
                            }
                            (*mtx)(m, 0) = val + 1;
                            (*mtx)(m, 1) = val + static_cast<int>(mtch[i].length());
                        }
                        if (hastok)
                        {
                            int n = 0;
                            if (!cell)
                            {
                                cell.reset(EvaluatorInterface::allocateCellArray(1, 1));
                            }
                            else
                            {
                                n = cell->N();
                                utils.CheckMathStatus(eval, cell->Resize(1, n + 1));
                            }
                            (*cell)(n) = mtch[i].str();
                        }
                        if (hasnames && !namesdata.empty())
                        {
                            std::string val(mtch.str(i));
                            int n = 0;
                            std::string name = (i - 1 < fieldnames.size()) ? fieldnames[i - 1] : "";
                            HML_CELLARRAY* cell = namesdata[name];
                            if (!cell)
                            {
                                cell = EvaluatorInterface::allocateCellArray(1, 1);                               
                            }
                            else
                            {
                                n = cell->N();
                                utils.CheckMathStatus(eval, cell->Resize(1, n + 1));
                            }
                            (*cell)(n) = val;
                            namesdata[name] = cell; 
                        }
                    }

                }
                if (hasext)
                {
                    int n = 0;
                    if (!extcell)
                    {
                        extcell.reset(EvaluatorInterface::allocateCellArray(1, 1));
                    }
                    else
                    {
                        n = extcell->N();
                        utils.CheckMathStatus(eval, extcell->Resize(1, n + 1));
                    }
                    if (!mtx)
                    {
                        mtx.reset(EvaluatorInterface::allocateMatrix(0, 2, hwMatrix::REAL));
                    }
                    (*extcell)(0, n) = mtx.release();
                }
                if (hastok)
                {
                    int n = 0;
                    if (!tokcell)
                    {
                        tokcell.reset(EvaluatorInterface::allocateCellArray(1, 1));
                    }
                    else
                    {
                        n = tokcell->N();
                        utils.CheckMathStatus(eval, tokcell->Resize(1, n + 1));
                    }
                    if (!cell)
                    {
                        cell.reset(EvaluatorInterface::allocateCellArray(1, 0));
                    }
                    (*tokcell)(n) = cell.release();
                }

				if (once && mtch.length() >= 1)
				{
					break;
				}
                index += (int)(matchpos + matchstr.length());

                std::string updatedSearchStr (mtch.suffix().str());

                // Fix for memory leak if input contains '\n'
                if (updatedSearchStr == str)
                {
                    break;
                }
                str = updatedSearchStr;
            }
        }

        if (hassplit)
        {
            int n = 0;
            if (!splitcell)
            {
                splitcell.reset(EvaluatorInterface::allocateCellArray(
                    1, 1));
            }
            else
            {
                n = splitcell->N();
                utils.CheckMathStatus(eval, splitcell->Resize(1, n + 1));
            }
            (*splitcell)(0, n) = str;
        }
    }
    catch (std::regex_error& err)
    {
        BuiltInFuncsUtils utils;
        utils.ThrowRegexError(err.code());
    }
    catch (const OML_Error& e)
    {
        throw e;
    }
    catch (...)
    {
        throw OML_Error("Error: unknown error from regexp");
    }


    // Set outputs with user specified order
    Currency startcur;
    if (hasstart)
    {
        if (!startmtx)
        {
            startmtx.reset(EvaluatorInterface::allocateMatrix(1, 0, hwMatrix::REAL));
        }
        startcur = startmtx.release();
    }

    Currency endcur;
    if (hasend)
    {
        if (!endmtx)
        {
            endmtx.reset(EvaluatorInterface::allocateMatrix(1, 0, hwMatrix::REAL));
        }
        endcur = endmtx.release();
    }

    Currency matchcur;
    if (hasmatch)
    {
        if (!matchcell)
        {
            matchcell.reset(EvaluatorInterface::allocateCellArray());
        }
        matchcur = matchcell.release();
    }

    Currency split;
    if (hassplit)
    {
        if (!splitcell)
        {
            splitcell.reset(EvaluatorInterface::allocateCellArray());
        }
        split = splitcell.release();
    }

    Currency ext;
    if (hasext)
    {
        if (!extcell)
        {
            extcell.reset(EvaluatorInterface::allocateCellArray());
        }
        ext = extcell.release();
    }

    Currency tok;
    if (hastok)
    {
        if (!tokcell)
        {
            tokcell.reset(EvaluatorInterface::allocateCellArray());
        }
        tok = tokcell.release();
    }

    Currency nmcur;
    if (hasnames)
    {
        std::map<std::string, HML_CELLARRAY* > ::const_iterator itr = namesdata.begin();

        for (; itr != namesdata.end(); ++itr)
        {
            std::string name(itr->first);
            HML_CELLARRAY* tmp = itr->second;
            if (!tmp || tmp->IsEmpty())
            {
                nm->SetValue(0, 0, itr->first, "");
            }
            else if (tmp->Size() == 1)
            {
                nm->SetValue(0, 0, itr->first, (*tmp)(0).StringVal());
            }
            else
            {
                nm->SetValue(0, 0, itr->first, tmp);
            }
        }
        nmcur = nm.release();
    }
    for (std::vector<OMLREGEXP>::const_iterator itr = options.begin();
         itr != options.end(); ++itr)
    {
        OMLREGEXP option (*itr);

        switch (option)
        {
            case OMLREGEXP_START: outputs.push_back(startcur); break;
            case OMLREGEXP_END:   outputs.push_back(endcur);   break;
            case OMLREGEXP_EXT:   outputs.push_back(ext);      break;
            case OMLREGEXP_MATCH: outputs.push_back(matchcur); break;
            case OMLREGEXP_TOK:   outputs.push_back(tok);      break;
            case OMLREGEXP_NAMES: outputs.push_back(nmcur);    break;
            case OMLREGEXP_SPLIT: outputs.push_back(split);    break;
            default: break;
        }
    }

    return outputs;
}
//------------------------------------------------------------------------------
// Does regular expression string matching
//------------------------------------------------------------------------------
void GetRegExpOutput( EvaluatorInterface              eval,
                      const Currency&                 searchCur,
                      const Currency&                 patternCur,
                      const std::vector<OMLREGEXP>&   options,
	                  const std::vector<std::string>& flags,
                      std::vector<Currency>&          outputs)
{
    HML_CELLARRAY *cell1 = 0;
    HML_CELLARRAY *cell2 = 0;
    std::string str1, str2;

    Currency cur1 = toCurrencyStr(eval, searchCur, false, true);
    if (cur1.IsString())
    {
        str1 = cur1.StringVal();
    }
    else if (cur1.IsCellArray())
    {
        cell1 = cur1.CellArray();
    }
    else
    {
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 1, OML_VAR_TYPE);
    }
       
    Currency cur2 = toCurrencyStr(eval, patternCur, false, true);
    if (cur2.IsString())
    {
        str2 = readString(cur2);
    }
    else if (cur2.IsCellArray())
    {
        cell2 = cur2.CellArray();
    }
    else
    {
        throw OML_Error(OML_ERR_STRING_STRINGCELL, 2, OML_VAR_TYPE);
    }

    if (cell2 && cell2->Size() == 1)
    {
        Currency temp = Unnest(cur2, OML_ERR_STRING_STRINGCELL, 2);
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
        outputs = DoRegExp(eval, str1, str2, options, flags);
        return;
    }

    if (cell1)
    {
        if (cell2)
        {
            int cell2size = cell2->Size();
            int m = std::max(cell1->M(), cell2->M());
            int n = std::max(cell1->N(), cell2->N());
            for (int k = 0; k < cell2size; ++k)
            {
                Currency elem2 = Unnest((*cell2)(k), OML_ERR_STRING_STRINGCELL, 2);
                if (!elem2.IsString())
                {
                    throw OML_Error(OML_ERR_STRING_STRINGCELL, 2);
                }
                DoRegExp(eval, cell1, readString(elem2), options, flags, k, m, n, outputs);
            }

        }
        else
        {
            DoRegExp(eval, cell1, str2, options, flags, -1, cell1->M(), cell1->N(), outputs);
            return;
        }
    }
    else if (cell2)
    {
        for (int i = 0; i < cell2->Size(); i++)
        {
            Currency elem2 = unnest((*cell2)(i), HW_ERROR_INPUTSTRINGCELLARRAY);
            if (!elem2.IsString())
                throw OML_Error(HW_ERROR_INPUTSTRINGCELLARRAY);

            std::vector<Currency> outvec = DoRegExp(eval, str1, readString(elem2), options, flags);

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
                size_t curdimssize = curdims.size();
                size_t dimssize    = dims.size();
                for (int i = 0; i < numDims; ++i)
                {
                    if (i >= static_cast<int>(curdimssize))
                    {
                        throw OML_Error(OML_ERR_ARRAYCATDIM, 2, i + 1);
                    }
                    if (i == dim-1)
                    {
                        dims[dim-1] += curdims[i];
                    }
                    else
                    {
                        if (i >= static_cast<int>(dimssize) || 
                            dims[i] != curdims[i])
                        {
                            throw OML_Error(OML_ERR_ARRAYCATDIM, 2, i + 1);
                        }
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
                temp.Convert2DtoND(*matrix, false);
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
// Helper function for unique command
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
// Helper function for unique command for 2D matrices or strings
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
// Helper function for unique command for cell arrays
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
// Throws an error without stack information [errormsgonly]
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
//------------------------------------------------------------------------------
// Sets a warning without stack information [warningmsgonly]
//------------------------------------------------------------------------------
bool oml_warningmsgonly(EvaluatorInterface           eval, 
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
    if (lower.find("warning: ", 0, 9) == std::string::npos)
    {
        msg.insert(0, "Warning: ");
    }
    Currency warn(msg);
    warn.DispOutput();
    eval.PrintResult(warn);
    EvaluatorInterface::SetLastWarning(msg);

    return true;
}
//------------------------------------------------------------------------------
// Gets ordered vector of output options to display for regexp
//------------------------------------------------------------------------------
std::vector<OMLREGEXP> GetRegexOptions(const std::vector<Currency>& inputs,
                                       int                          nargout,
	                                   std::vector<std::string>&    flags)    
{
    assert(!inputs.empty());

    bool hasstart = false;
    bool hasend = false;
    bool hasext = false;
    bool hasmatch = false;
    bool hastok = false;
    bool hasnames = false;
    bool hassplit = false;

    size_t numoutputs = (nargout == 0) ? 1 : nargout; // At least one output
    std::vector<OMLREGEXP> options;
    options.reserve(nargout);
    
    if (inputs.size() > 2)
    {
        // First currency in input is search string
        // Second currency in input is pattern string
        // Rest are options for outputs specified by user
        int i      = 3;
        int outidx = 0;
        for (std::vector<Currency>::const_iterator itr = inputs.begin() + 2;
			itr != inputs.end(); ++itr, ++i)
        {
            if (!(*itr).IsString())
            {
                throw OML_Error(OML_ERR_STRING, i);
            }
            std::string opt = (*itr).StringVal();
            if (opt.empty())
            {
                throw OML_Error(OML_ERR_NONEMPTY_STR, i);
            }

            std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);

            if (opt == "start")
            {
                hasstart = true;
				if (outidx < nargout)
				{
					options.push_back(OMLREGEXP_START);
					outidx++;
				}
                continue;
            }
            else if (opt == "end")
            {
                hasend = true;
				if (outidx < nargout)
				{
					options.push_back(OMLREGEXP_END);
					outidx++;
				}
                continue;
            }
            else if (opt == "tokenextents")
            {
                hasext = true;
				if (outidx < nargout)
				{
					options.push_back(OMLREGEXP_EXT);
					outidx++;
				}
                continue;
            }
            else if (opt == "match")
            {
                hasmatch = true;
				if (outidx < nargout)
				{
					options.push_back(OMLREGEXP_MATCH);
					outidx++;
				}
                continue;
            }
            else if (opt == "tokens")
            {
                hastok = true;
				if (outidx < nargout)
				{
					options.push_back(OMLREGEXP_TOK);
					outidx++;
				}
                continue;
            }
            else if (opt == "names")
            {
                hasnames = true;
				if (outidx < nargout)
				{
					options.push_back(OMLREGEXP_NAMES);
					outidx++;
				}
                continue;
            }
            else if (opt == "split")
            {
                hassplit = true;
				if (outidx < nargout)
				{
					options.push_back(OMLREGEXP_SPLIT);
					outidx++;
				}
                continue;
            }
			else if (opt == "ignorecase")
			{
				if (flags.empty() ||
					std::find(flags.begin(), flags.end(), "ignorecase") == flags.end())
				{
					flags.push_back("ignorecase");
				}
				continue;
			}
			else if (opt == "matchcase")
			{
				if (flags.empty())
				{
					continue;  // Nothing to change as case sensitive comparison is default
				}
				std::vector<std::string>::iterator fitr = std::find(
					flags.begin(), flags.end(), "ignorecase");
				if (fitr != flags.end())
				{
					flags.erase(fitr);
				}
				continue;
			}
			else if (opt == "once")
			{
				flags.push_back("once");
				continue;
			}
			throw OML_Error(OML_ERR_OPT_UNSUPPORTED, i);
        }
    }

    // Now push in other options    
    if (!hasstart && options.size() < numoutputs)   // Start indices of each matching substring
    {
        options.push_back(OMLREGEXP_START);
    }
    if (!hasend && options.size() < numoutputs)    // End indices of each matching substring
    {
        options.push_back(OMLREGEXP_END);
    }
    if (!hasext && options.size() < numoutputs)    // Extents of matching tokens surrounded by (..)
    {
        options.push_back(OMLREGEXP_EXT);
    }
    if (!hasmatch && options.size() < numoutputs)  // Cell array of the text of each matching string
    {
        options.push_back(OMLREGEXP_MATCH);
    }
    if (!hastok && options.size() < numoutputs)    // Cell array of each token matched
    {
        options.push_back(OMLREGEXP_TOK);
    }
    if (!hasnames && options.size() < numoutputs)  // Text of each matched named token
    {
        options.push_back(OMLREGEXP_NAMES);
    }
    if (!hassplit && options.size() < numoutputs) // Cell array of what is not returned by match
    {
        options.push_back(OMLREGEXP_SPLIT);
    }
    return options;
}

//------------------------------------------------------------------------------
// Function placeholder - implementation is on client side
//------------------------------------------------------------------------------
bool oml_getargc(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs)
{
    outputs.push_back(1);
    return true;
}
//------------------------------------------------------------------------------
// Function placeholder - implementation is on client side
//------------------------------------------------------------------------------
bool oml_getargv(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs)
{
    outputs.push_back("");
    return true;
}
//------------------------------------------------------------------------------
// Returns true, gets the days, with Jan 1 0000 considered as day 1 [now]
//------------------------------------------------------------------------------
bool oml_now(EvaluatorInterface           eval,
             const std::vector<Currency>& inputs,
             std::vector<Currency>&       outputs)
{
    if (!inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    time_t rawtime = time(nullptr);
    struct tm* curtime = localtime(&rawtime);

    double year = curtime->tm_year + 1900;
    double mon  =  curtime->tm_mon;
    double day  = curtime->tm_mday;

    year += static_cast<int>(mon/12.0);
    mon = rem(mon, 12.0);

    // deal with leap years
    double result = year-- * 365;
    result += floor(year / 4.0) - floor(year / 100.0) + floor(year / 400.0) + 1;

    //months
    result += daysInMonths((int)mon, (int)(year + 1.0));

    // days
    result += day;

    // Integer part of result will give the number of days and the fractional
    // part of results should give the number of seconds
    struct tm start;
    start.tm_year = curtime->tm_year;
    start.tm_mon = curtime->tm_mon;
    start.tm_mday = curtime->tm_mday;
    start.tm_hour = 0;
    start.tm_min = 0;
    start.tm_sec = 0;
    start.tm_isdst = -1;

    std::chrono::system_clock::time_point today =
        std::chrono::system_clock::from_time_t(std::mktime(&start));
    std::chrono::system_clock::duration d = std::chrono::system_clock::now() - today;

    std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(d);
    long long secs = s.count();
    int numdigits  = (secs < 1) ? 1 : static_cast<int>(log10(secs) + 1);

    result += secs / pow(10, numdigits);
    
    outputs.push_back(result);    
   
    return true;
}

bool oml_parcluster(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error(OML_ERR_NUMARGIN);

	Currency val = inputs[0];

	if (!val.IsScalar())
		throw OML_Error(OML_ERR_SCALAR);

	if (!val.IsPositiveInteger())
		throw OML_Error(OML_ERR_POSINTEGER);

	eval.SetNumberOfThreads((int)val.Scalar());

	return true;
}

//------------------------------------------------------------------------------
// Reads from file [fread]
//------------------------------------------------------------------------------
bool oml_fread(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    BuiltInFuncsUtils utils;
    int fileID = utils.GetFileId(eval, inputs[0], 1);
    utils.CheckFileIndex(eval, fileID, 1, true);

    int blockSize = 1;
    int maxLoops = -1; // unlimited
    int nrows = -1;
    int ncols = -1;
    int skip = 0;

    bool        signedOutput = false;
    DataType    dtype = Char;
    size_t      size = sizeof(unsigned char);
    size_t      nargin = inputs.size();
    std::string origPrec;

    try
    {
        if (nargin > 1)
        {
            const Currency& input2 = inputs[1];
            if (nargin == 2 && input2.IsString())
            {
                origPrec = input2.StringVal();
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
                        maxLoops = (int)round(temp);
                }
                else if (input2.IsMatrix())
                {
                    const hwMatrix* m = input2.Matrix();
                    if (m->IsReal())
                    {
                        if (m->Size() == 2)
                        {
                            nrows = (int)round((*m)(0));
                            ncols = (int)round((*m)(1));
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
                {
                    std::string msg("Error: invalid input in argument 2; must be a positive integer,");
                    msg += " a vector of positive integers or infinity";
                    throw OML_Error(msg);
                }
            }
            if (nargin > 2)
            {
                const Currency& input3 = inputs[2];
                origPrec = readOption(eval, unnest(input3));
                Precision p = getPrecision(eval, input3);
                signedOutput = p.sign;
                blockSize = p.blockSize;
                size = p.numBytes;
                dtype = p.dtype;

                if (nargin > 3)
                {
                    skip = (int)inputs[3].Scalar();

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
            maxLoops = (int)ceil(maxLoops / (double)blockSize);
        }

        if (!nrows || !ncols || !blockSize || !file)
        {
            outputs.push_back(EvaluatorInterface::allocateMatrix());
            outputs.push_back(0);
            return true;
        }

        outputs = FreadData(eval, dtype, signedOutput, maxLoops, size,
            blockSize, skip, file, nrows, ncols);
    }
    catch (const OML_Error & e)
    {
        eval.CloseFile(fileID);  // \todo: Should not be closing files
        throw e;
    }

    if (!origPrec.empty() && !outputs.empty())
    {
        std::transform(origPrec.begin(), origPrec.end(), origPrec.begin(), ::tolower);
        if (origPrec.find("=>char") != std::string::npos)
        {
            outputs[0].SetMask(Currency::MASK_STRING);
        }
    }

    return true;
}
//------------------------------------------------------------------------------
// Helper method to do regular expression matching
//------------------------------------------------------------------------------
void DoRegExp(EvaluatorInterface&             eval,
              HML_CELLARRAY*                  cell,
              const std::string&              pattern,
              const std::vector<OMLREGEXP>&   options,
              const std::vector<std::string>& flags,
              int                             idx,
              int                             m,
              int                             n,
              std::vector<Currency>&          outputs)
{
    if (!cell)
    {
        return;
    }

    int cellsize = cell->Size();

    for (int i = 0; i < cellsize; ++i)
    {
        Currency elem1 = Unnest((*cell)(i), OML_ERR_STRING_STRINGCELL, 1);
        if (!elem1.IsString())
        {
            throw OML_Error(OML_ERR_STRING_STRINGCELL, 1);
        }
        std::vector<Currency> outvec (DoRegExp(eval, readString(elem1), pattern, 
            options, flags));

        if (!outvec.empty())
        {
            int numoutputs = static_cast<int>(outvec.size());
            if (i == 0 && outputs.size() < static_cast<size_t>(numoutputs))
            {
                outputs.reserve(numoutputs);
                for (int j = 0; j < numoutputs; ++j)
                {
                    outputs.push_back(EvaluatorInterface::allocateCellArray(
                        m, n));
                }
            }
            for (int j = 0; j < numoutputs; ++j)
            {
                int elem = (idx == -1) ? i : idx;
                if (!outputs[j].IsCellArray() || 
                    (outputs[j].CellArray())->Size() < elem)
                {
                    continue;
                }
                (*outputs[j].CellArray())(elem) = outvec[j];
            }

        }
    }
}
//------------------------------------------------------------------------------
// Returns the first element in a nested cell array
//------------------------------------------------------------------------------
Currency Unnest(Currency nested, omlMathErrCode err, int idx)
{
    while (nested.IsCellArray())
    {
        HML_CELLARRAY* cell = nested.CellArray();
        if (cell->IsEmpty())
        {
            if (idx > 0)
            {
                throw OML_Error(err, idx);
            }
            else
            {
                throw OML_Error(err);
            }
        }
        nested = (*cell)(0);
    }
    return nested;
}
//------------------------------------------------------------------------------
// Helper method for file read operations
//------------------------------------------------------------------------------
std::vector<Currency> FreadData(EvaluatorInterface eval, DataType type, 
    bool sign, int nLoops, size_t size, int count, int skip, 
    std::FILE* file, int nrows, int ncols)
{
    // If size is not specified, everything will be pushed to a 
    // column vector and it will be reshaped/resized later, based on the data
    std::unique_ptr<hwMatrix> mtx(EvaluatorInterface::allocateMatrix(
        1, 1, hwMatrix::REAL));
    bool grow = (nrows < 0 || ncols < 0);
    if (!grow)
    {
        mtx.reset(EvaluatorInterface::allocateMatrix(
            nrows, ncols, hwMatrix::REAL));
    }
    mtx->SetElements(0);

    int total     = 0;
    int origncols = ncols;
    int idx       = 0;

    for (int i = 0; i != nLoops && !feof(file); ++i)
    {
        bool result = false;
        switch (type)
        {
            case Double: 
                result = ReadDoubleBlock(eval, file, grow, size, count, mtx.get(), total, idx);
                break;

            case Float:
                result = ReadFloatBlock(eval, file, grow, size, count, mtx.get(), total, idx);
                break;

            case Int: 
                result = ReadIntBlock(eval, file, sign, grow, size,
                    count, mtx.get(), total, idx);
                break;

            case Int8:
                result = ReadInt8Block(eval, file, sign, grow, size,
                    count, mtx.get(), total, idx);
                break;

            case Int16:
                result = ReadInt16Block(eval, file, sign, grow, size,
                    count, mtx.get(), total, idx);
                break;

            case Int32:
                result = ReadInt32Block(eval, file, sign, grow, size,
                    count, mtx.get(), total, idx);
                break;

            case Int64:
                result = ReadInt64Block(eval, file, sign, grow, size, count, 
                    mtx.get(), total, idx);
                break;

            case Long:
                result = ReadLongBlock(eval, file, sign, grow, size, count,
                    mtx.get(), total, idx);
                break;

            case LongLong:
                result = ReadLongLongBlock(eval, file, sign, grow, size, count,
                    mtx.get(), total, idx);
                break;

            case Short:
                result = ReadShortBlock(eval, file, sign, grow, size,
                    count, mtx.get(), total, idx);
                break;

            case Char:
                result = ReadCharBlock(eval, file, sign, grow, size,
                    count, mtx.get(), total, idx);
                break;

            default:
                throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR)); 
                break;
        }
        if (!result)
        {
            break;
        }

        if (skip)
        {
            fseek(file, skip, SEEK_CUR);
        }
    }

    if (ferror(file))
    {
        total = -1;
    }

    std::vector<Currency> outputs;
    if (total < 1)
    {
        outputs.push_back(EvaluatorInterface::allocateMatrix());
    }
    else
    {
        if (grow)
        {
            if (nrows > 0)
            {
                ncols = static_cast<int>(ceil(total / static_cast<double>(nrows)));
            }
            else if (origncols > 0)
            {
                nrows = static_cast<int>(ceil(total / static_cast<double>(origncols)));
            }
            mtx->Reshape(nrows, ncols);
        }

        // Prevents unnecessary columns from being created if data does not fill
        // all columns specified by user
        if (nrows > -1)
        {
            ncols = (int)min((double)origncols, ceil((double)(total / (double)nrows)));
        }
        if (ncols < origncols && ncols > 0)
        {
            mtx->Resize(mtx->M(), ncols, true);
        }

        outputs.push_back(mtx.release());
    }
    outputs.push_back(total);
    return outputs;
}
//------------------------------------------------------------------------------
// Returns result after reading an unsigned/signed char block of data
//------------------------------------------------------------------------------
bool ReadCharBlock(EvaluatorInterface eval, std::FILE* file, bool sign, 
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx) 
{
    int msize = mtx->Size();

    if (!sign)
    {
        std::unique_ptr<unsigned char[]> buf(new unsigned char[count + 1]);
        size_t read = fread(buf.get(), size, count, file);
        if (read <= 0 || read > count)
        {
            return false;
        }
        total += static_cast<int>(read);

        if (grow && msize < total)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
            msize = mtx->Size();
        }

        for (size_t i = 0; i < read && idx < msize; ++i)
        {
            (*mtx)(idx++) = (buf.get())[i];
        }
        return true;
    }

    std::unique_ptr<signed char[]> buf(new signed char[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = (buf.get())[i];
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns result after reading a double block of data
//------------------------------------------------------------------------------
bool ReadDoubleBlock(EvaluatorInterface eval, std::FILE* file, 
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx)
{
    int msize = mtx->Size();
    std::unique_ptr<double[]> buf(new double[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = (buf.get())[i];
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns result after reading a float block of data
//------------------------------------------------------------------------------
bool ReadFloatBlock(EvaluatorInterface eval, std::FILE* file, 
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx)
{
    int msize = mtx->Size();
    std::unique_ptr<float[]> buf(new float[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = (buf.get())[i];
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns result after reading an unsigned/signed int block of data
//------------------------------------------------------------------------------
bool ReadIntBlock(EvaluatorInterface eval, std::FILE* file, bool sign,
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx)
{
    int msize = mtx->Size();
    if (!sign)
    {
        std::unique_ptr<unsigned int[]> buf(new unsigned int[count]);
        size_t read = fread(buf.get(), size, count, file);
        if (read <= 0 || read > count)
        {
            return false;
        }

        total += static_cast<int>(read);

        if (grow && msize < total)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
            msize = mtx->Size();
        }

        for (size_t i = 0; i < read && idx < msize; ++i)
        {
            (*mtx)(idx++) = (buf.get())[i];
        }
        return true;
    }

    std::unique_ptr<signed int[]> buf(new signed int[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = (buf.get())[i];
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns result after reading an unsigned/signed short block of data
//------------------------------------------------------------------------------
bool ReadShortBlock(EvaluatorInterface eval, std::FILE* file, bool sign,
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx)
{
    int msize = mtx->Size();
    if (!sign)
    {
        std::unique_ptr<unsigned short[]> buf(new unsigned short[count]);
        size_t read = fread(buf.get(), size, count, file);
        if (read <= 0 || read > count)
        {
            return false;
        }

        total += static_cast<int>(read);

        if (grow && msize < total)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
            msize = mtx->Size();
        }

        for (size_t i = 0; i < read && idx < msize; ++i)
        {
            (*mtx)(idx++) = (buf.get())[i];
        }
        return true;
    }

    std::unique_ptr<signed short[]> buf(new signed short[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = (buf.get())[i];
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns result after reading an unsigned/signed long block of data
//------------------------------------------------------------------------------
bool ReadLongBlock(EvaluatorInterface eval, std::FILE* file, bool sign,
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx)
{
    int msize = mtx->Size();
    if (!sign)
    {
        std::unique_ptr<unsigned long[]> buf(new unsigned long[count]);
        size_t read = fread(buf.get(), size, count, file);
        if (read <= 0 || read > count)
        {
            return false;
        }

        total += static_cast<int>(read);

        if (grow && msize < total)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
            msize = mtx->Size();
        }

        for (size_t i = 0; i < read && idx < msize; ++i)
        {
            (*mtx)(idx++) = (buf.get())[i];
        }
        return true;
    }

    std::unique_ptr<signed long[]> buf(new signed long[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = (buf.get())[i];
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns result after reading an unsigned/signed long long block of data
//------------------------------------------------------------------------------
bool ReadLongLongBlock(EvaluatorInterface eval, std::FILE* file, bool sign,
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx)
{
    int msize = mtx->Size();
    if (!sign)
    {
        std::unique_ptr<unsigned long long[]> buf(new unsigned long long[count]);
        size_t read = fread(buf.get(), size, count, file);
        if (read <= 0 || read > count)
        {
            return false;
        }

        total += static_cast<int>(read);

        if (grow && msize < total)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
            msize = mtx->Size();
        }

        for (size_t i = 0; i < read && idx < msize; ++i)
        {
            (*mtx)(idx++) = static_cast<double>((buf.get())[i]);
        }
        return true;
    }

    std::unique_ptr<signed long long[] > buf(new signed long long[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = static_cast<double>((buf.get())[i]);
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns result after reading an unsigned/signed int8_t block of data
//------------------------------------------------------------------------------
bool ReadInt8Block(EvaluatorInterface eval, std::FILE* file, bool sign,
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx)
{
    int msize = mtx->Size();
    if (!sign)
    {
        std::unique_ptr<uint8_t[]> buf(new uint8_t[count]);
        size_t read = fread(buf.get(), size, count, file);
        if (read <= 0 || read > count)
        {
            return false;
        }

        total += static_cast<int>(read);

        if (grow && msize < total)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
            msize = mtx->Size();
        }

        for (size_t i = 0; i < read && idx < msize; ++i)
        {
            (*mtx)(idx++) = (buf.get())[i];
        }
        return true;
    }

    std::unique_ptr<int8_t[]> buf(new int8_t[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = (buf.get())[i];
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns result after reading an unsigned/signed int16_t block of data
//------------------------------------------------------------------------------
bool ReadInt16Block(EvaluatorInterface eval, std::FILE* file, bool sign,
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx)
{
    int msize = mtx->Size();
    if (!sign)
    {
        std::unique_ptr<uint16_t[]> buf(new uint16_t[count]);
        size_t read = fread(buf.get(), size, count, file);
        if (read <= 0 || read > count)
        {
            return false;
        }

        total += static_cast<int>(read);

        if (grow && msize < total)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
            msize = mtx->Size();
        }

        for (size_t i = 0; i < read && idx < msize; ++i)
        {
            (*mtx)(idx++) = (buf.get())[i];
        }
        return true;
    }

    std::unique_ptr<int16_t[]> buf(new int16_t[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = (buf.get())[i];
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns result after reading an unsigned/signed int32_t block of data
//------------------------------------------------------------------------------
bool ReadInt32Block(EvaluatorInterface eval, std::FILE* file, bool sign,
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx)
{
    int msize = mtx->Size();
    if (!sign)
    {
        std::unique_ptr<uint32_t[]> buf(new uint32_t[count]);
        size_t read = fread(buf.get(), size, count, file);
        if (read <= 0 || read > count)
        {
            return false;
        }

        total += static_cast<int>(read);

        if (grow && msize < total)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
            msize = mtx->Size();
        }

        for (size_t i = 0; i < read && idx < msize; ++i)
        {
            (*mtx)(idx++) = (buf.get())[i];
        }
        return true;
    }

    std::unique_ptr<int32_t[]> buf(new int32_t[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = (buf.get())[i];
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns result after reading an unsigned/signed int64_t block of data
//------------------------------------------------------------------------------
bool ReadInt64Block(EvaluatorInterface eval, std::FILE* file, bool sign,
    bool grow, size_t size, int count, hwMatrix* mtx, int& total, int& idx)
{
    int msize = mtx->Size();
    if (!sign)
    {
        std::unique_ptr<uint64_t[]> buf(new uint64_t[count]);
        size_t read = fread(buf.get(), size, count, file);
        if (read <= 0 || read > count)
        {
            return false;
        }

        total += static_cast<int>(read);

        if (grow && msize < total)
        {
            BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
            msize = mtx->Size();
        }

        for (size_t i = 0; i < read && idx < msize; ++i)
        {
            (*mtx)(idx++) = static_cast<double>((buf.get())[i]);
        }
        return true;
    }

    std::unique_ptr<int64_t[]> buf(new int64_t[count]);
    size_t read = fread(buf.get(), size, count, file);
    if (read <= 0 || read > count)
    {
        return false;
    }

    total += static_cast<int>(read);

    if (grow && msize < total)
    {
        BuiltInFuncsUtils::CheckMathStatus(eval, mtx->Resize(total, 1, true));
        msize = mtx->Size();
    }

    for (size_t i = 0; i < read && idx < msize; ++i)
    {
        (*mtx)(idx++) = static_cast<double>((buf.get())[i]);
    }
    return true;
}

// If we keep all of these functions to define an AST in memory, it probably makes sense to try memory pooling the OMLTree
// since we're allocating a bunch of them.
bool oml_p_definefunction(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 3)
		throw OML_Error("Illegal function definition");

	Currency func_name = inputs[0];
	Currency func_inputs = inputs[1];
	Currency func_outputs = inputs[2];

	if (!func_name.IsString())
		throw OML_Error("Illegal function definition");

	if (!func_inputs.IsCellArray())
		throw OML_Error("Illegal function definition");

	if (!func_outputs.IsCellArray())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(FUNC_DEF, "FUNC_DEF", NULL, 0, 5);
	my_tree->AddChild(new OMLTree(IDENT, func_name.StringVal(), NULL, 0, 0));

	HML_CELLARRAY* cells = func_inputs.CellArray();

	OMLTree* input_tree = new OMLTree(ID_LIST, "ID_LIST", NULL, 0, cells->Size());

	for (int j = 0; j < cells->Size(); j++)
	{
		Currency temp = (*cells)(j);

		if (!temp.IsString())
			throw OML_Error("Illegal function definition");

		input_tree->AddChild(new OMLTree(IDENT, temp.StringVal(), NULL, 0, 0));
	}

	cells = func_outputs.CellArray();

	OMLTree* output_tree = new OMLTree(ID_LIST, "ID_LIST", NULL, 0, cells->Size());

	for (int j = 0; j < cells->Size(); j++)
	{
		Currency temp = (*cells)(j);

		if (!temp.IsString())
			throw OML_Error("Illegal function definition");

		output_tree->AddChild(new OMLTree(IDENT, temp.StringVal(), NULL, 0, 0));
	}

	my_tree->AddChild(input_tree);
	my_tree->AddChild(output_tree);

	my_tree->AddChild(new OMLTree(STATEMENT_LIST, "STATEMENT_LIST", NULL, 0, 0));
	my_tree->AddChild(new OMLTree(END, "end", NULL, 0, 0));

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_registerfunction(EvaluatorInterface ei, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal function definition");

	Currency func_def = inputs[0];

	if (!func_def.IsTree())
		throw OML_Error("Illegal function definition");

	ei.RegisterFunc(func_def.Tree());

	return true;
}

bool oml_p_assign(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	rhs.ConvertToTree();

	if (!lhs.IsTree() && !lhs.IsString())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* assign_tree = new OMLTree(ASSIGN, "=", NULL, 0, 2);

	if (lhs.IsString())
		assign_tree->AddChild(new OMLTree(IDENT, lhs.StringVal(), NULL, 0, 0));
	else
		assign_tree->AddChild(lhs.Tree());

	assign_tree->AddChild(rhs.Tree());

	outputs.push_back(assign_tree);

	return true;
}

bool oml_p_multassign(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	rhs.ConvertToTree();

	if (!lhs.IsCellArray())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* assign_tree = new OMLTree(MR_FUNC, "MR_FUNC", NULL, 0, 2);

	HML_CELLARRAY* cells = lhs.CellArray();
	OMLTree* params_tree = new OMLTree(ID_LIST, "ID_LIST", NULL, 0, cells->Size());

	for (int j = 0; j < cells->Size(); ++j)
	{
		Currency temp = (*cells)(j);

		if (temp.IsString())
		{
			OMLTree* temp_tree = new OMLTree(IDENT, temp.StringVal(), NULL, 0, 1);
			params_tree->AddChild(temp_tree);
		}
	}

	assign_tree->AddChild(params_tree);
	assign_tree->AddChild(rhs.Tree());

	outputs.push_back(assign_tree);

	return true;
}

bool oml_p_struct(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() < 2)
		throw OML_Error("Illegal function definition");

	OMLTree* assign_tree = NULL;

	for (int j = 1; j<inputs.size(); ++j)
	{
		if (!assign_tree)
		{
			Currency temp = inputs[j-1];

			if (!temp.IsString())
				throw OML_Error("Illegal struct definition");

			assign_tree = new OMLTree(STRUCT, "STRUCT", NULL, 0, 2);
			assign_tree->AddChild(new OMLTree(IDENT, temp.StringVal(), NULL, 0, 0));

			Currency temp2 = inputs[j];

			if (!temp2.IsString())
				throw OML_Error("Illegal struct definition");

			assign_tree->AddChild(new OMLTree(IDENT, temp2.StringVal(), NULL, 0, 0));
		}
		else
		{
			Currency temp = inputs[j];

			OMLTree* new_tree = new OMLTree(STRUCT, "STRUCT", NULL, 0, 2);

			new_tree->AddChild(assign_tree);
			new_tree->AddChild(new OMLTree(IDENT, temp.StringVal(), NULL, 0, 0));

			assign_tree = new_tree;
		}
	}

	outputs.push_back(assign_tree);

	return true;
}

bool oml_p_func(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() == 0)
		throw OML_Error("Illegal function definition");

	Currency func_name = inputs[0];

	OMLTree* func_tree = func_name.ConvertToTree();
		
	if (!func_tree)
		throw OML_Error("Illegal function definition");

	int numargs = (int)inputs.size() - 1;

	OMLTree* my_tree = new OMLTree(FUNC, "FUNC", NULL, 0, 1);

	my_tree->AddChild(func_tree);

	OMLTree* param_tree = new OMLTree(PARAM_LIST, "PARAM_LIST", NULL, 0, numargs);
	my_tree->AddChild(param_tree);

	for (int j = 0; j < numargs; ++j)
	{
		Currency arg = inputs[j + 1];
		arg.ConvertToTree();

		if (!arg.IsTree())		
			throw OML_Error("Illegal function definition");

		param_tree->AddChild(arg.Tree());
	}

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_add(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(PLUS, "+", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_statement(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if ((inputs.size() != 2) && (inputs.size() != 3))
		throw OML_Error("Illegal function definition");

	Currency parent = inputs[0];
	Currency stmt = inputs[1];

	int semicolon = 0;
	
	if (inputs.size() == 3)
	{
		Currency semic = inputs[2];

		if (semic.IsScalar())
		{
			if (semic.Scalar())
				semicolon = 1;
		}
		else
		{
			throw OML_Error("Illegal function definition");
		}
	}

	if (!parent.IsTree())
		throw OML_Error("Illegal function definition");

	if (!stmt.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* cur_tree = parent.Tree();
	OMLTree* stmts = NULL;
	
	if (parent.Tree()->GetType() == FUNC_DEF)
		stmts = cur_tree->GetChild(3);
	else if (parent.Tree()->GetType() == IF)
		stmts = cur_tree->GetChild(1);
	else if (parent.Tree()->GetType() == ELSEIF)
		stmts = cur_tree->GetChild(1);
	else if (parent.Tree()->GetType() == ELSE)
		stmts = cur_tree->GetChild(0);
	else if (parent.Tree()->GetType() == FOR)
		stmts = cur_tree->GetChild(2);
	else if (parent.Tree()->GetType() == WHILE)
		stmts = cur_tree->GetChild(1);
	else if (parent.Tree()->GetType() == CASE)
		stmts = cur_tree->GetChild(1);
	else if (parent.Tree()->GetType() == OTHERWISE)
		stmts = cur_tree->GetChild(0);
	else 
		stmts = parent.Tree();

	int num_children = 1;

	if (semicolon)
		num_children = 2;

	if (stmt.Tree()->GetType() == FUNC_DEF) 
	{
		stmts->AddChild(stmt.Tree());
	}
	else if (stmt.Tree()->GetType() == CONDITIONAL)
	{
		stmts->AddChild(stmt.Tree());
	}
	else if (stmt.Tree()->GetType() == SWITCH)
	{
		stmts->AddChild(stmt.Tree());
	}
	else if (stmt.Tree()->GetType() == FOR)
	{
		stmts->AddChild(stmt.Tree());
	}
	else if (stmt.Tree()->GetType() == WHILE)
	{
		stmts->AddChild(stmt.Tree());
	}
	else if (stmt.Tree()->GetType() == TRY)
	{
		stmts->AddChild(stmt.Tree());
	}
	else
	{
		OMLTree* stmt_tree = new OMLTree(STMT, "STMT", NULL, 0, num_children);
		stmt_tree->AddChild(stmt.Tree());

		if (num_children == 2)
			stmt_tree->AddChild(new OMLTree(SEMIC, ";", NULL, 0, 0));

		stmts->AddChild(stmt_tree);
	}

	return true;
}

bool oml_p_handle(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal handle definition");

	Currency handle_def = inputs[0];

	if (!handle_def.IsString())
		throw OML_Error("Illegal handle definition");

	OMLTree* new_tree = new OMLTree(FUNC_HANDLE, "FUNC_HANDLE", NULL, 0, 2);
	new_tree->AddChild(new OMLTree(IDENT, handle_def.StringVal(), NULL, 0, 0));
	
	outputs.push_back(new_tree);

	return true;
}

bool oml_p_string(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal string definition");

	Currency string_def = inputs[0];

	if (!string_def.IsString())
		throw OML_Error("Illegal handle definition");

	std::string new_string = "'" + string_def.StringVal() + "'";

	OMLTree* new_tree = new OMLTree(HML_STRING, "HML_STRING", NULL, 0, 1);
	new_tree->AddChild(new OMLTree(IDENT, new_string, NULL, 0, 0));

	outputs.push_back(new_tree);

	return true;
}

bool oml_p_equal(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(EQUAL, "==", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
	return true;
}

bool oml_p_if(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal if statement");

	Currency if_test = inputs[0];

	if (!if_test.IsTree())
		throw OML_Error("Illegal if statement");

	OMLTree* my_tree = new OMLTree(IF, "if", NULL, 0, 2);
	my_tree->AddChild(if_test.Tree());
	my_tree->AddChild(new OMLTree(STATEMENT_LIST, "STATEMENT_LIST", NULL, 0, 0));

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_elseif(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal elseif statement");

	Currency elseif_test = inputs[0];

	if (!elseif_test.IsTree())
		throw OML_Error("Illegal elseif statement");

	OMLTree* my_tree = new OMLTree(IF, "IF", NULL, 0, 2);
	my_tree->AddChild(elseif_test.Tree());
	my_tree->AddChild(new OMLTree(STATEMENT_LIST, "STATEMENT_LIST", NULL, 0, 0));

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_else(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 0)
		throw OML_Error("Illegal elseif statement");

	OMLTree* my_tree = new OMLTree(ELSE, "ELSE", NULL, 0, 2);
	my_tree->AddChild(new OMLTree(STATEMENT_LIST, "STATEMENT_LIST", NULL, 0, 0));

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_conditional(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	OMLTree* my_tree = new OMLTree(CONDITIONAL, "CONDITIONAL", NULL, 0, 2);

	for (int j = 0; j < inputs.size(); ++j)
	{
		Currency branch = inputs[j];

		if (!branch.IsTree())
			throw OML_Error("Illegal conditional branch");

		my_tree->AddChild(branch.Tree());
	}

	outputs.push_back(my_tree);

	return true;
}
bool oml_p_divide(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(DIV, "/", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_bitand(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(AND, "&", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_entrypow(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(DOTPOW, ".^", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_entrydivide(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(EDIV, "./", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_entryleftdivide(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(ELDIV, ".\\", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_greaterequal(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(GEQ, ">=", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_greaterthan(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(GTHAN, ">", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_logand(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(LAND, "&&", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_leftdivide(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(LDIV, "\\", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_lessequal(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(LEQ, "<=", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_logor(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(LOR, "||", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_lessthan(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(LTHAN, "<", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_subtract(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(MINUS, "-", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_notequal(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(NEQUAL, "~=", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_or(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(OR, "|", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_pow(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(POW, "^", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_multiply(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal function definition");

	Currency lhs = inputs[0];
	Currency rhs = inputs[1];

	lhs.ConvertToTree();
	rhs.ConvertToTree();

	if (!lhs.IsTree())
		throw OML_Error("Illegal function definition");

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(TIMES, "*", NULL, 0, 2);
	my_tree->AddChild(lhs.Tree());
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_negate(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal function definition");

	Currency rhs = inputs[0];

	rhs.ConvertToTree();

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(NEGATE, "~", NULL, 0, 1);
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_uminus(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal function definition");

	Currency rhs = inputs[0];

	rhs.ConvertToTree();

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(UMINUS, "UMINUS", NULL, 0, 1);
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_transpose(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal function definition");

	Currency rhs = inputs[0];

	rhs.ConvertToTree();

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(TRANSP, "TRANSP", NULL, 0, 1);
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_ctranspose(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal function definition");

	Currency rhs = inputs[0];

	rhs.ConvertToTree();

	if (!rhs.IsTree())
		throw OML_Error("Illegal function definition");

	OMLTree* my_tree = new OMLTree(CTRANSP, "CTRANSP", NULL, 0, 1);
	my_tree->AddChild(rhs.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_global(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal function definition");

	Currency rhs = inputs[0];

	if (!rhs.IsCellArray())
		throw OML_Error("Illegal function definition");

	HML_CELLARRAY* cells = rhs.CellArray();

	int num_children = cells->Size();

	OMLTree* my_tree = new OMLTree(GLOBAL, "global", NULL, 0, num_children);

	for (int j = 0; j < num_children; j++)
	{
		Currency temp = (*cells)(j);

		if (!temp.IsString())
			throw OML_Error("Illegal global definition");

		my_tree->AddChild(new OMLTree(IDENT, temp.StringVal(), NULL, 0, 0));
	}

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_persistent(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal function definition");

	Currency rhs = inputs[0];

	if (!rhs.IsCellArray())
		throw OML_Error("Illegal function definition");

	HML_CELLARRAY* cells = rhs.CellArray();

	int num_children = cells->Size();

	OMLTree* my_tree = new OMLTree(PERSISTENT, "persistent", NULL, 0, num_children);

	for (int j = 0; j < num_children; j++)
	{
		Currency temp = (*cells)(j);

		if (!temp.IsString())
			throw OML_Error("Illegal global definition");

		my_tree->AddChild(new OMLTree(IDENT, temp.StringVal(), NULL, 0, 0));
	}

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_range(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if ((inputs.size() != 2) && (inputs.size() != 3))
		throw OML_Error("Illegal range definition");

	Currency cur1 = inputs[0];
	Currency cur2 = inputs[1];

	cur1.ConvertToTree();
	cur2.ConvertToTree();

	if (!cur1.IsTree())
		throw OML_Error("Illegal range definition");

	if (!cur2.IsTree())
		throw OML_Error("Illegal range definition");

	OMLTree* my_tree = new OMLTree(COLON, ":", NULL, 0, 2);
	my_tree->AddChild(cur1.Tree());
	my_tree->AddChild(cur2.Tree());

	if (inputs.size() == 3)
	{
		Currency cur3 = inputs[2];

		cur3.ConvertToTree();

		if (!cur3.IsTree())
			throw OML_Error("Illegal range definition");

		OMLTree* new_tree = new OMLTree(COLON, ":", NULL, 0, 2); 

		new_tree->AddChild(my_tree);
		new_tree->AddChild(cur3.Tree());
		
		my_tree = new_tree;
	}

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_forloop(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal range definition");

	Currency cur1 = inputs[0];
	Currency cur2 = inputs[1];

	cur1.ConvertToTree();
	cur2.ConvertToTree();

	if (!cur1.IsTree())
		throw OML_Error("Illegal for definition");

	if (!cur2.IsTree())
		throw OML_Error("Illegal for definition");

	OMLTree* my_tree = new OMLTree(FOR, "for", NULL, 0, 3);
	my_tree->AddChild(cur1.Tree());
	my_tree->AddChild(cur2.Tree());
	my_tree->AddChild(new OMLTree(STATEMENT_LIST, "STATEMENT_LIST", NULL, 0, 0));

	outputs.push_back(my_tree);

	return true;
}
bool oml_p_trycatch(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 3)
		throw OML_Error("Illegal trycatch definition");

	Currency cur1 = inputs[0];
	Currency cur2 = inputs[1];
	Currency cur3 = inputs[2];

	cur3.ConvertToTree();

	if (!cur1.IsTree())
		throw OML_Error("Illegal trycatch definition");

	if (!cur2.IsTree())
		throw OML_Error("Illegal trycatch definition");

	if (!cur3.IsTree())
		throw OML_Error("Illegal trycatch definition");

	OMLTree* my_tree = new OMLTree(TRY, "try", NULL, 0, 3);
	my_tree->AddChild(cur1.Tree());
	my_tree->AddChild(cur2.Tree());
	my_tree->AddChild(cur3.Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_try(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 0)
		throw OML_Error("Illegal try statement");

	OMLTree* my_tree = new OMLTree(STATEMENT_LIST, "STATEMENT_LIST", NULL, 0, 0);

	outputs.push_back(my_tree);

	return true;
}
bool oml_p_catch(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 0)
		throw OML_Error("Illegal try statement");

	OMLTree* my_tree = new OMLTree(STATEMENT_LIST, "STATEMENT_LIST", NULL, 0, 0);

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_cellvalue(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 2)
		throw OML_Error("Illegal cell dereference");

	Currency cur1 = inputs[0];
	Currency cur2 = inputs[1];
	
	cur1.ConvertToTree();

	if (!cur1.IsTree())
		throw OML_Error("Illegal cell dereference");

	if (!cur2.IsCellArray())
		throw OML_Error("Illegal cell dereference");

	HML_CELLARRAY* cells = cur2.CellArray();

	OMLTree* my_tree = new OMLTree(CELL_VAL, "CELL_VAL", NULL, 0, 2);

	my_tree->AddChild(cur1.Tree());
	
	OMLTree* param_tree = new OMLTree(PARAM_LIST, "PARAM_LIST", NULL, 0, cells->Size());

	for (int j = 0; j < cells->Size(); ++j)
	{
		Currency temp = (*cells)(j);
		temp.ConvertToTree();

		if (!temp.IsTree())
			throw OML_Error("Illegal cell dereference");

		param_tree->AddChild(temp.Tree());
	}

	my_tree->AddChild(param_tree);

	outputs.push_back(my_tree);

	return true;
}
bool oml_p_whileloop(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal while definition");

	Currency cur1 = inputs[0];

	cur1.ConvertToTree();

	if (!cur1.IsTree())
		throw OML_Error("Illegal for definition");

	OMLTree* my_tree = new OMLTree(WHILE, "while", NULL, 0, 3);
	my_tree->AddChild(cur1.Tree());
	my_tree->AddChild(new OMLTree(STATEMENT_LIST, "STATEMENT_LIST", NULL, 0, 0));

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_case(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 1)
		throw OML_Error("Illegal case statement");

	Currency case_test = inputs[0];

	case_test.ConvertToTree();

	if (!case_test.IsTree())
		throw OML_Error("Illegal case statement");

	OMLTree* my_tree = new OMLTree(CASE, "case", NULL, 0, 2);
	my_tree->AddChild(case_test.Tree());
	my_tree->AddChild(new OMLTree(STATEMENT_LIST, "STATEMENT_LIST", NULL, 0, 0));

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_otherwise(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 0)
		throw OML_Error("Illegal otherwise statement");

	OMLTree* my_tree = new OMLTree(OTHERWISE, "otherwise", NULL, 0, 2);
	my_tree->AddChild(new OMLTree(STATEMENT_LIST, "STATEMENT_LIST", NULL, 0, 0));

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_switch(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	int num_inputs = (int)inputs.size();

	OMLTree* my_tree = new OMLTree(SWITCH, "switch", NULL, 0, num_inputs);

	Currency case_test = inputs[0];

	case_test.ConvertToTree();

	if (!case_test.IsTree())
		throw OML_Error("Illegal switch statement");

	my_tree->AddChild(case_test.Tree());

	for (int j = 1; j < num_inputs; ++j)
		my_tree->AddChild(inputs[j].Tree());

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_matrix(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	int num_inputs = (int)inputs.size();

	OMLTree* my_tree = new OMLTree(MATRIX, "MATRIX", NULL, 0, num_inputs);

	for (int j = 0; j < num_inputs; ++j)
	{
		Currency cur = inputs[j];
		cur.ConvertToTree();

		if (!cur.IsTree())
			throw OML_Error("Illegal vector");

		// might want to check that all the trees are vectors ...

		my_tree->AddChild(cur.Tree());
	}

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_vector(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	int num_inputs = (int)inputs.size();

	OMLTree* my_tree = new OMLTree(VECTOR, "VECTOR", NULL, 0, num_inputs);

	for (int j = 0; j < num_inputs; ++j)
	{
		Currency cur = inputs[j];
		cur.ConvertToTree();

		if (!cur.IsTree())
			throw OML_Error("Illegal vector");

		my_tree->AddChild(cur.Tree());
	}

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_cellarray(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	int num_inputs = (int)inputs.size();

	OMLTree* my_tree = new OMLTree(CELL_ARRAY, "CELL_ARRAY", NULL, 0, num_inputs);

	for (int j = 0; j < num_inputs; ++j)
	{
		Currency cur = inputs[j];
		cur.ConvertToTree();

		if (!cur.IsTree())
			throw OML_Error("Illegal vector");

		// might want to check that all the trees are vectors ...

		my_tree->AddChild(cur.Tree());
	}

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_colon(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 0)
		throw OML_Error("Illegal use of colon");

	OMLTree* my_tree = new OMLTree(COLON, ":", NULL, 0, 2);

	outputs.push_back(my_tree);

	return true;
}

bool oml_p_return(EvaluatorInterface, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() != 0)
		throw OML_Error("Illegal use of return");

	OMLTree* my_tree = new OMLTree(RETURN, "return", NULL, 0, 2);

	outputs.push_back(my_tree);

	return true;
}
bool oml_getprofiledata(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	outputs.push_back(eval.GetProfileData());
	return true;
}

bool oml_clearprofiledata(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	eval.ClearProfileData();
	return true;
}

bool oml_profile(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (inputs.size() == 1)
	{
		Currency temp = inputs[0];

		if (temp.IsLogical() && (temp.Scalar() == 1.0))
			eval.Profile(true);
		else
			eval.Profile(false);
	}

	return true;
}

bool oml_properties(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (!inputs.size())
	{
		throw OML_Error(OML_ERR_INPUT_EMPTY);
	}

	Currency cur = inputs[0];

	std::string classname;

	if (cur.IsString())
	{
		classname = cur.StringVal();
	}
	else if (cur.IsObject())
	{
		classname = cur.GetClassname();
	}

	std::vector<std::string> properties = eval.GetProperties(classname);

	outputs.push_back(properties);

	return true;
}

bool oml_methods(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
	if (!inputs.size())
	{
		throw OML_Error(OML_ERR_INPUT_EMPTY);
	}

	Currency cur = inputs[0];

	std::string classname;

	if (cur.IsString())
	{
		classname = cur.StringVal();
	}
	else if (cur.IsObject())
	{
		classname = cur.GetClassname();
	}

	std::vector<std::string> methods = eval.GetMethods(classname);

	outputs.push_back(methods);

	return true;
}

bool oml_inputname(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs)
{
    if (inputs.size() != 1)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency val = inputs[0];
    if (!val.IsPositiveInteger())
        throw OML_Error(OML_ERR_POSINTEGER);

    int idx = (int)val.Scalar();

    outputs.push_back(eval.GetFunctionArgumentName(idx));

    return true;
}