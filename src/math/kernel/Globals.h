/**
* @file Globals.h
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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
#ifndef _MathCore_Globals_h
#define _MathCore_Globals_h

#include <string>
#include <MathCoreExports.h>

using std::string;

//////////////////////////////////////////////////////////////////////
// Function pointer default

typedef void (*hMathFuncPtr) ();

//////////////////////////////////////////////////////////////////////
// Constants

#ifdef OS_WIN
   #define NULL 0
#endif

/*
Cephes Math Library, Release 2.8: June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier

The Cephes Library is public domain, available from http://www.netlib.org/cephes.
These constants were extracted and reformatted by Brian Shier.
*/

static const double PI = 3.14159265358979323846;
static const double MACHEP  = 1.11022302462515654042E-16;   // 2**-53
static const double MACHEP2 = 2.22044604925031308084E-16;   // 2**-52
static const double MAXNUM  = 1.7976931348623158E+308;      // 2**1024*(1-MACHEP)
static const double MAXLOG  = 7.09782712893383996732E+2;    // log(MAXNUM)
static const double MINLOG  = -7.08396418532264106224E+2;   // log(2**-1022)
static const double MAXGAM  = 171.624376956302725;

//////////////////////////////////////////////////////////////////////
// Warning messages
#define HW_MATH_MSG_SUCCESS         "Status: Success"
// general warnings
#define HW_MATH_MSG_NONNONNEGINT_W  "Warning: invalid value; must be nonnegative integer(s)"
// matrix warnings
#define HW_MATH_MSG_SINGMATRIX      "Warning: singular matrix (to within machine precision)"
#define HW_MATH_MSG_SINGMATRIXDIV   "Warning: singular matrix divisor (to within machine precision)"
#define HW_MATH_MSG_MTXNOTFULLRANK_W "Warning: rank deficient matrix; matrix does not have full rank"
#define HW_MATH_MSG_MATRIXDEPENDCOL "Warning: rank deficient matrix; columns are not linearly independent"
// algorithm warnings
#define HW_MATH_MSG_MAXITERATE      "Warning: maximum iterations exceeded, returned last estimate"
#define HW_MATH_MSG_MAXFUNCEVAL     "Warning: maximum number of function evaluations exceeded, returned last estimate"
#define HW_MATH_MSG_NOTCONVERGE_W   "Warning: algorithm did not converge, returned best estimate"
#define HW_MATH_MSG_LOCALMIN        "Warning: algorithm converged at local minimum"
#define HW_MATH_MSG_NOSOLUTION      "Warning: no solution, returned best fit"
#define HW_MATH_MSG_TOOFEWPOINTS_W  "Warning: does not contain enough data points, no unique solution"
#define HW_MATH_MSG_NOPEAKS         "Warning: invalid input; has no peaks"
#define HW_MATH_MSG_GAMAXITERS_W    "Warning: maximum iterations will be capped at 1000"
#define HW_MATH_MSG_CONVIOLATEBEST  "Warning: the best reported design violates the constraint"
#define HW_MATH_MSG_POORPOLYFIT     "Warning: if poor fit quality, consult the help document"
// signal processing warnings
#define HW_MATH_MSG_FILTERSIMPLOW   "Warning: syntax simplifies to low pass filter"
#define HW_MATH_MSG_FILTERSIMPHIGH  "Warning: syntax simplifies to high pass filter"
#define HW_MATH_MSG_FILTERSPEC_W    "Warning: filter specification is too extreme, results cannot be trusted"
#define HW_MATH_MSG_NYQUIST         "Warning: data exceeds Nyquist frequency"
#define HW_MATH_MSG_UNEQUALTIMESAMP "Warning: samples are not equally spaced"
#define HW_MATH_MSG_FILTERCFC       "Warning: invalid frequency; filter should not be used for CFC > 180"
// statistics warnings
#define HW_MATH_MSG_STATTESTALPHA_W "Warning: invalid alpha; should be in the interval (0.0, 0.5)" 
#define HW_MATH_MSG_SMALLVARIANCE   "Warning: data has minor variation; reported correlation may not be meaningful"

//////////////////////////////////////////////////////////////////////
// Error messages
// general errors
#define HW_MATH_MSG_NULLPOINTER     "Error: null pointer"
#define HW_MATH_MSG_INVALIDINPUT    "Error: invalid input"
#define HW_MATH_MSG_NOTSTRING       "Error: invalid type; must be string"
#define HW_MATH_MSG_INVALIDINTERVAL "Error: invalid interval"
#define HW_MATH_MSG_NOTALLOWED      "Error: operation not allowed"
#define HW_MATH_MSG_NOTIMPLEMENT    "Error: operation not implemented"
#define HW_MATH_MSG_COMPLEX         "Error: invalid type; must be real"
#define HW_MATH_MSG_NEEDCOMPLEX     "Error: invalid type; must be complex"
#define HW_MATH_MSG_COMPLEXSUPPORT  "Error: unsupported operation for complex data; must be real"
#define HW_MATH_MSG_ZERO            "Error: invalid value; must not be zero"
//#define HW_MATH_MSG_NONZERO         "Error: invalid value; must be zero"
#define HW_MATH_MSG_DIVIDEZERO      "Error: division by zero"
#define HW_MATH_MSG_NONFINITEDATA   "Error: invalid value; must be finite"
#define HW_MATH_MSG_NEGATIVE        "Error: invalid value; must be nonnegative"
#define HW_MATH_MSG_NONPOSITIVE     "Error: invalid value; must be positive"
#define HW_MATH_MSG_NONINTEGER      "Error: invalid value; must be integer"
#define HW_MATH_MSG_NONPOSINT       "Error: invalid value; must be positive integer"
#define HW_MATH_MSG_NONNONNEGINT    "Error: invalid value; must be nonnegative integer(s)"
#define HW_MATH_MSG_NONINCREASE     "Error: invalid data; must be strictly increasing"
#define HW_MATH_MSG_MINMAXVALUES    "Error: invalid values; maximum must be >= minimum"
#define HW_MATH_MSG_BADRANGE        "Error: invalid data; is out of range"
#define HW_MATH_MSG_ZERORANGE       "Error: invalid data; has a range of zero"
#define HW_MATH_MSG_BADFILE         "Error: file missing or could not be read;"
// user function errors
#define HW_MATH_MSG_USERFUNCFAIL    "Error: user function failed; could not be evaluated"
#define HW_MATH_MSG_USERFUNCMATRIX  "Error: user function failed; must return a real matrix"
#define HW_MATH_MSG_USERFUNCREAL    "Error: user function failed; must return real data"
#define HW_MATH_MSG_USERFUNCREALNUM "Error: user function failed; must return a real number"
#define HW_MATH_MSG_USERFUNCREALMAT "Error: user function failed; must return a real matrix"
#define HW_MATH_MSG_USERFUNCSIZE    "Error: user function failed; has mismatched input/output dimensions"
#define HW_MATH_MSG_USERFUNCSIZE2   "Error: user function failed; has invalid output dimensions"
// matrix errors
#define HW_MATH_MSG_OUTOFMEMORY     "Error: available memory exhausted"
#define HW_MATH_MSG_ALLOCFAILED     "Error: allocation failure"
#define HW_MATH_MSG_NONEMPTYMATRIX  "Error: invalid matrix; must be []"
#define HW_MATH_MSG_EMPTYMATRIX     "Error: invalid matrix; must contain data"
#define HW_MATH_MSG_MATRIXRESIZE    "Error: invalid matrix; cannot be resized"
#define HW_MATH_MSG_MATRIXRESHAPE1  "Error: invalid request; only one dimension can be inferred"
#define HW_MATH_MSG_MATRIXRESHAPE2  "Error: invalid dimensions; number of elements must be preserved"
#define HW_MATH_MSG_PERMVEC1        "Error: invalid permutation vector; contains too few elements"
#define HW_MATH_MSG_PERMVEC2        "Error: invalid permutation vector; elements must be positive"
#define HW_MATH_MSG_PERMVEC3        "Error: invalid permutation vector; all dimensions must be included"
#define HW_MATH_MSG_PERMVEC4        "Error: invalid permutation vector; duplicate elements are not allowed"
#define HW_MATH_MSG_ARRAYDIM        "Error: invalid dimension; must be nonnegative integer"
#define HW_MATH_MSG_ARRAYSIZE       "Error: incompatible matrices; dimensions must be equal"
#define HW_MATH_MSG_ARRAYTOOFEWROWS "Error: invalid matrix; cannot have fewer rows than columns"
#define HW_MATH_MSG_ARRAYTYPE       "Error: incompatible matrices; data types must be the same"
#define HW_MATH_MSG_COLUMNDIM       "Error: incompatible matrices; column dimensions must be equal"
#define HW_MATH_MSG_COLUMNVEC       "Error: invalid matrix; must be column vector"
#define HW_MATH_MSG_INVALIDINDEX    "Error: invalid index; must be within allowed bounds"
#define HW_MATH_MSG_DATANOTFOUND    "Error: invalid data; not found"
#define HW_MATH_MSG_MTXNOTFULLRANK_E "Error: invalid matrix; must have full rank"
#define HW_MATH_MSG_MTXNOTSQUARE    "Error: invalid matrix; must be square"
#define HW_MATH_MSG_MTXNOTSYM       "Error: invalid matrix; must be symmetric"
#define HW_MATH_MSG_MTXNOTSPD       "Error: invalid matrix; must be symmetric positive definite"
#define HW_MATH_MSG_DECOMPFAIL      "Error: invalid matrix; could not be decomposed"
#define HW_MATH_MSG_DIAGFAIL        "Error: invalid matrix; could not be diagonalized"
#define HW_MATH_MSG_VECTOR12        "Error: invalid vector; length must be 1 or 2"
#define HW_MATH_MSG_VECTOR2         "Error: invalid vector; must have a length 2"
#define HW_MATH_MSG_VECTOR          "Error: invalid matrix; must be non-empty vector"
#define HW_MATH_MSG_NOVECTOR        "Error: invalid vector; must be a matrix"
#define HW_MATH_MSG_KLUD            "Error: invalid index; must be a positive integer < #columns"
// slice errors
#define HW_MATH_MSG_SLICE_INDEX     "Error: invalid slice; dimension is inconsistent or out of bounds"
#define HW_MATH_MSG_SLICE_NUMCOLON  "Error: invalid syntax; can have only one non-colon"
// algorithm errors
#define HW_MATH_MSG_NOTCONVERGE_E   "Error: algorithm did not converge"
#define HW_MATH_MSG_NOLOCALMIN      "Error: no local minimum near the initial condition"
#define HW_MATH_MSG_UNDERDETSYS_E   "Error: invalid system; must have at least as many equations as parameters"
#define HW_MATH_MSG_UNDERDETSYS_P   "Error: invalid system; must have at least as many points as parameters"
#define HW_MATH_MSG_TOOFEWPOINTS_E  "Error: invalid input; does not contain enough data points"
#define HW_MATH_MSG_DISTINCTVALS3   "Error: invalid input; must contain at least 3 distinct values"
#define HW_MATH_MSG_TIGHTTOL        "Error: invalid tolerance; may be too small"
#define HW_MATH_MSG_GAMAXITERS_E    "Error: requested maximum iterations exceeds 1000"
#define HW_MATH_MSG_POPSIZE         "Error: invalid population size; must be 0 or greater than 2"
#define HW_MATH_MSG_GLOBALSEARCH    "Error: invalid global search; must be 1, 2 or 3"
#define HW_MATH_MSG_TOOFEWVEHICLES  "Error: invalid vehicles; not enough available"
#define HW_MATH_MSG_GAUSSLEGENDRE   "Error: invalid number of points; must be from 1 to 10"
#define HW_MATH_MSG_GAUSSLOBATTO    "Error: invalid number of points; must be from 2 to 10"
#define HW_MATH_MSG_RK45STEPSIZE    "Error: invalid step size; unnecessarily small"
//#define HW_MATH_MSG_RK45WORKLOAD    "Error: the work load limit has been exceeded"
#define HW_MATH_MSG_QUADSTEPSIZE    "Error: minimum allowed step size has been reached"
#define HW_MATH_MSG_FUNCTIONCOUNT   "Error: maximum allowed number of function calls has been reached"
#define HW_MATH_MSG_NOUNIQUEROOT    "Error: invalid polynomial; does not have a unique root"
#define HW_MATH_MSG_INTERNALERROR   "Error: an internal error occurred"
// signal processing errors
#define HW_MATH_MSG_FILTERTYPE      "Error: invalid filter; must be a supported type"
#define HW_MATH_MSG_FILTERORDER     "Error: invalid filter; must be positive"
#define HW_MATH_MSG_FILTERORDERIRR  "Error: invalid filter; must be an integer from 1 to 14"
#define HW_MATH_MSG_FILTERORDERODD  "Error: invalid filter; requires even order"
#define HW_MATH_MSG_FILTERSPEC_E    "Error: invalid filter, specification is too extreme"
#define HW_MATH_MSG_FILTERFREQ_A    "Error: invalid filter; must be positive"
#define HW_MATH_MSG_FILTERFREQ_D    "Error: invalid filter; must be between 0 and 1"
#define HW_MATH_MSG_FILTERFREQS_EQ  "Error: invalid filter frequencies; cannot be equal"
#define HW_MATH_MSG_FILTERBANDCONF  "Error: invalid filter frequencies; bands conflict"
#define HW_MATH_MSG_FILTERCLASS     "Error: invalid filter; class must be positive"
#define HW_MATH_MSG_FREQCLASS       "Error: invalid sampling frequency; incorrect for this filter class"
#define HW_MATH_MSG_FILTERRIPPLE    "Error: invalid filter; stop band attenuation must exceed pass band"
#define HW_MATH_MSG_FILTERDENZERO   "Error: invalid filter first coefficient; cannot be zero"
#define HW_MATH_MSG_FILTFILTCOEFS   "Error: invalid filter coefficients; must have length greater that 1"
#define HW_MATH_MSG_FILTFILTDATA    "Error: invalid data vector; must be > 3 X filter order"
#define HW_MATH_MSG_FIRWINDOW       "Error: invalid window; length must be filter order + 1"
#define HW_MATH_MSG_DB_SIGN         "Error: invalid dB value; must be positive"
#define HW_MATH_MSG_OVERLAPPOINTS   "Error: invalid overlap; must be less than the block or window size"
#define HW_MATH_MSG_FTBLOCKSIZE     "Error: invalid block size; cannot exceed number of samples"
#define HW_MATH_MSG_FFTSIZE         "Error: invalid fft size; cannot be less than the block or window size"
#define HW_MATH_MSG_PADARG          "Error: invalid parameter; must be 'pad', the FFT size, or omitted"
#define HW_MATH_MSG_RESAMPOFFSET    "Error: invalid offset; must be positive integer less than period length"
// statistics errors
#define HW_MATH_MSG_INVALIDPROB     "Error: invalid probability; must be in the interval (0,1)" 
#define HW_MATH_MSG_NORMTESTPNTS    "Error: invalid data; must have at least 20 data points" 
#define HW_MATH_MSG_STATTESTALPHA_E "Error: invalid alpha; must be in the interval (0.0, 0.5)" 
#define HW_MATH_MSG_DOELEVEL        "Error: invalid DOE; must have more than one level"
// CVODE errors
#define HW_MATH_MSG_CVODE_WORKLOAD  "Error: ODE solver work load exceeded"
#define HW_MATH_MSG_CVODE_ACCURACY  "Error: ODE solver error; accuracy requirements too high"
#define HW_MATH_MSG_CVODE_ERR       "Error: ODE solver error; internal errors, possible step size related"
#define HW_MATH_MSG_CVODE_CONV      "Error: ODE solver error; convergence errors, possible step size related"
// IDA errors
#define HW_MATH_MSG_IDA_WORKLOAD    "Error: DAE solver work load exceeded"
#define HW_MATH_MSG_IDA_ACCURACY    "Error: DAE solver error; accuracy requirements too high"
#define HW_MATH_MSG_IDA_ERR         "Error: DAE solver error; internal errors, possibly step size related"
#define HW_MATH_MSG_IDA_CONV        "Error: DAE solver error; convergence errors, possibly step size related"
// unknown error
#define HW_MATH_MSG_UNKNOWN         "Error: contact Altair customer support"

enum hwMathMsgCode
{
    // the enum begins explicitly at 0
    // only positive values should be used
    // negative values are reserved for other modules
    // that may produce their own int error flags

    HW_MATH_ERR_NONE = 0,
    // general warnings
    HW_MATH_WARN_NONNONNEGINT,
    // matrix warnings
    HW_MATH_WARN_SINGMATRIX,
    HW_MATH_WARN_SINGMATRIXDIV,
    HW_MATH_WARN_MTXNOTFULLRANK,
    HW_MATH_WARN_MATRIXDEPENDCOL,
    // algorithm warnings
    HW_MATH_WARN_MAXITERATE,
    HW_MATH_WARN_MAXFUNCEVAL,
    HW_MATH_WARN_NOTCONVERGE,
    HW_MATH_WARN_LOCALMIN,
    HW_MATH_WARN_NOSOLUTION,
    HW_MATH_WARN_NOUNIQUESOL,
    HW_MATH_WARN_TOOFEWPOINTS,
    HW_MATH_WARN_NOPEAKS,
    HW_MATH_WARN_GAMAXITERS,
    HW_MATH_WARN_CONVIOLATEBEST,
    HW_MATH_WARN_POORPOLYFIT,
    // signal processing warnings
    HW_MATH_WARN_FILTERSIMPLOW,
    HW_MATH_WARN_FILTERSIMPHIGH,
    HW_MATH_WARN_FILTERSPEC_W,
    HW_MATH_WARN_NYQUIST,
    HW_MATH_WARN_UNEQUALTIMESAMP,
    HW_MATH_WARN_FILTERCFC,
    // statistics warnings
    HW_MATH_WARN_STATTESTALPHA,
    HW_MATH_WARN_SMALLVARIANCE,
    // end of warning list
    HW_MATH_WARN_EOL,
    // general errors
    HW_MATH_ERR_NULLPOINTER,
    HW_MATH_ERR_INVALIDINPUT,
    HW_MATH_ERR_NOTSTRING,
    HW_MATH_ERR_INVALIDINTERVAL,
    HW_MATH_ERR_NOTALLOWED,
    HW_MATH_ERR_NOTIMPLEMENT,
    HW_MATH_ERR_COMPLEX,
    HW_MATH_ERR_NEEDCOMPLEX,
    HW_MATH_ERR_COMPLEXSUPPORT,
    HW_MATH_ERR_ZERO,
    HW_MATH_ERR_NONZERO,
    HW_MATH_ERR_DIVIDEZERO,
    HW_MATH_ERR_NONFINITEDATA,
    HW_MATH_ERR_NEGATIVE,
    HW_MATH_ERR_NONPOSITIVE,
    HW_MATH_ERR_NONINTEGER,
    HW_MATH_ERR_NONPOSINT,
    HW_MATH_ERR_NONNONNEGINT,
    HW_MATH_ERR_NONINCREASE,
    HW_MATH_ERR_MINMAXVALUES,
    HW_MATH_ERR_BADRANGE,
    HW_MATH_ERR_ZERORANGE,
    HW_MATH_ERR_BADFILE,
    // user function errors
    HW_MATH_ERR_USERFUNCFAIL,
    HW_MATH_ERR_USERFUNCMATRIX,
    HW_MATH_ERR_USERFUNCREAL,
    HW_MATH_ERR_USERFUNCREALNUM,
    HW_MATH_ERR_USERFUNCREALMAT,
    HW_MATH_ERR_USERFUNCSIZE,
    HW_MATH_ERR_USERFUNCSIZE2,
    // matrix errors
    HW_MATH_ERR_OUTOFMEMORY,
    HW_MATH_ERR_ALLOCFAILED,
    HW_MATH_ERR_EMPTYMATRIX,
    HW_MATH_ERR_NONEMPTYMATRIX,
    HW_MATH_ERR_MATRIXRESIZE,
    HW_MATH_ERR_MATRIXRESHAPE1,
    HW_MATH_ERR_MATRIXRESHAPE2,
    HW_MATH_ERR_PERMVEC1,
    HW_MATH_ERR_PERMVEC2,
    HW_MATH_ERR_PERMVEC3,
    HW_MATH_ERR_PERMVEC4,
    HW_MATH_ERR_ARRAYDIM,
    HW_MATH_ERR_ARRAYSIZE,
    HW_MATH_ERR_ARRAYTOOFEWROWS,
    HW_MATH_ERR_ARRAYTYPE,
    HW_MATH_ERR_COLUMNDIM,
    HW_MATH_ERR_COLUMNVEC,
    HW_MATH_ERR_INVALIDINDEX,
    HW_MATH_ERR_DATANOTFOUND,
    HW_MATH_ERR_MTXNOTFULLRANK,
    HW_MATH_ERR_MTXNOTSQUARE,
    HW_MATH_ERR_MTXNOTSYM,
    HW_MATH_ERR_MTXNOTSPD,
    HW_MATH_ERR_DECOMPFAIL,
    HW_MATH_ERR_DIAGFAIL,
    HW_MATH_ERR_VECTOR12,
    HW_MATH_ERR_VECTOR2,
    HW_MATH_ERR_VECTOR,
    HW_MATH_ERR_NOVECTOR,
    HW_MATH_ERR_KLUD,
    // slice errors
    HW_MATH_ERR_SLICE_INDEX,
    HW_MATH_ERR_SLICE_NUMCOLON,
    // algorithm errors
    HW_MATH_ERR_NOTCONVERGE,
    HW_MATH_ERR_NOLOCALMIN,
    HW_MATH_ERR_UNDERDETSYS_E,
    HW_MATH_ERR_UNDERDETSYS_P,
    HW_MATH_ERR_TOOFEWPOINTS,
    HW_MATH_ERR_DISTINCTVALS3,
    HW_MATH_ERR_TIGHTTOL,
    HW_MATH_ERR_GAMAXITERS,
    HW_MATH_ERR_POPSIZE,
    HW_MATH_ERR_GLOBALSEARCH,
    HW_MATH_ERR_TOOFEWVEHICLES,
    HW_MATH_ERR_GAUSSLEGENDRE,
    HW_MATH_ERR_GAUSSLOBATTO,
    HW_MATH_ERR_RK45STEPSIZE,
    // HW_MATH_ERR_RK45WORKLOAD,
    HW_MATH_ERR_QUADSTEPSIZE,
    HW_MATH_ERR_FUNCTIONCOUNT,
    HW_MATH_ERR_NOUNIQUEROOT,
    HW_MATH_ERR_INTERNALERROR,
    // signal processing errors
    HW_MATH_ERR_FILTERTYPE,
    HW_MATH_ERR_FILTERORDER,
    HW_MATH_ERR_FILTERORDERIRR,
    HW_MATH_ERR_FILTERORDERODD,
    HW_MATH_ERR_FILTERSPEC_E,
    HW_MATH_ERR_FILTERFREQ_A,
    HW_MATH_ERR_FILTERFREQ_D,
    HW_MATH_ERR_FILTERFREQS_EQ,
    HW_MATH_ERR_FILTERBANDCONF,
    HW_MATH_ERR_FILTERCLASS,
    HW_MATH_ERR_FREQCLASS,
    HW_MATH_ERR_FILTERRIPPLE,
    HW_MATH_ERR_FILTERDENZERO,
    HW_MATH_ERR_FILTFILTCOEFS,
    HW_MATH_ERR_FILTFILTDATA,
    HW_MATH_ERR_FIRWINDOW,
    HW_MATH_ERR_DB_SIGN,
    HW_MATH_ERR_OVERLAPPOINTS,
    HW_MATH_ERR_FTBLOCKSIZE,
    HW_MATH_ERR_FFTSIZE,
    HW_MATH_ERR_PADARG,
    HW_MATH_ERR_RESAMPOFFSET,
    // statistics errors
    HW_MATH_ERR_INVALIDPROB,
    HW_MATH_ERR_NORMTESTPNTS,
    HW_MATH_ERR_STATTESTALPHA,
    HW_MATH_ERR_DOELEVEL,
    // CVODE errors
    HW_MATH_ERR_CVODE_WORKLOAD,
    HW_MATH_ERR_CVODE_ACCURACY,
    HW_MATH_ERR_CVODE_ERR,
    HW_MATH_ERR_CVODE_CONV,
    // IDA errors
    HW_MATH_ERR_IDA_WORKLOAD,
    HW_MATH_ERR_IDA_ACCURACY,
    HW_MATH_ERR_IDA_ERR,
    HW_MATH_ERR_IDA_CONV,
    // unknown
    HW_MATH_ERR_UNKNOWN
};

MATHCORE_DECLS std::string GetHMathErrMsg(hwMathMsgCode code);

#include <memory.h>

// memcpy wrapper
#ifndef OS_WIN
  #ifndef memcpy_s
    static void* memcpy_s_wrapper(void *s1, size_t dummy, const void *s2, size_t n)
    {
        return memcpy(s1, s2, n);
    }
    #define memcpy_s memcpy_s_wrapper
  #endif
#endif

#endif // _MathCore_Globals_h
