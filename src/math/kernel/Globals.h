/**
* @file Globals.h
* @date June 2007
* Copyright (C) 2007-2019 Altair Engineering, Inc.  
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
#ifndef _MathCore_Globals_h
#define _MathCore_Globals_h

#include <string>
#include "MathCoreExports.h"

using std::string;

//!
//! typedef hMathFuncPtr
//! Function pointer default
//!
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

static const double PI      = 3.14159265358979323846;
static const double MACHEP  = 1.11022302462515654042E-16;   // 2**-53
static const double MACHEP2 = 2.22044604925031308084E-16;   // 2**-52
static const double MAXNUM  = 1.7976931348623158E+308;      // 2**1024*(1-MACHEP)
static const double MAXLOG  = 7.09782712893383996732E+2;    // log(MAXNUM)
static const double MINLOG  = -7.08396418532264106224E+2;   // log(2**-1022)
static const double MAXGAM  = 171.624376956302725;

//!
//! \enum hwMathMsgCode
//! Enum begins explicitly at 0 and only positive values are used
//! Negative values are reserved for other modules
//!
enum hwMathMsgCode
{
    HW_MATH_ERR_NONE = 0,              // Start of enum, explicitly set to 0

    // General warnings
    HW_MATH_WARN_NONNONNEGINT,         // invalid value; must be nonnegative integer(s)
    // Matrix warnings
    HW_MATH_WARN_SINGMATRIX,           // singular matrix (to within machine precision)
    HW_MATH_WARN_SINGMATRIXDIV,        // singular matrix divisor (to within machine precision)
    HW_MATH_WARN_MTXNOTFULLRANK,       // rank deficient matrix; matrix does not have full rank
    HW_MATH_WARN_MATRIXDEPENDCOL,      // rank deficient matrix; columns are not linearly independent
    HW_MATH_WARN_NOSQRTM,              // singular matrix; the square root may not exist
    HW_MATH_WARN_MAXSQRTM,             // maximum number of allowed square roots has been exceeded
    HW_MATH_WARN_HANKEL,               // first element of arg2 ignored in favor of last element of arg1
    HW_MATH_WARN_TOEPLITZ,             // first element of arg2 ignored in favor of first element of arg1
    // Algorithm warnings
    HW_MATH_WARN_MAXITERATE,           // maximum iterations exceeded, returned last estimate
    HW_MATH_WARN_MAXFUNCEVAL,          // maximum number of function evaluations exceeded, returned last estimate
    HW_MATH_WARN_NOTCONVERGE,          // algorithm did not converge, returned best estimate
    HW_MATH_WARN_LOCALMIN,             // algorithm converged at local minimum
    HW_MATH_WARN_NOSOLUTION,           // no solution, returned best fit
    HW_MATH_WARN_NOUNIQUESOL,          // no unique solution, returned best fit
    HW_MATH_WARN_TOOFEWPOINTS,         // does not contain enough data points, no unique solution
    HW_MATH_WARN_NOPEAKS,              // invalid input; has no peaks
    HW_MATH_WARN_GAMAXITERS,           // maximum iterations will be capped at 1000
    HW_MATH_WARN_CONVIOLATEBEST,       // the best reported design violates the constraint
    HW_MATH_WARN_POORPOLYFIT,          // if poor fit quality, consult the help document
    HW_MATH_WARN_DOUBLEPRECLIMIT,      // double precision accuracy limits may have been exceeded
    // Signal processing warnings
    HW_MATH_WARN_FILTERORDERIRR,       // high filter order; numerical accuracy may be diminished
    HW_MATH_WARN_FILTERSIMPLOW,        // syntax simplifies to low pass filter
    HW_MATH_WARN_FILTERSIMPHIGH,       // syntax simplifies to high pass filter
    HW_MATH_WARN_FILTERSPEC_W,         // extreme filter specification; numerical accuracy may be diminished
    HW_MATH_WARN_NYQUIST,              // data exceeds Nyquist frequency
    HW_MATH_WARN_UNEQUALTIMESAMP,      // samples are not equally spaced
    HW_MATH_WARN_FILTERCFC,            // invalid frequency; filter should not be used for CFC > 180
    // Statistics warnings
    HW_MATH_WARN_STATTESTALPHA,        // invalid alpha; should be in the interval (0.0, 0.5)
    // CVODE warnings
    HW_MATH_WARN_CVODE_EVENT,          // ODE solver terminated by Events function
    // IDA warnings
    HW_MATH_WARN_IDA_EVENT,            // DAE solver terminated by Events function
    // Qhull warnings                 
    HW_MATH_WARN_QHULL_LEAK,           // Error: Qhull did not release all memory
    // end of warning list
    HW_MATH_WARN_EOL,                  // End of Warnings

    // Information codes
    HW_MATH_INFO_TOLFCONV,             // function converged to within tolFun specification
    HW_MATH_INFO_TOLFCONV_R,           // function converged to within relative tolFun specification
    HW_MATH_INFO_TOLXCONV,             // step size converged to within tolX specification
    HW_MATH_INFO_TOLXCONV_R,           // step size converged to within relative tolX specification
    HW_MATH_INFO_SMALLTRUST,           // aborted due to trust region becoming too small
    // end of information list
    HW_MATH_INFO_EOL,                  // End of information messages

    // General errors
    HW_MATH_ERR_NULLPOINTER,           // null pointer
    HW_MATH_ERR_INVALIDINPUT,          // invalid input
    HW_MATH_ERR_NOTSTRING,             // invalid type; must be string
    HW_MATH_ERR_INVALIDINTERVAL,       // invalid interval
    HW_MATH_ERR_NOTALLOWED,            // operation not allowed
    HW_MATH_ERR_NOTIMPLEMENT,          // operation not implemented
    HW_MATH_ERR_COMPLEX,               // invalid type; must be real
    HW_MATH_ERR_NEEDCOMPLEX,           // invalid type; must be complex
    HW_MATH_ERR_COMPLEXSUPPORT,        // unsupported operation for complex data; must be real
    HW_MATH_ERR_ZERO,                  // invalid value; must not be zero
    HW_MATH_ERR_DIVIDEZERO,            // division by zero
    HW_MATH_ERR_NONFINITEDATA,         // invalid value; must be finite
    HW_MATH_ERR_POSITIVE,              // invalid value; must be negative
    HW_MATH_ERR_NEGATIVE,              // invalid value; must be nonnegative
    HW_MATH_ERR_NONPOSITIVE,           // invalid value; must be positive
    HW_MATH_ERR_NONINTEGER,            // invalid value; must be integer
    HW_MATH_ERR_NONPOSINT,             // invalid value; must be positive integer
    HW_MATH_ERR_NONNONNEGINT,          // invalid value; must be nonnegative integer(s)
    HW_MATH_ERR_NONINCREASE,           // invalid data; must be strictly increasing
    HW_MATH_ERR_NONUNIQUE,             // invalid data; must contain all unique values
    HW_MATH_ERR_CONSECUTIVE3,          // invalid data; must contain at most 2 consecutive identical values
    HW_MATH_ERR_MINMAXVALUES,          // invalid values; maximum must be >= minimum
    HW_MATH_ERR_ARG1LTARG2,            // invalid values; arg1 must be >= arg2
    HW_MATH_ERR_BADRANGE,              // invalid data; value is out of range
    HW_MATH_ERR_ZERORANGE,             // invalid data; has a range of zero
    HW_MATH_ERR_BADFILE,               // file missing or could not be read
    // User function errors
    HW_MATH_ERR_USERFUNCFAIL,          // user function failed; could not be evaluated
    HW_MATH_ERR_USERFUNCMATRIX,        // user function failed; must return a real matrix
    HW_MATH_ERR_USERFUNCREAL,          // user function failed; must return real data
    HW_MATH_ERR_USERFUNCREALNUM,       // user function failed; must return a real number
    HW_MATH_ERR_USERFUNCREALMAT,       // user function failed; must return a real matrix
    HW_MATH_ERR_USERFUNCSIZE,          // user function failed; has mismatched input/output dimensions
    HW_MATH_ERR_USERFUNCSIZE2,         // user function failed; has invalid output dimensions
    // Matrix errors
    HW_MATH_ERR_OUTOFMEMORY,           // available memory exhausted
    HW_MATH_ERR_ALLOCFAILED,           // allocation failure
    HW_MATH_ERR_EMPTYMATRIX,           // invalid matrix; must contain data
    HW_MATH_ERR_NONEMPTYMATRIX,        // invalid matrix; must be []
    HW_MATH_ERR_MATRIXRESIZE,          // invalid matrix; cannot be resized
    HW_MATH_ERR_MATRIXRESHAPE1,        // invalid request; only one dimension can be inferred
    HW_MATH_ERR_MATRIXRESHAPE2,        // invalid dimensions; number of elements must be preserved
    HW_MATH_ERR_PERMVEC1,              // invalid permutation vector; contains too few elements
    HW_MATH_ERR_PERMVEC2,              // invalid permutation vector; elements must be positive
    HW_MATH_ERR_PERMVEC3,              // invalid permutation vector; all dimensions must be included
    HW_MATH_ERR_PERMVEC4,              // invalid permutation vector; duplicate elements are not allowed
    HW_MATH_ERR_ARRAYDIM,              // invalid dimension; must be nonnegative integer
    HW_MATH_ERR_ARRAYSIZE,             // incompatible matrices; dimensions must be consistent
    HW_MATH_ERR_ARRAYCOL3,             // invalid matrix; 3D data requires 3 columns"
    HW_MATH_ERR_ARRAYTOOFEWROWS,       // invalid matrix; cannot have fewer rows than columns
    HW_MATH_ERR_ARRAYTYPE,             // incompatible matrices; data types must be the same
    HW_MATH_ERR_COLUMNDIM,             // incompatible matrices; column dimensions must be equal
    HW_MATH_ERR_COLUMNVEC,             // invalid matrix; must be column vector
    HW_MATH_ERR_INVALIDINDEX,          // invalid index; must be within allowed bounds
    HW_MATH_ERR_DATANOTFOUND,          // invalid data; not found
    HW_MATH_ERR_SINGMATRIX,            // singular matrix (to within machine precision)
    HW_MATH_ERR_MTXNOTFULLRANK,        // invalid matrix; must have full rank
    HW_MATH_ERR_MTXNOTSQUARE,          // invalid matrix; must be square
    HW_MATH_ERR_MTXNOTSYM,             // invalid matrix; must be symmetric
    HW_MATH_ERR_MTXNOTSPD,             // invalid matrix; must be symmetric positive definite
    HW_MATH_ERR_DECOMPFAIL,            // invalid matrix; could not be decomposed
    HW_MATH_ERR_DIAGFAIL,              // invalid matrix; could not be diagonalized
    HW_MATH_ERR_VECTOR12,              // invalid vector; length must be 1 or 2
    HW_MATH_ERR_VECTOR2,               // invalid vector; must have a length 2
    HW_MATH_ERR_VECTOR2x,              // invalid vector; must have an even length
    HW_MATH_ERR_VECTOR,                // invalid matrix; must be non-empty vector
    HW_MATH_ERR_NOVECTOR,              // invalid vector; must be a matrix
    HW_MATH_ERR_KLUD,                  // invalid index; must be a positive integer < #columns
    // Slice errors
    HW_MATH_ERR_SLICE_INDEX,           // invalid slice; dimension is inconsistent or out of bounds
    HW_MATH_ERR_SLICE_NUMCOLON,        // invalid syntax; can have only one non-colon
    // Algorithm errors
    HW_MATH_ERR_NOTCONVERGE,           // algorithm did not converge
    HW_MATH_ERR_NOLOCALMIN,            // no local minimum near the initial condition
    HW_MATH_ERR_UNDERDETSYS_E,         // invalid system; must have at least as many equations as parameters
    HW_MATH_ERR_UNDERDETSYS_P,         // invalid system; must have at least as many points as parameters
    HW_MATH_ERR_TOOFEWPOINTS,          // invalid input; does not contain enough data points
    HW_MATH_ERR_TOOFEWPOINTS10,        // invalid input; must contain at least 10 data points
    HW_MATH_ERR_DISTINCTVALS3,         // invalid input; must contain at least 3 distinct values
    HW_MATH_ERR_TIGHTTOL,              // invalid tolerance; may be too small
    HW_MATH_ERR_GAMAXITERS,            // requested maximum iterations has been exceeded
    HW_MATH_ERR_POPSIZE,               // invalid population size; must be 0 or greater than 2
    HW_MATH_ERR_GLOBALSEARCH,          // invalid global search; must be 1, 2 or 3
    HW_MATH_ERR_TOOFEWVEHICLES,        // invalid vehicles; not enough available
    HW_MATH_ERR_GAUSSLEGENDRE,         // invalid number of points; must be from 1 to 10
    HW_MATH_ERR_GAUSSLOBATTO,          // invalid number of points; must be from 2 to 10
    HW_MATH_ERR_RK45STEPSIZE,          // invalid step size; unnecessarily small
    HW_MATH_ERR_QUADSTEPSIZE,          // minimum allowed step size has been reached
    HW_MATH_ERR_FUNCTIONCOUNT,         // maximum allowed number of function calls has been reached
    HW_MATH_ERR_NOUNIQUEROOT,          // invalid polynomial; does not have a unique root
    HW_MATH_ERR_NESTSUPPORT,           // nested calls are not supported
    HW_MATH_ERR_INTERNALERROR,         // an internal error occurred
    // Signal processing errors
    HW_MATH_ERR_FILTERTYPE,            // invalid filter; must be a supported type
    HW_MATH_ERR_FILTERORDER,           // invalid filter; must be positive
    HW_MATH_ERR_FILTERORDERODD,        // invalid filter; requires even order
    HW_MATH_ERR_FILTERSPEC_E,          // invalid filter, specification is too extreme
    HW_MATH_ERR_FILTERFREQ_A,          // invalid filter; must be positive
    HW_MATH_ERR_FILTERFREQ_D,          // invalid filter; must be between 0 and 1
    HW_MATH_ERR_FILTERFREQS_EQ,        // invalid filter frequencies; cannot be equal
    HW_MATH_ERR_FILTERBANDCONF,        // invalid filter frequencies; bands conflict
    HW_MATH_ERR_FILTERCLASS,           // invalid filter; class must be positive
    HW_MATH_ERR_FREQCLASS,             // invalid sampling frequency; incorrect for this filter class
    HW_MATH_ERR_FILTERRIPPLE,          // invalid filter; stop band attenuation must exceed pass band
    HW_MATH_ERR_FILTERDENZERO,         // invalid filter first coefficient; cannot be zero
    HW_MATH_ERR_FILTFILTCOEFS,         // invalid filter coefficients; must have length greater that 1
    HW_MATH_ERR_FILTFILTDATA,          // invalid data vector; must be > 3 X filter order
    HW_MATH_ERR_FIRWINDOW,             // invalid window; length must be filter order + 1
    HW_MATH_ERR_DB_SIGN,               // invalid dB value; must be positive
    HW_MATH_ERR_OVERLAPPOINTS,         // invalid overlap; must be less than the block or window size
    HW_MATH_ERR_FTBLOCKSIZE,           // invalid block size; cannot exceed number of samples
    HW_MATH_ERR_FFTSIZE,               // invalid fft size; cannot be less than the block or window size
    HW_MATH_ERR_PADARG,                // invalid parameter; must be 'pad', the FFT size, or omitted
    HW_MATH_ERR_RESAMPOFFSET,          // invalid offset; must be positive integer less than period length
    // Statistics errors
    HW_MATH_ERR_INVALIDPROB,           // invalid probability; must be in the interval (0,1)
    HW_MATH_ERR_NORMTESTPNTS,          // invalid data; must have at least 20 data points
    HW_MATH_ERR_STATTESTALPHA,         // invalid alpha; must be in the interval (0.0, 0.5)
    HW_MATH_ERR_DOELEVEL,              // invalid DOE; must have more than one level
    // CVODE errors
    HW_MATH_ERR_CVODE_WORKLOAD,        // ODE solver work load exceeded
    HW_MATH_ERR_CVODE_ACCURACY,        // ODE solver error; accuracy requirements too high
    HW_MATH_ERR_CVODE_ERR,             // ODE solver error; internal errors, possible step size related
    HW_MATH_ERR_CVODE_CONV,            // ODE solver error; convergence errors, possible step size related
    // IDA errors
    HW_MATH_ERR_IDA_WORKLOAD,          // DAE solver work load exceeded
    HW_MATH_ERR_IDA_ACCURACY,          // DAE solver error; accuracy requirements too high
    HW_MATH_ERR_IDA_ERR,               // DAE solver error; internal errors, possibly step size related
    HW_MATH_ERR_IDA_CONV,              // DAE solver error; convergence errors, possibly step size related
    // Qhull errors
    HW_MATH_ERR_QHULL,                 // Qhull function error
    HW_MATH_ERR_QHULL_PNTS,            // Qhull requires at least 3 input points
    HW_MATH_ERR_QHULL_DIMS,            // Qhull requires at least 2 dimensions
    HW_MATH_ERR_QHULL_DIMS23,          // invalid matrix; must have either 2 or 3 columns
    HW_MATH_ERR_QHULL_NS_FACET,        // Qhull returned non-simplicial facet; try other options
    // unknown
    HW_MATH_ERR_UNKNOWN                // contact Altair customer support
};

//!
//! Warning messages - other messages are in Globals.cxx
//!
#define HW_MATH_MSG_SUCCESS         "Status: Success"

//!
//! Gets error message
//! \param code Error code
//!
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
