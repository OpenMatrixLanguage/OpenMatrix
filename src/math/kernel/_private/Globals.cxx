/**
* @file Globals.cxx
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

#include <Globals.h>

// Warnings
// General warnings
#define HW_MATH_MSG_NONNONNEGINT_W   "Warning: invalid value; must be nonnegative integer(s)"
// Matrix warnings
#define HW_MATH_MSG_SINGMATRIX       "Warning: singular matrix (to within machine precision)"
#define HW_MATH_MSG_SINGMATRIXDIV    "Warning: singular matrix divisor (to within machine precision)"
#define HW_MATH_MSG_MTXNOTFULLRANK_W "Warning: rank deficient matrix; matrix does not have full rank"
#define HW_MATH_MSG_MATRIXDEPENDCOL  "Warning: rank deficient matrix; columns are not linearly independent"
#define HW_MATH_MSG_NOSQRTM          "Warning: singular matrix; it may not have a square root"
#define HW_MATH_MSG_MAXSQRTM         "Warning: maximum number of allowed square roots has been exceeded"
#define HW_MATH_MSG_HANKEL           "Warning: arg1(end) ~= arg2(0); arg2(0) is ignored"
#define HW_MATH_MSG_TOEPLITZ         "Warning: arg1(0) ~= arg2(0); arg2(0) is ignored"
// Algorithm warnings
#define HW_MATH_MSG_MAXITERATE       "Warning: maximum iterations exceeded, returned last estimate"
#define HW_MATH_MSG_MAXFUNCEVAL      "Warning: maximum number of function evaluations exceeded, returned last estimate"
#define HW_MATH_MSG_NOTCONVERGE_W    "Warning: algorithm did not converge, returned best estimate"
#define HW_MATH_MSG_LOCALMIN         "Warning: algorithm converged at local minimum"
#define HW_MATH_MSG_NOSOLUTION       "Warning: no solution, returned best fit"
#define HW_MATH_MSG_NOUNIQUESOL      "Warning: no unique solution, returned best fit"
#define HW_MATH_MSG_TOOFEWPOINTS_W   "Warning: does not contain enough data points, no unique solution"
#define HW_MATH_MSG_NOPEAKS          "Warning: invalid input; has no peaks"
#define HW_MATH_MSG_GAMAXITERS_W     "Warning: maximum iterations will be capped at 1000"
#define HW_MATH_MSG_CONVIOLATEBEST   "Warning: the best reported design violates the constraint"
#define HW_MATH_MSG_POORPOLYFIT      "Warning: if poor fit quality, consult the help document"
// Signal processing warnings
#define HW_MATH_MSG_FILTERORDERIRR  "Warning: high filter order; numerical accuracy may be diminished"
#define HW_MATH_MSG_FILTERSIMPLOW   "Warning: syntax simplifies to low pass filter"
#define HW_MATH_MSG_FILTERSIMPHIGH  "Warning: syntax simplifies to high pass filter"
#define HW_MATH_MSG_FILTERSPEC_W    "Warning: extreme filter specification; numerical accuracy may be diminished"
#define HW_MATH_MSG_NYQUIST         "Warning: data exceeds Nyquist frequency"
#define HW_MATH_MSG_UNEQUALTIMESAMP "Warning: samples are not equally spaced"
#define HW_MATH_MSG_FILTERCFC       "Warning: invalid frequency; filter should not be used for CFC > 180"
// Statistics warnings
#define HW_MATH_MSG_STATTESTALPHA_W "Warning: invalid alpha; should be in the interval (0.0, 0.5)" 
// Qhull warnings
#define HW_MATH_MSG_QHULL_LEAK      "Warning: Qhull did not release all memory"

// Information messages
#define HW_MATH_MSG_TOLFCONV        "Information: function converged to within tolFun specification"
#define HW_MATH_MSG_TOLFCONV_R      "Information: function converged to within relative tolFun specification"
#define HW_MATH_MSG_TOLXCONV        "Information: step size converged to within tolX specification"
#define HW_MATH_MSG_TOLXCONV_R      "Information: step size converged to within relative tolX specification"
#define HW_MATH_MSG_SMALLTRUST      "Information: aborted due to trust region becoming too small"

// Errors
// General errors
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
#define HW_MATH_MSG_DIVIDEZERO      "Error: division by zero"
#define HW_MATH_MSG_NONFINITEDATA   "Error: invalid value; must be finite"
#define HW_MATH_MSG_NEGATIVE        "Error: invalid value; must be nonnegative"
#define HW_MATH_MSG_NONPOSITIVE     "Error: invalid value; must be positive"
#define HW_MATH_MSG_NONINTEGER      "Error: invalid value; must be integer"
#define HW_MATH_MSG_NONPOSINT       "Error: invalid value; must be positive integer"
#define HW_MATH_MSG_NONNONNEGINT    "Error: invalid value; must be nonnegative integer(s)"
#define HW_MATH_MSG_NONINCREASE     "Error: invalid data; must be strictly increasing"
#define HW_MATH_MSG_NONUNIQUE       "Error: invalid data; must contain all unique values"
#define HW_MATH_MSG_CONSECUTIVE3    "Error: invalid data; must contain at most 2 consecutive identical values"
#define HW_MATH_MSG_MINMAXVALUES    "Error: invalid values; maximum must be >= minimum"
#define HW_MATH_MSG_BADRANGE        "Error: invalid data; value is out of range"
#define HW_MATH_MSG_ZERORANGE       "Error: invalid data; has a range of zero"
#define HW_MATH_MSG_BADFILE         "Error: file missing or could not be read;"
// User function errors
#define HW_MATH_MSG_USERFUNCFAIL    "Error: user function failed; could not be evaluated"
#define HW_MATH_MSG_USERFUNCMATRIX  "Error: user function failed; must return a real matrix"
#define HW_MATH_MSG_USERFUNCREAL    "Error: user function failed; must return real data"
#define HW_MATH_MSG_USERFUNCREALNUM "Error: user function failed; must return a real number"
#define HW_MATH_MSG_USERFUNCREALMAT "Error: user function failed; must return a real matrix"
#define HW_MATH_MSG_USERFUNCSIZE    "Error: user function failed; has mismatched input/output dimensions"
#define HW_MATH_MSG_USERFUNCSIZE2   "Error: user function failed; has invalid output dimensions"
// Matrix errors
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
#define HW_MATH_MSG_ARRAYSIZE       "Error: incompatible matrices; dimensions must be consistent"
#define HW_MATH_MSG_ARRAYCOL3       "Error: invalid matrix; 3D data requires 3 columns"
#define HW_MATH_MSG_ARRAYTOOFEWROWS "Error: invalid matrix; cannot have fewer rows than columns"
#define HW_MATH_MSG_ARRAYTYPE       "Error: incompatible matrices; data types must be the same"
#define HW_MATH_MSG_COLUMNDIM       "Error: incompatible matrices; column dimensions must be equal"
#define HW_MATH_MSG_COLUMNVEC       "Error: invalid matrix; must be column vector"
#define HW_MATH_MSG_INVALIDINDEX    "Error: invalid index; must be within allowed bounds"
#define HW_MATH_MSG_DATANOTFOUND    "Error: invalid data; not found"
#define HW_MATH_MSG_SINGMATRIX_E    "Error: singular matrix (to within machine precision)"
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
// Slice errors
#define HW_MATH_MSG_SLICE_INDEX     "Error: invalid slice; dimension is inconsistent or out of bounds"
#define HW_MATH_MSG_SLICE_NUMCOLON  "Error: invalid syntax; can have only one non-colon"
// Algorithm errors
#define HW_MATH_MSG_NOTCONVERGE_E   "Error: algorithm did not converge"
#define HW_MATH_MSG_NOLOCALMIN      "Error: no local minimum near the initial condition"
#define HW_MATH_MSG_UNDERDETSYS_E   "Error: invalid system; must have at least as many equations as parameters"
#define HW_MATH_MSG_UNDERDETSYS_P   "Error: invalid system; must have at least as many points as parameters"
#define HW_MATH_MSG_TOOFEWPOINTS_E  "Error: invalid input; does not contain enough data points"
#define HW_MATH_MSG_DISTINCTVALS3   "Error: invalid input; must contain at least 3 distinct values"
#define HW_MATH_MSG_TIGHTTOL        "Error: invalid tolerance; may be too small"
#define HW_MATH_MSG_GAMAXITERS_E    "Error: requested maximum iterations has been exceeded"
#define HW_MATH_MSG_POPSIZE         "Error: invalid population size; must be 0 or greater than 2"
#define HW_MATH_MSG_GLOBALSEARCH    "Error: invalid global search; must be 1, 2 or 3"
#define HW_MATH_MSG_TOOFEWVEHICLES  "Error: invalid vehicles; not enough available"
#define HW_MATH_MSG_GAUSSLEGENDRE   "Error: invalid number of points; must be from 1 to 10"
#define HW_MATH_MSG_GAUSSLOBATTO    "Error: invalid number of points; must be from 2 to 10"
#define HW_MATH_MSG_RK45STEPSIZE    "Error: invalid step size; unnecessarily small"
#define HW_MATH_MSG_QUADSTEPSIZE    "Error: minimum allowed step size has been reached"
#define HW_MATH_MSG_FUNCTIONCOUNT   "Error: maximum allowed number of function calls has been reached"
#define HW_MATH_MSG_NOUNIQUEROOT    "Error: invalid polynomial; does not have a unique root"
#define HW_MATH_MSG_NESTSUPPORT     "Error: nested calls are not supported"
#define HW_MATH_MSG_INTERNALERROR   "Error: an internal error occurred"
// Signal processing errors
#define HW_MATH_MSG_FILTERTYPE      "Error: invalid filter; must be a supported type"
#define HW_MATH_MSG_FILTERORDER     "Error: invalid filter; must be positive"
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
// Statistics errors
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
// Qhull errors
#define HW_MATH_MSG_QHULL           "Error: Qhull function error"
#define HW_MATH_MSG_QHULL_PNTS      "Error: Qhull requires at least 3 input points"
#define HW_MATH_MSG_QHULL_DIMS      "Error: Qhull requires at least 2 dimensions"
#define HW_MATH_MSG_QHULL_DIMS23    "Error: invalid matrix; must have either 2 or 3 columns"
#define HW_MATH_MSG_QHULL_NS_FACET  "Error: Qhull returned non-simplicial facet; try other options"
// Other errors
#define HW_MATH_MSG_UNKNOWN         "Error: contact Altair customer support"


//! Get the string corresponding to an error code.
//! The string will be edited later to insert details.
std::string GetHMathErrMsg(hwMathMsgCode code)
{
    std::string retStr;

    if (code < HW_MATH_WARN_EOL)   // warning
    {
        switch (code)
        {
            case HW_MATH_ERR_NONE:              retStr = HW_MATH_MSG_SUCCESS;           break;
            // general warnings
            case HW_MATH_WARN_NONNONNEGINT:     retStr = HW_MATH_MSG_NONNONNEGINT_W;    break;
            // matrix warnings
            case HW_MATH_WARN_SINGMATRIX:       retStr = HW_MATH_MSG_SINGMATRIX;        break;
            case HW_MATH_WARN_SINGMATRIXDIV:    retStr = HW_MATH_MSG_SINGMATRIXDIV;     break;
            case HW_MATH_WARN_MTXNOTFULLRANK:   retStr = HW_MATH_MSG_MTXNOTFULLRANK_W;  break;
            case HW_MATH_WARN_MATRIXDEPENDCOL:  retStr = HW_MATH_MSG_MATRIXDEPENDCOL;   break;
            case HW_MATH_WARN_NOSQRTM:          retStr = HW_MATH_MSG_NOSQRTM;           break;
            case HW_MATH_WARN_MAXSQRTM:         retStr = HW_MATH_MSG_MAXSQRTM;          break;
            case HW_MATH_WARN_HANKEL:           retStr = HW_MATH_MSG_HANKEL;            break;
            case HW_MATH_WARN_TOEPLITZ:         retStr = HW_MATH_MSG_TOEPLITZ;          break;
            // algorithm warnings
            case HW_MATH_WARN_MAXITERATE:       retStr = HW_MATH_MSG_MAXITERATE;        break;
            case HW_MATH_WARN_MAXFUNCEVAL:      retStr = HW_MATH_MSG_MAXFUNCEVAL;       break;
            case HW_MATH_WARN_NOTCONVERGE:      retStr = HW_MATH_MSG_NOTCONVERGE_W;     break;
            case HW_MATH_WARN_LOCALMIN:         retStr = HW_MATH_MSG_LOCALMIN;          break;
            case HW_MATH_WARN_NOSOLUTION:       retStr = HW_MATH_MSG_NOSOLUTION;        break;
			case HW_MATH_WARN_NOUNIQUESOL:      retStr = HW_MATH_MSG_NOUNIQUESOL;       break;
			case HW_MATH_WARN_TOOFEWPOINTS:     retStr = HW_MATH_MSG_TOOFEWPOINTS_W;    break;
            case HW_MATH_WARN_NOPEAKS:          retStr = HW_MATH_MSG_NOPEAKS;           break;
            case HW_MATH_WARN_GAMAXITERS:       retStr = HW_MATH_MSG_GAMAXITERS_W;      break;
            case HW_MATH_WARN_CONVIOLATEBEST:   retStr = HW_MATH_MSG_CONVIOLATEBEST;    break;
            case HW_MATH_WARN_POORPOLYFIT:      retStr = HW_MATH_MSG_POORPOLYFIT;       break;
            // signal processing warnings
            case HW_MATH_WARN_FILTERORDERIRR:   retStr = HW_MATH_MSG_FILTERORDERIRR;    break;
            case HW_MATH_WARN_FILTERSIMPLOW:    retStr = HW_MATH_MSG_FILTERSIMPLOW;     break;
            case HW_MATH_WARN_FILTERSIMPHIGH:   retStr = HW_MATH_MSG_FILTERSIMPHIGH;    break;
            case HW_MATH_WARN_FILTERSPEC_W:     retStr = HW_MATH_MSG_FILTERSPEC_W;      break;
            case HW_MATH_WARN_NYQUIST:          retStr = HW_MATH_MSG_NYQUIST;           break;
            case HW_MATH_WARN_UNEQUALTIMESAMP:  retStr = HW_MATH_MSG_UNEQUALTIMESAMP;   break;
            case HW_MATH_WARN_FILTERCFC:        retStr = HW_MATH_MSG_FILTERCFC;         break;
            // statistics warnings
            case HW_MATH_WARN_STATTESTALPHA:    retStr = HW_MATH_MSG_STATTESTALPHA_W;   break;
            // Qhull warnings
            case HW_MATH_WARN_QHULL_LEAK:       retStr = HW_MATH_MSG_QHULL_LEAK;        break;
            // unknown error
			default:                            retStr = HW_MATH_MSG_UNKNOWN;           break;
        }
    }
    else if (code < HW_MATH_INFO_EOL)   // information
    {
        switch (code)
        {
            // general errors
            case HW_MATH_INFO_TOLFCONV:         retStr = HW_MATH_MSG_TOLFCONV;          break;
            case HW_MATH_INFO_TOLFCONV_R:       retStr = HW_MATH_MSG_TOLFCONV_R;        break;
            case HW_MATH_INFO_TOLXCONV:         retStr = HW_MATH_MSG_TOLXCONV;          break;
            case HW_MATH_INFO_TOLXCONV_R:       retStr = HW_MATH_MSG_TOLXCONV_R;        break;
            case HW_MATH_INFO_SMALLTRUST:       retStr = HW_MATH_MSG_SMALLTRUST;        break;
        }
    }
    else    // error
    {
        switch (code)
        {
            // general errors
            case HW_MATH_ERR_NULLPOINTER:           retStr = HW_MATH_MSG_NULLPOINTER;       break;
            case HW_MATH_ERR_INVALIDINPUT:          retStr = HW_MATH_MSG_INVALIDINPUT;      break;
            case HW_MATH_ERR_NOTSTRING:             retStr = HW_MATH_MSG_NOTSTRING;         break;
            case HW_MATH_ERR_INVALIDINTERVAL:       retStr = HW_MATH_MSG_INVALIDINTERVAL;   break;
            case HW_MATH_ERR_NOTALLOWED:            retStr = HW_MATH_MSG_NOTALLOWED;        break;
            case HW_MATH_ERR_NOTIMPLEMENT:          retStr = HW_MATH_MSG_NOTIMPLEMENT;      break;
            case HW_MATH_ERR_COMPLEX:               retStr = HW_MATH_MSG_COMPLEX;           break;
            case HW_MATH_ERR_NEEDCOMPLEX:           retStr = HW_MATH_MSG_NEEDCOMPLEX;       break;
            case HW_MATH_ERR_COMPLEXSUPPORT:        retStr = HW_MATH_MSG_COMPLEXSUPPORT;    break;
            case HW_MATH_ERR_ZERO:                  retStr = HW_MATH_MSG_ZERO;              break;
            case HW_MATH_ERR_DIVIDEZERO:            retStr = HW_MATH_MSG_DIVIDEZERO;        break;
            case HW_MATH_ERR_NONFINITEDATA:         retStr = HW_MATH_MSG_NONFINITEDATA;     break;
            case HW_MATH_ERR_NEGATIVE:              retStr = HW_MATH_MSG_NEGATIVE;          break;
            case HW_MATH_ERR_NONPOSITIVE:           retStr = HW_MATH_MSG_NONPOSITIVE;       break;
            case HW_MATH_ERR_NONINTEGER:            retStr = HW_MATH_MSG_NONINTEGER;        break;
            case HW_MATH_ERR_NONPOSINT:             retStr = HW_MATH_MSG_NONPOSINT;         break;
            case HW_MATH_ERR_NONNONNEGINT:          retStr = HW_MATH_MSG_NONNONNEGINT;      break;
            case HW_MATH_ERR_NONINCREASE:           retStr = HW_MATH_MSG_NONINCREASE;       break;
            case HW_MATH_ERR_NONUNIQUE:             retStr = HW_MATH_MSG_NONUNIQUE;         break;
            case HW_MATH_ERR_CONSECUTIVE3:          retStr = HW_MATH_MSG_CONSECUTIVE3;      break;
            case HW_MATH_ERR_MINMAXVALUES:          retStr = HW_MATH_MSG_MINMAXVALUES;      break;
            case HW_MATH_ERR_BADRANGE:              retStr = HW_MATH_MSG_BADRANGE;          break;
            case HW_MATH_ERR_ZERORANGE:             retStr = HW_MATH_MSG_ZERORANGE;         break;
            case HW_MATH_ERR_BADFILE:               retStr = HW_MATH_MSG_BADFILE;           break;
            // user function errors
            case HW_MATH_ERR_USERFUNCFAIL:          retStr = HW_MATH_MSG_USERFUNCFAIL;      break;
            case HW_MATH_ERR_USERFUNCMATRIX:        retStr = HW_MATH_MSG_USERFUNCMATRIX;    break;
            case HW_MATH_ERR_USERFUNCREAL:          retStr = HW_MATH_MSG_USERFUNCREAL;      break;
            case HW_MATH_ERR_USERFUNCREALNUM:       retStr = HW_MATH_MSG_USERFUNCREALNUM;   break;
            case HW_MATH_ERR_USERFUNCREALMAT:       retStr = HW_MATH_MSG_USERFUNCREALMAT;   break;
            case HW_MATH_ERR_USERFUNCSIZE:          retStr = HW_MATH_MSG_USERFUNCSIZE;      break;
            case HW_MATH_ERR_USERFUNCSIZE2:         retStr = HW_MATH_MSG_USERFUNCSIZE2;     break;
            // matrix errors
            case HW_MATH_ERR_OUTOFMEMORY:           retStr = HW_MATH_MSG_OUTOFMEMORY;       break;
            case HW_MATH_ERR_ALLOCFAILED:           retStr = HW_MATH_MSG_ALLOCFAILED;       break;
            case HW_MATH_ERR_EMPTYMATRIX:           retStr = HW_MATH_MSG_EMPTYMATRIX;       break;
            case HW_MATH_ERR_NONEMPTYMATRIX:        retStr = HW_MATH_MSG_NONEMPTYMATRIX;    break;
            case HW_MATH_ERR_MATRIXRESIZE:          retStr = HW_MATH_MSG_MATRIXRESIZE;      break;
            case HW_MATH_ERR_MATRIXRESHAPE1:        retStr = HW_MATH_MSG_MATRIXRESHAPE1;    break;
            case HW_MATH_ERR_MATRIXRESHAPE2:        retStr = HW_MATH_MSG_MATRIXRESHAPE2;    break;
            case HW_MATH_ERR_PERMVEC1:              retStr = HW_MATH_MSG_PERMVEC1;          break;
            case HW_MATH_ERR_PERMVEC2:              retStr = HW_MATH_MSG_PERMVEC2;          break;
            case HW_MATH_ERR_PERMVEC3:              retStr = HW_MATH_MSG_PERMVEC3;          break;
            case HW_MATH_ERR_PERMVEC4:              retStr = HW_MATH_MSG_PERMVEC4;          break;
            case HW_MATH_ERR_ARRAYDIM:              retStr = HW_MATH_MSG_ARRAYDIM;          break;
            case HW_MATH_ERR_ARRAYSIZE:             retStr = HW_MATH_MSG_ARRAYSIZE;         break;
            case HW_MATH_ERR_ARRAYCOL3:             retStr = HW_MATH_MSG_ARRAYCOL3;         break;
            case HW_MATH_ERR_ARRAYTOOFEWROWS:       retStr = HW_MATH_MSG_ARRAYTOOFEWROWS;   break;
            case HW_MATH_ERR_ARRAYTYPE:             retStr = HW_MATH_MSG_ARRAYTYPE;         break;
            case HW_MATH_ERR_COLUMNDIM:             retStr = HW_MATH_MSG_COLUMNDIM;         break;
            case HW_MATH_ERR_COLUMNVEC:             retStr = HW_MATH_MSG_COLUMNVEC;         break;
            case HW_MATH_ERR_INVALIDINDEX:          retStr = HW_MATH_MSG_INVALIDINDEX;      break;
            case HW_MATH_ERR_DATANOTFOUND:          retStr = HW_MATH_MSG_DATANOTFOUND;      break;
            case HW_MATH_ERR_SINGMATRIX:            retStr = HW_MATH_MSG_SINGMATRIX_E;      break;
            case HW_MATH_ERR_MTXNOTFULLRANK:        retStr = HW_MATH_MSG_MTXNOTFULLRANK_E;  break;
            case HW_MATH_ERR_MTXNOTSQUARE:          retStr = HW_MATH_MSG_MTXNOTSQUARE;      break;
            case HW_MATH_ERR_MTXNOTSYM:             retStr = HW_MATH_MSG_MTXNOTSYM;         break;
            case HW_MATH_ERR_MTXNOTSPD:             retStr = HW_MATH_MSG_MTXNOTSPD;         break;
            case HW_MATH_ERR_DECOMPFAIL:            retStr = HW_MATH_MSG_DECOMPFAIL;        break;
            case HW_MATH_ERR_DIAGFAIL:              retStr = HW_MATH_MSG_DIAGFAIL;          break;
            case HW_MATH_ERR_VECTOR12:              retStr = HW_MATH_MSG_VECTOR12;          break;
            case HW_MATH_ERR_VECTOR2:               retStr = HW_MATH_MSG_VECTOR2;           break;
            case HW_MATH_ERR_VECTOR:                retStr = HW_MATH_MSG_VECTOR;            break;
            case HW_MATH_ERR_NOVECTOR:              retStr = HW_MATH_MSG_NOVECTOR;          break;
            case HW_MATH_ERR_KLUD:                  retStr = HW_MATH_MSG_KLUD;              break;
            // slice errors
            case HW_MATH_ERR_SLICE_INDEX:           retStr = HW_MATH_MSG_SLICE_INDEX;       break;
            case HW_MATH_ERR_SLICE_NUMCOLON:        retStr = HW_MATH_MSG_SLICE_NUMCOLON;    break;
            // algorithm errors
            case HW_MATH_ERR_NOTCONVERGE:           retStr = HW_MATH_MSG_NOTCONVERGE_E;     break;
            case HW_MATH_ERR_NOLOCALMIN:            retStr = HW_MATH_MSG_NOLOCALMIN;        break;
            case HW_MATH_ERR_UNDERDETSYS_E:         retStr = HW_MATH_MSG_UNDERDETSYS_E;     break;
            case HW_MATH_ERR_UNDERDETSYS_P:         retStr = HW_MATH_MSG_UNDERDETSYS_P;     break;
            case HW_MATH_ERR_TOOFEWPOINTS:          retStr = HW_MATH_MSG_TOOFEWPOINTS_E;    break;
            case HW_MATH_ERR_DISTINCTVALS3:         retStr = HW_MATH_MSG_DISTINCTVALS3;     break;
            case HW_MATH_ERR_GAMAXITERS:            retStr = HW_MATH_MSG_GAMAXITERS_E;      break;
            case HW_MATH_ERR_TIGHTTOL:              retStr = HW_MATH_MSG_TIGHTTOL;          break;
            case HW_MATH_ERR_POPSIZE:               retStr = HW_MATH_MSG_POPSIZE;           break;
            case HW_MATH_ERR_GLOBALSEARCH:          retStr = HW_MATH_MSG_GLOBALSEARCH;      break;
            case HW_MATH_ERR_TOOFEWVEHICLES:        retStr = HW_MATH_MSG_TOOFEWVEHICLES;    break;
            case HW_MATH_ERR_GAUSSLEGENDRE:         retStr = HW_MATH_MSG_GAUSSLEGENDRE;     break;
            case HW_MATH_ERR_GAUSSLOBATTO:          retStr = HW_MATH_MSG_GAUSSLOBATTO;      break;
            case HW_MATH_ERR_RK45STEPSIZE:          retStr = HW_MATH_MSG_RK45STEPSIZE;      break;
            case HW_MATH_ERR_QUADSTEPSIZE:          retStr = HW_MATH_MSG_QUADSTEPSIZE;      break;
            case HW_MATH_ERR_FUNCTIONCOUNT:         retStr = HW_MATH_MSG_FUNCTIONCOUNT;     break;
            case HW_MATH_ERR_NOUNIQUEROOT:          retStr = HW_MATH_MSG_NOUNIQUEROOT;      break;
            case HW_MATH_ERR_NESTSUPPORT:           retStr = HW_MATH_MSG_NESTSUPPORT;       break;
            case HW_MATH_ERR_INTERNALERROR:         retStr = HW_MATH_MSG_INTERNALERROR;     break;
            // signal processing errors
            case HW_MATH_ERR_FILTERTYPE:            retStr = HW_MATH_MSG_FILTERTYPE;        break;
            case HW_MATH_ERR_FILTERORDER:           retStr = HW_MATH_MSG_FILTERORDER;       break;
            case HW_MATH_ERR_FILTERORDERODD:        retStr = HW_MATH_MSG_FILTERORDERODD;    break;
            case HW_MATH_ERR_FILTERSPEC_E:          retStr = HW_MATH_MSG_FILTERSPEC_E;      break;
            case HW_MATH_ERR_FILTERFREQ_A:          retStr = HW_MATH_MSG_FILTERFREQ_A;      break;
            case HW_MATH_ERR_FILTERFREQ_D:          retStr = HW_MATH_MSG_FILTERFREQ_D;      break;
            case HW_MATH_ERR_FILTERFREQS_EQ:        retStr = HW_MATH_MSG_FILTERFREQS_EQ;    break;
            case HW_MATH_ERR_FILTERBANDCONF:        retStr = HW_MATH_MSG_FILTERBANDCONF;    break;
            case HW_MATH_ERR_FILTERCLASS:           retStr = HW_MATH_MSG_FILTERCLASS;       break;
            case HW_MATH_ERR_FREQCLASS:             retStr = HW_MATH_MSG_FREQCLASS;         break;
            case HW_MATH_ERR_FILTERRIPPLE:          retStr = HW_MATH_MSG_FILTERRIPPLE;      break;
            case HW_MATH_ERR_FILTERDENZERO:         retStr = HW_MATH_MSG_FILTERDENZERO;     break;
            case HW_MATH_ERR_FILTFILTCOEFS:         retStr = HW_MATH_MSG_FILTFILTCOEFS;     break;
            case HW_MATH_ERR_FILTFILTDATA:          retStr = HW_MATH_MSG_FILTFILTDATA;      break;
            case HW_MATH_ERR_FIRWINDOW:             retStr = HW_MATH_MSG_FIRWINDOW;         break;
            case HW_MATH_ERR_DB_SIGN:               retStr = HW_MATH_MSG_DB_SIGN;           break;
            case HW_MATH_ERR_OVERLAPPOINTS:         retStr = HW_MATH_MSG_OVERLAPPOINTS;     break;
            case HW_MATH_ERR_FTBLOCKSIZE:           retStr = HW_MATH_MSG_FTBLOCKSIZE;       break;
            case HW_MATH_ERR_FFTSIZE:               retStr = HW_MATH_MSG_FFTSIZE;           break;
            case HW_MATH_ERR_PADARG:                retStr = HW_MATH_MSG_PADARG;            break;
            case HW_MATH_ERR_RESAMPOFFSET:          retStr = HW_MATH_MSG_RESAMPOFFSET;      break;
            // statistics errors
            case HW_MATH_ERR_INVALIDPROB:           retStr = HW_MATH_MSG_INVALIDPROB;       break;
            case HW_MATH_ERR_NORMTESTPNTS:          retStr = HW_MATH_MSG_NORMTESTPNTS;      break;
            case HW_MATH_ERR_STATTESTALPHA:         retStr = HW_MATH_MSG_STATTESTALPHA_E;   break;
            case HW_MATH_ERR_DOELEVEL:              retStr = HW_MATH_MSG_DOELEVEL;          break;
            // CVODE errors
            case HW_MATH_ERR_CVODE_WORKLOAD:        retStr = HW_MATH_MSG_CVODE_WORKLOAD;    break;
            case HW_MATH_ERR_CVODE_ACCURACY:        retStr = HW_MATH_MSG_CVODE_ACCURACY;    break;
            case HW_MATH_ERR_CVODE_ERR:             retStr = HW_MATH_MSG_CVODE_ERR;         break;
            case HW_MATH_ERR_CVODE_CONV:            retStr = HW_MATH_MSG_CVODE_CONV;        break;
            // IDA errors
            case HW_MATH_ERR_IDA_WORKLOAD:          retStr = HW_MATH_MSG_IDA_WORKLOAD;      break;
            case HW_MATH_ERR_IDA_ACCURACY:          retStr = HW_MATH_MSG_IDA_ACCURACY;      break;
            case HW_MATH_ERR_IDA_ERR:               retStr = HW_MATH_MSG_IDA_ERR;           break;
            case HW_MATH_ERR_IDA_CONV:              retStr = HW_MATH_MSG_IDA_CONV;          break;
            // Qhull errors
            case HW_MATH_ERR_QHULL:                 retStr = HW_MATH_MSG_QHULL;             break;
            case HW_MATH_ERR_QHULL_PNTS:            retStr = HW_MATH_MSG_QHULL_PNTS;        break;
            case HW_MATH_ERR_QHULL_DIMS:            retStr = HW_MATH_MSG_QHULL_DIMS;        break;
            case HW_MATH_ERR_QHULL_DIMS23:          retStr = HW_MATH_MSG_QHULL_DIMS23;      break;
            case HW_MATH_ERR_QHULL_NS_FACET:        retStr = HW_MATH_MSG_QHULL_NS_FACET;    break;
            // unknown error
            default:                                retStr = HW_MATH_MSG_UNKNOWN;
        }
    }

    return retStr;
}
