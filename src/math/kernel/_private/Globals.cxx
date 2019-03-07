/**
* @file Globals.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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

//:---------------------------------------------------------------------------
//:Description
//
//  Global functions
//
//:---------------------------------------------------------------------------

#include <Globals.h>

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
            // algorithm warnings
            case HW_MATH_WARN_MAXITERATE:       retStr = HW_MATH_MSG_MAXITERATE;        break;
            case HW_MATH_WARN_MAXFUNCEVAL:      retStr = HW_MATH_MSG_MAXFUNCEVAL;       break;
            case HW_MATH_WARN_NOTCONVERGE:      retStr = HW_MATH_MSG_NOTCONVERGE_W;     break;
            case HW_MATH_WARN_LOCALMIN:         retStr = HW_MATH_MSG_LOCALMIN;          break;
            case HW_MATH_WARN_NOSOLUTION:       retStr = HW_MATH_MSG_NOSOLUTION;        break;
            case HW_MATH_WARN_TOOFEWPOINTS:     retStr = HW_MATH_MSG_TOOFEWPOINTS_W;    break;
            case HW_MATH_WARN_NOPEAKS:          retStr = HW_MATH_MSG_NOPEAKS;           break;
            case HW_MATH_WARN_GAMAXITERS:       retStr = HW_MATH_MSG_GAMAXITERS_W;      break;
            case HW_MATH_WARN_CONVIOLATEBEST:   retStr = HW_MATH_MSG_CONVIOLATEBEST;    break;
            case HW_MATH_WARN_POORPOLYFIT:      retStr = HW_MATH_MSG_POORPOLYFIT;       break;
            // signal processing warnings
            case HW_MATH_WARN_FILTERSIMPLOW:    retStr = HW_MATH_MSG_FILTERSIMPLOW;     break;
            case HW_MATH_WARN_FILTERSIMPHIGH:   retStr = HW_MATH_MSG_FILTERSIMPHIGH;    break;
            case HW_MATH_WARN_FILTERSPEC_W:     retStr = HW_MATH_MSG_FILTERSPEC_W;      break;
            case HW_MATH_WARN_NYQUIST:          retStr = HW_MATH_MSG_NYQUIST;           break;
            case HW_MATH_WARN_UNEQUALTIMESAMP:  retStr = HW_MATH_MSG_UNEQUALTIMESAMP;   break;
            case HW_MATH_WARN_FILTERCFC:        retStr = HW_MATH_MSG_FILTERCFC;         break;
            // statistics warnings
            case HW_MATH_WARN_STATTESTALPHA:    retStr = HW_MATH_MSG_STATTESTALPHA_W;   break;
            case HW_MATH_WARN_SMALLVARIANCE:    retStr = HW_MATH_MSG_SMALLVARIANCE;     break;
            // Qhull warnings
            case HW_MATH_WARN_QHULL_LEAK:       retStr = HW_MATH_MSG_QHULL_LEAK;        break;
            // unknown error
            default:                            retStr = HW_MATH_MSG_UNKNOWN;
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
            // case HW_MATH_ERR_NONZERO:               retStr = HW_MATH_MSG_NONZERO;           break;
            case HW_MATH_ERR_DIVIDEZERO:            retStr = HW_MATH_MSG_DIVIDEZERO;        break;
            case HW_MATH_ERR_NONFINITEDATA:         retStr = HW_MATH_MSG_NONFINITEDATA;     break;
            case HW_MATH_ERR_NEGATIVE:              retStr = HW_MATH_MSG_NEGATIVE;          break;
            case HW_MATH_ERR_NONPOSITIVE:           retStr = HW_MATH_MSG_NONPOSITIVE;       break;
            case HW_MATH_ERR_NONINTEGER:            retStr = HW_MATH_MSG_NONINTEGER;        break;
            case HW_MATH_ERR_NONPOSINT:             retStr = HW_MATH_MSG_NONPOSINT;         break;
            case HW_MATH_ERR_NONNONNEGINT:          retStr = HW_MATH_MSG_NONNONNEGINT;      break;
            case HW_MATH_ERR_NONINCREASE:           retStr = HW_MATH_MSG_NONINCREASE;       break;
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
            case HW_MATH_ERR_ARRAYTOOFEWROWS:       retStr = HW_MATH_MSG_ARRAYTOOFEWROWS;   break;
            case HW_MATH_ERR_ARRAYTYPE:             retStr = HW_MATH_MSG_ARRAYTYPE;         break;
            case HW_MATH_ERR_COLUMNDIM:             retStr = HW_MATH_MSG_COLUMNDIM;         break;
            case HW_MATH_ERR_COLUMNVEC:             retStr = HW_MATH_MSG_COLUMNVEC;         break;
            case HW_MATH_ERR_INVALIDINDEX:          retStr = HW_MATH_MSG_INVALIDINDEX;      break;
            case HW_MATH_ERR_DATANOTFOUND:          retStr = HW_MATH_MSG_DATANOTFOUND;      break;
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
            // case HW_MATH_ERR_RK45WORKLOAD:          retStr = HW_MATH_MSG_RK45WORKLOAD;      break;
            case HW_MATH_ERR_QUADSTEPSIZE:          retStr = HW_MATH_MSG_QUADSTEPSIZE;      break;
            case HW_MATH_ERR_FUNCTIONCOUNT:         retStr = HW_MATH_MSG_FUNCTIONCOUNT;     break;
            case HW_MATH_ERR_NOUNIQUEROOT:          retStr = HW_MATH_MSG_NOUNIQUEROOT;      break;
            case HW_MATH_ERR_NESTSUPPORT:           retStr = HW_MATH_MSG_NESTSUPPORT;       break;
            case HW_MATH_ERR_INTERNALERROR:         retStr = HW_MATH_MSG_INTERNALERROR;     break;
            // signal processing errors
            case HW_MATH_ERR_FILTERTYPE:            retStr = HW_MATH_MSG_FILTERTYPE;        break;
            case HW_MATH_ERR_FILTERORDER:           retStr = HW_MATH_MSG_FILTERORDER;       break;
            case HW_MATH_ERR_FILTERORDERIRR:        retStr = HW_MATH_MSG_FILTERORDERIRR;    break;
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
