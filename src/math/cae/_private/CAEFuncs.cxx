/**
* @file CAEFuncs.cxx
* @date February 2017
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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

#include <vector>
#include <algorithm>
#include "CAEFuncs.h"

#include "hwMatrix.h"
#include "RainFlow.h"
#include "hwISO6487.h"
#include "hwSAE.h"
#include "hwSAE1995.h"

//------------------------------------------------------------------------------
// ISO 6487 filter
//------------------------------------------------------------------------------
hwMathStatus ISO6487(const hwMatrix& inSignal, 
                     double          sampFreq, 
                     double          cfc,
                     hwMatrix&       outSignal)
{
    hwISO6487    iso6487;
    hwMathStatus status = iso6487.Filter(inSignal, sampFreq, cfc, outSignal);

    return status;
}
//------------------------------------------------------------------------------
// SAE filter
//------------------------------------------------------------------------------
hwMathStatus SAEFilter(const hwMatrix& inSignal, 
                       double          sampFreq,
                       double          SAE_class, 
                       hwMatrix&       outSignal,
                       int             fftSize)
{
    if (!inSignal.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!inSignal.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (fftSize < 0)
    {
        hwMathStatus status(HW_MATH_ERR_NONNONNEGINT, 5);
    }
    
    int          numPnts = inSignal.Size();
    hwMathStatus status  = outSignal.Dimension(numPnts, hwMatrix::REAL);
    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(4);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    hwSAEFilter sae(SAE_class, sampFreq, fftSize);

    status = sae.Status();

    if (!status.IsOk())
    {
        switch(status.GetArg1())
        {
            case 1: status.SetArg1(3);   break;
            case 2: break;
            case 3: status.SetArg1(5);   break;
            default: status.ResetArgs(); break;
        }
        return status;
    }

    status = sae.Compute(inSignal, outSignal);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 2)
        {
            status.SetArg1(4);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    return status;
}
//------------------------------------------------------------------------------
// SAE filter 1995
//------------------------------------------------------------------------------
hwMathStatus SAEFilt95(const hwMatrix& inSignal, 
                       double          sampFreq, 
                       double          cfc,
                       int             stdpad, 
                       int             direction, 
                       hwMatrix&       outSignal)
{
    hwSAE1995    sae95;
    hwMathStatus status = sae95.Filter(inSignal, sampFreq, cfc, stdpad, 
                                       direction, outSignal);

    return status;
}
//------------------------------------------------------------------------------
// Rainflow commmon processing
//------------------------------------------------------------------------------
static
hwMathStatus RainFlowFunc(const hwMatrix& signal,
                          int             numBins,
                          int             hysteresis,
                          RainFlow&       rainflow,
                          hwMatrix&       index)
{
    // check inputs
    if (signal.IsEmpty())
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);

    if (!signal.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEXSUPPORT, 1);

    if (!signal.IsEmptyOrVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

    if (numBins < 1)
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);

    if (hysteresis != 0 && hysteresis != 1)
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 3);

    // pass data to rainflow
    int NoofPoints = signal.Size();
    SignalAttributes CurrentSignalData;

    CurrentSignalData.start = 0;
    CurrentSignalData.nChans = 1;
    CurrentSignalData.nVals = NoofPoints;

    rainflow.setsigstart(CurrentSignalData.start);
    rainflow.setsigNC(CurrentSignalData.nChans);
    rainflow.setsigNV(CurrentSignalData.nVals);

    hwMathStatus status = index.Dimension(NoofPoints, 1, hwMatrix::REAL);

    if (!status.IsOk())
        return status;

    double* _xchanneldata = index.GetRealData();
    const double* _ychanneldata = signal.GetRealData();

    for (int i = 0; i < CurrentSignalData.nVals; i++)
    {
        _xchanneldata[i] = static_cast<double> (i);
    }

    rainflow.setXchanneldata(_xchanneldata);
    rainflow.setYchanneldata(_ychanneldata);

    rainflow.setUser_NBins(numBins);
    rainflow.setHysteresis(hysteresis);

    return status;
}

//------------------------------------------------------------------------------
// Rainflow
//------------------------------------------------------------------------------
hwMathStatus RainFlowFunc(const hwMatrix& signal,
                          int             numBins,
                          int             hysteresis,
                          hwMatrix&       countMatrix,
                          hwMatrix&       damageMatrix,
                          hwMatrix&       rangeBinVec,
                          hwMatrix&       meanBinVec)
{
    RainFlow rainflow;
    hwMatrix index;

    // perform common processing
    hwMathStatus status = RainFlowFunc(signal, numBins, hysteresis,
                                       rainflow, index);

    if (!status.IsOk())
        return status;

    // run algorithm
    rainflow.Evaluate();

    // get outputs
    rangeBinVec  = *rainflow.getrangebins();
    meanBinVec   = *rainflow.getmeanbins();
    damageMatrix = *rainflow.getmeanrangematrix();
    countMatrix  = *rainflow.getcountmatrix();

    return status;
}
//------------------------------------------------------------------------------
// Rainflow
//------------------------------------------------------------------------------
hwMathStatus RainFlowFunc(const hwMatrix& signal,
                          int             numBins,
                          int             hysteresis,
                          const hwMatrix* countMatrix,
                          hwMatrix&       damageMatrix,
                          hwMatrix&       rangeBinVec,
                          hwMatrix&       meanBinVec)
{
    RainFlow rainflow;
    hwMatrix index;

    // perform common processing
    hwMathStatus status = RainFlowFunc(signal, numBins, hysteresis,
                                       rainflow, index);

    if (!status.IsOk())
        return status;

    // run algorithm
    int retv = rainflow.Evaluate(countMatrix);

    if (retv != 0)
        return status(HW_MATH_ERR_INTERNALERROR);

    // get outputs
    rangeBinVec = *rainflow.getrangebins();
    meanBinVec = *rainflow.getmeanbins();
    damageMatrix = *rainflow.getmeanrangematrix();

    return status;
}
