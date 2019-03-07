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
// Rainflow
//------------------------------------------------------------------------------
hwMathStatus RainFlowFunc(const hwMatrix& time, int numBins, double minRange, double maxRange,
                          int output, int hysteresis, hwMatrix& result)
{
    // check inputs
    if (time.IsEmpty())
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);

    if (!time.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEXSUPPORT, 1);

    if (!time.IsEmptyOrVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

    if (numBins < 1)
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);

    if (maxRange < minRange)
        return hwMathStatus(HW_MATH_ERR_MINMAXVALUES, 3, 4);

    if (output < 0 || output > 4)
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 5);

    if (hysteresis != 0 && hysteresis != 1)
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 6);

    // call algorithm
    RainFlow rainflow;

    int NoofPoints = time.Size();
    SignalAttributes CurrentSignalData;

    CurrentSignalData.start=0;
    CurrentSignalData.nChans=1;
    CurrentSignalData.nVals=NoofPoints;

    rainflow.setsigstart(CurrentSignalData.start);
    rainflow.setsigNC(CurrentSignalData.nChans);
    rainflow.setsigNV(CurrentSignalData.nVals);

    double* _xchanneldata  = new double[NoofPoints];
    double* _ychanneldata  = new double[NoofPoints];
    rainflow.setXchanneldata(_xchanneldata);
    rainflow.setYchanneldata(_ychanneldata);

    for (int i = 0;i<CurrentSignalData.nVals;i++)
    {
        _xchanneldata[i]=i;
        _ychanneldata[i]=time(i);
    }

    rainflow.setUser_NBins(numBins);
    rainflow.setminpoint(&minRange);
    rainflow.setmaxpoint(&maxRange);
	rainflow.setHysteresis(hysteresis);
    hwMathStatus status;

    if (output == 0)
    {
		rainflow.ComputeRangeBins();
        const STD::vector<double>& rangebins = rainflow.getrangebins();
        status = result.Dimension(1, (int) rangebins.size(), hwMatrix::REAL);

		for (int i = 0; i < (int) rangebins.size(); i++)
			result(i) = rangebins[i];
    }
    else if (output == 1)
    {
		rainflow.ComputeMeanBins();
        const STD::vector<double>& meanbins = rainflow.getmeanbins();
        status = result.Dimension(1, (int) meanbins.size(), hwMatrix::REAL);

		for (int i = 0; i < (int) meanbins.size(); i++)
			result(i) = meanbins[i];
    }
    else
    {
		rainflow.Evaluate();
        const STD::vector< STD::vector<double> >& meanrangematrix = rainflow.getmeanrangematrix();

        if (output == 2)
        {
            status = result.Dimension(numBins, numBins, hwMatrix::REAL);

            for (int ival=0; ival<numBins; ival++)
            {
                for (int jval=0; jval<numBins; jval++)
                {
                    result(ival, jval) = (meanrangematrix[ival])[jval]; 					
                }
            }
        }
        else if (output == 3)
        {
            status = result.Dimension(1, numBins, hwMatrix::REAL);

            result.SetElements(0.0);

            for (int ival=0; ival<numBins; ival++)
            {
                for (int jval=0; jval<numBins; jval++)
                {
                    result(ival) += (meanrangematrix[jval])[ival];
                }
            }
        }
        else // (output == 4)
        {
            status = result.Dimension(1, numBins, hwMatrix::REAL);

            result.SetElements(0.0);

            for (int ival=0; ival<numBins; ival++)
            {
                for (int jval=0; jval<numBins; jval++)
                {
                     result(ival) += (meanrangematrix[ival])[jval];
                }
            }
        }
    }

    return status;
}
