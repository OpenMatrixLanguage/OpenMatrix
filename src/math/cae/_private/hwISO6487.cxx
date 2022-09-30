/**
* @file  hwISO6487.cxx
* @date June 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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
#include "hwISO6487.h"

#include "hwMatrix.h"
#include "hwMathStatus.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwISO6487::hwISO6487()
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwISO6487::~hwISO6487()
{
}
//------------------------------------------------------------------------------
// Filters the signal and returns the status
//------------------------------------------------------------------------------
hwMathStatus hwISO6487::Filter(const hwMatrix& response, 
                               double          sampFreq, 
                               double          cfc,
                               hwMatrix&       output)
{
    hwMathStatus status;
    if (!response.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!response.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }
    if (sampFreq <= 0.0 || !IsFinite_T(sampFreq))
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (cfc <= 0.0)
    {
        return status(HW_MATH_ERR_FILTERCLASS, 3);
    }

    hwMathStatus status2;
    if (cfc > 180.0)
    {
        status2(HW_MATH_WARN_FILTERCFC, 3);
    }

    // allocate padded input
    double delta      = 1.0 / sampFreq;
    int    vectorsize = response.Size();

    if (vectorsize < 10)
    {
        return status(HW_MATH_ERR_TOOFEWPOINTS10, 1);
    }

    int numberOfAdditionalPoints = static_cast<int>(0.01 / delta);
    if (numberOfAdditionalPoints < 100)
    {
        numberOfAdditionalPoints = 100;
    }

    int padcnt = (numberOfAdditionalPoints > vectorsize-1) ?
                 vectorsize-1 : numberOfAdditionalPoints;

    hwMatrix pad_resp;

    status = FillIt(response, padcnt, pad_resp);
    if (!status.IsOk())
    {
        return status;
    }

    // filter forward, and then backward
    Filter(pad_resp, delta, cfc, 0);
    Filter(pad_resp, delta, cfc, 1);

    // crop the filtered data back down to the original length
    status = output.Dimension(vectorsize, hwMatrix::REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    double* filt_resp = pad_resp.GetRealData();
    double* ret_vec   = output.GetRealData();
    int i = 0;

    while (i < vectorsize)
    {
        ret_vec[i] = filt_resp[i+padcnt];
        i++;
    }

    if (!status.IsOk())
    {
        return status;
    }
    return status2;
}
//------------------------------------------------------------------------------
// Pads the signal and returns status
//------------------------------------------------------------------------------
hwMathStatus hwISO6487::FillIt(const hwMatrix& response, 
                               int             padcnt, 
                               hwMatrix&       pad_resp)
{
    // pad the provided vector with mirror strategy
    int size          = response.Size();
    int vectorsize    = response.Size();
    int newvectorsize = vectorsize + 2*padcnt;

    hwMathStatus status = pad_resp.Dimension(newvectorsize, hwMatrix::REAL);
    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    int i = padcnt;
    const double* src_vec = response.GetRealData();
    double*       ret_vec = pad_resp.GetRealData();

    while (i < vectorsize+padcnt)
    {
        ret_vec[i] = src_vec[i-padcnt];
        i++;
    }

    i = 0;

    while (i < padcnt)
    {
        ret_vec[padcnt - i - 1] = 2 * src_vec[0] - src_vec[i+1];
        ret_vec[vectorsize + padcnt + i] = 2 * src_vec[vectorsize - 1] - src_vec[vectorsize - i - 2];
        i++;
    }

    return status;
}
//------------------------------------------------------------------------------
// Filters the signal in place and returns status
//------------------------------------------------------------------------------
void hwISO6487::Filter(hwMatrix& response,
                       double    delta,
                       double    cfc,
                       int       dir)
{
    // define some constants needed for filtering
    double wd = 2.0*PI*cfc*25/12;
    double wa = sin(wd*delta/2.0)/cos(wd*delta/2.0);
    double b0 = wa*wa/(1.0+sqrt(2.0)*wa+wa*wa);
    double b1 = 2.0*b0;
    double b2 = b0;
    double a1 = -2.0*(wa*wa-1.0)/(1.0+sqrt(2.0)*wa+wa*wa);
    double a2 = (-1.0+sqrt(2.0)*wa-wa*wa)/(1.0+sqrt(2.0)*wa+wa*wa);
    double y1 = 0.0;
    double x0 = 0.0;
    double x1 = 0.0;
    double x2 = 0.0;
    double* resp_vec = response.GetRealData();
    int i;
    int vectorsize = response.Size();

    // perform the actual filtering
    if (dir == 0) // go forward
    {
        i = 0;

        while (i < 10)
        {
            y1 += resp_vec[i];
            i++;
        }

        y1 /= 10.0;
        x1 = resp_vec[0];
        x0 = resp_vec[1];
        resp_vec[0] = resp_vec[1] = y1;
        i = 2; // start at the beginning (almost) and go the end

        while (i < vectorsize)
        {
            x2 = x1;
            x1 = x0;
            x0 = resp_vec[i];
            resp_vec[i] = b0*x0 + a1*resp_vec[i-1] + a2*resp_vec[i-2] + b1*x1 + b2*x2;
            i++;
        }
    }
    else // go backward
    {
        i = vectorsize - 10;

        while (i < vectorsize)
        {
            y1 += resp_vec[i];
            i++;
        }

        y1 /= 10.0;
        x1 = resp_vec[vectorsize - 1];
        x0 = resp_vec[vectorsize - 2];
        resp_vec[vectorsize - 1] = resp_vec[vectorsize - 2] = y1;
        i = vectorsize - 3;

        while (i >= 0)
        {
            x2 = x1;
            x1 = x0;
            x0 = resp_vec[i];
            resp_vec[i] = b0*x0 + b1*x1 + b2*x2 + a1*resp_vec[i+1] + a2*resp_vec[i+2];
            i--;
        }
    }
}
