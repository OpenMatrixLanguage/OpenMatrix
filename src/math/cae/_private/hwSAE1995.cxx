/**
* @file hwSAE1995.cxx
* @date October 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
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
#include "hwSAE1995.h"

#include "hwMatrix.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwSAE1995::hwSAE1995()
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwSAE1995::~hwSAE1995()
{
}
//------------------------------------------------------------------------------
// Manages filter options, calls filter function and returns status
//------------------------------------------------------------------------------
hwMathStatus hwSAE1995::Filter(const hwMatrix& response, 
                               double          sampFreq,
                               double          cfc, 
                               int             stdpad,
                               int             direction, 
                               hwMatrix&       output)
{
    if (!response.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!response.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (sampFreq <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }
    if (cfc <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_FILTERCLASS, 3);
    }

    hwMathStatus status;
    hwMathStatus status2;
    if (cfc > 180.0)
    {
        // cfc = channel frequency class
        status2(HW_MATH_WARN_FILTERCFC, 3);
    }

    hwMatrix pad_resp;
    hwMatrix temp_resp;

    // pad input vectors in order to smooth out tails
    if (stdpad == 0)        // no padding
    {
        status = PadVec(response, 0, stdpad, pad_resp);
    }
    else if (stdpad < 0)    // pad with zeros
    {
        status = PadVec(response, 1, -stdpad, pad_resp);
    }
    else // (stdpad > 0)    // pad using mirror
    {
        status = PadVec(response, 3, stdpad, pad_resp);
    }
    if (!status.IsOk())
    {
        if (status.GetArg1() == 3)
        {
            status.SetArg1(4);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    // filter the padded data
    int size = 2 * abs(stdpad) + response.Size();
    hwMatrix filt_resp(size, hwMatrix::REAL);

    switch (direction) 
    {
        case 1:     // forward only
            Filter1D(pad_resp, sampFreq, cfc, 0, filt_resp);
            break;
        case 2:     // backward only
            Filter1D(pad_resp, sampFreq, cfc, 1, filt_resp);
            break;
        case 3:    // forward, then backward
            temp_resp.Dimension(size, hwMatrix::REAL);
            Filter1D(pad_resp, sampFreq, cfc, 0, temp_resp);
            Filter1D(temp_resp, sampFreq, cfc, 1, filt_resp);
            break;
        case 4:     // backward, then forward
            temp_resp.Dimension(size, hwMatrix::REAL);
            Filter1D(pad_resp, sampFreq, cfc, 1, temp_resp);
            Filter1D(temp_resp, sampFreq, cfc, 0, filt_resp);
            break;
        default:
            return status(HW_MATH_ERR_INVALIDINPUT, 5);
            break;
    }

    // crop the filtered data back down to the original length
    size   = response.Size();
    status = output.Dimension(size, hwMatrix::REAL);

    if (!status.IsOk())
    {
        return status;
    }

    for (int i = 0; i < size; ++i)
    {
        output(i) = filt_resp(i+abs(stdpad));
    }

    if (!status.IsOk())
    {
        return status;
    }
    return status2;
}
//------------------------------------------------------------------------------
// Pads the source vector
//------------------------------------------------------------------------------
hwMathStatus hwSAE1995::PadVec(const hwMatrix& source, 
                               int             padmode,
                               int             padcnt, 
                               hwMatrix&       source_pad)
{
    // padmode = 0 -> no padding, just reference the vector.
    // padmode = 1 -> pre/append "padcnt" zeros to vector.
    // padmode = 2 -> pre/append "padcnt" intervals to vector (time).
    // padmode = 3 -> pre/append "padcnt" mirrored values to vector.

    int          size   = source.Size();
    hwMathStatus status = source_pad.Dimension(2*padcnt + size, hwMatrix::REAL);
    if (!status.IsOk())
    {
        return status;
    }

    // create padded vector based on type and count specified
    int i = 0;
    const double* src_vec = source.GetRealData();
    double*       ret_vec = source_pad.GetRealData();
    double        delta   = (src_vec[size-1] - src_vec[0]) / (size-1);

    switch (padmode)
    {
        case 1: // pad with zeros, could also pad with constants
            while (i < padcnt)
            {
                // ret_vec[i] = src_vec[0];
                ret_vec[i] = 0.0;
                ++i;
            }

            while (i < padcnt + size)
            {
                ret_vec[i] = src_vec[i-padcnt];
                ++i;
            }

            while (i < padcnt + size + padcnt)
            {
                // ret_vec[i] = src_vec[size-1];
                ret_vec[i] = 0.0;
                ++i;
            }
            break;

        case 2: // pad with intervals (time)
            while (i < padcnt)
            {
                ret_vec[i] = -delta*(padcnt-i) + src_vec[0];
                ++i;
            }

            while (i < padcnt + size)
            {
                ret_vec[i] = src_vec[i-padcnt];
                ++i;
            }

            while (i < padcnt + size + padcnt)
            {
                ret_vec[i] = src_vec[size-1] + delta*(i-padcnt-size+1);
                ++i;
            }   
            break;

        case 3: // pad with mirror
            if (padcnt > size-1)
            {
                return status(HW_MATH_ERR_INVALIDINPUT, 3);
            }

            while (i < padcnt)
            {
                ret_vec[i] = 2.0*src_vec[0] - src_vec[padcnt-i];
                ++i;
            }

            while (i < padcnt+size)
            {
                ret_vec[i] = src_vec[i-padcnt];
                ++i;
            }

            while (i < padcnt+size+padcnt)
            {
                ret_vec[i] = 2.0*src_vec[size-1] - src_vec[2*size+padcnt-2-i];
                ++i;
            }
            break;

        default: // no padding
            source_pad = source;
            break;
    }

    return status;
}
//------------------------------------------------------------------------------
// Performs the actual filtering
//------------------------------------------------------------------------------
void hwSAE1995::Filter1D(const hwMatrix& response, 
                         double          sampFreq,
                         double          cfc, 
                         int             direction,
                         hwMatrix&       output)
{
    double wd = 2.0*PI*cfc*2.0775;  // 2.0775 = 5/3 * two-pass compensation for -6dB
    double wa = tan(wd/(2.0*sampFreq));
    double a0 = wa*wa/(1.0+sqrt(2.0)*wa+wa*wa);
    double a1 = 2.0*a0;
    double a2 = a0;
    double b1 = -2.0*(wa*wa-1.0)/(1.0+sqrt(2.0)*wa+wa*wa);
    double b2 = (-1.0+sqrt(2.0)*wa-wa*wa)/(1.0+sqrt(2.0)*wa+wa*wa);

    int size = response.Size();
    int i;
    const double* resp_vec = response.GetRealData();
    double*       ret_vec  = output.GetRealData();

    if (direction == 0)     // go forward
    {
        ret_vec[0] = resp_vec[0];
        ret_vec[1] = resp_vec[0];   // repeat previous point to avoid initial feedback spike
        i = 2;      // start at the beginning (almost) and go the end

        while (i < size)
        {
            ret_vec[i] = a0*resp_vec[i]
            + a1*resp_vec[i-1]
            + a2*resp_vec[i-2]
            + b1*ret_vec[i-1]
            + b2*ret_vec[i-2];
            ++i;
        }
    }
    else                    // go backward
    {
        ret_vec[size-1] = resp_vec[size-1];
        ret_vec[size-2] = resp_vec[size-1];   // repeat previous point to avoid initial feedback spike
        i = size-2; // start at the tail (almost) and go to the beginning

        while (i > 0)
        {
            --i;
            ret_vec[i] = a0*resp_vec[i]
            + a1*resp_vec[i+1]
            + a2*resp_vec[i+2]
            + b1*ret_vec[i+1]
            + b2*ret_vec[i+2];
        }
    }
}
