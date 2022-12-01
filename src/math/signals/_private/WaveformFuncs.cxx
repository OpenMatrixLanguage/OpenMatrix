/**
* @file * @file WaveformFuncs.cxx
* @date February 2022
* Copyright (C) 2007-2022 Altair Engineering, Inc.
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

#include "WaveformFuncs.h"
#include "hwMatrix.h"
#include "MKLutilities.h"
#include "PolynomFuncs.h"

#define MKLuD MKLutilitiesD

//------------------------------------------------------------------------------
// Computes a rectangular pulse and returns status
//------------------------------------------------------------------------------
hwMathStatus RectPulse(const hwMatrix& time,
                       double          width,
                       hwMatrix&       waveform)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!time.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (width <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    hwMathStatus status = waveform.Dimension(time.M(), time.N(), hwMatrix::REAL);
    int size = time.Size();
    double edge = width / 2.0;

    for (int i = 0; i < size; ++i)
    {
        if (time(i) >= -edge && time(i) < edge)
            waveform(i) = 1.0;
        else
            waveform(i) = 0.0;
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes a triangular pulse and returns status
//------------------------------------------------------------------------------
hwMathStatus TriPulse(const hwMatrix& time,
                      double          width,
                      double          skew,
                      hwMatrix&       waveform)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!time.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (width <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    if (skew < -1.0 || skew > 1.0)
    {
        return hwMathStatus(HW_MATH_ERR_BADRANGE, 3);   // improve this
    }

    hwMathStatus status = waveform.Dimension(time.M(), time.N(), hwMatrix::REAL);
    int size = time.Size();
    double edge = width / 2.0;
    double peak = skew * edge;

    waveform.SetElements(0.0);

    for (int i = 0; i < size; ++i)
    {
        if (time(i) >= -edge && time(i) <= peak)
            waveform(i) = (time(i) + edge) / (peak + edge);
        else if (time(i) > peak && time(i) < edge)
            waveform(i) = (time(i) - edge) / (peak - edge);
    }
 
    return status;
}
//------------------------------------------------------------------------------
// Computes a Gaussian pulse and returns status
//------------------------------------------------------------------------------
hwMathStatus GausPulse(const hwMatrix& time,
                       double          fc,
                       double          bw,
                       double          bwr,
                       hwMatrix&       waveform)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!time.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (fc < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 2);
    }

    if (bw <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 3);
    }

    if (bwr >= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_POSITIVE, 4);
    }

    // definition: bw = (fh - fl) / fc
    //             bw = ((fh - fc) + (fc - fl)) / fc
    //             bw = (delta_f + delta_f) / fc
    //             deltaF = bw * fc / 2
    // choose sigmaF such that:
    //             e^(-deltaF^2 / (2 * sigmaF^2)) = 10^(bwr / 20)
    //             -deltaF^2 / (2 * sigmaF^2) = (bwr / 20) * ln(10)
    //             sigmaF^2 = deltaF^2 / (-2 * (bwr / 20) * ln(10))
    //             sigmaF^2 = (bw * fc)^2 / (-8 * (bwr / 20) * ln(10))
    // convert to time domain:
    //             sigmaT^2 = 1 / (2 * pi * sigmaF)^2
    //             1 / sigmaT^2 = (2 * pi * sigmaF)^2

    double sigmaFsq = (bw * fc) * (bw * fc) / (-8.0 * (bwr / 20.0) * log(10.0));
    double sigmaTsqR = 4.0 * PI * PI * sigmaFsq;
    hwMatrix arg1 = MKLuD::PowerByElems(time, 2.0) * (-0.5 * sigmaTsqR);
    hwMatrix arg2 = (2.0 * PI * fc) * time;
    waveform = MKLuD::MultByElems(MKLuD::Exp(arg1), MKLuD::Cos(arg2));

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes a Gaussian pulse cutoff time and returns status
//------------------------------------------------------------------------------
hwMathStatus GausPulse(double  fc,
                       double  bw,
                       double  bwr,
                       double  tpr,
                       double& tc)
{
    if (fc < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 1);
    }

    if (bw <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    if (bwr >= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 3);
    }

    if (tpr >= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 4);
    }

    // see previous overload for derivation
    double sigmaFsq = (bw * fc) * (bw * fc) / (-8.0 * (bwr / 20.0) * log(10.0));
    double sigmaTsq = 1.0 / (4.0 * PI * PI * sigmaFsq);
    tc = sqrt(-2.0 * sigmaTsq * (tpr / 20.0) * log(10.0));

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes a Pulse train and returns status
//------------------------------------------------------------------------------
hwMathStatus PulsTran(const hwMatrix&    time,
                      const hwMatrix&    delay,
                      const std::string& func,
                      const hwMatrix&    args,
                      hwMatrix&          waveform)
{
    hwMathStatus status;
    hwMatrix pulswave;
    int argSize = args.Size();
    int numDelay;
    int indexN;

    if (!time.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!delay.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 2);
    }

    if (!args.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 4);
    }

    if (delay.IsVector())
    {
        numDelay = delay.Size();
        indexN = 0;
    }
    else
    {
        numDelay = delay.M();
        indexN = delay.N() - 1;

        if (delay.N() != 2)
        {
            return status(HW_MATH_ERR_INVALIDINPUT, 2);
        }
    }
   
    status = waveform.Dimension(time.M(), time.N(), hwMatrix::REAL);

    if (waveform.IsEmpty())
    {
        return status;
    }

    waveform.SetElements(0.0);

    hwMatrix tcopy(time);
    double* t = tcopy.GetRealData();
    double* w = waveform.GetRealData();
    int numTimes = time.Size();
    double ts = (time(numTimes - 1) - time(0)) / (numTimes - 1.0);
    double fs = 1.0 / ts;

    if (func == "rectpuls")
    {
        double width = 1.0;

        if (argSize == 1)
        {
            width = args(0);
        }
        else if (argSize > 1)
        {
            return status(HW_MATH_ERR_INVALIDINPUT, 4);
        }

        // pulse width processing
        int samplesPerPulse = static_cast<int> (ceil(width * fs - 1e-10));
        double w2 = width / 2.0;

        // process pulses
        for (int i = 0; i < numDelay; ++i)
        {
            int start_index = static_cast<int> (floor((delay(i) - w2) * fs + 1e-10));
            int stop_index = start_index + samplesPerPulse;
            start_index = _max(0, start_index);
            stop_index = _min(stop_index, numTimes - 1);

            if (start_index > stop_index)
                continue;

            double start_time = time(start_index);
            double stop_time = time(stop_index);

            // adjust for rounding error or possible non-constant sampling rate
            while (start_index > 0 && delay(i) - start_time < w2)
            {
                start_index--;
                start_time = time(start_index);
            }

            while (stop_index < numTimes - 1 && stop_time - delay(i) < w2)
            {
                stop_index++;
                stop_time = time(stop_index);
            }

            double* t_offset = t + start_index;
            double* w_offset = w + start_index;
            int samples = stop_index - start_index + 1;
            hwMatrix timeDelay(samples, 1, t_offset, hwMatrix::REAL);
            hwMatrix waveDelay(samples, 1, w_offset, hwMatrix::REAL);

            timeDelay -= delay(i, 0);

            status = RectPulse(timeDelay, width, pulswave);

            if (!status.IsOk())
            {
                int arg1 = status.GetArg1();

                if (arg1 > 1)
                {
                    status.SetArg1(arg1 + 2);
                }

                return status;
            }

            if (!indexN)
                waveDelay += pulswave;
            else
                waveDelay += pulswave * delay(i, 1);

            timeDelay += delay(i, 0);

            if (i == numDelay - 1)
                break;
        }
    }
    else if (func == "tripuls")
    {
        double width = 1.0;
        double skew = 0.0;

        if (argSize)
        {
            width = args(0);
        }

        if (argSize > 1)
        {
            skew = args(1);
        }

        if (argSize > 2)
        {
            return status(HW_MATH_ERR_INVALIDINPUT, 4);
        }

        // pulse width processing
        int samplesPerPulse = static_cast<int> (ceil(width * fs - 1e-10));
        double w2 = width / 2.0;

        // process pulses
        for (int i = 0; i < numDelay; ++i)
        {
            int start_index = static_cast<int> (floor((delay(i) - w2) * fs + 1e-10));
            int stop_index = start_index + samplesPerPulse;
            start_index = _max(0, start_index);
            stop_index = _min(stop_index, numTimes - 1);

            if (start_index > stop_index)
                continue;

            double start_time = time(start_index);
            double stop_time = time(stop_index);

            // adjust for rounding error or possible non-constant sampling rate
            while (start_index > 0 && delay(i) - start_time < w2)
            {
                start_index--;
                start_time = time(start_index);
            }

            while (stop_index < numTimes - 1 && stop_time - delay(i) < w2)
            {
                stop_index++;
                stop_time = time(stop_index);
            }

            double* t_offset = t + start_index;
            double* w_offset = w + start_index;
            int samples = stop_index - start_index + 1;
            hwMatrix timeDelay(samples, 1, t_offset, hwMatrix::REAL);
            hwMatrix waveDelay(samples, 1, w_offset, hwMatrix::REAL);

            timeDelay -= delay(i, 0);

            status = TriPulse(timeDelay, width, skew, pulswave);

            if (!status.IsOk())
            {
                int arg1 = status.GetArg1();

                if (arg1 > 1)
                {
                    status.SetArg1(arg1 + 2);
                }

                return status;
            }

            if (!indexN)
                waveDelay += pulswave;
            else
                waveDelay += pulswave * delay(i, 1);

            timeDelay += delay(i, 0);

            if (i == numDelay - 1)
                break;
        }
    }
    else if (func == "gauspuls")
    {
        double fc = 1000.0;
        double bw = 0.5;
        double bwr = -6.0;
        double tpr = -60.0;

        if (argSize)
        {
            fc = args(0);
        }

        if (argSize > 1)
        {
            bw = args(1);
        }

        if (argSize > 2)
        {
            bwr = args(2);
        }

        if (argSize > 3)
        {
            tpr = args(3);
        }

        if (argSize > 4)
        {
            return status(HW_MATH_ERR_INVALIDINPUT, 4);
        }

        // pulse width processing
        double w2;
        status = GausPulse(fc, bw, bwr, tpr, w2);
        
        if (!status.IsOk())
        {
            int arg1 = status.GetArg1();

            if (arg1 > 1)
            {
                status.SetArg1(arg1 + 3);
            }

            return status;
        }

        double width = 2.0 * w2;
        int samplesPerPulse = static_cast<int> (ceil(width * fs - 1e-10));

        // process pulses
        for (int i = 0; i < numDelay; ++i)
        {
            int start_index = static_cast<int> (floor((delay(i) - w2) * fs + 1e-10));
            int stop_index = start_index + samplesPerPulse;
            start_index = _max(0, start_index);
            stop_index = _min(stop_index, numTimes - 1);

            if (start_index > stop_index)
                continue;

            double start_time = time(start_index);
            double stop_time = time(stop_index);

            // adjust for rounding error or possible non-constant sampling rate
            while (start_index > 0 && delay(i) - start_time < w2)
            {
                start_index--;
                start_time = time(start_index);
            }

            while (stop_index < numTimes - 1 && stop_time - delay(i) < w2)
            {
                stop_index++;
                stop_time = time(stop_index);
            }

            double* t_offset = t + start_index;
            double* w_offset = w + start_index;
            int samples = stop_index - start_index + 1;
            hwMatrix timeDelay(samples, 1, t_offset, hwMatrix::REAL);
            hwMatrix waveDelay(samples, 1, w_offset, hwMatrix::REAL);

            timeDelay -= delay(i, 0);

            status = GausPulse(timeDelay, fc, bw, bwr, pulswave);

            if (!status.IsOk())
            {
                int arg1 = status.GetArg1();

                if (arg1 > 1)
                {
                    status.SetArg1(arg1 + 2);
                }

                return status;
            }

            if (!indexN)
                waveDelay += pulswave;
            else
                waveDelay += pulswave * delay(i, 1);

            timeDelay += delay(i, 0);

            if (i == numDelay - 1)
                break;
        }
    }
    else
    {
        return status(HW_MATH_ERR_INVALIDINPUT, 3);
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes a Pulse train and returns status
//------------------------------------------------------------------------------
hwMathStatus PulsTran(const hwMatrix&    time,
                      const hwMatrix&    delay,
                      const hwMatrix&    pulse,
                      double             fs_pulse,
                      const std::string& method,
                      hwMatrix&          waveform)
{
    hwMathStatus status;
    hwMatrix pulswave;
    int numDelay;
    int indexN;

    if (!time.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!delay.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 2);
    }

    if (delay.IsVector())
    {
        numDelay = delay.Size();
        indexN = 0;
    }
    else
    {
        numDelay = delay.M();
        indexN = delay.N() - 1;

        if (delay.N() != 2)
        {
            return status(HW_MATH_ERR_INVALIDINPUT, 2);
        }
    }

    if (!pulse.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 3);
    }

    if (fs_pulse <= 0.0)
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 4);
    }

    status = waveform.Dimension(time.M(), time.N(), hwMatrix::REAL);
    waveform.SetElements(0.0);

    hwMatrix tcopy(time);
    double* t = tcopy.GetRealData();
    double* w = waveform.GetRealData();
    int numTimes = time.Size();
    double ts = (time(numTimes - 1) - time(0)) / (numTimes - 1.0);
    double fs = 1.0 / ts;

    // process pulse variables
    double ts_pulse = 1.0 / fs_pulse;
    int pulseSize = pulse.Size();
    double width = (pulseSize - 1) * ts_pulse;
    double w2 = width / 2.0;
    int samplesPerPulse = static_cast<int> (ceil(width * fs - 1e-10));
    hwMatrix pulseTime(pulseSize, 1, hwMatrix::REAL);

    for (int i = 0; i < pulseSize; ++i)
    {
        pulseTime(i) = i * ts_pulse;
    }

    // process pulses
    for (int i = 0; i < numDelay; ++i)
    {
        int start_index = static_cast<int> (floor(delay(i) * fs + 1e-10));
        int stop_index = start_index + samplesPerPulse;
        start_index = _max(0, start_index);
        stop_index = _min(stop_index, numTimes - 1);

        if (start_index > stop_index)
            continue;

        double* t_offset = t + start_index;
        double* w_offset = w + start_index;
        int samples = stop_index - start_index + 1;
        hwMatrix timeDelay(samples, 1, t_offset, hwMatrix::REAL);
        hwMatrix waveDelay(samples, 1, w_offset, hwMatrix::REAL);

        timeDelay -= delay(i, 0);

        // compute pulswave
        if (method == "linear")
        {
            status = LinearInterp(pulseTime, pulse, timeDelay, pulswave, true, true);
        }
        else if (method == "pchip")
        {
            status = PchipInterp(pulseTime, pulse, timeDelay, pulswave, true);
        }
        else if (method == "spline")
        {
            status = Spline(pulseTime, pulse, timeDelay, pulswave, true);
        }
        else
        {
            status(HW_MATH_ERR_INVALIDINPUT, 5);
            break;
        }

        if (!indexN)
            waveDelay += pulswave;
        else
            waveDelay += pulswave * delay(i, 1);

        timeDelay += delay(i, 0);

        if (i == numDelay - 1)
            break;
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes a square pulse and returns status
//------------------------------------------------------------------------------
hwMathStatus SquarePulse(const hwMatrix& time,
                         double          duty,
                         hwMatrix&       waveform)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!time.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (duty <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    duty /= 100.0;

    hwMatrix tnorm = time / (2.0 * PI);
    waveform = MKLuD::Subtr(tnorm, MKLuD::Floor(tnorm));
    int size = waveform.Size();

    for (int i = 0; i < size; ++i)
    {
        if (waveform(i) < duty)
            waveform(i) = 1.0;
        else
            waveform(i) = -1.0;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes a sawtooth pulse and returns status
//------------------------------------------------------------------------------
hwMathStatus SawToothPulse(const hwMatrix& time,
                           double          width,
                           hwMatrix&       waveform)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!time.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (width < 0.0 || width > 1)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    hwMatrix tnorm = time / (2.0 * PI);
    hwMatrix tmod = MKLuD::Subtr(tnorm, MKLuD::Floor(tnorm));
    int size = tmod.Size();

    hwMathStatus status = waveform.Dimension(time.M(), time.N(), hwMatrix::REAL);
    waveform.SetElements(0.0);

    if (width != 0.0)
    {
        for (int i = 0; i < size; ++i)
        {
            if (tmod(i) < width)
            {
                waveform(i) = 2.0 * tmod(i) / width - 1.0;
            }
        }
    }

    if (width != 1.0)
    {
        for (int i = 0; i < size; ++i)
        {
            if (tmod(i) >= width)
            {
                waveform(i) = -2.0 * (tmod(i) - width) / (1.0 - width) + 1.0;
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes the Dirichlet function and returns status
//------------------------------------------------------------------------------
hwMathStatus Diric(const hwMatrix& time,
                   int             n,
                   hwMatrix&       waveform)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!time.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (n < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }

    double nn = static_cast<double> (n);

    waveform = MKLuD::DivideByElems(MKLuD::Sin(time * (nn / 2.0)),
                                    MKLuD::Sin(time / 2.0) * nn);
    double twoPI = 2.0 * PI;
    hwMatrix m = MKLuD::Mod(time, twoPI);
    int size = time.Size();

    for (int i = 0; i < size; ++i)
    {
        if (m(i) == 0.0)
        {
            waveform(i) = pow(-1.0, (nn - 1.0) * floor(time(i) / twoPI + 0.5));
        }
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes a chirp pulse and returns status
//------------------------------------------------------------------------------
hwMathStatus ChirpPulse(const hwMatrix& time,
                        double          f0,
                        double          t1,
                        double          f1,
                        const char*     shape,
                        double          phase,
                        hwMatrix&       waveform)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!time.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    phase *= (PI / 180.0);

    if (!strcmp(shape, "linear"))
    {
        double a = PI * (f1 - f0) / t1;
        double b = 2 * PI * f0;
        waveform = MKLuD::Cos(a * MKLuD::PowerByElems(time, 2.0) + b * time + phase);
    }
    else if (!strcmp(shape, "quadratic"))
    {
        double a = 2.0 / 3.0 * PI * (f1 - f0) / (t1 * t1);
        double b = 2 * PI * f0;
        waveform = MKLuD::Cos(a * MKLuD::PowerByElems(time, 3.0) + b * time + phase);
    }
    else if (!strcmp(shape, "logarithmic"))
    {
        double a = 2 * PI * f0 * t1 / log(f1 / f0);
        hwMatrix base(time.M(), time.N(), hwMatrix::REAL);
        base.SetElements(pow(f1 / f0, 1 / t1));
        waveform = MKLuD::Cos(a * MKLuD::PowerByElems(base, time) + phase);
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 5);
    }

    return hwMathStatus();
}
