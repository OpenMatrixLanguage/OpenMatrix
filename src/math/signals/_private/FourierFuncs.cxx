/**
* @file FourierFuncs.cxx
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

#include "FourierFuncs.h"

#include "hwMatrix.h"
#include "hwPSD.h"
#include "hwCSD.h"
#include "hwWindowFunc.h"

//------------------------------------------------------------------------------
// Get the mean sample rate from a time vector and compute the std deviation
//------------------------------------------------------------------------------
hwMathStatus SampleRate(const hwMatrix& time, 
                        double&         meanSampRate,
                        double&         stdDevRate, 
                        double          scale)
{
    if (!time.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!time.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int n = time.Size();
    if (n < 2)
    {
        return hwMathStatus(HW_MATH_ERR_TOOFEWPOINTS, 1);
    }

    if (scale <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 4);
    }

    double sum = 0.0;
    double sumSq = 0.0;
    double value;

    for (int i = 1; i < n; ++i)
    {
        value = (time(i) - time(i-1));

        if (value <= 0.0)
        {
            return hwMathStatus(HW_MATH_ERR_NONINCREASE, 1);
        }

        value = scale / value;    // use reciprocals as of 11.0 SA130

        sum += value;
        sumSq += value * value;
    }

    --n;
    meanSampRate = sum / n;
    value = (n * sumSq - sum * sum) / (double) (n * n);

    stdDevRate = (value > 0.0) ? sqrt(value) : 0.0;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Get the sampling frequency from a time vector and returns the status
//------------------------------------------------------------------------------
hwMathStatus SampleFreq(const hwMatrix& time,
                        double&         sampFreq)
{
    hwMathStatus status;

    if (!time.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!time.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }
    int numPnts = time.Size();

    if (numPnts < 2)
    {
        return status(HW_MATH_ERR_TOOFEWPOINTS, 1);
    }
    double delta_t0 = time(1) - time(0);

    if (delta_t0 <= 0.0)
    {
        return status(HW_MATH_ERR_NONINCREASE, 1);
    }
    if (numPnts == 2)
    {
        sampFreq = 1 / delta_t0;
        return status;
    }

    //// check for constant sampling rate
    // closely compare first and second halves
    int half = numPnts / 2;
    double delta_t1 = (time(half) - time(0)) / half;
    double delta_t2 = (time(numPnts-1) - time(half)) / (numPnts-1-half);

    if (delta_t1 <= 0.0)
    {
        return status(HW_MATH_ERR_NONINCREASE, 1);
    }
    if (delta_t2 <= 0.0)
    {
        return status(HW_MATH_ERR_NONINCREASE, 1);
    }
    if (fabs(delta_t1 / delta_t2 - 1.0) > 1.0e-5)
    {
        status(HW_MATH_WARN_UNEQUALTIMESAMP, 1);
    }

    // loosely compare first sample and full interval
    double delta_t = (time(numPnts-1) - time(0)) / (numPnts-1);

    if (fabs(delta_t0 / delta_t - 1.0) > 1.0e-2)
    {
        status(HW_MATH_WARN_UNEQUALTIMESAMP, 1);
    }
    // compute sampling rate
    sampFreq = 1 / delta_t;

    return status;
}
//------------------------------------------------------------------------------
// Create a frequency vector from a sampling frequency
//------------------------------------------------------------------------------
hwMathStatus Freq(int         numPnts,
                  double      sampFreq,
                  hwMatrix&   freq,
                  const char* option)
{
    if (numPnts < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);
    }

    if (sampFreq < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    hwMathStatus status = freq.Dimension(numPnts, hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(3);
        }
        return status;
    }

    if (!strcmp(option, "onesided"))
    {
        int n = 2 * (numPnts - 1) + numPnts % 2;
        double delta_f = sampFreq / n;
        freq(0) = 0.0;

        for (int i = 1; i < numPnts; ++i)
        {
            freq(i) = freq(i - 1) + delta_f;
        }
    }
    else if (!strcmp(option, "twosided"))
    {
        double delta_f = sampFreq / numPnts;
        freq(0) = 0.0;

        for (int i = 1; i < numPnts; ++i)
        {
            freq(i) = freq(i - 1) + delta_f;
        }
    }
    else if (!strcmp(option, "shift"))
    {
        double delta_f = sampFreq / numPnts;
        freq(0) = -(numPnts / 2) * delta_f;

        for (int i = 1; i < numPnts; ++i)
        {
            freq(i) = freq(i - 1) + delta_f;
        }
    }
    else
    {
        status(HW_MATH_ERR_INVALIDINPUT, 4);
    }

    return status;
}
//------------------------------------------------------------------------------
// Create a frequency vector from a time vector
//------------------------------------------------------------------------------
hwMathStatus Freq(const hwMatrix& time,
                  hwMatrix&       freq,
                  const char*     option)
{
    double sampFreq;
    hwMathStatus status_sf = SampleFreq(time, sampFreq);

    if (!status_sf.IsOk())
    {
        if (status_sf != HW_MATH_WARN_UNEQUALTIMESAMP)
        {
            return status_sf;
        }
    }

    hwMathStatus status = Freq(time.Size(), sampFreq, freq, option);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 3)
        {
            status.SetArg1(2);
        }
        else if (status.GetArg1() == 2)
        {
            status.SetArg1(1);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    return status_sf;
}
//------------------------------------------------------------------------------
// Create a 1-sided response vector by folding a 2-sided magnitude response vector
// A complex input is assumed to have conjugate symmetry
//------------------------------------------------------------------------------
hwMathStatus Fold(const hwMatrix& mag_twosided,
                  hwMatrix&       mag_onesided)
{
    if (!mag_twosided.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int numPnts2 = mag_twosided.Size();
    int numPnts1 = numPnts2/2 + 1;
    hwMathStatus status;

    if (mag_twosided.M() == 1)
    {
        status = mag_onesided.Dimension(1, numPnts1, mag_twosided.Type());
    }
    else if (mag_twosided.N() == 1)
    {
        status = mag_onesided.Dimension(numPnts1, 1, mag_twosided.Type());
    }
    else
    {
        status = mag_onesided.Dimension(0, 0, hwMatrix::REAL);
    }

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(2);
        }
        return status;
    }

    if (numPnts1 < 1)
    {
        return status;
    }

	if (mag_onesided.IsReal())
	{
		mag_onesided(0) = mag_twosided(0);

		for (int i = 1; i < numPnts1-1; ++i)
		{
			mag_onesided(i) = 2.0 * mag_twosided(i);
		}

		if (numPnts2%2 == 0)
		{
			mag_onesided(numPnts1-1) = mag_twosided(numPnts1-1);
		}
		else
		{
			mag_onesided(numPnts1-1) = 2.0 * mag_twosided(numPnts1-1);
		}
	}
	else
	{
		mag_onesided.z(0) = mag_twosided.z(0);

		for (int i = 1; i < numPnts1-1; ++i)
		{
			mag_onesided.z(i) = 2.0 * mag_twosided.z(i);
		}

		if (numPnts2%2 == 0)
		{
			mag_onesided.z(numPnts1-1) = mag_twosided.z(numPnts1-1);
		}
		else
		{
			mag_onesided.z(numPnts1-1) = 2.0 * mag_twosided.z(numPnts1-1);
		}
	}

    return status;
}
//------------------------------------------------------------------------------
// 2D FFT of a real or complex signal
//------------------------------------------------------------------------------
hwMathStatus Fft2(const hwMatrix& signal,
                  hwMatrix&       freqRes)
{
    int m = signal.M();
    int n = signal.N();

    hwMathStatus status = freqRes.Dimension(m, n, hwMatrix::COMPLEX);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(2);
        }
        else
        {
            status.ResetArgs();
        }

        return status;
    }

    if (freqRes.IsEmpty())
    {
        return status;
    }

    if (signal.IsReal())
    {
        int mm = m/2 + 1;
        int nn = n/2 + 1;
        hwMatrix* temp = (hwMatrix*) &signal;
        double* sr = temp->GetRealData();
        fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n * mm);     

        if (!out)
        {
            return status(HW_MATH_ERR_OUTOFMEMORY);
        }

        unsigned flags = FFTW_ESTIMATE | FFTW_PRESERVE_INPUT;
        fftw_plan p    = fftw_plan_dft_r2c_2d(n, m, sr, out, flags);

        if (!p)
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }

        fftw_execute(p);

        int i;
        int j;
        int ii;
        int jj;     // jj will replace index
        int count = 0;
        double real;
        double imag;

        if (m%2 == 0)
        {
            --mm;
        }

        if (n%2 == 0)
        {
            --nn;
        }

        // process DC column
        if (true)
        {
            // process DC row
            freqRes.z(0, 0).Set(out[count][0], out[count][1]);
            ++count;

            // process positive and negative frequency rows up to last odd jj
            ii = m - 1;

            for (i = 1; i < mm; ++i)
            {
                real = out[count][0];
                imag = out[count][1];
                freqRes.z(i, 0).Set(real, imag);
                freqRes.z(ii, 0).Set(real, -imag);    // phase of FT of real signal is an odd function
                ++count;
                --ii;
            }

            // process max row frequency if an even number
            if (m%2 == 0)
            {
                freqRes.z(i, 0).Set(out[count][0], out[count][1]);
                ++count;
            }
        }

        // process positive frequency columns up to max
        for (j = 1; j < nn; ++j)
        {
            // process DC row
            freqRes.z(0, j).Set(out[count][0], out[count][1]);
            ++count;

            // process positive and negative frequency rows up to last odd jj
            ii = m - 1;
            jj = n - j;

            for (i = 1; i < mm; ++i)
            {
                real = out[count][0];
                imag = out[count][1];
                freqRes.z(i, j).Set(real, imag);
                freqRes.z(ii, jj).Set(real, -imag);   // phase of FT of real signal is an odd function
                ++count;
                --ii;
            }

            // process max row frequency if an even number
            if (m%2 == 0)
            {
                freqRes.z(i, j).Set(out[count][0], out[count][1]);
                ++count;
            }
        }

        // process max column frequency if an even number
        if (n%2 == 0)
        {
            // process DC row
            freqRes.z(0, j).Set(out[count][0], out[count][1]);
            ++count;

            // process positive and negative frequency rows up to last odd jj
            ii = m - 1;

            for (i = 1; i < mm; ++i)
            {
                real = out[count][0];
                imag = out[count][1];
                freqRes.z(i, j).Set(real, imag);
                freqRes.z(ii, j).Set(real, -imag);   // phase of FT of real signal is an odd function
                ++count;
                --ii;
            }

            // process max row frequency if an even number
            if (m%2 == 0)
            {
                freqRes.z(i, j).Set(out[count][0], out[count][1]);
                ++count;
            }

            ++j;
        }

        // process negative frequency columns down to DC
        for ( ; j < n; ++j)
        {
            // process DC row
            freqRes.z(0, j).Set(out[count][0], out[count][1]);
            ++count;

            // process positive and negative frequency rows up to last odd jj
            ii = m - 1;
            jj = n - j;

            for (i = 1; i < mm; ++i)
            {
                real = out[count][0];
                imag = out[count][1];
                freqRes.z(i, j).Set(real, imag);
                freqRes.z(ii, jj).Set(real, -imag);   // phase of FT of real signal is an odd function
                ++count;
                --ii;
            }

            // process max row frequency if an even number
            if (m%2 == 0)
            {
                freqRes.z(i, j).Set(out[count][0], out[count][1]);
                ++count;
            }
        }

        fftw_destroy_plan(p); 
        fftw_free(out);
    }
    else
    {
        hwMatrix* temp = (hwMatrix*) &signal;
        fftw_complex* sc  = (fftw_complex*) temp->GetComplexData();
        fftw_complex* out = (fftw_complex*) freqRes.GetComplexData();     

        if (!out)
        {
            return status(HW_MATH_ERR_OUTOFMEMORY);
        }

        unsigned flags = FFTW_ESTIMATE | FFTW_PRESERVE_INPUT;

        fftw_plan p = fftw_plan_dft_2d(n, m, sc, out, FFTW_FORWARD, flags);
        if (!p)
        {
            return status(HW_MATH_ERR_ALLOCFAILED);
        }
        fftw_execute(p);
        fftw_destroy_plan(p); 
    }

    return status;
}
//------------------------------------------------------------------------------
// 2D FFT of a real or complex signal
//------------------------------------------------------------------------------
hwMathStatus Fft2(const hwMatrix& signal,
                  int             m,
                  int             n,
                  hwMatrix&       freqRes)
{
    hwMatrix copy(signal);

    hwMathStatus status = copy.Resize(m, n, true);
    if (!status.IsOk())
    {
        if (status.GetArg1() == 1)
        {
            status.SetArg1(2);
        }
        else if (status.GetArg1() == 2)
        {
            status.SetArg1(3);
        }
        else
        {
            status.ResetArgs();
        }
    }

    status = Fft2(copy, freqRes);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 2)
        {
            status.SetArg1(3);
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// FFT of a real or complex signal
//------------------------------------------------------------------------------
hwMathStatus Fft(const hwMatrix& signal,
                 hwMatrix&       freqRes,
                 int             fftSize)
{
    hwFFT_f fft(fftSize);

    hwMathStatus status = fft.Status();

    if (!status.IsOk())
    {
        status.SetArg1(3);
        return status;
    }

    status = fft.Compute(signal, freqRes);

    return status;
}
//------------------------------------------------------------------------------
// FFT of a real or complex signal
//------------------------------------------------------------------------------
hwMathStatus Fft(const hwMatrix& signal,
                 hwMatrix&       freqRes,
                 int             dim,
                 int             fftSize)   // ignored for now
{
    hwFFT_f fft;

    hwMathStatus status = fft.Status();

    if (!status.IsOk())
    {
        status.SetArg1(3);
        return status;
    }

    status = fft.Compute(signal, dim, freqRes);

    return status;
}
//------------------------------------------------------------------------------
// FFT of a real or complex signal
//------------------------------------------------------------------------------
hwMathStatus Fft(const hwMatrixN& signal,
                 hwMatrixN&       freqRes,
                 int              dim,
                 int              fftSize)   // ignored for now
{
    hwFFT_f fft;

    hwMathStatus status = fft.Status();

    if (!status.IsOk())
    {
        status.SetArg1(3);
        return status;
    }

    status = fft.Compute(signal, dim, freqRes);

    return status;
}
//------------------------------------------------------------------------------
// Inverse FFT of a real or complex spectrum
//------------------------------------------------------------------------------
hwMathStatus Ifft(const hwMatrix& freqRes,
                  hwMatrix&       signal,
                  int             fftSize)
{
    hwFFT_r ifft(fftSize, false);

    hwMathStatus status = ifft.Status();
    if (!status.IsOk())
    {
        status.SetArg1(3);
        return status;
    }

    status = ifft.Compute(freqRes, signal);

    return status;
}
//------------------------------------------------------------------------------
// Inverse FFT of a real or complex spectrum
//------------------------------------------------------------------------------
hwMathStatus Ifft(const hwMatrix& freqRes,
                  hwMatrix&       signal,
                  int             dim,
                  int             fftSize)   // ignored for now
{
    hwFFT_r ifft(0, false);

    hwMathStatus status = ifft.Status();
    if (!status.IsOk())
    {
        status.SetArg1(3);
        return status;
    }

    status = ifft.Compute(freqRes, dim, signal);

    return status;
}
//------------------------------------------------------------------------------
// Inverse FFT of a real or complex spectrum
//------------------------------------------------------------------------------
hwMathStatus Ifft(const hwMatrixN& freqRes,
                  hwMatrixN&       signal,
                  int              dim,
                  int              fftSize)   // ignored for now
{
    hwFFT_r ifft(0, false);

    hwMathStatus status = ifft.Status();
    if (!status.IsOk())
    {
        status.SetArg1(3);
        return status;
    }

    status = ifft.Compute(freqRes, dim, signal);

    return status;
}
//------------------------------------------------------------------------------
// Circular convolution of a real signal pair for periodic signals 
//------------------------------------------------------------------------------
hwMathStatus ConvCirc(const hwMatrix& signal1, 
                      const hwMatrix& signal2,
                      hwMatrix&       conv, 
                      int             fftSize)
{
    if (!signal1.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!signal1.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (!signal2.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }

    if (!signal2.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }

    if (signal2.Size() != signal1.Size())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    // compute fft of each signal
    hwFFT_f fft(fftSize);   

    hwMathStatus status = fft.Status();
    if (!status.IsOk())
    {
        status.SetArg1(4);
        return status;
    }

    hwMatrix freqRes1;
    status = fft.Compute(signal1, freqRes1);
    if (!status.IsOk())
    {
        if (status.GetArg1() != 1)
        {
            status.ResetArgs();
        }
        return status;
    }

    hwMatrix freqRes2;
    status = fft.Compute(signal2, freqRes2);
    if (!status.IsOk())
    {
        if (status.GetArg1() == 1)
        {
            status.SetArg1(2);
        }
        else
        {
            status.ResetArgs();
        }

        return status;
    }

    // compute convolution
    hwFFT_r ifft(fftSize);
    hwMatrix product;

    if (freqRes2.M() != freqRes1.M())
    {
        status = freqRes2.Transpose();
    }

    status = product.MultByElems(freqRes1, freqRes2);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = ifft.Compute(product, conv);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 2)
        {
            status.SetArg1(3);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    double factor = 1.0 / signal1.Size();
    conv *= factor;

    return status;
}
//------------------------------------------------------------------------------
// Circular correlation of a real signal pair for periodic signals 
//------------------------------------------------------------------------------
hwMathStatus CorrCirc(const hwMatrix& signal1, 
                      const hwMatrix& signal2,
                      hwMatrix&       corr, 
                      int             fftSize)
{
    if (!signal1.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!signal1.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (!signal2.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }

    if (!signal2.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }

    if (signal2.Size() != signal1.Size())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    // compute fft of each signal
    hwFFT_f fft(fftSize);
    hwMathStatus status = fft.Status();
    if (!status.IsOk())
    {
        status.SetArg1(4);
        return status;
    }

    hwMatrix freqRes1;
    status = fft.Compute(signal1, freqRes1);
    if (!status.IsOk())
    {
        return status;
    }

    hwMatrix freqRes2;
    status = fft.Compute(signal2, freqRes2);
    if (!status.IsOk())
    {
        return status;
    }

    // compute correlation
    hwFFT_r ifft(fftSize);
    hwMatrix freqResSq1;
    hwMatrix freqResSq2;
    hwMatrix auto1;
    hwMatrix auto2;
    hwMatrix cross;

    status = freqResSq1.AbsSq(freqRes1);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    freqResSq1.MakeComplex();

    status = ifft.Compute(freqResSq1, auto1);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = freqResSq2.AbsSq(freqRes2);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    freqResSq2.MakeComplex();    // should not have to do this

    status = ifft.Compute(freqResSq2, auto2);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    freqRes2.Conjugate();

    if (freqRes2.M() != freqRes1.M())
    {
        status = freqRes2.Transpose();
    }

    status = cross.MultByElems(freqRes1, freqRes2);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    status = ifft.Compute(cross, corr);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    double factor = 1.0 / sqrt(auto1(0) * auto2(0));
    corr *= factor;

    return status;
}
//------------------------------------------------------------------------------
// Coherence of a real signal pair with a window function for periodic signals 
//------------------------------------------------------------------------------
hwMathStatus Coherence(const hwMatrix& sysInput, 
                       const hwMatrix& sysOutput,
                       const hwMatrix& window, 
                       int             num_overlap_points,
                       hwMatrix&       cohere, 
                       int             fftSize)
{
    if (!sysInput.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!sysInput.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (!sysOutput.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!sysOutput.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    if (sysOutput.Size() != sysInput.Size())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    if (!window.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);
    }
    if (!window.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 3);
    }

    int numPnts   = sysInput.Size();
    int blockSize = window.Size();
    hwWindowFunc windowFunc(window);

    if (blockSize > numPnts)
    {
        blockSize = numPnts;
    }

    if (num_overlap_points >= blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_OVERLAPPOINTS, 4);
    }

    if (num_overlap_points < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 4);
    }

    if (fftSize == 0)
    {
        fftSize = blockSize;
    }
    else if (fftSize < blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_FFTSIZE, 3, 6);
    }

    int num_ffts = (int)floor((double)(numPnts-blockSize)/
                   (double)(blockSize-num_overlap_points)) + 1;

    if (num_ffts < 1)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 3);
    }

    hwMathStatus status = cohere.Dimension(fftSize, hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(5);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    hwFFT_f fft(fftSize);
    hwMatrix temp_real_in(fftSize, hwMatrix::REAL);
    hwMatrix temp_real_out(fftSize, hwMatrix::REAL);
    hwMatrix input_spec(fftSize, hwMatrix::COMPLEX);
    hwMatrix output_spec(fftSize, hwMatrix::COMPLEX);
    hwMatrix cross_spec(fftSize, hwMatrix::COMPLEX);
    hwMatrix input_mag(fftSize, hwMatrix::REAL);
    hwMatrix output_mag(fftSize, hwMatrix::REAL);

    cross_spec.SetElements(hwComplex());
    input_mag.SetElements(0.0);
    output_mag.SetElements(0.0);

    int index;

    for (int k = 0; k < num_ffts; ++k)
    {
        temp_real_in.SetElements(0.0);
        temp_real_out.SetElements(0.0);
        index = k * (blockSize-num_overlap_points);

        windowFunc.ApplyWindow(sysInput, index, temp_real_in, true);
        windowFunc.ApplyWindow(sysOutput, index, temp_real_out, true);

        fft.Compute(temp_real_in, input_spec);
        fft.Compute(temp_real_out, output_spec);

        input_spec.Conjugate();

        for (int j = 0; j < fftSize; ++j)
        {
            cross_spec.z(j) += input_spec.z(j) * output_spec.z(j);
            input_mag(j) += input_spec.z(j).MagSq();
            output_mag(j) += output_spec.z(j).MagSq();
        }
    }

    double numer;
    double denom;

    for (int k = 0; k < fftSize; ++k)
    {
        numer = cross_spec.z(k).MagSq();
        denom = input_mag(k) * output_mag(k);

        cohere(k) = (denom == 0.0) ? 1.0 : numer / denom;
    }

    return status;
}
//------------------------------------------------------------------------------
// Coherence of a real signal pair for periodic signals 
//------------------------------------------------------------------------------
hwMathStatus Coherence(const hwMatrix& sysInput, 
                       const hwMatrix& sysOutput,
                       int             blockSize,
                       int             num_overlap_points,
                       hwMatrix&       cohere, 
                       int             fftSize)
{
    if (!sysInput.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!sysInput.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (!sysOutput.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!sysOutput.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    if (sysOutput.Size() != sysInput.Size())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    int numPnts = sysInput.Size();
    if (blockSize > numPnts)
    {
        blockSize = numPnts;
    }

    if (blockSize < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 3);
    }

    if (num_overlap_points >= blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_OVERLAPPOINTS, 4);
    }

    if (num_overlap_points < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 4);
    }

    if (fftSize == 0)
    {
        fftSize = blockSize;
    }
    else if (fftSize < blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_FFTSIZE, 3, 6);
    }

    int num_ffts = (int)floor((double)(numPnts-blockSize)/
                   (double)(blockSize-num_overlap_points)) + 1;

    if (num_ffts < 1)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 3);
    }
    hwMathStatus status = cohere.Dimension(fftSize, hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(5);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    hwFFT_f fft(fftSize);

    hwMatrix temp_real_in(blockSize, hwMatrix::REAL);
    hwMatrix temp_real_out(blockSize, hwMatrix::REAL);
    hwMatrix input_spec(fftSize, hwMatrix::COMPLEX);
    hwMatrix output_spec(fftSize, hwMatrix::COMPLEX);
    hwMatrix cross_spec(fftSize, hwMatrix::COMPLEX);
    hwMatrix input_mag(fftSize, hwMatrix::REAL);
    hwMatrix output_mag(fftSize, hwMatrix::REAL);

    cross_spec.SetElements(hwComplex());
    input_mag.SetElements(0.0);
    output_mag.SetElements(0.0);

    bool discard = false;
    int max_idx = sysInput.Size()-1;
    long index;

    for (long k = 0; k < num_ffts; ++k)
    {
        temp_real_in.SetElements(0.0);
        temp_real_out.SetElements(0.0);
        index = k * (blockSize-num_overlap_points);

        for (long j = 0; j < blockSize; ++j)
        {
            temp_real_in(j) = sysInput(index);
            temp_real_out(j) = sysOutput(index);

            if (index == max_idx)
            {
                if (j < blockSize-1)
                {
                    discard = true;
                    break;
                }
            }

            ++index;
        }

        if (discard == true)
        {
            break;
        }

        fft.Compute(temp_real_in, input_spec);
        fft.Compute(temp_real_out, output_spec);

        input_spec.Conjugate();

        for (long j = 0; j < fftSize; ++j)
        {
            cross_spec.z(j) += input_spec.z(j) * output_spec.z(j);
            input_mag(j) += input_spec.z(j).MagSq();
            output_mag(j) += output_spec.z(j).MagSq();
        }
    }

    double numer;
    double denom;

    for (long k = 0; k < fftSize; ++k)
    {
        numer = cross_spec.z(k).MagSq();
        denom = input_mag(k) * output_mag(k);

        cohere(k) = (denom == 0.0) ? 1.0 : numer / denom;
    }

    return status;
}
//------------------------------------------------------------------------------
// Power spectral density of a real signal 
//------------------------------------------------------------------------------
hwMathStatus PSD(const hwMatrix& signal, 
                 double          sampFreq,
                 hwMatrix&       density, 
                 int             fftSize)
{
    int numPnts = signal.Size();
    
    if (!signal.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!signal.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (fftSize < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 4);
    }

    hwMathStatus status = density.Dimension(numPnts, hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(3);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    hwPSD psd(sampFreq, fftSize);

    status = psd.Compute(signal, density);
    if (!status.IsOk() && status.GetArg1() == 2)
    {
        status.SetArg1(3);
    }

    return status;
}
//------------------------------------------------------------------------------
// Cross power spectral density of a real signal pair
//------------------------------------------------------------------------------
hwMathStatus CPSD(const hwMatrix& signal1, 
                  const hwMatrix& signal2,
                  double          sampFreq,
                  hwMatrix&       density,
                  int             fftSize)
{
    if (!signal1.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!signal1.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    
    if (!signal2.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!signal2.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    if (signal2.Size() != signal1.Size())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    if (fftSize < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 5);
    }

    int numPnts = signal1.Size();
    hwMathStatus status = density.Dimension(numPnts, hwMatrix::REAL);

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

    hwCSD csd(sampFreq, fftSize);

    status = csd.Compute(signal1, signal2, density);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 3)
        {
            status.SetArg1(4);
        }
        return status;
    }

    return status;
}
//------------------------------------------------------------------------------
// Block power spectral density of a real signal with a window function
//------------------------------------------------------------------------------
hwMathStatus BlockPSD(const hwMatrix& signal, 
                      const hwMatrix& window,
                      int             num_overlap_points, 
                      double          sampFreq,
                      hwMatrix&       density, 
                      int             fftSize)
{
    if (!signal.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!signal.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (!window.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!window.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }

    int blockSize = window.Size();
    hwWindowFunc winFunc(window);

    int numPnts = signal.Size();

    if (blockSize > numPnts)
    {
        return hwMathStatus(HW_MATH_ERR_FTBLOCKSIZE, 1, 2);
    }
    if (num_overlap_points >= blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_OVERLAPPOINTS, 3);
    }
    if (num_overlap_points < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 3);
    }
    if (sampFreq <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 4);
    }

    if (fftSize == 0)
    {
        fftSize = blockSize;
    }
    else if (fftSize < blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_FFTSIZE, 2, 6);
    }

    int num_ffts = (int)floor((double)(numPnts - blockSize) /
               (double)(blockSize - num_overlap_points)) + 1;

    if (num_ffts < 1)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    hwMathStatus status = density.Dimension(fftSize, hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(5);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    hwPSD    psd(sampFreq, fftSize);
    hwMatrix sum(fftSize, hwMatrix::REAL);
    sum.SetElements(0.0);

    int index;
    double value;
    hwMatrix temp(fftSize, hwMatrix::REAL);
    hwMatrix resp(fftSize, hwMatrix::REAL);

    for (int k = 0; k < num_ffts; ++k)
    {
        temp.SetElements(0.0);
        index = k * (blockSize-num_overlap_points);

        winFunc.ApplyWindow(signal, index, temp, true);

        status = psd.Compute(temp, resp);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        sum += resp;
    }

    // compute output
    value = 1.0 / (double) num_ffts;
    density = sum * value;

    return status;
}
//------------------------------------------------------------------------------
// Block power spectral density of a real signal
//------------------------------------------------------------------------------
hwMathStatus BlockPSD(const hwMatrix& signal, 
                      int             blockSize,
                      int             num_overlap_points, 
                      double          sampFreq,
                      hwMatrix&       density, 
                      int             fftSize)
{
    if (!signal.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!signal.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    int numPnts = signal.Size();
    if (blockSize > numPnts)
    {
        return hwMathStatus(HW_MATH_ERR_FTBLOCKSIZE, 1, 2);
    }
    if (blockSize < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 2);
    }

    if (num_overlap_points >= blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_OVERLAPPOINTS, 3);
    }
    if (num_overlap_points < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 3);
    }

    if (sampFreq <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 4);
    }

    if (fftSize == 0)
    {
        fftSize = blockSize;
    }
    else if (fftSize < blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_FFTSIZE, 2, 6);
    }  

    int num_ffts = (int)floor((double)(numPnts - blockSize) /
                   (double)(blockSize - num_overlap_points)) + 1;
    if (num_ffts < 1)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    hwMathStatus status = density.Dimension(fftSize, hwMatrix::REAL);
    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(5);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    hwPSD    psd(sampFreq, fftSize);
    hwMatrix sum(fftSize, hwMatrix::REAL);
    sum.SetElements(0.0);

    bool discard = false;
    int max_idx = signal.Size()-1;
    int index;
    hwMatrix temp(fftSize, hwMatrix::REAL);
    hwMatrix resp(fftSize, hwMatrix::REAL);

    for (int k = 0; k < num_ffts; ++k)
    {
        temp.SetElements(0.0);
        index = k * (blockSize-num_overlap_points);

        for (int j = 0; j < blockSize; ++j)
        {
            temp(j) = signal(index);

            if (index == max_idx)
            {
                if (j < blockSize-1)
                {
                    discard = true;
                    break;
                }
            }

            ++index;
        }

        if (discard == true)
            break;

        status = psd.Compute(temp, resp);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        sum += resp;
    }

    // compute output
    double value = 1.0 / (double) num_ffts;
    density = sum * value;

    return status;
}
//------------------------------------------------------------------------------
// Block crosss power spectral density of a real signal pair with a window function
//------------------------------------------------------------------------------
hwMathStatus BlockCPSD(const hwMatrix& signal1, 
                       const hwMatrix& signal2,
                       const hwMatrix& window, 
                       int             num_overlap_points,
                       double          sampFreq, 
                       hwMatrix&       density, 
                       int             fftSize)
{
    if (!signal1.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!signal1.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    
    if (!signal2.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!signal2.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    if (signal2.Size() != signal1.Size())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    if (!window.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);
    }
    if (!window.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 3);
    }

    int numPnts   = signal1.Size();
    int blockSize = window.Size();
    hwWindowFunc winFunc(window);

    if (blockSize > numPnts)
    {
        return hwMathStatus(HW_MATH_ERR_FTBLOCKSIZE, 1, 3);
    }

    if (num_overlap_points >= blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_OVERLAPPOINTS, 4);
    }
    if (num_overlap_points < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 4);
    }
    if (sampFreq <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 5);
    }

    if (fftSize == 0)
    {
        fftSize = blockSize;
    }
    else if (fftSize < blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_FFTSIZE, 3, 7);
    }

    int index;
    hwMatrix temp1(fftSize, hwMatrix::REAL);
    hwMatrix temp2(fftSize, hwMatrix::REAL);
    hwMatrix resp;

    int num_ffts = (int)floor((double)(numPnts - blockSize) /
                   (double)(blockSize - num_overlap_points)) + 1;
    if (num_ffts < 1)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 3);
    }

    hwMathStatus status = density.Dimension(fftSize, hwMatrix::REAL);
    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(6);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    hwCSD cpsd(sampFreq, fftSize);
    hwMatrix sum(fftSize, hwMatrix::REAL);
    sum.SetElements(0.0);

    for (int k = 0; k < num_ffts; ++k)
    {
        temp1.SetElements(0.0);
        temp2.SetElements(0.0);
        index = k * (blockSize-num_overlap_points);

        winFunc.ApplyWindow(signal1, index, temp1,  true);
        winFunc.ApplyWindow(signal2, index, temp2,  true);

        status = cpsd.Compute(temp1, temp2, resp);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        sum += resp;
    }

    // compute output
    double value = 1.0 / (double) num_ffts;
    density = sum * value;

    return status;
}
//------------------------------------------------------------------------------
// Block crosss power spectral density of a real signal pair
//------------------------------------------------------------------------------
hwMathStatus BlockCPSD(const hwMatrix& signal1, 
                       const hwMatrix& signal2,
                       int             blockSize, 
                       int             num_overlap_points,
                       double          sampFreq, 
                       hwMatrix&       density, 
                       int             fftSize)
{   
    if (!signal1.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!signal1.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    
    if (!signal2.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!signal2.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }
    if (signal2.Size() != signal1.Size())
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    int numPnts = signal1.Size();
    if (blockSize > numPnts)
    {
        return hwMathStatus(HW_MATH_ERR_FTBLOCKSIZE, 1, 3);
    }
    if (blockSize < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 3);
    }
    if (num_overlap_points >= blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_OVERLAPPOINTS, 4);
    }
    if (num_overlap_points < 0)
    {
        return hwMathStatus(HW_MATH_ERR_NONNONNEGINT, 4);
    }
    if (sampFreq <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 5);
    }

    if (fftSize == 0)
    {
        fftSize = blockSize;
    }
    else if (fftSize < blockSize)
    {
        return hwMathStatus(HW_MATH_ERR_FFTSIZE, 3, 7);
    }

    bool discard = false;
    int max_idx = numPnts-1;
    int index;
    hwMatrix temp_in1(fftSize, hwMatrix::REAL);
    hwMatrix temp_in2(fftSize, hwMatrix::REAL);
    hwMatrix temp_out;

    int num_ffts = (int)floor((double)(numPnts - blockSize) /
                   (double)(blockSize - num_overlap_points)) + 1;
    if (num_ffts < 1)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 3);
    }

    hwMathStatus status = density.Dimension(fftSize, hwMatrix::REAL);
    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(6);
        }
        else
        {
            status.ResetArgs();
        }
        return status;
    }

    hwCSD cpsd(sampFreq, fftSize);
    hwMatrix sum(fftSize, hwMatrix::REAL);
    sum.SetElements(0.0);

    for (int k = 0; k < num_ffts; ++k)
    {
        temp_in1.SetElements(0.0);
        temp_in2.SetElements(0.0);
        index = k * (blockSize-num_overlap_points);

        for (int j = 0; j < blockSize; ++j)
        {
            temp_in1(j) = signal1(index);
            temp_in2(j) = signal2(index);

            if (index == max_idx)
            {
                if (j < blockSize-1)
                {
                    discard = true;
                    break;
                }
            }

            ++index;
        }

        if (discard == true)
            break;

        status = cpsd.Compute(temp_in1, temp_in2, temp_out);

        if (!status.IsOk())
        {
            if (status.GetArg1() == 3)
                status.SetArg1(6);
            else
                status.ResetArgs();

            return status;
        }

        sum += temp_out;
    }

    // compute output
    double value = 1.0 / (double) num_ffts;
    density = sum * value;

    return status;
}
