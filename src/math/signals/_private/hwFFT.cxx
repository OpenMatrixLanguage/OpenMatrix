/**
* @file  hwFFT.cxx
* @date July 2011
* Copyright (C) 2011-2018 Altair Engineering, Inc.  
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
#include "hwFFT.h"

#include <math.h>
#include "hwMatrix.h"
#include "hwMatrixN.h"

// store FFTW plans and info for reuse
bool initializedFFTW    = false;
fftw_plan file_planR2C  = NULL;
fftw_plan file_planC2R  = NULL;
fftw_plan file_planC2C  = NULL;
int file_plan_fftsize   = -1;
unsigned file_plan_dir  = 0;   // only affects FFT_C2C

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwFFTW::hwFFTW(unsigned direction, int fftSize)
    : m_type          (FFT_NOTSET)
    , m_realInput     (nullptr)
    , m_realOutput    (nullptr)
    , m_complexInput  (nullptr)
    , m_complexOutput (nullptr)
    , m_fftSize       (-1)
{
    if (fftSize < 0)
    {
        m_status(HW_MATH_ERR_NONNONNEGINT, 2);
    }

    m_fftSize = fftSize;

    if (file_plan_dir != direction)
    {
        if (file_planC2C)
        {
            fftw_destroy_plan(file_planC2C);
            file_planC2C = NULL;
        }

        file_plan_dir = direction;
    }

    if (!initializedFFTW)
    {
        int nthreads = 6;
        int val      = fftw_init_threads();
        fftw_plan_with_nthreads(nthreads);
        initializedFFTW = true;
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwFFTW::~hwFFTW()
{
    if (m_type == FFT_R2C)
    {
        // only m_direction == FFTW_FORWARD for now
        if (m_fftSize > m_dataSize)
        {
            if (m_realInput)
            {
                delete [] m_realInput;
            }
        }
    }
    else if (m_type == FFT_C2R)
    {
        // only m_direction == FFTW_BACKWARD for now
        if (m_fftSize != m_dataSize)
        {
            if (m_complexInput)
            {
                delete [] m_complexInput;
            }
        }
    }
    else if (m_type == FFT_C2C)
    {
        if (file_plan_dir == FFTW_FORWARD)
        {
            if (m_fftSize > m_dataSize)
            {
                if (m_complexInput)
                {
                    delete [] m_complexInput;
                }
            }
        }
        else    // m_direction == FFTW_BACKWARD
        {
            if (m_fftSize != m_dataSize)
            {
                if (m_complexInput)
                {
                    delete [] m_complexInput;
                }
            }
        }
    }
}
//------------------------------------------------------------------------------
// Determine FFTW plan type
//------------------------------------------------------------------------------
void hwFFTW::DetermineType(const hwMatrixN& input)
{
    if (input.IsReal())
    {
        m_type = FFT_R2C;
    }
    else
    {
        m_type = FFT_C2C;
    }
}
//------------------------------------------------------------------------------
// Allocated memory based on FFT type
//------------------------------------------------------------------------------
void hwFFTW::SetupFFT(const hwMatrix& input)
{
    if (m_fftSize < 0)
    {
        m_dataSize = 0;
        m_status(HW_MATH_ERR_NONNONNEGINT);
        return;
    }

    m_dataSize = input.Size();

    if (m_fftSize == 0)
    {
        m_fftSize = m_dataSize;
    }

    // validate size of file scope FFTW plans
    if (file_plan_fftsize != m_fftSize)
    {
        if (file_planR2C)
        {
            fftw_destroy_plan(file_planR2C);
            file_planR2C = NULL;
        }

        if (file_planC2R)
        {
            fftw_destroy_plan(file_planC2R);
            file_planC2R = NULL;
        }

        if (file_planC2C)
        {
            fftw_destroy_plan(file_planC2C);
            file_planC2C = NULL;
        }

        fftw_cleanup_threads();
        file_plan_fftsize = m_fftSize;
    }

    DetermineType(input);

    if (m_type == FFT_R2C)
    {
        // only m_direction == FFTW_FORWARD for now
        if (m_fftSize > m_dataSize)
        {
            try
            {
                m_realInput = new double[m_fftSize];
            }
            catch (std::bad_alloc)
            {
                m_status(HW_MATH_ERR_OUTOFMEMORY);
                return;
            }
        }
    }
    else if (m_type == FFT_C2R)
    {
        // only m_direction == FFTW_BACKWARD for now
        if (m_fftSize != m_dataSize)
        {
            int n2 = m_fftSize/2 + 1;

            try
            {
                m_complexInput = new fftw_complex[n2];
            }
            catch (std::bad_alloc)
            {
                m_status(HW_MATH_ERR_OUTOFMEMORY);
                return;
            }
        }
    }
    else if (m_type == FFT_C2C)
    {
        if (file_plan_dir == FFTW_FORWARD)
        {
            if (m_fftSize > m_dataSize)
            {
                try
                {
                    m_complexInput = new fftw_complex[m_fftSize];
                }
                catch (std::bad_alloc)
                {
                    m_status(HW_MATH_ERR_OUTOFMEMORY);
                    return;
                }
            }
        }
        else    // m_direction == FFTW_BACKWARD
        {
            if (m_fftSize != m_dataSize)
            {
                try
                {
                    m_complexInput = new fftw_complex[m_fftSize];
                }
                catch (std::bad_alloc)
                {
                    m_status(HW_MATH_ERR_OUTOFMEMORY);
                    return;
                }
            }
        }
    }
}
//------------------------------------------------------------------------------
// Call FFTW for const input
//------------------------------------------------------------------------------
hwMathStatus hwFFTW::Compute(const hwMatrix& input, hwMatrix& output)
{
    if (m_type == FFT_NOTSET)
    {
        SetupFFT(input);
    }
    if (!m_status.IsOk())
    {
        return m_status;
    }

    unsigned options;

    if (m_type == FFT_R2C || m_type == FFT_C2C)
    {
        if (m_fftSize > m_dataSize)
        {
            options = FFTW_ESTIMATE;
        }
        else
        {
            options = FFTW_ESTIMATE|FFTW_PRESERVE_INPUT;
        }
    }
    else if (m_type == FFT_C2R)
    {
        if (m_dataSize != m_fftSize)
        {
            options = FFTW_ESTIMATE;
        }
        else
        {
            options = FFTW_ESTIMATE|FFTW_PRESERVE_INPUT;
        }
    }
    else
    {
        if (m_dataSize != 0)
        {
            return m_status(HW_MATH_ERR_NOTALLOWED);
        }
        return m_status;     // empty input
    }

    hwMatrix* input_nc = (hwMatrix*) (&input);      // override const with preserved input
    return Compute(*input_nc, output, options);
}
//------------------------------------------------------------------------------
// Call FFTW for non-const input
//------------------------------------------------------------------------------
hwMathStatus hwFFTW::Compute(hwMatrix& input, hwMatrix& output)
{
    if (m_type == FFT_NOTSET)
    {
        SetupFFT(input);
    }
    if (!m_status.IsOk())
    {
        return m_status;
    }
    unsigned options = FFTW_ESTIMATE;

    return Compute(input, output, options);
}
//------------------------------------------------------------------------------
//! Call FFTW, selecting the appropriate plan
//------------------------------------------------------------------------------
hwMathStatus hwFFTW::Compute(hwMatrix&       input, 
                             hwMatrix&       output, 
                             const unsigned& options)
{
    if (!input.IsEmptyOrVector())
    {
        return m_status(HW_MATH_ERR_VECTOR, 1);
    }

    if (input.Size() != m_dataSize)
    {
        return m_status(HW_MATH_ERR_ARRAYSIZE, 1);
    }

    // dimension output
    if (input.M() == 1)
    {
        switch (m_type)
        {
            case FFT_R2C:
                m_status = output.Dimension(1, m_fftSize, hwMatrix::COMPLEX);
                break;
            case FFT_C2R:
                m_status = output.Dimension(1, m_fftSize, hwMatrix::REAL);
                break;
            case FFT_C2C:
                m_status = output.Dimension(1, m_fftSize, hwMatrix::COMPLEX);
                break;
            default: break;
        }
    }
    else if (input.N() == 1)
    {
        switch (m_type)
        {
            case FFT_R2C:  
                m_status = output.Dimension(m_fftSize, 1, hwMatrix::COMPLEX);
                break;
            case FFT_C2R:
                m_status = output.Dimension(m_fftSize, 1, hwMatrix::REAL);
                break;
            case FFT_C2C:
                m_status = output.Dimension(m_fftSize, 1, hwMatrix::COMPLEX);
                break;
            default: break;
        }
    }
    else if (input.M() == 0)
    {
        int size = (m_fftSize == 0) ? input.N() : m_fftSize;

        switch (m_type)
        {
            case FFT_R2C:
                m_status = output.Dimension(0, size, hwMatrix::COMPLEX);
                break;
            case FFT_C2R:
                m_status = output.Dimension(0, size, hwMatrix::REAL);
                break;
            case FFT_C2C:
                m_status = output.Dimension(0, size, hwMatrix::COMPLEX);
                break;
            default: break;
        }
    }
    else // input.N() = 0
    {
        int size = (m_fftSize == 0) ? input.M() : m_fftSize;

        switch (m_type)
        {
            case FFT_R2C:
                m_status = output.Dimension(size, 0, hwMatrix::COMPLEX);
                break;
            case FFT_C2R:
                m_status = output.Dimension(size, 0, hwMatrix::REAL);
                break;
            case FFT_C2C:
                m_status = output.Dimension(size, 0, hwMatrix::COMPLEX);
                break;
            default: break;
        }
    }

    if (!m_status.IsOk())
    {
        m_status.SetArg1(2);
        return m_status;
    }

    if (output.IsEmpty())
    {
        return m_status;
    }

    // prepare FFTW plan based on the data type
    if (m_type == FFT_R2C)
    {
        // only m_direction == FFTW_FORWARD for now
        if (!input.IsReal())
        {
            return m_status(HW_MATH_ERR_ARRAYTYPE, 1);
        }

        if (m_fftSize <= m_dataSize)
        {
            m_realInput = input.GetRealData();
        }
        else
        {
            const double* tempData = input.GetRealData();

            memcpy_s(m_realInput, m_dataSize * sizeof(double), tempData, m_dataSize * sizeof(double));
            memset(m_realInput + m_dataSize, 0, (m_fftSize-m_dataSize) * sizeof(double));
        }

        m_complexOutput = (fftw_complex*) output.GetComplexData();

        if (!file_planR2C)
        {
            file_planR2C = fftw_plan_dft_r2c_1d(m_fftSize, m_realInput, m_complexOutput, options);
        }

        fftw_execute_dft_r2c(file_planR2C, m_realInput, m_complexOutput);
    }
    else if (m_type == FFT_C2R)
    {
        // only m_direction == FFTW_BACKWARD for now
        if (input.IsReal())
        {
            return m_status(HW_MATH_ERR_ARRAYTYPE, 1);
        }

        int n2;
        int pad = 0;
        
        if (m_fftSize == m_dataSize)
        {
            m_complexInput = (fftw_complex*) input.GetComplexData();
        }
        else if (m_fftSize > m_dataSize)
        {
            // pad symmetrically about the Nyquist frequency, not at the end
            const double* tempData = (const double*) input.GetComplexData();

            if (m_dataSize%2 == 0)
            {
                // copy data
                n2 = m_dataSize/2;
                memcpy_s(m_complexInput, 2*n2 * sizeof(double), tempData, 2*n2 * sizeof(double));
                // split Nyquist component
                m_complexInput[n2][0] = 0.5 * tempData[2*n2];
                m_complexInput[n2][1] = 0.5 * tempData[2*n2+1];
                // pad zeros    
                pad = m_fftSize/2 - n2;
                memset(m_complexInput + (n2+1), 0, 2*pad * sizeof(double));
            }
            else
            {
                // copy data
                n2 = (m_dataSize+1) / 2;
                memcpy_s(m_complexInput, 2*n2 * sizeof(double), tempData, 2*n2 * sizeof(double));
                // pad zeros    
                pad = m_fftSize/2 + 1 - n2;
                memset(m_complexInput + n2, 0, 2*pad * sizeof(double));
            }
        }
        else    // m_fftSize < m_dataSize
        {
            // truncate symmetrically about the Nyquist frequency, not at the end
            const double* tempData = (const double*) input.GetComplexData();

            if (m_fftSize%2 == 0)
            {
                // copy data
                n2 = m_fftSize/2;
                memcpy_s(m_complexInput, 2*n2 * sizeof(double), tempData, 2*n2 * sizeof(double));
                // add symmetric components to produce new Nyquist component. This works because the
                // imaginary component of a sampled Nyquist frequency phasor is always zero.
                m_complexInput[n2][0] = tempData[2*n2] + tempData[2*(m_dataSize-n2)];
                m_complexInput[n2][1] = 0.0;    // = tempData[2*n2+1] + tempData[2*(m_dataSize-n2)+1];
            }
            else
            {
                // copy data
                n2 = (m_fftSize+1)/2;
                memcpy_s(m_complexInput, 2*n2 * sizeof(double), tempData, 2*n2 * sizeof(double));
            }
        }

        m_realOutput = output.GetRealData();

        if (!file_planC2R)
        {
            file_planC2R = fftw_plan_dft_c2r_1d(m_fftSize, m_complexInput, m_realOutput, options);
        }

        fftw_execute_dft_c2r(file_planC2R, m_complexInput, m_realOutput);
    }
    else if (m_type == FFT_C2C)
    {
        if (input.IsReal())
        {
            return m_status(HW_MATH_ERR_ARRAYTYPE, 1);
        }
        if (file_plan_dir == FFTW_FORWARD)
        {
            if (m_fftSize <= m_dataSize)
            {
                m_complexInput = (fftw_complex*) input.GetComplexData();
            }
            else
            {
                const double* tempData = (double*) input.GetComplexData();
                // copy data
                memcpy_s(m_complexInput, 2*m_dataSize * sizeof(double), tempData, 2*m_dataSize * sizeof(double));
                // pad zeros    
                memset(m_complexInput + m_dataSize, 0, 2*(m_fftSize-m_dataSize) * sizeof(double));
            }

            m_complexOutput = (fftw_complex*) output.GetComplexData();
        }
        else    // m_direction == FFTW_BACKWARD
        {
            if (m_fftSize == m_dataSize)
            {
                m_complexInput = (fftw_complex*) input.GetComplexData();
            }
            else if (m_fftSize > m_dataSize)
            {
                // pad symmetrically about the Nyquist frequency, not at the end
                const double* tempData = (const double*) input.GetComplexData();
                int n2;
                int pad;

                if (m_dataSize%2 == 0)
                {
                    // copy data
                    n2 = m_dataSize/2;
                    memcpy_s(m_complexInput, 2*n2 * sizeof(double), tempData, 2*n2 * sizeof(double));
                    // split Nyquist component
                    m_complexInput[n2][0] = 0.5 * tempData[2*n2];
                    m_complexInput[n2][1] = 0.5 * tempData[2*n2+1];
                    // pad zeros    
                    pad = m_fftSize - m_dataSize - 1;
                    memset(m_complexInput + (n2+1), 0, 2*pad * sizeof(double));
                    // split Nyquist component
                    m_complexInput[n2+pad+1][0] = 0.5 * tempData[2*n2];
                    m_complexInput[n2+pad+1][1] = 0.5 * tempData[2*n2+1];
                    // copy data
                    memcpy_s(m_complexInput + (n2+pad+2), 2*(n2-1) * sizeof(double), tempData + 2*(n2+1), 2*(n2-1) * sizeof(double));
                }
                else
                {
                    // copy data
                    int n2 = (m_dataSize+1) / 2;
                    memcpy_s(m_complexInput, 2*n2 * sizeof(double), tempData, 2*n2 * sizeof(double));
                    // pad zeros
                    pad = m_fftSize - m_dataSize;
                    memset(m_complexInput + n2, 0, 2*pad * sizeof(double));
                    // copy data
                    memcpy_s(m_complexInput + (n2+pad), 2*(n2-1) * sizeof(double), tempData + 2*n2, 2*(n2-1) * sizeof(double));
                }
            }
            else    // m_fftSize < m_dataSize
            {
                // truncate symmetrically about the Nyquist frequency, not at the end
                const double* tempData = (const double*) input.GetComplexData();
                int n2;
                int omit;

                if (m_fftSize%2 == 0)
                {
                    // copy data
                    n2 = m_fftSize/2;
                    memcpy_s(m_complexInput, 2*n2 * sizeof(double), tempData, 2*n2 * sizeof(double));
                    // add symmetric components to produce new Nyquist component. This works because the
                    // imaginary component of a sampled Nyquist frequency phasor is always zero for a real
                    // signal, and real component of a sampled Nyquist frequency phasor is always zero for a
                    // pure imaginary signal. These two properties can be combined to produce the following.
                    m_complexInput[n2][0] = tempData[2*n2] + tempData[2*(m_dataSize-n2)];
                    m_complexInput[n2][1] = tempData[2*n2+1] + tempData[2*(m_dataSize-n2)+1];
                    omit = m_dataSize - (2*n2-1);
                    memcpy_s(m_complexInput + n2+1, 2*n2 * sizeof(double), tempData + 2*(n2+omit), 2*(n2-1) * sizeof(double));
                }
                else
                {
                    n2 = (m_fftSize+1)/2;
                    memcpy_s(m_complexInput, 2*n2 * sizeof(double), tempData, 2*n2 * sizeof(double));
                    omit = m_dataSize - (2*n2-1);
                    memcpy_s(m_complexInput + n2, 2*(n2-1) * sizeof(double), tempData + 2*(n2+omit), 2*(n2-1) * sizeof(double));
                }
            }

            m_complexOutput = (fftw_complex*) output.GetComplexData();
        }

        if (!file_planC2C)
        {
            file_planC2C = fftw_plan_dft_1d(m_fftSize, m_complexInput, m_complexOutput, file_plan_dir, options);
        }

        fftw_execute_dft(file_planC2C, m_complexInput, m_complexOutput);
    }

    // prepare output
    if (m_type == FFT_R2C)
    {
        // populate conjugate symmetrix half of spectrum
        int n2 = (m_fftSize+1) / 2;

        for (int i = 1; i < n2; ++i) 
        {
            m_complexOutput[m_fftSize-i][0] =  m_complexOutput[i][0];
            m_complexOutput[m_fftSize-i][1] = -m_complexOutput[i][1];
        }
    }
    else if (m_type == FFT_C2R)
    {
        double value = 1.0 / (double) m_dataSize;
        
        for (int i = 0; i < m_fftSize; ++i) 
        {
            m_realOutput[i] *= value;
        }
    }
    else if (m_type == FFT_C2C)
    {
        if (file_plan_dir == FFTW_FORWARD)
        {
            // nothing to do
        }
        else if (file_plan_dir == FFTW_BACKWARD)
        {
            double value = 1.0 / (double) m_dataSize;

            for (int i = 0; i < m_fftSize; ++i) 
            {
                m_complexOutput[i][0] *= value;
                m_complexOutput[i][1] *= value;
            }
        }
    }

    return m_status;
}
//------------------------------------------------------------------------------
// Call FFTW for a specified dimension for const input and return status
//------------------------------------------------------------------------------
hwMathStatus hwFFTW::Compute(const hwMatrix& input, int dim, hwMatrix& output)
{
    // clean up
    if (file_planR2C)
    {
        fftw_destroy_plan(file_planR2C);
        file_planR2C = NULL;
    }

    if (file_planC2R)
    {
        fftw_destroy_plan(file_planC2R);
        file_planC2R = NULL;
    }

    if (file_planC2C)
    {
        fftw_destroy_plan(file_planC2C);
        file_planC2C = NULL;
    }

    fftw_cleanup_threads();
    file_plan_fftsize = -1;
    file_plan_dir = 0;

    // start over
    int rank = 1;
    int howmany;
    int stride;
    int dist;
    unsigned options = FFTW_ESTIMATE;

    if (dim == 0)       // FFT of columns
    {
        stride = 1;
        howmany = input.N();
        dist = input.M();
        m_fftSize = input.M();
    }
    else if (dim == 1)  // FFT of rows
    {
        stride = input.M();
        howmany = input.M();
        dist = 1;
        m_fftSize = input.N();
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    m_dataSize = m_fftSize;

    DetermineType(input);

    if (file_plan_dir == FFTW_FORWARD)
    {
        DetermineType(input);
    }
    else
    {
        DetermineType(input, howmany, dist, stride);
    }

    if (m_type == FFT_R2C)
    {
        m_realInput = const_cast<double*> (input.GetRealData());
        hwMathStatus status = output.Dimension(input.M(), input.N(), hwMatrix::COMPLEX);
        m_complexOutput = reinterpret_cast<fftw_complex*> (output.GetComplexData());

        file_planR2C = fftw_plan_many_dft_r2c(rank, &m_fftSize, howmany,
            m_realInput, 0, stride, dist, m_complexOutput, 0, stride, dist, options);

        fftw_execute_dft_r2c(file_planR2C, m_realInput, m_complexOutput);

        // populate conjugate symmetrix half of spectrum
        int n2 = (m_fftSize + 1) / 2;

        for (int j = 0; j < howmany; ++j)
        {
            for (int k = stride; k < n2 * stride; k += stride)
            {
                m_complexOutput[stride * m_fftSize - k][0] = m_complexOutput[k][0];
                m_complexOutput[stride * m_fftSize - k][1] = -m_complexOutput[k][1];
            }

            m_complexOutput += dist;
        }

        fftw_destroy_plan(file_planR2C);
        file_planR2C = NULL;
    }
    else if (m_type == FFT_C2R)
    {
        m_complexInput = (fftw_complex*)input.GetComplexData();
        output.Dimension(input.M(), input.N(), hwMatrix::REAL);
        m_realOutput = output.GetRealData();

        file_planC2R = fftw_plan_many_dft_c2r(rank, &m_fftSize, howmany,
            m_complexInput, 0, stride, dist, m_realOutput, 0, stride, dist, options);

        fftw_execute_dft_c2r(file_planC2R, m_complexInput, m_realOutput);

        for (int k = 0; k < howmany * m_fftSize; ++k)
        {
            m_realOutput[k] /= m_fftSize;
        }

        fftw_destroy_plan(file_planC2R);
        file_planC2R = NULL;
    }
    else  // m_type == FFT_C2C
    {
        m_complexInput = (fftw_complex*)input.GetComplexData();
        hwMathStatus status = output.Dimension(input.M(), input.N(), hwMatrix::COMPLEX);
        m_complexOutput = reinterpret_cast<fftw_complex*> (output.GetComplexData());

        file_planC2C = fftw_plan_many_dft(rank, &m_fftSize, howmany,
            m_complexInput, 0, stride, dist, m_complexOutput, 0, stride, dist, file_plan_dir, options);

        fftw_execute_dft(file_planC2C, m_complexInput, m_complexOutput);

        if (file_plan_dir == FFTW_BACKWARD)
        {
            for (int k = 0; k < howmany * m_fftSize; ++k)
            {
                m_complexOutput[k][0] /= m_fftSize;
                m_complexOutput[k][1] /= m_fftSize;
            }
        }

        fftw_destroy_plan(file_planC2C);
        file_planC2C = NULL;
    }

    fftw_cleanup_threads();
        
    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Call FFTW for a specified dimension for const input and return status
//------------------------------------------------------------------------------
hwMathStatus hwFFTW::Compute(const hwMatrixN& input, int dim, hwMatrixN& output)
{
    // clean up
    if (file_planR2C)
    {
        fftw_destroy_plan(file_planR2C);
        file_planR2C = NULL;
    }

    if (file_planC2R)
    {
        fftw_destroy_plan(file_planC2R);
        file_planC2R = NULL;
    }

    if (file_planC2C)
    {
        fftw_destroy_plan(file_planC2C);
        file_planC2C = NULL;
    }

    fftw_cleanup_threads();
    file_plan_fftsize = -1;

    // start over
    unsigned options = FFTW_ESTIMATE;
    m_fftSize = input.Dimensions()[dim];
    std::vector<int> index(input.Dimensions().size());
    ++index[dim];
    int stride = input.Index(index);
    int rank = 1;
    int howmany;                            // number of vector per stride group
    int dist;                               // distance between start of stride groups
    int numVecs = input.Size() / m_fftSize; // number of vectors in the dimension
    int numGroups;                          // number of stride groups

    if (stride == 1)
    {
        numGroups = 1;
        howmany = numVecs;
        dist = m_fftSize;
    }
    else
    {
        numGroups = numVecs / stride;
        howmany = stride;
        dist = 1;
    }

    m_dataSize = m_fftSize;

    if (file_plan_dir == FFTW_FORWARD)
    {
        DetermineType(input);
    }
    else
    {
        DetermineType(input, numGroups, howmany, dist, stride);
    }

    if (m_type == FFT_R2C)
    {
        m_realInput = const_cast<double*> (input.GetRealData());
        output.Dimension(input.Dimensions(), hwMatrixN::COMPLEX);
        m_complexOutput = reinterpret_cast<fftw_complex*> (output.GetComplexData());

        file_planR2C = fftw_plan_many_dft_r2c(rank, &m_fftSize, howmany,
            m_realInput, 0, stride, dist, m_complexOutput, 0, stride, dist, options);

        for (int i = 0; i < numGroups; ++i)
        {
            fftw_execute_dft_r2c(file_planR2C, m_realInput, m_complexOutput);
            m_realInput += stride * m_fftSize;  // advance to next stride group
            m_complexOutput += stride * m_fftSize;
        }

        // populate conjugate symmetrix half of spectrum
        int n2 = (m_fftSize + 1) / 2;

        m_complexOutput = reinterpret_cast<fftw_complex*> (output.GetComplexData());

        for (int i = 0; i < numGroups; ++i)
        {
            for (int j = 0; j < howmany; ++j)
            {
                for (int k = stride; k < n2 * stride; k += stride)
                {
                    m_complexOutput[stride * m_fftSize - k][0] = m_complexOutput[k][0];
                    m_complexOutput[stride * m_fftSize - k][1] = -m_complexOutput[k][1];
                }

                m_complexOutput += dist;
            }

            m_complexOutput += stride * m_fftSize - howmany * dist;
        }

        fftw_destroy_plan(file_planR2C);
        file_planR2C = NULL;
    }
    else if (m_type == FFT_C2R)
    {
        m_complexInput = (fftw_complex*)input.GetComplexData();
        output.Dimension(input.Dimensions(), hwMatrixN::REAL);
        m_realOutput = output.GetRealData();

        file_planC2R = fftw_plan_many_dft_c2r(rank, &m_fftSize, howmany,
            m_complexInput, 0, stride, dist, m_realOutput, 0, stride, dist, options);

        for (int i = 0; i < numGroups; ++i)
        {
            fftw_execute_dft_c2r(file_planC2R, m_complexInput, m_realOutput);

            for (int k = 0; k < howmany * m_fftSize; ++k)
            {
                m_realOutput[k] /= m_fftSize;
            }

            m_complexInput += stride * m_fftSize;  // advance to next stride group
            m_realOutput += stride * m_fftSize;
        }

        fftw_destroy_plan(file_planC2R);
        file_planC2R = NULL;
    }
    else  // m_type == FFT_C2C
    {
        m_complexInput = (fftw_complex*)input.GetComplexData();
        output.Dimension(input.Dimensions(), hwMatrixN::COMPLEX);
        m_complexOutput = reinterpret_cast<fftw_complex*> (output.GetComplexData());

        file_planC2C = fftw_plan_many_dft(rank, &m_fftSize, howmany,
            m_complexInput, 0, stride, dist, m_complexOutput, 0, stride, dist, file_plan_dir, options);

        for (int i = 0; i < numGroups; ++i)
        {
            fftw_execute_dft(file_planC2C, m_complexInput, m_complexOutput);

            if (file_plan_dir == FFTW_BACKWARD)
            {
                for (int k = 0; k < howmany * m_fftSize; ++k)
                {
                    m_complexOutput[k][0] /= m_fftSize;
                    m_complexOutput[k][1] /= m_fftSize;
                }
            }

            m_complexInput += stride * m_fftSize;  // advance to next stride group
            m_complexOutput += stride * m_fftSize;
        }

        fftw_destroy_plan(file_planC2C);
        file_planC2C = NULL;
    }

    fftw_cleanup_threads();

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwFFT_f::hwFFT_f(int fftSize)
    : hwFFTW(FFTW_FORWARD, fftSize)
{
    if (!m_status.IsOk())
    {
        m_status.SetArg1(1);
    }
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwFFT_f::~hwFFT_f()
{
    if (!m_status.IsOk())
    {
        m_status.SetArg1(1);
    }
}
//------------------------------------------------------------------------------
// Determine FFTW plan type
//------------------------------------------------------------------------------
void hwFFT_f::DetermineType(const hwMatrix& input)
{
    if (input.IsReal())
    {
        m_type = FFT_R2C;
    }
    else
    {
        m_type = FFT_C2C;
    }
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwFFT_r::hwFFT_r(int fftSize, bool assumeConjSym)
       : hwFFTW(FFTW_BACKWARD, fftSize)
{
    m_assumeSymmetry = assumeConjSym;
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwFFT_r::~hwFFT_r()
{
}
//------------------------------------------------------------------------------
// Determine FFTW plan type
//------------------------------------------------------------------------------
void hwFFT_r::DetermineType(const hwMatrix& input)
{
    bool conjugateSymmetric = true;

    if (!m_assumeSymmetry)
    {
        if (input.IsReal())
        {
            // This section is not in use. It would become case FFT_R2R.
            // The autocorrelation of a real signal could utilize it.
            const double* realInput = input.GetRealData();

            for (int i = 1; i < m_dataSize/2; ++i)
            {
                if (realInput[i] != realInput[m_dataSize-i])
                {
                    conjugateSymmetric = false;
                    break;
                }
            }
        }
        else
        {
            const hwComplex* complexInput = input.GetComplexData();

            if (m_dataSize)
            {
                if (complexInput[0].Imag() != 0.0)
                {
                    conjugateSymmetric = false;
                }
                else
                {
                    for (int i = 1; i < m_dataSize/2; ++i)
                    {
                        if (complexInput[i] != complexInput[m_dataSize-i].Conjugate())
                        {
                            conjugateSymmetric = false;
                            break;
                        }
                    }

                    if (conjugateSymmetric && m_dataSize % 2 == 0)
                    {
                        if (complexInput[m_dataSize / 2].Imag() != 0.0)
                        {
                            conjugateSymmetric = false;
                        }
                    }
                }
            }
        }
    }

    if (conjugateSymmetric)
    {
        m_type = FFT_C2R;
    }
    else
    {
        m_type = FFT_C2C;
    }
}
//------------------------------------------------------------------------------
// Determine FFTW plan type
//------------------------------------------------------------------------------
void hwFFT_r::DetermineType(const hwMatrix& input, int howmany,
                            int dist, int stride)
{
    bool conjugateSymmetric = true;

    if (!m_assumeSymmetry)
    {
        int n2 = (m_fftSize + 1) / 2;

        if (input.IsReal())
        {
            // This section is not in use. It would become case FFT_R2R.
            // The autocorrelation of a real signal could utilize it.
            const double* realInput = input.GetRealData();

            for (int j = 0; j < howmany; ++j)
            {
                for (int k = stride; k < n2 * stride; k += stride)
                {
                    if (realInput[k] != realInput[stride * m_fftSize - k])
                    {
                        conjugateSymmetric = false;
                        break;
                    }
                }

                if (!conjugateSymmetric)
                    break;

                realInput += dist;
            }

            realInput += stride * m_fftSize - howmany * dist;
        }
        else
        {
            const hwComplex* complexInput = input.GetComplexData();

            for (int j = 0; j < howmany; ++j)
            {
                if (complexInput[0].Imag() != 0.0)
                {
                    conjugateSymmetric = false;
                }
                else
                {
                    for (int k = stride; k < n2 * stride; k += stride)
                    {
                        if (complexInput[k] != complexInput[stride * m_fftSize - k].Conjugate())
                        {
                            conjugateSymmetric = false;
                            break;
                        }
                    }

                    complexInput += dist;
                }

                if (!conjugateSymmetric)
                    break;
            }

            complexInput += stride * m_fftSize - howmany * dist;
        }
    }

    if (conjugateSymmetric)
    {
        m_type = FFT_C2R;
    }
    else
    {
        m_type = FFT_C2C;
    }
}
//------------------------------------------------------------------------------
// Determine FFTW plan type
//------------------------------------------------------------------------------
void hwFFT_r::DetermineType(const hwMatrixN& input, int numGroups, int howmany,
                            int dist, int stride)
{
    bool conjugateSymmetric = true;

    if (!m_assumeSymmetry)
    {
        int n2 = (m_fftSize + 1) / 2;

        if (input.IsReal())
        {
            // This section is not in use. It would become case FFT_R2R.
            // The autocorrelation of a real signal could utilize it.
            const double* realInput = input.GetRealData();

            for (int i = 0; i < numGroups; ++i)
            {
                for (int j = 0; j < howmany; ++j)
                {
                    for (int k = stride; k < n2 * stride; k += stride)
                    {
                        if (realInput[k] != realInput[stride * m_fftSize - k])
                        {
                            conjugateSymmetric = false;
                            break;
                        }
                    }

                    if (!conjugateSymmetric)
                        break;

                    realInput += dist;
                }

                if (!conjugateSymmetric)
                    break;

                realInput += stride * m_fftSize - howmany * dist;
            }
        }
        else
        {
            const hwComplex* complexInput = input.GetComplexData();

            for (int i = 0; i < numGroups; ++i)
            {
                for (int j = 0; j < howmany; ++j)
                {
                    if (complexInput[0].Imag() != 0.0)
                    {
                        conjugateSymmetric = false;
                    }
                    else
                    {
                        for (int k = stride; k < n2 * stride; k += stride)
                        {
                            if (complexInput[k] != complexInput[stride * m_fftSize - k].Conjugate())
                            {
                                conjugateSymmetric = false;
                                break;
                            }
                        }

                        complexInput += dist;
                    }

                    if (!conjugateSymmetric)
                        break;
                }

                if (!conjugateSymmetric)
                    break;

                complexInput += stride * m_fftSize - howmany * dist;
            }
        }
    }

    if (conjugateSymmetric)
    {
        m_type = FFT_C2R;
    }
    else
    {
        m_type = FFT_C2C;
    }
}
