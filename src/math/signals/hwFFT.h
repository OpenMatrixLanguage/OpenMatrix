/**
* @file  hwFFT.h
* @date July 2011
* Copyright (C) 2011-2018 Altair Engineering, Inc.  
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
#ifndef _Signals_FFT_h
#define _Signals_FFT_h

#include "fftw3.h"
#include "hwMathStatus.h"

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

//!
//! \enum FFT_TYPE
//!
enum FFT_TYPE
{ 
    FFT_NOTSET, 
    FFT_R2C, 
    FFT_C2R, 
    FFT_C2C
};

//------------------------------------------------------------------------------
//!
//! \class hwFFTW
//! \brief Fast fourier transform classes
//!
//------------------------------------------------------------------------------
class hwFFTW
{
public:
    //!
    //! Returns size
    //!
    int Size() const { return m_fftSize; }
    //!
    //! Returns status
    //!
    const hwMathStatus& Status() const { return m_status;}

    //!
    //! Call FFTW for const input and return status
    //! \param input
    //! \param output
    //!
    hwMathStatus Compute(const hwMatrix& input, 
                         hwMatrix&       output);
    //!
    //! Call FFTW for non-const input and return status
    //! \param input
    //! \param output
    //!
    hwMathStatus Compute(hwMatrix& input, 
                         hwMatrix& output);

protected:
    hwFFTW(const unsigned& direction, int fftSize = 0);
    virtual ~hwFFTW();

    FFT_TYPE m_type;
    int m_dataSize;
    int m_fftSize;
    unsigned m_direction;
    double* m_realInput;
    double* m_realOutput;
    fftw_complex* m_complexInput;
    fftw_complex* m_complexOutput;
    fftw_plan m_plan;
    hwMathStatus m_status;

    virtual void DetermineType(const hwMatrix& input) = 0;

private:
    //!
    //! Allocate memory based on FFT type
    //! \param input
    //!
    void SetupFFT(const hwMatrix& input);
    //!
    //! Call FFTW, selecting the appropriate plan
    //! \param input
    //! \param output
    //! \param options
    //!

    hwMathStatus Compute(hwMatrix&       input, 
                         hwMatrix&       output, 
                         const unsigned& options);
};
//------------------------------------------------------------------------------
//!
//! \class hwFFT_f
//! \brief Fast fourier transform classes
//!
//------------------------------------------------------------------------------
class hwFFT_f : public hwFFTW
{
public:
    //!
    //! Constructor
    //! \param fftSize
    //!
    hwFFT_f(int fftSize = 0);
    //!
    //! Destructor
    //!
    virtual ~hwFFT_f();

protected:
    //!
    //! Determines the type
    //! \param input
    //!
    void DetermineType(const hwMatrix& input);
};
//------------------------------------------------------------------------------
//!
//! \class hwFFT_r
//! \brief Fast fourier transform classes
//!
//------------------------------------------------------------------------------
class hwFFT_r : public hwFFTW
{
public:
    //!
    //! Constructor
    //! \param fftSize
    //! \param assumeSymmetry 
    //!
    hwFFT_r(int  fftSize = 0, 
            bool assumeSymmetry = true);
    //!
    //! Destructor
    //!
    ~hwFFT_r();

protected:
    //!
    //! Determines the type
    //! \param input
    //!
    void DetermineType(const hwMatrix& input);

private:
    bool m_assumeSymmetry; //!< 

};

#endif // _Signals_FFT_h
