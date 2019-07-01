/**
* @file  hwFilterManager.cxx
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
#include "hwMatrix.h"
#include "hwFilterManager.h"

#ifndef OS_WIN
#define _stricmp strcmp
#endif

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwFilterManager::hwFilterManager()
    : m_pFilter(NULL)
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwFilterManager::~hwFilterManager()
{
    if (m_pFilter)
    {
        delete m_pFilter;
    }
}
//------------------------------------------------------------------------------
// Returns status after creating filter with transfer function polynomials
//------------------------------------------------------------------------------
hwMathStatus hwFilterManager::CreateFilter(const hwMatrix& numerCoef, 
                                           const hwMatrix* denomCoef)
{
    if (m_pFilter)
    {
        delete m_pFilter;
        m_pFilter = nullptr;
    }

    m_pFilter = new hwDigitalFilter(numerCoef, denomCoef);
    return m_pFilter->Status();
}
//------------------------------------------------------------------------------
// Returns status after applying filter to signal data
//------------------------------------------------------------------------------
hwMathStatus hwFilterManager::ApplyFilter(const hwMatrix& inSignal, 
                                          hwMatrix&       outSignal,
                                          const hwMatrix* initCond)
{
    if (!m_pFilter)
    {
        return hwMathStatus(HW_MATH_ERR_NULLPOINTER);
    }
    if (!inSignal.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!inSignal.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (initCond)
    {
        if (!initCond->IsReal())
        {
            return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);
        }
        if (!initCond->IsVector())
        {
            return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 3);
        }
        if (initCond->Size() != m_pFilter->GetDenomCoefs()->Size() - 1)
        {
            return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 3);
        }
    }

    hwMathStatus status = outSignal.Dimension(inSignal.M(), inSignal.N(), 
                                              hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(2);
        }
        return status;
    }

    int max;
    int size = inSignal.Size();

    const double* pInSignal  = inSignal.GetRealData();
    double*       pOutSignal = outSignal.GetRealData();

    int     numNumerCoefs = m_pFilter->GetNumerCoefs()->Size();
    int     numDenomCoefs = 0;
    double* numerCoef     = m_pFilter->GetNumerCoefs()->GetRealData();
    double* denomCoef     = nullptr;
    if (m_pFilter->GetDenomCoefs())
    {
        denomCoef     = m_pFilter->GetDenomCoefs()->GetRealData();
        numDenomCoefs = m_pFilter->GetDenomCoefs()->Size();
    }

    if (denomCoef)
    {
        if (IsZero(denomCoef[0], 1.0e-12))
        {
            return status(HW_MATH_ERR_FILTERDENZERO);
        }
    }

    for (int i = 0; i < size; ++i)
    {
        max = _min(i + 1, numNumerCoefs);
        outSignal(i) = 0.0;

        for (int j = 0; j < max; ++j)
        {
            pOutSignal[i] += numerCoef[j] * pInSignal[i-j];
        }

        if (denomCoef)     // used with IIR filters, but not FIR
        {
            max = _min(i + 1, numDenomCoefs);

            for (int j = 1; j < max; ++j)
            {
                pOutSignal[i] -= denomCoef[j] * pOutSignal[i-j];
            }
            if (denomCoef[0] != 1.0)
            {
                pOutSignal[i] /= denomCoef[0];
            }
            // initCond values are for the state vector, not for the output signal.
            if (initCond && i < numDenomCoefs - 1)
            {
                pOutSignal[i] += (*initCond)(i);
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Returns status after applying filter to signal data in reverse
//------------------------------------------------------------------------------
hwMathStatus hwFilterManager::ApplyReverseFilter(const hwMatrix& inSignal,
                                                 hwMatrix&       outSignal,
                                                 const hwMatrix* initCond)
{
    if (!m_pFilter)
    {
        return hwMathStatus(HW_MATH_ERR_NULLPOINTER);
    }
    if (!inSignal.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!inSignal.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    if (initCond)
    {
        if (!initCond->IsReal())
        {
            return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);
        }
        if (!initCond->IsVector())
        {
            return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 3);
        }
        if (initCond->Size() != m_pFilter->GetDenomCoefs()->Size() - 1)
        {
            return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 3);
        }
    }

    hwMathStatus status = outSignal.Dimension(inSignal.M(), inSignal.N(), 
                                              hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(2);
        }
        return status;
    }

    int max;
    int size = inSignal.Size();
    int numNumerCoefs = m_pFilter->GetNumerCoefs()->Size();
    int numDenomCoefs = 0;
    double* numerCoef = m_pFilter->GetNumerCoefs()->GetRealData();
    double* denomCoef = nullptr;
    const double* pInSignal  = inSignal.GetRealData();
    double*       pOutSignal = outSignal.GetRealData();

    if (m_pFilter->GetDenomCoefs())
    {
        denomCoef = m_pFilter->GetDenomCoefs()->GetRealData();
        numDenomCoefs = m_pFilter->GetDenomCoefs()->Size();
    }
    if (denomCoef)
    {
        if (IsZero(denomCoef[0], 1.0e-12))
        {
            return status(HW_MATH_ERR_FILTERDENZERO);
        }
    }

    int index = size-1;

    for (int i = 0; i < size; ++i)
    {
        max = _min(i + 1, numNumerCoefs);

        outSignal(index) = 0.0;

        for (int j = 0; j < max; ++j)
        {
            pOutSignal[index] += numerCoef[j] * pInSignal[index+j];
        }

        if (denomCoef)     // used with IIR filters, but not FIR
        {
            max = _min(i + 1, numDenomCoefs);

            for (int j = 1; j < max; ++j)
            {
                pOutSignal[index] -= denomCoef[j] * pOutSignal[index+j];
            }

            if (denomCoef[0] != 1.0)
            {
                pOutSignal[index] /= denomCoef[0];
            }
            // initCond values are for the state vector, not for the output signal.
            if (initCond && i < numDenomCoefs - 1)
            {
                pOutSignal[index] += (*initCond)(i);
            }
        }

        index--;
    }

    return status;
}
//------------------------------------------------------------------------------
// Returns status after applying zero phase filter to signal data
//------------------------------------------------------------------------------
hwMathStatus hwFilterManager::ApplyZeroPhaseFilter(const hwMatrix& inSignal, 
                                                   hwMatrix&       outSignal)
{
    if (!m_pFilter)
    {
        return hwMathStatus(HW_MATH_ERR_NULLPOINTER);
    }
    if (!inSignal.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!inSignal.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }
    hwMathStatus status = outSignal.Dimension(inSignal.M(), inSignal.N(), 
                                              hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(2);
        }
        return status;
    }

    // zero pad if needed
    double*          numerCoef  = m_pFilter->GetNumerCoefs()->GetRealData();
    double*          denomCoef  = nullptr;
    hwMatrix*        temp       = nullptr;
    hwDigitalFilter* tempFilter = nullptr;

    int numNumerCoefs = m_pFilter->GetNumerCoefs()->Size();
    int numDenomCoefs = 0;

    if (m_pFilter->GetDenomCoefs())
    {
        denomCoef     = m_pFilter->GetDenomCoefs()->GetRealData();
        numDenomCoefs = m_pFilter->GetDenomCoefs()->Size();
    }

    int numCoef = _max(numNumerCoefs, numDenomCoefs);
    int n       = inSignal.Size();
    int nfact   = 3 * (numCoef - 1);

    if (numCoef < 2)
    {
        return hwMathStatus(HW_MATH_ERR_FILTFILTCOEFS);
    }

    if (n <= nfact)
    {
        return hwMathStatus(HW_MATH_ERR_FILTFILTDATA, 1);
    }
    if (numNumerCoefs < numCoef)
    {
        temp = new hwMatrix(numCoef, hwMatrix::REAL);

        for (int i = 0; i < numNumerCoefs; ++i)
        {
            (*temp)(i) = numerCoef[i];
        }

        for (int i = numNumerCoefs; i < numCoef; ++i)
        {
            (*temp)(i) = 0.0;
        }

        tempFilter = new hwDigitalFilter(*temp, m_pFilter->GetDenomCoefs());
        delete m_pFilter;
        m_pFilter = tempFilter;
        numerCoef = m_pFilter->GetNumerCoefs()->GetRealData();
    }

    if (numDenomCoefs < numCoef)
    {
        temp = new hwMatrix(numCoef, hwMatrix::REAL);

        if (numDenomCoefs == 0)     // FIR condition
        {
            (*temp)(0) = 1.0;
            ++numDenomCoefs;
        }
        else
            (*temp)(0) = denomCoef[0];

        for (int i = 1; i < numDenomCoefs; ++i)
        {
            (*temp)(i) = denomCoef[i];
        }
        for (int i = numDenomCoefs; i < numCoef; ++i)
        {
            (*temp)(i) = 0.0;
        }
        tempFilter = new hwDigitalFilter(*m_pFilter->GetNumerCoefs(), temp);
        delete m_pFilter;
        m_pFilter = tempFilter;
        numerCoef = m_pFilter->GetNumerCoefs()->GetRealData();
        denomCoef = m_pFilter->GetDenomCoefs()->GetRealData();
    }

    if (temp)
    {
        delete temp;
        temp = nullptr;
    }

    // set working matrix A
    hwMatrix A(numCoef-1, numCoef-1, hwMatrix::REAL);

    // set A: first row
    A(0, 0) = 1.0 + denomCoef[1];

    if (numCoef > 2)
    {
        A(0, 1) = -1.0;
    }

    for (int j = 2; j < numCoef-1; ++j)
    {
        A(0, j) = 0.0;
    }

    // set A: next numCoef-3 rows
    for (int i = 1; i < numCoef-2; ++i)
    {
        A(i, 0) = denomCoef[i+1];

        for (int j = 1; j < i; ++j)
        {
            A(i, j) = 0.0;
        }
        A(i, i) = 1.0;

        if (i + 1 < numCoef-1)
        {
            A(i, i+1) = -1.0;
        }
        for (int j = i + 2; j < numCoef-1; ++j)
        {
            A(i, j) = 0.0;
        }
    }

    // set A: last row
    if (numCoef > 2)
    {
        A(numCoef-2, 0) = denomCoef[numCoef-1];

        for (int j = 1; j < numCoef-2; ++j)
        {
            A(numCoef-2, j) = 0.0;
        }
        A(numCoef-2, numCoef-2) = 1.0;
    }

    // extrapolate beginning and end of data using "reflection method"
    int           n_out      = n + 2*nfact;
    const double* pInSignal  = inSignal.GetRealData();
    double*       pOutSignal = outSignal.GetRealData();

    hwMatrix tempSignal_1(n_out, hwMatrix::REAL);
    double*  pTempSignal_1 = tempSignal_1.GetRealData();

    for (int i = 0; i < nfact; ++i)
    {
        pTempSignal_1[i] = 2.0 * pInSignal[0] - pInSignal[nfact - i];
    }
    memcpy(pTempSignal_1 + nfact, pInSignal, n * sizeof(double));

    for (int i = n + nfact; i < n_out; ++i)
    {
        pTempSignal_1[i] = 2.0 * pInSignal[n - 1] - pInSignal[2*n + nfact - i - 2];
    }

    // compute filter state vector steady state (DC) values
    hwMatrix a2(numCoef-1, 1, hwMatrix::REAL);
    hwMatrix b2(numCoef-1, 1, hwMatrix::REAL);
    hwMatrix ssCond(numCoef-1, 1, hwMatrix::REAL);    // steady state condition

    for (int i = 0; i < numCoef-1; ++i)
    {
        a2(i) = denomCoef[i+1];
        b2(i) = numerCoef[i+1];
    }

    a2 *= numerCoef[0];

    status = ssCond.LSolve(A, b2 - a2);

    if (!status.IsOk())
    {
        status.ResetArgs();

        if (!status.IsWarning())
        {
            return status;
        }
        status = hwMathStatus();    // suppress warning for now
    }

    // filter the signal for zero phase
    hwMatrix initCond(numCoef-1, hwMatrix::REAL);     // state vector initial conditions
    hwMatrix tempSignal_2(n_out, hwMatrix::REAL);
    double* pTempSignal_2 = tempSignal_2.GetRealData();

    initCond = ssCond * pTempSignal_1[0];
    status = ApplyFilter(tempSignal_1, tempSignal_2, &initCond);     // filter

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    initCond = ssCond * pTempSignal_2[n_out-1];
    ApplyReverseFilter(tempSignal_2, tempSignal_1, &initCond);      // reverse filter

    // remove extrapolated segments from signal
    for (int i = 0; i < n; ++i)
    {
        pOutSignal[i] = pTempSignal_1[nfact + i];
    }
    return status;
}
//------------------------------------------------------------------------------
// Returns status and ups sample by a given factor
//------------------------------------------------------------------------------
hwMathStatus hwFilterManager::UpSample(const hwMatrix& inSignal, 
                                       hwMatrix&       outSignal,
                                       int             k, 
                                       int             phase)
{
    if (!inSignal.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!inSignal.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    // compute upsample size
    if (k < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 3);
    }
    if (phase < 0 || phase > k - 1)
    {
        return hwMathStatus(HW_MATH_ERR_RESAMPOFFSET, 4);
    }
    int m = inSignal.Size() * k;  // Length of output signal
    hwMathStatus status = outSignal.Dimension(m, hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(2);
        }
        return status;
    }

    // up sample
    const double* pInSignal  = inSignal.GetRealData();
    double*       pOutSignal = outSignal.GetRealData();
    
    int count = 0;
    for (int i = 0; i < m; ++i)
    {
        if ((i - phase)%k == 0)
        {
            pOutSignal[i] = pInSignal[count++];
        }
        else
        {
            pOutSignal[i] = 0;
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Returns status and downs sample by a given factor
//------------------------------------------------------------------------------
hwMathStatus hwFilterManager::DownSample(const hwMatrix& inSignal, 
                                         hwMatrix&       outSignal,
                                         int             k, 
                                         int             phase)
{
    if (!inSignal.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!inSignal.IsEmptyOrVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    // compute downsample size
    if (k < 1)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSINT, 3);
    }
    if (phase < 0 || phase > k - 1)
    {
        return hwMathStatus(HW_MATH_ERR_RESAMPOFFSET, 4);
    }
    int m = (inSignal.Size() - phase - 1) / k + 1; // Length of outsignal
    hwMathStatus status = outSignal.Dimension(m, hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(2);
        }
        return status;
    }

    // down sample
    const double* pInSignal  = inSignal.GetRealData();
    double*       pOutSignal = outSignal.GetRealData();

    int count = 0;
    for (int i = 0; i < m; ++i)
    {
        pOutSignal[i] = pInSignal[count + phase];
        count += k;
    }

    return hwMathStatus();
}
