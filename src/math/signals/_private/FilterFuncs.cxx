/**
* @file * @file FilterFuncs.cxx
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

#include "FilterFuncs.h"

#include "hwMatrix.h"
#include "PolynomFuncs.h"
#include "MathUtilsFuncs.h"
#include "StatisticsFuncs.h"
#include "hwFilterManager.h"

// digital
#include "hwBessel_z.h"
#include "hwButterworth_z.h"
#include "hwChebyshev_I_z.h"
#include "hwChebyshev_II_z.h"
#include "hwElliptic_z.h"
#include "hwDigitalFilter_FIR.h"
// analog
#include "hwBessel_s.h"
#include "hwButterworth_s.h"
#include "hwChebyshev_I_s.h"
#include "hwChebyshev_II_s.h"
#include "hwElliptic_s.h"
// design
#include "hwBessel_Design.h"
#include "hwButterworth_Design.h"
#include "hwChebyshev_I_Design.h"
#include "hwChebyshev_II_Design.h"
#include "hwElliptic_Design.h"
// windows
#include "hwChebyshev.h"
#include "hwHamming.h"
#include "hwHanning.h"
#include "hwBartlettHann.h"
#include "hwBlackman.h"
#include "hwWelch.h"
#include "hwParzen.h"
#include "hwKaiserBessel.h"

//------------------------------------------------------------------------------
// Computes filter transfer function coefficients and returns status
//------------------------------------------------------------------------------
hwMathStatus Besself(int         order, 
                     double      lowCutoffFreq, 
                     double      highCutoffFreq,
                     hwMatrix&   numerCoef,
                     hwMatrix&   denomCoef, 
                     const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 6);
    }

    hwFilter* pFilter = nullptr;
    if (!strcmp(type, "z"))
    {
        pFilter = new hwBessel_z(order, lowCutoffFreq, highCutoffFreq, "besself");
    }
    else if (!strcmp(type, "s"))
    {
        pFilter = new hwBessel_s(order, lowCutoffFreq, highCutoffFreq, "besself");
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_FILTERTYPE, 6);
    }

    hwMathStatus status = pFilter->Status();
    if (!status.IsOk())
	{
        if (status.GetArg1() == 4)
        {
            status.SetArg1(6);
        }

        if (!status.IsWarning())
        {
            delete pFilter;
            return status;
        }
    }

    numerCoef = *(pFilter->GetNumerCoefs());
    denomCoef = *(pFilter->GetDenomCoefs());

    delete pFilter;
    return status;
}
//------------------------------------------------------------------------------
// Computes filter transfer function coefficients and returns status
//------------------------------------------------------------------------------
hwMathStatus Besself3(int         order, 
                      double      lowCutoffFreq, 
                      double      highCutoffFreq,
                      hwMatrix&   numerCoef, 
                      hwMatrix&   denomCoef, 
                      const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 6);
    }

    hwFilter* pFilter = nullptr; 
    if (!strcmp(type, "z"))
    {
        pFilter = new hwBessel_z(order, lowCutoffFreq, highCutoffFreq, "besself3");
    }
    else if (!strcmp(type, "s"))
    {
        pFilter = new hwBessel_s(order, lowCutoffFreq, highCutoffFreq, "besself3");
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_FILTERTYPE, 6);
    }

    hwMathStatus status = pFilter->Status();
    if (!status.IsOk())
	{
        if (status.GetArg1() == 4)
        {
            status.SetArg1(6);
        }
        if (!status.IsWarning())
        {
            delete pFilter;
            return status;
        }
    }

    numerCoef = *(pFilter->GetNumerCoefs());
    denomCoef = *(pFilter->GetDenomCoefs());

    delete pFilter;
    return status;
}
//------------------------------------------------------------------------------
// Computes Butterworth filter transfer function coefficients and returns status
//------------------------------------------------------------------------------
hwMathStatus Butter(int         order, 
                    double      lowCutoffFreq, 
                    double      highCutoffFreq,
                    hwMatrix&   numerCoef, 
                    hwMatrix&   denomCoef, 
                    const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 6);
    }

    hwFilter* pFilter = nullptr;
    if (!strcmp(type, "z"))
    {
        pFilter = new hwButterworth_z(order, lowCutoffFreq, highCutoffFreq);
    }
    else if (!strcmp(type, "s"))
    {
        pFilter = new hwButterworth_s(order, lowCutoffFreq, highCutoffFreq);
    }
    else
    {
        return hwMathStatus (HW_MATH_ERR_FILTERTYPE, 6);
    }

    hwMathStatus status = pFilter->Status();
    if (!status.IsOk() && !status.IsWarning())
    {
        delete pFilter;
        return status;
    }

    numerCoef = *(pFilter->GetNumerCoefs());
    denomCoef = *(pFilter->GetDenomCoefs());

    delete pFilter;

    return status;
}
//------------------------------------------------------------------------------
// Computes Chebyshev type I filter transfer function coefficients
//------------------------------------------------------------------------------
hwMathStatus Cheby1(int         order, 
                    double      lowCutoffFreq, 
                    double      highCutoffFreq,
                    double      passEdgeDb, 
                    hwMatrix&   numerCoef, 
                    hwMatrix&   denomCoef,
                    const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 7);
    }

    hwFilter* pFilter = nullptr;
    if (!strcmp(type, "z"))
    {
        pFilter = new hwChebyshev_I_z(order, lowCutoffFreq, highCutoffFreq, passEdgeDb);
    }
    else if (!strcmp(type, "s"))
    {
        pFilter = new hwChebyshev_I_s(order, lowCutoffFreq, highCutoffFreq, passEdgeDb);
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_FILTERTYPE, 7);
    }

    hwMathStatus status = pFilter->Status();
    if (!status.IsOk() && !status.IsWarning())
    {
        delete pFilter;
        return status;
    }

    numerCoef = *(pFilter->GetNumerCoefs());
    denomCoef = *(pFilter->GetDenomCoefs());

    delete pFilter;

    return status;
}
//------------------------------------------------------------------------------
// Computes Chebyshev type II filter transfer function coefficients
//------------------------------------------------------------------------------
hwMathStatus Cheby2(int         order, 
                    double      lowCutoffFreq, 
                    double      highCutoffFreq,
                    double      stopEdgeDb, 
                    hwMatrix&   numerCoef, 
                    hwMatrix&   denomCoef,
                    const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 7);
    }

    hwFilter* pFilter = nullptr;
    if (!strcmp(type, "z"))
    {
        pFilter = new hwChebyshev_II_z(order, lowCutoffFreq, highCutoffFreq, stopEdgeDb);
    }
    else if (!strcmp(type, "s"))
    {
        pFilter = new hwChebyshev_II_s(order, lowCutoffFreq, highCutoffFreq, stopEdgeDb);
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_FILTERTYPE, 7);
    }

    hwMathStatus status = pFilter->Status();
    if (!status.IsOk() && !status.IsWarning())
    {
        delete pFilter;
        return status;
    }

    numerCoef = *(pFilter->GetNumerCoefs());
    denomCoef = *(pFilter->GetDenomCoefs());

    delete pFilter;

    return status;
}
//------------------------------------------------------------------------------
// Computes Elliptic filter transfer function coefficients and returns status
//------------------------------------------------------------------------------
hwMathStatus Ellip(int         order,
                   double      lowCutoffFreq,
                   double      highCutoffFreq,
                   double      passEdgeDb,
                   double      stopEdgeDb,
                   hwMatrix&   numerCoef,
                   hwMatrix&   denomCoef,
                   const char* type)
{
    if (!type)
    {
        return hwMathStatus (HW_MATH_ERR_NOTSTRING, 8);
    }

    hwFilter* pFilter = nullptr;
    if (!strcmp(type, "z"))
    {
        pFilter = new hwElliptic_z(order, lowCutoffFreq, highCutoffFreq, 
                                   passEdgeDb, stopEdgeDb);
    }
    else if (!strcmp(type, "s"))
    {
        pFilter = new hwElliptic_s(order, lowCutoffFreq, highCutoffFreq, 
                                   passEdgeDb, stopEdgeDb);
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_FILTERTYPE, 8);
    }

    hwMathStatus status = pFilter->Status();
    if (!status.IsOk() && !status.IsWarning())
    {
        delete pFilter;
        return status;
    }

    numerCoef = *(pFilter->GetNumerCoefs());
    denomCoef = *(pFilter->GetDenomCoefs());

    delete pFilter;

    return status;
}
//------------------------------------------------------------------------------
// Computes FIR filter transfer function coefficients and returns status
//------------------------------------------------------------------------------
hwMathStatus Fir(int             order,
                 double          lowCutoffFreq,
                 double          highCutoffFreq,
                 const hwMatrix* window,
                 hwMatrix&       numerCoef,
                 bool            normalize)
{
    hwDigitalFilter_FIR filter(order, lowCutoffFreq, highCutoffFreq, window, normalize);
    hwMathStatus status = filter.Status();

    if (!status.IsOk())
        return status;

    numerCoef = (*filter.GetNumerCoefs());

    return status;
}
//------------------------------------------------------------------------------
// Computes FIR filter transfer function coefficients for multiband filters
//------------------------------------------------------------------------------
hwMathStatus Fir(int             order, 
                 const hwMatrix& cutoffFreq, 
                 const hwMatrix* window, 
                 hwMatrix&       numerCoef,
                 const char*     type, 
                 bool            normalize)
{
    hwDigitalFilter_FIR filter(order, cutoffFreq, window, type, normalize);

    hwMathStatus status = filter.Status();
    if (!status.IsOk())
    {
        return status;
    }

    numerCoef = (*filter.GetNumerCoefs());

    return status;
}
//------------------------------------------------------------------------------
// Computes FIR filter transfer function coefficients for multiband filters
//------------------------------------------------------------------------------
hwMathStatus FirLS(int             order,
                   const hwMatrix& freq,
                   const hwMatrix& mag,
                   const hwMatrix* weight,
                   hwMatrix&       filterCoef)
{
    // check inputs
    if (order < 1)
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 1);

    int n = order + order % 2;     // force even order
    int m = n / 2;

    if (!freq.IsVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);

    if (!freq.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);

    if (!mag.IsVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 3);

    if (!mag.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);

    int numPts = freq.Size();
    int numBands = numPts / 2;

    if (numPts%2 != 0)
        return hwMathStatus(HW_MATH_ERR_VECTOR2x, 2);

    if (mag.Size() != numPts)
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 2, 3);

    hwMatrix wvec;
    hwMathStatus status;

    if (weight)
    {
        if (!weight->IsVector())
            return hwMathStatus(HW_MATH_ERR_VECTOR, 4);

        if (!weight->IsReal())
            return hwMathStatus(HW_MATH_ERR_COMPLEX, 4);

        if (2 * weight->Size() != numPts)
            return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 2, 4);       // need new message

        wvec = (*weight);

        status = wvec.Reshape(numBands, 1);
    }
    else
    {
        status = wvec.Dimension(numBands, 1, hwMatrix::REAL);
        wvec.SetElements(1.0);
    }

    hwMatrix wk;    // weights 
    hwMatrix pm(2, 1, hwMatrix::REAL);  // plus / minus

    pm(0) = -1.0;
    pm(1) = 1.0;

    wk.Kronecker(wvec, pm);

    // set up intermediate variables
    hwMatrix omegaR(1, numPts, hwMatrix::REAL);   // row of angular frequency band limits
    hwMatrix omegaL(1, numBands, hwMatrix::REAL); // row of lower angular frequency band limits
    hwMatrix omegaU(1, numBands, hwMatrix::REAL); // row of upper angular frequency band limits
    hwMatrix ampL(1, numBands, hwMatrix::REAL);   // row of lower angular frequency band amplitudes
    hwMatrix ampU(1, numBands, hwMatrix::REAL);   // row of upper angular frequency band amplitudes

    for (int i = 0; i < numBands; ++i)
    {
        int ii = 2 * i;
        omegaR(ii) = omegaL(i) = freq(ii) * PI;
        omegaR(ii + 1) = omegaU(i) = freq(ii + 1) * PI;
        ampL(i) = mag(ii);
        ampU(i) = mag(ii + 1);
    }

    // construct Toeplitz + Hankel LHS matrix
    hwMatrix idxN(n, 1, hwMatrix::REAL);
    hwMatrix idxN2(n + 1, 1, hwMatrix::REAL);

    idxN2(0) = 1.0;

    for (int i = 1; i < n + 1; ++i)
    {
        idxN(i - 1) = static_cast<double> (i);
        idxN2(i)    = 1.0 / static_cast<double> (i);
    }

    hwMatrix temp1(idxN * omegaR);
    hwMatrix lhs_ints(n + 1, numPts, hwMatrix::REAL);
    status = lhs_ints.WriteRow(0, omegaR);

    for (int j = 0; j < numPts; ++j)
    {
        for (int i = 1; i < n + 1; ++i)
        {
            lhs_ints(i, j) = sin(temp1(i - 1, j));
        }
    }
 
    hwMatrix DiagElems;
    status = DiagElems.MultByElems(idxN2, lhs_ints * wk);
    hwMatrix ToeplitzElems(1, m + 1, (void*) DiagElems.GetRealData(), hwMatrix::REAL);
    hwMatrix HankelElems(1, m + 1, (void*) (DiagElems.GetRealData() + m), hwMatrix::REAL);
    hwMatrix T;
    hwMatrix H;

    T.Toeplitz(ToeplitzElems);
    H.Hankel(ToeplitzElems, &HankelElems);
    hwMatrix C(T + H);

    // construct RHS matrix
    hwMatrix omegaL2;
    hwMatrix omegaU2;
    hwMatrix idxM(m, 1, hwMatrix::REAL);
    hwMatrix idxM2(m + 1, 1, hwMatrix::REAL);

    status = omegaL2.MultByElems(omegaL, omegaL);
    status = omegaU2.MultByElems(omegaU, omegaU);

    idxM2(0) = 2.0;

    for (int i = 1; i < m + 1; ++i)
    {
        idxM(i - 1) = static_cast<double> (i);
        idxM2(i)    = static_cast<double> (i);
    }

    hwMatrix temp2(idxM * omegaU);
    hwMatrix temp3(idxM * omegaL);
    hwMatrix temp4(m + 1, numBands, hwMatrix::REAL);

    status = temp4.WriteRow(0, omegaL2 - omegaU2);

    for (int j = 0; j < numBands; ++j)
    {
        for (int i = 1; i < m + 1; ++i)
        {
            temp4(i, j) = cos(temp2(i - 1, j)) - cos(temp3(i - 1, j));
        }
    }

    hwMatrix temp5(idxM2 * (omegaU - omegaL));
    hwMatrix rhs_ints;
    
    status = rhs_ints.DivideByElems(temp4, temp5);

    // construct d
    hwMatrix d(2, numBands, hwMatrix::REAL);

    status = wvec.Reshape(1, numBands);
    status = temp2.MultByElems(-wvec, ampL);
    status = temp3.MultByElems( wvec, ampU);
    status = d.WriteSubmatrix(0, 0, temp2);
    status = d.WriteSubmatrix(1, 0, temp3);
    status = d.Reshape(d.Size(), 1);

    // construct b
    idxM2(0) = 1.0;

    for (int i = 1; i < m + 1; ++i)
    {
        idxM2(i) = 1.0 / static_cast<double> (i);
    }

    hwMatrix temp6;
    hwMatrix temp7;
    hwMatrix b;
    hwMatrix a;
    hwMatrix pp(1, 2, hwMatrix::REAL);
    pp.SetElements(1.0);
    status = temp6.Kronecker(rhs_ints, pp);
    status = temp7.ReadSubmatrix(0, 0, m + 1, numPts, lhs_ints);
    status = b.MultByElems(idxM2, (temp6 + temp7) * d);
    status = a.QRSolve(C, b);

    filterCoef.Dimension(n + 1, 1, hwMatrix::REAL);

    for (int i = 0; i < m; ++i)
    {
        filterCoef(i) = a(m - i);
        filterCoef(n - i) = filterCoef(i);
    }

    filterCoef(m) = 2.0 * a(0);

    return status;
}
//------------------------------------------------------------------------------
// Designs a Bessel filter and returns the status
//------------------------------------------------------------------------------
hwMathStatus BesselOrd(double      passBandFreq, 
                       double      stopBandFreq,
                       double      passEdgeDb, 
                       double      stopEdgeDb, 
                       int&        order,
                       hwMatrix&   freqC, 
                       const char* type)
{
    hwBessel_Design bessel;
    hwMathStatus status;

    if (!type)
    {
        status(HW_MATH_ERR_NOTSTRING, 7);
    }
    else if (!strcmp(type, "z"))
    {
        status = bessel.Digital(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                                order, freqC);
    }
    else if (!strcmp(type, "s"))
    {
        status = bessel.Analog(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                               order, freqC);
    }
    else
    {
        status(HW_MATH_ERR_FILTERTYPE, 7);
    }

    return status;
}
//------------------------------------------------------------------------------
// Designs a Bessel filter and returns the status
//------------------------------------------------------------------------------
hwMathStatus BesselOrd(const hwMatrix& passBandFreq, 
                       const hwMatrix& stopBandFreq,
                       double          passEdgeDb, 
                       double          stopEdgeDb, 
                       int&            order,
                       hwMatrix&       freqC1, 
                       hwMatrix&       freqC2, 
                       const char*     type)
{
    hwBessel_Design bessel;
    hwMathStatus status;

    if (!type)
    {
        status(HW_MATH_ERR_NOTSTRING, 8);
    }
    else if (!strcmp(type, "z"))
    {
        status = bessel.Digital(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                                order, freqC1, freqC2);
    }
    else if (!strcmp(type, "s"))
    {
        status = bessel.Analog(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                               order, freqC1, freqC2);
    }
    else
    {
        status(HW_MATH_ERR_FILTERTYPE, 8);
    }
    return status;
}
//------------------------------------------------------------------------------
// Designs a Butterworth filter and returns the status
//------------------------------------------------------------------------------
hwMathStatus ButterOrd(double      passBandFreq, 
                       double      stopBandFreq,
                       double      passEdgeDb, 
                       double      stopEdgeDb, 
                       int&        order, 
                       hwMatrix&   freqC, 
                       const char* type)
{
    hwButterworth_Design butter;
    hwMathStatus status;

    if (!type)
    {
        status(HW_MATH_ERR_NOTSTRING, 7);
    }
    else if (!strcmp(type, "z"))
    {
        status = butter.Digital(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                                order, freqC);
    }
    else if (!strcmp(type, "s"))
    {
        status = butter.Analog(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                               order, freqC);
    }
    else
    {
        status(HW_MATH_ERR_FILTERTYPE, 7);
    }

    return status;
}
//------------------------------------------------------------------------------
// Designs a Butterworth filter and returns the status
//------------------------------------------------------------------------------
hwMathStatus ButterOrd(const hwMatrix& passBandFreq, 
                       const hwMatrix& stopBandFreq,
                       double          passEdgeDb, 
                       double          stopEdgeDb,
                       int&            order,
                       hwMatrix&       freqC1, 
                       hwMatrix&       freqC2, 
                       const char*     type)
{
    hwButterworth_Design butter;
    hwMathStatus status;

    if (!type)
    {
        status(HW_MATH_ERR_NOTSTRING, 8);
    }
    else if (!strcmp(type, "z"))
    {
        status = butter.Digital(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                                order, freqC1, freqC2);
    }
    else if (!strcmp(type, "s"))
    {
        status = butter.Analog(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                               order, freqC1, freqC2);
    }
    else
    {
        status(HW_MATH_ERR_FILTERTYPE, 8);
    }

    return status;
}
//------------------------------------------------------------------------------
// Designs a Chebyshev type I filter and returns the status
//------------------------------------------------------------------------------
hwMathStatus Cheby1Ord(double      passBandFreq, 
                       double      stopBandFreq,
                       double      passEdgeDb,
                       double      stopEdgeDb,
                       int&        order,
                       hwMatrix&   freqC, 
                       const char* type)
{
    hwChebyshev_I_Design chebyI;
    hwMathStatus         status;

    if (!type)
    {
        status(HW_MATH_ERR_NOTSTRING, 7);
    }
    else if (!strcmp(type, "z"))
    {
        status = chebyI.Digital(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                                order, freqC);
    }
    else if (!strcmp(type, "s"))
    {
        status = chebyI.Analog(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                               order, freqC);
    }
    else
    {
        status(HW_MATH_ERR_FILTERTYPE, 7);
    }
    return status;
}
//------------------------------------------------------------------------------
// Designs a Chebyshev type I filter and returns the status
//------------------------------------------------------------------------------
hwMathStatus Cheby1Ord(const hwMatrix& passBandFreq, 
                       const hwMatrix& stopBandFreq,
                       double          passEdgeDb, 
                       double          stopEdgeDb, 
                       int&            order,
                       hwMatrix&       freqC1, 
                       hwMatrix&       freqC2, 
                       const char*     type)
{
    hwChebyshev_I_Design chebyI;
    hwMathStatus status;

    if (!type)
    {
        status(HW_MATH_ERR_NOTSTRING, 8);
    }
    else if (!strcmp(type, "z"))
    {
        status = chebyI.Digital(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                                order, freqC1, freqC2);
    }
    else if (!strcmp(type, "s"))
    {
        status = chebyI.Analog(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                               order, freqC1, freqC2);
    }
    else
    {
        status(HW_MATH_ERR_FILTERTYPE, 8);
    }
    return status;
}
//------------------------------------------------------------------------------
// Designs a Chebyshev type II filter and returns the status
//------------------------------------------------------------------------------
hwMathStatus Cheby2Ord(double      passBandFreq, 
                       double      stopBandFreq,
                       double      passEdgeDb, 
                       double      stopEdgeDb,
                       int&        order, 
                       hwMatrix&   freqC, 
                       const char* type)
{
    hwChebyshev_II_Design chebyII;
    hwMathStatus status;

    if (!type)
    {
        status(HW_MATH_ERR_NOTSTRING, 7);
    }
    else if (!strcmp(type, "z"))
    {
        status = chebyII.Digital(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                                 order, freqC);
    }
    else if (!strcmp(type, "s"))
    {
        status = chebyII.Analog(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                                order, freqC);
    }
    else
    {
        status(HW_MATH_ERR_FILTERTYPE, 7);
    }
    return status;
}
//------------------------------------------------------------------------------
// Designs a Chebyshev type II filter and returns the status
//------------------------------------------------------------------------------
hwMathStatus Cheby2Ord(const hwMatrix& passBandFreq, 
                       const hwMatrix& stopBandFreq,
                       double          passEdgeDb, 
                       double          stopEdgeDb, 
                       int&            order,
                       hwMatrix&       freqC1, 
                       hwMatrix&       freqC2, 
                       const char*     type)
{
    hwChebyshev_II_Design chebyII;
    hwMathStatus status;

    if (!type)
    {
        status(HW_MATH_ERR_NOTSTRING, 8);
    }
    else if (!strcmp(type, "z"))
    {
        status = chebyII.Digital(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                                 order, freqC1, freqC2, true);
    }
    else if (!strcmp(type, "s"))
    {
        status = chebyII.Analog(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                                order, freqC1, freqC2, true);
    }
    else
    {
        status(HW_MATH_ERR_FILTERTYPE, 8);
    }

    return status;
}
//------------------------------------------------------------------------------
// Designs an Elliptic filter and returns the status
//------------------------------------------------------------------------------
hwMathStatus EllipOrd(double      passBandFreq, 
                      double      stopBandFreq,
                      double      passEdgeDb, 
                      double      stopEdgeDb,
                      int&        order, 
                      hwMatrix&   freqC, 
                      const char* type)
{
    hwElliptic_Design ellip;
    hwMathStatus status;

    if (!type)
    {
        status(HW_MATH_ERR_NOTSTRING, 7);
    }
    else if (!strcmp(type, "z"))
    {
        status = ellip.Digital(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                               order, freqC);
    }
    else if (!strcmp(type, "s"))
    {
        status = ellip.Analog(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                              order, freqC);
    }
    else
    {
        status(HW_MATH_ERR_FILTERTYPE, 7);
    }

    return status;
}
//------------------------------------------------------------------------------
// Designs an Elliptic filter and returns the status
//------------------------------------------------------------------------------
hwMathStatus EllipOrd(const hwMatrix& passBandFreq, 
                      const hwMatrix& stopBandFreq,
                      double          passEdgeDb, 
                      double          stopEdgeDb, 
                      int&            order,
                      hwMatrix&       freqC1,
                      hwMatrix&       freqC2, 
                      const char*     type)
{
    hwElliptic_Design ellip;
    hwMathStatus status;

    if (!type)
    {
        status(HW_MATH_ERR_NOTSTRING, 8);
    }
    else if (!strcmp(type, "z"))
    {
        status = ellip.Digital(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                               order, freqC1, freqC2);
    }
    else if (!strcmp(type, "s"))
    {
        status = ellip.Analog(passBandFreq, stopBandFreq, passEdgeDb, stopEdgeDb,
                              order, freqC1, freqC2);
    }
    else
    {
        status(HW_MATH_ERR_FILTERTYPE, 8);
    }

    return status;
}
//------------------------------------------------------------------------------
// Filter a signal using the transfer function and return the status
//------------------------------------------------------------------------------
hwMathStatus Filter(const hwMatrix& numerCoef, 
                    const hwMatrix* denomCoef,
                    const hwMatrix& inSignal, 
                    hwMatrix&       outSignal)
{
    hwFilterManager manager;
    hwMathStatus    status = manager.CreateFilter(numerCoef, denomCoef);
    if (!status.IsOk())
    {
        return status;
    }

    status = manager.ApplyFilter(inSignal, outSignal);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 1)
        {
            status.SetArg1(3);
        }
        else if (status.GetArg1() == 2)
        {
            status.SetArg1(4);
        }
        else if (status == HW_MATH_ERR_FILTERDENZERO)
        {
            status.SetArg1(2);
        }
        else
        {
            status.ResetArgs();
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Filter a signal using the transfer function, making two passes to produce a
// zero phase shift due to two-pass transer function
//------------------------------------------------------------------------------
hwMathStatus FiltFilt(const hwMatrix& numerCoef, 
                      const hwMatrix* denomCoef,
                      const hwMatrix& inSignal,
                      hwMatrix&       outSignal)
{
    hwFilterManager manager;
    hwMathStatus    status = manager.CreateFilter(numerCoef, denomCoef);
    if (!status.IsOk())
    {
        return status;
    }

    status = manager.ApplyZeroPhaseFilter(inSignal, outSignal);
    if (!status.IsOk())
    {
        if (status.GetArg1() == 1)
        {
            status.SetArg1(3);
        }
        else if (status.GetArg1() == 2)
        {
            status.SetArg1(4);
        }
        else if (status == HW_MATH_ERR_FILTERDENZERO)
        {
            status.SetArg1(2);
        }
        else
        {
            status.ResetArgs();
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Down sample a signal by an integer factor
//------------------------------------------------------------------------------
hwMathStatus DownSample(const hwMatrix& inSignal,
                        hwMatrix&       outSignal,
                        int             k,
                        int             phase)
{
    hwFilterManager manager;

    return manager.DownSample(inSignal, outSignal, k, phase);
}
//------------------------------------------------------------------------------
// Up sample a signal by an integer factor
//------------------------------------------------------------------------------
hwMathStatus UpSample(const hwMatrix& inSignal,
                      hwMatrix&       outSignal,
                      int             k,
                      int             phase)
{
    hwFilterManager manager;

    return manager.UpSample(inSignal, outSignal, k, phase);
}
//------------------------------------------------------------------------------
// Compute the response of a filter at specified frequencies using the transfer 
// function
//------------------------------------------------------------------------------
hwMathStatus Response(const hwMatrix& numerCoef, 
                      const hwMatrix* denomCoef,
                      const hwMatrix& freq, 
                      hwMatrix&       response, 
                      const double*   sampFreq)
{
    // Note: freq is assumed to be in Hz for a digital filter and in radians/sec
    // for an analog filter. This was done to support freqz and freqs for hml2. It
    // is a change from the convention in MagRes and PhaseRes. 
    if (!freq.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);
    }

    hwMathStatus status = response.Dimension(freq.M(), freq.N(), hwMatrix::COMPLEX);
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
    }

    int numFreq = freq.Size();

    if (sampFreq)   // digital
    {
        if (*sampFreq < 1.0e-10)
        {
            return status(HW_MATH_ERR_NONPOSITIVE, 5);
        }

        double value = 2.0 * PI / (*sampFreq);      // see hwFilterSpecs for PI vs 2 PI
        hwDigitalFilter filter(numerCoef, denomCoef);

        status = filter.Status();
        if (!status.IsOk())
        {
            return status;
        }

        for (int i = 0; i < numFreq; ++i)
        {
            filter.Response(value * freq(i), response.z(i));

            if (filter.Status() == HW_MATH_WARN_NYQUIST)
            {
                status(HW_MATH_WARN_NYQUIST, 3, 5);
            }
        }
    }
    else    // analog
    {
        hwAnalogFilter filter(numerCoef, denomCoef);

        status = filter.Status();

        if (!status.IsOk())
        {
            return status;
        }

        for (int i = 0; i < numFreq; ++i)
        {
            filter.Response(freq(i), response.z(i));
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute the magnitude response of a filter at specified frequencies using 
// the transfer function
//------------------------------------------------------------------------------
hwMathStatus MagRes(const hwMatrix& numerCoef, 
                    const hwMatrix* denomCoef,
                    const hwMatrix& freq, 
                    hwMatrix&       mag, 
                    const double*   sampFreq)
{
    if (!freq.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);
    }

    hwMathStatus status = mag.Dimension(freq.M(), freq.N(), hwMatrix::REAL);
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
    }

    int numFreq = freq.Size();
    double value;

    if (sampFreq)   // digital
    {
        if (*sampFreq < 1.0e-10)
        {
            return status(HW_MATH_ERR_NONPOSITIVE, 5);
        }

        value = 2.0 * PI / (*sampFreq);      // see hwFilterSpecs for PI vs 2 PI
        hwDigitalFilter filter(numerCoef, denomCoef);

        status = filter.Status();

        if (!status.IsOk())
        {
            return status;
        }

        for (int i = 0; i < numFreq; ++i)
        {
            filter.RespMag(value * freq(i), mag(i));

            if (filter.Status() == HW_MATH_WARN_NYQUIST)
            {
                status(HW_MATH_WARN_NYQUIST, 3, 5);
            }
        }
    }
    else    // analog
    {
        value = 2.0 * PI;
        hwAnalogFilter filter(numerCoef, denomCoef);

        status = filter.Status();

        if (!status.IsOk())
        {
            return status;
        }

        for (int i = 0; i < numFreq; ++i)
        {
            filter.RespMag(value * freq(i), mag(i));
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute the phase response of a filter at specified frequencies using the
// transfer function
//------------------------------------------------------------------------------
hwMathStatus PhaseRes(const hwMatrix& numerCoef,
                      const hwMatrix* denomCoef,
                      const hwMatrix& freq,
                      hwMatrix&       phase,
                      const double*   sampFreq)
{
    if (!freq.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);
    }

    hwMathStatus status = phase.Dimension(freq.M(), freq.N(), hwMatrix::REAL);
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
    }

    int numFreq = freq.Size();
    double value;

    if (sampFreq)   // digital
    {
        if (*sampFreq < 1.0e-10)
        {
            return status(HW_MATH_ERR_NONPOSITIVE, 5);
        }

        value = 2.0 * PI / (*sampFreq);      // see hwFilterSpecs for PI vs 2 PI
        hwDigitalFilter filter(numerCoef, denomCoef);

        status = filter.Status();
        if (!status.IsOk())
        {
            return status;
        }

        for (int i = 0; i < numFreq; ++i)
        {
            filter.RespPhase(value * freq(i), phase(i));

            if (filter.Status() == HW_MATH_WARN_NYQUIST)
            {
                status(HW_MATH_WARN_NYQUIST, 3, 5);
            }
        }
    }
    else    // analog
    {
        value = 2.0 * PI;
        hwAnalogFilter filter(numerCoef, denomCoef);

        status = filter.Status();
        if (!status.IsOk())
        {
            return status;
        }

        for (int i = 0; i < numFreq; ++i)
        {
            filter.RespPhase(value * freq(i), phase(i));
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Compute the impulse response of a filter at specified frequencies using the
// transfer function
//------------------------------------------------------------------------------
hwMathStatus ImpulseRes(const hwMatrix& numerCoef, 
                        const hwMatrix* denomCoef,
                        hwMatrix&       impulseRes, 
                        hwMatrix*       time,
                        int             numPnts, 
                        double          sampFreq)
{
    hwFilterManager manager;
    hwMathStatus    status = manager.CreateFilter(numerCoef, denomCoef);
    if (!status.IsOk())
    {
        return status;
    }

    if (numPnts < 0)
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 5);
    }

    if (sampFreq < 0.0)
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 6);
    }

    // compute number of samples to report
    if (numPnts == 0)
    {
        if (denomCoef && denomCoef->Size() > 1)
        {
            hwMatrix roots;

            status = PolyRoots(*denomCoef, roots);
            if (!status.IsOk())
            {
                status.ResetArgs();
                return status;
            }

            // find pole with maximum magnitude
            int    numRoots = roots.Size();
            int    index;
            double angFreq;     // angular frequency
            double temp;
            double mag       = -999.0;
            double precision = 1e-6;

            for (int i = 0; i < numRoots; ++i)
            {
                temp = (roots.IsReal()) ? fabs(roots(i)) : roots.z(i).Mag();
                if (temp > mag)
                {
                    index = i;
                    mag   = temp;
                }
            }

            // find angular frequency of max magnitude pole
            if (roots.IsReal())
            {
                angFreq = (roots(index) < 0.0) ? PI : 0.0;
            }
            else
            { 
                angFreq = roots.z(index).Arg();
            }

            // set number of output samples
            if (mag > 1.0 + precision)
            {
                // unstable filter - stop at 120 dB amplification
                if (fabs(angFreq) < precision)
                {
                    numPnts = (int) floor(6.0/log10(mag));  // dc
                }
                else if (angFreq > PI - precision)
                {
                    numPnts = (int) floor(6.0/log10(mag));  // Nyquist
                }
                else
                {
                    numPnts = (int) floor(3.0/log10(mag));  // complex conjugate pair
                }
            }
            else if (mag < 1.0 - precision)
            {
                // stable filter - stop at -120 dB attenuation
                if (fabs(angFreq) < precision)
                {
                    numPnts = (int) floor(-6.0/log10(mag));  // dc
                }
                else if (angFreq > PI - precision)
                {
                    numPnts = (int) floor(-6.0/log10(mag));  // Nyquist
                }
                else
                {
                    numPnts = (int) floor(-3.0/log10(mag));  // complex conjugate pair
                }
            }
            else
            {
                // stop after five periods of lowest frequency
                double angFreq_min = PI;

                for (int i = 0; i < numRoots; ++i)
                {
                    if (roots.IsReal())
                    {
                        if (fabs(roots(i) - mag) < precision)
                        {
                            if (fabs(roots(i)) < precision)
                            {
                                angFreq_min = 0.0;
                            }
                        }
                    }
                    else
                    {
                        if (fabs(roots.z(i).Mag() - mag) < precision)
                        {
                            if (fabs(roots.z(i).Arg()) < angFreq_min)
                            {
                                angFreq_min = fabs(roots.z(i).Arg());
                            }
                        }
                    }
                }

                numPnts = (int) ceil(10.0*PI/angFreq_min);
            }
        }

        numPnts += numerCoef.Size();
    }

    // filter the impulse and add the time vector
    hwMatrix signal(numPnts, hwMatrix::REAL);

    signal.SetElements(0.0);
    signal(0) = 1.0;

    status = Filter(numerCoef, denomCoef, signal, impulseRes);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    if (time)
    {
        status = time->Dimension(impulseRes.Size(), hwMatrix::REAL);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        double ts = 1.0 / sampFreq;

        for (int i = 0; i < time->Size(); ++i)
        {
            (*time)(i) = i * ts;
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Estimate IIR transfer function coefficients from frequency response information
//------------------------------------------------------------------------------
hwMathStatus FilterID(const hwMatrix& response, 
                      const hwMatrix& freq,
                      hwMatrix&       numerCoef, 
                      hwMatrix&       denomCoef,
                      const hwMatrix* weight, 
                      const char*     type)
{
    if (response.IsReal())
    {
        hwMatrix temp;
        temp.PackComplex(response);
        return FilterID(temp, freq, numerCoef, denomCoef, weight, type);
    }

    hwMathStatus status;
    int m = response.Size();

    if (!response.IsVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }

    if (!freq.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 2);
    }

    if (!freq.IsVector())
    {
        return status(HW_MATH_ERR_VECTOR, 2);
    }

    if (freq.Size() != m)
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    if (weight)
    {
        if (weight->Size() != m)
        {
            return status(HW_MATH_ERR_ARRAYSIZE, 1, 5);
        }
        if (!weight->IsReal())
        {
            return status(HW_MATH_ERR_COMPLEX, 5);
        }
    }

    if (strcmp(type, "s") && strcmp(type, "z"))
    {
        return status(HW_MATH_ERR_FILTERTYPE, 6);
    }

    int nb = numerCoef.Size();
    int na = denomCoef.Size();

    if (nb < 1)
    {
        return status(HW_MATH_ERR_FILTERORDER, 3);
    }

    if (na < 1)
    {
        return status(HW_MATH_ERR_FILTERORDER, 4);
    }

    int n = (na - 1) + nb;
    hwMatrix LHS(2*m, n, hwMatrix::REAL);
    hwMatrix RHS(2*m, 1, hwMatrix::REAL);
    hwMatrix Ab(2*m, n+1, hwMatrix::REAL);   // augmented matrix
    hwMatrix Q;
    hwMatrix R;
    hwComplex zk;

    if (weight)
    {
        for (int i = 0; i < m; ++i)
        {
            hwComplex c(0.0, freq(i));
            hwComplex h = response.z(i);
            double w    = sqrt((*weight)(i));

            for (int j = 0; j < nb; ++j)
            {
                if (!strcmp(type, "s"))
                {
                    zk = hwComplex::pow(c, nb-1-j);
                }
                else    // type == "z"
                {
                    zk = hwComplex::exp(c * (nb-1-j));
                }
                
                Ab(i, j) = zk.Real() * w;
                Ab(i+m, j) = zk.Imag() * w;
            }

            for (int j = nb; j < n; ++j)
            {
                if (!strcmp(type, "s"))
                {
                    zk = hwComplex::pow(c, n-1-j+1) * h;
                }
                else    // type == "z"
                {
                    zk = hwComplex::exp(c * (n-1-j+1)) * h;
                }
                Ab(i, j) = -zk.Real() * w;
                Ab(i+m, j) = -zk.Imag() * w;
            }

            Ab(i, n)   = h.Real() * w;
            Ab(i+m, n) = h.Imag() * w;
        }
    }
    else
    {
        for (int i = 0; i < m; ++i)
        {
            hwComplex c(0.0, freq(i));
            hwComplex h = response.z(i);

            for (int j = 0; j < nb; ++j)
            {
                if (!strcmp(type, "s"))
                {
                    zk = hwComplex::pow(c, nb-1-j);
                }
                else    // type == "z"
                {
                    zk = hwComplex::exp(c * (nb-1-j));
                }
                Ab(i, j) = zk.Real();
                Ab(i+m, j) = zk.Imag();
            }

            for (int j = nb; j < n; ++j)
            {
                if (!strcmp(type, "s"))
                {
                    zk = hwComplex::pow(c, n-1-j+1) * h;
                }
                else    // type == "z"
                {
                    zk = hwComplex::exp(c * (n-1-j+1)) * h;
                }
                Ab(i, j)   = -zk.Real();
                Ab(i+m, j) = -zk.Imag();
            }

            Ab(i, n)   = h.Real();
            Ab(i+m, n) = h.Imag();
        }
    }

    // Minimize sum|response(z) - numer(z)/denom(z)|^2, ultimately with Steiglitz-McBride.
    // The first step is minimizing an approximation sum|denom(z) - response(z)*numer(z)|^2. 
    // This is converted into a linear least squares problem Ax=b, solved via QR.
    // Instead of computing A = QR, compute qr([A b]) = [Q q][R r; 0 e]. 
    // Assuming that b is in the column space of A we can compute x = R^(-1)r 
    // without Q. Substituting the solved x into QRx = Qr + qe reveals e to be 
    // the error in the fit. See "A Note on the MGS Approach" in Golub and 
    // Van Loan, 5.3.5.

    status = Ab.QR(Q, R);   // R(n,n) = e
    if (!status.IsOk())
    {
        status.ResetArgs();

        if (!status.IsWarning())
        {
            return status();
        }
    }

    hwMatrix RR(n, n, hwMatrix::REAL);
    hwMatrix r(n, 1, hwMatrix::REAL);
    hwMatrix X;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            RR(i, j) = R(i, j);
        }
        r(i, 0) = R(i, n);
    }

    status = X.LSolve(RR, r);   // needs a trangular solver

    if (!status.IsOk())
    {
        status.ResetArgs();

        if (!status.IsWarning())
        {
            return status();
        }
        status = hwMathStatus();    // suppress warning for now
    }

    for (int i = 0; i < nb; ++i)
    {
        numerCoef(i) = X(i);
    }

    for (int i = 1; i < na; ++i)
    {
        denomCoef(i-1) = X(i-1+nb);
    }
    denomCoef(na-1) = 1.0;

    // normalize to denomCoef(0) to match butter, cheb1, etc
    for (int i = 0; i < nb; ++i)
    {
        numerCoef(i) /= denomCoef(0);
    }

    for (int i = 1; i < na; ++i)
    {
        denomCoef(i) /= denomCoef(0);
    }
    denomCoef(0) = 1.0;

    return status;
}
//------------------------------------------------------------------------------
// Compute Chebyshev window weights and return status
//------------------------------------------------------------------------------
hwMathStatus ChebyWin(hwMatrix&   weight,
                      double      sideLobe,
                      const char* type)
{
    if (sideLobe < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_DB_SIGN, 2);
    }
    
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 3);
    }
    
    bool periodic;
    if (!strcmp(type, "symmetric"))
    {
        periodic = false;
    }
    else if (!strcmp(type, "periodic"))
    {
        periodic = true;
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 3);
    }

    hwChebyshev window(sideLobe, periodic);
    return window.ComputeWeights(weight);

}
//------------------------------------------------------------------------------
// Compute Hamming window weights and return status
//------------------------------------------------------------------------------
hwMathStatus HammWin(hwMatrix&   weight,
                     const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 2);
    }

    bool periodic;
    if (!strcmp(type, "symmetric"))
    {
        periodic = false;
    }
    else if (!strcmp(type, "periodic"))
    {
        periodic = true;
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    hwHamming window(periodic);
    return window.ComputeWeights(weight);
}
//------------------------------------------------------------------------------
// Compute Hann window weights and return status
//------------------------------------------------------------------------------
hwMathStatus HannWin(hwMatrix&   weight,
                     const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 2);
    }

    bool periodic;
    if (!strcmp(type, "symmetric"))
    {
        periodic = false;
    }
    else if (!strcmp(type, "periodic"))
    {
        periodic = true;
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    hwHanning window(periodic);
    return window.ComputeWeights(weight);
}
//------------------------------------------------------------------------------
// Compute Bartlett-Hann window weights and return status
//------------------------------------------------------------------------------
hwMathStatus BartHannWin(hwMatrix&   weight,
                         const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 2);
    }

    bool periodic;
    if (!strcmp(type, "symmetric"))
    {
        periodic = false;
    }
    else if (!strcmp(type, "periodic"))
    {
        periodic = true;
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    hwBartlettHann window(periodic);
    return window.ComputeWeights(weight);
}
//------------------------------------------------------------------------------
// Compute Blackman window weights and return status
//------------------------------------------------------------------------------
hwMathStatus BlackmanWin(hwMatrix&   weight,
                         const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 2);
    }

    bool periodic;
    if (!strcmp(type, "symmetric"))
    {
        periodic = false;
    }
    else if (!strcmp(type, "periodic"))
    {
        periodic = true;
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    hwBlackman window(periodic);
    return window.ComputeWeights(weight);
}
//------------------------------------------------------------------------------
// Compute Welch window weights and return status
//------------------------------------------------------------------------------
hwMathStatus WelchWin(hwMatrix&   weight,
                      const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 2);
    }

    bool periodic;
    if (!strcmp(type, "symmetric"))
    {
        periodic = false;
    }
    else if (!strcmp(type, "periodic"))
    {
        periodic = true;
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    hwWelch window(periodic);
    return window.ComputeWeights(weight);
}
//------------------------------------------------------------------------------
// Compute Parzen window weights and return status
//------------------------------------------------------------------------------
hwMathStatus ParzenWin(hwMatrix&   weight,
                       const char* type)
{
    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 2);
    }

    bool periodic;
    if (!strcmp(type, "symmetric"))
    {
        periodic = false;
    }
    else if (!strcmp(type, "periodic"))
    {
        periodic = true;
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2);
    }

    hwParzen window(periodic);
    return window.ComputeWeights(weight);
}
//------------------------------------------------------------------------------
// Compute Kaiser-Bessel window weights and return status
//------------------------------------------------------------------------------
hwMathStatus KaiserBesselWin(hwMatrix&   weight,
                             double      beta,
                             const char* type)
{
    if (beta < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 2);
    }

    if (!type)
    {
        return hwMathStatus(HW_MATH_ERR_NOTSTRING, 3);
    }

    bool periodic;
    if (!strcmp(type, "symmetric"))
    {
        periodic = false;
    }
    else if (!strcmp(type, "periodic"))
    {
        periodic = true;
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 3);
    }

    hwKaiserBessel window(beta, periodic);
    return window.ComputeWeights(weight);
}
//------------------------------------------------------------------------------
// Acoustic A weight magnitude function
//------------------------------------------------------------------------------
hwMathStatus dBa(const hwMatrix& freq, 
                 const hwMatrix& mag_in,
                 hwMatrix&       mag_out,
                 double          reference)
{
    hwMathStatus status;

    if (!freq.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!freq.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }

    if (!mag_in.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 2);
    }

    if (!mag_in.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 2);
    }

    if (freq.Size() != mag_in.Size())
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    status = mag_out.Dimension(mag_in.M(), mag_in.N(), hwMatrix::REAL);
    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(3);
        }
        return status;
    }

    if (reference <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 4);
    }

    double f4;
    double f2;
    const double* x          = freq.GetRealData();
    const double* y          = mag_in.GetRealData();
    double*       result     = mag_out.GetRealData();
    long          num_points = static_cast<long>(freq.Size());

    double n1 = 12200.0 * 12200.0;
    double n2 = 20.6 * 20.6;
    double n3 = 107.7 * 107.7;
    double n4 = 737.9 *737.9;

    for (long j = 0; j < num_points; ++j)
    {
        f2 = x[j] * x[j];
        f4 = f2 * f2;
 
        result[j] = (n1*f4)/((f2+n2)*(f2+n1)*sqrt(f2+n3)*sqrt(f2+n4));
        result[j] /= 0.79434639580229505 * reference;
        result[j] *= y[j];
    }

    return status;
}
//------------------------------------------------------------------------------
// Acoustic B weight magnitude function
//------------------------------------------------------------------------------
hwMathStatus dBb(const hwMatrix& freq, 
                 const hwMatrix& mag_in, 
                 hwMatrix&       mag_out,
                 double          reference)
{
    hwMathStatus status;

    if (!freq.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!freq.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }

    if (!mag_in.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 2);
    }

    if (!mag_in.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 2);
    }

    if (freq.Size() != mag_in.Size())
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    status = mag_out.Dimension(mag_in.M(), mag_in.N(), hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(3);
        }
        return status;
    }

    if (reference <= 0.0)
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 4);
    }

    double f3;
    double f2;
    const double* x = freq.GetRealData();
    const double* y = mag_in.GetRealData();
    double* result     = mag_out.GetRealData();
    long    num_points = static_cast<long>(freq.Size());

    double n1 = 12200.0 * 12200.0;
    double n2 = 20.6 * 20.6;
    double n3 = 158.5 * 158.5;

    for (long j = 0; j < num_points; ++j)
    {
        f2 = x[j] * x[j];
        f3 = f2   * x[j];
 
        result[j] = (n1*f3)/((f2+n2)*(f2+n1)*sqrt(f2+n3));
        result[j] /= 0.98066304307189356 * reference;
        result[j] *= y[j];
    }

    return status;
}
//------------------------------------------------------------------------------
// Acoustic C weight magnitude function
//------------------------------------------------------------------------------
hwMathStatus dBc(const hwMatrix& freq, 
                 const hwMatrix& mag_in, 
                 hwMatrix&       mag_out,
                 double          reference)
{
    hwMathStatus status;

    if (!freq.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!freq.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }

    if (!mag_in.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 2);
    }

    if (!mag_in.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 2);
    }

    if (freq.Size() != mag_in.Size())
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    status = mag_out.Dimension(mag_in.M(), mag_in.N(), hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(3);
        }
        return status;
    }

    if (reference <= 0.0)
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 4);
    }

    double f2;
    const double* x = freq.GetRealData();
    const double* y = mag_in.GetRealData();
    double* result     = mag_out.GetRealData();
    long    num_points = static_cast<long>(freq.Size());

    double n1 = 12200.0 * 12200.0;
    double n2 = 20.6 * 20.6;

    for (long j = 0; j < num_points; ++j)
    {
        f2 = x[j] * x[j];
 
        result[j] = (n1*f2)/((f2+n2)*(f2+n1));
        result[j] /= 0.9929048655202054 * reference;
        result[j] *= y[j];
    }
    return status;
}
//------------------------------------------------------------------------------
// Acoustic U weight magnitude function
//------------------------------------------------------------------------------
hwMathStatus dBu(const hwMatrix& freq, 
                 const hwMatrix& mag_in, 
                 hwMatrix&       mag_out,
                 double          reference)
{
    hwMathStatus status;

    if (!freq.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!freq.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }

    if (!mag_in.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 2);
    }

    if (!mag_in.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 2);
    }

    if (freq.Size() != mag_in.Size())
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);
    }

    status = mag_out.Dimension(mag_in.M(), mag_in.N(), hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
        {
            status.SetArg1(3);
        }
        return status;
    }

    if (reference <= 0.0)
    {
        return status(HW_MATH_ERR_NONPOSITIVE, 4);
    }

    double f2;
    double f;

    hwMatrix* X = (hwMatrix*) &freq;    // override const
    hwMatrix* Y = (hwMatrix*) &mag_in;

    const double* x = X->GetRealData();
    const double* y = Y->GetRealData();
    double* result     = mag_out.GetRealData();
    long    num_points = static_cast<long>(freq.Size());

    double n1 = 12200.0 * 12200.0;
    double n3 = 8800.0;
    double n4 = 7850.0 * 7850.0;
    double n5 = 12150.0;
    double n6 = 2900.0 * 2900.0;

    for (long j = 0; j < num_points; ++j)
    {
        f = x[j];
        f2 = x[j] * x[j];
 
        result[j] = n1 / (f2 + n1);
        result[j] /= sqrt((f + n3) * (f + n3) + n4) * 
                        sqrt((f - n3) * (f - n3) + n4);
        result[j] /= sqrt((f + n5) * (f + n5) + n6) * 
                        sqrt((f - n5) * (f - n5) + n6);
        result[j] /= 4.6078641892082647e-17 * reference;
        result[j] *= y[j];
    }

    return status;
}
//------------------------------------------------------------------------------
// Find local extrema of a signal
//------------------------------------------------------------------------------
hwMathStatus FindPeaks(const hwMatrix& signal,
                       bool            twoSided,
                       double          minPeakHeight,
                       int             minPeakDistance,
                       int             minPeakWidth,
                       hwMatrix&       peaks,
                       hwMatrixI&      index,
                       int             indexOrigin,
                       PeakInfo*       extra)
{
    if (!signal.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);

    if (!signal.IsEmptyOrVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

    hwMathStatus status;
    hwMatrix signalD;

    // remove mean
    if (twoSided)
    {
        status = Detrend(signal, "constant", signalD);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }
    }

    // find peak locations
    int n = signal.Size();
    int numPeaks = 0;
    hwMatrixI indexRaw(0, 1, hwMatrixI::REAL);

    for (int i = 1; i < n-1; ++i)
    {
        if ((signal(i-1) <= signal(i) && signal(i) > signal(i+1)) ||
            (twoSided && (signal(i-1) >= signal(i) && signal(i) < signal(i+1))))
        {
            ++numPeaks;

            if (numPeaks == 1)
                indexRaw.Dimension(1, 1, hwMatrixI::REAL);
            else
                indexRaw.Resize(numPeaks, 1);

            indexRaw(numPeaks-1) = i;
        }
    }

    // apply minumum height criterion
    if (minPeakHeight == -1.0)
    {
        minPeakHeight = MACHEP2;
    }

    if (minPeakHeight > 0)
    {
        if (twoSided)
        {
            for (int i = 0; i < numPeaks; ++i)
            {
                if (fabs(signalD(indexRaw(i))) < minPeakHeight)
                {
                    indexRaw.DeleteRows(i);
                    --i;
                    --numPeaks;
                }
            }
        }
        else
        {
            for (int i = 0; i < numPeaks; ++i)
            {
                if (signal(indexRaw(i)) < minPeakHeight)
                {
                    indexRaw.DeleteRows(i);
                    --i;
                    --numPeaks;
                }
            }
        }
    }

    // apply minumum distance criterion
    if (minPeakDistance == -1)
        minPeakDistance = 4;

    if (numPeaks && minPeakDistance > 0)
    {
        hwMatrix signalPeaks(numPeaks, 1, hwMatrix::REAL);
        hwMatrix dummy;
        hwMatrixI indexSorted;

        for (int i = 0; i < numPeaks; ++i)
            signalPeaks(i) = signal(indexRaw(i));

        status = Sort(signalPeaks, dummy, &indexSorted, false);

        for (int i = 0; i < numPeaks-1; ++i)
        {
            int peakIndex = indexRaw(indexSorted(i));

            for (int j = i + 1; j < numPeaks; ++j)
            {
                if (abs(indexRaw(indexSorted(j)) - peakIndex) < minPeakDistance)
                {
                    status = indexSorted.DeleteRows(j);
                    --j;
                    --numPeaks;

                    if (!status.IsOk())
                    {
                        status.ResetArgs();
                        return status;
                    }
                }
            }
        }

        // unsort indices with respect to peak values
        hwMatrix unsort(numPeaks, hwMatrix::REAL);
        hwMatrix sorted(numPeaks, hwMatrix::REAL);

        for (int i = 0; i < numPeaks; ++i)
            sorted(i) = static_cast<double>(indexRaw(indexSorted(i)));

        status = Sort(sorted, unsort, nullptr, true);   // templatize to avoid copy/casts

        status = index.Dimension(numPeaks, 1, hwMatrixI::REAL);

        for (int i = 0; i < numPeaks; ++i)
            index(i) = static_cast<int>(unsort(i));
    }
    else
    {
        index = indexRaw;
    }

    // apply minumum width criterion
    if (minPeakWidth == -1)
        minPeakWidth = 2;

    if (extra)
    {
        extra->parabol_pp = new hwMatrix(0, 3, hwMatrix::REAL);
        extra->parabol_x = new hwMatrix(0, 2, hwMatrix::REAL);
        extra->height = new hwMatrix(1, 0, hwMatrix::REAL);
        extra->baseline = new hwMatrix(1, 0, hwMatrix::REAL);
        extra->roots = new hwMatrix(0, 2, hwMatrix::REAL);
    }

    if (numPeaks && minPeakWidth > 0)
    {
        for (int i = 0; i < numPeaks; ++i)
        {
            int left  = _max(index(i) - (minPeakDistance+1)/2, 0);
            int right = _min(index(i) + (minPeakDistance+1)/2, n-1);
            int size = right - left + 1;
            hwMatrix indexD(size, hwMatrix::REAL);
            hwMatrix neighbors(size, hwMatrix::REAL);
            hwMatrix coef;
            int idx = left;

            for (int j = 0; j < size; ++j)
            {
                indexD(j) = static_cast<double>(idx + indexOrigin);

                if (twoSided)
                    neighbors(j) = signalD(idx++);
                else
                    neighbors(j) = signal(idx++);
            }

            status = PolyCurveFit(indexD, neighbors, 2, coef);    // ascending coef

            double a = coef(2);
            double b = coef(1);
            double c = coef(0);
            double height = c - b*b/(4.0*a);
            double zero1;
            double zero2;
            bool isOk;

            if (height > minPeakHeight) // maxima peak for a good parabola fit
                isOk = quadraticRoots(a, b, c-(height+minPeakHeight)/2.0, zero1, zero2);
            else                        // minima peak for a good parabola fit
                isOk = quadraticRoots(a, b, c+(height-minPeakHeight)/2.0, zero1, zero2);

            double width = abs(zero2 - zero1);

            if (!isOk || width < minPeakWidth)
            {
                status = index.DeleteRows(i);
                --i;
                --numPeaks;

                if (!status.IsOk())
                {
                    status.ResetArgs();
                    return status;
                }
            }
            else if (extra)
            {
                if (i == 0)
                {
                    extra->parabol_pp->Dimension(1, 3, hwMatrix::REAL);
                    extra->parabol_x->Dimension(1, 2, hwMatrix::REAL);
                    extra->height->Dimension(1, 1, hwMatrix::REAL);
                    extra->baseline->Dimension(1, 1, hwMatrix::REAL);
                    extra->roots->Dimension(1, 2, hwMatrix::REAL);
                }
                else
                {
                    extra->parabol_pp->Resize(i+1, 3);
                    extra->parabol_x->Resize(i+1, 2);
                    extra->height->Resize(1, i+1);
                    extra->baseline->Resize(1, i+1);
                    extra->roots->Resize(i+1, 2);
                }

                (*extra->parabol_pp)(i, 0) = a;
                (*extra->parabol_pp)(i, 1) = b;
                (*extra->parabol_pp)(i, 2) = c;
                (*extra->parabol_x)(i, 0) = left + indexOrigin;
                (*extra->parabol_x)(i, 1) = right + indexOrigin;
                (*extra->height)(0, i) = height;
                (*extra->baseline)(0, i) = (height+minPeakHeight) / 2.0;

                if (zero1 < zero2)
                {
                    (*extra->roots)(i, 0) = zero1;
                    (*extra->roots)(i, 1) = zero2;
                }
                else
                {
                    (*extra->roots)(i, 0) = zero2;
                    (*extra->roots)(i, 1) = zero1;
                }
            }
        }
    }

    // populate peaks argument
    status = peaks.Dimension(numPeaks, 1, hwMatrix::REAL);

    for (int i = 0; i < numPeaks; ++i)
        peaks(i) = signal(index(i));

    if (indexOrigin)
        index += indexOrigin;

    if (signal.M() == 1 && signal.N() != 1)
    {
        peaks.Transpose();
        index.Transpose();
    }

    return status;
}