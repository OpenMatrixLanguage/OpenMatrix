/**
* @file PolynomFuncs.cxx
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
#include "PolynomFuncs.h"

#include "MathUtilsFuncs.h"
#include "hwMatrix.h"
#include "hwMatrixN.h"
//#include "hwMathStatus.h"

//------------------------------------------------------------------------------
// Computes polynomial roots and returns status
//------------------------------------------------------------------------------
hwMathStatus PolyRoots(const hwMatrix& coef,
                       hwMatrix&       roots)
{
	// check inputs
	int i, j;
	int numRoots = coef.Size()-1;
	int numLeadZeros = 0;
	hwMathStatus status;

	if (!coef.IsEmptyOrVector())
		return status(HW_MATH_ERR_VECTOR, 1);

    // count leading zeros, if any
	if (coef.IsReal())
	{
		for (i = 0; i <= numRoots; ++i)
		{
			if (coef(i) == 0.0)
				++numLeadZeros;
			else
				break;
		}
	}
	else    // !coef.IsReal()
	{
		for (i = 0; i <= numRoots; ++i)
		{
			if (coef.z(i) == 0.0)
				++numLeadZeros;
			else
				break;
		}
	}

    if (numLeadZeros >= numRoots)
    {
        if (numRoots == -1)
            roots.Dimension(0, 0, hwMatrix::REAL);
        else if (coef.IsReal() && coef(numRoots) != 0.0)
            roots.Dimension(0, 1, hwMatrix::REAL);
        else if (!coef.IsReal() && coef.z(numRoots) != 0.0)
            roots.Dimension(0, 1, hwMatrix::REAL);
        else    // infinite number of roots
            roots.Dimension(0, 0, hwMatrix::REAL);

        return status;
    }

	numRoots -= numLeadZeros;

	// setup the problem as an Eigen problem
	if (coef.IsReal())
	{
		hwMatrix A(numRoots, numRoots, hwMatrix::REAL);

		if (A.IsEmpty())
			return status(HW_MATH_ERR_OUTOFMEMORY);     // allocation failure

		// find the roots
		A.SetElements(0.0);

		for (i = 1; i < numRoots; i++)
			A(i, i-1) = 1.0;

		for (j = 0; j < numRoots; j++)
			A(0, j) = -coef(j+numLeadZeros+1) / coef(numLeadZeros);

		status = A.Eigen(true, NULL, roots);
	}
	else    // !coef.IsReal()
	{
		hwMatrix A(numRoots, numRoots, hwMatrix::COMPLEX);

		if (A.IsEmpty())
			return status(HW_MATH_ERR_OUTOFMEMORY);     // allocation failure

		// find the roots
		A.SetElements(0.0);

		for (i = 1; i < numRoots; i++)
			A.z(i, i-1) = 1.0;

		for (j = 0; j < numRoots; j++)
			A.z(0, j) = -coef.z(j+numLeadZeros+1) / coef.z(numLeadZeros);

		status = A.Eigen(true, NULL, roots);
	}

	return status;
}
//------------------------------------------------------------------------------
// Evaluates a polynomial and returns status
//------------------------------------------------------------------------------
hwMathStatus PolyVal(const hwMatrix& coef,
                     double          x,
                     double&         y)
{
	if (!coef.IsReal())
		return hwMathStatus(HW_MATH_ERR_COMPLEXSUPPORT, 1);

	if (!coef.IsEmptyOrVector())
		return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

	int numCoefs = coef.Size();
	if (numCoefs == 0) 
	{ 
		y = 0.; 
	} 
	else
	{
		int i = coef.Size() - 1;
		y = coef(i);

		while (i)
		{
			i--;
			y = y * x + coef(i);
		}
	}
	return hwMathStatus();
}
//------------------------------------------------------------------------------
// Evaluates a polynomial and returns status
//------------------------------------------------------------------------------
hwMathStatus PolyVal(const hwMatrix&  coef,
                     const hwComplex& x,
                     hwComplex&       y)
{
	if (!coef.IsEmptyOrVector())
		return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

	int i;
	int numCoefs = coef.Size();
	if (numCoefs == 0) 
	{
		// PolyVal([],i) returns 0
		y = 0;
		return hwMathStatus();
	}
	else
	{
		i = coef.Size() - 1;
	}

	if (coef.IsReal())
	{
		y = coef(i);

		while (i)
		{
			i--;
			y = y * x + coef(i);
		}
	}
	else
	{
		y = coef.z(i);

		while (i)
		{
			i--;
			y = y * x + coef.z(i);
		}
	}

	return hwMathStatus();
}
//------------------------------------------------------------------------------
// Evaluates a polynomial and returns status
//------------------------------------------------------------------------------
hwMathStatus PolyVal(const hwMatrix& coef,
                     const hwMatrix& x,
                     hwMatrix&       y)
{
	int i;
	int j;
	int numPnts = x.Size();
	int numCoefs = coef.Size();
	hwMathStatus status;

	if (!coef.IsEmptyOrVector())
		return status(HW_MATH_ERR_VECTOR, 1);

	if (!x.IsEmptyOrVector())
		return status(HW_MATH_ERR_VECTOR, 2);

	if (coef.IsReal() && x.IsReal())
		status = y.Dimension(x.M(), x.N(), hwMatrix::REAL);
	else
		status = y.Dimension(x.M(), x.N(), hwMatrix::COMPLEX);

	if (!status.IsOk())
	{
		if (status.GetArg1() == 0)
			status.SetArg1(3);
		else
			status.ResetArgs();

		return status;
	}

	if (coef.IsEmpty())
	{
		for (i = 0; i < numPnts; i++)
		{
			if (x.IsReal())
			{
				y(i) = 0.; 
			}
			else
			{
				y.z(i) = 0.;
			}
		}
	}
	else if (coef.IsReal())
	{
		if (x.IsReal())
		{
			double xx, yy;

			for (i = 0; i < numPnts; i++)
			{
				xx = x(i);
				j = numCoefs - 1;
				yy = coef(j);
				while (j)
				{
					j--;
					yy = yy * xx + coef(j);
				}

				y(i) = yy;
			}
		}
		else // !x.IsReal()
		{
			hwComplex xx, yy;

			for (i = 0; i < numPnts; i++)
			{
				xx = x.z(i);
				j = numCoefs - 1;
				yy = coef(j);
				while (j)
				{
					j--;
					yy = yy * xx + coef(j);
				}

				y.z(i) = yy;
			}
		}
	}
	else // !coef.IsReal()
	{
		if (x.IsReal())
		{
			double xx;
			hwComplex yy;

			for (i = 0; i < numPnts; i++)
			{
				xx = x(i);
				j = numCoefs - 1;
				yy = coef.z(j);
				while (j)
				{
					j--;
					yy = yy * xx + coef.z(j);
				}

				y.z(i) = yy;
			}
		}
		else // !x.IsReal()
		{
			hwComplex xx, yy;

			for (i = 0; i < numPnts; i++)
			{
				xx = x.z(i);
				j = numCoefs - 1;
				yy = coef.z(j);
				while (j)
				{
					j--;
					yy = yy * xx + coef.z(j);
				}

				y.z(i) = yy;
			}
		}
	}

	return status;
}
//------------------------------------------------------------------------------
// Performs polynomial division or deconvolution and returns status
//------------------------------------------------------------------------------
hwMathStatus PolyDivide(const hwMatrix& A,
                        const hwMatrix& B,
	                    hwMatrix&       Q,
                        hwMatrix&       R)
{
    // see deconv
    // return Q and R such that A = ConvLin(Q,B) + R
    hwMathStatus status;

    if (!A.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!B.IsVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    int asize = A.Size();
    int bsize = B.Size();

    // check B for leading zeros
    int numLeadZeros = 0;

    if (B.IsReal())
    {
        for (int i = 0; i < bsize; i++)
        {
            if (B(i) == 0.0)
                ++numLeadZeros;
            else
                break;
        }
    }
    else
    {
        for (int i = 0; i < bsize; i++)
        {
            if (B.z(i) == 0.0)
                ++numLeadZeros;
            else
                break;
        }
    }

    if (numLeadZeros == bsize)
        return status(HW_MATH_ERR_DIVIDEZERO, 2);

    bsize -= numLeadZeros;

    // generate Q and R
    R = A;  // initialize remainder to use as working memory

    if (bsize > asize)
    {
        status = Q.Dimension(1, hwMatrix::REAL);
        Q(0) = 0.0;
        return status;
    }

    int qn = asize-bsize+1;

    if (A.IsReal() && B.IsReal())
    {
        if (A.M() == 1)
            status = Q.Dimension(1, qn, hwMatrix::REAL);
        else
            status = Q.Dimension(qn, 1, hwMatrix::REAL);

        if (!status.IsOk())
        {
            status.SetArg1(3);
            return status;
        }

        for (int i = 0; i < qn; ++i)
        {
            Q(i) = R(i) / B(numLeadZeros);
            R(i) = 0.0;

            for (int j = 1; j < bsize; ++j)
                R(i+j) -= Q(i) * B(numLeadZeros+j);
        }
    }
    else if (!A.IsReal() && !B.IsReal())
    {
        if (A.M() == 1)
            status = Q.Dimension(1, qn, hwMatrix::COMPLEX);
        else
            status = Q.Dimension(qn, 1, hwMatrix::COMPLEX);

        if (!status.IsOk())
        {
            status.SetArg1(3);
            return status;
        }

        for (int i = 0; i < qn; ++i)
        {
            Q.z(i) = R.z(i) / B.z(numLeadZeros);
            R.z(i) = 0.0;

            for (int j = 1; j < bsize; ++j)
                R.z(i+j) -= Q.z(i) * B.z(numLeadZeros+j);
        }
    }
    else if (A.IsReal() && !B.IsReal())
    {
        hwMatrix temp;
        temp.PackComplex(A, nullptr);

        return PolyDivide(temp, B, Q, R);
    }
    else if (!A.IsReal() && B.IsReal())
    {
        hwMatrix temp;
        temp.PackComplex(B, nullptr);

        return PolyDivide(A, temp, Q, R);
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes the derivative of a polynomial and returns status
//------------------------------------------------------------------------------
hwMathStatus PolyDer(const hwMatrix& P,
                     hwMatrix&       D)
{
    int i;
    int size = P.Size() - 1;
	hwMathStatus status;

    if (!P.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!size)
    {
        status = D.Dimension(1, hwMatrix::REAL);
        D(0) = 0.0;
        return status;
    }

    if (P.IsReal())
        status = D.Dimension(size, hwMatrix::REAL);
    else
        status = D.Dimension(size, hwMatrix::COMPLEX);

	if (!status.IsOk())
	{
		status.SetArg1(2);
		return status;
	}

    if (P.IsReal())
    {
        for (i = 0; i < size; ++i)
            D(i) = (double) (size-i) * P(i);
    }
    else
    {
        for (i = 0; i < size; ++i)
            D.z(i) = (double) (size-i) * P.z(i);
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes the derivative of a polynomial product and returns status
//------------------------------------------------------------------------------
hwMathStatus PolyDer(const hwMatrix& A,
                     const hwMatrix& B,
                     hwMatrix&       D)
{
	hwMathStatus status;
    hwMatrix P;

    status = P.ConvLin(A, B);   // polynomial multiplication

	if (!status.IsOk())
    {
        status.ResetArgs();
		return status;
    }

    status = PolyDer(P, D);

    return status;
}
//------------------------------------------------------------------------------
// Computes the derivative of a quotient and returns status
//------------------------------------------------------------------------------
hwMathStatus PolyDer(const hwMatrix& A,
                     const hwMatrix& B,
                     hwMatrix&       P,
                     hwMatrix&       Q)
{
	hwMathStatus status;
    hwMatrix dAdx;
    hwMatrix dBdx;
    hwMatrix temp;

    status = PolyDer(A, dAdx);

	if (!status.IsOk())
    {
        if (status.GetArg1() != 1)
            status.ResetArgs();

		return status;
    }

    status = PolyDer(B, dBdx);

	if (!status.IsOk())
    {
        if (status.GetArg1() == 1)
            status.SetArg1(2);
        else
            status.ResetArgs();

		return status;
    }

    // numerator
    status = P.ConvLin(dAdx, B);   // polynomial multiplication

	if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(3);

		return status;
    }

    status = temp.ConvLin(A, dBdx);   // polynomial multiplication

	if (!status.IsOk())
    {
        status.ResetArgs();
		return status;
    }

    if (P.N() != 1)
        P.Transpose();

    if (temp.N() != 1)
        temp.Transpose();

    P -= temp;

    // denominator
    status = Q.ConvLin(B, B);   // polynomial multiplication

	if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(4);
    }

    if (Q.N() != 1)
        Q.Transpose();

    if (Q.IsReal())
    {
        P /= Q(0);
        Q /= Q(0);
    }
    else
    {
        P /= Q.z(0);
        Q /= Q.z(0);
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes the integral of a polynomial and returns status
//------------------------------------------------------------------------------
hwMathStatus PolyInt(const hwMatrix& D,
                     hwMatrix&       P,
                     double          k)
{
    int i;
    int size = D.Size() + 1;
	hwMathStatus status;

    if (!D.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (D.IsReal())
        status = P.Dimension(size, hwMatrix::REAL);
    else
        status = P.Dimension(size, hwMatrix::COMPLEX);

	if (!status.IsOk())
	{
		status.SetArg1(2);
		return status;
	}

    if (D.IsReal())
    {
        for (i = 0; i < size-1; ++i)
            P(i) = D(i) / (double) (size-1-i);

        P(size-1) = k;
    }
    else
    {
        for (i = 0; i < size-1; ++i)
            P.z(i) = D.z(i) / (double) (size-1-i);

        P.z(size-1) = k;
    }

    return status;
}
//------------------------------------------------------------------------------
// Performs a bisection search on a reverse sorted data array and returns status
//------------------------------------------------------------------------------
int BinarySearchR(const double* data,
                  int           npts,
                  double        find)
{
    if (!data)
    {
        return 0;
    }

    if (find > data[0])
    {
        return -1;
    }

    int num_points = npts;
    int temp_idx = num_points / 2;
    int max_idx = num_points - 1;
    int min_idx = 0;

    if (find < data[max_idx])
    {
        return max_idx;
    }

    for (int cnt = 0; cnt < num_points; cnt++)  // allow for bail-out
    {
        if (data[temp_idx] == find)
        {
            //if (temp_idx + 1 <= max_idx && data[temp_idx + 1] == find)
            //{
            //    return temp_idx + 1;
            //}
            if (temp_idx - 1 >= min_idx && data[temp_idx - 1] == find)
            {
                return temp_idx - 1;
            }
            else
            {
                return temp_idx;
            }
        }
        else if (data[temp_idx] > find)
        {
            if (temp_idx == max_idx)
            {
                return temp_idx;
            }

            if (data[temp_idx + 1] < find)
            {
                return temp_idx;
            }
            else if (data[temp_idx + 1] == find)
            {
                return temp_idx + 1;
            }

            min_idx = temp_idx;
        }
        else
        {
            if (temp_idx == min_idx)
            {
                return temp_idx;
            }

            if (data[temp_idx - 1] > find)
            {
                return temp_idx - 1;
            }

            max_idx = temp_idx;
        }

        temp_idx = (max_idx + min_idx) / 2;
    }

    return 0; // bail-out case
}
//------------------------------------------------------------------------------
// Performs linear interpolation and returns status
//------------------------------------------------------------------------------
hwMathStatus LinearInterp(const hwMatrix& x_old,
                          const hwMatrix& y_old,
                          const hwMatrix& x_new,
                          hwMatrix&       y_new,
                          bool            requireUniqueX,
                          int             extrap,
                          double          extrapVal)
{
    hwMathStatus status;

    if (!x_old.IsReal())
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 1);

    if (!x_old.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!y_old.IsReal())
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 2);

    if (!y_old.IsVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    if (!x_new.IsReal())
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 3);

    if (!x_new.IsEmptyOrVector())
        return status(HW_MATH_ERR_VECTOR, 3);

    if (x_old.Size() < 2)
        return status(HW_MATH_ERR_TOOFEWPOINTS, 1);

    if (x_old.Size() != y_old.Size())
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    status = y_new.Dimension(x_new.M(), x_new.N(), hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(4);
        else
            status.ResetArgs();

        return status;
    }

    int i;
    int n = x_old.Size();
    int nn = x_new.Size();
    long idx;
    bool ascendingX;

    if (x_old(0) < x_old(n - 1))
        ascendingX = true;
    else
        ascendingX = false;

    // check for multiple x_old values
    for (i = 0; i < n - 2; i++)
    {
        if (requireUniqueX && x_old(i) == x_old(i + 1))
            return status(HW_MATH_ERR_NONUNIQUE, 1);

        if (x_old(i) == x_old(i+2))
            return status(HW_MATH_ERR_CONSECUTIVE3, 1);
    }

    if (requireUniqueX && n > 1 && x_old(n - 2) == x_old(n - 1))
        return status(HW_MATH_ERR_NONUNIQUE, 1);

    for (i = 0; i < nn; i++)
    {
        if (ascendingX)
            idx = BinarySearch(x_old.GetRealData(), n, x_new(i));
        else
            idx = BinarySearchR(x_old.GetRealData(), n, x_new(i));

        if (extrap == 1)
        {
            // extrapolate
            if (idx < 0)
                idx = 0;
            else if (idx >= n - 1)
                idx = n - 2;

            if (x_old(idx) == x_old(idx + 1))
            {
                if (requireUniqueX)
                    return status(HW_MATH_ERR_NONUNIQUE, 1);

                y_new(i) = y_old(idx + 1);
            }
            else
            {
                y_new(i) = (y_old(idx + 1) - y_old(idx)) / (x_old(idx + 1) - x_old(idx)) *
                           (x_new(i) - x_old(idx)) + y_old(idx);
            }
        }
        else if (extrap == 0)
        {
            // no extrapolation, assign fill value
            if (idx < 0)
            {
                y_new(i) = extrapVal;
            }
            else if (idx == n - 1)
            {
                if (fabs(x_new(i) - x_old(idx)) < 1.0e-10)
                    y_new(i) = y_old(idx);
                else
                    y_new(i) = extrapVal;
            }
            else if (x_old(idx) == x_old(idx + 1))
            {
                y_new(i) = y_old(idx + 1);
            }
            else
            {
                y_new(i) = (y_old(idx + 1) - y_old(idx)) / (x_old(idx + 1) - x_old(idx)) *
                           (x_new(i) - x_old(idx)) + y_old(idx);
            }
        }
        else
        {
            // no extrapolation, return error
            if (idx < 0)
            {
                return status(HW_MATH_ERR_BADRANGE, 3);
            }
            else if (idx == n - 1)
            {
                if (fabs(x_new(i) - x_old(idx)) < 1.0e-10)
                    y_new(i) = y_old(idx);
                else
                    return status(HW_MATH_ERR_BADRANGE, 3);
            }
            else if (x_old(idx) == x_old(idx + 1))
            {
                y_new(i) = y_old(idx + 1);
            }
            else
            {
                y_new(i) = (y_old(idx + 1) - y_old(idx)) / (x_old(idx + 1) - x_old(idx)) *
                           (x_new(i) - x_old(idx)) + y_old(idx);
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes piecewise linear interpolation coefficients and returns status
//------------------------------------------------------------------------------
hwMathStatus LinearInterp(const hwMatrix& x_old,
                          const hwMatrix& y_old,
                          hwMatrix&       coefs)
{
    hwMathStatus status;

    if (!x_old.IsReal())
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 1);

    if (!x_old.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!y_old.IsReal())
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 2);

    if (!y_old.IsVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    if (x_old.Size() < 2)
        return status(HW_MATH_ERR_TOOFEWPOINTS, 1);

    if (x_old.Size() != y_old.Size())
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    // compute coefficients
    int n = x_old.Size();
    status = coefs.Dimension(n - 1, 2, hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(3);
        else
            status.ResetArgs();

        return status;
    }

    for (int i = 0; i < n - 1; ++i)
    {
        coefs(i, 0) = (y_old(i + 1) - y_old(i)) / (x_old(i + 1) - x_old(i));
        coefs(i, 1) = y_old(i);
    }

    return status;
}
//------------------------------------------------------------------------------
// Performs piecewise cubic hermite interpolation and returns status
//------------------------------------------------------------------------------
hwMathStatus PchipInterp(const hwMatrix& x_old,
                         const hwMatrix& y_old,
                         const hwMatrix& x_new,
                         hwMatrix&       y_new,
                         int             extrap,
                         double          extrapVal)
{
    hwMathStatus status;

    if (!x_old.IsReal())
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 1);

    if (!x_old.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!y_old.IsReal())
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 2);

    if (!y_old.IsVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    if (!x_new.IsReal())
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 3);

    if (!x_new.IsEmptyOrVector())
        return status(HW_MATH_ERR_VECTOR, 3);

    if (x_old.Size() < 2)
        return status(HW_MATH_ERR_TOOFEWPOINTS, 1);

    if (x_old.Size() != y_old.Size())
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    status = y_new.Dimension(x_new.M(), x_new.N(), hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(4);
        else
            status.ResetArgs();

        return status;
    }

    // compute derivatives
    int n = x_old.Size();

    hwMatrix d(n, hwMatrix::REAL);
    hwMatrix h(n-1, hwMatrix::REAL);
    hwMatrix delta(n-1, hwMatrix::REAL);

    d.SetElements(0.0);

    for (int i = 0; i < n-1; ++i)
    {
        h(i) = x_old(i+1) - x_old(i);

        if (h(i) == 0.0)
            return status(HW_MATH_ERR_NONUNIQUE, 1);

        delta(i) = (y_old(i+1) - y_old(i)) / h(i);

        if (i > 0 && (delta(i)*delta(i-1) > 0.0))
        {
            double w1 = 2.0 * h(i) + h(i-1);
            double w2 = h(i) + 2.0 * h(i-1);
            d(i) = (w1+w2) / (w1/delta(i-1) + w2/delta(i));
        }
    }

    d(0) = ((2.0*h(0)+h(1))*delta(0) - h(0)*delta(1)) / (h(0)+h(1));

    if (d(0) * delta(0) < 0.0)
        d(0) = 0.0;
    else if (delta(0) * delta(1) < 0.0 && fabs(d(0)) > fabs(3.0*delta(0)))
        d(0) = 3.0 * delta(0);

    d(n-1) = ((2.0*h(n-2)+h(n-3))*delta(n-2) - h(n-2)*delta(n-3)) / (h(n-2)+h(n-3));

    if (d(n-1) * delta(n-2) < 0.0)
        d(n-1) = 0.0;
    else if (delta(n-2) * delta(n-3) < 0.0 && fabs(d(n-1)) > fabs(3.0*delta(n-2)))
        d(n-1) = 3.0 * delta(n-2);

    // interpolate
    int nn = x_new.Size();
    long idx;
    bool ascendingX;

    if (x_old(0) < x_old(n - 1))
        ascendingX = true;
    else
        ascendingX = false;    // the sign changes take care of themselves

    for (int i = 0; i < nn; i++)
    {
        if (ascendingX)
            idx = BinarySearch(x_old.GetRealData(), n, x_new(i));
        else
            idx = BinarySearchR(x_old.GetRealData(), n, x_new(i));

        if (extrap == 1)
        {
            // extrapolate
            if (idx < 0)
            {
                idx = 0;
            }
            else if (idx >= n - 1)
            {
                idx = n - 2;
            }

            double c = (3.0 * delta(idx) - 2.0 * d(idx) - d(idx + 1)) / h(idx);
            double b = (d(idx) - 2.0 * delta(idx) + d(idx + 1)) / (h(idx) * h(idx));
            double s = x_new(i) - x_old(idx);

            y_new(i) = y_old(idx) + s * (d(idx) + s * (c + s * b));
        }
        else if (extrap == 0)
        {
            // no extrapolation, assign fill value
            if (idx < 0)
            {
                y_new(i) = extrapVal;
            }
            else if (idx == n - 1)
            {
                if (fabs(x_new(i) - x_old(idx)) < 1.0e-10)
                {
                    y_new(i) = y_old(idx);
                }
                else
                {
                    y_new(i) = extrapVal;
                }
            }
            else
            {
                double c = (3.0 * delta(idx) - 2.0 * d(idx) - d(idx + 1)) / h(idx);
                double b = (d(idx) - 2.0 * delta(idx) + d(idx + 1)) / (h(idx) * h(idx));
                double s = x_new(i) - x_old(idx);

                y_new(i) = y_old(idx) + s * (d(idx) + s * (c + s * b));
            }
        }
        else
        {
            // no extrapolation, return error
            if (idx < 0)
            {
                return status(HW_MATH_ERR_BADRANGE, 3);
            }
            else if (idx == n - 1)
            {
                if (fabs(x_new(i) - x_old(idx)) < 1.0e-10)
                {
                    y_new(i) = y_old(idx);
                }
                else
                {
                    return status(HW_MATH_ERR_BADRANGE, 3);
                }
            }
            else
            {
                double c = (3.0 * delta(idx) - 2.0 * d(idx) - d(idx + 1)) / h(idx);
                double b = (d(idx) - 2.0 * delta(idx) + d(idx + 1)) / (h(idx) * h(idx));
                double s = x_new(i) - x_old(idx);

                y_new(i) = y_old(idx) + s * (d(idx) + s * (c + s * b));
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes piecewise cubic hermite coefficients and returns status
//------------------------------------------------------------------------------
hwMathStatus PchipInterp(const hwMatrix& x_old,
                         const hwMatrix& y_old,
                         hwMatrix&       coefs)
{
    hwMathStatus status;

    if (!x_old.IsReal())
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 1);

    if (!x_old.IsVector())
        return status(HW_MATH_ERR_VECTOR, 1);

    if (!y_old.IsReal())
        return status(HW_MATH_ERR_COMPLEXSUPPORT, 2);

    if (!y_old.IsVector())
        return status(HW_MATH_ERR_VECTOR, 2);

    if (x_old.Size() < 2)
        return status(HW_MATH_ERR_TOOFEWPOINTS, 1);

    if (x_old.Size() != y_old.Size())
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

    // compute derivatives
    int n = x_old.Size();

    hwMatrix d(n, hwMatrix::REAL);
    hwMatrix h(n - 1, hwMatrix::REAL);
    hwMatrix delta(n - 1, hwMatrix::REAL);

    d.SetElements(0.0);

    for (int i = 0; i < n - 1; ++i)
    {
        h(i) = x_old(i + 1) - x_old(i);

        if (h(i) == 0.0)
            return status(HW_MATH_ERR_NONUNIQUE, 1);

        delta(i) = (y_old(i + 1) - y_old(i)) / h(i);

        if (i > 0 && (delta(i) * delta(i - 1) > 0.0))
        {
            double w1 = 2.0 * h(i) + h(i - 1);
            double w2 = h(i) + 2.0 * h(i - 1);
            d(i) = (w1 + w2) / (w1 / delta(i - 1) + w2 / delta(i));
        }
    }

    d(0) = ((2.0 * h(0) + h(1)) * delta(0) - h(0) * delta(1)) / (h(0) + h(1));

    if (d(0) * delta(0) < 0.0)
        d(0) = 0.0;
    else if (delta(0) * delta(1) < 0.0 && fabs(d(0)) > fabs(3.0 * delta(0)))
        d(0) = 3.0 * delta(0);

    d(n - 1) = ((2.0 * h(n - 2) + h(n - 3)) * delta(n - 2) - h(n - 2) * delta(n - 3)) / (h(n - 2) + h(n - 3));

    if (d(n - 1) * delta(n - 2) < 0.0)
        d(n - 1) = 0.0;
    else if (delta(n - 2) * delta(n - 3) < 0.0 && fabs(d(n - 1)) > fabs(3.0 * delta(n - 2)))
        d(n - 1) = 3.0 * delta(n - 2);

    // compute coefficients
    status = coefs.Dimension(n - 1, 4, hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(3);
        else
            status.ResetArgs();

        return status;
    }

    for (int i = 0; i < n - 1; ++i)
    {
        double c = (3.0 * delta(i) - 2.0 * d(i) - d(i+1)) / h(i);
        double b = (d(i) - 2.0 * delta(i) + d(i + 1)) / (h(i) * h(i));
        coefs(i, 0) = b;
        coefs(i, 1) = c;
        coefs(i, 2) = d(i);
        coefs(i, 3) = y_old(i);
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes knot-a-not cubic spline second derivatives and returns status
//------------------------------------------------------------------------------
static hwMathStatus SplineDerivatives2(const hwMatrix& x_old, 
                                       const hwMatrix& y_old,
	                                   hwMatrix&       deriv_2)
{
	if (!x_old.IsReal())
		return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);

	if (!x_old.IsVector())
		return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

	if (!y_old.IsReal())
		return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);

	if (!y_old.IsVector())
		return hwMathStatus(HW_MATH_ERR_VECTOR, 2);

	int n = x_old.Size();

	if (n < 4)
		return hwMathStatus(HW_MATH_ERR_TOOFEWPOINTS, 1);

	if (y_old.Size() != n)
		return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);

	hwMathStatus status;

	status = deriv_2.Dimension(n, 1, hwMatrix::REAL);

	if (!status.IsOk())
	{
		if (status.GetArg1() == 0)
			status.SetArg1(3);
		else
			status.ResetArgs();

		return status;
	}

	// calculate 2nd derivatives
	int i;
	hwMatrix a(n, hwMatrix::REAL);
	hwMatrix b(n, hwMatrix::REAL);
	hwMatrix c(n, hwMatrix::REAL);
	hwMatrix d(n, hwMatrix::REAL);
/*
    // natural spline (based on book from Yogesh Jaluria)
    //
    // Note: This code is disabled in favor of the not-a-knot
    // spline and is left here for historical value. The natural
    // spline code uses 1 based indexing. If the natural spline
    // is re-enabled it would be better to commonize with the
    // not-a-knot code. See the comments in the not-a-knot code
    // below for details.

    for (i = 1; i < n-1; i++)
	{
		a(i) = x_old(i) - x_old(i-1);

		if (a(i) == 0.0)
			return status(HW_MATH_ERR_NONUNIQUE, 1);

		b(i) = 2.0 * (x_old(i+1) - x_old(i-1));
		c(i) = x_old(i+1) - x_old(i);
		d(i) = 6.0 * ((y_old(i+1) - y_old(i)) / c(i) -
			          (y_old(i) - y_old(i-1)) / a(i));
	}

    // solve tridiagonal system and compute 2nd derivatives
	for (i = 2; i < n-1; i++)
	{
		b(i) -= a(i) * c(i-1) / b(i-1);
		d(i) -= a(i) * d(i-1) / b(i-1);
	}

	deriv_2(0) = 0.0;
	deriv_2(n-1) = 0.0;
	deriv_2(n-2) = d(n-2) / b(n-2);

	for (i = 1; i < n-2; i++)
	{
		idx = n - i - 2;
		deriv_2(idx) = (d(idx) - c(idx) * deriv_2(idx+1)) / b(idx);
	}
    // end natural spline
*/
    // populate spline matrix
	hwMatrix h(n, hwMatrix::REAL);

    for (i = 0; i < n-1; i++)
    {
        h(i) = x_old(i+1) - x_old(i);

		if (h(i) == 0.0)
			return status(HW_MATH_ERR_NONUNIQUE, 1);
    }

    for (i = 0; i < n-3; i++)
    {
        a(i) = h(i+1);
        b(i) = 2.0 * (h(i) + h(i+1));
        c(i) = h(i+1);
        d(i) = 6.0 * ((y_old(i+2) - y_old(i+1)) / h(i+1) -
                      (y_old(i+1) - y_old(i)) / h(i));
    }

    b(n-3) = 2.0 * (h(n-3) + h(n-2));
    d(n-3) = 6.0 * ((y_old(n-1) - y_old(n-2)) / h(n-2) -
                    (y_old(n-2) - y_old(n-3)) / h(n-3));

    // set not-a-knot end conditions
    // disable these lines for the natural spline
    a(n-4) -= (h(n-2)*h(n-2) / h(n-3));
    b(0) += (h(0) + h(0)*h(0) / h(1));
    b(n-3) += (h(n-2) + h(n-2)*h(n-2) / h(n-3));
    c(0) -= (h(0)*h(0) / h(1));

    // solve tridiagonal system using Thomas algorithm adjusted
    // for system of n-2 equations
    d(0) /= b(0);

    for (i = 0; i < n-3; ++i)
    {
        c(i) /= b(i);
        b(i+1) -= a(i)*c(i);
    }

    for (i = 1; i < n-2; ++i)
        d(i) = (d(i) - a(i-1) * d(i-1)) / b(i);

    for (i = n-4; i >= 0; --i)
        d(i) -= d(i+1) * c(i);

    // compute 2nd derivatives
    for (i = n-3; i >= 0; --i)
        deriv_2(i+1) = d(i);

    // set not-a-knot conditions
    // set these values equal to 0.0 for the natural spline
    deriv_2(0) = (1.0 + h(0) / h(1)) * deriv_2(1) - h(0) / h(1) * deriv_2(2);
    deriv_2(n-1) = (1.0 + h(n-2) / h(n-3)) * deriv_2(n-2) - h(n-2) / h(n-3) * deriv_2(n-3);

	return status;
}
//------------------------------------------------------------------------------
// Performs knot-a-not cubic spline interpolation and returns status
//------------------------------------------------------------------------------
hwMathStatus Spline(const hwMatrix& x_old,
                    const hwMatrix& y_old,
	                const hwMatrix& x_new,
                    hwMatrix&       y_new,
                    int             extrap,
                    double          extrapVal)
{
    hwMathStatus status;
	hwMatrix deriv_2;
    
    // get second derivatives
    status = SplineDerivatives2(x_old, y_old, deriv_2);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 4)
            status.ResetArgs();

        return status;
    }

	if (!x_new.IsReal())
		return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);

	if (!x_new.IsEmptyOrVector())
		return hwMathStatus(HW_MATH_ERR_VECTOR, 3);

	status = y_new.Dimension(x_new.M(), x_new.N(), hwMatrix::REAL);

	if (!status.IsOk())
	{
		if (status.GetArg1() == 0)
			status.SetArg1(4);
		else
			status.ResetArgs();

		return status;
	}

    // compute interpolated/extrapolated points
	int i;
	int n = x_old.Size();
	int nn = x_new.Size();
    long idx;
	const double* x_start = x_old.GetRealData();
    bool ascendingX;

    if (x_old(0) < x_old(n - 1))
        ascendingX = true;
    else
        ascendingX = false;    // the sign changes take care of themselves

    for (i = 0; i < nn; i++)
    {
        if (ascendingX)
            idx = BinarySearch(x_start, n, x_new(i));
        else
            idx = BinarySearchR(x_start, n, x_new(i));

        if (extrap == 1)
        {
            // extrapolate
            if (idx < 0)
            {
                idx = 0;
            }
            else if (idx >= n - 1)
            {
                idx = n - 2;
            }

            double s1 = x_old(idx + 1) - x_old(idx);
            double s2 = x_new(i) - x_old(idx);
            double s3 = x_old(idx + 1) - x_new(i);

            y_new(i) = deriv_2(idx) * s3 * (s3 * s3 / s1 - s1) / 6.0
                     + deriv_2(idx + 1) * s2 * (s2 * s2 / s1 - s1) / 6.0
                     + y_old(idx) * s3 / s1 + y_old(idx + 1) * s2 / s1;
        }
        else if (extrap == 0)
        {
            // no extrapolation, assign fill value
            if (idx < 0)
            {
                y_new(i) = extrapVal;
            }
            else if (idx == n - 1)
            {
                if (fabs(x_new(i) - x_old(idx)) < 1.0e-10)
                {
                    y_new(i) = y_old(idx);
                }
                else
                {
                    y_new(i) = extrapVal;
                }
            }
            else
            {
                double s1 = x_old(idx + 1) - x_old(idx);
                double s2 = x_new(i) - x_old(idx);
                double s3 = x_old(idx + 1) - x_new(i);

                y_new(i) = deriv_2(idx) * s3 * (s3 * s3 / s1 - s1) / 6.0
                    + deriv_2(idx + 1) * s2 * (s2 * s2 / s1 - s1) / 6.0
                    + y_old(idx) * s3 / s1 + y_old(idx + 1) * s2 / s1;
            }
        }
        else
        {
            // no extrapolation, return error
            if (idx < 0)
            {
                return status(HW_MATH_ERR_BADRANGE, 3);
            }
            else if (idx == n - 1)
            {
                if (fabs(x_new(i) - x_old(idx)) < 1.0e-10)
                {
                    y_new(i) = y_old(idx);
                }
                else
                {
                    return status(HW_MATH_ERR_BADRANGE, 3);
                }
            }
            else
            {
                double s1 = x_old(idx + 1) - x_old(idx);
                double s2 = x_new(i) - x_old(idx);
                double s3 = x_old(idx + 1) - x_new(i);

                y_new(i) = deriv_2(idx) * s3 * (s3 * s3 / s1 - s1) / 6.0
                         + deriv_2(idx + 1) * s2 * (s2 * s2 / s1 - s1) / 6.0
                         + y_old(idx) * s3 / s1 + y_old(idx + 1) * s2 / s1;
            }
        }
    }

	return status;
}
//------------------------------------------------------------------------------
// Computes knot-a-not cubic spline coefficients and returns status
//------------------------------------------------------------------------------
hwMathStatus Spline(const hwMatrix& x_old,
                    const hwMatrix& y_old,
                    hwMatrix&       coefs)
{
    hwMathStatus status;
	int n = x_old.Size();
    double h;
	hwMatrix deriv_2;
    
    // get second derivatives
    status = SplineDerivatives2(x_old, y_old, deriv_2);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 4)
            status.ResetArgs();

        return status;
    }

    // compute coefficients
	status = coefs.Dimension(n-1, 4, hwMatrix::REAL);

	if (!status.IsOk())
	{
		if (status.GetArg1() == 0)
			status.SetArg1(3);
		else
			status.ResetArgs();

		return status;
	}

    for (int i = 0; i < n-1; ++i)
    {
        h = x_old(i+1) - x_old(i);
        coefs(i, 0) = (deriv_2(i+1)-deriv_2(i)) / (6.0*h);
        coefs(i, 1) = deriv_2(i) / 2.0;
        coefs(i, 2) = (y_old(i+1) - y_old(i))/h - (deriv_2(i+1) + 2.0*deriv_2(i))*(h/6.0);
        coefs(i, 3) = y_old(i);
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes clamped cubic spline second derivatives and returns status
//------------------------------------------------------------------------------
static
hwMathStatus SplineDerivatives2(const hwMatrix& x_old,
                                const hwMatrix& y_old,
	                            double          fp1,
                                double          fp2,
                                hwMatrix&       deriv_2)
{
	if (!x_old.IsReal())
		return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);

	if (!x_old.IsVector())
		return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

	if (!y_old.IsReal())
		return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);

	if (!y_old.IsVector())
		return hwMathStatus(HW_MATH_ERR_VECTOR, 2);

	int n = x_old.Size();

	if (n < 4)
		return hwMathStatus(HW_MATH_ERR_TOOFEWPOINTS, 1);

	if (y_old.Size() != n)
		return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);

	hwMathStatus status;

	status = deriv_2.Dimension(n, 1, hwMatrix::REAL);

	if (!status.IsOk())
	{
		if (status.GetArg1() == 0)
			status.SetArg1(5);
		else
			status.ResetArgs();

		return status;
	}

	// calculate 2nd derivatives
	int i;
	hwMatrix a(n-1, hwMatrix::REAL);   // diagonal
	hwMatrix b(n, hwMatrix::REAL); // super-diagonal
	hwMatrix c(n-1, hwMatrix::REAL);   // sub-diagonal
	hwMatrix& d = deriv_2;

    // populate spline matrix
	hwMatrix h(n-1, hwMatrix::REAL);

    for (i = 0; i < n-1; i++)
    {
        h(i) = x_old(i+1) - x_old(i);

		if (h(i) <= 0.0)
			return status(HW_MATH_ERR_NONINCREASE, 1);
    }

    a(0) = h(0);
    b(0) = 2.0 * h(0);
    d(0) = 6.0 * ((y_old(1) - y_old(0)) / h(0) - fp1);

    for (i = 1; i < n-1; i++)
    {
        a(i) = h(i);
        b(i) = 2.0 * (h(i-1) + h(i));
        c(i-1) = h(i-1);
        d(i) = 6.0 * ((y_old(i+1) - y_old(i)) / h(i) -
                      (y_old(i) - y_old(i-1)) / h(i-1));
    }

    b(n-1) = 2.0 * h(n-2);
    c(n-2) = h(n-2);
    d(n-1) = 6.0 * (fp2 - (y_old(n-1) - y_old(n-2)) / h(n-2));

    // solve tridiagonal system using Thomas algorithm
    d(0) /= b(0);

    for (i = 0; i < n-1; ++i)
    {
        c(i) /= b(i);
        b(i+1) -= a(i)*c(i);
    }

    for (i = 1; i < n; ++i)
        d(i) = (d(i) - a(i-1) * d(i-1)) / b(i);

    for (i = n-2; i >= 0; --i)
        d(i) -= d(i+1) * c(i);

    return status;
}
//------------------------------------------------------------------------------
// Performs clamped cubic spline interpolation and returns status
//------------------------------------------------------------------------------
hwMathStatus Spline(const hwMatrix& x_old,
                    const hwMatrix& y_old,
                    double          fp1,
                    double          fp2,
                    const hwMatrix& x_new,
                    hwMatrix&       y_new,
                    int             extrap,
                    double          extrapVal)
{
    hwMathStatus status;
    hwMatrix deriv_2;

    // get second derivatives
    status = SplineDerivatives2(x_old, y_old, fp1, fp2, deriv_2);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 5)
            status.ResetArgs();

        return status;
    }

    if (!x_new.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);

    if (!x_new.IsEmptyOrVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 3);

    status = y_new.Dimension(x_new.M(), x_new.N(), hwMatrix::REAL);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 0)
            status.SetArg1(6);
        else
            status.ResetArgs();

        return status;
    }

    // compute interpolated/extrapolated points
    int i;
    int idx;
    int n = x_old.Size();
    int nn = x_new.Size();
    const double* x_start = x_old.GetRealData();

    for (i = 0; i < nn; i++)
    {
        idx = BinarySearch(x_start, n, x_new(i));

        if (extrap == 1)
        {
            // extrapolate
            if (idx < 0)
            {
                idx = 0;
            }
            else if (idx >= n - 1)
            {
                idx = n - 2;
            }

            double s1 = x_old(idx + 1) - x_old(idx);
            double s2 = x_new(i) - x_old(idx);
            double s3 = x_old(idx + 1) - x_new(i);

            y_new(i) = deriv_2(idx) * s3 * (s3 * s3 / s1 - s1) / 6.0
                     + deriv_2(idx + 1) * s2 * (s2 * s2 / s1 - s1) / 6.0
                     + y_old(idx) * s3 / s1 + y_old(idx + 1) * s2 / s1;
        }
        else if (extrap == 0)
        {
            // no extrapolation, assign fill value
            if (idx < 0)
            {
                y_new(i) = extrapVal;
            }
            else if (idx == n - 1)
            {
                if (fabs(x_new(i) - x_old(idx)) < 1.0e-10)
                {
                    y_new(i) = y_old(idx);
                }
                else
                {
                    y_new(i) = extrapVal;
                }
            }
            else
            {
                double s1 = x_old(idx + 1) - x_old(idx);
                double s2 = x_new(i) - x_old(idx);
                double s3 = x_old(idx + 1) - x_new(i);

                y_new(i) = deriv_2(idx) * s3 * (s3 * s3 / s1 - s1) / 6.0
                         + deriv_2(idx + 1) * s2 * (s2 * s2 / s1 - s1) / 6.0
                         + y_old(idx) * s3 / s1 + y_old(idx + 1) * s2 / s1;
            }
        }
        else
        {
            // no extrapolation, return error
            if (idx < 0)
            {
                return status(HW_MATH_ERR_BADRANGE, 3);
            }
            else if (idx == n - 1)
            {
                if (fabs(x_new(i) - x_old(idx)) < 1.0e-10)
                {
                    y_new(i) = y_old(idx);
                }
                else
                {
                    return status(HW_MATH_ERR_BADRANGE, 3);
                }
            }
            else
            {
                double s1 = x_old(idx + 1) - x_old(idx);
                double s2 = x_new(i) - x_old(idx);
                double s3 = x_old(idx + 1) - x_new(i);

                y_new(i) = deriv_2(idx) * s3 * (s3 * s3 / s1 - s1) / 6.0
                         + deriv_2(idx + 1) * s2 * (s2 * s2 / s1 - s1) / 6.0
                         + y_old(idx) * s3 / s1 + y_old(idx + 1) * s2 / s1;
            }
        }
    }

	return status;
}
//------------------------------------------------------------------------------
// Computes clamped cubic spline coefficients and returns status
//------------------------------------------------------------------------------
hwMathStatus Spline(const hwMatrix& x_old,
                    const hwMatrix& y_old,
                    double          fp1,
                    double          fp2,
                    hwMatrix&       coefs)
{
    hwMathStatus status;
	int n = x_old.Size();
    double h;
	hwMatrix deriv_2;
    
    // get second derivatives
    status = SplineDerivatives2(x_old, y_old, fp1, fp2, deriv_2);

    if (!status.IsOk())
    {
        if (status.GetArg1() == 4)
            status.ResetArgs();

        return status;
    }

    // compute coefficients
    status = coefs.Dimension(n-1, 4, hwMatrix::REAL);

	if (!status.IsOk())
	{
		if (status.GetArg1() == 0)
			status.SetArg1(5);
		else
			status.ResetArgs();

		return status;
	}

    for (int i = 0; i < n-1; ++i)
    {
        h = x_old(i+1) - x_old(i);
        coefs(i, 0) = (deriv_2(i+1)-deriv_2(i)) / (6.0*h);
        coefs(i, 1) = deriv_2(i) / 2.0;
        coefs(i, 2) = (y_old(i+1) - y_old(i))/h - (deriv_2(i+1) + 2.0*deriv_2(i))*(h/6.0);
        coefs(i, 3) = y_old(i);
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Performs bilinear interpolation on a patch
// https://en.wikipedia.org/wiki/Bilinear_interpolation
//------------------------------------------------------------------------------
static double BilinearPatch(const hwMatrix& x_old, 
                            const hwMatrix& y_old, 
                            const hwMatrix& z_old,
                            int             idxr, 
                            int             idxc, 
                            double          x_new, 
                            double          y_new)

{
    // scale coordinates to unit cube
    double x = (x_new - x_old(idxc)) / (x_old(idxc+1) - x_old(idxc));
    double y = (y_new - y_old(idxr)) / (y_old(idxr+1) - y_old(idxr));
    double a00 = z_old(idxr, idxc);
    double a10 = z_old(idxr, idxc+1) - a00;
    double a01 = z_old(idxr+1, idxc) - a00;
    double a11 = z_old(idxr+1, idxc+1) + a00 - (z_old(idxr, idxc+1) + z_old(idxr+1, idxc));

    return a00 + a10*x + a01*y + a11*x*y;
}
//------------------------------------------------------------------------------
// Performs bilinear interpolation and returns status
//------------------------------------------------------------------------------
hwMathStatus BilinearInterp(const hwMatrix& x_old,
                            const hwMatrix& y_old,
                            const hwMatrix& z_old,
                            const hwMatrix& x_new,
                            const hwMatrix& y_new,
                            hwMatrix&       z_new,
                            int             extrap,
                            double          extrapVal)
{
    hwMathStatus status;

    if (!x_old.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (!y_old.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 2);

    if (!z_old.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 3);

    if (!x_new.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 4);

    if (!y_new.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 5);

    int zm_old = z_old.M();
    int zn_old = z_old.N();
    int zm_new;
    int zn_new;

    if (x_old.M() == y_old.M() && x_old.N() == y_old.N() &&
        x_old.M() == zm_old    && x_old.N() == zn_old)
    {
        // extract single vectors from x_old and y_old, assuming that the
        // rows of x_old are identical and columns of y_old are also
        int m = y_old.M();

        if (m)
        {
            const double* y_old_col = y_old.GetRealData();
            hwMatrix y_temp(m, (void*)y_old_col, hwMatrix::REAL);
            hwMatrix x_temp;
            status = x_old.ReadRow(0, x_temp);
            return BilinearInterp(x_temp, y_temp, z_old, x_new, y_new, z_new, extrap, extrapVal);
        }
    }

    if (!x_old.IsVector() || !y_old.IsVector())
    {
        if (x_old.M() != y_old.M() || x_old.N() != y_old.N())
            return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        if (x_old.M() != z_old.M() || x_old.N() != z_old.N())
            return status(HW_MATH_ERR_ARRAYSIZE, 1, 3);
    }
    else if (x_old.Size() != z_old.N())
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 3);
    }
    else if (y_old.Size() != z_old.M())
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 2, 3);
    }

    // interpolate/extrapolate each x,y location using bilinear patch 
    int idxc;
    int idxr;
    const double* x_start = x_old.GetRealData();
    const double* y_start = y_old.GetRealData();

    if (x_new.M() == y_new.M() && x_new.N() == y_new.N())
    {
        zm_new = x_new.M();
        zn_new = x_new.N();

        status = z_new.Dimension(zm_new, zn_new, hwMatrix::REAL);

        if (!status.IsOk())
        {
            status.SetArg1(6);
            return status;
        }

        int zsize_old = zm_old * zn_old;
        int zsize_new = zm_new * zn_new;

        // find row and column indices
        for (int i = 0; i < zsize_new; ++i)
        {
            idxr = BinarySearch(y_start, zm_old, y_new(i));
            idxc = BinarySearch(x_start, zn_old, x_new(i));

            if (extrap == 1)
            {
                // extrapolate
                if (idxr < 0)
                {
                    idxr = 0;
                }
                else if (idxr >= zm_old - 1)
                {
                    idxr = zm_old - 2;
                }

                if (idxc < 0)
                {
                    idxc = 0;
                }
                else if (idxc >= zn_old - 1)
                {
                    idxc = zn_old - 2;
                }
            }
            else if (extrap == 0)
            {
                // no extrapolation, assign fill value
                if (idxr < 0)
                {
                    z_new(i) = extrapVal;
                    continue;
                }
                else if (idxr == zm_old - 1)
                {
                    if (fabs(y_new(i) - y_old(idxr)) > 1.0e-10)
                    {
                        z_new(i) = extrapVal;
                        continue;
                    }

                    --idxr;
                }

                if (idxc < 0)
                {
                    z_new(i) = extrapVal;
                    continue;
                }
                else if (idxc == zn_old - 1)
                {
                    if (fabs(x_new(i) - x_old(idxc)) > 1.0e-10)
                    {
                        z_new(i) = extrapVal;
                        continue;
                    }

                    --idxc;
                }
            }
            else
            {
                // no extrapolation, return error
                if (idxr < 0)
                {
                    return status(HW_MATH_ERR_BADRANGE, 4);
                }
                else if (idxr == zm_old - 1)
                {
                    if (fabs(y_new(i) - y_old(idxr)) > 1.0e-10)
                        return status(HW_MATH_ERR_BADRANGE, 4);

                    --idxr;
                }

                if (idxc < 0)
                {
                    return status(HW_MATH_ERR_BADRANGE, 5);
                }
                else if (idxc == zn_old - 1)
                {
                    if (fabs(x_new(i) - x_old(idxc)) > 1.0e-10)
                        return status(HW_MATH_ERR_BADRANGE, 5);

                    --idxc;
                }
            }

            z_new(i) = BilinearPatch(x_old, y_old, z_old, idxr, idxc,
                                     x_new(i), y_new(i));
        }
    }
    else if (x_new.IsVector() && y_new.IsVector())
    {
        zm_new = y_new.Size();
        zn_new = x_new.Size();

        status = z_new.Dimension(zm_new, zn_new, hwMatrix::REAL);

        if (!status.IsOk())
        {
            status.SetArg1(6);
            return status;
        }

        hwMatrixI idxr_vec(zm_new, hwMatrixI::REAL);

        // find row indices
        for (int i = 0; i < zm_new; ++i)
        {
            idxr = BinarySearch(y_start, zm_old, y_new(i));

            if (extrap == 1)
            {
                if (idxr < 0)
                {
                    idxr = 0;
                }
                else if (idxr >= zm_old - 1)
                {
                    idxr = zm_old - 2;
                }
            }
            else if (extrap != 0)
            {
                if (idxr < 0)
                {
                    return status(HW_MATH_ERR_BADRANGE, 4);
                }
                else if (idxr == zm_old - 1)
                {
                    if (fabs(y_new(i) - y_old(idxr)) > 1.0e-10)
                        return status(HW_MATH_ERR_BADRANGE, 4);

                    --idxr;
                }
            }

            idxr_vec(i) = idxr;
        }

        // find column indices
        for (int j = 0; j < zn_new; ++j)
        {
            idxc = BinarySearch(x_start, zn_old, x_new(j));

            if (extrap == 1)
            {
                // extrapolate
                if (idxc < 0)
                {
                    idxc = 0;
                }
                else if (idxc >= zn_old - 1)
                {
                    idxc = zn_old - 2;
                }

                // interpolate each patch
                for (int i = 0; i < zm_new; ++i)
                {
                    z_new(i, j) = BilinearPatch(x_old, y_old, z_old, idxr_vec(i), idxc,
                                                x_new(j), y_new(i));
                }
            }
            else if (extrap == 0)
            {
                // no extrapolation, assign fill value
                // interpolate each patch
                for (int i = 0; i < zm_new; ++i)
                {
                    if (idxr_vec(i) < 0)
                    {
                        z_new(i, j) = extrapVal;
                        continue;
                    }
                    else if (idxr_vec(i) == zm_old - 1)
                    {
                        if (fabs(y_new(i) - y_old(idxr_vec(i))) > 1.0e-10)
                        {
                            z_new(i, j) = extrapVal;
                            continue;
                        }

                        --idxr_vec(i);
                    }

                    if (idxc < 0)
                    {
                        z_new(i, j) = extrapVal;
                        continue;
                    }
                    else if (idxc == zn_old - 1)
                    {
                        if (fabs(x_new(j) - x_old(idxc)) > 1.0e-10)
                        {
                            z_new(i, j) = extrapVal;
                            continue;
                        }

                        --idxc;
                    }

                    z_new(i, j) = BilinearPatch(x_old, y_old, z_old, idxr_vec(i), idxc,
                                                x_new(j), y_new(i));
                }
            }
            else
            {
                // no extrapolation, return error
                if (idxc < 0)
                {
                    return status(HW_MATH_ERR_BADRANGE, 5);
                }
                else if (idxc == zn_old - 1)
                {
                    if (fabs(x_new(j) - x_old(idxc)) > 1.0e-10)
                        return status(HW_MATH_ERR_BADRANGE, 5);

                    --idxc;
                }

                // interpolate each patch
                for (int i = 0; i < zm_new; ++i)
                {
                    z_new(i, j) = BilinearPatch(x_old, y_old, z_old, idxr_vec(i), idxc,
                                                x_new(j), y_new(i));
                }
            }
        }
    }
    else
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 4, 5);
    }

    return status;
}
//------------------------------------------------------------------------------
// Performs bicubic interpolation on patch
// https://en.wikipedia.org/wiki/Bicubic_interpolation
//------------------------------------------------------------------------------
static
double BicubicPatch(const hwMatrix& x_old,
                    const hwMatrix& y_old,
                    const hwMatrix& z_old,
                    int             idxr,
                    int             idxc,
                    double          x_new,
                    double          y_new,
                    const hwMatrix& dzdx,
                    const hwMatrix& dzdy,
                    const hwMatrix& d2zdxdy)

{
    // scale coordinates and derivatives to unit cube
    double deltaX = x_old(idxc+1) - x_old(idxc);
    double deltaY = y_old(idxr+1) - y_old(idxr);
    double deltaXY = deltaX * deltaY;
    double x = (x_new - x_old(idxc)) / deltaX;
    double y = (y_new - y_old(idxr)) / deltaY;

    hwMatrix a(4, 4, hwMatrix::REAL);
    hwMatrix f(4, 4, hwMatrix::REAL);
    hwMatrix m1(4, 4, hwMatrix::REAL);
    hwMatrix m2(4, 4, hwMatrix::REAL);
    hwMatrix vx(1, 4, hwMatrix::REAL);
    hwMatrix vy(4, 1, hwMatrix::REAL);
    hwMatrix p(1, 1, hwMatrix::REAL);

    f(0) = z_old(idxr, idxc);         f(4) = z_old(idxr+1, idxc);          f(8) = dzdy(idxr, idxc)*deltaY;         f(12) = dzdy(idxr+1, idxc)*deltaY;
    f(1) = z_old(idxr, idxc+1);       f(5) = z_old(idxr+1, idxc+1);        f(9) = dzdy(idxr, idxc+1)*deltaY;       f(13) = dzdy(idxr+1, idxc+1)*deltaY;
    f(2) = dzdx(idxr, idxc)*deltaX;   f(6) = dzdx(idxr+1, idxc)*deltaX;    f(10) = d2zdxdy(idxr, idxc)*deltaXY;    f(14) = d2zdxdy(idxr+1, idxc)*deltaXY;
    f(3) = dzdx(idxr, idxc+1)*deltaX; f(7) = dzdx(idxr+1, idxc+1)*deltaX;  f(11) = d2zdxdy(idxr, idxc+1)*deltaXY;  f(15) = d2zdxdy(idxr+1, idxc+1)*deltaXY;

    m1(0) = 1.0;        m1(4) = 0.0;        m1(8) = 0.0;        m1(12) = 0.0;
    m1(1) = 0.0;        m1(5) = 0.0;        m1(9) = 1.0;        m1(13) = 0.0;
    m1(2) = -3.0;       m1(6) = 3.0;        m1(10) = -2.0;      m1(14) = -1.0;
    m1(3) = 2.0;        m1(7) = -2.0;       m1(11) = 1.0;       m1(15) = 1.0;

    m2(0) = 1.0;        m2(4) = 0.0;        m2(8) = -3.0;       m2(12) = 2.0;
    m2(1) = 0.0;        m2(5) = 0.0;        m2(9) = 3.0;        m2(13) = -2.0;
    m2(2) = 0.0;        m2(6) = 1.0;        m2(10) = -2.0;      m2(14) = 1.0;
    m2(3) = 0.0;        m2(7) = 0.0;        m2(11) = -1.0;      m2(15) = 1.0;

    a = m1 * f * m2;

    vx(0) = 1.0;        vx(1) = x;          vx(2) = x*x;        vx(3) = x*x*x;
    vy(0) = 1.0;        vy(1) = y;          vy(2) = y*y;        vy(3) = y*y*y;

    p = vx * a * vy;

    return p(0);
}
//------------------------------------------------------------------------------
// Performs spline surface interpolation and returns status
//------------------------------------------------------------------------------
hwMathStatus Spline2D(const hwMatrix& x_old,
                      const hwMatrix& y_old,
                      const hwMatrix& z_old,
                      const hwMatrix& x_new,
                      const hwMatrix& y_new,
                      hwMatrix&       z_new,
                      int             extrap,
                      double          extrapVal)
{
    hwMathStatus status;

    if (!x_old.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 1);

    if (!y_old.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 2);

    if (!z_old.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 3);

    if (!x_new.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 4);

    if (!y_new.IsReal())
        return status(HW_MATH_ERR_COMPLEX, 5);

    int zm_old = z_old.M();
    int zn_old = z_old.N();
    int zm_new;
    int zn_new;

    if (x_old.M() == y_old.M() && x_old.N() == y_old.N() &&
        x_old.M() == zm_old && x_old.N() == zn_old)
    {
        // extract single vectors from x_old and y_old, assuming that the
        // rows of x_old are identical and columns of y_old are also
        int m = y_old.M();

        if (m)
        {
            const double* y_old_col = y_old.GetRealData();
            hwMatrix y_temp(m, (void*)y_old_col, hwMatrix::REAL);
            hwMatrix x_temp;
            status = x_old.ReadRow(0, x_temp);
            return Spline2D(x_temp, y_temp, z_old, x_new, y_new, z_new, extrap, extrapVal);
        }
    }

    if (!x_old.IsVector() || !y_old.IsVector())
    {
        if (x_old.M() != y_old.M() || x_old.N() != y_old.N())
            return status(HW_MATH_ERR_ARRAYSIZE, 1, 2);

        if (x_old.M() != z_old.M() || x_old.N() != z_old.N())
            return status(HW_MATH_ERR_ARRAYSIZE, 1, 3);
    }
    else if (x_old.Size() != z_old.N())
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 1, 3);
    }
    else if (y_old.Size() != z_old.M())
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 2, 3);
    }

    // compute spline partial derivatives
    hwMatrix dzdx(zm_old, zn_old, hwMatrix::REAL);
    hwMatrix dzdy(zm_old, zn_old, hwMatrix::REAL);
    hwMatrix d2zdxdy(zm_old, zn_old, hwMatrix::REAL);
    hwMatrix coefs;
    hwMatrix polyDeriv(3, hwMatrix::REAL);

    if (x_new.IsVector() && y_new.IsVector())
    {
        // (y_old, x_old) contain the (row, col) domain values
        zm_new = y_new.Size();
        zn_new = x_new.Size();
    }
    else if (x_new.M() == y_new.M() && x_new.N() == y_new.N())
    {
        zm_new = x_new.M();
        zn_new = x_new.N();
    }
    else
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 4, 5);
    }

    status = z_new.Dimension(zm_new, zn_new, hwMatrix::REAL);

    if (!status.IsOk())
    {
        status.SetArg1(6);
        return status;
    }

    // compute spline 1st partial derivatives w.r.t. y for each column of z
    const double* z_old_ptr = z_old.GetRealData();

    for (int j = 0; j < zn_old; ++j)
    {
        hwMatrix z_old_col(zm_old, 1, (void*) z_old_ptr, hwMatrix::REAL);

        // get spline coefficients (stored in rows)
        status = Spline(y_old, z_old_col, coefs);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        // compute derivatives
        for (int i = 0; i < zm_old-1; ++i)
        {
            polyDeriv(2) = 3.0 * coefs(i, 0);
            polyDeriv(1) = 2.0 * coefs(i, 1);
            polyDeriv(0) = coefs(i, 2);
            PolyVal(polyDeriv, 0.0, dzdy(i, j));    // piecewise polynomials have offsets
        }

        PolyVal(polyDeriv, y_old(zm_old-1)-y_old(zm_old-2), dzdy(zm_old-1, j));
        z_old_ptr += zm_old;
    }

    // compute spline 1st partial derivatives w.r.t. x for each row of z
    hwMatrix z_old_row;

    for (int i = 0; i < zm_old; ++i)
    {
        // get spline coefficients (stored in rows)
        status = z_old.ReadRow(i, z_old_row);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        status = Spline(x_old, z_old_row, coefs);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        // compute derivatives
        for (int j = 0; j < zn_old-1; ++j)
        {
            polyDeriv(2) = 3.0 * coefs(j, 0);
            polyDeriv(1) = 2.0 * coefs(j, 1);
            polyDeriv(0) = coefs(j, 2);
            PolyVal(polyDeriv, 0.0, dzdx(i, j));    // piecewise polynomials have offsets
        }

        PolyVal(polyDeriv, x_old(zn_old-1)-x_old(zn_old-2), dzdx(i, zn_old-1));
    }

    // compute spline 2nd partial derivatives
    z_old_ptr = dzdx.GetRealData();

    for (int j = 0; j < zn_old; ++j)
    {
        hwMatrix z_old_col(zm_old, 1, (void*) z_old_ptr, hwMatrix::REAL);

        // get spline coefficients (stored in rows)
        status = Spline(y_old, z_old_col, coefs);

        if (!status.IsOk())
        {
            status.ResetArgs();
            return status;
        }

        // compute derivatives
        for (int i = 0; i < zm_old-1; ++i)
        {
            polyDeriv(2) = 3.0 * coefs(i, 0);
            polyDeriv(1) = 2.0 * coefs(i, 1);
            polyDeriv(0) = coefs(i, 2);
            PolyVal(polyDeriv, 0.0, d2zdxdy(i, j));    // piecewise polynomials have offsets
        }

        PolyVal(polyDeriv, y_old(zm_old-1)-y_old(zm_old-2), d2zdxdy(zm_old-1, j));
        z_old_ptr += zm_old;
    }

    // interpolate/extrapolate each x,y location using bicubic spline patch 
    int idxc;
    int idxr;
    const double* x_start = x_old.GetRealData();
    const double* y_start = y_old.GetRealData();

    if (x_new.M() == y_new.M() && x_new.N() == y_new.N())
    {
        zm_new = x_new.M();
        zn_new = x_new.N();

        status = z_new.Dimension(zm_new, zn_new, hwMatrix::REAL);

        if (!status.IsOk())
        {
            status.SetArg1(6);
            return status;
        }

        int zsize_old = zm_old * zn_old;
        int zsize_new = zm_new * zn_new;

        // find row and column indices
        for (int i = 0; i < zsize_new; ++i)
        {
            idxr = BinarySearch(y_start, zm_old, y_new(i));
            idxc = BinarySearch(x_start, zn_old, x_new(i));

            if (extrap == 1)
            {
                // extrapolate
                if (idxr < 0)
                {
                    idxr = 0;
                }
                else if (idxr >= zm_old - 1)
                {
                    idxr = zm_old - 2;
                }

                if (idxc < 0)
                {
                    idxc = 0;
                }
                else if (idxc >= zn_old - 1)
                {
                    idxc = zn_old - 2;
                }
            }
            else if (extrap == 0)
            {
                // no extrapolation, assign fill value
                if (idxr < 0)
                {
                    z_new(i) = extrapVal;
                    continue;
                }
                else if (idxr == zm_old - 1)
                {
                    if (fabs(y_new(i) - y_old(idxr)) > 1.0e-10)
                    {
                        z_new(i) = extrapVal;
                        continue;
                    }

                    --idxr;
                }

                if (idxc < 0)
                {
                    z_new(i) = extrapVal;
                    continue;
                }
                else if (idxc == zn_old - 1)
                {
                    if (fabs(x_new(i) - x_old(idxc)) > 1.0e-10)
                    {
                        z_new(i) = extrapVal;
                        continue;
                    }

                    --idxc;
                }
            }
            else
            {
                // no extrapolation, return error
                if (idxr < 0)
                {
                    return status(HW_MATH_ERR_BADRANGE, 4);
                }
                else if (idxr == zm_old - 1)
                {
                    if (fabs(y_new(i) - y_old(idxr)) > 1.0e-10)
                        return status(HW_MATH_ERR_BADRANGE, 4);

                    --idxr;
                }

                if (idxc < 0)
                {
                    return status(HW_MATH_ERR_BADRANGE, 5);
                }
                else if (idxc == zn_old - 1)
                {
                    if (fabs(x_new(i) - x_old(idxc)) > 1.0e-10)
                        return status(HW_MATH_ERR_BADRANGE, 5);

                    --idxc;
                }
            }

            z_new(i) = BicubicPatch(x_old, y_old, z_old, idxr, idxc,
                                    x_new(i), y_new(i), dzdx, dzdy, d2zdxdy);
        }
    }
    else if (x_new.IsVector() && y_new.IsVector())
    {
        zm_new = y_new.Size();
        zn_new = x_new.Size();

        status = z_new.Dimension(zm_new, zn_new, hwMatrix::REAL);

        if (!status.IsOk())
        {
            status.SetArg1(6);
            return status;
        }

        hwMatrixI idxr_vec(zm_new, hwMatrixI::REAL);

        // find row indices
        for (int i = 0; i < zm_new; ++i)
        {
            idxr = BinarySearch(y_start, zm_old, y_new(i));

            if (extrap == 1)
            {
                if (idxr < 0)
                {
                    idxr = 0;
                }
                else if (idxr >= zm_old - 1)
                {
                    idxr = zm_old - 2;
                }
            }
            else if (extrap != 0)
            {
                if (idxr < 0)
                {
                    return status(HW_MATH_ERR_BADRANGE, 4);
                }
                else if (idxr == zm_old - 1)
                {
                    if (fabs(y_new(i) - y_old(idxr)) > 1.0e-10)
                        return status(HW_MATH_ERR_BADRANGE, 4);

                    --idxr;
                }
            }

            idxr_vec(i) = idxr;
        }

        // find column indices
        for (int j = 0; j < zn_new; ++j)
        {
            idxc = BinarySearch(x_start, zn_old, x_new(j));

            if (extrap == 1)
            {
                // extrapolate
                if (idxc < 0)
                {
                    idxc = 0;
                }
                else if (idxc >= zn_old - 1)
                {
                    idxc = zn_old - 2;
                }

                // interpolate each patch
                for (int i = 0; i < zm_new; ++i)
                {
                    z_new(i, j) = BicubicPatch(x_old, y_old, z_old, idxr_vec(i), idxc,
                                               x_new(j), y_new(i), dzdx, dzdy, d2zdxdy);
                }
            }
            else if (extrap == 0)
            {
                // no extrapolation, assign fill value
                // interpolate each patch
                for (int i = 0; i < zm_new; ++i)
                {
                    if (idxr_vec(i) < 0)
                    {
                        z_new(i, j) = extrapVal;
                        continue;
                    }
                    else if (idxr_vec(i) == zm_old - 1)
                    {
                        if (fabs(y_new(i) - y_old(idxr_vec(i))) > 1.0e-10)
                        {
                            z_new(i, j) = extrapVal;
                            continue;
                        }

                        --idxr_vec(i);
                    }

                    if (idxc < 0)
                    {
                        z_new(i, j) = extrapVal;
                        continue;
                    }
                    else if (idxc == zn_old - 1)
                    {
                        if (fabs(x_new(j) - x_old(idxc)) > 1.0e-10)
                        {
                            z_new(i, j) = extrapVal;
                            continue;
                        }

                        --idxc;
                    }

                    z_new(i, j) = BicubicPatch(x_old, y_old, z_old, idxr_vec(i), idxc,
                                               x_new(j), y_new(i), dzdx, dzdy, d2zdxdy);
                }
            }
            else
            {
                // no extrapolation, return error
                if (idxc < 0)
                {
                    return status(HW_MATH_ERR_BADRANGE, 5);
                }
                else if (idxc == zn_old - 1)
                {
                    if (fabs(x_new(j) - x_old(idxc)) > 1.0e-10)
                        return status(HW_MATH_ERR_BADRANGE, 5);

                    --idxc;
                }

                // interpolate each patch
                for (int i = 0; i < zm_new; ++i)
                {
                    z_new(i, j) = BicubicPatch(x_old, y_old, z_old, idxr_vec(i), idxc,
                                               x_new(j), y_new(i), dzdx, dzdy, d2zdxdy);
                }
            }
        }
    }
    else
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 4, 5);
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes trilinear interpolation
// https://en.wikipedia.org/wiki/Trilinear_interpolation
//------------------------------------------------------------------------------
static void TrilinearInterp(const hwMatrix&  x_old,
                            const hwMatrix&  y_old,
                            const hwMatrix&  z_old,
                            const hwMatrixN& val_old,
                            const hwMatrix&  p_new,
                            hwMatrix&        val_new)
{
    int numPts = p_new.M();
    hwMatrix v000(numPts, hwMatrix::REAL);
    hwMatrix v100(numPts, hwMatrix::REAL);
    hwMatrix v010(numPts, hwMatrix::REAL);
    hwMatrix v001(numPts, hwMatrix::REAL);
    hwMatrix v011(numPts, hwMatrix::REAL);
    hwMatrix v101(numPts, hwMatrix::REAL);
    hwMatrix v110(numPts, hwMatrix::REAL);
    hwMatrix v111(numPts, hwMatrix::REAL);
    hwMatrix frac_x(numPts, hwMatrix::REAL);
    hwMatrix frac_y(numPts, hwMatrix::REAL);
    hwMatrix frac_z(numPts, hwMatrix::REAL);

    int    nx = x_old.Size();
    int    ny = y_old.Size();
    int    nz = z_old.Size();
    int    jx, jy, jz;
    double dx, dy, dz;
    std::vector<int> indx_old(3);
    const double* x_start = x_old.GetRealData();
    const double* y_start = y_old.GetRealData();
    const double* z_start = z_old.GetRealData();

    // compute interval sizes and fractions for unit cube scaling
    for (int i = 0; i < numPts; ++i)
    {
        jx = _max(BinarySearch(x_start, nx, p_new(i, 0)), 0);
        jy = _max(BinarySearch(y_start, ny, p_new(i, 1)), 0);
        jz = _max(BinarySearch(z_start, nz, p_new(i, 2)), 0);
        jx = _min(jx + 1, nx - 1);
        jy = _min(jy + 1, ny - 1);
        jz = _min(jz + 1, nz - 1);

        indx_old[1] = jx - 1;
        indx_old[0] = jy - 1;
        indx_old[2] = jz - 1;
        v000(i) = val_old(indx_old);
        indx_old[1] = jx;
        v100(i) = val_old(indx_old);
        indx_old[1] = jx - 1;
        indx_old[0] = jy;
        v010(i) = val_old(indx_old);
        indx_old[1] = jx;
        v110(i) = val_old(indx_old);
        indx_old[1] = jx - 1;
        indx_old[0] = jy - 1;
        indx_old[2] = jz;
        v001(i) = val_old(indx_old);
        indx_old[1] = jx;
        v101(i) = val_old(indx_old);
        indx_old[1] = jx - 1;
        indx_old[0] = jy;
        v011(i) = val_old(indx_old);
        indx_old[1] = jx;
        v111(i) = val_old(indx_old);

        dx = x_old(jx) - x_old(jx - 1);
        dy = y_old(jy) - y_old(jy - 1);
        dz = z_old(jz) - z_old(jz - 1);

        frac_x(i) = (p_new(i, 0) - x_old(jx - 1)) / dx;
        frac_y(i) = (p_new(i, 1) - y_old(jy - 1)) / dy;
        frac_z(i) = (p_new(i, 2) - z_old(jz - 1)) / dz;
    }

    // interpolation on unit cube
    // Vxyz = V000 (1 - x) (1 - y) (1 - z) +
    //        V100    x    (1 - y) (1 - z) +
    //        V010 (1 - x)    y    (1 - z) +
    //        V001 (1 - x) (1 - y)    z    +
    //        V101    x    (1 - y)    z    +
    //        V011 (1 - x)    y       z    +
    //        V110    x       y    (1 - z) +
    //        V111    x       y       z
    hwMathStatus status;
    hwMatrix temp1;
    hwMatrix temp2;
    hwMatrix temp3;
    hwMatrix temp4;
    hwMatrix frac_xc(1.0 - frac_x);
    hwMatrix frac_yc(1.0 - frac_y);
    hwMatrix frac_zc(1.0 - frac_z);

    status   = temp1.MultByElems(frac_yc, frac_zc);
    status   = temp2.MultByElems(v000, frac_xc);
    status   = temp3.MultByElems(v100, frac_x);
    status   = temp4.MultByElems(temp2 + temp3, temp1);
    val_new  = temp4;
    status   = temp1.MultByElems(frac_y, frac_zc);
    status   = temp2.MultByElems(v010, frac_xc);
    status   = temp3.MultByElems(v110, frac_x);
    status   = temp4.MultByElems(temp2 + temp3, temp1);
    val_new += temp4;
    status   = temp1.MultByElems(frac_yc, frac_z);
    status   = temp2.MultByElems(v001, frac_xc);
    status   = temp3.MultByElems(v101, frac_x);
    status   = temp4.MultByElems(temp2 + temp3, temp1);
    val_new += temp4;
    status   = temp1.MultByElems(frac_y, frac_z);
    status   = temp2.MultByElems(v011, frac_xc);
    status   = temp3.MultByElems(v111, frac_x);
    status   = temp4.MultByElems(temp2 + temp3, temp1);
    val_new += temp4;
}
//------------------------------------------------------------------------------
// Computes gradient of trilinear interpolation and returns status
// https://en.wikipedia.org/wiki/Trilinear_interpolation
//------------------------------------------------------------------------------
hwMathStatus TrilinearInterpGrad(const hwMatrix&  x_old,
                                 const hwMatrix&  y_old,
                                 const hwMatrix&  z_old,
                                 const hwMatrixN& val_old,
                                 const hwMatrix&  v_new,
                                 hwMatrix&        grad)
{
    // check inputs
    if (!x_old.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);

    if (!x_old.IsVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);

    if (!y_old.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);

    if (!y_old.IsVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);

    if (!z_old.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 3);

    if (!z_old.IsVector())
        return hwMathStatus(HW_MATH_ERR_VECTOR, 3);

    if (!val_old.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 4);

    int nx = x_old.Size();
    int ny = y_old.Size();
    int nz = z_old.Size();
    const std::vector<int> &dims_old = val_old.Dimensions();

    if (dims_old.size() != 3)
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 4);

    if (dims_old[1] != nx)
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 4);

    if (dims_old[0] != ny)
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 2, 4);

    if (dims_old[2] != nz)
        return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 3, 4);

    if (!v_new.IsReal())
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 5);

    if (v_new.N() != 3)
        return hwMathStatus(HW_MATH_ERR_ARRAYCOL3, 5);

    int numPts = v_new.M();

    // collect cube information
    // We could simply compute the gradient at each v_new(i,:) based on the 
    // grid element that contains it. However, this would be the equivalent
    // of a 3D forward difference.
    // To obtain a central difference, we construct a new grid element around
    // each v_new(i,:) of dimensions dx X dy X dz. The function value at each
    // vertex of the new grid element is interpolated. This new grid element
    // overlaps nodes of the original grid, and so captures second order
    // information for the central difference.
    // The gradient at each v_new(i,:) is computed on the new grid element.
    hwMatrix p000(numPts, 3, hwMatrix::REAL);
    hwMatrix p100(numPts, 3, hwMatrix::REAL);
    hwMatrix p010(numPts, 3, hwMatrix::REAL);
    hwMatrix p001(numPts, 3, hwMatrix::REAL);
    hwMatrix p011(numPts, 3, hwMatrix::REAL);
    hwMatrix p101(numPts, 3, hwMatrix::REAL);
    hwMatrix p110(numPts, 3, hwMatrix::REAL);
    hwMatrix p111(numPts, 3, hwMatrix::REAL);
    hwMatrix frac_x(numPts, hwMatrix::REAL);
    hwMatrix frac_y(numPts, hwMatrix::REAL);
    hwMatrix frac_z(numPts, hwMatrix::REAL);
    hwMatrix dx(numPts, hwMatrix::REAL);
    hwMatrix dy(numPts, hwMatrix::REAL);
    hwMatrix dz(numPts, hwMatrix::REAL);

    int jx, jy, jz;
    const double* x_start = x_old.GetRealData();
    const double* y_start = y_old.GetRealData();
    const double* z_start = z_old.GetRealData();

    for (int i = 0; i < numPts; ++i)
    {
        jx = _max(BinarySearch(x_start, nx, v_new(i, 0)), 0);
        jy = _max(BinarySearch(y_start, ny, v_new(i, 1)), 0);
        jz = _max(BinarySearch(z_start, nz, v_new(i, 2)), 0);
        jx = _min(jx + 1, nx - 1);
        jy = _min(jy + 1, ny - 1);
        jz = _min(jz + 1, nz - 1);

        dx(i) = (x_old(jx) - x_old(jx - 1)) / 2.0;
        dy(i) = (y_old(jy) - y_old(jy - 1)) / 2.0;
        dz(i) = (z_old(jz) - z_old(jz - 1)) / 2.0;

        p000(i, 0) = v_new(i, 0) - dx(i);
        p000(i, 1) = v_new(i, 1) - dy(i);
        p000(i, 2) = v_new(i, 2) - dz(i);

        p100(i, 0) = v_new(i, 0) + dx(i);
        p100(i, 1) = v_new(i, 1) - dy(i);
        p100(i, 2) = v_new(i, 2) - dz(i);

        p010(i, 0) = v_new(i, 0) - dx(i);
        p010(i, 1) = v_new(i, 1) + dy(i);
        p010(i, 2) = v_new(i, 2) - dz(i);

        p110(i, 0) = v_new(i, 0) + dx(i);
        p110(i, 1) = v_new(i, 1) + dy(i);
        p110(i, 2) = v_new(i, 2) - dz(i);

        p001(i, 0) = v_new(i, 0) - dx(i);
        p001(i, 1) = v_new(i, 1) - dy(i);
        p001(i, 2) = v_new(i, 2) + dz(i);

        p101(i, 0) = v_new(i, 0) + dx(i);
        p101(i, 1) = v_new(i, 1) - dy(i);
        p101(i, 2) = v_new(i, 2) + dz(i);

        p011(i, 0) = v_new(i, 0) - dx(i);
        p011(i, 1) = v_new(i, 1) + dy(i);
        p011(i, 2) = v_new(i, 2) + dz(i);

        p111(i, 0) = v_new(i, 0) + dx(i);
        p111(i, 1) = v_new(i, 1) + dy(i);
        p111(i, 2) = v_new(i, 2) + dz(i);
    }

    hwMatrix v000(numPts, hwMatrix::REAL);
    hwMatrix v100(numPts, hwMatrix::REAL);
    hwMatrix v010(numPts, hwMatrix::REAL);
    hwMatrix v001(numPts, hwMatrix::REAL);
    hwMatrix v011(numPts, hwMatrix::REAL);
    hwMatrix v101(numPts, hwMatrix::REAL);
    hwMatrix v110(numPts, hwMatrix::REAL);
    hwMatrix v111(numPts, hwMatrix::REAL);

    TrilinearInterp(x_old, y_old, z_old, val_old, p000, v000);
    TrilinearInterp(x_old, y_old, z_old, val_old, p100, v100);
    TrilinearInterp(x_old, y_old, z_old, val_old, p010, v010);
    TrilinearInterp(x_old, y_old, z_old, val_old, p110, v110);
    TrilinearInterp(x_old, y_old, z_old, val_old, p001, v001);
    TrilinearInterp(x_old, y_old, z_old, val_old, p101, v101);
    TrilinearInterp(x_old, y_old, z_old, val_old, p011, v011);
    TrilinearInterp(x_old, y_old, z_old, val_old, p111, v111);

    hwMathStatus status;
    hwMatrix temp1;

    // compute dVdx
    // note: the factors of 8 below are because dx,dy,dz were divided by 2
    // above, and having each v_new(i,:) at the center of its new grid element
    // introduces two factors of 2 into each term of the gradient of Vxyz.
    hwMatrix dVdx;
    dVdx  = v100 - v000;
    dVdx += v110 - v010;
    dVdx += v101 - v001;
    dVdx += v111 - v011;
    temp1 = dVdx;
    status = dVdx.DivideByElems(temp1, 8.0 * dx);

    // compute dVdy
    hwMatrix dVdy;
    dVdy  = v010 - v000;
    dVdy += v110 - v100;
    dVdy += v011 - v001;
    dVdy += v111 - v101;
    temp1 = dVdy;
    status = dVdy.DivideByElems(temp1, 8.0 * dy);

    // compute dVdz
    hwMatrix dVdz;
    dVdz  = v001 - v000;
    dVdz += v101 - v100;
    dVdz += v011 - v010;
    dVdz += v111 - v110;
    temp1 = dVdz;
    status = dVdz.DivideByElems(temp1, 8.0 * dz);

    // prepare output
    status = grad.Dimension(numPts, 3, hwMatrix::REAL);

    status = grad.WriteColumn(0, dVdx);
    status = grad.WriteColumn(1, dVdy);
    status = grad.WriteColumn(2, dVdz);

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Performs multilinear interpolation
//------------------------------------------------------------------------------
hwMathStatus MultilinearInterp(const hwMatrix** x_old,
                               const hwMatrixN& val_old,
                               const hwMatrix** x_new,
                               hwMatrixN&       val_new,
                               int              extrap,
                               double           extrapVal)
{
    // check inputs
    // assumption: x_old = new hwMatrix** [numDims];
    // assumption: x_new = new hwMatrix** [numDims];
    const std::vector<int>& dims = val_old.Dimensions();
    int numDims = static_cast<int> (dims.size());

    if (!val_old.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1, numDims + 1);
    }

    if (numDims == 2 && dims[1] == 1)
    {
        numDims = 1;    // vector val_old case
    }

    for (int i = 0; i < numDims; ++i)
    {
        if (x_old[i]->Size() != dims[i])
            return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, i + 1, numDims);
    }

    hwMatrixI** index_x_old = new hwMatrixI* [numDims];
    hwMatrix** coef1 = new hwMatrix* [numDims];
    hwMatrix** coef2 = new hwMatrix* [numDims];
    std::vector<int> newDims(numDims);

    // allocate work space and outputs
    for (int i = 0; i < numDims; ++i)
    {
        newDims[i] = x_new[i]->Size();
        index_x_old[i] = new hwMatrixI(newDims[i], hwMatrixI::REAL);
        coef1[i] = new hwMatrix(newDims[i], hwMatrix::REAL);
        coef2[i] = new hwMatrix(newDims[i], hwMatrix::REAL);
    }

    val_new.Dimension(newDims, hwMatrixN::REAL);

    // look up the index of each x_new[i](j) value in x_old[i]
    // and compute interval fractions
    for (int i = 0; i < numDims; ++i)
    {
        int xoldsz = dims[i];
        int xnewsz = newDims[i];

        if (!x_old[i]->IsReal())
        {
            return hwMathStatus(HW_MATH_ERR_COMPLEX, i + 1);
        }

        if (!x_new[i]->IsReal())
        {
            return hwMathStatus(HW_MATH_ERR_COMPLEX, numDims + i + 1);
        }

        const double* xold = x_old[i]->GetRealData();
        const double* xnew = x_new[i]->GetRealData();
        int* xoldIdx = index_x_old[i]->GetRealData();
        double* coef1Vals = coef1[i]->GetRealData();
        double* coef2Vals = coef2[i]->GetRealData();

        for (int j = 0; j < xnewsz; ++j)
        {
            int k = BinarySearch(xold, xoldsz, xnew[j]);

            if (k >= 0 && k < xoldsz - 1)
            {
                xoldIdx[j] = k;
                coef1Vals[j] = (xnew[j] - xold[k]) / (xold[k + 1] - xold[k]);
                coef2Vals[j] = 1.0 - coef1Vals[j];
            }
            else if (k == xoldsz - 1 && xnew[j] == xold[xoldsz - 1])
            {
                xoldIdx[j] = k - 1;
                coef1Vals[j] = 1.0;
                coef2Vals[j] = 0.0;
            }
            else if (extrap == 1)
            {
                if (k == -1)
                {
                    k = 0;
                }

                if (k == xoldsz - 1)
                {
                    --k;
                }

                xoldIdx[j] = k;
                coef1Vals[j] = (xnew[j] - xold[k]) / (xold[k + 1] - xold[k]);
                coef2Vals[j] = 1.0 - coef1Vals[j];
            }
            else
            {
                xoldIdx[j] = k;
                coef1Vals[j] = std::numeric_limits<double>::quiet_NaN();
                coef2Vals[j] = std::numeric_limits<double>::quiet_NaN();
            }
        }
    }

    // allocate a hypercube
    int numHCvertices = 1 << numDims;     // = 2 ^ numDims
    std::vector<int> HCdims(numDims);

    for (int i = 0; i < numDims; ++i)
        HCdims[i] = 2;

    hwMatrixN HCvals(HCdims, hwMatrixN::REAL);

    // interpolate each new point over its hypercube
    int numPnts = val_new.Size();
    std::vector<int> HCvertex(numDims);     // relative to val_old
    std::vector<int> HCindices(numDims);    // relative to vertex 0
    std::vector<int> val_new_indices(numDims);

    for (int i = 0; i < numPnts; ++i)
    {
        // locate hypercube vertex 0 for this point
        bool outOfRange = false;

        for (int j = 0; j < numDims; ++j)
        {
            int* xoldIdx = index_x_old[j]->GetRealData();
            int k = xoldIdx[val_new_indices[j]];

            if (k < 0)
            {
                outOfRange = true;
                k = 0;
            }
            else if (k > dims[j] - 2)
            {
                outOfRange = true;
                k = dims[j] - 2;
            }

            HCvertex[j] = k;
        }

        if (outOfRange)
        {
            if (extrap == 0)
            {
                val_new(val_new_indices) = extrapVal;
            }
            else if (extrap == -1)
            {
                // no extrapolation, return error
                return hwMathStatus(HW_MATH_ERR_BADRANGE);
            }
        }

        if (!outOfRange || extrap == 1)
        {
            // populate the hypercube
            for (int j = 0; j < numHCvertices; ++j)
            {
                HCvals(HCindices) = val_old(HCvertex);

                // advance indices to next hypercube vertex
                for (int k = 0; k < numDims; ++k)
                {
                    // increment index j if possible
                    if (HCindices[k] == 0)
                    {
                        HCindices[k] = 1;
                        ++HCvertex[k];
                        break;
                    }

                    // index j is maxed out, so reset and continue to j+1
                    HCindices[k] = 0;
                    --HCvertex[k];
                }
            }

            // interpolate the hypercube one dimension at a time
            for (int j = 0; j < numDims; ++j)
            {
                int* xoldIdx = index_x_old[j]->GetRealData();
                double* coef1Vals = coef1[j]->GetRealData();
                double* coef2Vals = coef2[j]->GetRealData();

                for (int k = 0; k < numHCvertices >> (j + 1); ++k)
                {
                    double val1 = HCvals(HCindices);
                    HCindices[j] = 1;
                    double val2 = HCvals(HCindices);
                    HCindices[j] = 0;

                    // the interpolated values are overwritten back onto
                    // location 0 of that hypercube dimension, ignoring
                    // location 1 thereafter, so that in concept, the
                    // dimensions successively collapse back onto vertex 0
                    HCvals(HCindices) = coef1Vals[val_new_indices[j]] * val2 +
                                        coef2Vals[val_new_indices[j]] * val1;

                    // advance indices to next hypercube vertex
                    for (int ii = j + 1; ii < numDims; ++ii)
                    {
                        // increment index j if possible
                        if (HCindices[ii] == 0)
                        {
                            HCindices[ii] = 1;
                            break;
                        }

                        // index ii is maxed out, so reset and continue to ii+1
                        HCindices[ii] = 0;
                    }
                }
            }

            val_new(val_new_indices) = HCvals(0);
        }

        // advance to next interpolation point
        for (int j = 0; j < numDims; ++j)
        {
            // increment index j if possible
            if (val_new_indices[j] < newDims[j] - 1)
            {
                ++val_new_indices[j];
                break;
            }

            // index j is maxed out, so reset and continue to j+1
            val_new_indices[j] = 0;
        }
    }

    // clean up work space
    for (int i = 0; i < numDims; ++i)
    {
        delete index_x_old[i];
        delete coef1[i];
        delete coef2[i];
    }

    delete [] index_x_old;
    delete [] coef1;
    delete [] coef2;

    return hwMathStatus();
}
