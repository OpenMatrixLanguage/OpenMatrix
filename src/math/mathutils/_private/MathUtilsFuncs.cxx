/**
* @file MathUtilsFuncs.cxx
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
#include "MathUtilsFuncs.h"
#include "hwMatrix.h"
#include "GeneralFuncs.h"

#include <vector>
#include <algorithm>

static const double*    data_priv_mag;
static const double*    data_priv_real;
static const hwComplex* data_priv_complex;

//------------------------------------------------------------------------------
// Returns true if element at position elem1 is lesser than that at elem2
//------------------------------------------------------------------------------
static bool AscendReal(int elem1,
                       int elem2)
{
    return data_priv_real[elem1] < data_priv_real[elem2];
}
//------------------------------------------------------------------------------
// Returns true if element at position elem1 is greater than that at elem2
//------------------------------------------------------------------------------
static bool DescendReal(int elem1,
                        int elem2)
{
    return data_priv_real[elem1] > data_priv_real[elem2];
}
//------------------------------------------------------------------------------
// Returns true if element at position elem1 is lesser than that at elem2
//------------------------------------------------------------------------------
static bool AscendComplex(int elem1,
                          int elem2)
{
    if (data_priv_mag[elem1] < data_priv_mag[elem2])
    {
        return true;
    }

    if (data_priv_mag[elem1] > data_priv_mag[elem2])
    {
        return false;
    }

    if (data_priv_complex[elem1].Arg() < data_priv_complex[elem2].Arg())
    {
        return true;
    }
    return false;
}
//------------------------------------------------------------------------------
// Returns true if element at position elem1 is greater than that at elem2
//------------------------------------------------------------------------------
static bool DescendComplex(int elem1,
                           int elem2)
{
    if (data_priv_mag[elem1] > data_priv_mag[elem2])
    {
        return true;
    }
    if (data_priv_mag[elem1] < data_priv_mag[elem2])
    {
        return false;
    }
    if (data_priv_complex[elem1].Arg() > data_priv_complex[elem2].Arg())
    {
        return true;
    }
    return false;
}
//------------------------------------------------------------------------------
// Computes the minimum value in a real matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus Min(const hwMatrix& A,
                 double*         value,
                 int*            index)
{
    if (A.IsEmpty())
    {
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);
    }
    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int indx = 0;
    int n    = A.Size();

    double min = A(0);

    for (int i = 1; i < n; i++)
    {
        if (A(i) < min)
        {
            min = A(i);
            indx = i;
        }
    }

    if (value)
    {
        (*value) = min;
    }
    if (index)
    {
        (*index) = indx;
    }
    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes the minimum value in a complex matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus Min(const hwMatrix& A,
                 hwComplex*      value,
                 int*            index)
{
    if (A.IsEmpty())
    {
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);
    }

    if (A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 1);
    }

    int indx = 0;
    int n    = A.Size();

    double mag;
    double min = A.z(0).Mag();

    for (int i = 1; i < n; i++)
    {
        mag = A.z(i).Mag();

        if (mag < min)
        {
            min = mag;
            indx = i;
        }
        else if (mag == min)
        {
            if (A.z(i).Arg() < A.z(indx).Arg())
            {
                indx = i;
            }
        }
    }

    if (value)
    {
        (*value) = A.z(indx);
    }

    if (index)
    {
        (*index) = indx;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes the minimum value from each column of a matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus Min(const hwMatrix& A,
                 hwMatrix*       mag,
                 hwMatrixI*      row)
{
    hwMathStatus status;

    int n = A.N();

    if (mag)
    {
        if (A.IsReal())
        {
            if (A.IsEmpty())
            {
                status = (A.M() == 0) ? 
                            mag->Dimension(0, n, hwMatrix::REAL) :
                            mag->Dimension(1, 0, hwMatrix::REAL);
            }
            else
            {
                status = mag->Dimension(1, n, hwMatrix::REAL);
            }
        }
        else
        {
            status = mag->Dimension(1, n, hwMatrix::COMPLEX);
        }

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
    }

    if (row)
    {
        if (A.IsEmpty())
        {
            status = (A.M() == 0) ?
                        row->Dimension(0, n, hwMatrixI::REAL):
                        row->Dimension(1, 0, hwMatrixI::REAL);
        }
        else
        {
            status = row->Dimension(1, n, hwMatrixI::REAL);
        }

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
    }

    int    index;
    double min;

    if (A.IsReal())
    {
        if (A.IsEmpty())
        {
            return status;
        }

        for (int j = 0; j < n; j++)
        {
            min   = A(0, j);
            index = 0;

            for (int i = 1; i < A.M(); i++)
            {
                if (A(i, j) < min)
                {
                    min   = A(i, j);
                    index = i;
                }
            }

            if (mag)
            {
                (*mag)(0, j) = min;
            }
            if (row)
            {
                (*row)(0, j) = index;
            }
        }
        return status;
    }

    // A is complex
    double mag_ij;

    for (int j = 0; j < n; j++)
    {
        min   = A.z(0, j).Mag();
        index = 0;

        for (int i = 1; i < A.M(); i++)
        {
            mag_ij = A.z(i, j).Mag();

            if (mag_ij < min)
            {
                min   = mag_ij;
                index = i;
            }
            else if (mag_ij == min)
            {
                if (A.z(i, j).Arg() < A.z(index, j).Arg())
                {
                    index = i;
                }
            }
        }

        if (mag)
        {
            (*mag).z(0, j) = A.z(index, j);
        }

        if (row)
        {
            (*row)(0, j) = index;
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes the minimum value from each vector of a matrix along a specified
// dimension and returns status
//------------------------------------------------------------------------------
hwMathStatus Min(const hwMatrix& A,
                 int             dim,
                 hwMatrix*       mag,
                 hwMatrixI*      index)
{
    hwMathStatus status;

    if (dim == 1)
    {
        return Min(A, mag, index);
    }
    else if (dim != 2)
    {
        return status(HW_MATH_ERR_INVALIDINPUT, 4);
    }

    // dim is 2
    int m = A.M();
    int n = A.N();

    if (mag)
    {
        if (A.IsReal())
        {
            if (A.IsEmpty())
            {
                status = (A.N() == 0) ?
                            mag->Dimension(m, 0, hwMatrix::REAL) :
                            mag->Dimension(0, 1, hwMatrix::REAL);
            }
            else
            {
                status = mag->Dimension(m, 1, hwMatrix::REAL);
            }
        }
        else
        {
            status = mag->Dimension(m, 1, hwMatrix::COMPLEX);
        }

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
    }

    if (index)
    {
        if (A.IsEmpty())
        {
            status = (A.N() == 0) ?
                        index->Dimension(m, 0, hwMatrixI::REAL):
                        index->Dimension(0, 1, hwMatrixI::REAL);
        }
        else
        {
            status = index->Dimension(m, 1, hwMatrixI::REAL);
        }

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
    }

    int    indx;
    double min;

    if (A.IsReal())
    {
        if (A.IsEmpty())
        {
            return status;
        }

        for (int i = 0; i < m; i++)
        {
            min  = A(i, 0);
            indx = 0;

            for (int j = 1; j < n; j++)
            {
                if (A(i, j) < min)
                {
                    min  = A(i, j);
                    indx = j;
                }
            }

            if (mag)
            {
                (*mag)(i) = min;
            }

            if (index)
            {
                (*index)(i) = indx;
            }
        }
    }
    else
    {
        int    indx;
        double mag_ij;

        for (int i = 0; i < m; i++)
        {
            min = A.z(i, 0).Mag();
            indx = 0;

            for (int j = 1; j < n; j++)
            {
                mag_ij = A.z(i, j).Mag();

                if (mag_ij < min)
                {
                    min  = mag_ij;
                    indx = j;
                }
                else if (mag_ij == min)
                {
                    if (A.z(i, j).Arg() < A.z(i, indx).Arg())
                    {
                        min  = mag_ij;
                        indx = j;
                    }
                }
            }

            if (mag)
            {
                (*mag).z(i) = A.z(i, indx);
            }
            if (index)
            {
                (*index)(i) = indx;
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes a matrix whose values are the minima of the arguments and
// returns status
//------------------------------------------------------------------------------
hwMathStatus Min(const hwMatrix& A,
                 double          b,
                 hwMatrix&       mag)
{
    hwMathStatus status;

    int m = A.M();
    int n = A.N();

    if (A.IsReal())
    {
        status = mag.Dimension(m, n, hwMatrix::REAL);
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

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                mag(i, j) = (A(i, j) < b) ? A(i, j) :  b;
            }
        }
    }
    else
    {
        status = mag.Dimension(m, n, hwMatrix::COMPLEX);
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

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (A.z(i, j).Mag() < b)
                {
                    mag.z(i, j) = A.z(i, j);
                }
                else if (A.z(i, j).Mag() > b)
                {
                    mag.z(i, j) = b;
                }
                else if (A.z(i, j).Imag() < 0.0)
                {
                    mag.z(i, j) = A.z(i, j);
                }
                else
                {
                    mag.z(i, j) = b;
                }
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes a matrix whose values are the minima of the arguments and
// returns status
//------------------------------------------------------------------------------
hwMathStatus Min(const hwMatrix&  A,
                 const hwComplex& b,
                 hwMatrix&        mag)
{
    int m = A.M();
    int n = A.N();

    hwMathStatus status = mag.Dimension(m, n, hwMatrix::COMPLEX);

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

    if (A.IsReal())
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (A(i, j) < b.Mag())
                {
                    mag.z(i, j) = A(i, j);
                }
                else if (A(i, j) > b.Mag())
                {
                    mag.z(i, j) = b;
                }
                else if (b.Imag() > 0.0)
                {
                    mag.z(i, j) = A(i, j);
                }
                else
                {
                    mag.z(i, j) = b;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (A.z(i, j).Mag() < b.Mag())
                {
                    mag.z(i, j) = A.z(i, j);
                }
                else if (A.z(i, j).Mag() > b.Mag())
                {
                    mag.z(i, j) = b;
                }
                else if (A.z(i, j).Arg() < b.Arg())
                {
                    mag.z(i, j) = A.z(i, j);
                }
                else
                {
                    mag.z(i, j) = b;
                }
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes a matrix whose values are the minima of the arguments and
// returns status
//------------------------------------------------------------------------------
hwMathStatus Min(const hwMatrix& A,
                 const hwMatrix& B,
                 hwMatrix&       mag)
{
    if (A.M() != B.M() || A.N() != B.N())
    {
        if (A.M() == 1 && A.N() == 1)
        {
            if (A.IsReal())
            {
                return Min(B, A(0), mag);
            }
            return Min(B, A.z(0), mag);
        }
        else if (B.M() == 1 && B.N() == 1)
        {
            if (B.IsReal())
            {
                return Min(A, B(0), mag);
            }
            return Min(A, B.z(0), mag);
        }
        else
        {
            return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
        }
    }

    int m = A.M();
    int n = A.N();

    hwMathStatus status;

    if (A.IsReal())
    {
        if (B.IsReal())
        {
            status = mag.Dimension(m, n, hwMatrix::REAL);

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

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    mag(i, j) = (A(i, j) < B(i, j)) ? A(i, j) : B(i, j);
                }
            }
        }
        else
        {
            double magB;

            status = mag.Dimension(m, n, hwMatrix::COMPLEX);
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

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    magB = B.z(i, j).Mag();

                    if (A(i, j) < magB)
                    {
                        mag.z(i, j).Set(A(i, j), 0.0);
                    }
                    else if (A(i, j) > magB)
                    {
                        mag.z(i, j) = B.z(i, j);
                    }
                    else if (B.z(i, j).Imag() > 0.0)
                    {
                        mag.z(i, j).Set(A(i, j), 0.0);
                    }
                    else
                    {
                        mag.z(i, j) = B.z(i, j);
                    }
                }
            }
        }
    }
    else
    {
        double magA;

        status = mag.Dimension(m, n, hwMatrix::COMPLEX);

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

        if (B.IsReal())
        {
            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    magA = A.z(i, j).Mag();

                    if (magA > B(i, j))
                    {
                        mag.z(i, j).Set(B(i, j), 0.0);
                    }
                    else if (magA < B(i, j))
                    {
                        mag.z(i, j) = A.z(i, j);
                    }
                    else if (A.z(i, j).Imag() > 0.0)
                    {
                        mag.z(i, j).Set(B(i, j), 0.0);
                    }
                    else
                    {
                        mag.z(i, j) = A.z(i, j);
                    }
                }
            }
        }
        else
        {
            double magB;

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    magA = A.z(i, j).Mag();
                    magB = B.z(i, j).Mag();

                    if (magA < magB)
                    {
                        mag.z(i, j) = A.z(i, j);
                    }
                    else if (magA > magB)
                    {
                        mag.z(i, j) = B.z(i, j);
                    }
                    else if (A.z(i, j).Arg() < B.z(i, j).Arg())
                    {
                        mag.z(i, j) = A.z(i, j);
                    }
                    else
                    {
                        mag.z(i, j) = B.z(i, j);
                    }
                }
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes the maximum value in a real matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus Max(const hwMatrix& A,
                 double*         value,
                 int*            index)
{
    if (A.IsEmpty())
    {
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);
    }

    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int indx = 0;
    int n    = A.Size();

    double max = A(0);

    for (int i = 1; i < n; i++)
    {
        if (A(i) > max)
        {
            max = A(i);
            indx = i;
        }
    }

    if (value)
    {
        (*value) = max;
    }

    if (index)
    {
        (*index) = indx;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes the maximum value in a complex matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus Max(const hwMatrix& A,
                 hwComplex*      value,
                 int*            index)
{
    if (A.IsEmpty())
    {
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);
    }

    if (A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 1);
    }

    int indx = 0;
    int n    = A.Size();

    double mag;
    double max = A.z(0).Mag();

    for (int i = 1; i < n; i++)
    {
        mag = A.z(i).Mag();

        if (mag > max)
        {
            max  = mag;
            indx = i;
        }
        else if (mag == max)
        {
            if (A.z(i).Arg() > A.z(indx).Arg())
            {
                indx = i;
            }
        }
    }

    if (value)
    {
        (*value) = A.z(indx);
    }

    if (index)
    {
        (*index) = indx;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes the maximum value from each column of a matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus Max(const hwMatrix& A,
                 hwMatrix*       mag,
                 hwMatrixI*      row)
{
    hwMathStatus status;

    int n = A.N();

    if (mag)
    {
        if (A.IsReal())
        {
            if (A.IsEmpty())
            {
                status = (A.M() == 0) ?
                            mag->Dimension(0, n, hwMatrix::REAL):
                            mag->Dimension(1, 0, hwMatrix::REAL);
            }
            else
            {
                status = mag->Dimension(1, n, hwMatrix::REAL);
            }
        }
        else
        {
            status = mag->Dimension(1, n, hwMatrix::COMPLEX);
        }

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
    }

    if (row)
    {
        if (A.IsEmpty())
        {
            status = (A.M() == 0) ? 
                        row->Dimension(0, n, hwMatrixI::REAL):
                        row->Dimension(1, 0, hwMatrixI::REAL);
        }
        else
        {
            status = row->Dimension(1, n, hwMatrixI::REAL);
        }
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
    }

    int    index;
    double max;

    if (A.IsReal())
    {
        if (A.IsEmpty())
        {
            return status;
        }
        for (int j = 0; j < n; j++)
        {
            max   = A(0, j);
            index = 0;

            for (int i = 1; i < A.M(); i++)
            {
                if (A(i, j) > max)
                {
                    max   = A(i, j);
                    index = i;
                }
            }

            if (mag)
            {
                (*mag)(0, j) = max;
            }
            if (row)
            {
                (*row)(0, j) = index;
            }
        }
    }
    else
    {
        double mag_ij;

        for (int j = 0; j < n; j++)
        {
            max   = A.z(0, j).Mag();
            index = 0;

            for (int i = 1; i < A.M(); i++)
            {
                mag_ij = A.z(i, j).Mag();

                if (mag_ij > max)
                {
                    max   = mag_ij;
                    index = i;
                }
                else if (mag_ij == max)
                {
                    if (A.z(i, j).Arg() > A.z(index, j).Arg())
                    {
                        index = i;
                    }
                }
            }

            if (mag)
            {
                (*mag).z(0, j) = A.z(index, j);
            }
            if (row)
            {
                (*row)(0, j) = index;
            }
        }
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes the maximum value from each vector of a matrix along a specified
// dimension and returns status
//------------------------------------------------------------------------------
hwMathStatus Max(const hwMatrix& A,
                 int             dim,
                 hwMatrix*       mag,
                 hwMatrixI*      index)
{
    if (dim == 1)
    {
        return Max(A, mag, index);
    }
    else if (dim != 2)
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 4);
    }

    // Dim is 2
    hwMathStatus status;

    int m = A.M();
    int n = A.N();

    if (mag)
    {
        if (A.IsReal())
        {
            if (A.IsEmpty())
            {
                status = (A.N() == 0) ? 
                         mag->Dimension(m, 0, hwMatrix::REAL):
                         mag->Dimension(0, 1, hwMatrix::REAL);
            }
            else
            {
                status = mag->Dimension(m, 1, hwMatrix::REAL);
            }
        }
        else
        {
            status = mag->Dimension(m, 1, hwMatrix::COMPLEX);
        }

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
    }

    if (index)
    {
        if (A.IsEmpty())
        {
            status = (A.N() == 0) ?
                     index->Dimension(m, 0, hwMatrixI::REAL) :
                     index->Dimension(0, 1, hwMatrixI::REAL);
        }
        else
        {
            status = index->Dimension(m, 1, hwMatrixI::REAL);
        }

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
    }

    int    indx;
    double max;

    if (A.IsReal())
    {
        if (A.IsEmpty())
        {
            return status;
        }

        for (int i = 0; i < m; i++)
        {
            max  = A(i, 0);
            indx = 0;

            for (int j = 1; j < n; j++)
            {
                if (A(i, j) > max)
                {
                    max  = A(i, j);
                    indx = j;
                }
            }

            if (mag)
            {
                (*mag)(i) = max;
            }

            if (index)
            {
                (*index)(i) = indx;
            }
        }
    }
    else
    {
        int    indx;
        double mag_ij;

        for (int i = 0; i < m; i++)
        {
            max  = A.z(i, 0).Mag();
            indx = 0;

            for (int j = 1; j < n; j++)
            {
                mag_ij = A.z(i, j).Mag();

                if (mag_ij > max)
                {
                    max  = mag_ij;
                    indx = j;
                }
                else if (mag_ij == max)
                {
                    if (A.z(i, j).Arg() > A.z(i, indx).Arg())
                    {
                        max  = mag_ij;
                        indx = j;
                    }
                }
            }

            if (mag)
            {
                (*mag).z(i) = A.z(i, indx);
            }

            if (index)
            {
                (*index)(i) = indx;
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes a matrix whose values are the maxima of the arguments and
// returns status
//------------------------------------------------------------------------------
hwMathStatus Max(const hwMatrix& A,
                 double          b,
                 hwMatrix&       mag)
{
    hwMathStatus status;

    int m = A.M();
    int n = A.N();

    if (A.IsReal())
    {
        status = mag.Dimension(m, n, hwMatrix::REAL);
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

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                mag(i, j) = (A(i, j) > b) ? A(i, j) : b;
            }
        }
        return status;
    }
    
    // A is complex
    status = mag.Dimension(m, n, hwMatrix::COMPLEX);

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

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (A.z(i, j).Mag() > b)
            {
                mag.z(i, j) = A.z(i, j);
            }
            else if (A.z(i, j).Mag() < b)
            {
                mag.z(i, j) = b;
            }
            else if (A.z(i, j).Imag() > 0.0)
            {
                mag.z(i, j) = A.z(i, j);
            }
            else
            {
                mag.z(i, j) = b;
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes a matrix whose values are the maxima of the arguments and
// returns status
//------------------------------------------------------------------------------
hwMathStatus Max(const hwMatrix&  A,
                 const hwComplex& b,
                 hwMatrix&        mag)
{
    int m = A.M();
    int n = A.N();

    hwMathStatus status = mag.Dimension(m, n, hwMatrix::COMPLEX);

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

    if (A.IsReal())
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (A(i, j) > b.Mag())
                {
                    mag.z(i, j) = A(i, j);
                }
                else if (A(i, j) < b.Mag())
                {
                    mag.z(i, j) = b;
                }
                else if (b.Imag() < 0.0)
                {
                    mag.z(i, j) = A(i, j);
                }
                else
                {
                    mag.z(i, j) = b;
                }
            }
        }
        return status;
    }

    // A is comples

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (A.z(i, j).Mag() > b.Mag())
            {
                mag.z(i, j) = A.z(i, j);
            }
            else if (A.z(i, j).Mag() < b.Mag())
            {
                mag.z(i, j) = b;
            }
            else if (A.z(i, j).Arg() > b.Arg())
            {
                mag.z(i, j) = A.z(i, j);
            }
            else
            {
                mag.z(i, j) = b;
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes a matrix whose values are the maxima of the arguments and
// returns status
//------------------------------------------------------------------------------
hwMathStatus Max(const hwMatrix& A,
                 const hwMatrix& B,
                 hwMatrix&       mag)
{
    if (A.M() != B.M() || A.N() != B.N())
    {
        if (A.M() == 1 && A.N() == 1)
        {
            if (A.IsReal())
            {
                return Max(B, A(0), mag);
            }
            return Max(B, A.z(0), mag);
        }
        else if (B.M() == 1 && B.N() == 1)
        {
            if (B.IsReal())
            {
                return Max(A, B(0), mag);
            }
            return Max(A, B.z(0), mag);
        }
        else
        {
            return hwMathStatus(HW_MATH_ERR_ARRAYSIZE, 1, 2);
        }
    }

    hwMathStatus status;

    int m = A.M();
    int n = A.N();

    if (A.IsReal())
    {
        if (B.IsReal())
        {
            status = mag.Dimension(m, n, hwMatrix::REAL);

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

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    mag(i, j) = (A(i, j) > B(i, j)) ? A(i, j) : B(i, j);
                }
            } 
        }
        else
        {
            double magB;

            status = mag.Dimension(m, n, hwMatrix::COMPLEX);

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

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    magB = B.z(i, j).Mag();

                    if (A(i, j) > magB)
                    {
                        mag.z(i, j).Set(A(i, j), 0.0);
                    }
                    else if (A(i, j) < magB)
                    {
                        mag.z(i, j) = B.z(i, j);
                    }
                    else if (B.z(i, j).Imag() < 0.0)
                    {
                        mag.z(i, j).Set(A(i, j), 0.0);
                    }
                    else
                    {
                        mag.z(i, j) = B.z(i, j);
                    }
                }
            }
        }
        return status;
    }

    // A is complex
    double magA;

    status = mag.Dimension(m, n, hwMatrix::COMPLEX);
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

    if (B.IsReal())
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                magA = A.z(i, j).Mag();

                if (magA < B(i, j))
                {
                    mag.z(i, j).Set(B(i, j), 0.0);
                }
                else if (magA > B(i, j))
                {
                    mag.z(i, j) = A.z(i, j);
                }
                else if (A.z(i, j).Imag() < 0.0)
                {
                    mag.z(i, j).Set(B(i, j), 0.0);
                }
                else
                {
                    mag.z(i, j) = A.z(i, j);
                }
            }
        }
    }
    else
    {
        double magB;

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                magA = A.z(i, j).Mag();
                magB = B.z(i, j).Mag();

                if (magA > magB)
                {
                    mag.z(i, j) = A.z(i, j);
                }
                else if (magA < magB)
                {
                    mag.z(i, j) = B.z(i, j);
                }
                else if (A.z(i, j).Arg() > B.z(i, j).Arg())
                {
                    mag.z(i, j) = A.z(i, j);
                }
                else
                {
                    mag.z(i, j) = B.z(i, j);
                }
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes the minimum absolute value in a real matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus AbsMin(const hwMatrix& A,
                    double*         value,
                    int*            index)
{
    if (A.IsEmpty())
    {
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);
    }

    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int indx = 0;
    int n    = A.Size();

    double min = fabs(A(0));
    double value1;

    for (int i = 1; i < n; i++)
    {
        value1 = fabs(A(i));

        if (value1 < min)
        {
            min  = value1;
            indx = i;
        }
    }

    if (value)
    {
        (*value) = min;
    }

    if (index)
    {
        (*index) = indx;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes the minimum absolute value from each column of a matrix and
// returns status
//------------------------------------------------------------------------------
hwMathStatus AbsMin(const hwMatrix& A,
                    hwMatrix*       mag,
                    hwMatrixI*      row)
{
    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int n = A.N();

    hwMathStatus status;
    if (mag)
    {
        status = (n > 0) ? 
                    mag->Dimension(1, n, hwMatrix::REAL):
                    mag->Dimension(0, 0, hwMatrix::REAL);

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
    }

    if (row)
    {
        status = (n > 0) ? 
                 row->Dimension(1, n, hwMatrixI::REAL):
                 row->Dimension(0, 0, hwMatrixI::REAL);

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
    }

    int    index;
    double min;
    double value;

    for (int j = 0; j < n; j++)
    {
        min   = fabs(A(0, j));
        index = 0;

        for (int i = 1; i < A.M(); i++)
        {
            value = fabs(A(i, j));

            if (value < min)
            {
                min   = value;
                index = i;
            }
        }

        if (mag)
        {
            (*mag)(0, j) = min;
        }

        if (row)
        {
            (*row)(0, j) = index;
        }
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes the maximum absolute value in a real matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus AbsMax(const hwMatrix& A,
                    double*         value,
                    int*            index)
{
    if (A.IsEmpty())
    {
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);
    }

    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int indx = 0;
    int n    = A.Size();

    double max = fabs(A(0));
    double value1;

    for (int i = 1; i < n; i++)
    {
        value1 = fabs(A(i));

        if (value1 > max)
        {
            max  = value1;
            indx = i;
        }
    }

    if (value)
    {
        (*value) = max;
    }

    if (index)
    {
        (*index) = indx;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Computes the maximum absolute value from each column of a matrix and
// returns status
//------------------------------------------------------------------------------
hwMathStatus AbsMax(const hwMatrix& A,
                    hwMatrix*       mag,
                    hwMatrixI*      row)
{
    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    int n = A.N();

    hwMathStatus status;

    if (mag)
    {
        status = (n > 0) ? mag->Dimension(1, n, hwMatrix::REAL) :
                           mag->Dimension(0, 0, hwMatrix::REAL);

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
    }

    if (row)
    {
        status = (n > 0) ? row->Dimension(1, n, hwMatrixI::REAL):
                           row->Dimension(0, 0, hwMatrixI::REAL);

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
    }

    int    index;
    double max;
    double value;

    for (int j = 0; j < n; j++)
    {
        max   = fabs(A(0, j));
        index = 0;

        for (int i = 1; i < A.M(); i++)
        {
            value = fabs(A(i, j));

            if (value > max)
            {
                max   = value;
                index = i;
            }
        }

        if (mag)
        {
            (*mag)(0, j) = max;
        }

        if (row)
        {
            (*row)(0, j) = index;
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes the peaks from a vector and returns status
//------------------------------------------------------------------------------
hwMathStatus Peaks(const hwMatrix& A,
                   hwMatrix*       mag,
                   hwMatrixI*      index,
                   int             option)
{
    hwMathStatus status;

    if (!A.IsReal())
    {
        return status(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!A.IsVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }

    int n = A.Size();
    if (n < 3)
    {
        return status(HW_MATH_ERR_ARRAYSIZE, 1);
    }

    bool getMax = false;
    bool getMin = true;
    if (option == 0.0)
    {
        getMax = true;
        getMin = true;
    }
    else if (option > 0.0)
    {
        getMax = true;
        getMin = false;
    }

    int numPeaks = 0;
    
    double value3;

    // first pass - count the peaks
    double value1 = A(0);
    double value2 = A(1);

    for (int i = 2; i < n; i++)
    {
        value3 = A(i);

        // check for maxima
        if (getMax)
        {
            if (value1 <= value2 && value2 >= value3)
            {
                numPeaks++;
            }
        }

        // check for minima
        if (getMin)
        {
            if (value1 >= value2 && value2 <= value3)
            {
                numPeaks++;
            }
        }

        value1 = value2;
        value2 = value3;
    }

    if (numPeaks == 0)
    {
        return status(HW_MATH_WARN_NOPEAKS, 1);
    }

    // second pass - record the peaks
    if (mag)
    {
        status = mag->Dimension(numPeaks, hwMatrix::REAL);

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
    }

    if (index)
    {
        status = index->Dimension(numPeaks, hwMatrixI::REAL);

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
    }

    numPeaks = 0;
    value1 = A(0);
    value2 = A(1);

    for (int i = 2; i < n; i++)
    {
        value3 = A(i);

        // check for maxima
        if (getMax)
        {
            if (value1 <= value2)
            {
                if (value2 >= value3)
                {
                    if (mag)
                    {
                        (*mag)(numPeaks) = value2;
                    }

                    if (index)
                    {
                        (*index)(numPeaks) = i - 1;
                    }
                    numPeaks++;
                }
            }
        }

        // check for minima
        if (getMin)
        {
            if (value1 >= value2)
            {
                if (value2 <= value3)
                {
                    if (mag)
                    {
                        (*mag)(numPeaks) = value2;
                    }
                    if (index)
                    {
                        (*index)(numPeaks) = i - 1;
                    }
                    numPeaks++;
                }
            }
        }

        value1 = value2;
        value2 = value3;
    }

    return status;
}
//------------------------------------------------------------------------------
// Finds the first index where A(index) = value and returns status
//------------------------------------------------------------------------------
hwMathStatus GetIndex(const hwMatrix& A,
                      double          value,
                      int&            index,
                      double          tol)
{
    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (!A.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (tol < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 4);
    }

    int size = A.Size();

    for (int i = 0; i < size; i++)
    {
        if (IsZero(A(i) - value, tol))
        {
            index = i;
            return hwMathStatus();
        }
    }

    return hwMathStatus(HW_MATH_ERR_DATANOTFOUND, 1, 2);
}
//------------------------------------------------------------------------------
// Finds the first pair, (i, j), where A(i,j) = value and returns status
//------------------------------------------------------------------------------
hwMathStatus GetIndexPair(const hwMatrix& A, 
                          double          value, 
                          int&            ii, 
                          int&            jj,
                          double          tol, 
                          bool            byCol)
{
    if (!A.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }

    if (tol < 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NEGATIVE, 5);
    }

    int m = A.M();
    int n = A.N();

    if (byCol)
    {
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < m; i++)
            {
                if (IsZero(A(i, j) - value, tol))
                {
                    ii = i;
                    jj = j;
                    return hwMathStatus();
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (IsZero(A(i, j) - value, tol))
                {
                    ii = i;
                    jj = j;
                    return hwMathStatus();
                }
            }
        }
    }

    return hwMathStatus(HW_MATH_ERR_DATANOTFOUND, 1, 2);
}
//------------------------------------------------------------------------------
// Performs a bisection search on a sorted data array and returns status
//------------------------------------------------------------------------------
int BinarySearch(const double *data,
                 int          npts,
                 double       find)
{
    if (!data)
    {
        return 0;
    }

    if (find < data[0])
    {
        return -1;
    }

    int num_points = npts;
    int temp_idx = num_points / 2;
    int max_idx = num_points - 1;
    int min_idx = 0;

    if (find > data[max_idx])
    {
        return max_idx;
    }

    for (int cnt = 0; cnt < num_points; cnt++)  // allow for bail-out
    {
        if (data[temp_idx] == find)
        {
            if (temp_idx - 1 >= min_idx && data[temp_idx - 1] == find)
            {
                return temp_idx - 1;
            }
            else
            {
                return temp_idx;
            }
        }
        else if (data[temp_idx] < find)
        {
            if (temp_idx == max_idx)
            {
                return temp_idx;
            }

            if (data[temp_idx + 1] > find)
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

            if (data[temp_idx - 1] < find)
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
// Computes the bin locations to create a histogram and returns status
//------------------------------------------------------------------------------
hwMathStatus Bins(const hwMatrix& data,
                  hwMatrix&       bin)
{
    if (data.IsEmpty())
    {
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);
    }
    if (!data.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 1);
    }
    if (!data.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 1);
    }

    if (!bin.IsReal())
    {
        return hwMathStatus(HW_MATH_ERR_COMPLEX, 2);
    }
    if (!bin.IsVector())
    {
        return hwMathStatus(HW_MATH_ERR_VECTOR, 2);
    }

    int n       = data.Size();
    int numBins = bin.Size();

    double temp_max = data(0);
    double temp_min = data(0);

    for (int i = 0; i < n; i++)
    {
        if (data(i) > temp_max)
        {
            temp_max = data(i);
        }
        else if (data(i) < temp_min)
        {
            temp_min = data(i);
        }
    }

    double bin_size = (temp_max - temp_min) / static_cast<double>(numBins);

    for (int i = 0; i < numBins; i++)
    {
        bin(i) = (i + 0.5) * bin_size + temp_min;
    }

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Sorts a real or complex vector in the specified direction, returning the 
// sorted values, indices and status
//------------------------------------------------------------------------------
static hwMathStatus SortVector(const hwMatrix& unsorted, 
                               hwMatrix&       sorted, 
                               hwMatrixI*      index,
                               bool            ascend)
{
    hwMathStatus status = (unsorted.IsReal()) ?
        sorted.Dimension(unsorted.M(), unsorted.N(), hwMatrix::REAL):
        sorted.Dimension(unsorted.M(), unsorted.N(), hwMatrix::COMPLEX);

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

    if (index)
    {
        status = index->Dimension(unsorted.M(), unsorted.N(), hwMatrixI::REAL);

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
    }

    if (!unsorted.IsEmptyOrVector())
    {
        return status(HW_MATH_ERR_VECTOR, 1);
    }

    int n = unsorted.Size();
    std::vector<int> temp(n);

    for (int i = 0; i < n; i++)
    {
        temp[i] = i;
    }

    hwMatrix* unsorted_nc = (hwMatrix*) &unsorted;

    data_priv_real = unsorted_nc->GetRealData();

    if (unsorted.IsReal())
    {
        if (ascend)
        {
            sort(temp.begin(), temp.end(), AscendReal);
        }
        else
        {
            sort(temp.begin(), temp.end(), DescendReal);
        }

        if (index)
        {
            for (int i = 0; i < n; i++)
            {
                (*index)(i) = temp[i];
                sorted(i)   = unsorted(temp[i]);
            }
        }
        else
        {
            // TODO: no need to use "temp" if indices not wanted
            for (int i = 0; i < n; i++)
            {
                sorted(i) = unsorted(temp[i]);
            }
        }
    }
    else
    {
        hwMatrix mag(unsorted.M(), unsorted.N(), hwMatrix::REAL);

        data_priv_mag     = mag.GetRealData();
        data_priv_complex = unsorted_nc->GetComplexData();

        for (int i = 0; i < n; i++)
        {
            mag(i) = data_priv_complex[i].Mag();
        }

        if (ascend)
        {
            sort(temp.begin(), temp.end(), AscendComplex);
        }
        else
        {
            sort(temp.begin(), temp.end(), DescendComplex);
        }

        if (index)
        {
            for (int i = 0; i < n; i++)
            {
                (*index)(i) = temp[i];
                sorted.z(i) = unsorted.z(temp[i]);
            }
        }
        else
        {
            // TODO: no need to use "temp" if indices not wanted
            for (int i = 0; i < n; i++)
            {
                sorted.z(i) = unsorted.z(temp[i]);
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Sorts a matrix by column in a given direction, returning the sorted values,
// indices and status
//------------------------------------------------------------------------------
hwMathStatus Sort(const hwMatrix& unsorted, 
                  hwMatrix&       sorted,
                  hwMatrixI*      index,
                  bool            ascend)
{
    if (unsorted.IsEmptyOrVector())
    {
        return SortVector(unsorted, sorted, index, ascend);
    }

	hwMathStatus status;
    if (unsorted.IsReal())
    {
        status = sorted.Dimension(unsorted.M(), unsorted.N(), hwMatrix::REAL);
    }
    else
    {
        status = sorted.Dimension(unsorted.M(), unsorted.N(), hwMatrix::COMPLEX);
    }

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

    if (index)
    {
        status = index->Dimension(unsorted.M(), unsorted.N(), hwMatrixI::REAL);

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
    }

    for (int k = 0; k < unsorted.N(); ++k)
    {
        if (unsorted.IsReal())
        {
            const double* unsortedColumnPtr = &(unsorted(0, k));
            double* sortedColumnPtr = &(sorted(0, k));

            const hwMatrix unsortedColumn(unsorted.M(), 
                (void*) unsortedColumnPtr, hwMatrix::REAL);
            hwMatrix sortedColumn(sorted.M(), 
                sortedColumnPtr, hwMatrix::REAL);

            if (index)
            {
                int* indexColumnPtr = &((*index)(0, k));
                hwMatrixI indexColumn(index->M(), indexColumnPtr, hwMatrixI::REAL);
                status = SortVector(unsortedColumn, sortedColumn, &indexColumn, ascend);
            }
            else
            {
                status = SortVector(unsortedColumn, sortedColumn, NULL, ascend);
            }
        }
        else    // complex
        {
            const hwComplex* unsortedColumnPtr = &(unsorted.z(0, k));
            hwComplex* sortedColumnPtr = &(sorted.z(0, k));

            const hwMatrix unsortedColumn(unsorted.M(), 
                (void*) unsortedColumnPtr, hwMatrix::COMPLEX);
            hwMatrix sortedColumn(sorted.M(), 
                sortedColumnPtr, hwMatrix::COMPLEX);

            if (index)
            {
                int* indexColumnPtr = &((*index)(0, k));
                hwMatrixI indexColumn(index->M(), indexColumnPtr, hwMatrixI::REAL);
                status = SortVector(unsortedColumn, sortedColumn, &indexColumn, ascend);
            }
            else
            {
                status = SortVector(unsortedColumn, sortedColumn, NULL, ascend);
            }
        }
    }

    return status;
}
//------------------------------------------------------------------------------
// Computes the continued fraction approximation of a value and returns
// the status
//------------------------------------------------------------------------------
hwMathStatus ContFrac(double     value, 
                      double     tol, 
                      double&    num,
                      double&    den,
                      hwMatrix&  cf)
{
    if (tol <= 0.0)
    {
        return hwMathStatus(HW_MATH_ERR_NONPOSITIVE, 2);
    }

    hwMathStatus status = cf.Dimension(1, hwMatrix::REAL);

    int    numTerms = 1;
    double intgr    = RoundT(value);
    double frac     = value - intgr;
    double n_old2   = 1.0;
    double d_old2   = 0.0;

    num   = intgr;
    den   = 1;
    cf(0) = intgr;

    while (fabs(value - num/den) >= tol)
    {
        numTerms++;
        status = cf.Resize(numTerms);

        double recip = 1.0 / frac;
        intgr = RoundT(recip);

        cf(numTerms-1) = intgr;
        frac = recip - intgr;

        double n_old1 = num;
        double d_old1 = den;
        num = num * cf(numTerms-1) + n_old2;
        den = den * cf(numTerms-1) + d_old2;
        n_old2 = n_old1;
        d_old2 = d_old1;
    }

    if (den < 0)
    {
        num = -num;
        den = -den;
    }

    return hwMathStatus();
}
