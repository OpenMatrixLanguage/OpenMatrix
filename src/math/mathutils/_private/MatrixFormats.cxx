/**
* @file MatrixFormats.cxx
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
#include <MatrixFormats.h>
#include <hwMathStatus.h>

// The following functions are used to store matrices in special formats
// that are used by LAPACK and Sundials third party libraries

//------------------------------------------------------------------------------
// Sizes a dense matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus SizeDenseMatrix(int       m,
                             int       n,
                             hwMatrix& A,
                             bool      initZero)
{
    // prepare A to store a dense matrix by column

    hwMathStatus status = A.Dimension(m, n, hwMatrix::REAL);

    if (!status.IsOk())
    {
        switch (status.GetArg1())
        {
            case 0:  status.SetArg1(3);  break;
            case 1:  status.SetArg1(2);  break;
            case 2:  status.SetArg1(1);  break;
            default: status.ResetArgs(); break;
        }
        return status;
    }

    if (initZero)
    {
        A.SetElements(0.0);
    }
    return status;
}
//------------------------------------------------------------------------------
// Sizes a band matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus SizeBandMatrix(int       n,
                            int       kl,
                            int       ku,
                            hwMatrix& A,
                            bool      initZero)
{
    if (n < 1)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 1);
    }

    if (kl < 0 || kl > n - 1)
    {
        return hwMathStatus(HW_MATH_ERR_KLUD, 2);
    }

    if (ku < 0 || ku > n - 1)
    {
        return hwMathStatus(HW_MATH_ERR_KLUD, 3);
    }

    // prepare A to store a CDS matrix by column
    int ldab = 2 * kl + ku + 1;

    hwMathStatus status = A.Dimension(ldab, n, hwMatrix::REAL);
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

    if (initZero)
    {
        A.SetElements(0.0);
    }

    return status;
}
//------------------------------------------------------------------------------
// Sizes a symmetric band matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus SizeSymBandMatrix(int       n,
                               int       kd,
                               hwMatrix& A,
                               bool      initZero)
{
    if (n < 1)
    {
        return hwMathStatus(HW_MATH_ERR_ARRAYDIM, 1);
    }

    if (kd < 0 || kd > n - 1)
    {
        return hwMathStatus(HW_MATH_ERR_KLUD, 2);
    }

    // prepare A to store a symmetric CDS matrix by column
    int ldab = kd + 1;

    hwMathStatus status = A.Dimension(ldab, n, hwMatrix::REAL);
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

    if (initZero)
    {
        A.SetElements(0.0);
    }
    return status;
}
//------------------------------------------------------------------------------
// Sets an element in a dense matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus SetDenseMatrixElem(int       i,
                                int       j,
                                double    value,
                                hwMatrix& A)
{
    int row = i;
    int col = j;

    A(row, col) = value;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Sets an element in a band matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus SetBandMatrixElem(int       i, 
                               int       j, 
                               double    value, 
                               hwMatrix& A,
                               int       kl, 
                               int       ku)
{
    if (kl < 0 || kl > A.M() - 1)
    {
        return hwMathStatus(HW_MATH_ERR_KLUD, 5);
    }
    if (ku < 0 || ku > A.M() - 1)
    {
        return hwMathStatus(HW_MATH_ERR_KLUD, 6);
    }

    int row = kl + ku + i - j;
    int col = j;

    A(row, col) = value;

    return hwMathStatus();
}
//------------------------------------------------------------------------------
// Sets an element in a symmetric band matrix and returns status
//------------------------------------------------------------------------------
hwMathStatus SetSymBandMatrixElem(int       i, 
                                  int       j,
                                  double    value, 
                                  hwMatrix& A,
                                  int       kd)
{
    if (kd < 0 || kd > A.M() - 1)
    {
        return hwMathStatus(HW_MATH_ERR_KLUD, 5);
    }

    int row = i - j;
    int col = j;

    A(row, col) = value;

    return hwMathStatus();
}


