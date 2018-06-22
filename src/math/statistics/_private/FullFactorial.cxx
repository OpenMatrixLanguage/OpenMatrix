/**
* @file FullFactorial.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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
#include "FullFactorial.h"

#include "hwMatrix.h"

//------------------------------------------------------------------------------
// Construct a full factorial design matrix and return status
//------------------------------------------------------------------------------
hwMathStatus FullFact(const hwMatrixI& runs, hwMatrixI& matrix)
{
    int n = runs.Size();  // Number of factors
    if (n < 1)
    {
        return hwMathStatus(HW_MATH_ERR_EMPTYMATRIX, 1);
    }

    int m = runs(0);
    if (m < 2)
    {
        return hwMathStatus(HW_MATH_ERR_DOELEVEL, 1);
    }

    for (int i = 1; i < n; i++)
    {
        if (runs(i) < 2)
        {
            return hwMathStatus(HW_MATH_ERR_DOELEVEL, 1);
        }
        m *= runs(i);
    }

    hwMathStatus status = matrix.Dimension(m, n, hwMatrixI::REAL);
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

    for (int j = 0; j < n; j++)     // first row is all ones
    {
        matrix(0, j) = 1;
    }

    for (int i = 1; i < m; i++)
    {
        matrix(i, 0) = matrix(i-1, 0) % runs(0) + 1;                // cycle first column

        for (int j = 1; j < n; j++)
        {
            if (matrix(i, j-1) == 1 && matrix(i-1, j-1) != 1)       // prev column switched to 1
            {
                matrix(i, j) = matrix(i-1, j) % runs(j) + 1;        // cycle the column
            }
            else
            {
                matrix(i, j) = matrix(i-1, j);                      // carry prev value down
            }
        }
    }

    return status;
}
