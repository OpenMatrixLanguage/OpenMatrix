/**
* @file BoxBehnken.cxx
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
#include "BoxBehnken.h"

#include "hwMatrix.h"
#include "hwMathStatus.h"

//------------------------------------------------------------------------------
// Returns status after constructing a Box-Behnken design matrix
//------------------------------------------------------------------------------
hwMathStatus BoxBehnken(int n, hwMatrixI& A)
{
    int m;
    switch (n)
    {
        case 3: m = 15; break;
        
        case 4: m = 27; break;

        case 5: m = 46; break;
    
        case 6: m = 54; break;
        
        case 7: m = 62; break;
        
        default: return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 1); break;
    };

    hwMathStatus status = A.Dimension(m, n, hwMatrixI::REAL);

    if (!status.IsOk())
    {
        status.ResetArgs();
        return status;
    }

    switch (n)
    {
    case 3:
        A(0,0) = -1;    A(0,1) = -1;    A(0,2) =  0;
        A(1,0) = -1;    A(1,1) =  1;    A(1,2) =  0;
        A(2,0) =  1;    A(2,1) = -1;    A(2,2) =  0;
        A(3,0) =  1;    A(3,1) =  1;    A(3,2) =  0;
        A(4,0) = -1;    A(4,1) =  0;    A(4,2) = -1;
        A(5,0) = -1;    A(5,1) =  0;    A(5,2) =  1;
        A(6,0) =  1;    A(6,1) =  0;    A(6,2) = -1;
        A(7,0) =  1;    A(7,1) =  0;    A(7,2) =  1;
        A(8,0) =  0;    A(8,1) = -1;    A(8,2) = -1;
        A(9,0) =  0;    A(9,1) = -1;    A(9,2) =  1;
        A(10,0) = 0;    A(10,1) =  1;   A(10,2) = -1;
        A(11,0) = 0;    A(11,1) =  1;   A(11,2) =  1;
        A(12,0) = 0;    A(12,1) =  0;   A(12,2) =  0;
        A(13,0) = 0;    A(13,1) =  0;   A(13,2) =  0;
        A(14,0) = 0;    A(14,1) =  0;   A(14,2) =  0;
        break;
    case 4:
        A(0,0) = -1;    A(0,1) = -1;    A(0,2) =  0;    A(0,3) =  0;
        A(1,0) = -1;    A(1,1) =  1;    A(1,2) =  0;    A(1,3) =  0;
        A(2,0) =  1;    A(2,1) = -1;    A(2,2) =  0;    A(2,3) =  0;
        A(3,0) =  1;    A(3,1) =  1;    A(3,2) =  0;    A(3,3) =  0;
        A(4,0) =  0;    A(4,1) =  0;    A(4,2) = -1;    A(4,3) = -1;
        A(5,0) =  0;    A(5,1) =  0;    A(5,2) = -1;    A(5,3) =  1;
        A(6,0) =  0;    A(6,1) =  0;    A(6,2) =  1;    A(6,3) = -1;
        A(7,0) =  0;    A(7,1) =  0;    A(7,2) =  1;    A(7,3) =  1;
        A(8,0) = -1;    A(8,1) =  0;    A(8,2) =  0;    A(8,3) = -1;
        A(9,0) = -1;    A(9,1) =  0;    A(9,2) =  0;    A(9,3) =  1;
        A(10,0) =  1;   A(10,1) =  0;   A(10,2) =  0;   A(10,3) = -1;
        A(11,0) =  1;   A(11,1) =  0;   A(11,2) =  0;   A(11,3) =  1;
        A(12,0) =  0;   A(12,1) = -1;   A(12,2) = -1;   A(12,3) =  0;
        A(13,0) =  0;   A(13,1) = -1;   A(13,2) =  1;   A(13,3) =  0;
        A(14,0) =  0;   A(14,1) =  1;   A(14,2) = -1;   A(14,3) =  0;
        A(15,0) =  0;   A(15,1) =  1;   A(15,2) =  1;   A(15,3) =  0;
        A(16,0) = -1;   A(16,1) =  0;   A(16,2) = -1;   A(16,3) =  0;
        A(17,0) = -1;   A(17,1) =  0;   A(17,2) =  1;   A(17,3) =  0;
        A(18,0) =  1;   A(18,1) =  0;   A(18,2) = -1;   A(18,3) =  0;
        A(19,0) =  1;   A(19,1) =  0;   A(19,2) =  1;   A(19,3) =  0;
        A(20,0) =  0;   A(20,1) = -1;   A(20,2) =  0;   A(20,3) = -1;
        A(21,0) =  0;   A(21,1) = -1;   A(21,2) =  0;   A(21,3) =  1;
        A(22,0) =  0;   A(22,1) =  1;   A(22,2) =  0;   A(22,3) = -1;
        A(23,0) =  0;   A(23,1) =  1;   A(23,2) =  0;   A(23,3) =  1;
        A(24,0) =  0;   A(24,1) =  0;   A(24,2) =  0;   A(24,3) =  0;
        A(25,0) =  0;   A(25,1) =  0;   A(25,2) =  0;   A(25,3) =  0;
        A(26,0) =  0;   A(26,1) =  0;   A(26,2) =  0;   A(26,3) =  0;
        break;
    case 5:
        A(0,0) = -1;    A(0,1) = -1;    A(0,2) =  0;    A(0,3) =  0;    A(0,4) =  0;
        A(1,0) = -1;    A(1,1) =  1;    A(1,2) =  0;    A(1,3) =  0;    A(1,4) =  0;
        A(2,0) =  1;    A(2,1) = -1;    A(2,2) =  0;    A(2,3) =  0;    A(2,4) =  0;
        A(3,0) =  1;    A(3,1) =  1;    A(3,2) =  0;    A(3,3) =  0;    A(3,4) =  0;
        A(4,0) =  0;    A(4,1) =  0;    A(4,2) = -1;    A(4,3) = -1;    A(4,4) =  0;
        A(5,0) =  0;    A(5,1) =  0;    A(5,2) = -1;    A(5,3) =  1;    A(5,4) =  0;
        A(6,0) =  0;    A(6,1) =  0;    A(6,2) =  1;    A(6,3) = -1;    A(6,4) =  0;
        A(7,0) =  0;    A(7,1) =  0;    A(7,2) =  1;    A(7,3) =  1;    A(7,4) =  0;
        A(8,0) =  0;    A(8,1) = -1;    A(8,2) =  0;    A(8,3) =  0;    A(8,4) = -1;
        A(9,0) =  0;    A(9,1) = -1;    A(9,2) =  0;    A(9,3) =  0;    A(9,4) =  1;
        A(10,0) =  0;   A(10,1) =  1;   A(10,2) =  0;   A(10,3) =  0;   A(10,4) = -1;
        A(11,0) =  0;   A(11,1) =  1;   A(11,2) =  0;   A(11,3) =  0;   A(11,4) =  1;
        A(12,0) = -1;   A(12,1) =  0;   A(12,2) = -1;   A(12,3) =  0;   A(12,4) =  0;
        A(13,0) = -1;   A(13,1) =  0;   A(13,2) =  1;   A(13,3) =  0;   A(13,4) =  0;
        A(14,0) =  1;   A(14,1) =  0;   A(14,2) = -1;   A(14,3) =  0;   A(14,4) =  0;
        A(15,0) =  1;   A(15,1) =  0;   A(15,2) =  1;   A(15,3) =  0;   A(15,4) =  0;
        A(16,0) =  0;   A(16,1) =  0;   A(16,2) =  0;   A(16,3) = -1;   A(16,4) = -1;
        A(17,0) =  0;   A(17,1) =  0;   A(17,2) =  0;   A(17,3) = -1;   A(17,4) =  1;
        A(18,0) =  0;   A(18,1) =  0;   A(18,2) =  0;   A(18,3) =  1;   A(18,4) = -1;
        A(19,0) =  0;   A(19,1) =  0;   A(19,2) =  0;   A(19,3) =  1;   A(19,4) =  1;
        A(20,0) =  0;   A(20,1) = -1;   A(20,2) = -1;   A(20,3) =  0;   A(20,4) =  0;
        A(21,0) =  0;   A(21,1) = -1;   A(21,2) =  1;   A(21,3) =  0;   A(21,4) =  0;
        A(22,0) =  0;   A(22,1) =  1;   A(22,2) = -1;   A(22,3) =  0;   A(22,4) =  0;
        A(23,0) =  0;   A(23,1) =  1;   A(23,2) =  1;   A(23,3) =  0;   A(23,4) =  0;
        A(24,0) = -1;   A(24,1) =  0;   A(24,2) =  0;   A(24,3) = -1;   A(24,4) =  0;
        A(25,0) = -1;   A(25,1) =  0;   A(25,2) =  0;   A(25,3) =  1;   A(25,4) =  0;
        A(26,0) =  1;   A(26,1) =  0;   A(26,2) =  0;   A(26,3) = -1;   A(26,4) =  0;
        A(27,0) =  1;   A(27,1) =  0;   A(27,2) =  0;   A(27,3) =  1;   A(27,4) =  0;
        A(28,0) =  0;   A(28,1) =  0;   A(28,2) = -1;   A(28,3) =  0;   A(28,4) = -1;
        A(29,0) =  0;   A(29,1) =  0;   A(29,2) = -1;   A(29,3) =  0;   A(29,4) =  1;
        A(30,0) =  0;   A(30,1) =  0;   A(30,2) =  1;   A(30,3) =  0;   A(30,4) = -1;
        A(31,0) =  0;   A(31,1) =  0;   A(31,2) =  1;   A(31,3) =  0;   A(31,4) =  1;
        A(32,0) = -1;   A(32,1) =  0;   A(32,2) =  0;   A(32,3) =  0;   A(32,4) = -1;
        A(33,0) = -1;   A(33,1) =  0;   A(33,2) =  0;   A(33,3) =  0;   A(33,4) =  1;
        A(34,0) =  1;   A(34,1) =  0;   A(34,2) =  0;   A(34,3) =  0;   A(34,4) = -1;
        A(35,0) =  1;   A(35,1) =  0;   A(35,2) =  0;   A(35,3) =  0;   A(35,4) =  1;
        A(36,0) =  0;   A(36,1) = -1;   A(36,2) =  0;   A(36,3) = -1;   A(36,4) =  0;
        A(37,0) =  0;   A(37,1) = -1;   A(37,2) =  0;   A(37,3) =  1;   A(37,4) =  0;
        A(38,0) =  0;   A(38,1) =  1;   A(38,2) =  0;   A(38,3) = -1;   A(38,4) =  0;
        A(39,0) =  0;   A(39,1) =  1;   A(39,2) =  0;   A(39,3) =  1;   A(39,4) =  0;
        A(40,0) =  0;   A(40,1) =  0;   A(40,2) =  0;   A(40,3) =  0;   A(40,4) =  0;
        A(41,0) =  0;   A(41,1) =  0;   A(41,2) =  0;   A(41,3) =  0;   A(41,4) =  0;
        A(42,0) =  0;   A(42,1) =  0;   A(42,2) =  0;   A(42,3) =  0;   A(42,4) =  0;
        A(43,0) =  0;   A(43,1) =  0;   A(43,2) =  0;   A(43,3) =  0;   A(43,4) =  0;
        A(44,0) =  0;   A(44,1) =  0;   A(44,2) =  0;   A(44,3) =  0;   A(44,4) =  0;
        A(45,0) =  0;   A(45,1) =  0;   A(45,2) =  0;   A(45,3) =  0;   A(45,4) =  0;
        break;
    case 6:
        A(0,0) = -1;    A(0,1) = -1;    A(0,2) =  0;    A(0,3) = -1;    A(0,4) =  0;    A(0,5) =  0;
        A(1,0) = -1;    A(1,1) = -1;    A(1,2) =  0;    A(1,3) =  1;    A(1,4) =  0;    A(1,5) =  0;
        A(2,0) = -1;    A(2,1) =  1;    A(2,2) =  0;    A(2,3) = -1;    A(2,4) =  0;    A(2,5) =  0;
        A(3,0) = -1;    A(3,1) =  1;    A(3,2) =  0;    A(3,3) =  1;    A(3,4) =  0;    A(3,5) =  0;
        A(4,0) =  1;    A(4,1) = -1;    A(4,2) =  0;    A(4,3) = -1;    A(4,4) =  0;    A(4,5) =  0;
        A(5,0) =  1;    A(5,1) = -1;    A(5,2) =  0;    A(5,3) =  1;    A(5,4) =  0;    A(5,5) =  0;
        A(6,0) =  1;    A(6,1) =  1;    A(6,2) =  0;    A(6,3) = -1;    A(6,4) =  0;    A(6,5) =  0;
        A(7,0) =  1;    A(7,1) =  1;    A(7,2) =  0;    A(7,3) =  1;    A(7,4) =  0;    A(7,5) =  0;
        A(8,0) =  0;    A(8,1) = -1;    A(8,2) = -1;    A(8,3) =  0;    A(8,4) = -1;    A(8,5) =  0;
        A(9,0) =  0;    A(9,1) = -1;    A(9,2) = -1;    A(9,3) =  0;    A(9,4) =  1;    A(9,5) =  0;
        A(10,0) =  0;   A(10,1) = -1;   A(10,2) =  1;   A(10,3) =  0;   A(10,4) = -1;   A(10,5) =  0;
        A(11,0) =  0;   A(11,1) = -1;   A(11,2) =  1;   A(11,3) =  0;   A(11,4) =  1;   A(11,5) =  0;
        A(12,0) =  0;   A(12,1) =  1;   A(12,2) = -1;   A(12,3) =  0;   A(12,4) = -1;   A(12,5) =  0;
        A(13,0) =  0;   A(13,1) =  1;   A(13,2) = -1;   A(13,3) =  0;   A(13,4) =  1;   A(13,5) =  0;
        A(14,0) =  0;   A(14,1) =  1;   A(14,2) =  1;   A(14,3) =  0;   A(14,4) = -1;   A(14,5) =  0;
        A(15,0) =  0;   A(15,1) =  1;   A(15,2) =  1;   A(15,3) =  0;   A(15,4) =  1;   A(15,5) =  0;
        A(16,0) =  0;   A(16,1) =  0;   A(16,2) = -1;   A(16,3) = -1;   A(16,4) =  0;   A(16,5) = -1;
        A(17,0) =  0;   A(17,1) =  0;   A(17,2) = -1;   A(17,3) = -1;   A(17,4) =  0;   A(17,5) =  1;
        A(18,0) =  0;   A(18,1) =  0;   A(18,2) = -1;   A(18,3) =  1;   A(18,4) =  0;   A(18,5) = -1;
        A(19,0) =  0;   A(19,1) =  0;   A(19,2) = -1;   A(19,3) =  1;   A(19,4) =  0;   A(19,5) =  1;
        A(20,0) =  0;   A(20,1) =  0;   A(20,2) =  1;   A(20,3) = -1;   A(20,4) =  0;   A(20,5) = -1;
        A(21,0) =  0;   A(21,1) =  0;   A(21,2) =  1;   A(21,3) = -1;   A(21,4) =  0;   A(21,5) =  1;
        A(22,0) =  0;   A(22,1) =  0;   A(22,2) =  1;   A(22,3) =  1;   A(22,4) =  0;   A(22,5) = -1;
        A(23,0) =  0;   A(23,1) =  0;   A(23,2) =  1;   A(23,3) =  1;   A(23,4) =  0;   A(23,5) =  1;
        A(24,0) = -1;   A(24,1) =  0;   A(24,2) =  0;   A(24,3) = -1;   A(24,4) = -1;   A(24,5) =  0;
        A(25,0) = -1;   A(25,1) =  0;   A(25,2) =  0;   A(25,3) = -1;   A(25,4) =  1;   A(25,5) =  0;
        A(26,0) = -1;   A(26,1) =  0;   A(26,2) =  0;   A(26,3) =  1;   A(26,4) = -1;   A(26,5) =  0;
        A(27,0) = -1;   A(27,1) =  0;   A(27,2) =  0;   A(27,3) =  1;   A(27,4) =  1;   A(27,5) =  0;
        A(28,0) =  1;   A(28,1) =  0;   A(28,2) =  0;   A(28,3) = -1;   A(28,4) = -1;   A(28,5) =  0;
        A(29,0) =  1;   A(29,1) =  0;   A(29,2) =  0;   A(29,3) = -1;   A(29,4) =  1;   A(29,5) =  0;
        A(30,0) =  1;   A(30,1) =  0;   A(30,2) =  0;   A(30,3) =  1;   A(30,4) = -1;   A(30,5) =  0;
        A(31,0) =  1;   A(31,1) =  0;   A(31,2) =  0;   A(31,3) =  1;   A(31,4) =  1;   A(31,5) =  0;
        A(32,0) =  0;   A(32,1) = -1;   A(32,2) =  0;   A(32,3) =  0;   A(32,4) = -1;   A(32,5) = -1;
        A(33,0) =  0;   A(33,1) = -1;   A(33,2) =  0;   A(33,3) =  0;   A(33,4) = -1;   A(33,5) =  1;
        A(34,0) =  0;   A(34,1) = -1;   A(34,2) =  0;   A(34,3) =  0;   A(34,4) =  1;   A(34,5) = -1;
        A(35,0) =  0;   A(35,1) = -1;   A(35,2) =  0;   A(35,3) =  0;   A(35,4) =  1;   A(35,5) =  1;
        A(36,0) =  0;   A(36,1) =  1;   A(36,2) =  0;   A(36,3) =  0;   A(36,4) = -1;   A(36,5) = -1;
        A(37,0) =  0;   A(37,1) =  1;   A(37,2) =  0;   A(37,3) =  0;   A(37,4) = -1;   A(37,5) =  1;
        A(38,0) =  0;   A(38,1) =  1;   A(38,2) =  0;   A(38,3) =  0;   A(38,4) =  1;   A(38,5) = -1;
        A(39,0) =  0;   A(39,1) =  1;   A(39,2) =  0;   A(39,3) =  0;   A(39,4) =  1;   A(39,5) =  1;
        A(40,0) = -1;   A(40,1) =  0;   A(40,2) = -1;   A(40,3) =  0;   A(40,4) =  0;   A(40,5) = -1;
        A(41,0) = -1;   A(41,1) =  0;   A(41,2) = -1;   A(41,3) =  0;   A(41,4) =  0;   A(41,5) =  1;
        A(42,0) = -1;   A(42,1) =  0;   A(42,2) =  1;   A(42,3) =  0;   A(42,4) =  0;   A(42,5) = -1;
        A(43,0) = -1;   A(43,1) =  0;   A(43,2) =  1;   A(43,3) =  0;   A(43,4) =  0;   A(43,5) =  1;
        A(44,0) =  1;   A(44,1) =  0;   A(44,2) = -1;   A(44,3) =  0;   A(44,4) =  0;   A(44,5) = -1;
        A(45,0) =  1;   A(45,1) =  0;   A(45,2) = -1;   A(45,3) =  0;   A(45,4) =  0;   A(45,5) =  1;
        A(46,0) =  1;   A(46,1) =  0;   A(46,2) =  1;   A(46,3) =  0;   A(46,4) =  0;   A(46,5) = -1;
        A(47,0) =  1;   A(47,1) =  0;   A(47,2) =  1;   A(47,3) =  0;   A(47,4) =  0;   A(47,5) =  1;
        A(48,0) =  0;   A(48,1) =  0;   A(48,2) =  0;   A(48,3) =  0;   A(48,4) =  0;   A(48,5) =  0;
        A(49,0) =  0;   A(49,1) =  0;   A(49,2) =  0;   A(49,3) =  0;   A(49,4) =  0;   A(49,5) =  0;
        A(50,0) =  0;   A(50,1) =  0;   A(50,2) =  0;   A(50,3) =  0;   A(50,4) =  0;   A(50,5) =  0;
        A(51,0) =  0;   A(51,1) =  0;   A(51,2) =  0;   A(51,3) =  0;   A(51,4) =  0;   A(51,5) =  0;
        A(52,0) =  0;   A(52,1) =  0;   A(52,2) =  0;   A(52,3) =  0;   A(52,4) =  0;   A(52,5) =  0;
        A(53,0) =  0;   A(53,1) =  0;   A(53,2) =  0;   A(53,3) =  0;   A(53,4) =  0;   A(53,5) =  0;
        break;
    case 7:
        A(0,0) =  0;    A(0,1) =  0;    A(0,2) =  0;    A(0,3) = -1;    A(0,4) = -1;    A(0,5) = -1;    A(0,6) =  0;
        A(1,0) =  0;    A(1,1) =  0;    A(1,2) =  0;    A(1,3) = -1;    A(1,4) = -1;    A(1,5) =  1;    A(1,6) =  0;
        A(2,0) =  0;    A(2,1) =  0;    A(2,2) =  0;    A(2,3) = -1;    A(2,4) =  1;    A(2,5) = -1;    A(2,6) =  0;
        A(3,0) =  0;    A(3,1) =  0;    A(3,2) =  0;    A(3,3) = -1;    A(3,4) =  1;    A(3,5) =  1;    A(3,6) =  0;
        A(4,0) =  0;    A(4,1) =  0;    A(4,2) =  0;    A(4,3) =  1;    A(4,4) = -1;    A(4,5) = -1;    A(4,6) =  0;
        A(5,0) =  0;    A(5,1) =  0;    A(5,2) =  0;    A(5,3) =  1;    A(5,4) = -1;    A(5,5) =  1;    A(5,6) =  0;
        A(6,0) =  0;    A(6,1) =  0;    A(6,2) =  0;    A(6,3) =  1;    A(6,4) =  1;    A(6,5) = -1;    A(6,6) =  0;
        A(7,0) =  0;    A(7,1) =  0;    A(7,2) =  0;    A(7,3) =  1;    A(7,4) =  1;    A(7,5) =  1;    A(7,6) =  0;
        A(8,0) = -1;    A(8,1) =  0;    A(8,2) =  0;    A(8,3) =  0;    A(8,4) =  0;    A(8,5) = -1;    A(8,6) = -1;
        A(9,0) = -1;    A(9,1) =  0;    A(9,2) =  0;    A(9,3) =  0;    A(9,4) =  0;    A(9,5) = -1;    A(9,6) =  1;
        A(10,0) = -1;   A(10,1) =  0;   A(10,2) =  0;   A(10,3) =  0;   A(10,4) =  0;   A(10,5) =  1;   A(10,6) = -1;
        A(11,0) = -1;   A(11,1) =  0;   A(11,2) =  0;   A(11,3) =  0;   A(11,4) =  0;   A(11,5) =  1;   A(11,6) =  1;
        A(12,0) =  1;   A(12,1) =  0;   A(12,2) =  0;   A(12,3) =  0;   A(12,4) =  0;   A(12,5) = -1;   A(12,6) = -1;
        A(13,0) =  1;   A(13,1) =  0;   A(13,2) =  0;   A(13,3) =  0;   A(13,4) =  0;   A(13,5) = -1;   A(13,6) =  1;
        A(14,0) =  1;   A(14,1) =  0;   A(14,2) =  0;   A(14,3) =  0;   A(14,4) =  0;   A(14,5) =  1;   A(14,6) = -1;
        A(15,0) =  1;   A(15,1) =  0;   A(15,2) =  0;   A(15,3) =  0;   A(15,4) =  0;   A(15,5) =  1;   A(15,6) =  1;
        A(16,0) =  0;   A(16,1) = -1;   A(16,2) =  0;   A(16,3) =  0;   A(16,4) = -1;   A(16,5) =  0;   A(16,6) = -1;
        A(17,0) =  0;   A(17,1) = -1;   A(17,2) =  0;   A(17,3) =  0;   A(17,4) = -1;   A(17,5) =  0;   A(17,6) =  1;
        A(18,0) =  0;   A(18,1) = -1;   A(18,2) =  0;   A(18,3) =  0;   A(18,4) =  1;   A(18,5) =  0;   A(18,6) = -1;
        A(19,0) =  0;   A(19,1) = -1;   A(19,2) =  0;   A(19,3) =  0;   A(19,4) =  1;   A(19,5) =  0;   A(19,6) =  1;
        A(20,0) =  0;   A(20,1) =  1;   A(20,2) =  0;   A(20,3) =  0;   A(20,4) = -1;   A(20,5) =  0;   A(20,6) = -1;
        A(21,0) =  0;   A(21,1) =  1;   A(21,2) =  0;   A(21,3) =  0;   A(21,4) = -1;   A(21,5) =  0;   A(21,6) =  1;
        A(22,0) =  0;   A(22,1) =  1;   A(22,2) =  0;   A(22,3) =  0;   A(22,4) =  1;   A(22,5) =  0;   A(22,6) = -1;
        A(23,0) =  0;   A(23,1) =  1;   A(23,2) =  0;   A(23,3) =  0;   A(23,4) =  1;   A(23,5) =  0;   A(23,6) =  1;
        A(24,0) = -1;   A(24,1) = -1;   A(24,2) =  0;   A(24,3) = -1;   A(24,4) =  0;   A(24,5) =  0;   A(24,6) =  0;
        A(25,0) = -1;   A(25,1) = -1;   A(25,2) =  0;   A(25,3) =  1;   A(25,4) =  0;   A(25,5) =  0;   A(25,6) =  0;
        A(26,0) = -1;   A(26,1) =  1;   A(26,2) =  0;   A(26,3) = -1;   A(26,4) =  0;   A(26,5) =  0;   A(26,6) =  0;
        A(27,0) = -1;   A(27,1) =  1;   A(27,2) =  0;   A(27,3) =  1;   A(27,4) =  0;   A(27,5) =  0;   A(27,6) =  0;
        A(28,0) =  1;   A(28,1) = -1;   A(28,2) =  0;   A(28,3) = -1;   A(28,4) =  0;   A(28,5) =  0;   A(28,6) =  0;
        A(29,0) =  1;   A(29,1) = -1;   A(29,2) =  0;   A(29,3) =  1;   A(29,4) =  0;   A(29,5) =  0;   A(29,6) =  0;
        A(30,0) =  1;   A(30,1) =  1;   A(30,2) =  0;   A(30,3) = -1;   A(30,4) =  0;   A(30,5) =  0;   A(30,6) =  0;
        A(31,0) =  1;   A(31,1) =  1;   A(31,2) =  0;   A(31,3) =  1;   A(31,4) =  0;   A(31,5) =  0;   A(31,6) =  0;
        A(32,0) =  0;   A(32,1) =  0;   A(32,2) = -1;   A(32,3) = -1;   A(32,4) =  0;   A(32,5) =  0;   A(32,6) = -1;
        A(33,0) =  0;   A(33,1) =  0;   A(33,2) = -1;   A(33,3) = -1;   A(33,4) =  0;   A(33,5) =  0;   A(33,6) =  1;
        A(34,0) =  0;   A(34,1) =  0;   A(34,2) = -1;   A(34,3) =  1;   A(34,4) =  0;   A(34,5) =  0;   A(34,6) = -1;
        A(35,0) =  0;   A(35,1) =  0;   A(35,2) = -1;   A(35,3) =  1;   A(35,4) =  0;   A(35,5) =  0;   A(35,6) =  1;
        A(36,0) =  0;   A(36,1) =  0;   A(36,2) =  1;   A(36,3) = -1;   A(36,4) =  0;   A(36,5) =  0;   A(36,6) = -1;
        A(37,0) =  0;   A(37,1) =  0;   A(37,2) =  1;   A(37,3) = -1;   A(37,4) =  0;   A(37,5) =  0;   A(37,6) =  1;
        A(38,0) =  0;   A(38,1) =  0;   A(38,2) =  1;   A(38,3) =  1;   A(38,4) =  0;   A(38,5) =  0;   A(38,6) = -1;
        A(39,0) =  0;   A(39,1) =  0;   A(39,2) =  1;   A(39,3) =  1;   A(39,4) =  0;   A(39,5) =  0;   A(39,6) =  1;
        A(40,0) = -1;   A(40,1) =  0;   A(40,2) = -1;   A(40,3) =  0;   A(40,4) = -1;   A(40,5) =  0;   A(40,6) =  0;
        A(41,0) = -1;   A(41,1) =  0;   A(41,2) = -1;   A(41,3) =  0;   A(41,4) =  1;   A(41,5) =  0;   A(41,6) =  0;
        A(42,0) = -1;   A(42,1) =  0;   A(42,2) =  1;   A(42,3) =  0;   A(42,4) = -1;   A(42,5) =  0;   A(42,6) =  0;
        A(43,0) = -1;   A(43,1) =  0;   A(43,2) =  1;   A(43,3) =  0;   A(43,4) =  1;   A(43,5) =  0;   A(43,6) =  0;
        A(44,0) =  1;   A(44,1) =  0;   A(44,2) = -1;   A(44,3) =  0;   A(44,4) = -1;   A(44,5) =  0;   A(44,6) =  0;
        A(45,0) =  1;   A(45,1) =  0;   A(45,2) = -1;   A(45,3) =  0;   A(45,4) =  1;   A(45,5) =  0;   A(45,6) =  0;
        A(46,0) =  1;   A(46,1) =  0;   A(46,2) =  1;   A(46,3) =  0;   A(46,4) = -1;   A(46,5) =  0;   A(46,6) =  0;
        A(47,0) =  1;   A(47,1) =  0;   A(47,2) =  1;   A(47,3) =  0;   A(47,4) =  1;   A(47,5) =  0;   A(47,6) =  0;
        A(48,0) =  0;   A(48,1) = -1;   A(48,2) = -1;   A(48,3) =  0;   A(48,4) =  0;   A(48,5) = -1;   A(48,6) =  0;
        A(49,0) =  0;   A(49,1) = -1;   A(49,2) = -1;   A(49,3) =  0;   A(49,4) =  0;   A(49,5) =  1;   A(49,6) =  0;
        A(50,0) =  0;   A(50,1) = -1;   A(50,2) =  1;   A(50,3) =  0;   A(50,4) =  0;   A(50,5) = -1;   A(50,6) =  0;
        A(51,0) =  0;   A(51,1) = -1;   A(51,2) =  1;   A(51,3) =  0;   A(51,4) =  0;   A(51,5) =  1;   A(51,6) =  0;
        A(52,0) =  0;   A(52,1) =  1;   A(52,2) = -1;   A(52,3) =  0;   A(52,4) =  0;   A(52,5) = -1;   A(52,6) =  0;
        A(53,0) =  0;   A(53,1) =  1;   A(53,2) = -1;   A(53,3) =  0;   A(53,4) =  0;   A(53,5) =  1;   A(53,6) =  0;
        A(54,0) =  0;   A(54,1) =  1;   A(54,2) =  1;   A(54,3) =  0;   A(54,4) =  0;   A(54,5) = -1;   A(54,6) =  0;
        A(55,0) =  0;   A(55,1) =  1;   A(55,2) =  1;   A(55,3) =  0;   A(55,4) =  0;   A(55,5) =  1;   A(55,6) =  0;
        A(56,0) =  0;   A(56,1) =  0;   A(56,2) =  0;   A(56,3) =  0;   A(56,4) =  0;   A(56,5) =  0;   A(56,6) =  0;
        A(57,0) =  0;   A(57,1) =  0;   A(57,2) =  0;   A(57,3) =  0;   A(57,4) =  0;   A(57,5) =  0;   A(57,6) =  0;
        A(58,0) =  0;   A(58,1) =  0;   A(58,2) =  0;   A(58,3) =  0;   A(58,4) =  0;   A(58,5) =  0;   A(58,6) =  0;
        A(59,0) =  0;   A(59,1) =  0;   A(59,2) =  0;   A(59,3) =  0;   A(59,4) =  0;   A(59,5) =  0;   A(59,6) =  0;
        A(60,0) =  0;   A(60,1) =  0;   A(60,2) =  0;   A(60,3) =  0;   A(60,4) =  0;   A(60,5) =  0;   A(60,6) =  0;
        A(61,0) =  0;   A(61,1) =  0;   A(61,2) =  0;   A(61,3) =  0;   A(61,4) =  0;   A(61,5) =  0;   A(61,6) =  0;
        break;
    default:
        return status(HW_MATH_ERR_INVALIDINPUT, 1);
    };

    return status;
}
