/**
* @file GeometryFuncs.h
* @date December, 2018
* Copyright (C) 2017-2018 Altair Engineering, Inc.
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

#ifndef _Geometry_Funcs_h
#define _Geometry_Funcs_h

#include <string>
#include "GeometryExports.h"

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;
typedef hwTMatrix<int, hwTComplex<int> > hwMatrixI;

//------------------------------------------------------------------------------
//!
//! \brief Geometry functions
//!
//------------------------------------------------------------------------------

//! 
//! Computes 2D convex hull and returns status
//! \param x vector
//! \param y vector
//! \param options QHull options
//! \param hull convex hull vertices vector
//! \param area convex hull area
//! \param errfile file for writing QHull error information
//!
GEOMETRY_DECLS hwMathStatus ConvexHull(const hwMatrix&    x,
                                       const hwMatrix&    y,
                                       const std::string& options,
                                       hwMatrixI&         hull,
                                       double&            area,
                                       FILE*              errfile);

//! 
//! Computes ND convex hull and returns status
//! \param P points matrix
//! \param options QHull options
//! \param hull convex hull vertices matrix
//! \param convex hull volume
//! \param errfile file for writing QHull error information
//!
GEOMETRY_DECLS hwMathStatus ConvexHulln(const hwMatrix&    P,
                                        const std::string& options,
                                        hwMatrixI&         hull,
                                        double&            volume,
                                        FILE*              errfile);
//! 
//! Computes ND dealaunay tringulation and returns status
//! \param options QHull options
//! \param P points matrix
//! \param triang matrix of triangulation points
//! \param errfile file for writing QHull error information
//!
GEOMETRY_DECLS hwMathStatus Delaunayn(const hwMatrix&    P,
                                      const std::string& options,
                                      hwMatrixI&         triang,
                                      FILE*              errfile);

#endif // _Geometry_Funcs_h