/**
* @file EllipticFuncs.h
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
#ifndef _MathUtils_EllipticFuncs_h
#define _MathUtils_EllipticFuncs_h

#include <MathUtilsExports.h>

//------------------------------------------------------------------------------
//!
//! \brief Elliptic functions 
//!
//------------------------------------------------------------------------------

//!
//! Returns the complete elliptic integral of the first kind
//! \param x Input
//!
MATHUTILS_DECLS double ellpk(double x);
//!
//! Returns the complete elliptic integral of the second kind
//! \param x Input
//!
MATHUTILS_DECLS double ellpe(double x);
//!
//! Returns the incomplete elliptic integral of the first kind
//! \param phi Given amplitude
//! \param m   Given modulus
//!
MATHUTILS_DECLS double ellik(double phi, 
                             double m);
//!
//! Returns 0 if successful in evaluating the Jacobian elliptic functions 
//! sn(u|m), cn(u|m), and dn(u|m) of parameter m between 0 and 1, and real 
//! argument u. Returns -1 in case of failure.
//! \param u  Real argument
//! \param m  Input parameter
//! \param sn Output
//! \param cn Output
//! \param dn Output
//! \param ph Output
//!
MATHUTILS_DECLS int ellpj (double u, 
                           double m,
                           double&       sn, 
                           double&       cn, 
                           double&       dn, 
                           double&       ph);
//!
//! Finds parameter corresponding to given nome by expansion in theta functions
//! \param q Input
//!
MATHUTILS_DECLS double cay(double q);

#endif // _MathUtils_EllipticFuncs_h
