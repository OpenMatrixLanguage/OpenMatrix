/**
* @file SpecialFuncs.h
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
#ifndef _MathUtils_SpecialFuncs_h
#define _MathUtils_SpecialFuncs_h

#include <MathUtilsExports.h>

//------------------------------------------------------------------------------
//!
//! \brief Contains mathematical utility functions
//!
//------------------------------------------------------------------------------

// Gamma related
//! 
//! Returns the gamma function of the given argument
//! \param x Input argument
//!
MATHUTILS_DECLS double GammaFunc(double x);
//! 
//! Returns the base e (2.718...) logarithm of the absolute value of the gamma 
//! function of the input argument
//! \param x Input argument
//!
MATHUTILS_DECLS double GammaLog(double x);
//!
//! Returns the Psi (digamma) value which is the logarithmic derivative of the 
//! gamma function
//! \param x Input argument
//!
MATHUTILS_DECLS double Digamma(double x);
//! 
//! Trigamma function
//! \param x Input argument
//!
MATHUTILS_DECLS double Trigamma(double x);

// Beta related
//! 
//! Beta function
//! \param x Input argument
//! \param y Input argument
//!
MATHUTILS_DECLS double BetaFunc(double x, 
                                double y);
//! 
//! Log-Beta function
//! \param x Input argument
//! \param y Input argument
//!
MATHUTILS_DECLS double BetaLog(double x, 
                               double y);

// Bessel
//! 
//! Zero order modified Bessel function of the first kind (Io)
//! \param alpha Input argument
//!
MATHUTILS_DECLS double Bessel(double alpha);

#endif // _MathUtils_SpecialFuncs_h

