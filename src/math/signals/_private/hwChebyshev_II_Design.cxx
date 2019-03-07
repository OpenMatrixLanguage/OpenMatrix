/**
* @file hwChebyshev_II_Design.h
* @date April 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
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

#include <cmath>

#include "hwChebyshev_II_Design.h"


//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwChebyshev_II_Design::hwChebyshev_II_Design()
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwChebyshev_II_Design::~hwChebyshev_II_Design()
{
}
//------------------------------------------------------------------------------
// Design a low pass prototype filter
//------------------------------------------------------------------------------
void hwChebyshev_II_Design::DesignPrototype(double  epsilon, 
                                            double  delta,
                                            double  omegaS, 
                                            int&    order,
                                            double& omegaS1, 
                                            double& omegaS2)
{
    // design a low pass prototye filter with normalized corner frequencies
    // compute order and range of possible normalized omegaS values
    // omegaP1 = default pass band corner > omegaP2, not computed
    // omegaP2 = alternate pass band corner = 1.0
    // omegaS1 = default stop band corner, the input upper bound = omegaS
    // omegaS2 = alternate stop band corner < omegaS1
    // omegaSP = omegaS2 / omegaP2 = omegaS1 / omegaP1, adjusted for integer filter size
    double value = delta / epsilon;
    double orderFp = acosh(value) / acosh(omegaS);  // floating point filter order
    order = static_cast<int>(ceil(orderFp - 1.0e-12));
    double omegaSP = cosh(acosh(value) / static_cast<double>(order));
    omegaS1 = omegaS;  // hit stop band requirement
    omegaS2 = omegaSP; // hit pass band requirement
}
