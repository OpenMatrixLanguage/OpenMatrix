/**
* @file hwButterworth_Design.cxx
* @date April 2009
* Copyright (C) 2009-2018 Altair Engineering, Inc.  
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
#include "hwButterworth_Design.h"

#include <math.h>

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwButterworth_Design::hwButterworth_Design()
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwButterworth_Design::~hwButterworth_Design()
{
}
//------------------------------------------------------------------------------
// Design a low pass prototype filter with normalized -3db cutoff frequencies
//------------------------------------------------------------------------------
void hwButterworth_Design::DesignPrototype(double  epsilon, 
                                           double  delta,
                                           double  omegaSP,
                                           int&    order,
                                           double& omegaC1, 
                                           double& omegaC2)
{
    // compute order and range of possible normalized omegaC values
    // omegaP1 = default pass band corner = 1.0
    // omegaP2 = alternate pass band corner > omegaP1
    // omegaS1 = default stop band corner < omegaS2, not computed
    // omegaS2 = alternate stop band corner, the upper bound
    // omegaSP = omegaS2 / omegaP2 = omegaS1 / omegaP1, adjusted for integer filter size
    double orderFp = log(delta/epsilon) / log(omegaSP); // floating point filter order
    order = (int) ceil(orderFp - 1.0e-12);
    omegaC1 = exp(-log(epsilon) / static_cast<double>(order));         // hit pass band requirement
    omegaC2 = omegaSP * exp(-log(delta) / static_cast<double>(order)); // hit stop band requirement
}
