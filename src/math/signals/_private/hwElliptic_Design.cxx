/**
* @file  hwElliptic_Design.cxx
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
#include "hwElliptic_Design.h"

#include "EllipticFuncs.h"

#ifndef OS_WIN
    #include <math.h>
#endif

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwElliptic_Design::hwElliptic_Design()
{
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwElliptic_Design::~hwElliptic_Design()
{
}
//------------------------------------------------------------------------------
// Design a low pass prototype filter
//------------------------------------------------------------------------------
void hwElliptic_Design::DesignPrototype(double  epsilon, 
                                        double  delta, 
                                        double  omegaS2,
                                        int&    order, 
                                        double& omegaP1, 
                                        double& omegaP2)
{
    // See James A. Moorer, "The Manifold Joys of Conformal Mapping"
    // design a low pass prototye filter with normalized corner frequencies
    // compute order and range of possible normalized omegaP values
    // omegaP1 = default pass band corner = 1.0
    // omegaP2 = alternate pass band corner > omegaP1
    // omegaS1 = default stop band corner < omegaS2, not computed
    // omegaS2 = alternate stop band corner, the input upper bound
    // omegaPS = omegaP2 / omegaS2 = omegaP1 / omegaS1, adjusted for integer filter size

    // Note: ellipk() used here is not the standard form. K(x) = ellpk(1-x^2).
    //       Likewise for ellipe()
    double value1 = epsilon / delta;
    double value2 = 1.0 / omegaS2;

    value1 *= value1;
    value2 *= value2;

    double value3 = ellpk(1.0 - value2);
    double value4 = ellpk(value1);
    double value5 = ellpk(1.0 - value1);
    double value6 = ellpk(value2);
    double orderFp = (value3 * value4) / (value5 * value6); // floating point filter order
    order = static_cast<int>(ceil(orderFp - 1.0e-12));
    double nome_hat = exp(-PI * value4 / (value5 * static_cast<double>(order)));
    double omegaPS = cay(nome_hat);
    omegaP1 = 1.0;               // hit pass band requirement
    omegaP2 = omegaS2 * omegaPS; // hit stop band requirement
}
