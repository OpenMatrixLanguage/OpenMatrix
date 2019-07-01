/**
* @file hwGaussQuadrature.cxx
* @date June 2007
* Copyright (C) 2007-2019 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language ("OpenMatrix") software.
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

#include "hwGaussQuadrature.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwGaussQuadrature::hwGaussQuadrature()
{
    n = 0;
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwGaussQuadrature::~hwGaussQuadrature()
{
}
//------------------------------------------------------------------------------
// Returns status and gets area after integrating pFunc(z) from a to b
//------------------------------------------------------------------------------
hwMathStatus hwGaussQuadrature::Compute(const QuadFunc1 pFunc,
                                        double          a,
                                        double          b, 
                                        double&         area)
{
    if (!pFunc)
    {
        return hwMathStatus(HW_MATH_ERR_NULLPOINTER, 1);
    }

    double c = 0.5 * (b + a);
    double d = 0.5 * (b - a);
    hwMatrix z = X * d + c;
    hwMatrix fz;

    hwMathStatus status = pFunc(z, fz);

    if (!status.IsOk())
    {
        return status;
    }

    area = 0.0;

    for (int i = 0; i < n; ++i)
    {
        area += W(i) * fz(i);
    }

    area *= d;

    return status;
}
//------------------------------------------------------------------------------
// Returns status and gets area after integrating from a to b
//------------------------------------------------------------------------------
hwMathStatus hwGaussQuadrature::ComputeRLog(const QuadFunc1 pFunc, 
                                            double          a,
                                            double          b, 
                                            double&         area)
{
    // apply transform
    bool rightSide;

    if (a > 0.0 && b > 0.0)
    {
        rightSide = true;
        a = 1.0 / log(a+1);
        b = 1.0 / log(b+1);
    }
    else if (a < 0.0 && b < 0.0)
    {
        rightSide = false;
        a = 1.0 / log(-a+1);
        b = 1.0 / log(-b+1);
    }
    else
    {
        return hwMathStatus(HW_MATH_ERR_INVALIDINPUT, 2, 3);
    }

    double c = 0.5 * (b + a);
    double d = 0.5 * (b - a);
    hwMatrix z = X * d + c;
    hwMatrix u(X.M(), X.N(), hwMatrix::REAL);
    hwMatrix fz;
    hwMathStatus status;

    for (int i = 0; i < n; ++i)
        u(i) = exp(1.0 / z(i));

    if (rightSide)
        status = pFunc(u - 1.0, fz);
    else
        status = pFunc(1.0 - u, fz);

    if (!status.IsOk())
        return status;

    area = 0.0;

    for (int i = 0; i < n; ++i)
    {
        if (fz(i))
        {
            area += W(i) * fz(i) * u(i) * 1.0/(z(i)*z(i));
        }
    }

    if (rightSide)
        area = -area;

    area *= d;

    return status;
}
//------------------------------------------------------------------------------
// Returns status and gets area after integrating from a to b
//------------------------------------------------------------------------------
hwMathStatus hwGaussQuadrature::ComputeSqrt1(const QuadFunc1 pFunc, 
                                             double          a,
                                             double          b, 
                                             double&         area)
{
    // apply transform
    double aa = 0.0;
    double bb = sqrt(b-a);
    double c = 0.5 * (bb + aa);
    double d = 0.5 * (bb - aa);
    hwMatrix z = X * d + c;
    hwMatrix u;
    hwMatrix fz;
    hwMathStatus status;

    status = u.MultByElems(z, z);
    status = pFunc(u + a, fz);

    if (!status.IsOk())
        return status;

    area = 0.0;

    for (int i = 0; i < n; ++i)
        area += W(i) * fz(i) * (2.0*z(i));

    area *= d;

    return status;
}
//------------------------------------------------------------------------------
//! Returns status and gets area after integrating from a to b
//------------------------------------------------------------------------------
hwMathStatus hwGaussQuadrature::ComputeSqrt2(const QuadFunc1 pFunc, 
                                             double          a,
                                             double          b, 
                                             double&         area)
{
    // apply transform
    double aa = sqrt(b-a);
    double bb = 0.0;
    double c = 0.5 * (bb + aa);
    double d = 0.5 * (bb - aa);
    hwMatrix z = X * d + c;
    hwMatrix u;
    hwMatrix fz;
    hwMathStatus status;

    status = u.MultByElems(z, z);
    status = pFunc(b - u, fz);

    if (!status.IsOk())
        return status;

    area = 0.0;

    for (int i = 0; i < n; ++i)
        area -= W(i) * fz(i) * (2.0*z(i));

    area *= d;

    return status;
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwGaussLegendre::hwGaussLegendre(int n_)
{
    SetPnts(n_);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwGaussLegendre::~hwGaussLegendre()
{
}
//------------------------------------------------------------------------------
// Sets points
//------------------------------------------------------------------------------
void hwGaussLegendre::SetPnts(int n_)
{
    if (n == n_)
    {
        return;
    }

    if (n_ > 0 && n_ < 11)
    {
        n = n_;
    }
    else
    {
        m_status(HW_MATH_ERR_GAUSSLEGENDRE, 1);
        n = 0;
    }

    hwMathStatus status = X.Dimension(n, hwMatrix::REAL);
    if (!status.IsOk())
    {
        n = 0;
        m_status = status;
        return;
    }

    status = W.Dimension(n, hwMatrix::REAL);
    if (!status.IsOk())
    {
        n = 0;
        m_status = status;
        return;
    }

    GetLocations();
    GetWeights();
}
//------------------------------------------------------------------------------
// Populates locations
//------------------------------------------------------------------------------
void hwGaussLegendre::GetLocations()
{
    switch(n)
    {
        case 1:
            X(0) = 0.000000000000000;
            break;
        case 2:
            X(0) = -0.577350269189626;  
            X(1) =  0.577350269189626;
            break;
        case 3:
            X(0) = -0.774596669241483;
            X(1) =  0.000000000000000;
            X(2) =  0.774596669241483;
            break;
        case 4:
            X(0) = -0.861136311594053;
            X(1) = -0.339981043584856;
            X(2) =  0.339981043584856;
            X(3) =  0.861136311594053;
            break;
        case 5:
            X(0) = -0.906179845938664;
            X(1) = -0.538469310105683;
            X(2) =  0.000000000000000;
            X(3) =  0.538469310105683;
            X(4) =  0.906179845938664;
            break;
        case 6:
            X(0) = -0.932469514203152;
            X(1) = -0.661209386466265;
            X(2) = -0.238619186083197;
            X(3) =  0.238619186083197;
            X(4) =  0.661209386466265;
            X(5) =  0.932469514203152;
            break;
        case 7:
            X(0) = -0.949107912342759;
            X(1) = -0.741531185599394;
            X(2) = -0.405845151377397;
            X(3) =  0.000000000000000;
            X(4) =  0.405845151377397;
            X(5) =  0.741531185599394;
            X(6) =  0.949107912342759;
            break;
        case 8:
            X(0) = -0.960289856497536;
            X(1) = -0.796666477413627;
            X(2) = -0.525532409916329;
            X(3) = -0.183434642495650;
            X(4) =  0.183434642495650;
            X(5) =  0.525532409916329;
            X(6) =  0.796666477413627;
            X(7) =  0.960289856497536;
            break;
        case 9:
            X(0) = -0.968160239507626;
            X(1)= -0.836031107326636;
            X(2)= -0.613371432700590;
            X(3)= -0.324253423403809;
            X(4)=  0.000000000000000;
            X(5)=  0.324253423403809;
            X(6)=  0.613371432700590;
            X(7)=  0.836031107326636;
            X(8)=  0.968160239507626;
            break;
        case 10:
            X(0)= -0.973906528517172;
            X(1)= -0.865063366688985;
            X(2)= -0.679409568299024;
            X(3)= -0.433395394129247;
            X(4)= -0.148874338981631;
            X(5)=  0.148874338981631;
            X(6)=  0.433395394129247;
            X(7)=  0.679409568299024;
            X(8)=  0.865063366688985;
            X(9)=  0.973906528517172;
            break;
        default: break;
    }
}
//------------------------------------------------------------------------------
// Populates weights
//------------------------------------------------------------------------------
void hwGaussLegendre::GetWeights()
{
    switch(n)
    {
        case 1:
            W(0) = 2.000000000000000;
            break;
        case 2:
            W(0) = 1.000000000000000;  
            W(1) = 1.000000000000000;
            break;
        case 3:
            W(0) = 0.555555555555556;
            W(1) = 0.888888888888889;
            W(2) = 0.555555555555556;
            break;
        case 4:
            W(0) = 0.347854845137454;
            W(1) = 0.652145154862546;
            W(2) = 0.652145154862546;
            W(3) = 0.347854845137454;
            break;
        case 5:
            W(0) = 0.236926885056189;
            W(1) = 0.478628670499366;
            W(2) = 0.568888888888889;
            W(3) = 0.478628670499366;
            W(4) = 0.236926885056189;
            break;
        case 6:
            W(0) = 0.171324492379170;
            W(1) = 0.360761573048139;
            W(2) = 0.467913934572691;
            W(3) = 0.467913934572691;
            W(4) = 0.360761573048139;
            W(5) = 0.171324492379170;
            break;
        case 7:
            W(0) = 0.129484966168870;
            W(1) = 0.279705391489277;
            W(2) = 0.381830050505119;
            W(3) = 0.417959183673469;
            W(4) = 0.381830050505119;
            W(5) = 0.279705391489277;
            W(6) = 0.129484966168870;
            break;
        case 8:
            W(0) = 0.101228536290376;
            W(1) = 0.222381034453374;
            W(2) = 0.313706645877887;
            W(3) = 0.362683783378362;
            W(4) = 0.362683783378362;
            W(5) = 0.313706645877887;
            W(6) = 0.222381034453374;
            W(7) = 0.101228536290376;
            break;
        case 9:
            W(0) = 0.081274388361574;
            W(1) = 0.180648160694857;
            W(2) = 0.260610696402935;
            W(3) = 0.312347077040003;
            W(4) = 0.330239355001260;
            W(5) = 0.312347077040003;
            W(6) = 0.260610696402935;
            W(7) = 0.180648160694857;
            W(8) = 0.081274388361574;
            break;
        case 10:
            W(0) = 0.066671344308688;
            W(1) = 0.149451349150581;
            W(2) = 0.219086362515982;
            W(3) = 0.269266719309996;
            W(4) = 0.295524224714753;
            W(5) = 0.295524224714753;
            W(6) = 0.269266719309996;
            W(7) = 0.219086362515982;
            W(8) = 0.149451349150581;
            W(9) = 0.066671344308688;
            break;
        default: break;
    }
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwGaussLobatto::hwGaussLobatto(int n_)
{
    SetPnts(n_);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwGaussLobatto::~hwGaussLobatto()
{
}
//------------------------------------------------------------------------------
// Sets points
//------------------------------------------------------------------------------
void hwGaussLobatto::SetPnts(int n_)
{
    if (n == n_)
    {
        return;
    }

    if (n_ > 1 && n_ < 11)
    {
        n = n_;
    }
    else
    {
        m_status(HW_MATH_ERR_GAUSSLOBATTO, 1);
        n = 0;
    }

    hwMathStatus status = X.Dimension(n, hwMatrix::REAL);
    if (!status.IsOk())
    {
        n = 0;
        m_status = status;
        return;
    }

    status = W.Dimension(n, hwMatrix::REAL);
    if (!status.IsOk())
    {
        n = 0;
        m_status = status;
        return;
    }

    GetLocations();
    GetWeights();
}
//------------------------------------------------------------------------------
// Populates locations
//------------------------------------------------------------------------------
void hwGaussLobatto::GetLocations()
{
    switch(n)
    {
        case 2:
            X(0) = -1.0;  
            X(1) =  1.0;
            break;
        case 3:
            X(0) = -1.0;
            X(1) =  0.0;
            X(2) =  1.0;
            break;
        case 4:
            X(0) = -1.0;
            X(1) = -sqrt(0.2);
            X(2) =  -X(1);
            X(3) =  1.0;
            break;
        case 5:
            X(0) = -1.0;
            X(1) = -sqrt(3.0 / 7.0);
            X(2) =  0.000000000000000;
            X(3) = -X(1);
            X(4) =  1.0;
            break;
        case 6:
            X(0) = -1.0;
            X(1) = -0.76505532392946;
            X(2) = -0.28523151648065;
            X(3) =  0.28523151648065;
            X(4) =  0.76505532392946;
            X(5) =  1.0;
            break;
        case 7:
            X(0) = -1.0;
            X(1) = -0.83022389627857;
            X(2) = -0.46884879347071;
            X(3) =  0.0;
            X(4) =  0.46884879347071;
            X(5) =  0.83022389627857;
            X(6) =  1.0;
            break;
        case 8:
            X(0) = -1.0;
            X(1) = -0.871740148509607;
            X(2) = -0.591700181433142;
            X(3) = -0.209299217902479;
            X(4) =  0.209299217902479;
            X(5) =  0.591700181433142;
            X(6) =  0.871740148509607;
            X(7) =  1.0;
            break;
        case 9:
            X(0) = -1.0;
            X(1) = -0.899757995411460;
            X(2) = -0.677186279510738;
            X(3) = -0.363117463826178;
            X(4) =  0.0;
            X(5) =  0.363117463826178;
            X(6) =  0.677186279510738;
            X(7) =  0.899757995411460;
            X(8) =  1.0;
            break;
        case 10:
            X(0) = -1.0;
            X(1) = -0.919533908166459;
            X(2) = -0.738773865105505;
            X(3) = -0.477924949810444;
            X(4) = -0.165278957666387;
            X(5) =  0.165278957666387;
            X(6) =  0.477924949810444;
            X(7) =  0.738773865105505;
            X(8) =  0.919533908166459;
            X(9) =  1.0;
            break;
        default: break;
    }
}
//------------------------------------------------------------------------------
// Populates weights
//------------------------------------------------------------------------------
void hwGaussLobatto::GetWeights()
{
    switch(n)
    {
        case 2:
            W(0) = 1.000000000000000;  
            W(1) = 1.000000000000000;
            break;
        case 3:
            W(0) = 1.0 / 3.0;
            W(1) = 4.0 / 3.0;
            W(2) = 1.0 / 3.0;
            break;
        case 4:
            W(0) = 0.1666666666666667;
            W(1) = 0.8333333333333333;
            W(2) = 0.8333333333333333;
            W(3) = 0.1666666666666667;
            break;
        case 5:
            W(0) = 0.1;
            W(1) = 49.0 / 90.0;
            W(2) = 32.0 / 45.0;
            W(3) = 49.0 / 90.0;
            W(4) = 0.1;
            break;
        case 6:
            W(0) = 0.066666666666667;
            W(1) = 0.378474956297847;
            W(2) = 0.554858377035486;
            W(3) = 0.554858377035486;
            W(4) = 0.378474956297847;
            W(5) = 0.066666666666667;
            break;
        case 7:
            W(0) = 0.047619047619048;
            W(1) = 0.276826047361566;
            W(2) = 0.431745381209863;
            W(3) = 0.487619047619048;
            W(4) = 0.431745381209863;
            W(5) = 0.276826047361566;
            W(6) = 0.047619047619048;
            break;
        case 8:
            W(0) = 0.035714285714286;
            W(1) = 0.210704227143506;
            W(2) = 0.341122692483504;
            W(3) = 0.412458794658704;
            W(4) = 0.412458794658704;
            W(5) = 0.341122692483504;
            W(6) = 0.210704227143506;
            W(7) = 0.035714285714286;
            break;
        case 9:
            W(0) = 0.027777777777778;
            W(1) = 0.165495361560806;
            W(2) = 0.274538712500162;
            W(3) = 0.346428510973043;
            W(4) = 0.371519274376417;
            W(5) = 0.346428510973043;
            W(6) = 0.274538712500162;
            W(7) = 0.165495361560806;
            W(8) = 0.027777777777778;
            break;
        case 10:
            W(0)= 0.0222222222222222;
            W(1)= 0.1333059908510701;
            W(2)= 0.2248893420631265;
            W(3)= 0.2920426836796838;
            W(4)= 0.3275397611838975;
            W(5)= 0.3275397611838975;
            W(6)= 0.2920426836796838;
            W(7)= 0.2248893420631265;
            W(8)= 0.1333059908510701;
            W(9)= 0.0222222222222222;
            break;
        default: break;
    }
}
