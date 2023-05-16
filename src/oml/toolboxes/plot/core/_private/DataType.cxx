/**
* @file DataType.cxx
* @date May 2018
* Copyright (C) 2018-2023 Altair Engineering, Inc.  
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

#include "DataType.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include "OML_Error.h"

namespace omlplot{

    LineData::LineData()
        :index(0), xcolcount(1), ycolcount(1), zcolcount(1), bar3PerRow(false)
    {
    }

    LineData::~LineData(){
    }

    GetData::GetData(){
    }

    GetData::~GetData(){
    }

    SetData::SetData(){
    }

    SetData::~SetData(){
    }

    Color::Color(){
        m_vComponent.push_back(0);
        m_vComponent.push_back(0);
        m_vComponent.push_back(0);
    }

    Color::Color(Currency c){
        if (c.IsString()){
            init(c.StringVal());
        } else if (c.IsMatrix()){
            init(c.Matrix());
        }
    }

    void Color::init(std::string name){
        m_vComponent.resize(3);
        if (name == "blue" || name == "b"){
            m_vComponent[0] = 0;
            m_vComponent[1] = 0;
            m_vComponent[2] = 255;
        } else if (name == "red" || name == "r"){
            m_vComponent[0] = 255;
            m_vComponent[1] = 0;
            m_vComponent[2] = 0;
        } else if (name == "green" || name == "g"){
            m_vComponent[0] = 0;
            m_vComponent[1] = 255;
            m_vComponent[2] = 0;
        } else if (name == "cyan" || name == "c"){
            m_vComponent[0] = 0;
            m_vComponent[1] = 255;
            m_vComponent[2] = 255;
        } else if (name == "yellow" || name == "y"){
            m_vComponent[0] = 255;
            m_vComponent[1] = 255;
            m_vComponent[2] = 0;
        } else if (name == "magenta" || name == "m"){
            m_vComponent[0] = 255;
            m_vComponent[1] = 0;
            m_vComponent[2] = 255;
        } else if (name == "white" || name == "w"){
            m_vComponent[0] = 255;
            m_vComponent[1] = 255;
            m_vComponent[2] = 255;
        } else if (name == "black" || name == "k"){
            m_vComponent[0] = 0;
            m_vComponent[1] = 0;
            m_vComponent[2] = 0;
        } else {
            throw OML_Error(OML_ERR_PLOT_INVALID_COLOR);
        }
    }

    void Color::init(const hwMatrix *mat){
        if (! (mat->IsVector() && mat->Size() == 3)){
            throw OML_Error(OML_ERR_PLOT_INVALID_COLOR);
        }
        double r = (*mat)(0), g = (*mat)(1), b = (*mat)(2);
        if (r <=1 && g <= 1 && b <= 1){
            r *= 255;
            g *= 255;
            b *= 255;
        }
        m_vComponent.push_back(int(r));
        m_vComponent.push_back(int(g));
        m_vComponent.push_back(int(b));
    }

    Color::Color(vector<int> c){
        if (c.size() != 3){
            throw;
        }
        if (!( ((0 <= c[0]) && (c[0] <= 255)) &&
               ((0 <= c[1]) && (c[1] <= 255)) &&
               ((0 <= c[2]) && (c[2] <= 255)) )){
            throw;
        }
        m_vComponent = c;
    }

    vector<double> Color::getComponent(){
        vector<double> res;
        res.push_back(m_vComponent[0]);
        res.push_back(m_vComponent[1]);
        res.push_back(m_vComponent[2]);
        return res;
    }

    std::string Color::getString(){
        stringstream ss;
        ss << std::hex << '#' 
           << std::setw(2) << std::setfill('0') << m_vComponent[0]
           << std::setw(2) << std::setfill('0') << m_vComponent[1]
           << std::setw(2) << std::setfill('0') << m_vComponent[2];
        return ss.str();
    }
} // namespace omlplot
