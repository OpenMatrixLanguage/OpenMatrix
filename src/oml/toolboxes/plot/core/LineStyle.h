/**
* @file LineStyle.h
* @date May 2017
* Copyright (C) 2017-2018 Altair Engineering, Inc.  
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

#ifndef _LINESTYLE_H_
#define _LINESTYLE_H_

#include "OmlPlotExport.h"
#include "DataType.h"
#include <string>
#include <vector>

using namespace std;

namespace omlplot{
    class OMLPLOT_EXPORT LineStyle{
    public:
        explicit LineStyle(const LineData&);

    private:
        void SplitFormat(const string fmt, string &style, string &legend);
        bool isAbbrColor(const char &c);
        bool isAbbrLine(const char &l);
        bool isAbbrMarker(const char &m);

    public:
        string m_lineStyle;
        Color m_lineColor;
        int m_lineWidth;
        string m_markerStyle;
        Color m_markerColor;
        int m_markerSize;
        string m_legend;

        static bool helper;
        static bool initHelper();
        static std::vector<int> c1;
        static std::vector<int> c2;
        static std::vector<int> c3;
        static std::vector<int> c4;
        static std::vector<int> c5;
    };
}

#endif
