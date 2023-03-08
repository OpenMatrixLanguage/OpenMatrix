/**
* @file DataType.h
* @date May 2017
* Copyright (C) 2017-2023 Altair Engineering, Inc.  
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

#ifndef _PLOT_TYPE_H_
#define _PLOT_TYPE_H_

#include "OmlPlotExport.h"
#include <memory>
#include <vector>
#include <string>
#include "Currency.h"
#include "hwMatrix.h"

using namespace std;

namespace omlplot{

    // All data type should independent of Currency or hwMatrix,
    // so save data as vectores.
    struct OMLPLOT_EXPORT LineData{
        LineData();
        ~LineData();

        double parent;
        vector<double> x;
        vector<double> y;
        vector<double> z;
        vector<double> u;
        vector<double> v;
        vector<double> colorMat;
        string style;
        vector<string> properties;
        vector<Currency> values;
        // TODO: Currency contourLevel;
        double areaLevel;
        int index;
        int xcolcount;
        int ycolcount;
        int zcolcount;
        Currency xCategories;
        bool bar3PerRow;
    };

    struct OMLPLOT_EXPORT SurfData{
        double parent;
        vector<double> x;
        vector<double> y;
        vector<double> z;
        vector<string> properties;
        vector<Currency> values;
    };

    struct OMLPLOT_EXPORT GetData{
        GetData();
        ~GetData();
        vector<double> handles;
        vector<string> properties;
    };

    struct OMLPLOT_EXPORT SetData{
        SetData();
        ~SetData();
        vector<double> handles;
        vector<string> properties;
        vector<Currency> values;
    };

    struct OMLPLOT_EXPORT FigureData{
        vector<double> fids;
        vector<string> properties;
        vector<Currency> values;
    };

    struct OMLPLOT_EXPORT AxesData{
        vector<double> handles;
        vector<string> properties;
        vector<Currency> values;
    };

    struct OMLPLOT_EXPORT LegendData{
        std::vector<std::string> legends;
        std::vector<double> handles;
        vector<string> properties;
        vector<Currency> values;
    };

    struct OMLPLOT_EXPORT ColorbarData {
        ColorbarData()
            : toggleVisibility(false), visible(false)
        {

        }
        ~ColorbarData() {}
        std::vector<string> properties;
        std::vector<Currency> values;
        bool toggleVisibility;
        bool visible;
    };

    struct OMLPLOT_EXPORT TextData{
        std::vector<double> xpos;
        std::vector<double> ypos;
        std::vector<double> zpos;
        std::vector<std::string> text;
        vector<string> properties;
        vector<Currency> values;
    };

    struct OMLPLOT_EXPORT LimData{
        LimData() : handle(0){}
        double handle;
        std::vector<double> limits;
    };

    class OMLPLOT_EXPORT Color{
    public:
        Color();
        Color(Currency c);
        Color(vector<int>);

        void init(std::string);
        void init(const hwMatrix *mat);
        vector<double> getComponent();
        std::string getString();
    private:
        vector<int> m_vComponent;
    };

    struct OMLPLOT_EXPORT QueryData {
        enum SEARCH_OP { OP_INV, OP_AND, OP_OR, OP_XOR, OP_NOT };
        QueryData()
            :m_onlyProperty(false), m_depth(-1)
        {
        }
        ~QueryData() {}

        std::vector<double> handles;
        std::vector<std::string> properties;
        std::vector<Currency> values;
        std::vector<int> ops;
        bool m_onlyProperty;
        int m_depth;
    };

}	  // namespace omlplot

#endif	// _PLOT_TYPE_H_
