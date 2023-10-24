/**
* @file DataStateMachine.h
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

#ifndef _DATA_STATE_MACHINE_H_
#define _DATA_STATE_MACHINE_H_

#include "DataType.h"
#include "CoreMain.h"
#include <Currency.h>
#include <memory>
#include <algorithm>
#include "OML_Error.h"
#include "hwMatrix.h"
#include "hwMatrixN.h"

namespace omlplot{

const std::string DEFAULT_STYLE = ""; // blue, point, solid line

class OMLPLOT_EXPORT DataStateMachine
{
public:
    DataStateMachine();
    ~DataStateMachine();

    std::vector<LineData> getLineData(const std::vector<Currency> &inputs);
    std::vector<LineData> getBarData(const std::vector<Currency> &inputs);
    std::vector<LineData> getFillData(const std::vector<Currency> &inputs);
    std::vector<LineData> getSurfData(const std::vector<Currency> &inputs);
    std::vector<LineData> getHistData(const std::vector<Currency>& inputs);
    std::vector<LineData> getXLineData(const std::vector<Currency> &inputs);
    std::unique_ptr<GetData> getGetData(const std::vector<Currency> &inputs);
    std::unique_ptr<SetData> getSetData(const std::vector<Currency> &inputs);
    std::unique_ptr<FigureData> getFigureData(const std::vector<Currency> &inputs);
    std::unique_ptr<AxesData> getAxesData(const std::vector<Currency> &inputs);
    std::unique_ptr<LegendData> getLegendData(const std::vector<Currency> &inputs);
    std::unique_ptr<TextData> getTextData(const std::vector<Currency> &inputs);
    std::unique_ptr<LimData> getLimData(const std::vector<Currency> &inputs);
    std::unique_ptr<QueryData> getQueryData(const std::vector<Currency>& inputs);
    std::vector<LineData> getShapeData(const std::vector<Currency>& inputs, const std::string& shape);
    std::vector<LineData> getPatchData(const std::vector<Currency>& inputs);
    std::unique_ptr<ColorbarData> getColorbarData(const std::vector<Currency>& inputs);
    std::vector<LineData> getQuiverData(const std::vector<Currency>& inputs);
    std::vector<LineData> getBar3Data(const std::vector<Currency>& inputs);
    std::vector<LineData> getHist3Data(const std::vector<Currency>& inputs, std::vector<Currency>& hist3dData);
    std::vector<LineData> getTriplotData(const std::vector<Currency>& inputs);
    std::vector<LineData> getTrisurfData(const std::vector<Currency>& inputs);

private:
    enum State { START, GET_X, GET_Y, GET_Z, // general state
                 GET_PROP_VALUE, GET_FMT, STATE_ERROR, // general state
                 GET_X_SCALAR, GET_X_VECTOR, GET_WIDTH, // state in bar
                 GET_LEVEL,     // state in area
                 GET_COMPLEX,   // state in polar
                 GET_COLOR_MAT, // state in surf
                 GET_MARKER_SIZE, GET_MARKER_COLOR, // state in scatter
                 GET_MARKER_STYLE, GET_MARKER_FILLED, // state in scatter
                 GET_CONTOUR_LEVEL, // state in contour
                 GET_FID,           // state in figure
                 GET_HANDLE,        // state in set
                 GET_LIM,           // state in lim
                 GET_LEGEND,        // state in legend
                 GET_PATH, GET_MAT, // state in image
                 GET_TEXT,
                 GET_SINGLE_PROP,
                 GET_SEARCH_DEPTH,
                 GET_U, GET_V,
                 DONE
    };
    enum GetDir{
        BY_ROW, BY_COL
    };

private:
    CoreMain *cm;
    double _temp_parent;

    bool isColor(const string &s);
    bool isMarker(const string &s);
	Currency getNextInput(const std::vector<Currency> &inputs, const int pos);
    std::shared_ptr<hwMatrix> GetMatrix(const Currency &c);
    vector<double> MatrixToVector(const hwMatrix *const);
    std::shared_ptr<hwMatrix> GetMatrixIndex(std::shared_ptr<hwMatrix> mat, GetDir dir);
    void GetHistData(const hwMatrix& data, const int numBins, hwMatrix &bin, hwMatrix &freq);
    void ExtractData(const hwMatrix * const x, const hwMatrix * const y,
                     vector<LineData> &result /*out*/,
                     int &count               /*out*/);

};

}     // namespace omlplot

#endif
