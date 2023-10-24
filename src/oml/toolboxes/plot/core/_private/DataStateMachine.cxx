/**
* @file DataStateMachine.cxx
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

#include "DataStateMachine.h"
#include "CoreMain.h"
#include <memory>
#include "EvaluatorInt.h"

using namespace std;

namespace omlplot{

    DataStateMachine::DataStateMachine() {
        cm = CoreMain::getInstance();
    }

    DataStateMachine::~DataStateMachine() {
        CoreMain::releaseInstance();
        cm = nullptr;
    }

    unique_ptr<GetData> DataStateMachine::getGetData(const vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();

        Currency input;
        unique_ptr<GetData> gd(new GetData);
        State state = START;
        do {
            if (pos == inputSize){
                if (pos == 0){
                    throw OML_Error(OML_ERR_NUMARGIN);
                }
                break;
            }

            input = getNextInput(inputs, pos);
            switch (state){
            case START: {
                if (input.IsScalar()){
                    double h = input.Scalar();
                    if (cm->isHandle(h)){
                        gd->handles.push_back(h);
                    } else {
                        // TODO: throw
                    }
                } else if (input.IsVector()){
                    gd->handles = input.Vector();
                } else {
                    throw OML_Error(OML_ERR_SCALARVECTOR, pos);
                }
                state = GET_HANDLE;
                ++pos;
                break;
            }
            case GET_HANDLE: {
                if (input.IsString() ){
                    gd->properties.push_back(input.StringVal());
                } else if (input.IsCellArray()){
                    const HML_CELLARRAY *ca = input.CellArray();
                    int row = ca->M();
                    int col = ca->N();
                    for (int j = 0; j < col; j++){
                        for (int i = 0; i < row; i++){
                            Currency c = (*ca)(i, j);
                            if (c.IsString()){
                                gd->properties.push_back(c.StringVal());
                            } else {
                                throw OML_Error(OML_ERR_STRING_STRINGCELL, pos);
                            }
                        }
                    }
                } else {
                    throw OML_Error(OML_ERR_STRING_STRINGCELL, pos);
                }
                state = DONE;
                ++pos;
                break;
            }
            case DONE:
                throw OML_Error(OML_ERR_NUMARGIN);
                break;
            }
        } while (true);
            
        return gd;
    }

    unique_ptr<SetData> DataStateMachine::getSetData(const vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();

        Currency input;
        unique_ptr<SetData> data(new SetData);
        State state = START;
        do {
            if (pos == inputSize){
                if (pos == 0){
                    throw OML_Error(OML_ERR_NUMARGIN);
                }
                break;
            }

            input = getNextInput(inputs, pos);
            switch (state){
            case START: {
                if (input.IsScalar()){
                    double h = input.Scalar();
                    if (cm->isHandle(h)){
                        data->handles.push_back(h);
                    } else {
                        // TODO: throw
                    }
                } else if (input.IsVector()){
                    data->handles = input.Vector();
                } else {
                    throw OML_Error(OML_ERR_SCALARVECTOR, pos);
                }
                state = GET_HANDLE;
                ++pos;
                break;
            }
            case GET_HANDLE:
            case GET_PROP_VALUE: {
                if (input.IsString() ){
                    data->properties.push_back(input.StringVal());
                    ++pos;
                    if (pos < inputSize){
                        data->values.push_back(inputs[pos]);
                        ++pos;
                        state = GET_PROP_VALUE;
                    } else {
                        state = STATE_ERROR;
                    }
                } else {
                    throw OML_Error(OML_ERR_STRING, pos);
                }
                break;
            }
            case DONE:
                throw OML_Error(OML_ERR_NUMARGIN);
                break;
            }
        } while (true);
            
        return data;
    }

    unique_ptr<FigureData> DataStateMachine::getFigureData(const vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();

        Currency input;
        unique_ptr<FigureData> fd(new FigureData);
        State state = START;
        do {
            if (pos == inputSize){
                break;
            }

            input = getNextInput(inputs, pos);
            switch (state){
            case START: {
                if (input.IsScalar()){
                    int h = int(input.Scalar());
                    if (h < 1){
                        throw OML_Error(OML_ERR_POSINTEGER, pos + 1);
                    }
                    fd->fids.push_back(h);
                    ++pos;
                } else if (input.IsString()){
                    fd->fids.push_back(-1);
                    
                } else {
                    state = STATE_ERROR;
                }
                state = GET_HANDLE;
                break;
            }
            case GET_HANDLE:
            case GET_PROP_VALUE: {
                if (input.IsString() ){
                    string prop = input.StringVal();
                    if (cm->isObjectPropertyName<Figure>(prop)){
                        if ((pos + 1) >= inputSize){
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        fd->properties.push_back(prop);
                        fd->values.push_back(inputs[++pos]);
                    }                    
                } else {
                    throw OML_Error(OML_ERR_STRING, pos);
                }
                state = GET_PROP_VALUE;
                ++pos;
                break;
            }
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);
            
        return fd;
    }

    unique_ptr<AxesData> DataStateMachine::getAxesData(const vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();

        Currency input;
        unique_ptr<AxesData> data(new AxesData);
        return data;
        State state = START;
        do {
            if (pos == inputSize){
                break;
            }

            input = getNextInput(inputs, pos);
            switch (state){
            case START: {
                if (input.IsScalar()){
                    double h = input.Scalar();
                    data->handles.push_back(h);
                } else {
                    throw OML_Error(OML_ERR_SCALAR, pos);
                }
                state = GET_HANDLE;
                ++pos;
                break;
            }
            case GET_HANDLE:
            case GET_PROP_VALUE: {
                if (input.IsString() ){
                    string prop = input.StringVal();
                    if (cm->isObjectPropertyName<Axes>(prop)){
                        data->properties.push_back(prop);
                    }
                    
                } else {
                    throw OML_Error(OML_ERR_STRING, pos);
                }
                state = DONE;
                ++pos;
                break;
            }
            case DONE:
                throw OML_Error(OML_ERR_NUMARGIN);
                break;
            }
        } while (true);
            
        return data;
    }

    vector<LineData> DataStateMachine::getLineData(const vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData pd;
        vector<LineData> res;
        State state = START;
        std::shared_ptr<hwMatrix> x;
        std::shared_ptr<hwMatrix> y;
        int count = 0;

        do {
            // this is the only state to quit the machine.
            if (pos == inputSize){  // reach the end of input
                if (state == GET_X){
                    ExtractData(x.get(), y.get(), res, count);
                }
                break;
            }
            input = getNextInput(inputs, pos);

            switch (state){
            case START:
                if (input.IsScalar()){
                    double h = input.Scalar();
                    if (cm->isAxes(h)){
                        pd.parent = h;
                        _temp_parent = h;
                        ++pos;
                    }
                }
                state = GET_HANDLE;
                break;
            case GET_HANDLE:
                if (input.IsScalar() || input.IsComplex() ||
                    input.IsMatrix() || input.IsNDMatrix() ){
                    x = GetMatrix(input);
                    pos += 1;
                    state = GET_X;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_X:
                if (input.IsScalar() || input.IsComplex() ||
                    input.IsMatrix() || input.IsNDMatrix() ){
                    y = GetMatrix(input);
                    ExtractData(x.get(), y.get(), res, count);
                    
                    pos += 1;
                    state = GET_Y;
                } else if (input.IsString()){
                    // no y input
                    ExtractData(x.get(), y.get(), res, count);
                    state = GET_Y;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_Y:
                if (input.IsScalar() || input.IsComplex() ||
                    input.IsMatrix() || input.IsNDMatrix() ){
                    x = GetMatrix(input);
                    pos += 1;
                    state = GET_X;
                } else if (input.IsString()){
                    string strInput = input.StringVal();
                    while (true){   // get every string inputs
                        if (cm->isObjectPropertyName<Line>(strInput)){
                            if ( (pos + 1) >= inputSize){
                                // out of range
                                throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                            }
                            string propertyName = strInput;
                            Currency propertyValue = inputs[++pos];

                            size_t end = res.size(), start = end - count;
                            for (size_t i = start; i < end; i++){
                                res[i].properties.push_back(propertyName);
                                res[i].values.push_back(propertyValue);
                            }
                        } else if (true){//isPlotStyle(strInput)){
                            size_t end = res.size(), start = end - count;
                            for (size_t i = start; i < end; i++){
                                res[i].style = strInput;
                            }
                        } else {
                            throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                        }

                        // there are more string input
                        if ( ((pos + 1) < inputSize) &&
                             inputs[pos + 1].IsString()){
                            ++pos;
                            strInput = inputs[pos].StringVal();
                        } else {    // no more string inputs
                            break;  // exit while loop
                        }
                    }

                    // line info saved in res, use pd again to save next info
                    pos += 1;
                    state = GET_FMT;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_FMT:             // connected to one state: GET_X
                if (input.IsScalar() || input.IsComplex() ||
                    input.IsMatrix() || input.IsNDMatrix() ){
                    x = GetMatrix(input);
                    pos += 1;
                    state = GET_X;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while(true);

        _temp_parent = 0;
        return res;
    }

    vector<LineData> DataStateMachine::getBarData(const vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData pd;
        vector<LineData> res;
        State state = START;
        std::shared_ptr<hwMatrix> x;
        std::shared_ptr<hwMatrix> y;
        int count = 0;
        Currency xCategories;

        do {
            // this is the only state to quit the machine.
            if (pos == inputSize){  // reach the end of input
                if (state == GET_X){
                    if (x.get())
                    ExtractData(x.get(), y.get(), res, count);
                    else
                        throw OML_Error(OML_ERR_PLOT_MATXY_NOT_MATCH);
                }
                break;
            }
            input = getNextInput(inputs, pos);

            switch (state){
            case START:
                if (input.IsScalar()){
                    double h = input.Scalar();
                    if (cm->isAxes(h)){
                        pd.parent = h;
                        _temp_parent = h;
                        ++pos;
                    }
                }
                state = GET_HANDLE;
                break;
            case GET_HANDLE:
                if (input.IsScalar() || input.IsComplex() ||
                    input.IsMatrix() || input.IsNDMatrix() ){
                    x = GetMatrix(input);
                    ++pos;
                    state = GET_X;
                } else if (input.IsCellArray()) {
                    HML_CELLARRAY* c = input.CellArray();
                    int count = c->Size();
                    for (int i = 0; i < count; ++i) {
                        Currency cur = (*c)(i);
                        if (!cur.IsString()) {
                            throw OML_Error("Error: invalid input; must be a cell array of strings");
                        }
                    }
                    xCategories = input;
                    state = GET_X;
                    ++pos;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_X:
                if (input.IsScalar() || input.IsComplex() ||
                    input.IsMatrix() || input.IsNDMatrix() ){
                    
                    if (input.IsScalar() && x.get() && x->Size() > 1){
                        //this input is width
                        ExtractData(x.get(), y.get(), res, count);

                        size_t end = res.size(), start = end - count;
                        for (size_t i = start; i < end; i++){
                            res[i].properties.push_back("barwidth");
                            res[i].values.push_back(input);
                        }
                        ++pos;
                        state = GET_WIDTH;
                    }
                    else {
                        if (xCategories.IsCellArray()) {
                            y = GetMatrix(input);
                            if (y->N() != xCategories.CellArray()->Size())
                                throw OML_Error("Error: length of categories must match the number of columns of y");
                            y->Transpose();
                            ExtractData(y.get(), x.get(), res, count);
                            size_t end = res.size(), start = end - count;
                            for (size_t i = start; i < end; i++) {
                                res[i].xCategories = xCategories;
                            }
                    }
                        else {
                    y = GetMatrix(input);
                    ExtractData(x.get(), y.get(), res, count);
                        }
                        ++pos;
                    state = GET_Y;
                    }
                } else {
                    // no y input
                    if (x.get()) {
                    ExtractData(x.get(), y.get(), res, count);
                    state = GET_Y;
                }
                    else {
                        state = STATE_ERROR;
                    }

                }
				break;
            case GET_Y:
            case GET_WIDTH:
            case GET_FMT:
            case GET_PROP_VALUE:
                if (input.IsScalar()){ // width
                    size_t end = res.size(), start = end - count;
                    for (size_t i = start; i < end; i++){
                        res[i].properties.push_back("barwidth");
                        res[i].values.push_back(input);
                    }
                    ++pos;
                    state = GET_WIDTH;
                } else if (input.IsString()){ // prop/value or style
                    string strInput = input.StringVal();

                    if (cm->isObjectPropertyName<HggroupBar>(strInput)){
                        if ( (pos + 1) >= inputSize){
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        string propertyName = strInput;
                        Currency propertyValue = inputs[++pos];

                        size_t end = res.size(), start = end - count;
                        for (size_t i = start; i < end; i++){
                            res[i].properties.push_back(propertyName);
                            res[i].values.push_back(propertyValue);
                        }
                        ++pos;
                        state = GET_PROP_VALUE;
                    } else if (strInput == "grouped" ||
                               strInput == "stacked"){
                        size_t end = res.size(), start = end - count;
                        for (size_t i = start; i < end; i++){
                            res[i].properties.push_back("barlayout");
                            res[i].values.push_back(input);
                        }
                        ++pos;
                        state = GET_FMT;
                        break;
                    } else {
                        size_t end = res.size(), start = end - count;
                        for (size_t i = start; i < end; i++){
                            res[i].style = strInput;
                        }
                        state = GET_FMT;
                        ++pos;
                    }
                } else {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while(true);

        _temp_parent = 0;
        return res;
    }

    vector<LineData> DataStateMachine::getFillData(const vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData pd;
        vector<LineData> res;
        State state = START;
        std::shared_ptr<hwMatrix> x;
        std::shared_ptr<hwMatrix> y;
        std::shared_ptr<hwMatrix> colorMat;
        int count = 0;

        do {
            // this is the only state to quit the machine.
            if (pos == inputSize){  // reach the end of input
                if (state == GET_COLOR_MAT ||
                    state == GET_PROP_VALUE){
                    break;
                } else {
                    throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                }
            }
            input = getNextInput(inputs, pos);

            switch (state){
            case START:
                if (input.IsScalar()){
                    double h = input.Scalar();
                    if (cm->isAxes(h)){
                        pd.parent = h;
                        _temp_parent = h;
                        ++pos;
                    }
                }
                state = GET_HANDLE;
                break;
            case GET_HANDLE:
                if (input.IsScalar() || input.IsMatrix()){
                    x = GetMatrix(input);
                    pos += 1;
                    state = GET_X;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_X:
                if (input.IsScalar() || input.IsMatrix() ){
                    y = GetMatrix(input);
                    ExtractData(x.get(), y.get(), res, count);
                    pos += 1;
                    state = GET_Y;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_Y:
                if (input.IsMatrix() ){
                    colorMat = GetMatrix(input);
                    size_t end = res.size(), start = end - count;
                    for (size_t i = start; i < end; i++){
                        //res[i].style = strInput;
                    }
                    pos += 1;
                    state = GET_COLOR_MAT;
                } else if (input.IsString()){
                    string strInput = input.StringVal();
                    size_t end = res.size(), start = end - count;
                    for (size_t i = start; i < end; i++){
                        res[i].style = strInput;
                    }
                    pos += 1;
                    state = GET_COLOR_MAT;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_COLOR_MAT:
                if (input.IsScalar() || input.IsMatrix() ){
                    x = GetMatrix(input);
                    pos += 1;
                    state = GET_X;
                } else if (input.IsString()){
                    string strInput = input.StringVal();
                    if (cm->isObjectPropertyName<Fill>(strInput)){
                        if ( (pos + 1) >= inputSize){
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        string propertyName = strInput;
                        Currency propertyValue = inputs[++pos];

                        size_t end = res.size(), start = end - count;
                        for (size_t i = start; i < end; i++){
                            res[i].properties.push_back(propertyName);
                            res[i].values.push_back(propertyValue);
                        }
                        ++pos;
                        state = GET_PROP_VALUE;
                    } else {
                        throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                    }
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_PROP_VALUE:
                if (input.IsString()){
                    string strInput = input.StringVal();
                    if (cm->isObjectPropertyName<Fill>(strInput)){
                        if ( (pos + 1) >= inputSize){
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        string propertyName = strInput;
                        Currency propertyValue = inputs[++pos];

                        size_t end = res.size(), start = end - count;
                        for (size_t i = start; i < end; i++){
                            res[i].properties.push_back(propertyName);
                            res[i].values.push_back(propertyValue);
                        }
                        ++pos;
                        state = GET_PROP_VALUE;
                    } else {
                        throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                    }
                } else {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while(true);

        _temp_parent = 0;
        return res;
    }

    vector<LineData> DataStateMachine::getSurfData(const vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData data;
        vector<LineData> res;
        State state = START;
        std::shared_ptr<hwMatrix> x;
        std::shared_ptr<hwMatrix> y;
        std::shared_ptr<hwMatrix> z;

        do {
            // this is the only state to quit the machine.
            if (pos == inputSize){  // reach the end of input
                if (state == GET_Z || state == GET_PROP_VALUE){
                    break;
                } else {
                    throw OML_Error(OML_ERR_OPT_UNSUPPORTED, pos + 1);
                }
            }
            input = getNextInput(inputs, pos);

            switch (state){
            case START:
                if (input.IsScalar() || input.IsMatrix()){
                    if (inputSize >= 3 &&
                        (inputs[1].IsScalar() ||
                         inputs[1].IsMatrix()) &&
                        (inputs[2].IsScalar() ||
                         inputs[2].IsMatrix()) ){ // Z
                        x = GetMatrix(input);
                        data.x = MatrixToVector(x.get());
                        if (x->IsVector()){
                            data.xcolcount = 1;
                        } else {
                            data.xcolcount = x->N();
                        }
                        ++pos;
                        state = GET_X;
                    } else {
                        z = GetMatrix(input);
                        data.z = MatrixToVector(z.get());
                        data.zcolcount = z->N();
                        x = GetMatrixIndex(z, BY_COL);
                        data.x = MatrixToVector(x.get());
                        data.xcolcount = 1;
                        y = GetMatrixIndex(z, BY_ROW);
                        data.y = MatrixToVector(y.get());
                        data.ycolcount = 1;
                        ++pos;
                        state = GET_Z;
                    }
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_X:
                if (input.IsScalar() || input.IsMatrix()){
                    y = GetMatrix(input);
                    data.y = MatrixToVector(y.get());
                    if (y->IsVector()){
                        data.ycolcount = 1;
                    } else {
                        data.ycolcount = y->N();
                    }
                    ++pos;
                    state = GET_Y;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_Y:
                if (input.IsScalar() || input.IsMatrix()){
                    z = GetMatrix(input);
                    data.z = MatrixToVector(z.get());
                    if (z->IsVector()){
                        data.zcolcount = 1;
                    } else {
                        data.zcolcount = z->N();
                    }
                    ++pos;
                    state = GET_Z;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_Z:
                if (input.IsScalar() || input.IsMatrix()){
                    throw OML_Error(OML_ERR_OPT_UNSUPPORTED, pos + 1);
                    //data.colorMat = GetMatrix(input);
                    ++pos;
                    state = GET_COLOR_MAT;
                } else if (input.IsString()){
                    string strInput = input.StringVal();
                    if (strInput == "curves" || strInput == "text" ||
                        cm->isObjectPropertyName<Surface>(strInput)){
                        if ((pos + 1) >= inputSize){
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        data.properties.push_back(strInput);
                        data.values.push_back(inputs[++pos]);
                        ++pos;
                        state = GET_PROP_VALUE;                    
                    } else {
                        data.style = strInput;
                        ++pos;
                    }
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_COLOR_MAT:
            case GET_PROP_VALUE:
                if (input.IsString()){
                    string strInput = input.StringVal();
                    if (strInput == "curves" || strInput == "text" || 
                        cm->isObjectPropertyName<Surface>(strInput)){
                        if ((pos + 1) >= inputSize){
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        data.properties.push_back(strInput);
                        data.values.push_back(inputs[++pos]);
                        ++pos;
                        state = GET_PROP_VALUE;                    
                    } else {
                        data.style = strInput;
                        ++pos;
                    }
                } else {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while(true);
        res.push_back(std::move(data));
        return res;
    }

    std::unique_ptr<LegendData> DataStateMachine::getLegendData(const std::vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;

        std::unique_ptr<LegendData> ld(new LegendData);
        State state = START;

        do{
            if (pos == inputSize){  // reach the end of input
                break;
            }
            input = getNextInput(inputs, pos);

            switch (state){
            case START:
                if (input.IsString()){
                    string str = input.StringVal();
                    if(cm->isObjectPropertyName<Legend>(str)) {
                        state = GET_PROP_VALUE;;
                    }
                    else {
                    ld->legends.push_back(input.StringVal());
                    ++pos;
                    state = GET_LEGEND;
                    }
                } else if (input.IsMatrix()){
                    if ((pos + 1) >= inputSize) {
                        // out of range
                        throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                    }
                    std::shared_ptr<hwMatrix> h = GetMatrix(input);
                    int matSize = h->Size();
                    ++pos;
                    input = getNextInput(inputs, pos);
                    if (! input.IsCellArray()){
                        throw OML_Error(OML_ERR_OPTIONVAL, pos + 1,
                                        OML_VAR_PARAMETER);
                    }
                    HML_CELLARRAY *c = input.CellArray();
                    int cellSize = c->Size();
                    if (matSize != cellSize){
                        throw OML_Error(OML_ERR_PLOT_DIM_NOT_MATCH, pos,
                                        OML_VAR_PARAMETER);
                    }
                    for (int i = 0; i < matSize; ++i){
                        ld->handles.push_back((*h)(i));
                        ld->legends.push_back((*c)(i).StringVal());
                    }
                    ++pos;
                    state = GET_HANDLE;
                } else if (input.IsCellArray()){
                    HML_CELLARRAY *c = input.CellArray();
                    int count = c->Size();
                    for (int i = 0; i < count; ++i){
                        ld->legends.push_back((*c)(i).StringVal());
                    }
                    ++pos;
                    state = GET_LEGEND;
                } else {
                    throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                }
                break;
            case GET_HANDLE:
                if (input.IsString()){
                    state = GET_PROP_VALUE;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_LEGEND:
                if (input.IsString()){
                    string str = input.StringVal();
                    if (cm->isObjectPropertyName<Legend>(str)) {
                        state = GET_PROP_VALUE;
                    }
                    else {
                    ld->legends.push_back(input.StringVal());
                    ++pos;
                    state = GET_LEGEND;
                    }
                } else if (input.IsCellArray()){
                    HML_CELLARRAY *c = input.CellArray();
                    int count = c->Size();
                    for (int i = 0; i < count; ++i){
                        ld->legends.push_back((*c)(i).StringVal());
                    }
                    ++pos;
                    state = GET_LEGEND;
                } else if (input.IsString()){
                    state = GET_PROP_VALUE;
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_PROP_VALUE:
                if (input.IsString()) {
                    string str = input.StringVal();
                    if (cm->isObjectPropertyName<Legend>(str)) {
                        if ((pos + 1) >= inputSize) {
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        ld->properties.push_back(str);
                        ld->values.push_back(inputs[++pos]);
                ++pos;
                        state = GET_PROP_VALUE;
                    }
                    else {
                        state = STATE_ERROR;
                    }
                }
                else {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);
        return ld;
    }

    std::unique_ptr<ColorbarData> DataStateMachine::getColorbarData(const std::vector<Currency>& inputs)
    {
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        std::unique_ptr<ColorbarData> data(new ColorbarData);

        if (pos == inputSize)
        {
            // reach the end of input
            data->toggleVisibility = true;
            return data;
        }

        State state = START;

        do {
            // this is the only state to quit the machine.
            if (pos == inputSize)
            {
                // reach the end of input
                break;
            }
            input = getNextInput(inputs, pos);

            switch (state)
            {
            case START:
            {
                if (input.IsScalar())
                {
                    // axis handle -> go to next input
                    ++pos;
                    break;
                }
                else if (input.IsString())
                {
                    std::string val = input.StringVal();
                    if (val == "on" || val == "off")
                    {
                        data->properties.push_back(std::string("visible"));
                        data->values.push_back(input);
                        pos++;
                    }
                    state = GET_PROP_VALUE;
                    break;
                }
                else
                {
                    if (!input.IsVector() || !input.IsRealVector() ||
                        (input.IsVector() && input.Vector().size() != 2))
                    {
                        throw OML_Error(OML_ERR_VECTOR2, pos + 1);
                    }
                    data->properties.push_back(std::string("visible"));
                    data->values.push_back(Currency("on"));
                    data->properties.push_back(std::string("clim"));
                    data->values.push_back(input);
                    state = GET_PROP_VALUE;
                    pos += 1;
                    break;
                }
                state = STATE_ERROR;
                pos += 1;
                break;
            }
            case GET_PROP_VALUE:
                if (input.IsString())
                {
                    std::string str = input.StringVal();
                    if ((pos + 1) >= inputSize)
                    {
                        throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                    }
                    data->properties.push_back(str);
                    data->values.push_back(inputs[++pos]);
                    ++pos;
                    state = GET_PROP_VALUE;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            default:
                throw OML_Error(OML_ERR_NUMARGIN);
            }
        } while (true);
        return data;
    }

    std::vector<LineData> DataStateMachine::getQuiverData(const std::vector<Currency>& inputs)
    {
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData pd;
        vector<LineData> res;
        State state = START;
        std::shared_ptr<hwMatrix> x;
        std::shared_ptr<hwMatrix> y;
        std::shared_ptr<hwMatrix> u;
        std::shared_ptr<hwMatrix> v;
        int count = 0;

        bool hasXYValues = false;
        int inputS = static_cast<int>(inputs.size());
        do {
            // this is the only state to quit the machine.
            if (pos == inputSize){
                if (pd.x.empty())
                    throw OML_Error(OML_ERR_NUMARGIN);
                break;
            }
                

            input = getNextInput(inputs, pos);

            switch (state)
            {
            case START:
                if (input.IsScalar() && cm->isAxes(input.Scalar())) {
                    pd.parent = input.Scalar();
                    ++pos;
                    hasXYValues = inputS >= 5;
                }
                else {
                    hasXYValues = inputS >= 4;
                }
                state = GET_HANDLE;
                break;
            case GET_HANDLE:
                if (input.IsScalar() || input.IsMatrix())
                {
                    if (hasXYValues)
                    {
                        x = GetMatrix(input);
                        state = GET_X;
                    }
                    else
                    {
                        u = GetMatrix(input);
                        state = GET_U;
                    }
                    pos += 1;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_X:
                if (input.IsScalar() || input.IsMatrix() && hasXYValues)
                {
                    y = GetMatrix(input);
                    pos += 1;
                    state = GET_Y;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_Y:
                if (input.IsScalar() || input.IsMatrix())
                {
                    u = GetMatrix(input);
                    pos += 1;
                    state = GET_U;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_U:
                if (input.IsScalar() || input.IsMatrix())
                {
                    v = GetMatrix(input);
                    pos += 1;
                    state = GET_V;

                    // end of x, y, u, v input - validate
                    int M = u->M();
                    int N = u->N();

                    if (M != v->M() || N != v->N())
                        throw OML_Error(OML_ERR_PLOT_DIM_NOT_MATCH);

                    pd.u = MatrixToVector(u.get());
                    pd.v = MatrixToVector(v.get());

                    if (hasXYValues) {
                        if (x->M() != M || x->N() != N ||
                            y->M() != M || y->N() != N)
                            throw OML_Error(OML_ERR_PLOT_DIM_NOT_MATCH);

                        pd.x = MatrixToVector(x.get());
                        pd.y = MatrixToVector(y.get());
                    }
                    else {
                        Currency newx(new hwMatrix(M, N, hwMatrix::REAL));
                        Currency newy(new hwMatrix(M, N, hwMatrix::REAL));
                        hwMatrix* newxMat = newx.GetWritableMatrix();
                        hwMatrix* newyMat = newy.GetWritableMatrix();
                        for (int i = 0; i < M; ++i)
                        {
                            for (int j = 0; j < N; ++j)
                            {
                                (*newxMat)(i, j) = j + 1;
                                (*newyMat)(i, j) = i + 1;
                            }
                        }
                        pd.x = MatrixToVector(newxMat);
                        pd.y = MatrixToVector(newyMat);
                    }
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_V:
                if (input.IsScalar())
                {
                    pd.properties.push_back("autoscalefactor");
                    pd.values.push_back(input);
                    ++pos;
                }
                else if (input.IsString())
                {
                    pd.style = input.StringVal();
                    ++pos;
                }
                else
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);
        res.push_back(pd);
        return res;
    }

    std::vector<LineData> DataStateMachine::getBar3Data(const std::vector<Currency>& inputs)
    {
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData pd;
        vector<LineData> res;
        State state = START;
        int count = 0;

        bool hasXYValues = false;
        int inputS = static_cast<int>(inputs.size());
        do {
            // this is the only state to quit the machine.
            if (pos == inputSize) {
                if (pd.z.empty())
                    throw OML_Error(OML_ERR_NUMARGIN);
                break;
            }


            input = getNextInput(inputs, pos);

            switch (state)
            {
            case START:
                if (input.IsScalar() && cm->isAxes(input.Scalar())) 
                {
                    pd.parent = input.Scalar();
                    ++pos;
                }
                state = GET_HANDLE;
                break;
            case GET_HANDLE:
                if (input.IsMatrix())
                {
                    std::shared_ptr<hwMatrix> z = GetMatrix(input);
                    // zcolcount is the number of items per bar row
                    pd.zcolcount = z->IsVector() ? 1 : z->N();
                    int numRows = z->IsVector() ? z->Size(): z->M();
                    z->Transpose();
                    pd.z = MatrixToVector(z.get());
                    
                    pd.x.resize(pd.zcolcount);
                    for (int i = 0; i < pd.zcolcount; ++i)
                    {
                        pd.x[i] = i + 1;
                    }

                    pd.y.resize(numRows);
                    for (int i = 0; i < numRows; ++i)
                    {
                        pd.y[i] = i + 1;
                    }

                    state = GET_Z;
                    pos += 1;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_Z:
                if (input.IsString())
                {
                    std::string str = input.StringVal();
                    if (str == "rowhandle")
                    {
                        pd.bar3PerRow = true;
                        ++pos;
                        break;
                    }
                    else
                    {
                        pd.style = str;
                        state = GET_COLOR_MAT;
                        ++pos;
                        break;
                    }
                }
                state = STATE_ERROR;
                break;
            case GET_COLOR_MAT:
                if (input.IsString() && input.StringVal() == "rowhandle")
                {
                    pd.bar3PerRow = true;
                    ++pos;
                    break;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);
        res.push_back(pd);
        return res;
    }

    void CreateEquallySpacedBins(hwMatrix* bins, double min, double max, int numBins)
    {
        if (min == max)
        {
            double startBin = min - numBins / 2;
            for (int i = 0; i < numBins; i++)
            {
                (*bins)(i) = startBin + i;
            }
        }
        else
        {
            double binXSize = (max - min) / (double)numBins;
            for (int i = 0; i < numBins; i++)
            {
                (*bins)(i) = (i + 0.5) * binXSize + min;
            }
        }
    }

    Currency Get3DHistogramBins(const Currency& inputMat, int nBinsX, int nBinsY)
    {
        if (!inputMat.IsMatrix())
            throw OML_Error(OML_ERR_MATRIX);

        const hwMatrix* mat = inputMat.Matrix();
        if (mat->N() != 2)
            throw OML_Error("Input must be a Nx2 matrix");

        bool isReal = mat->IsReal();
        double initVal1 = isReal ? (*mat)(0, 0) : mat->z(0, 0).Real();
        double initVal2 = isReal ? (*mat)(0, 1) : mat->z(0, 1).Real();
        // find min/max per column
        double min1 = initVal1;
        double max1 = initVal1;
        double min2 = initVal2;
        double max2 = initVal2;
        for (int ii = 0; ii < mat->M(); ++ii)
        {
            double tmp1 = isReal ? (*mat)(ii, 0) : mat->z(ii, 0).Real();
            if (!isnan(tmp1))
            {
                if (tmp1 < min1)
                    min1 = tmp1;
                if (tmp1 > max1)
                    max1 = tmp1;
            }
            double tmp2 = isReal ? (*mat)(ii, 1) : mat->z(ii, 1).Real();
            if (!isnan(tmp2))
            {
                if (tmp2 < min2)
                    min2 = tmp2;
                if (tmp2 > max2)
                    max2 = tmp2;
            }
        }

        hwMatrix* binX = EvaluatorInterface::allocateMatrix(1, nBinsX, true);
        CreateEquallySpacedBins(binX, min1, max1, nBinsX);

        hwMatrix* binY = EvaluatorInterface::allocateMatrix(1, nBinsY, true);
        CreateEquallySpacedBins(binY, min2, max2, nBinsY);

        HML_CELLARRAY* cell = EvaluatorInterface::allocateCellArray(1, 2);

        (*cell)(0) = Currency(binX);
        (*cell)(1) = Currency(binY);

        return Currency(cell);
    }

    std::vector<Currency> Get3DHistogramBinsAndDists(const Currency& inputMat,
        const Currency& binsInput, int binsType)
    {
        // Validate inputs
        if (!inputMat.IsMatrix())
            throw OML_Error(OML_ERR_MATRIX);

        const hwMatrix* mat = inputMat.Matrix();
        if (mat->N() != 2)
            throw OML_Error("must be a Nx2 matrix");

        bool isBinsN = binsType == 0;
        bool isCenters = binsType == 1;
        bool isEdges = binsType == 2;

        if (!(isBinsN || isCenters || isEdges))
            throw OML_Error(OML_ERR_OPTION);

        if (isCenters || isEdges)
        {
            if (!binsInput.IsCellArray())
                throw OML_Error(OML_ERR_CELLARRAY);

            HML_CELLARRAY* bc = binsInput.CellArray();
            if (bc->Size() != 2 || !(*bc)(0).IsRealVector() || !(*bc)(1).IsRealVector())
                throw OML_Error(OML_Error(OML_ERR_CELLARRAY).GetErrorMessage() + " that contains 2 vectors");

            // monotonous increasing values
            std::vector<double> tmp = (*bc)(0).Vector();
            for (int i = 1; i < static_cast<int>(tmp.size()); ++i)
            {
                if (tmp[i] <= tmp[i - 1])
                    throw OML_Error("values must be monotonously increasing");
            }

            tmp = (*bc)(1).Vector();
            for (int i = 1; i < static_cast<int>(tmp.size()); ++i)
            {
                if (tmp[i] <= tmp[i - 1])
                    throw OML_Error("values must be monotonously increasing");
            }
        }

        // Create "center" and "edges" matrices
        Currency binCenters; // hist3 output
        hwMatrix* edgeX = nullptr;
        hwMatrix* edgeY = nullptr;

        if (isBinsN || isCenters)
        {
            if (isBinsN)
            {
                if (binsInput.IsVector() && binsInput.Vector().size() == 2)
                {
                    std::vector<double> nbins = binsInput.Vector();
                    binCenters = Get3DHistogramBins(inputMat, (int)nbins[0], (int)nbins[1]);
                }
                else
                {
                    binCenters = Get3DHistogramBins(inputMat, 10, 10);
                }
            }
            else
            {
                binCenters = binsInput;
            }

            // Calculate the bin edges
            HML_CELLARRAY* bc = binCenters.CellArray();
            const hwMatrix* cntr1 = (*bc)(0).Matrix();
            const hwMatrix* cntr2 = (*bc)(1).Matrix();

            edgeX = EvaluatorInterface::allocateMatrix(1, cntr1->Size() - 1, true);
            for (int i = 1; i < cntr1->Size(); i++)
            {
                double tmpEdge = (*cntr1)(i - 1) + ((*cntr1)(i) - (*cntr1)(i - 1)) / 2.0;
                (*edgeX)(i - 1) = tmpEdge;
            }

            edgeY = EvaluatorInterface::allocateMatrix(1, cntr2->Size() - 1, true);
            for (int i = 1; i < cntr2->Size(); i++)
            {
                double tmpEdge = (*cntr2)(i - 1) + ((*cntr2)(i) - (*cntr2)(i - 1)) / 2.0;
                (*edgeY)(i - 1) = tmpEdge;
            }
        }
        else
        {
            HML_CELLARRAY* bc = binsInput.CellArray();
            edgeX = EvaluatorInterface::allocateMatrix((*bc)(0).Matrix());
            edgeY = EvaluatorInterface::allocateMatrix((*bc)(1).Matrix());

            // calculate the center of the bins
            hwMatrix* centerX = EvaluatorInterface::allocateMatrix(1, edgeX->Size(), true);
            double tmpD = 0;
            int i = 1;
            for (; i < edgeX->Size(); i++)
            {
                tmpD = ((*edgeX)(i) - (*edgeX)(i - 1)) / 2.0;
                (*centerX)(i - 1) = (*edgeX)(i - 1) + tmpD;
            }
            (*centerX)(i - 1) = (*edgeX)(i - 1) + tmpD;

            hwMatrix* centerY = EvaluatorInterface::allocateMatrix(1, edgeY->Size(), true);
            i = 1;
            for (; i < edgeY->Size(); i++)
            {
                tmpD = ((*edgeY)(i) - (*edgeY)(i - 1)) / 2.0;
                (*centerY)(i - 1) = (*edgeY)(i - 1) + tmpD;
            }
            (*centerY)(i - 1) = (*edgeY)(i - 1) + tmpD;

            HML_CELLARRAY* cell = EvaluatorInterface::allocateCellArray(1, 2);

            (*cell)(0) = Currency(centerX);
            (*cell)(1) = Currency(centerY);

            binCenters = Currency(cell);
        }

        if (!edgeX || !edgeY)
            throw OML_Error(OML_ERR_PLOT_UNKNOWN_ERROR);

        int dM = isEdges ? edgeX->Size() : edgeX->Size() + 1;
        int dN = isEdges ? edgeY->Size() : edgeY->Size() + 1;
        hwMatrix* dists = EvaluatorInterface::allocateMatrix(dM, dN, 0.0);

        bool isReal = mat->IsReal();
        // count elements in each bin
        for (int ii = 0; ii < mat->M(); ++ii)
        {
            double tmp1 = isReal ? (*mat)(ii, 0) : mat->z(ii, 0).Real();
            double tmp2 = isReal ? (*mat)(ii, 1) : mat->z(ii, 1).Real();

            if (isnan(tmp1) || isnan(tmp2))
                continue;

            bool skipRow = false;
            int j = 0;
            for (; j < edgeX->Size(); ++j)
            {
                if (tmp1 < (*edgeX)(j))
                {
                    if (isEdges && j == 0)
                        skipRow = true; // if 'edges' are given skip items less than min
                    break;
                }
            }
            if (skipRow)
                continue;

            int k = 0;
            for (; k < edgeY->Size(); ++k)
            {
                if (tmp2 < (*edgeY)(k))
                {
                    if (isEdges && k == 0)
                        skipRow = true; // if 'edges' are given skip items less than min
                    break;
                }
            }
            if (skipRow)
                continue;

            if (isEdges)
            {
                // if 'edges' skip items in the last bin except if they are equal to the edge
                if (j == edgeX->Size() && tmp1 != (*edgeX)(edgeX->Size() - 1))
                    continue;
                if (k == edgeY->Size() && tmp2 != (*edgeY)(edgeY->Size() - 1))
                    continue;
            }
            if (isEdges)
            {
                --j;
                --k;
            }
            if (j < dists->M() && k < dists->N())
            {
                (*dists)(j, k)++;
            }
        }

        delete edgeX;
        delete edgeY;

        std::vector<Currency> ret;
        ret.push_back(dists);
        ret.push_back(binCenters);

        return ret;
    }

    std::vector<LineData> DataStateMachine::getHist3Data(const std::vector<Currency>& inputs, std::vector<Currency>& hist3dData)
    {
        int binType = -1; // 0: nbins, 1: centers, 2: edges
        Currency binsInput;
        Currency colorInput;
        if (!inputs[0].IsMatrix())
            throw OML_Error(OML_ERR_MATRIX, 1);

        const hwMatrix* mtxInput = inputs[0].Matrix();
        for (int ii = 1; ii < static_cast<int>(inputs.size()); ++ii)
        {
            if (inputs[ii].IsString())
            {
                // check if a next argument exists
                if ((ii + 1) >= inputs.size())
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, ii + 2, OML_VAR_PARAMETER);
                }

                // check if one of the options - 'nbins', 'ctrs', 'edges', or one of the bar3 properties
                std::string opt = inputs[ii].StringVal();
                std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);
                if (opt == "nbins")
                {
                    if (inputs[ii + 1].IsVector() && inputs[ii + 1].Vector().size() == 2)
                    {
                        binsInput = inputs[++ii];
                        binType = 0;
                    }
                    else
                    {
                        throw OML_Error(OML_ERR_VECTOR2, ii + 2, OML_VAR_PARAMETER);
                    }
                }
                else if (opt == "ctrs")
                {
                    if (inputs[ii + 1].IsCellArray())
                    {
                        binsInput = inputs[++ii];
                        binType = 1;
                    }
                    else
                    {
                        throw OML_Error(OML_ERR_CELLARRAY, ii + 2, OML_VAR_PARAMETER);
                    }
                }
                else if (opt == "edges")
                {
                    if (inputs[ii + 1].IsCellArray())
                    {
                        binsInput = inputs[++ii];
                        binType = 2;
                    }
                    else
                    {
                        throw OML_Error(OML_ERR_CELLARRAY, ii + 2, OML_VAR_PARAMETER);
                    }
                }
                else if (opt == "color") // only property of a bar3
                {
                    colorInput = inputs[++ii];
                }
                else
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, ii + 1, OML_VAR_PARAMETER);
                }

            }
            else if (inputs[ii].IsVector() && inputs[ii].Vector().size() == 2)
            {
                // nbins
                binsInput = inputs[ii];
                binType = 0;
            }
            else if (inputs[ii].IsCellArray())
            {
                // centers
                binsInput = inputs[ii];
                binType = 1;
            }
            else
            {
                throw OML_Error(OML_ERR_OPTIONVAL, ii + 1, OML_VAR_PARAMETER);
            }
        }

        // no bin input given, set the default [10, 10]
        if (binType == -1)
        {
            std::vector<double> nbins = { 10.0, 10.0 };
            binsInput = Currency(nbins);
            binType = 0;
        }

        // calculate histogram bins and dists
        hist3dData = Get3DHistogramBinsAndDists(inputs[0], binsInput, binType);

        // bar3 elements 
        hwMatrix* mat = hist3dData[0].GetWritableMatrix();

        LineData pd;
        //mat->Transpose();
        pd.z = MatrixToVector(mat);
        pd.zcolcount = mat->M();

        HML_CELLARRAY* cell = hist3dData[1].CellArray();
        pd.x = MatrixToVector((*cell)(0).Matrix());
        pd.y = MatrixToVector((*cell)(1).Matrix());

        pd.properties.push_back("color");
        if (!colorInput.IsEmpty())
        {
            pd.values.push_back(colorInput);
        }
        else
        {
            pd.values.push_back("b");
        }

        std::vector<LineData> ret;
        ret.push_back(pd);
        return ret;
    }

    std::vector<LineData> DataStateMachine::getTriplotData(const std::vector<Currency>& inputs)
    {
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData data;
        std::vector<LineData> res;
        State state = START;
        
        int inputS = static_cast<int>(inputs.size());
        do {
            // this is the only state to quit the machine.
            if (pos == inputSize)
            {
                if (data.x.empty() || data.y.empty())
                    throw OML_Error(OML_ERR_NUMARGIN);
                break;
            }

            input = getNextInput(inputs, pos);

            switch (state)
            {
            case START:
                if (input.IsScalar()) {
                    double h = input.Scalar();
                    if (cm->isAxes(h)) {
                        data.parent = h;
                        _temp_parent = h;
                        ++pos;
                    }
                }
                state = GET_HANDLE;
                break;
            case GET_HANDLE:
                if (input.IsMatrix())
                {
                    std::shared_ptr<hwMatrix> tri = GetMatrix(input);
                    if (tri->N() != 3)
                    {
                        state = STATE_ERROR;
                        break;
                    }
                    if (!tri->IsVector())
                        tri->Transpose();

                    data.tri = MatrixToVector(tri.get());
                    data.triCount = tri->IsVector() ? 1 : tri->N();
                    state = GET_X;
                    pos += 1;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_X:
                if (input.IsMatrix())
                {
                    std::shared_ptr<hwMatrix> x = GetMatrix(input);
                    data.x = MatrixToVector(x.get());
                    pos += 1;
                    state = GET_Y;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_Y:
                if (input.IsMatrix())
                {
                    std::shared_ptr<hwMatrix> y =GetMatrix(input);
                    data.y = MatrixToVector(y.get());
                    pos += 1;
                    state = GET_FMT;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_FMT:
                if (input.IsString())
                {
                    std::string strInput = input.StringVal();
                    while (true)
                    {
                        if (cm->isObjectPropertyName<Line>(strInput))
                        {
                            // out of range
                            if ((pos + 1) >= inputSize)
                                throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);

                            data.properties.push_back(strInput);
                            data.values.push_back(inputs[++pos]);
                        }
                        else
                        {
                            if (data.style.empty())
                                data.style = strInput;
                        }

                        // there are more string input
                        if (((pos + 1) < inputSize) && inputs[pos + 1].IsString())
                        {
                            ++pos;
                            strInput = inputs[pos].StringVal();
                        }
                        else
                        {    // no more string inputs, exit while loop
                            break;
                        }
                    }
                    pos += 1;
                }
                else
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);
        res.push_back(std::move(data));
        return res;
    }

    std::vector<LineData> DataStateMachine::getTrisurfData(const std::vector<Currency>& inputs)
    {
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData data;
        std::vector<LineData> res;
        State state = START;

        int inputS = static_cast<int>(inputs.size());
        do {
            // this is the only state to quit the machine.
            if (pos == inputSize)
            {
                if (data.x.empty() || data.y.empty() || data.z.empty())
                    throw OML_Error(OML_ERR_NUMARGIN);
                break;
            }

            input = getNextInput(inputs, pos);

            switch (state)
            {
            case START:
                if (input.IsScalar()) {
                    double h = input.Scalar();
                    if (cm->isAxes(h)) {
                        data.parent = h;
                        _temp_parent = h;
                        ++pos;
                    }
                }
                state = GET_HANDLE;
                break;
            case GET_HANDLE:
                if (input.IsMatrix())
                {
                    std::shared_ptr<hwMatrix> tri = GetMatrix(input);
                    if (tri->N() != 3)
                    {
                        state = STATE_ERROR;
                        break;
                    }
                    if (!tri->IsVector())
                        tri->Transpose();

                    data.tri = MatrixToVector(tri.get());
                    data.triCount = tri->IsVector() ? 1 : tri->N();
                    state = GET_X;
                    pos += 1;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_X:
                if (input.IsMatrix())
                {
                    std::shared_ptr<hwMatrix> x = GetMatrix(input);
                    data.x = MatrixToVector(x.get());
                    pos += 1;
                    state = GET_Y;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_Y:
                if (input.IsMatrix())
                {
                    std::shared_ptr<hwMatrix> y = GetMatrix(input);
                    data.y = MatrixToVector(y.get());
                    pos += 1;
                    state = GET_Z;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_Z:
                if (input.IsMatrix())
                {
                    std::shared_ptr<hwMatrix> z = GetMatrix(input);
                    data.z = MatrixToVector(z.get());
                    pos += 1;
                    state = GET_FMT;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_FMT:
                if (input.IsMatrix() || input.IsNDMatrix())
                {
                    if (data.cData.IsEmpty())
                    {
                        data.cData = input;
                        pos += 1;
                        state = GET_FMT;
                    }
                    else
                    {
                        state = STATE_ERROR;
                    }
                }
                else if (input.IsString())
                {
                    std::string strInput = input.StringVal();
                    while (true)
                    {
                        if (cm->isObjectPropertyName<Line>(strInput))
                        {
                            // out of range
                            if ((pos + 1) >= inputSize)
                                throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);

                            data.properties.push_back(strInput);
                            data.values.push_back(inputs[++pos]);
                        }
                        else
                        {
                            if (data.style.empty())
                                data.style = strInput;
                        }

                        // there are more string input
                        if (((pos + 1) < inputSize) && inputs[pos + 1].IsString())
                        {
                            ++pos;
                            strInput = inputs[pos].StringVal();
                        }
                        else
                        {    // no more string inputs, exit while loop
                            break;
                        }
                    }
                    pos += 1;
                }
                else
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);
        res.push_back(std::move(data));
        return res;
    }

    std::unique_ptr<TextData> DataStateMachine::getTextData(const std::vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;

        std::unique_ptr<TextData> td(new TextData);
        State state = START;

        do{
            if (pos == inputSize){  // reach the end of input
                break;
            }
            input = getNextInput(inputs, pos);

            switch (state){
            case START:
                if (input.IsScalar()){
                    double x = input.Scalar();
                    td->xpos.push_back(x);
                } else if(input.IsVector()){
                    td->xpos = input.Vector();
                } else {
                    throw OML_Error(OML_ERR_REALVECTOR, pos + 1, OML_VAR_PARAMETER);
                }
                ++pos;
                state = GET_X;
                break;
            case GET_X:
                if (input.IsScalar()){
                    if (td->xpos.size() > 1){
                        throw OML_Error(OML_ERR_ARRAYSIZE, 2);
                    }
                    double y = input.Scalar();
                    td->ypos.push_back(y);
                } else if(input.IsVector()){
                    if (td->xpos.size() != input.Vector().size()){
                        throw OML_Error(OML_ERR_ARRAYSIZE, 2);
                    }
                    td->ypos = input.Vector();
                } else {
                    throw OML_Error(OML_ERR_REALVECTOR, pos + 1, OML_VAR_PARAMETER);
                }
                ++pos;
                state = GET_Y;
                break;
            case GET_Y:
                if (input.IsScalar()){
                    if (td->ypos.size() > 1){
                        throw OML_Error(OML_ERR_ARRAYSIZE, 2);
                    }
                    double z = input.Scalar();
                    td->zpos.push_back(z);
                } else if(input.IsVector()){
                    if (td->ypos.size() != input.Vector().size()){
                        throw OML_Error(OML_ERR_ARRAYSIZE, 2);
                    }
                    td->zpos = input.Vector();
                } else {
                    for (int i = 0; i < td->ypos.size(); i++){
                        td->zpos.push_back(0);
                    }                    
                    state = GET_Z;
                    break;
                }
                ++pos;
                state = GET_Z;
                break;
            case GET_Z:
                if (input.IsString()){
                    string s = input.StringVal();
                    for (int i = 0; i < td->zpos.size(); i++){
                        td->text.push_back(s);
                    }                    
                } else if (input.IsCellArray()){
                    HML_CELLARRAY *strCell = input.CellArray();
                    int strSize = strCell->Size();
                    size_t xSize = td->xpos.size();
                    if (strSize == 1){
                        string s = (*strCell)(0).StringVal();
                        for (int i = 0; i < xSize; i++){
                            td->text.push_back(s);
                        }                    
                    } else if (strSize > 1 && xSize == 1){
                        string s = (*strCell)(0).StringVal();
                        for (int i = 1; i < strCell->Size(); i++){
                            s += "\\n";
                            s += (*strCell)(i).StringVal();
                        }
                        for (int i = 0; i < xSize; i++){
                            td->text.push_back(s);
                        }                    
                    } else if (strSize > 1 && strSize == xSize){
                        for (int i = 0; i < strSize; i++){
                            string s = (*strCell)(i).StringVal();
                            td->text.push_back(s);
                        }
                    } else {
                        throw OML_Error(OML_ERR_CELLSIZE, 3);
                    }
                }
                ++pos;
                state = GET_TEXT;
                break;
            case GET_TEXT:
            case GET_PROP_VALUE:
                if (input.IsString()) {
                    string strInput = input.StringVal();
                    if ((pos + 1) >= inputSize) {
                        // out of range
                        throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                    }
                    td->properties.push_back(strInput);
                    td->values.push_back(inputs[++pos]);
                    state = GET_PROP_VALUE;
                ++pos;
                }
                else {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);
        return td;
    }

    vector<LineData> DataStateMachine::getHistData(const vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData pd;
        vector<LineData> res;
        State state = START;
        std::shared_ptr<hwMatrix> x;
        std::shared_ptr<hwMatrix> y;
        int count = 0;

        do {
            // this is the only state to quit the machine.
            if (pos == inputSize){  // reach the end of input
                if (y == 0){
                    y = GetMatrix(Currency(10));
                    //GetHistData();
                }
                break;
            }
            input = getNextInput(inputs, pos);

            switch (state){
            case START:
                if ( input.IsVector() ||
                     input.IsMatrix() ){
                    state = GET_X;
                    x = GetMatrix(input);
                    pos += 1;
                } else if ( input.IsScalar() ){
                    throw OML_Error(OML_ERR_REALVECTOR, pos + 1);
                }else {
                    state = STATE_ERROR;
                }
                break;
            case GET_X:
                if (input.IsScalar() || input.IsVector() ){
                    y = GetMatrix(input);
                    pos += 1;
                    state = GET_Y;
                } else if (input.IsString()){
                    y = GetMatrix(Currency(10));

                    string strInput = input.StringVal();

                    if (cm->isObjectPropertyName<HggroupBar>(strInput)){
                        if ( (pos + 1) >= inputSize){
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        pd.properties.push_back(strInput);
                        pd.values.push_back(inputs[++pos]);
                        ++pos;
                        state = GET_PROP_VALUE;
                    } else {
                        throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                    }
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_Y:
                if (input.IsString()){ // prop/value or style
                    string strInput = input.StringVal();

                    if (cm->isObjectPropertyName<HggroupBar>(strInput)){
                        if ( (pos + 1) >= inputSize){
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        pd.properties.push_back(strInput);
                        pd.values.push_back(inputs[++pos]);
                        ++pos;
                        state = GET_PROP_VALUE;
                    } else {
                        throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                    }
                } else {
                    state = STATE_ERROR;
                }
                break;
            case GET_PROP_VALUE:
                if (input.IsString()){
                    string strInput = input.StringVal();

                    if (cm->isObjectPropertyName<HggroupBar>(strInput)){
                        if ( (pos + 1) >= inputSize){
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        pd.properties.push_back(strInput);
                        pd.values.push_back(inputs[++pos]);
                        ++pos;
                        state = GET_PROP_VALUE;
                    } else {
                        throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                    }
                } else {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);

        std::shared_ptr<hwMatrix> data = x;
        if (! data->IsRealData()){
            throw OML_Error(OML_ERR_REALMATRIX, 1, OML_VAR_MATRIX);
        }
        std::shared_ptr<hwMatrix> binData = y;
        int numBins = (int)(*binData)(0);
        if (numBins <= 0){
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_PARAMETER);
        }

        std::unique_ptr<hwMatrix> bin, pdf;
        int row = (int)data->M();
        int col = (int)data->N();

        bin = std::unique_ptr<hwMatrix>(new hwMatrix(1, numBins, hwMatrix::REAL));
        pdf = std::unique_ptr<hwMatrix>(new hwMatrix(numBins, (row == 1 ? 1 : col), hwMatrix::REAL));
        
        if (row == 1){
            hwMatrix yPoints;
            data->ReadRow(0, yPoints);
            yPoints.Transpose();
            GetHistData(yPoints, numBins, *bin, *pdf);
        } else {
            GetHistData(*data, numBins, *bin, *pdf);
        }

        if (row == 1){         // vector data
            pd.x = MatrixToVector(bin.get());
            pd.y = MatrixToVector(pdf.get());
            res.push_back(pd);
        } else {                    // matrix data
            hwMatrix yPoints;
            for (int i = 0; i < col; i++){ // density on each column
                pdf->ReadColumn(i, yPoints);
                pd.x = MatrixToVector(bin.get());
                pd.y = MatrixToVector(&yPoints);
                res.push_back(pd);
                pd = LineData();
            }
        }
        return res;
    }

    std::vector<LineData> DataStateMachine::getXLineData(const std::vector<Currency>& inputs)
    {
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData pd;
        vector<LineData> res;
        State state = START;
        std::shared_ptr<hwMatrix> x;
        std::shared_ptr<hwMatrix> y;
        int count = 0;

        do {
            // this is the only state to quit the machine.
            if (pos == inputSize)
                break;
            
            input = getNextInput(inputs, pos);

            switch (state) {
            case START:
                if (input.IsScalar()) {
                    double h = input.Scalar();
                    if (cm->isAxes(h)) {
                        pd.parent = h;
                        _temp_parent = h;
                        ++pos;
                    }
                }
                state = GET_HANDLE;
                break;
            case GET_HANDLE:
                if (input.IsScalar() || input.IsComplex() ||
                    input.IsMatrix() || input.IsNDMatrix()) {
                    x = GetMatrix(input);
                    pos += 1;
                    ExtractData(x.get(), y.get(), res, count);
                    state = GET_X;
                }
                else {
                    state = STATE_ERROR;
                }
                break;
            case GET_X:
                if (input.IsString()) {
                    string strInput = input.StringVal();
                    while (true) {   // get every string inputs
                        if (cm->isObjectPropertyName<Line>(strInput)) {
                            if ((pos + 1) >= inputSize) {
                                // out of range
                                throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                            }
                            string propertyName = strInput;
                            Currency propertyValue = inputs[++pos];

                            size_t end = res.size(), start = end - count;
                            for (size_t i = start; i < end; i++) {
                                res[i].properties.push_back(propertyName);
                                res[i].values.push_back(propertyValue);
                            }
                        }
                        else {
                            throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                        }

                        // there are more string input
                        if (((pos + 1) < inputSize) &&
                            inputs[pos + 1].IsString()) {
                            ++pos;
                            strInput = inputs[pos].StringVal();
                        }
                        else {    // no more string inputs
                            break;  // exit while loop
                        }
                    }

                    // line info saved in res, use pd again to save next info
                    pos += 1;
                }
                else {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);

        _temp_parent = 0;
        return res;
    }

    void DataStateMachine::GetHistData(const hwMatrix& data, const int numBins, hwMatrix &bin, hwMatrix &freq){
        double temp_max = data(0);
        double temp_min = data(0);
        int n = data.Size();
        for (int i = 0; i < n; i++) {
            if (data(i) > temp_max){
                temp_max = data(i);
            } else if (data(i) < temp_min) {
                temp_min = data(i);
            }
        }
        double bin_size = (temp_max - temp_min) / (double) numBins;
        for (int i = 0; i < numBins; i++){
            bin(i) = (i + 0.5) * bin_size + temp_min;
        }

        if (bin_size < 1.0e-12){    // bin size too small
            if (numBins == 1){      // one bin is ok
                freq(0) = n;
                return;
            } else {
                throw OML_Error(OML_ERR_PLOT_ZERORANGE, 1, OML_VAR_MATRIX);
            }
        }

        freq.SetElements(0);

        int index;
        for (int col = 0; col < data.N(); ++col){
            for (int row = 0; row < data.M(); ++row){
                index = (int) floor( (data(row, col) - temp_min - bin_size / 2) / bin_size);
                if (index < 0) {
                    index = 0;
                }
                freq(index, col) += 1;
            }
        }
    }

    std::unique_ptr<LimData> DataStateMachine::getLimData(const std::vector<Currency> &inputs){
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        State state = START;
        int count = 0;
        std::unique_ptr<LimData> ld(new LimData);

        do {
            // this is the only state to quit the machine.
            if (pos == inputSize){  // reach the end of input
                break;
            }
            input = getNextInput(inputs, pos);

            switch (state){
            case START:
                if (input.IsScalar()){
                    double h = input.Scalar();
                    if (cm->isAxes(h)){
                        ld->handle = h;
                        ++pos;
                    } else {
                        throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
                    }
                }
                state = GET_HANDLE;
                break;
            case GET_HANDLE:
                if (input.IsVector()) {
                    if (! input.IsRealVector()){
                        throw OML_Error(OML_ERR_REALVECTOR, pos + 1);
                    }
                    vector<double> lim = input.Vector();
                    if ((lim.size() != 2) &&
                        (lim.size() != 4) &&
                        (lim.size() != 6)  ){
                        throw OML_Error(OML_ERR_VECTOR, pos + 1);
                    }
                    ld->limits = lim;
                    state = GET_LIM;
                    pos += 1;
                } else {
                    throw OML_Error(OML_ERR_VECTOR, pos + 1);
                }
                break;
            case GET_LIM:
            default:
                throw OML_Error(OML_ERR_NUMARGIN);
            }
        } while (true);
        return ld;
    }

    QueryData::SEARCH_OP GetSearchOperator(const std::string& op)
    {
        if (op == "-and")
        {
            return QueryData::OP_AND;
        }
        else if (op == "-or")
        {
            return QueryData::OP_OR;
        }
        else if (op == "-xor")
        {
            return QueryData::OP_XOR;
        }
        else if (op == "-not")
        {
            return QueryData::OP_NOT;
        }
        return QueryData::OP_INV;
    }

    std::unique_ptr<QueryData> DataStateMachine::getQueryData(const std::vector<Currency>& inputs)
    {
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;

        std::unique_ptr<QueryData> data(new QueryData);
        State state = START;

        if (inputSize == 0)
            throw OML_Error(OML_ERR_NUMARGIN);

        do
        {
            if (pos == inputSize)
            {  // reach the end of input
                if (state != GET_PROP_VALUE && state != GET_SINGLE_PROP)
                    throw OML_Error(OML_ERR_NUMARGIN);
                break;
            }
            input = getNextInput(inputs, pos);

            switch (state)
            {
            case START:
                if (input.IsScalar())
                {
                    data->handles.push_back(input.Scalar());
                    ++pos;
                    state = GET_HANDLE;
                }
                else if (input.IsVector())
                {
                    std::vector<double> xx = input.Vector();
                    for (int i = 0; i < xx.size(); i++)
                        data->handles.push_back(xx[i]);
                    ++pos;
                    state = GET_HANDLE;
                }
                else if (input.IsString())
                {
                    std::string val = input.StringVal();
                    if (val == "-property")
                    {
                        state = GET_SINGLE_PROP;
                    }
                    else if (val == "-depth")
                    {
                        state = GET_SEARCH_DEPTH;
                    }
                    else if (val == "flat")
                    {
                        data->m_depth = 0;
                        state = GET_HANDLE;
                        ++pos;
                    }
                    else
                    {
                        // check if operator
                        QueryData::SEARCH_OP op = GetSearchOperator(val);
                        if (op != QueryData::OP_INV)
                        {
                            data->ops.push_back(op);
                            state = GET_HANDLE;
                            ++pos;
                        }
                        else
                        {
                            state = GET_HANDLE;
                        }
                    }
                }
                break;
            case GET_HANDLE:
            case GET_PROP_VALUE:
                if (input.IsString())
                {
                    std::string str = input.StringVal();
                    if (str == "-property")
                    {
                        state = GET_SINGLE_PROP;
                        break;
                    }
                    else if (str == "-depth")
                    {
                        state = GET_SEARCH_DEPTH;
                        break;
                    }
                    else if (str == "flat")
                    {
                        data->m_depth = 0;
                        state = GET_HANDLE;
                        ++pos;
                    }
                    else
                    {
                        // check if operator
                        QueryData::SEARCH_OP op = GetSearchOperator(str);
                        if (op != QueryData::OP_INV)
                        {
                            data->ops.push_back(op);
                            state = GET_HANDLE;
                            ++pos;
                        }
                        else
                        {
                            if ((pos + 1) >= inputSize)
                            {
                                throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                            }
                            data->properties.push_back(str);
                            data->values.push_back(inputs[++pos]);
                            ++pos;
                            state = GET_PROP_VALUE;
                        }
                    }
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_SINGLE_PROP:
                // check there is a property
                if ((pos + 1) >= inputSize)
                {
                    throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                }
                ++pos;
                if (inputs[pos].IsString())
                {
                    data->m_onlyProperty = true;
                    data->properties.push_back(inputs[pos].StringVal());
                    state = GET_SINGLE_PROP;
                }
                else
                {
                    state = STATE_ERROR;
                }
                // if there are more arguments throw error
                if ((pos + 1) < inputSize)
                {
                    throw OML_Error(OML_ERR_NUMARGIN);
                }
                ++pos;
                break;
            case GET_SEARCH_DEPTH:
                // check there is a value
                if ((pos + 1) >= inputSize)
                {
                    throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                }
                ++pos;
                if (inputs[pos].IsScalar())
                {
                    int d = static_cast<int>(inputs[pos].Scalar());
                    if (d < 0)
                    {
                        throw OML_Error(OML_ERR_POSINTEGER, pos + 1);
                    }
                    data->m_depth = d;
                    state = GET_PROP_VALUE;
                }
                else
                {
                    state = STATE_ERROR;
                }
                ++pos;
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);

        if (!data->ops.empty() && data->ops.size() != ((int)data->properties.size() - 1))
            throw OML_Error(OML_ERR_NUMARGIN);

        return data;
    }

    std::vector<LineData> DataStateMachine::getShapeData(const std::vector<Currency>& inputs, 
        const std::string& shape) {
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;

        std::vector<LineData> vd;
        LineData data;// (new LineData);
        State state = START;

        if (inputSize == 0) {
            throw OML_Error(OML_ERR_NUMARGIN);
        }

        do {
            if (pos == inputSize) {  // reach the end of input
                break;
            }
            input = getNextInput(inputs, pos);

            switch (state) {
            case START:
                if (input.IsScalar()) {
                    data.parent = input.Scalar();
                    ++pos;
                }
                else {
                    data.parent = -1;
                }
                state = GET_HANDLE;
                break;
            case GET_HANDLE:
            case GET_PROP_VALUE:
                if (input.IsString()) {
                    string str = input.StringVal();
                    if ((shape == "ellipse" && cm->isObjectPropertyName<Ellipse>(str)) ||
                        (shape == "rectangle" && cm->isObjectPropertyName<Rectangle>(str))) {
                        if ((pos + 1) >= inputSize) {
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        data.properties.push_back(str);
                        data.values.push_back(inputs[++pos]);
                        ++pos;
                        state = GET_PROP_VALUE;
                    }
                    else {
                        state = STATE_ERROR;
                    }
                }
                else {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);
        vd.push_back(data);
        return vd;
    }

    std::vector<LineData> DataStateMachine::getPatchData(const std::vector<Currency>& inputs)
    {
        int pos = 0;
        size_t inputSize = inputs.size();
        Currency input;
        LineData data;
        std::vector<LineData> res;
        State state = START;
        std::shared_ptr<hwMatrix> x;
        std::shared_ptr<hwMatrix> y;
        std::shared_ptr<hwMatrix> z;

        do {
            // this is the only state to quit the machine.
            if (pos == inputSize) {  // reach the end of input
                if (state == GET_COLOR_MAT || state == GET_PROP_VALUE) {
                    break;
                }
                else {
                    throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                }
            }
            input = getNextInput(inputs, pos);

            switch (state) {
            case START:
                if (input.IsScalar() || input.IsMatrix()) {
                    x = GetMatrix(input);
                    data.x = MatrixToVector(x.get());
                    if (x->IsVector())
                        data.xcolcount = 1;
                    else
                        data.xcolcount = x->N();
                    pos += 1;
                    state = GET_X;
                }
                else if (input.IsString())
                {
                    state = GET_COLOR_MAT;
                }
                else
                {
                    state = STATE_ERROR;
                }
                break;
            case GET_X:
                if (input.IsScalar() || input.IsMatrix()) {
                    y = GetMatrix(input);
                    data.y = MatrixToVector(y.get());
                    if (y->IsVector())
                        data.ycolcount = 1;
                    else
                        data.ycolcount = y->N();
                    pos += 1;
                    state = GET_Y;
                }
                else {
                    state = STATE_ERROR;
                }
                break;
            case GET_Y:
                if (input.IsMatrix()) {
                    z = GetMatrix(input);
                    data.z = MatrixToVector(z.get());
                    if (z->IsVector())
                        data.zcolcount = 1;
                    else
                        data.zcolcount = z->N();
                    pos += 1;
                    state = GET_Z;
                }
                else if (input.IsString()) {
                    string strInput = input.StringVal();
                    data.style = strInput;
                    pos += 1;
                    state = GET_COLOR_MAT;
                }
                else {
                    state = STATE_ERROR;
                }
                break;
            case GET_Z:
                if (input.IsString()) {
                    string strInput = input.StringVal();
                    data.style = strInput;
                    pos += 1;
                    state = GET_COLOR_MAT;
                }
                else {
                    state = STATE_ERROR;
                }
                break;
            case GET_COLOR_MAT:
                if (input.IsScalar() || input.IsMatrix()) {
                    res.push_back(data);

                    data = LineData();
                    x = GetMatrix(input);
                    data.x = MatrixToVector(x.get());
                    if (x->IsVector())
                        data.xcolcount = 1;
                    else
                        data.xcolcount = x->N();
                    pos += 1;
                    state = GET_X;
                }
                else if (input.IsString()) {
                    string strInput = input.StringVal();
                    if (cm->isObjectPropertyName<Patch>(strInput) || 
                        strInput =="vertices" || strInput == "faces") {
                        if ((pos + 1) >= inputSize) {
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        data.properties.push_back(strInput);
                        data.values.push_back(inputs[++pos]);
                        ++pos;
                        state = GET_PROP_VALUE;
                    }
                    else {
                        throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                    }
                }
                else {
                    state = STATE_ERROR;
                }
                break;
            case GET_PROP_VALUE:
                if (input.IsString()) {
                    string strInput = input.StringVal();
                    if (cm->isObjectPropertyName<Patch>(strInput) ||
                        strInput == "vertices" || strInput == "faces") {
                        if ((pos + 1) >= inputSize) {
                            // out of range
                            throw OML_Error(OML_ERR_PLOT_MISSING_VALUE, pos + 1);
                        }
                        data.properties.push_back(strInput);
                        data.values.push_back(inputs[++pos]);
                        ++pos;
                        state = GET_PROP_VALUE;
                    }
                    else {
                        throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                    }
                }
                else {
                    state = STATE_ERROR;
                }
                break;
            case STATE_ERROR:
                throw OML_Error(OML_ERR_OPTIONVAL, pos + 1, OML_VAR_PARAMETER);
                break;
            }
        } while (true);
        res.push_back(data);
        return res;
    }

    std::shared_ptr<hwMatrix> DataStateMachine::GetMatrix(const Currency &c){
        std::shared_ptr<hwMatrix> mat;
        if (c.IsScalar()){
            mat.reset(new hwMatrix(1, 1, hwMatrix::REAL));
            (*mat)(0) = c.Scalar();
        } else if (c.IsComplex()){
            mat.reset(new hwMatrix(1, 1, hwMatrix::COMPLEX));
            mat->z(0) = c.Complex();
        } else if (c.IsMatrix() ){
            const hwMatrix *data = c.Matrix();
            mat.reset(new hwMatrix(*data));
        } else if (c.IsNDMatrix()){
            const hwMatrixN *data = c.MatrixN();
            const vector<int> dim = data->Dimensions();
            if (dim.size() != 2){
                throw OML_Error(OML_ERR_OPTIONVAL);
            }
            int m = dim[0], n = dim[1];
            mat.reset(new hwMatrix(m, n,
                                   data->IsReal()? hwMatrix::REAL : hwMatrix::COMPLEX));
            for (int i = 0; i < m; i++){
                for (int j = 0; j < n; j ++){
                    double x = (*data)(i + j * m);
                    (*mat)(i + j *m) = x;
                }
            }
        } else {
            throw OML_Error(OML_ERR_OPTIONVAL);
        }
        return mat;
    }

    vector<double> DataStateMachine::MatrixToVector(const hwMatrix * const mat){
        vector<double> v;
        int count = mat->Size();
		bool isReal = mat->IsReal();
        for (int i = 0; i < count; i++){
            v.push_back(isReal?(*mat)(i):mat->z(i).Real());
        }
        return v;
    }

    std::shared_ptr<hwMatrix> DataStateMachine::GetMatrixIndex(std::shared_ptr<hwMatrix> __m, GetDir dir){
        std::shared_ptr<hwMatrix> mat;
        int col = 0;
        if (dir == BY_ROW){
            col = __m->M();
        } else{
            col = __m->N();
        }
        mat = make_shared<hwMatrix>(1, col, hwMatrix::REAL);
        for (int i = 0; i < col; i++){
            (*mat)(i) = i + 1;
        }
        return mat;
    }

    Currency DataStateMachine::getNextInput(const vector<Currency> &inputs, const int pos){
        return inputs[pos];
    }

    void DataStateMachine::ExtractData(const hwMatrix * const x,
                                       const hwMatrix * const y,
                                       vector<LineData> &result /*out*/,
                                       int &count               /*out*/){
        LineData pd;
        pd.parent = _temp_parent;
        if (y == nullptr){
            // no y input
            if (x->IsReal()){
                if (x->IsVector()){
                    pd.y = MatrixToVector(x);
                    pd.x.resize(pd.y.size());
                    for (int i = 0; i < pd.y.size(); i++){
                        pd.x[i] = i + 1;
                    }
                    count = 1;
                    result.push_back(pd);
                } else {
                    int size = x->M();
                    for (int col = 0; col < x->N(); col++){
                        pd.x.resize(size);
                        pd.y.resize(size);
                        for (int row = 0; row < size; row++){
                            pd.x[row] = row + 1;
                            pd.y[row] = x->operator()(row,col);
                        }
                        result.push_back(pd);
                        pd = LineData();
                        pd.parent = _temp_parent;
                    }
                    count = x->N();
                }
            } else {
                // complex data
                hwMatrix xmat(x->M(), x->N(), hwMatrix::REAL);
                hwMatrix ymat(x->M(), x->N(), hwMatrix::REAL);
                x->UnpackComplex(&xmat, &ymat);
                ExtractData(&xmat, &ymat, result, count); // recursive call
            }
            return;
        }
        if (x->IsVector()){
            if (y->IsVector()){
                pd.x = MatrixToVector(x);
                pd.y = MatrixToVector(y);
                result.push_back(pd);
                count = 1;
            } else {            // y is matrix
                const int xsize = x->Size();
                const int yrow = y->M();
                const int ycol = y->N();
                if (yrow == xsize){
                    hwMatrix columnData;
                    for (int i = 0; i < ycol; i++){
                        y->ReadColumn(i, columnData);
                        pd.x = MatrixToVector(x);
                        pd.y = MatrixToVector(&columnData);
                        result.push_back(pd);
                        pd = LineData();
                        pd.parent = _temp_parent;
                    }
                    count = ycol;
                } else if (ycol == xsize){
                    hwMatrix rowData;
                    for (int i = 0; i < yrow; i++){
                        y->ReadRow(i, rowData);
                        pd.x = MatrixToVector(x);
                        pd.y = MatrixToVector(&rowData);
                        result.push_back(pd);
                        pd = LineData();
                        pd.parent = _temp_parent;
                    }
                    count = yrow;
                } else {
                    throw;
                }
            }
        } else {                // x is matrix
            if (y->IsVector()){
                const int ysize = y->Size();
                const int xrow = x->M();
                const int xcol = x->N();
                if (xrow == ysize){
                    hwMatrix columnData;
                    for (int i = 0; i < xcol; i++){
                        x->ReadColumn(i, columnData);
                        pd.x = MatrixToVector(&columnData);
                        pd.y = MatrixToVector(y);
                        result.push_back(pd);
                        pd = LineData();
                        pd.parent = _temp_parent;
                    }
                    count = xcol;
                } else if (xcol == ysize) {
                    hwMatrix rowData;
                    for (int i = 0; i < xrow; i++){
                        x->ReadRow(i, rowData);
                        pd.x = MatrixToVector(&rowData);
                        pd.y = MatrixToVector(y);
                        result.push_back(pd);
                        pd = LineData();
                        pd.parent = _temp_parent;
                    }
                    count = xrow;
                } else {
                    throw;
                }
            } else {            // y is matrix
                const int xrow = x->M();
                const int xcol = x->N();
                const int yrow = y->M();
                const int ycol = y->N();
                if ((xrow == yrow) &&
                    (xcol == ycol)) {
                    for (int i = 0; i < xcol; i++){
                        hwMatrix xColumnData;
                        hwMatrix yColumnData;
                        x->ReadColumn(i, xColumnData);
                        y->ReadColumn(i, yColumnData);
                        pd.x = MatrixToVector(&xColumnData);
                        pd.y = MatrixToVector(&yColumnData);
                        result.push_back(pd);
                        pd = LineData();
                        pd.parent = _temp_parent;
                    }
                    count = xcol;
                } else {
                    throw;
                }
            }
        }
        _temp_parent = 0;
    }

}     // namespace omlplot
