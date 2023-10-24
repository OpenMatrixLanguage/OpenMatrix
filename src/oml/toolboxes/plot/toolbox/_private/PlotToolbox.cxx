/**
* @file PlotToolbox.cxx
* @date March 2017
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

#include "OmlPlotExport.h"
#include "EvaluatorInt.h" 
#include "Currency.h"
#include "StructData.h"
#include "core/DataStateMachine.h"
#include "core/CoreMain.h"
#include "core/DataType.h"
#include "core/Object.h"
#include <BuiltInFuncsUtils.h>
#include <BuiltInFuncsSystem.h>
#include <memory>

using namespace std;

#define CATCH_ALL_PLOT_ERROR throw OML_Error(OML_ERR_PLOT_UNKNOWN_ERROR);

namespace omlplot{
    void UpdateDateticksForAxis(EvaluatorInterface eval, double axesHandle, const string& axis);

    //FILE *pout;
    DataStateMachine dsm;
    CoreMain *cm = CoreMain::getInstance();

    static Currency castValue(CurrencyAndColor p){
        if (p.isCurrency()){
            Currency c = p.getCurrency();
            if (c.IsBoundObject()){
                Object *o = (Object *)c.BoundObject();
                return o->getHandle();
            } else {
                return c;
            }
        } else if (p.isColor()){
            return p.getColor();
        } else {
            return Currency();
        }
    }

    bool plot(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getLineData(inputs);
        vector<double> hs = cm->plot(vld);
        outputs.push_back(hs);
        return true;
    }

    bool bar(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getBarData(inputs);
        vector<double> hs = cm->bar(vld);
        outputs.push_back(hs);
        return true;
    }

    bool hist(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getHistData(inputs);
        vector<double> hs = cm->hist(vld);
        outputs.push_back(hs);
        return true;
    }

    bool area(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getLineData(inputs);
        vector<double> hs = cm->area(vld);
        outputs.push_back(hs);
        return true;
    }

    bool polar(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getLineData(inputs);
        vector<double> hs = cm->polar(vld);
        outputs.push_back(hs);
        return true;
    }

    bool line(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getLineData(inputs);
        vector<double> hs = cm->line(vld);
        outputs.push_back(hs);
        return true;
    }

    bool scatter(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getLineData(inputs);
        vector<double> hs = cm->scatter(vld);
        outputs.push_back(hs);
        return true;
    }

    bool plot3(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getSurfData(inputs);
        vector<double> hs = cm->plot3(vld);
        outputs.push_back(hs);
        return true;
    }

    bool fill(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getFillData(inputs);
        vector<double> hs = cm->fill(vld);
        outputs.push_back(hs);
        return true;
    }

    bool scatter3(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getSurfData(inputs);
        vector<double> hs = cm->scatter3(vld);
        outputs.push_back(hs);
        return true;
    }

    bool surf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getSurfData(inputs);
        vector<double> hs = cm->surf(vld);
        outputs.push_back(hs);
        return true;
    }

    bool mesh(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getSurfData(inputs);
        vector<double> hs = cm->mesh(vld);
        outputs.push_back(hs);
        return true;
    }

    bool contour3(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getSurfData(inputs);
        vector<double> hs = cm->contour3(vld);
        outputs.push_back(hs);
        return true;
    }

    bool contour(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getSurfData(inputs);
        vector<double> hs = cm->contour(vld);
        outputs.push_back(hs);
        return true;
    }

    bool stem(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getLineData(inputs);
        vector<double> hs = cm->stem(vld);
        outputs.push_back(hs);
        return true;
    }

    bool loglog(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getLineData(inputs);
        vector<double> hs = cm->loglog(vld);
        outputs.push_back(hs);
        return true;
    }

    bool semilogx(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getLineData(inputs);
        vector<double> hs = cm->semilogx(vld);
        outputs.push_back(hs);
        return true;
    }

    bool semilogy(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        vector<LineData> vld = dsm.getLineData(inputs);
        vector<double> hs = cm->semilogy(vld);
        outputs.push_back(hs);
        return true;
    }

    bool get(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        bool res = false;
        unique_ptr<GetData> gd = dsm.getGetData(inputs);
        if (gd->handles.size() == 1){
            double h = gd->handles[0];
            vector<string> props;
            if (gd->properties.size() == 0){
                props = cm->getObjectPropertyNames(h);
            } else {
                props = gd->properties;
            }
            if (props.size() == 1){
                if (!cm->isPropertySupported(h, props[0]))
				{
					BuiltInFuncsUtils::SetWarning(eval, "Property ["+props[0]+"] is not supported in OpenMatrix");
					return false;
				}
                outputs.push_back(castValue(cm->getObjectPropertyValue(h, props[0])));
            } else {
                StructData *sd = new StructData();
                vector<string>::iterator it = props.begin();
                for (; it != props.end(); ++it){
                    string pname = *it;
                    if (!cm->isPropertySupported(h, pname))
						continue;

                    sd->SetValue(0, 0, pname, castValue(cm->getObjectPropertyValue(h, pname)));
                }
                outputs.push_back(sd);
            }
            res = true;
        } else {
            res = false;
        }
        return res;
    }

    bool set(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        unique_ptr<SetData> data = dsm.getSetData(inputs);
		
		// temp solution for setting a warning for not supported properties
		vector<string> notSupportedProperties;
        cm->set(data, notSupportedProperties);
		if (!notSupportedProperties.empty())
		{
			for (vector<string>::const_iterator it = notSupportedProperties.cbegin(); it != notSupportedProperties.cend();++it)
				BuiltInFuncsUtils::SetWarning(eval, "Property [" + (*it) + "] is not supported in OpenMatrix");

		}
        string wrn = cm->GetWarningString();
        if (!wrn.empty())
            BuiltInFuncsUtils::SetWarning(eval, wrn);

        std::vector<std::pair<double, string>> dt = cm->getUpdateDatetickFlag();
        std::vector<std::pair<double, string>>::const_iterator it = dt.cbegin();
        for (; it != dt.cend(); ++it) {
            UpdateDateticksForAxis(eval, it->first, it->second);
        }
        return true;
    }

    bool figure(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        unique_ptr<FigureData> fd = dsm.getFigureData(inputs);
        double h = cm->figure(fd);
        outputs.push_back(h);
        return true;
    }

    bool close(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        double f = 0;
        if (inputs.size() == 1){
            Currency c = inputs[0];
            if (c.IsScalar()){
                f = c.Scalar();
            } else if (c.IsString()){
                string s = c.StringVal();
                if (s != "all"){
                    // throw
                }
                f = -1;
            }
        } else if (inputs.size() == 0){
            f = cm->gcf();
        } else {
            // throw
        }
        cm->close((int)f);
        return true;
    }

    bool gcf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        double h = cm->gcf();
        outputs.push_back(h);
        return true;
    }

    bool clf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        double f = 0;
        if (inputs.size() == 1){
            Currency c = inputs[0];
            if (c.IsScalar()){
                f = c.Scalar();
            } else {
                // throw
            }
        } else if (inputs.size() == 0){
            f = cm->gcf();
        } else {
            // throw
        }
        int h = cm->clf((int)f);
        outputs.push_back(h);
        return true;
    }

    bool axes(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        unique_ptr<AxesData> data = dsm.getAxesData(inputs);
        double h = cm->axes(data);
        outputs.push_back(h);
        return true;
    }

    bool gca(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        //unique_ptr<FigureData> fd = dsm.getFigureData(inputs);
        double h = cm->gca();
        outputs.push_back(h);
        return true;
    }

    bool cla(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        double a = 0;
        if (inputs.size() == 1){
            Currency c = inputs[0];
            if (c.IsScalar()){
                a = c.Scalar();
            } else {
                // throw
            }
        } else if (inputs.size() == 0){
            a = cm->gca();
        } else {
            // throw
        }
        double h = cm->cla(a);
        outputs.push_back(h);
        return true;
    }

    bool hold(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        // TODO: check input first
        if (inputs.size() == 0){
            double axes = cm->gca();
            cm->hold(axes, (! cm->ishold(axes)) );
        } else if(inputs.size() == 1){
            if (inputs[0].IsString()) {
                string state = inputs[0].StringVal();
                double axes = cm->gca();
                if (state == "on"){
                    cm->hold(axes, true);
                } else if (state == "off"){
                    cm->hold(axes, false);
                } else {
                    // throw
                }
            } else if (inputs[0].IsScalar()){
                double axes = inputs[0].Scalar();
                cm->hold(axes, cm->ishold(axes));
            } else {
                // throw
            }
        } else if (inputs.size() == 2){
            double axes = inputs[0].Scalar();
            string state = inputs[1].StringVal();
            if (state == "on"){
                cm->hold(axes, true);
            } else if (state == "off"){
                cm->hold(axes, false);
            } else {
                // throw
            }
        } else {
            // TODO: throw
        }
        return true;
    }

    bool ishold(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        bool hold = false;
        if (inputs.size() == 0){
            hold = cm->ishold(cm->gca());
        } else if (inputs.size() == 1){
            hold = cm->ishold(inputs[0].Scalar());
        } else {
            throw OML_Error(OML_ERR_NUMARGIN);
        }
        outputs.push_back(hold);
        return true;
    }

    bool ishandle(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        if (inputs.size() != 1) {
            throw OML_Error(OML_ERR_NUMARGIN);
        }

        bool res = false;
        if (inputs[0].IsScalar()){
            double h = inputs[0].Scalar();
            res = cm->isHandle(h);
        }
        outputs.push_back(res);
        return true;
    }

    bool isfigure(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        bool res = false;
        if (inputs.size() != 1) {
            throw OML_Error(OML_ERR_NUMARGIN);
        }
        if (! inputs[0].IsScalar()){
            throw OML_Error(OML_ERR_SCALAR, 1, OML_VAR_INPUT);
        }
        double h = inputs[0].Scalar();
        res = cm->isFigure(h);

        outputs.push_back(res);
        return true;
    }

    bool isaxes(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        bool res = false;
        if (inputs.size() != 1) {
            throw OML_Error(OML_ERR_NUMARGIN);
        }
        if (! inputs[0].IsScalar()){
            throw OML_Error(OML_ERR_SCALAR, 1);
        }
        double h = inputs[0].Scalar();
        res = cm->isAxes(h);

        outputs.push_back(res);
        return true;
    }

    bool subplot(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        if (inputs.size() == 3){
            if (!inputs[0].IsPositiveInteger())
                throw OML_Error(OML_ERR_SCALAR, 1);
            
            if (!inputs[1].IsPositiveInteger())
                throw OML_Error(OML_ERR_SCALAR, 2);
            
            if (!inputs[2].IsPositiveInteger() && !inputs[2].IsVector())
                throw OML_Error(OML_ERR_SCALAR, 3);
            
        } else if ( inputs.size() == 1 ){
            if (!inputs[0].IsPositiveInteger())
                throw OML_Error(OML_ERR_SCALAR, 1);
            
        } else {
            throw OML_Error(OML_ERR_NUMARGIN);
        }
        int row = 1, col = 1;
        std::vector<int> active;
        size_t size = inputs.size();
        if (size == 3){
            row = (int)inputs[0].Scalar();
            col = (int)inputs[1].Scalar();
            if (inputs[2].IsScalar())
            {
                active.push_back((int)inputs[2].Scalar());
            }
            else
            {
                std::vector<double> idxD = inputs[2].Vector();
                std::vector<double>::const_iterator it = idxD.cbegin();
                for (; it != idxD.cend(); ++it)
                    active.push_back((int)(*it));
            }
        } else if (size == 1){
            int rcn = (int)inputs[0].Scalar();
            if (rcn < 111){
                throw OML_Error(OML_ERR_OPTION);
            }
            row = rcn / 100;
            col = (rcn / 10) % 10;
            active.push_back(rcn % 10);
        }
        double h = cm->subplot(row, col, active);
        outputs.push_back(h);
        return true;
    }

    bool grid(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        if (inputs.size() == 0) {
            cm->grid(cm->gca());
        } else if (inputs.size() == 1) {
            if (inputs[0].IsString()) {
                string state = inputs[0].StringVal();
                if (state != "on" && state != "off" && state != "minor"){
                    throw OML_Error(OML_ERR_FUNCSWITCH);
                }
                if (state == "minor"){
                    cm->minorgrid(cm->gca());
                } else {
                    cm->grid(cm->gca(), state);
                }
            } else if (inputs[0].IsScalar()) {
                double axes = inputs[0].Scalar();
                cm->grid(axes);
            } else {
                throw OML_Error(OML_ERR_OPTIONVAL, 1, OML_VAR_PARAMETER);
            }
        } else if (inputs.size() == 2) {
            if (inputs[0].IsScalar()) {
                double axes = inputs[0].Scalar();
                string state = inputs[1].StringVal();
                if (state != "on" && state != "off"){
                    throw OML_Error(OML_ERR_FUNCSWITCH, 2);
                }
                cm->grid(axes, state);
            } else if(inputs[0].IsString()){
                string option = inputs[0].StringVal();
                if (option != string("minor")){
                    throw OML_Error(OML_ERR_OPTION);
                }
                string state = inputs[1].StringVal();
                if (state != "on" && state != "off"){
                    throw OML_Error(OML_ERR_FUNCSWITCH, 2);
                }
                cm->minorgrid(cm->gca(), state);
            } else {
                throw OML_Error(OML_ERR_OPTION);
            }
        } else {
            throw OML_Error(OML_ERR_NUMARGIN);
        }
        return true;
    }

    bool title(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        if(inputs.empty() || inputs.size() > 2)
            throw OML_Error(OML_ERR_NUMARGIN);
        else if(inputs.size() == 2) {
            if ( ! inputs[0].IsScalar() ){
                throw OML_Error(OML_ERR_SCALAR, 1);
            }
            if ( ! inputs[1].IsString() ){
                throw OML_Error(OML_ERR_STRING, 2);
            }
        }
        else if(inputs.size() == 1 && ! inputs[0].IsString()){
            throw OML_Error(OML_ERR_STRING, 1);
        }

        double axes = 0;
        string strTitle = "";
        if (inputs.size() == 1) {
            axes = cm->gca();
            strTitle = inputs[0].StringVal();
        } else if (inputs.size() == 2) {
            axes = inputs[0].Scalar();
            strTitle = inputs[1].StringVal();
        }
        outputs.push_back(cm->title(axes, strTitle));
        return true;
    }

    bool xlabel(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        if(inputs.empty() || inputs.size() > 2)
            throw OML_Error(OML_ERR_NUMARGIN);
        else if(inputs.size() == 2) {
            if ( ! inputs[0].IsScalar() ){
                throw OML_Error(OML_ERR_SCALAR, 1);
            }
            if ( ! inputs[1].IsString() ){
                throw OML_Error(OML_ERR_STRING, 2);
            }
        }
        else if(inputs.size() == 1 && ! inputs[0].IsString()){
            throw OML_Error(OML_ERR_STRING, 1);
        }

        double axes = 0;
        string str = "";
        if (inputs.size() == 1) {
            axes = cm->gca();
            str = inputs[0].StringVal();
        } else if (inputs.size() == 2) {
            axes = inputs[0].Scalar();
            str = inputs[1].StringVal();
        }
        outputs.push_back(cm->xlabel(axes, str));
        return true;
    }

    bool ylabel(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        if(inputs.empty() || inputs.size() > 2)
            throw OML_Error(OML_ERR_NUMARGIN);
        else if(inputs.size() == 2) {
            if ( ! inputs[0].IsScalar() ){
                throw OML_Error(OML_ERR_SCALAR, 1);
            }
            if ( ! inputs[1].IsString() ){
                throw OML_Error(OML_ERR_STRING, 2);
            }
        }
        else if(inputs.size() == 1 && ! inputs[0].IsString()){
            throw OML_Error(OML_ERR_STRING, 1);
        }

        double axes = 0;
        string str = "";
        if (inputs.size() == 1) {
            axes = cm->gca();
            str = inputs[0].StringVal();
        } else if (inputs.size() == 2) {
            axes = inputs[0].Scalar();
            str = inputs[1].StringVal();
        }
        outputs.push_back(cm->ylabel(axes, str));
        return true;
    }

    bool zlabel(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        if(inputs.empty() || inputs.size() > 2)
            throw OML_Error(OML_ERR_NUMARGIN);
        else if(inputs.size() == 2) {
            if ( ! inputs[0].IsScalar() ){
                throw OML_Error(OML_ERR_SCALAR, 1);
            }
            if ( ! inputs[1].IsString() ){
                throw OML_Error(OML_ERR_STRING, 2);
            }
        }
        else if(inputs.size() == 1 && ! inputs[0].IsString()){
            throw OML_Error(OML_ERR_STRING, 1);
        }

        double axes = 0;
        string str = "";
        if (inputs.size() == 1) {
            axes = cm->gca();
            str = inputs[0].StringVal();
        } else if (inputs.size() == 2) {
            axes = inputs[0].Scalar();
            str = inputs[1].StringVal();
        }
        outputs.push_back(cm->zlabel(axes, str));
        return true;
    }

    bool axis(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        if (inputs.size() == 1 && inputs.front().IsString())
        {
            std::string inVal = inputs.front().StringVal();
            cm->axis(inVal);
            string wrn = cm->GetWarningString();
            if (!wrn.empty())
                BuiltInFuncsUtils::SetWarning(eval, wrn);
                return true;
            }
        std::unique_ptr<LimData> ld = dsm.getLimData(inputs);
        double h = ld->handle;
        if (h == 0){
            h = cm->gca();
        }
        if (ld->limits.size() == 0){
            vector<double> res = cm->axis(h);
            outputs.push_back(res);
        } else {
            cm->axis(h, ld->limits);
        }
        return true;
    }

    bool xlim(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        std::unique_ptr<LimData> ld = dsm.getLimData(inputs);
        double h = ld->handle;
        if (h == 0){
            h = cm->gca();
        }
        if (ld->limits.size() == 0){
            vector<double> res = cm->xlim(h);
            outputs.push_back(res);
        } else {
            cm->xlim(h, ld->limits);
        }
        return true;
    }

    bool ylim(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        std::unique_ptr<LimData> ld = dsm.getLimData(inputs);
        double h = ld->handle;
        if (h == 0){
            h = cm->gca();
        }
        if (ld->limits.size() == 0){
            vector<double> res = cm->ylim(h);
            outputs.push_back(res);
        } else {
            cm->ylim(h, ld->limits);
        }
        return true;
    }

    bool zlim(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        std::unique_ptr<LimData> ld = dsm.getLimData(inputs);
        double h = ld->handle;
        if (h == 0){
            h = cm->gca();
        }
        if (ld->limits.size() == 0){
            vector<double> res = cm->zlim(h);
            outputs.push_back(res);
        } else {
            cm->zlim(h, ld->limits);
        }
        return true;
    }

    bool legend(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        std::unique_ptr<LegendData> ld = dsm.getLegendData(inputs);
        int nargout = eval.GetNargoutValue();
        if (nargout > 0 && inputs.empty())
        {
            ld->legends.push_back("on");
        }
        double lh = cm->legend(ld);

        if (nargout > 0)
        {
            outputs.push_back(lh);
        }
        return true;
    }

    bool text(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        std::unique_ptr<TextData> td = dsm.getTextData(inputs);
        vector<double> res = cm->text(td);
        outputs.push_back(res);
        std::string wrn = cm->GetWarningString();
        if (!wrn.empty()) {
            BuiltInFuncsUtils::SetWarning(eval, wrn);
        }
        return true;
    }

    bool saveas(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        if (inputs.size() < 2) {
            throw OML_Error(OML_ERR_NUMARGIN);
        }

        double handle = -1.0;
        string filename, fmt;
        int width = -1, height = -1;

        if (!inputs[0].IsScalar()) {
            throw OML_Error(OML_ERR_SCALAR, 1);
        }
        handle = inputs[0].Scalar();

        if (!inputs[1].IsString()) {
            throw OML_Error(OML_ERR_STRING, 2); 
        }
        filename = inputs[1].StringVal();

        if (inputs.size() == 3) {
            if (inputs[2].IsString()) {
                fmt = inputs[2].StringVal();
            }
            else if (inputs[2].IsVector()) {
                std::vector<double> s = inputs[2].Vector();
                if (s.size() == 2)
                {
                    width = (int)s[0];
                    height = (int)s[1];
                }
            }
            else {
                throw OML_Error(OML_ERR_STRINGVECTOR, 3);
            }
        }
        else if (inputs.size() > 3) {
            if (!inputs[2].IsString()) {
                throw OML_Error(OML_ERR_STRING, 3);
            }
            fmt = inputs[2].StringVal();

            if (!inputs[3].IsVector()) {
                throw OML_Error(OML_ERR_VECTOR2, 4);
            }
            std::vector<double> s = inputs[3].Vector();
            if (s.size() == 2)
            {
                width = (int)s[0];
                height = (int)s[1];
            }
        }

        BuiltInFuncsUtils utils;
        std::string dirPath = utils.GetBaseDir(filename);
        if (!utils.DoesPathExist(dirPath))
            throw OML_Error("invalid path; cannot create file [" + filename + "]");

        filename = utils.GetAbsolutePath(filename);
        try {
            cm->saveas(handle, filename, fmt, width, height);
            return true;
        } catch (const OML_Error& e){
            throw e;
        } catch (...){
            CATCH_ALL_PLOT_ERROR
        }
        return false;
    }

    bool dump(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        cm->dump();
        return true;
    }

    bool out(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        if (inputs.size() == 1) {
            if (inputs[0].IsString()) {
                string s = inputs[0].StringVal();
                cm->out(s);
            }
        } else {
            throw OML_Error(OML_ERR_STRING, 1);
        }
        return true;
    }

    bool oml_doNothing(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        static bool firstCall = true;
        if (firstCall){
            std::cout << "\nError: Gnuplot not found. Please install it and try again.\n" << std::endl;
            firstCall = false;
        }
        return true;
    }

	bool box(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        if (inputs.size() == 0) {
            cm->box(cm->gca());
        }
        else if (inputs.size() == 1) {
            if (inputs[0].IsString()) {
                string state = inputs[0].StringVal();
                if (state != "on" && state != "off" ) {
                    throw OML_Error(OML_ERR_FUNCSWITCH);
                }
                else {
                    cm->box(cm->gca(), state);
                }
            }
            else if (inputs[0].IsScalar()) {
                cm->box(inputs[0].Scalar());
            }
            else {
                throw OML_Error(OML_ERR_OPTIONVAL, 1, OML_VAR_PARAMETER);
	}
        }
        else if (inputs.size() == 2 && inputs[0].IsScalar() && inputs[1].IsString()) {
            cm->box(inputs[0].Scalar(), inputs[1].StringVal());
        }
        else {
            throw OML_Error(OML_ERR_NUMARGIN);
        }
		return true;
	}

	bool colorbar(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) 
    {
        std::unique_ptr<ColorbarData> cd = dsm.getColorbarData(inputs);
        double axesHandle = cm->gca();
        if (!inputs.empty() && inputs[0].IsScalar())
            axesHandle = inputs[0].Scalar();

        int nargout = eval.GetNargoutValue();
        if ((nargout > 0 && inputs.empty()) || !inputs.empty())
        {
            cd->toggleVisibility = false;
            cd->visible = true;
        }

        double lh = cm->colorbar(axesHandle, cd);
        
        if (nargout == 1)
        {
            outputs.push_back(lh);
        }
        else if (nargout == 2)
        {
            std::vector<double> r = cm->colorbarRange(axesHandle);
            outputs.insert(outputs.end(), r.begin(), r.end());
        }
        return true;
	}


    bool colormap(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        
        if (inputs.size() > 2)             
            throw OML_Error(OML_ERR_NUMARGIN);

        if (inputs.size() == 0) {
            double h = cm->gca();
            Currency cmap = cm->colormap(h);
            outputs.push_back(cmap);
        }
        else if (inputs.size() == 1) {
            if (inputs[0].IsScalar()) {
                double h = inputs[0].Scalar();
                Currency cmap = cm->colormap(h);
                outputs.push_back(cmap);
            }
            else if (inputs[0].IsMatrix()) {
                cm->colormap(cm->gca(), inputs[0]);
            }
            else
                throw OML_Error(OML_ERR_NUMARGIN);
        }
        else {
            if (inputs[0].IsScalar() && inputs[1].IsMatrix()) {
                cm->colormap(inputs[0].Scalar(), inputs[1]);
            }
            else
                throw OML_Error(OML_ERR_NUMARGIN);
        }
        return true;
    }

    bool plotyy(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        
        if (inputs.size() != 4)
            throw OML_Error(OML_ERR_NUMARGIN);

        vector<LineData> vld = dsm.getLineData(inputs);
        outputs = cm->plotyy(vld);
        
        return true;
    }

	bool drawnow(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        if (!inputs.empty())
            throw OML_Error(OML_ERR_NUMARGIN);
        cm->drawnow();
		return true;
	}

    bool findobj(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        std::unique_ptr<QueryData> data = dsm.getQueryData(inputs);
        std::vector<double> res = cm->findobj(data);
        outputs.push_back(res);
        return true;
    }

    bool deletefun(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        if (inputs.empty())
        {
            throw OML_Error(OML_ERR_NUMARGIN);
        }
        int nargin = static_cast<int>(inputs.size());

        std::vector<Currency> fileNames;      // These need to be deleted with BuiltInFuncsSystem::Delete
        std::vector<Currency> handles;    // Handle varnames to clear

        fileNames.reserve(nargin);
        handles.reserve(nargin);

        for (int i = 0; i < nargin; ++i)
        {
            const Currency& cur = inputs[i];
            if (cur.IsString())
            {
                fileNames.push_back(cur);
                continue;
            }
            else if (cur.IsScalar())
            {
                double handle = cur.Scalar();
                if (cm->deleteHandle(handle))
                {
                    handles.push_back(cur.GetOutputName());
                }
                else
                    fileNames.push_back(cur);
            }
        }
        if (!handles.empty())
        {
            eval.CallFunction("clear", handles);
        }
        if (!fileNames.empty())
        {
            BuiltInFuncsSystem::Delete(eval, fileNames, outputs);
        }
        return true;
    }

    bool view(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        if (inputs.size() == 0) {
            outputs.push_back(Currency());
            return true;
        }
        else {
            double axesHandle = -1;
            std::vector<double> viewVal;
            std::vector<Currency>::const_iterator it = inputs.cbegin();
            int pos = 1;
            for (; it != inputs.cend(); ++it)
            {
                if (it->IsScalar())
                {
                    if (cm->isAxes(it->Scalar()))
                        axesHandle = it->Scalar();
                    else
                        viewVal.push_back(it->Scalar());
                }
                else if (it->IsMatrix())
                {
                    const hwMatrix* mat = it->Matrix();
                    int c = mat->Size();
                    if (c == 2 || c == 3 || c == 16)
                    {
                        for (int i = 1; i <= c; i++)
                            viewVal.push_back((*mat)(i - 1));
                    }
                    else
                    {
                        throw OML_Error(OML_ERR_OPTIONVAL, pos, OML_VAR_PARAMETER);
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_OPTIONVAL, pos, OML_VAR_PARAMETER);
                }
                ++pos;
            }

            if (viewVal.empty()) {
                outputs.push_back(Currency());
            }
            else {
                cm->view(axesHandle, viewVal);
            }
            return true;
        }
        return true;
    }

    bool xline(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        if (inputs.size() == 0)
            throw OML_Error(OML_ERR_NUMARGIN);

        vector<LineData> vld = dsm.getXLineData(inputs);
        if (vld.empty())
            throw OML_Error(OML_ERR_PLOT_UNKNOWN_ERROR);

        double axis = (vld.front().parent > 0) ? vld.front().parent : axis = cm->gca();

        bool hold = cm->ishold(axis);
        cm->hold(axis, true);
        
        vector<double> hs = cm->xline(vld);
        outputs.push_back(hs);

        cm->hold(axis, hold);

        return true;
    }

    bool yline(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        if (inputs.size() == 0)
            throw OML_Error(OML_ERR_NUMARGIN);

        vector<LineData> vld = dsm.getXLineData(inputs);
        if (vld.empty())
            throw OML_Error(OML_ERR_PLOT_UNKNOWN_ERROR);

        double axis = (vld.front().parent > 0) ? vld.front().parent : axis = cm->gca();

        bool hold = cm->ishold(axis);
        cm->hold(axis, true);

        vector<double> hs = cm->yline(vld);
        outputs.push_back(hs);

        cm->hold(axis, hold);
        return true;
    }

    bool waterfall(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        vector<LineData> vld = dsm.getSurfData(inputs);
        vector<double> hs = cm->waterfall(vld);
        outputs.push_back(hs);
        return true;
    }

    bool ellipse(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        vector<LineData> vld = dsm.getShapeData(inputs, "ellipse");
        vector<double> hs = cm->ellipse(vld);
        outputs.push_back(hs);
        return true;
    }

    bool rectangle(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        vector<LineData> vld = dsm.getShapeData(inputs, "rectangle");
        vector<double> hs = cm->rectangle(vld);
        outputs.push_back(hs);
        return true;
    }

	bool fanplot(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        vector<LineData> vld = dsm.getSurfData(inputs);
        outputs = cm->fanplot(vld);
		return true;
	}

    bool pcolor(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        vector<LineData> vld = dsm.getSurfData(inputs);
        vector<double> hs = cm->pcolor(vld);
        outputs.push_back(hs);
        return true;
    }

    bool patch(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        vector<LineData> vld = dsm.getPatchData(inputs);
        vector<double> hs = cm->patch(vld);
        outputs.push_back(hs);
        return true;
    }

    bool stem3(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        vector<LineData> vld = dsm.getSurfData(inputs);
        vector<double> hs = cm->stem3(vld);
        outputs.push_back(hs);
        return true;
    }

    bool quiver(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        if (inputs.size() < 2)
            throw OML_Error(OML_ERR_NUMARGIN);
        std::vector<LineData> vld = dsm.getQuiverData(inputs);
        vector<double> hs = cm->quiver(vld);
        outputs.push_back(hs);
        return true;
    }

    bool findall(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        std::unique_ptr<QueryData> data = dsm.getQueryData(inputs);
        std::vector<double> res = cm->findobj(data, true);
        outputs.push_back(res);
        return true;
    }

    void UpdateDateticksForAxis(EvaluatorInterface eval, double axesHandle, const string& axis)
    {
        bool datetickEnabled = false;
        std::string datetickFMT;
        int datetickFMTIdx;
        cm->getAxisDatetickOptions(axesHandle, axis, datetickEnabled, datetickFMT, datetickFMTIdx);

        if (!datetickEnabled)
            return;
        
        Object* obj = cm->getObject(axesHandle);
        if (!obj->isAxes() && obj->getObjectType() != ObjectType::SECONDARYYAXIS)
            return;

        std::string prop("xtick");
        if (axis == "y" || obj->getObjectType() == ObjectType::SECONDARYYAXIS)
            prop = "ytick";

        Currency ticksCur = obj->getPropertyValue(prop).getCurrency();
        if (ticksCur.IsMatrix() && ticksCur.Matrix()->Size() > 0) {
            std::vector<Currency> ins, outs;
            ins.push_back(ticksCur);
            if (datetickFMT.empty())
                ins.push_back(datetickFMTIdx);
            else
                ins.push_back(datetickFMT);

            Currency out = eval.CallFunction("datestr", ins);
            // set the user ticks
            if (out.IsString())
            {
                std::string ret = out.StringVal();
                // if ret starts with '\n' remove it
                if (ret.find("\n") == 0)
                    ret = ret.substr(1);

                std::string split = "\n";
                std::string::size_type start = 0, end;
                std::vector<std::string> vs;
                do
                {
                    end = ret.find(split, start);
                    vs.push_back(ret.substr(start, end - start));
                    start = end + 1;
                } while (end != std::string::npos);
                if (!vs.empty() && vs.back() == "")
                    vs.pop_back();

                HML_CELLARRAY* cell = EvaluatorInterface::allocateCellArray(static_cast<int>(vs.size()), 1);
                for (int i = 0; i < vs.size(); i++)
                    (*cell)(i) = vs[i];

                unique_ptr<SetData> data(new SetData);
                data->handles.push_back(obj->getHandle());
                data->properties.push_back(prop + "label");
                data->values.push_back(cell);
                vector<string> wrn;
                cm->set(data, wrn);
            }
        }
        else
        {
            BuiltInFuncsUtils::SetWarning(eval, "property '" + prop + "' must be set to create date ticks;");
        }
    }

    bool datetick(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        
        double axesH = -1;
        bool firstInputIsAxis = false;
        if (inputs.empty() || inputs[0].IsInteger() || !inputs[0].IsScalar())
        {
            axesH = cm->gca();
        }
        else
        {
            axesH = inputs[0].Scalar();
            firstInputIsAxis = true;
        }
        std::string axis = "x";
        std::string dateFmt;
        int dateFmtIdx = 0;
        if (!inputs.empty())
        {
            int i = firstInputIsAxis ? 1 : 0;
            for (; i < inputs.size(); ++i)
            {
                if (inputs[i].IsInteger())
                {
                    dateFmtIdx = (int)inputs[i].Scalar();
                    dateFmt = "";
                    if (dateFmtIdx < 0 || dateFmtIdx>31)
                        throw OML_Error(OML_Error(OML_ERR_INTEGER).GetErrorMessage() + " in range [0, 31]", i + 1);
                }
                else if (inputs[i].IsString())
                {
                    std::string tmp = inputs[i].StringVal();
                    if (tmp == "x" || tmp == "y")
                    {
                        axis = tmp;
                    }
                    else if (tmp == "z")
                    {
                        throw OML_Error(OML_ERR_OPTIONVAL, i + 1);
                    }
                    else if (tmp == "keeplimits" || tmp == "keepticks")
                    {
                        BuiltInFuncsUtils::SetWarning(eval, "Option '" + tmp + "' is not supported");
                    }
                    else
                    {
                        dateFmt = tmp;
                        dateFmtIdx = -1;
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_STRING, i + 1);
                }
            }
        }

        // set the dateticks options
        cm->datetick(axesH, dateFmt, dateFmtIdx, axis);
        // update the tick labels
        UpdateDateticksForAxis(eval, axesH, axis);

        std::string warningMsg = cm->GetWarningString();
        if (!warningMsg.empty())
            BuiltInFuncsUtils::SetWarning(eval, warningMsg);
        return true;
    }

    bool bar3(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        vector<LineData> vld = dsm.getBar3Data(inputs);
        vector<double> hs = cm->bar3(vld);
        outputs.push_back(hs);
        return true;
    }

    bool copystyle(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [copystyle] is not supported in OpenMatrix");
        return false;
    }

    bool pastestyle(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [pastestyle] is not supported in OpenMatrix");
        return false;
    }

    bool hist3(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        
        if (inputs.empty())
            throw OML_Error(OML_ERR_NUMARGIN);

        if (!inputs[0].IsMatrix() ||
            (inputs[0].IsMatrix() && inputs[0].Matrix()->N() != 2))
        {
            throw OML_Error("must be a Nx2 matrix", 1);
        }
        
        
        // Create hist3 data
        std::vector<Currency> histOuts;
        vector<LineData> vld = dsm.getHist3Data(inputs, histOuts);

        int nargout = eval.GetNargoutValue();
        if (nargout == 1)
        {
            outputs.push_back(histOuts.front());
        }
        else if (nargout == 2)
        {
            outputs.push_back(histOuts[0]);
            outputs.push_back(histOuts[1]);
        }
        else if (nargout == 0 || nargout == 3)
        {
            // create the histogram plot
            std::vector<double> bar3Output = cm->bar3(vld);
            if (nargout == 3)
            {
                outputs.push_back(histOuts[0]);
                outputs.push_back(histOuts[1]);
                outputs.push_back(bar3Output.front());
            }

            // set the x categories and y categories
            HML_CELLARRAY* cell = histOuts[1].CellArray();

            unique_ptr<SetData> data(new SetData);
            data->handles.push_back(cm->gca());
            data->properties.push_back("xtick");
            data->values.push_back((*cell)(0));
            data->properties.push_back("ytick");
            data->values.push_back((*cell)(1));
            vector<string> wrn;
            cm->set(data, wrn);
        }
        return true;
    }

    bool triplot(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {

        if (inputs.empty())
            throw OML_Error(OML_ERR_NUMARGIN);

        if (inputs.size() < 3)
            throw OML_Error(OML_ERR_NUMARGIN);
        if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
            throw OML_Error(OML_ERR_SCALARMATRIX, 1);
        if (!inputs[1].IsMatrix())
            throw OML_Error(OML_ERR_MATRIX, 2);
        if (!inputs[2].IsMatrix())
            throw OML_Error(OML_ERR_MATRIX, 3);

        std::vector<LineData> vld = dsm.getTriplotData(inputs);
        outputs.push_back(cm->triplot(vld));
        return true;
    }

    bool trimesh(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {

        if (inputs.size() < 3)
            throw OML_Error(OML_ERR_NUMARGIN);
        if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
            throw OML_Error(OML_ERR_SCALARMATRIX, 1);
        if (!inputs[1].IsMatrix())
            throw OML_Error(OML_ERR_MATRIX, 2);
        if (!inputs[2].IsMatrix())
            throw OML_Error(OML_ERR_MATRIX, 3);

        std::vector<LineData> vld = dsm.getTrisurfData(inputs);
        // trimesh without z data falls back to triplot!
        if (vld.front().z.empty())
            outputs.push_back(cm->triplot(vld));
        else
            outputs.push_back(cm->trimesh(vld));
        return true;
    }

    bool trisurf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        
        if (inputs.size() < 3)
            throw OML_Error(OML_ERR_NUMARGIN);
        if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
            throw OML_Error(OML_ERR_SCALARMATRIX, 1);
        if (!inputs[1].IsMatrix())
            throw OML_Error(OML_ERR_MATRIX, 2);
        if (!inputs[2].IsMatrix())
            throw OML_Error(OML_ERR_MATRIX, 3);

        std::vector<LineData> vld = dsm.getTrisurfData(inputs);
        // trimesh without z data falls back to triplot!
        if (vld.front().z.empty())
            throw OML_Error(OML_ERR_NUMARGIN);
        outputs.push_back(cm->trisurf(vld));
        return true;
    }

    bool autumn(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [autumn] is not supported in OpenMatrix");
        return false;
    }

    bool bone(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [bone] is not supported in OpenMatrix");
        return false;
    }

    bool cividis(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [cividis] is not supported in OpenMatrix");
        return false;
    }

    bool cool(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [cool] is not supported in OpenMatrix");
        return false;
    }

    bool copper(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [copper] is not supported in OpenMatrix");
        return false;
    }

    bool cubehelix(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [cubehelix] is not supported in OpenMatrix");
        return false;
    }

    bool flag(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [flag] is not supported in OpenMatrix");
        return false;
    }

    bool gray(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [gray] is not supported in OpenMatrix");
        return false;
    }

    bool hot(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [hot] is not supported in OpenMatrix");
        return false;
    }

    bool hsv(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [hsv] is not supported in OpenMatrix");
        return false;
    }

    bool inferno(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [inferno] is not supported in OpenMatrix");
        return false;
    }

    bool jet(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [jet] is not supported in OpenMatrix");
        return false;
    }

    bool lines(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [lines] is not supported in OpenMatrix");
        return false;
    }

    bool magma(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [magma] is not supported in OpenMatrix");
        return false;
    }

    bool ocean(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [ocean] is not supported in OpenMatrix");
        return false;
    }

    bool pink(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [pink] is not supported in OpenMatrix");
        return false;
    }

    bool plasma(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [plasma] is not supported in OpenMatrix");
        return false;
    }

    bool prism(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [prism] is not supported in OpenMatrix");
        return false;
    }

    bool rainbow(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [rainbow] is not supported in OpenMatrix");
        return false;
    }

    bool spring(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [spring] is not supported in OpenMatrix");
        return false;
    }

    bool summer(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [summer] is not supported in OpenMatrix");
        return false;
    }

    bool twilight(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [twilight] is not supported in OpenMatrix");
        return false;
    }

    bool viridis(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [viridis] is not supported in OpenMatrix");
        return false;
    }

    bool white(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [white] is not supported in OpenMatrix");
        return false;
    }

    bool winter(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
        BuiltInFuncsUtils::SetWarning(eval, "Command [winter] is not supported in OpenMatrix");
        return false;
    }

#define TBOXVERSION 1.13
    extern "C" OMLPLOT_EXPORT
    double GetToolboxVersion(EvaluatorInterface eval){
        return TBOXVERSION;
    }

    extern "C" OMLPLOT_EXPORT
    int InitDll(EvaluatorInterface evl)
    {
        if (! GnuplotOutput::isGnuplotReady()){
            evl.RegisterBuiltInFunction("plot", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("bar", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("hist", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("area", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("polar", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("line", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("scatter", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("plot3", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("fill", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("stem", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("scatter3", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("surf", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("mesh", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("contour3", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("contour", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("loglog", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("semilogx", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("semilogy", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("figure", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("subplot", oml_doNothing, FunctionMetaData(3, 1, "Plotting"));
            evl.RegisterBuiltInFunction("axes", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("get", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("set", oml_doNothing, FunctionMetaData(3, 0, "Plotting"));
            evl.RegisterBuiltInFunction("gcf", oml_doNothing, FunctionMetaData(0, 1, "Plotting"));
            evl.RegisterBuiltInFunction("gca", oml_doNothing, FunctionMetaData(0, 1, "Plotting"));
            evl.RegisterBuiltInFunction("clf", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("cla", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("close", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("hold", oml_doNothing, FunctionMetaData(-1, 0, "Plotting"));
            evl.RegisterBuiltInFunction("ishold", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("ishandle", oml_doNothing, FunctionMetaData(1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("isfigure", oml_doNothing, FunctionMetaData(1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("isaxes", oml_doNothing, FunctionMetaData(1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("grid", oml_doNothing, FunctionMetaData(-1, 0, "Plotting"));
            evl.RegisterBuiltInFunction("title", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("xlabel", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("ylabel", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("zlabel", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("axis", oml_doNothing, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("xlim", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("ylim", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("zlim", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("legend", oml_doNothing, FunctionMetaData(-1, 0, "Plotting"));
            evl.RegisterBuiltInFunction("text", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("saveas", oml_doNothing, FunctionMetaData(-1, 0, "Plotting"));
            evl.RegisterBuiltInFunction("view", oml_doNothing, FunctionMetaData(-1, 0, "Plotting"));
            evl.RegisterBuiltInFunction("colorbar", oml_doNothing, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("plotyy", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("colormap", oml_doNothing, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("waterfall", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("box", oml_doNothing, FunctionMetaData(-1, 0, "Plotting"));
            evl.RegisterBuiltInFunction("drawnow", oml_doNothing, FunctionMetaData(0, 0, "Plotting"));
            evl.RegisterBuiltInFunction("fanplot", oml_doNothing, FunctionMetaData(3, -1, "Plotting"));
            evl.RegisterBuiltInFunction("findobj", oml_doNothing, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("xline", oml_doNothing, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("yline", oml_doNothing, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("ellipse", oml_doNothing, FunctionMetaData(-1, 0, "Plotting"));
            evl.RegisterBuiltInFunction("rectangle", oml_doNothing, FunctionMetaData(-1, 0, "Plotting"));
            evl.RegisterBuiltInFunction("patch", oml_doNothing, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("pcolor", oml_doNothing, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("stem3", oml_doNothing, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("quiver", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("findall", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("datetick", oml_doNothing, FunctionMetaData(-4, 0, "Plotting"));
            evl.RegisterBuiltInFunction("bar3", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("hist3", oml_doNothing, FunctionMetaData(-2, 1, "Plotting"));
            evl.RegisterBuiltInFunction("triplot", oml_doNothing, FunctionMetaData(-3, -1, "Plotting"));
            evl.RegisterBuiltInFunction("trimesh", oml_doNothing, FunctionMetaData(-3, -1, "Plotting"));
            evl.RegisterBuiltInFunction("trisurf", oml_doNothing, FunctionMetaData(-3, -1, "Plotting"));

            // Not yet supported commands
            evl.RegisterBuiltInFunction("copystyle", copystyle, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("pastestyle", pastestyle, FunctionMetaData(-1, -1, "Plotting"));
            evl.RegisterBuiltInFunction("autumn", autumn, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("bone", bone, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("cividis", cividis, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("cool", cool, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("copper", copper, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("cubehelix", cubehelix, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("flag", flag, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("gray", gray, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("hot", hot, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("hsv", hsv, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("inferno", inferno, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("jet", jet, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("lines", lines, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("magma", magma, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("ocean", ocean, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("pink", pink, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("plasma", plasma, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("prism", prism, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("rainbow", rainbow, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("spring", spring, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("summer", summer, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("twilight", twilight, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("viridis", viridis, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("white", white, FunctionMetaData(-1, 1, "Plotting"));
            evl.RegisterBuiltInFunction("winter", winter, FunctionMetaData(-1, 1, "Plotting"));

            return 0;
        }

        evl.RegisterBuiltInFunction("plot", plot, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("bar", bar, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("hist", hist, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("area", area, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("polar", polar, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("line", line, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("scatter", scatter, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("plot3", plot3, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("fill", fill, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("stem", stem, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("scatter3", scatter3, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("surf", surf, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("mesh", mesh, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("contour3", contour3, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("contour", contour, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("loglog", loglog, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("semilogx", semilogx, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("semilogy", semilogy, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("figure", figure, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("subplot", subplot, FunctionMetaData(3, 1, "Plotting"));
        evl.RegisterBuiltInFunction("axes", axes, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("get", get, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("set", set, FunctionMetaData(3, 0, "Plotting"));
        evl.RegisterBuiltInFunction("gcf", gcf, FunctionMetaData(0, 1, "Plotting"));
        evl.RegisterBuiltInFunction("gca", gca, FunctionMetaData(0, 1, "Plotting"));
        evl.RegisterBuiltInFunction("clf", clf, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("cla", cla, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("close", close, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("hold", hold, FunctionMetaData(-1, 0, "Plotting"));
        evl.RegisterBuiltInFunction("ishold", ishold, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("ishandle", ishandle, FunctionMetaData(1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("isfigure", isfigure, FunctionMetaData(1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("isaxes", isaxes, FunctionMetaData(1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("grid", grid, FunctionMetaData(-1, 0, "Plotting"));
        evl.RegisterBuiltInFunction("title", title, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("xlabel", xlabel, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("ylabel", ylabel, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("zlabel", zlabel, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("axis", axis, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("xlim", xlim, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("ylim", ylim, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("zlim", zlim, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("legend", legend, FunctionMetaData(-1, 0, "Plotting"));
        evl.RegisterBuiltInFunction("text", text, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("saveas", saveas, FunctionMetaData(-1, 0, "Plotting"));
        evl.RegisterBuiltInFunction("view", view, FunctionMetaData(-1, 0, "Plotting"));
        evl.RegisterBuiltInFunction("box", box, FunctionMetaData(-1, 0, "Plotting"));
        evl.RegisterBuiltInFunction("colorbar", colorbar, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("plotyy", plotyy, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("colormap", colormap, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("delete", deletefun, FunctionMetaData(-1, 0, "System"));
        evl.RegisterBuiltInFunction("ellipse", ellipse, FunctionMetaData(-1, 0, "Plotting"));
        evl.RegisterBuiltInFunction("rectangle", rectangle, FunctionMetaData(-1, 0, "Plotting"));
        evl.RegisterBuiltInFunction("waterfall", waterfall, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("patch", patch, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("pcolor", pcolor, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("drawnow", drawnow, FunctionMetaData(0, 0, "Plotting"));
        evl.RegisterBuiltInFunction("fanplot", fanplot, FunctionMetaData(3, -1, "Plotting"));
        evl.RegisterBuiltInFunction("findobj", findobj, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("xline", xline, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("yline", yline, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("stem3", stem3, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("quiver", quiver, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("findall", findall, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("datetick", datetick, FunctionMetaData(-4, 0, "Plotting"));
        evl.RegisterBuiltInFunction("bar3", bar3, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("hist3", hist3, FunctionMetaData(-2, 1, "Plotting"));
        evl.RegisterBuiltInFunction("triplot", triplot, FunctionMetaData(-3, -1, "Plotting"));
        evl.RegisterBuiltInFunction("trimesh", trimesh, FunctionMetaData(-3, -1, "Plotting"));
        evl.RegisterBuiltInFunction("trisurf", trisurf, FunctionMetaData(-3, -1, "Plotting"));

#ifdef _DEBUG
        evl.RegisterBuiltInFunction("dump", dump, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("out", out, FunctionMetaData(-1, 1, "Plotting"));
#endif
		// Not yet supported commands
        evl.RegisterBuiltInFunction("copystyle", copystyle, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("pastestyle", pastestyle, FunctionMetaData(-1, -1, "Plotting"));
        evl.RegisterBuiltInFunction("autumn", autumn, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("bone", bone, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("cividis", cividis, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("cool", cool, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("copper", copper, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("cubehelix", cubehelix, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("flag", flag, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("gray", gray, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("hot", hot, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("hsv", hsv, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("inferno", inferno, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("jet", jet, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("lines", lines, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("magma", magma, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("ocean", ocean, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("pink", pink, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("plasma", plasma, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("prism", prism, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("rainbow", rainbow, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("spring", spring, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("summer", summer, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("twilight", twilight, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("viridis", viridis, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("white", white, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("winter", winter, FunctionMetaData(-1, 1, "Plotting"));

        return 0;
    }

}
