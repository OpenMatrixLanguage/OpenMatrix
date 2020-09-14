/**
* @file PlotToolbox.cxx
* @date March 2017
* Copyright (C) 2017-2020 Altair Engineering, Inc.  
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
#include <memory>

using namespace std;

#define CATCH_ALL_PLOT_ERROR throw OML_Error(OML_ERR_PLOT_UNKNOWN_ERROR);

namespace omlplot{

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

    bool waterfall(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        /*
        vector<LineData> vld = dsm.getSurfData(inputs);
        vector<double> hs = cm->waterfall(vld);
        outputs.push_back(hs);
        */
        outputs.push_back(1.2);
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
                Property p = cm->getObjectProperty(h, props[0]);
				if (p.getType() == UNSUPPORTED)
				{
					BuiltInFuncsUtils::SetWarning(eval, "Property ["+props[0]+"] is not supported in OpenMatrix");
					return false;
				}
				outputs.push_back(castValue(p.getValue()));
            } else {
                StructData *sd = new StructData();
                vector<string>::iterator it = props.begin();
                for (; it != props.end(); ++it){
                    string pname = *it;
                    Property p = cm->getObjectProperty(h, pname);
					if (p.getType() == UNSUPPORTED)
						continue;

                    sd->SetValue(0, 0, pname,
                                 castValue(p.getValue()) );
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
        return true;
    }

    bool figure(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        unique_ptr<FigureData> fd = dsm.getFigureData(inputs);
        double h = cm->figure(fd);
        outputs.push_back(h);
        return true;
    }

    bool close(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        int f = 0;
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
        cm->close(f);
        return true;
    }

    bool gcf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        double h = cm->gcf();
        outputs.push_back(h);
        return true;
    }

    bool clf(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        int f = 0;
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
        int h = cm->clf(f);
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
            if ( ! inputs[0].IsPositiveInteger() ){
                throw OML_Error(OML_ERR_SCALAR, 1);
            }
            if ( ! inputs[1].IsPositiveInteger() ){
                throw OML_Error(OML_ERR_SCALAR, 2);
            }
            if ( ! inputs[2].IsPositiveInteger() ){
                throw OML_Error(OML_ERR_SCALAR, 3);
            }
        } else if ( inputs.size() == 1 ){
            if ( ! inputs[0].IsPositiveInteger() ){
                throw OML_Error(OML_ERR_SCALAR, 1);
            }
        } else {
            throw OML_Error(OML_ERR_NUMARGIN);
        }
        int row = 1, col = 1, active = 1;
        size_t size = inputs.size();
        if (size == 3){
            row = (int)inputs[0].Scalar();
            col = (int)inputs[1].Scalar();
            active = (int)inputs[2].Scalar();
        } else if (size == 1){
            int rcn = (int)inputs[0].Scalar();
            if (rcn < 111){
                throw OML_Error(OML_ERR_OPTION);
            }
            row = rcn / 100;
            col = (rcn / 10) % 10;
            active = rcn % 10;
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
        cm->legend(ld);
        return true;
    }

    bool text(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        std::unique_ptr<TextData> td = dsm.getTextData(inputs);
        vector<double> res = cm->text(td);
        outputs.push_back(res);
        return true;
    }

    bool saveas(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs){
        if (inputs.size() < 2) {
            throw OML_Error(OML_ERR_NUMARGIN);
        }
        if ( ! inputs[0].IsScalar() ){
            throw OML_Error(OML_ERR_SCALAR, 1);
        }
        if ( ! inputs[1].IsString() ){
            throw OML_Error(OML_ERR_STRING, 2);
        }
        if ( inputs.size() == 3 &&
             (! inputs[2].IsString()) ){
            throw OML_Error(OML_ERR_STRING, 3);
        }

        double handle = inputs[0].Scalar();
        string filename, fmt;
        filename = inputs[1].StringVal();
        if (inputs.size() == 3){
            fmt = inputs[2].StringVal();
        }
        try {
            cm->saveas(handle, filename, fmt);
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
		BuiltInFuncsUtils::SetWarning(eval, "Command [box] is not supported in OpenMatrix");
		return false;
	}

	bool colorbar(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
		BuiltInFuncsUtils::SetWarning(eval, "Command [colorbar] is not supported in OpenMatrix");
		return false;
	}

	bool drawnow(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
		BuiltInFuncsUtils::SetWarning(eval, "Command [drawnow] is not supported in OpenMatrix");
		return false;
	}

	bool fanplot(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
		BuiltInFuncsUtils::SetWarning(eval, "Command [fanplot] is not supported in OpenMatrix");
		return false;
	}

	bool getmousepos(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
		BuiltInFuncsUtils::SetWarning(eval, "Command [getmousepos] is not supported in OpenMatrix");
		return false;
	}

	bool imagesc(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
		BuiltInFuncsUtils::SetWarning(eval, "Command [imagesc] is not supported in OpenMatrix");
		return false;
	}

	bool plotyy(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
		BuiltInFuncsUtils::SetWarning(eval, "Command [plotyy] is not supported in OpenMatrix");
		return false;
	}

	bool uiresume(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
		BuiltInFuncsUtils::SetWarning(eval, "Command [uiresume] is not supported in OpenMatrix");
		return false;
	}

	bool uiwait(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs) {
		BuiltInFuncsUtils::SetWarning(eval, "Command [uiwait] is not supported in OpenMatrix");
		return false;
	}

#define TBOXVERSION 2019.0
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
            evl.RegisterBuiltInFunction("waterfall", oml_doNothing, FunctionMetaData(-1, 1, "Plotting"));
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
			// Not yet supported commands
			evl.RegisterBuiltInFunction("box", box, FunctionMetaData(-1, 0, "Plotting"));
			evl.RegisterBuiltInFunction("colorbar", colorbar, FunctionMetaData(-1, -1, "Plotting"));
			evl.RegisterBuiltInFunction("drawnow", drawnow, FunctionMetaData(0, 0, "Plotting"));
			evl.RegisterBuiltInFunction("fanplot", fanplot, FunctionMetaData(3, -1, "Plotting"));
			evl.RegisterBuiltInFunction("getmousepos", getmousepos, FunctionMetaData(0, 2, "Plotting"));
			evl.RegisterBuiltInFunction("imagesc", imagesc, FunctionMetaData(-1, 1, "Plotting"));
			evl.RegisterBuiltInFunction("plotyy", plotyy, FunctionMetaData(-1, 1, "Plotting"));
			evl.RegisterBuiltInFunction("uiresume", uiresume, FunctionMetaData(-1, 1, "Plotting"));
			evl.RegisterBuiltInFunction("uiwait", uiwait, FunctionMetaData(-1, 1, "Plotting"));

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
        evl.RegisterBuiltInFunction("waterfall", waterfall, FunctionMetaData(-1, 1, "Plotting"));
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
        evl.RegisterBuiltInFunction("view", oml_doNothing, FunctionMetaData(-1, 0, "Plotting"));
#ifdef _DEBUG
        evl.RegisterBuiltInFunction("dump", dump, FunctionMetaData(-1, 1, "Plotting"));
        evl.RegisterBuiltInFunction("out", out, FunctionMetaData(-1, 1, "Plotting"));
#endif
		// Not yet supported commands
		evl.RegisterBuiltInFunction("box", box, FunctionMetaData(-1, 0, "Plotting"));
		evl.RegisterBuiltInFunction("colorbar", colorbar, FunctionMetaData(-1, -1, "Plotting"));
		evl.RegisterBuiltInFunction("drawnow", drawnow, FunctionMetaData(0, 0, "Plotting"));
		evl.RegisterBuiltInFunction("fanplot", fanplot, FunctionMetaData(3, -1, "Plotting"));
		evl.RegisterBuiltInFunction("getmousepos", getmousepos, FunctionMetaData(0, 2, "Plotting"));
		evl.RegisterBuiltInFunction("imagesc", imagesc, FunctionMetaData(-1, 1, "Plotting"));
		evl.RegisterBuiltInFunction("plotyy", plotyy, FunctionMetaData(-1, 1, "Plotting"));
		evl.RegisterBuiltInFunction("uiresume", uiresume, FunctionMetaData(-1, 1, "Plotting"));
		evl.RegisterBuiltInFunction("uiwait", uiwait, FunctionMetaData(-1, 1, "Plotting"));
        return 0;
    }

}
