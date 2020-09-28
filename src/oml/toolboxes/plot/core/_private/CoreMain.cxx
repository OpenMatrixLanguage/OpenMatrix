/**
* @file CoreMain.cxx
* @date May 2018
* Copyright (C) 2018-2020 Altair Engineering, Inc.  
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

#include "CoreMain.h"
#include <algorithm>
#include "OML_Error.h"

namespace omlplot{

    CoreMain *CoreMain::m_instance = nullptr;

    CoreMain::CoreMain()
        :root(new Root){
    }

    CoreMain::~CoreMain(){
        root.reset(nullptr);
    }

    CoreMain* CoreMain::getInstance(){
        if (! m_instance){
            m_instance = new CoreMain();
        }
        return m_instance;
    }

    void CoreMain::releaseInstance(){
        delete m_instance;
    }

    bool CoreMain::isHandle(double h){
        return root->isHandle(h);
    }

    bool CoreMain::isFigure(double h){
        return root->isFigure(h);
    }

    bool CoreMain::isAxes(double h){
        return root->isAxes(h);
    }


    Object *CoreMain::getObject(double h){
        return root->getObject(h);
    }

    Property CoreMain::getObjectProperty(double h, string name){
        Object *o = getObject(h);
        return o->getProperty(name);
    }

    vector<string> CoreMain::getObjectPropertyNames(double h){
        Object *o = getObject(h);
        return o->getPropertyNames();
    }

    double CoreMain::figure(unique_ptr<FigureData> &fd){
        size_t idCount = fd->fids.size();

        // get/alloc the figure first
        Figure *f = nullptr;
        double fHandle = 0;
        if (idCount == 0){
            f = allocObject<Figure>();
            root->addChild(f);
        } else if (idCount == 1){
            int fid = int(fd->fids[0]);
            if (isFigure(fid) ){
                f = dynamic_cast<Figure *>(getObject(fid));
            } else {
                if (fid == -1){
                    f = allocObject<Figure>();
                } else {
                    f = allocObject<Figure>(fid);
                }
                root->addChild(f);
            }
        } else {
        }

        // set properties
        if (! f){
            // TODO throw
        }
        fHandle = f->getHandle();
        root->setPropertyValue("currentfigure", fHandle);
        return fHandle;
    }

    double CoreMain::axes(unique_ptr<AxesData> &data){
        Object *figure = getObject(gcf());
        size_t count = data->handles.size();
        Axes *a = nullptr;
        if (count == 0){
            a = allocObject<Axes>();
            figure->addChild(a);
        } else if (count == 1) {
            double handle = data->handles[0];
            if (isAxes(handle)){
                a = dynamic_cast<Axes *>(getObject(handle));
            } else {
                throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
            }
        } else {
        }
        if (! a){
            throw;
        }

        double h = a->getHandle();
        figure->setPropertyValue("currentaxes", h);
        root->repaint();
        return h;
    }

    double CoreMain::subplot(int row, int col, int active){
        Object *figure = getObject(gcf());
        Axes *a = allocObject<Axes>();

        int r = (active - 1) / col;
        int c = (active - 1) % col;
        int rr = row - 1 - r;
        double x = double(c) / col;
        double y = double(rr) / row;
        double width = 1.0 / col;
        double height = 1.0 / row;
        vector<double> pos;
        pos.push_back(x); pos.push_back(y);
        pos.push_back(width); pos.push_back(height);
        a->setPropertyValue("pos", pos);

        figure->addChild(a);
        double h = a->getHandle();
        figure->setPropertyValue("currentaxes", h);
        root->repaint();
        return h;
    }

    vector<double> CoreMain::plot(vector<LineData> &ldVec){
        return _T_2D_PLOT<Line>(ldVec);
    }

    vector<double> CoreMain::bar(vector<LineData> &ldVec){
        return _T_2D_PLOT<HggroupBar>(ldVec);        
    }

    vector<double> CoreMain::fill(vector<LineData> &ldVec){
        return _T_2D_PLOT<Fill>(ldVec);
    }
    vector<double> CoreMain::hist(vector<LineData> &ldVec){
        return _T_2D_PLOT<Hist>(ldVec);        
    }

    vector<double> CoreMain::area(vector<LineData> &ldVec){
        return _T_2D_PLOT<Area>(ldVec);        
    }

    vector<double> CoreMain::polar(vector<LineData> &ldVec){
        return _T_2D_PLOT<Polar>(ldVec);        
    }

    vector<double> CoreMain::line(vector<LineData> &ldVec){
        double ah = 0;
        if (ldVec.size() > 0){
            ah = ldVec[0].parent;
        }
        if (! isAxes(ah)){
            ah = gca();
        }
        Axes *axes = dynamic_cast<Axes *>(getObject(ah));
        vector<double> res;
        Drawable *pLine = nullptr;
        vector<LineData>::iterator it = ldVec.begin();
        for (; it != ldVec.end(); ++it){
            LineData ld = *it;
            try {
                pLine = allocObject<Line>();
                pLine->init(ld);
                axes->addChild(pLine);
                res.push_back(pLine->getHandle());
            } catch (OML_Error &e){
                delete pLine;
                throw e;
            }
        }
        axes->repaint();
        return res;
    }

    vector<double> CoreMain::scatter(vector<LineData> &ldVec){
        return _T_2D_PLOT<HggroupScatter>(ldVec);        
    }

    vector<double> CoreMain::plot3(vector<LineData> &ldVec){
        return _T_3D_PLOT<Line3>(ldVec);        
    }

    vector<double> CoreMain::scatter3(vector<LineData> &ldVec){
        return _T_3D_PLOT<HggroupScatter3>(ldVec);
    }
    vector<double> CoreMain::surf(vector<LineData> &ldVec){
        return _T_3D_PLOT<Surface>(ldVec);        
    }

    vector<double> CoreMain::mesh(vector<LineData> &ldVec){
        return _T_3D_PLOT<Mesh>(ldVec);        
    }

    vector<double> CoreMain::contour3(vector<LineData> &ldVec){
        return _T_3D_PLOT<Contour3>(ldVec);        
    }

    vector<double> CoreMain::contour(vector<LineData> &ldVec){
        return _T_3D_PLOT<Contour>(ldVec);        
    }

    vector<double> CoreMain::stem(vector<LineData> &ldVec){
        return _T_2D_PLOT<Stem>(ldVec);
    }

    vector<double> CoreMain::loglog(vector<LineData> &ldVec){
        return _T_2D_PLOT<Loglog>(ldVec);
    }

    vector<double> CoreMain::semilogx(vector<LineData> &ldVec){
        return _T_2D_PLOT<Semilogx>(ldVec);
    }

    vector<double> CoreMain::semilogy(vector<LineData> &ldVec){
        return _T_2D_PLOT<Semilogy>(ldVec);
    }

    void CoreMain::set(unique_ptr<SetData> &data, vector<string>& notSupported)
	{
        vector<double> hs =data->handles;
        vector<double>::iterator it = hs.begin();
        for ( ; it != hs.end(); ++it)
		{
            double handle = *it;
            Object *po = getObject(handle);
            size_t propCount = data->properties.size();
            for (int index = 0; index < propCount; index++)
			{
                bool ret = po->setPropertyValue(data->properties[index], data->values[index]);
				if (!ret)
					notSupported.push_back(data->properties[index]);
            }
        }
        root->update();
    }

    double CoreMain::gcf(){
        double h = 0;
        h = root->getPropertyValue("currentfigure").Scalar();
        if (h == 0) {
            unique_ptr<FigureData> fd(new FigureData);
            h = figure(fd);
        }
        return h;
            
    }

    double CoreMain::gca(){
        double h = 0;
        Object *figure = getObject(gcf() );
        h = figure->getPropertyValue("currentaxes").Scalar();
        if (h == 0) {
            unique_ptr<AxesData> data(new AxesData);
            h = axes(data);
        }        
        return h;
    }

    void CoreMain::close(int f){
        if (! isFigure(f)){
            if (! (f == -1)){
                throw OML_Error(OML_ERR_PLOT_INVALID_FIGURE_HANDLE);
            }
        }
        if (f == -1){
            vector<double> fs = root->getPropertyValue("children").Vector();
            vector<double>::iterator it = fs.begin();
            for (; it != fs.end(); ++it){
                Object *o = getObject(*it);
                delete o;
            }
            root->setPropertyValue("currentfigure", 0);
            return;
        }
        Object *o = getObject(f);
        delete o;
        vector<double> fs = root->getPropertyValue("children").Vector();
        size_t count = fs.size();
        if (count == 0){
            root->setPropertyValue("currentfigure", 0);
        } else {
            root->setPropertyValue("currentfigure", fs[count - 1]);
        }
    }

    int CoreMain::clf(int f){
        if (! isFigure(f)){
            throw OML_Error(OML_ERR_PLOT_INVALID_FIGURE_HANDLE);
        }
        Figure *fig = dynamic_cast<Figure *>(getObject(f));
        fig->clear();
        return (int)fig->getHandle();
    }

    double CoreMain::cla(double a){
        if (! isAxes(a)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *axes = dynamic_cast<Axes *>(getObject(a));
        axes->clear();
        return axes->getHandle();
    }

    void CoreMain::hold(double axes, bool hold){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *a = dynamic_cast<Axes *>(getObject(axes));
        a->hold(hold);
    }

    bool CoreMain::ishold(double axes){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *a = dynamic_cast<Axes *>(getObject(axes));
        return a->ishold();
    }

    void CoreMain::grid(double axes){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        string xgrid = pAxes->getPropertyValue("xgrid").StringVal();
        string ygrid = pAxes->getPropertyValue("ygrid").StringVal();
        string zgrid = pAxes->getPropertyValue("zgrid").StringVal();
        string status = "off";
        if (xgrid == "off" || ygrid == "off" || zgrid == "off"){
            status = "on";
        }
        pAxes->setPropertyValue("xgrid", status);
        pAxes->setPropertyValue("ygrid", status);
        pAxes->setPropertyValue("zgrid", status);
        pAxes->repaint();
    }

    void CoreMain::grid(double axes, string status){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        if (status != "on" && status != "off"){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        pAxes->setPropertyValue("xgrid", status);
        pAxes->setPropertyValue("ygrid", status);
        pAxes->setPropertyValue("zgrid", status);
        pAxes->repaint();
    }

    void CoreMain::minorgrid(double axes){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        string xgrid = pAxes->getPropertyValue("xminorgrid").StringVal();
        string ygrid = pAxes->getPropertyValue("yminorgrid").StringVal();
        string zgrid = pAxes->getPropertyValue("zminorgrid").StringVal();
        string status = "off";
        if (xgrid == "off" || ygrid == "off" || zgrid == "off"){
            status = "on";
        }
        pAxes->setPropertyValue("xminorgrid", status);
        pAxes->setPropertyValue("yminorgrid", status);
        pAxes->setPropertyValue("zminorgrid", status);
        pAxes->repaint();
    }

    void CoreMain::minorgrid(double axes, string status){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        if (status != "on" && status != "off"){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        pAxes->setPropertyValue("xminorgrid", status);
        pAxes->setPropertyValue("yminorgrid", status);
        pAxes->setPropertyValue("zminorgrid", status);
        pAxes->repaint();
    }

    double CoreMain::title(double axes, string str){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        VALUETYPE c = pAxes->getPropertyValue("title");
        Object *text = (Object *)c.BoundObject();
        text->setPropertyValue("string", str);
        pAxes->repaint();
        return text->getHandle();
    }

    double CoreMain::xlabel(double axes, string str){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        VALUETYPE c = pAxes->getPropertyValue("xlabel");
        Text *text = (Text *)c.BoundObject();
        text->setPropertyValue("string", str);
        pAxes->repaint();
        return text->getHandle();
    }

    double CoreMain::ylabel(double axes, string str){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        VALUETYPE c = pAxes->getPropertyValue("ylabel");
        Text *text = (Text *)c.BoundObject();
        text->setPropertyValue("string", str);
        pAxes->repaint();
        return text->getHandle();
    }

    double CoreMain::zlabel(double axes, string str){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        VALUETYPE c = pAxes->getPropertyValue("zlabel");
        Text *text = (Text *)c.BoundObject();
        text->setPropertyValue("string", str);
        pAxes->repaint();
        return text->getHandle();
    }

    vector<double> CoreMain::axis(double axes){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        vector<double> xlim = pAxes->getPropertyValue("xlim").Vector();
        vector<double> ylim = pAxes->getPropertyValue("ylim").Vector();
        vector<double> zlim = pAxes->getPropertyValue("zlim").Vector();
        vector<double> res;
        res.push_back(xlim[0]);res.push_back(xlim[1]);
        res.push_back(ylim[0]);res.push_back(ylim[1]);
        res.push_back(zlim[0]);res.push_back(zlim[1]);
        return res;
    }

    void CoreMain::axis(double axes, vector<double> data){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        vector<double> limits(2);
        limits[0] = data[0]; limits[1] = data[1];
        pAxes->setPropertyValue("xlim", limits);
        if (data.size() == 4){
            limits[0] = data[2]; limits[1] = data[3];
            pAxes->setPropertyValue("ylim", limits);
        }
        if (data.size() == 6){
            limits[0] = data[4]; limits[1] = data[5];
            pAxes->setPropertyValue("zlim", limits);
        }
        pAxes->setAxisNeedRepaint(true);
        pAxes->repaint();
    }

    vector<double> CoreMain::xlim(double axes){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        return pAxes->getPropertyValue("xlim").Vector();
    }

    void CoreMain::xlim(double axes, vector<double> limits){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        pAxes->setPropertyValue("xlim", limits);
        pAxes->setAxisNeedRepaint(true);
        pAxes->repaint();
    }

    vector<double> CoreMain::ylim(double axes){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        return pAxes->getPropertyValue("ylim").Vector();
    }

    void CoreMain::ylim(double axes, vector<double> limits){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        pAxes->setPropertyValue("ylim", limits);
        pAxes->setAxisNeedRepaint(true);
        pAxes->repaint();
    }

    vector<double> CoreMain::zlim(double axes){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        return pAxes->getPropertyValue("zlim").Vector();
    }

    void CoreMain::zlim(double axes, vector<double> limits){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        pAxes->setPropertyValue("zlim", limits);
        pAxes->setAxisNeedRepaint(true);
        pAxes->repaint();
    }

    void CoreMain::legend(unique_ptr<LegendData> &ld){
        Axes *pAxes = dynamic_cast<Axes *>(getObject(gca()));
        vector<string> ls = ld->legends;
        if (ls.size() == 0){
            string legend = pAxes->getPropertyValue("legend").StringVal();
            if (legend == "on"){
                pAxes->setPropertyValue("legend", string("off"));
            } else {
                pAxes->setPropertyValue("legend", string("on"));
            }
        } else if ((ls.size() == 1) && (ls[0] == "on" || ls[0] == "off")){
            if (ls[0] == "on"){
                pAxes->setPropertyValue("legend", string("on"));
            } else {            // off
                pAxes->setPropertyValue("legend", string("off"));
            }
        } else {
            vector<double> children = pAxes->getPropertyValue("children").Vector();
            int count = min((int)ls.size(), (int)children.size());
            Drawable *pLine = nullptr;
            for (int i = 0; i < count; i++){
                pLine = dynamic_cast<Drawable *>(getObject(children[i]));
                pLine->setPropertyValue("displayname", ls[i]);
            }
            pAxes->setPropertyValue("legend", string("on"));
        }
        pAxes->repaint();
    }

    vector<double> CoreMain::text(unique_ptr<TextData> &td){
        Axes *pAxes = dynamic_cast<Axes *>(getObject(gca()));
        vector<double> res;
        Text *pText = nullptr;

        size_t size = td->xpos.size();
        for (int i = 0; i < size; i++){
            try {
                pText = allocObject<Text>();
                pText->setPropertyValue("x", td->xpos[i]);
                pText->setPropertyValue("y", td->ypos[i]);
                pText->setPropertyValue("string", td->text[i]);
                pAxes->addChild(pText);
                res.push_back(pText->getHandle());
            } catch (OML_Error &e){
                delete pText;
                throw e;
            }
        }
        pAxes->repaint();
        return res;
    }

    void CoreMain::saveas(double h, string filename, string fmt){
        if (! isFigure(h)){
            throw OML_Error(OML_ERR_PLOT_INVALID_FIGURE_HANDLE);
        }

        Figure *figure = dynamic_cast<Figure*>(getObject(h));
        figure->saveas(filename, fmt);
        return;
    }


    void CoreMain::dump(){
        root->dump(root.get(), 0);
    }

    void CoreMain::out(string s){
        Figure *figure = dynamic_cast<Figure*>(getObject(gcf()));
        figure->out(s);
    }

    void CoreMain::box(double axes)
    {
        if (!isAxes(axes)) {
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }

        Axes* pAxes = dynamic_cast<Axes*>(getObject(axes));
        pAxes->setBorder(!pAxes->getBorder());
        pAxes->repaint();
    }

    void CoreMain::box(double axes, string state)
    {
        if (!isAxes(axes)) {
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }

        Axes* pAxes = dynamic_cast<Axes*>(getObject(axes));
        pAxes->setBorder(state == "on");
        pAxes->repaint();
    }

    void CoreMain::colorbar(double axes)
    {
        if (!isAxes(axes)) {
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }

        Axes* pAxes = dynamic_cast<Axes*>(getObject(axes));
        pAxes->setColorbarVisible(!pAxes->getColorbarVisible());
        pAxes->repaint();
    }

    void CoreMain::colorbar(double axes, std::string state, const std::vector<double>& range)
    {
        if (!isAxes(axes)) {
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }

        Axes* pAxes = dynamic_cast<Axes*>(getObject(axes));
        pAxes->setColorbarVisible(state == "on");
        if (range.size() == 2)
            pAxes->setColorbarRange(range);
        pAxes->repaint();
    }

    std::vector<double> CoreMain::colorbarRange(double axes)
    {
        if (!isAxes(axes)) {
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }

        Axes* pAxes = dynamic_cast<Axes*>(getObject(axes));
        return pAxes->getColorbarRange();
    }
}
