/**
* @file Object.cxx
* @date May 2018
* Copyright (C) 2018 Altair Engineering, Inc.  
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

#include "Object.h"
#include "LineStyle.h"
#include "OML_Error.h"
#include <cctype>
#include <algorithm>
#include <cfloat>
#include <iostream>
#include <cassert>
#include "GnuplotOutput.h"

using namespace std;

namespace omlplot{

    typedef Property::ValueType VALUETYPE;

    static bool isSameDouble(const double& a, const double& b){
        return ( abs(a - b) <= DBL_EPSILON );
    }

    HandlePool* Object::m_handlePool = new HandlePool;
    FigureHandlePool *Object::m_figureHandlePool = new FigureHandlePool;
    map<double, Object *> Object::m_objectMap;

    Object::Object()
    {
        m_ps.push_back(Property("parent", vector<double>(), PropertyType::VEC_DOUBLE) );
        m_ps.push_back(Property("children", vector<double>(), PropertyType::VEC_DOUBLE) );
    }

    Object::~Object(){
        double parent = getPropertyValue("parent").Scalar();
        Object *p = getObject(parent);
        p->removeChild(this);

        vector<double> children = getPropertyValue("children").Vector();
        Object *o = nullptr;
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it){
            double child = *it;
            o = getObject(child);
            delete o;
        }
        for (it = children.begin(); it != children.end(); ++it){
            double child = *it;
            m_objectMap.erase(child);
        }
        m_objectMap.erase(getHandle());
    }

    VALUETYPE Object::getPropertyValue(const string& name){
        Property &p = getProperty(name);
        return p.getValue();
    }

    vector<string> Object::getPropertyNames(){
        vector<string> names;
        vector<Property>::iterator it = m_ps.begin();
        for (; it != m_ps.end(); ++it){
            Property p = *it;
            names.push_back(p.getName());
        }
        return names;
    }

    void Object::setPropertyValue(const string& name, VALUETYPE value){
        Property &p = getProperty(name);
        if (name == "color" ||
            name == "xcolor" ||
            name == "ycolor" ||
            name == "zcolor" ||
            name == "facecolor"){
            if (value.isCurrency()){
                Currency c = value.getCurrency();
                Color clr(c);
                p.setValue(clr);
            } else if (value.isColor()){
                p.setValue(value);
            }
        } else {
            p.setValue(value);
        }
    }

    Property &Object::getProperty(const string& name){
        string n = name;
        int count = 0;
        // ignore the case
        std::transform(n.begin(), n.end(),
                       n.begin(), ::tolower);

        vector<Property>::iterator it, itret;
        // get all properties start with n
        for (it = m_ps.begin(); it != m_ps.end(); ++it){
            string pn = (*it).getName();
            if (pn.compare(0, n.length(), n) == 0 ){
                if (pn == n){
                    return *it;
                }
                ++count;
                itret = it;
            }            
        }

        // return the only match or throw error
        if (0 == count){
            throw OML_Error(OML_ERR_PLOT_AMBIGUOUS_PROPERTY);
        } else if (count > 1){
            throw OML_Error(OML_ERR_PLOT_AMBIGUOUS_PROPERTY);
        }
        return *itret;
    }

    double Object::getHandle(){
        VALUETYPE v = getPropertyValue("handle");
        return v.Scalar();
    }

    string Object::getType(){
        VALUETYPE v = getPropertyValue("type");
        return v.StringVal();
    }

    bool Object::isHandle(double h){
        return isSameDouble(h, 0) ||
            m_handlePool->has(h) || isFigure(h);
    }

    bool Object::isFigure(double h){
        return false;
    }

    bool Object::isAxes(double h){
        return false;
    }

    Object *Object::getObject(double h){
        map<double, Object *>::iterator it =
            m_objectMap.find(h);
        if (it != m_objectMap.end()){
            return it->second;
        }
        return nullptr;
    }

    void Object::setParent(Object *o){
        vector<double> parent;
        parent.push_back(o->getHandle());
        setPropertyValue("parent", parent);
    }

    void Object::addChild(Object *o){
        vector<double> children = getPropertyValue("children").Vector();
        children.push_back(o->getHandle());
        setPropertyValue("children", children);

        vector<double> parent;
        parent.push_back(this->getHandle());
        o->setPropertyValue("parent", parent);
    }

    void Object::removeChild(Object *o){
        vector<double> children = getPropertyValue("children").Vector();
        vector<double>::iterator it = children.begin();
        double ch = o->getHandle();
        for (; it != children.end(); ++it){
            if (*it == ch){
                break;
            }
        }
        if (it != children.end()){
            children.erase(it);
        }
        setPropertyValue("children", children);
    }

    void Object::repaint(){
        vector<double> parent = getPropertyValue("parent").Vector();
        Object *o = getObject(parent[0]);
        o->repaint();
    }

    void Object::dump(Object *o, int level){
        for (int i = 0; i < level; i++){
            cout << '\t';
        }
        cout << o->getType() << " ( " << o->getHandle() << " )" << endl;

        vector<double> children = o->getPropertyValue("children").Vector();
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it){
            double child = *it;
            dump(getObject(child), level+1);
        }
        if (level == 0){
            cout << "in the map:" << endl;
            auto it = m_objectMap.begin();
            for (; it != m_objectMap.end(); ++it){
                auto &o = *it;
                cout << o.second->getType() << " ( " << o.first << " ) " << endl;
            }
        }
    }

    Root::Root(){
        m_ps.push_back(Property("type", string("root"), PropertyType::STRING) );
        m_ps.push_back(Property("handle", double(0), PropertyType::DOUBLE) );
        m_ps.push_back(Property("currentfigure", double(0), PropertyType::DOUBLE) );
        m_objectMap[double(0)] = this;
    }

    bool Root::isFigure(double h){
        vector<double> children = getPropertyValue("children").Vector();
        Object *o = nullptr;
        vector<double>::iterator it = children.begin();
        for ( ; it != children.end(); ++it){
            double child = *it;
            o = getObject(child);
            if (o->isFigure(h) ){
                return true;
            }
        }
        return false;
    }

    bool Root::isAxes(double h){
        vector<double> children = getPropertyValue("children").Vector();
        Object *o = nullptr;
        vector<double>::iterator it = children.begin();
        for ( ; it != children.end(); ++it){
            double child = *it;
            o = getObject(child);
            if (o->isAxes(h) ){
                return true;
            }
        }
        return false;
    }

    void Root::repaint(){
        update(nullptr);
    }

    void Root::update(GnuplotOutput *out){
        vector<double> children = getPropertyValue("children").Vector();
        Figure *o = nullptr;
        vector<double>::iterator it = children.begin();
        for ( ; it != children.end(); ++it){
            double child = *it;
            o = dynamic_cast<Figure *>(getObject(child));
            o->update();
        }
    }

    Figure::Figure()
        :m_out(new GnuplotOutput)
    {
        double handle = m_figureHandlePool->allocHandle();
        init(handle);
    }

    Figure::Figure(int fh)
        :m_out(new GnuplotOutput)
    {
        if (! m_figureHandlePool->isFreeHandle(fh) ){
            assert(0);          // fh should always be free;
        }
        m_figureHandlePool->setHandleUsed(fh);
        init(fh);
    }

    void Figure::init(double handle){
        m_ps.push_back(Property("type", string("figure"), PropertyType::STRING) );
        m_ps.push_back(Property("handle", handle,PropertyType::DOUBLE) );
        m_ps.push_back(Property("currentaxes", double(0), PropertyType::DOUBLE) );
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("units", string("pixels"), PropertyType::STRING));
        vector<double> pos;
        pos.push_back(0);pos.push_back(0);
        pos.push_back(800);pos.push_back(600);
        m_ps.push_back(Property("position", pos, PropertyType::VEC_DOUBLE));

        m_objectMap[handle] = this;
        m_out->setWindowId(handle);
        m_out->printf("clear\n");
        m_out->flush();
    }

    Figure::~Figure(){
        int fh = int(getHandle());
        m_figureHandlePool->releaseHandle(fh);
    }

    bool Figure::isFigure(double h){
        const double thisH = getHandle();
        if (isSameDouble(thisH, h)) {
            return true;
        }
        return false;
    }

    bool Figure::isAxes(double h){
        vector<double> children = getPropertyValue("children").Vector();
        Object *o = nullptr;
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it){
            double child = *it;
            o = getObject(child);
            if (o->isAxes(h) ){
                return true;
            }
        }
        return false;
    }

    void Figure::repaint(){
        update(nullptr);
    }

    void Figure::update(GnuplotOutput *out){
        m_out->printf("clear\n");
        m_out->printf("set multiplot\n");
        m_out->flush();

        vector<double> children = getPropertyValue("children").Vector();
        Object *o = nullptr;
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it){
            double child = *it;
            o = getObject(child);
            o->update(m_out.get());
        }

        m_out->printf("unset multiplot\n\n\n\n");
        m_out->flush();
    }

    void Figure::clear(){
        vector<double> children = getPropertyValue("children").Vector();
        Object *o = nullptr;
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it){
            double child = *it;
            o = getObject(child);
            delete o;
        }
        
        for (it = children.begin(); it != children.end(); ++it){
            double child = *it;
            m_objectMap.erase(child);
        }
        children.clear();
        setPropertyValue("children", children);
        setPropertyValue("currentaxes", double(0));
        
        repaint();
    }

    void Figure::saveas(string _filename, string fmt){
        string filename = _filename;
        if (fmt == ""){
            fmt = string("png");
        } else {
            filename += fmt;
        }

        m_out->printf("set terminal push\n");
        m_out->printf("set terminal png\n");
        m_out->printf("set output \"%s\"\n", filename.c_str());
        m_out->printf("replot\n");
        m_out->printf("set terminal pop\n");
        m_out->printf("set output\n");
        m_out->flush();
    }

    void Figure::out(string s){
        m_out->printf("%s\n", s.c_str());
        m_out->flush();
    }

    Axes::Axes()
        :m_title(new Text),
         m_xlabel(new Text),
         m_ylabel(new Text),
         m_zlabel(new Text),
         _axisNeedRepaint(false)
    {
        double handle = m_handlePool->allocHandle();
        // all of the properties of axes are added here
        m_ps.push_back(Property("type", string("axes"), PropertyType::STRING) );
        m_ps.push_back(Property("handle", handle, PropertyType::DOUBLE) );
        m_ps.push_back(Property("title", Currency(m_title.get(), "text"), PropertyType::POINTER) );
        vector<double> lim;
        lim.push_back(0); lim.push_back(1);
        m_ps.push_back(Property("xlim", lim, PropertyType::VEC_DOUBLE) );
        m_ps.push_back(Property("ylim", lim, PropertyType::VEC_DOUBLE) );
        m_ps.push_back(Property("zlim", lim, PropertyType::VEC_DOUBLE) );
        m_ps.push_back(Property("xgrid", string("off"), PropertyType::STRING) );
        m_ps.push_back(Property("ygrid", string("off"), PropertyType::STRING) );
        m_ps.push_back(Property("zgrid", string("off"), PropertyType::STRING) );
        m_ps.push_back(Property("xlabel", Currency(m_xlabel.get(), "text"), PropertyType::POINTER) );
        m_ps.push_back(Property("ylabel", Currency(m_ylabel.get(), "text"), PropertyType::POINTER) );
        m_ps.push_back(Property("zlabel", Currency(m_zlabel.get(), "text"), PropertyType::POINTER) );
        m_ps.push_back(Property("xminorgrid", string("off"), PropertyType::STRING) );
        m_ps.push_back(Property("yminorgrid", string("off"), PropertyType::STRING) );
        m_ps.push_back(Property("zminorgrid", string("off"), PropertyType::STRING) );
        m_ps.push_back(Property("xscale", string("linear"), PropertyType::STRING) );
        m_ps.push_back(Property("yscale", string("linear"), PropertyType::STRING) );
        m_ps.push_back(Property("zscale", string("linear"), PropertyType::STRING) );
        m_ps.push_back(Property("xcolor", Color(string("black")), PropertyType::COLOR));
        m_ps.push_back(Property("ycolor", Color(string("black")), PropertyType::COLOR));
        m_ps.push_back(Property("zcolor", Color(string("black")), PropertyType::COLOR));
        m_ps.push_back(Property("nextplot", string("replace"), PropertyType::STRING) );
        m_ps.push_back(Property("units", string("normalized"), PropertyType::STRING) );
        vector<double> pos;
        pos.push_back(0);pos.push_back(0);
        pos.push_back(1);pos.push_back(1);
        m_ps.push_back(Property("position", pos, PropertyType::VEC_DOUBLE) );
        m_ps.push_back(Property("legend", string("off"), PropertyType::STRING) );
        m_ps.push_back(Property("color", Color(string("white")), PropertyType::COLOR));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));

        m_objectMap[handle] = this;
        m_title->setParent(this);
    }

    void Axes::clear(){
        vector<double> children = getPropertyValue("children").Vector();
        Object *o = nullptr;
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it){
            double child = *it;
            o = getObject(child);
            delete o;
        }

        for (it = children.begin(); it != children.end(); ++it){
            double child = *it;
            m_objectMap.erase(child);
        }

        children.clear();
        setPropertyValue("children", children);

        m_title->setPropertyValue("string", string(""));

        repaint();
    }

    void Axes::repaint(){
        vector<double> parent = getPropertyValue("parent").Vector();
        Figure *pFig = dynamic_cast<Figure *>(getObject(parent[0]));
        pFig->repaint();
    }        

    void Axes::update(GnuplotOutput *out){
        string v = this->getPropertyValue("visible").StringVal();
        if (v == "off"){
            return;
        }

        if (axisNeedRepaint()) {
            vector<double> xrange = getPropertyValue("xlim").Vector();
            out->printf("set xrange [%g:%g]\n", xrange[0], xrange[1]);
            vector<double> yrange = getPropertyValue("ylim").Vector();
            out->printf("set yrange [%g:%g]\n", yrange[0], yrange[1]);
            vector<double> zrange = getPropertyValue("zlim").Vector();
            out->printf("set zrange [%g:%g]\n", zrange[0], zrange[1]);
            setAxisNeedRepaint(false);
        }

        string bgcolor = getPropertyValue("color").ColorString();
        out->printf("set obj 1 rectangle from graph 0,0 to graph 1,1 behind fc rgb '%s'\n", bgcolor.c_str());        
        string xcolor = getPropertyValue("xcolor").ColorString();
        out->printf("set xtics textcolor rgb '%s'\n", xcolor.c_str());
        string ycolor = getPropertyValue("ycolor").ColorString();
        out->printf("set ytics textcolor rgb '%s'\n", ycolor.c_str());
        string zcolor = getPropertyValue("zcolor").ColorString();
        out->printf("set ztics textcolor rgb '%s'\n", zcolor.c_str());

        string xscale = getPropertyValue("xscale").ColorString();
        if (xscale == "log"){
            out->printf("set logscale x\n");
        }
        string yscale = getPropertyValue("yscale").ColorString();
        if (yscale == "log"){
            out->printf("set logscale y\n");
        }
        string zscale = getPropertyValue("zscale").ColorString();
        if (zscale == "log"){
            out->printf("set logscale z\n");
        }

        m_title->update(out, "title");
        m_xlabel->update(out, "xlabel");
        m_ylabel->update(out, "ylabel");
        m_zlabel->update(out, "zlabel");

        vector<double> pos = getPropertyValue("position").Vector();
        out->printf("set origin %f,%f\n", pos[0], pos[1]);
        out->printf("set size %f,%f\n", pos[2], pos[3]);

        {
            out->printf("unset key\nplot '-' with points ps 0\n0 0\n1 1\ne\n");
        }
        string xgrid = getPropertyValue("xgrid").StringVal();
        string xmgrid = getPropertyValue("xminorgrid").StringVal();
        if (xmgrid == "on"){
            out->printf("set mxtics 5\nset grid xtics mxtics\nrefresh\n");
        } else {
            if (xgrid == "on"){
                out->printf("set mxtics 1\nset grid xtics mxtics\nrefresh\n");
            } else {
                out->printf("set mxtics 1\nset grid noxtics\nrefresh\n");
            }
        }
        string ygrid = getPropertyValue("ygrid").StringVal();
        string ymgrid = getPropertyValue("yminorgrid").StringVal();
        if (ymgrid == "on"){
            out->printf("set mytics 5\nset grid ytics mytics\nrefresh\n");
        } else {
            if (ygrid == "on"){
                out->printf("set grid ytics\nrefresh\nrefresh\n");
            } else {
                out->printf("set mytics 1\nset grid noytics\nrefresh\n");
            }
        }
        string zgrid = getPropertyValue("zgrid").StringVal();
        if (zgrid == "on"){
            out->printf("set grid ztics\nrefresh\nrefresh\n");
        } else {
            out->printf("set grid noztics\nrefresh\n");
        }
        string lState = getPropertyValue("legend").StringVal();
        if (lState == "on"){
            out->printf("set key on\n");
        } else {
            out->printf("set key off\n");
        }

        vector<double> children = getPropertyValue("children").Vector();
        Drawable *pChild = nullptr;
        vector<string> linestyle;
        vector<string> usingclause;
        vector<string> legend;
        vector<string> withclause;
        vector<vector<double> > zdata;

        int lineId = 1;
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it){
            double child = *it;
            pChild = dynamic_cast<Drawable *>(getObject(child));
            if (pChild){
                linestyle.push_back(pChild->getLineStyle(lineId));
                usingclause.push_back(pChild->getUsingClause());
                legend.push_back(pChild->getLegend());
                withclause.push_back(pChild->getWithClause(lineId));
                zdata.push_back(pChild->getZdata());
            } else {
                Text *pChild = dynamic_cast<Text *>(getObject(child));
                if (pChild){
                    pChild->update(out);
                }
            }

            ++lineId;
        }

        bool is3D = false;
        if ( (zdata.size() > 0) &&
             (zdata[0].size() > 0) ){
            is3D = true;
        }

        int lineCount = usingclause.size();
        for (lineId = 0; lineId < lineCount; ++lineId){
            string style = linestyle[lineId];
            out->printf(style);
            out->flush();
        }

        if (lineCount > 0){
            out->printf("clear\n");
            stringstream ss;
            if (is3D){
                ss << "splot ";
            } else {
                ss << "plot ";
            }
            for (int i = 0; i < lineCount; i++){
                ss << usingclause[i] << " ";
                ss << legend[i] << " ";
                ss << withclause[i] << " ";
                if (i < (lineCount - 1) ){
                    ss << ", ";
                }
            }
            ss << "\n";
            string cmd = ss.str();
            out->printf(cmd.c_str());

            vector<double>::iterator it = children.begin();
            for (; it != children.end(); ++it){
                double child = *it;
                pChild = dynamic_cast<Drawable *>(getObject(child));
				if (pChild){
					pChild->putData(out);
                    out->printf("e\n");
				}
            }
        }
        it = children.begin();
        for (; it != children.end(); ++it){
            double child = *it;
            pChild = dynamic_cast<Drawable *>(getObject(child));
			if (pChild){
				pChild->cleanup(out);
			}
        }

        out->printf("set title ''\n");
        out->printf("set xlabel ''\n");
        out->printf("set ylabel ''\n");
        out->printf("set zlabel ''\n");
        out->printf("unset grid\n");
        out->printf("unset logscale xyz\n");
        out->printf("set key off\n");
        out->printf("\n\n");
    }

    bool Axes::isAxes(double h){
        const double thisH = getHandle();
        if (isSameDouble(thisH, h)){
            return true;
        }
        return false;
    }

    void Axes::hold(bool setHold){
        if (setHold){
            setPropertyValue("nextplot", "add");
        } else {
            setPropertyValue("nextplot", "replace");
        }
    }

    bool Axes::ishold(){
        string hold = getPropertyValue("nextplot").StringVal();
        return string("add") == hold;
    }

    void Axes::setAxisNeedRepaint(bool b){
        _axisNeedRepaint = b;
    }

    bool Axes::axisNeedRepaint(){
        return _axisNeedRepaint;
    }

    Drawable::Drawable(){
        double handle = m_handlePool->allocHandle();
        m_ps.push_back(Property("handle", handle , PropertyType::DOUBLE));
        vector<double> d;
        m_ps.push_back(Property("xdata", d, PropertyType::VEC_DOUBLE) );
        m_ps.push_back(Property("ydata", d, PropertyType::VEC_DOUBLE) );
        m_ps.push_back(Property("zdata", d, PropertyType::VEC_DOUBLE) );
        m_ps.push_back(Property("displayname", string(""), PropertyType::STRING) );

        m_objectMap[handle] = this;
    }

    void Drawable::init(const LineData &ld){
        setPropertyValue("xdata", ld.x);
        setPropertyValue("ydata", ld.y);
        setPropertyValue("zdata", ld.z);

        for (int i = 0; i < ld.properties.size(); i++){
            setPropertyValue(ld.properties[i], ld.values[i]);
        }

        _xcolcount = ld.xcolcount;
        _ycolcount = ld.ycolcount;
        _zcolcount = ld.zcolcount;
    }

    void Drawable::update(GnuplotOutput *){
    }

    vector<double> Drawable::getXdata(){
        return getPropertyValue("xdata").Vector();
    }

    vector<double> Drawable::getYdata(){
        return getPropertyValue("ydata").Vector();
    }

    vector<double> Drawable::getZdata(){
        return getPropertyValue("zdata").Vector();
    }

    void Drawable::putData(GnuplotOutput *out){
        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> z = getZdata();

        if (z.size() == 0){     // 2d data
            int count = min(x.size(), y.size());
            for (int j = 0; j < count; j++){
                out->printf("%g %g\n", x[j], y[j]);
            }
        } else {                // 3d data
            int count = min(x.size(), y.size());
            count = min(count, (int)z.size());
            for (int j = 0; j < count; j++){
                out->printf("%g %g %g\n", x[j], y[j], z[j]);
            }
        }
    }

    Line::Line(){
        m_ps.push_back(Property("type", string("line"), PropertyType::STRING) );
        m_ps.push_back(Property("linestyle", string("-"), PropertyType::STRING));
        m_ps.push_back(Property("color", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("linewidth", double(1), PropertyType::DOUBLE));
        m_ps.push_back(Property("marker", string(""), PropertyType::STRING));
        m_ps.push_back(Property("markerfacecolor", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("markersize", double(1), PropertyType::DOUBLE));
        m_ps.push_back(Property("basevalue", double(0), PropertyType::DOUBLE));
    }

    void Line::init(const LineData &ld){
        Drawable::init(ld);

        // update line-specific properties
        LineStyle ls = LineStyle(ld);
        setPropertyValue("linestyle", ls.m_lineStyle);
        setPropertyValue("color", ls.m_lineColor);
        setPropertyValue("linewidth", ls.m_lineWidth);
        setPropertyValue("marker", ls.m_markerStyle);
        setPropertyValue("markerfacecolor", ls.m_markerColor);
        setPropertyValue("markersize", ls.m_markerSize);
        setPropertyValue("displayname", ls.m_legend);

        int propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++){
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string Line::getUsingClause(){
        return "'-' using 1:2";
    }

    string Line::getLegend(){
        string legend = getPropertyValue("displayname").StringVal();
        return string("title '") + legend + string("'");
    }

    string Line::getWithClause(int lineId){
        string linestyle = getPropertyValue("linestyle").StringVal();
        stringstream ss;
        if (linestyle == ""){
            ss << "with points linestyle " << lineId;
        } else {
            ss << "with lp linestyle " << lineId;
        }
        return ss.str();
    }

    string Line::getLineStyle(int line){
        string linestyle = getPropertyValue("linestyle").StringVal();
        int dt = 1;
        if (linestyle == "-"){
            dt = 1;
        } else if (linestyle == "--") {
            dt = 2;
        } else if (linestyle == ":"){
            dt = 3;
        } else if (linestyle == "-."){
            dt = 4;
        } else if (linestyle == ""){
            dt = 0;
        }
        string lineColor = getPropertyValue("color").ColorString();
        int lineWidth = int(getPropertyValue("linewidth").Scalar());
        string marker = getPropertyValue("marker").StringVal();
        int pt = 0;
        if (marker == ""){
            pt = 0;
        } else if (marker == "s"){
            pt = 5;
        } else if (marker == "^"){
            pt = 9;
        } else if (marker == "v"){
            pt = 11;
        } else if (marker == "x"){
            pt = 2;
        } else if (marker == "o"){
            pt = 7;
        } else if (marker == "d"){
            pt = 13;
        } else if (marker == "+"){
            pt = 1;
        } else if (marker == "*"){
            pt = 3;
        } else if (marker == "."){
            pt = 0;
        } else {
            pt = 6;
        }
            
        string markerColor = getPropertyValue("markerfacecolor").ColorString();
        int markerSize = int(getPropertyValue("markersize").Scalar());

        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line << " linetype " << dt
           << " lc rgb '" << lineColor << "' lw " << lineWidth
           << " pt " << pt << " ps " << markerSize
            //<< " tc rgb '" << markerColor << "'"
           << "\n";
        return ss.str();
    }

    void Line::cleanup(GnuplotOutput *out){
    }

    string Polar::getLineStyle(int line) {
        string style = Line::getLineStyle(line);
        stringstream ss;
        ss << style << "set polar\n"
           << "set grid polar\n"
           << "unset border\n"
           << "unset xtics\n"
           << "unset ytics\n";
        return ss.str();
    }

    void Polar::cleanup(GnuplotOutput *out){
        out->printf("unset polar\n");
        out->printf("set border\n");
        out->printf("set xtics\n");
        out->printf("set ytics\n");
    }

    Fill::Fill(){
        m_ps.push_back(Property("color", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("type", string("patch"), PropertyType::STRING) );
    }

    void Fill::init(const LineData &ld){
        Drawable::init(ld);

        LineStyle ls = LineStyle(ld);
        setPropertyValue("color", ls.m_lineColor);

        int propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++){
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string Fill::getUsingClause(){
        return "'-' using 1:2";
    }

    string Fill::getLegend(){
        string legend = getPropertyValue("displayname").StringVal();
        return string("title '") + legend + string("'");
    }

    string Fill::getWithClause(int lineId){
        stringstream ss;
        ss << "with filledcurves linestyle " << lineId;
        return ss.str();
    }

    string Fill::getLineStyle(int line){
        int dt = 1;
        string lineColor = getPropertyValue("color").ColorString();

        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line << " linetype " << dt
           << " lc rgb '" << lineColor
           << "\n";
        return ss.str();
    }

    void Fill::cleanup(GnuplotOutput *out){
    }

    Area::Area(){
        m_ps.push_back(Property("type", string("hggroup"), PropertyType::STRING) );
        m_ps.push_back(Property("basevalue", double(0), PropertyType::DOUBLE));
        m_ps.push_back(Property("facecolor", Color(string("blue")), PropertyType::COLOR));
    }

    string Area::getUsingClause(){
        return "'-' using 1:2";
    }

    string Area::getLegend(){
        string legend = getPropertyValue("displayname").StringVal();
        return string("title '") + legend + string("'");
    }

    string Area::getWithClause(int lineId){
        stringstream ss;
        double basevalue = getPropertyValue("basevalue").Scalar();
        ss << "with filledcurve y1="
           << basevalue << " linestyle " << lineId;
        return ss.str();
    }

    string Area::getLineStyle(int line){
        string lineColor = getPropertyValue("facecolor").ColorString();

        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line
           << " lc rgb '" << lineColor
           << "\n";
        return ss.str();
    }

    void Area::cleanup(GnuplotOutput *out){
    }

    HggroupBar::HggroupBar(){
        m_ps.push_back(Property("type", string("hggroup"), PropertyType::STRING) );
        m_ps.push_back(Property("barwidth", double(0.8), PropertyType::DOUBLE));
        m_ps.push_back(Property("barlayout", string("grouped"), PropertyType::STRING));
        m_ps.push_back(Property("facecolor", Color(string("blue")), PropertyType::COLOR));
    }

    void HggroupBar::init(const LineData &ld){
        Drawable::init(ld);

        // update line-specific properties
        LineStyle ls = LineStyle(ld);
        setPropertyValue("facecolor", ls.m_lineColor);
        setPropertyValue("displayname", ls.m_legend);

        int propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++){
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string HggroupBar::getUsingClause(){
        string color = getPropertyValue("facecolor").ColorString();        
        stringstream ss;
        ss << "'-' using 2:xtic(1) lc rgb '" << color << "' ";
        return ss.str();
    }

    string HggroupBar::getLegend(){
        string legend = getPropertyValue("displayname").StringVal();
        return string("title '") + legend + string("'");
    }

    string HggroupBar::getWithClause(int lineId){
        return "";
    }

    string HggroupBar::getLineStyle(int line){
        double barwidth = getPropertyValue("barwidth").Scalar();
        string barlayout =  getPropertyValue("barlayout").StringVal();
        string color = getPropertyValue("facecolor").ColorString();        
        string layout = "cluster";
        if (barlayout == "stacked"){
            layout = "rowstacked";
        }
        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line
           << " lc rgb '" << color << "'\n";
        //ss << "set boxwidth " << barwidth/3 << " relative\n";
        ss << "set yrange [*<0:]\n";
        ss << "set style fill solid\n"
           << "set style data histogram\n"
           << "set style histogram " << layout <<"\n"
           << "set boxwidth 0.8\n";
        return ss.str();
    }

    void HggroupBar::cleanup(GnuplotOutput *out) {
    }

    Hist::Hist(){
        setPropertyValue("barwidth", 0.4);
    }

    void Hist::init(const LineData &ld){
        HggroupBar::init(ld);
    }

    string Hist::getLineStyle(int line){
        stringstream ss;
        ss << HggroupBar::getLineStyle(line);
        ss << "set style data histogram\n";
        ss << "set style histogram cluster gap 1\n";
        return ss.str();
    }

    HggroupScatter::HggroupScatter(){
        m_ps.push_back(Property("type", string("hggroup"), PropertyType::STRING) );
        m_ps.push_back(Property("marker", string("o"), PropertyType::STRING));
        m_ps.push_back(Property("markerfacecolor", string("blue"), PropertyType::COLOR));
        m_ps.push_back(Property("markersize", double(1), PropertyType::DOUBLE));
    }

    string HggroupScatter::getUsingClause(){
        return "'-' using 1:2";
    }

    string HggroupScatter::getLegend(){
        return "title ''";
    }

    string HggroupScatter::getWithClause(int lineId){
        stringstream ss;
        ss << "with points linestyle " << lineId;
        return ss.str();
    }

    string HggroupScatter::getLineStyle(int line){
        string marker = getPropertyValue("marker").StringVal();
        int pt = 0;
        if (marker == ""){
            pt = 0;
        } else if (marker == "s"){
            pt = 5;
        } else if (marker == "^"){
            pt = 9;
        } else if (marker == "v"){
            pt = 11;
        } else if (marker == "x"){
            pt = 2;
        } else if (marker == "o"){
            pt = 7;
        } else if (marker == "d"){
            pt = 13;
        } else if (marker == "+"){
            pt = 1;
        } else if (marker == "*"){
            pt = 3;
        } else if (marker == "."){
            pt = 0;
        } else {
            pt = 6;
        }
            
        string markerColor = getPropertyValue("markerfacecolor").ColorString();
        int markerSize = int(getPropertyValue("markersize").Scalar());

        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line << " dt " << 0
           << " lc rgb '" << markerColor << "'"
           << " pt " << pt << " ps " << markerSize
            //<< " tc rgb '" << markerColor << "'"
           << "\n";
        return ss.str();
    }

    void HggroupScatter::cleanup(GnuplotOutput *out) {
    }

    string HggroupScatter3::getUsingClause(){
        return "'-' using 1:2:3";
    }

    Surface::Surface(){
        m_ps.push_back(Property("type", string("surface"), PropertyType::STRING) );
    }

    string Surface::getUsingClause(){
        return "'-' using 1:2:3 with lines";
    }

    string Surface::getLegend(){
        return "title ''";
    }

    string Surface::getWithClause(int lineId){
        return " ";
    }

    string Surface::getLineStyle(int line){
        return "set palette defined ( 0 \"blue\", 3 \"green\", 6 \"yellow\", 10 \"red\")\n"
            "set pm3d\n"
            "unset colorbox\n"
            "set view 60, 315,1,1\n";
    }

    void Surface::cleanup(GnuplotOutput *out) {
    }

    void Surface::putData(GnuplotOutput *out){
        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> z = getZdata();

        if (_xcolcount > 1 && _ycolcount > 1){
            int xsize = x.size();
            int ysize = y.size();
            int zsize = z.size();
            if (! ((xsize == ysize) &&
                   (xsize == zsize)) ){
                throw;
            }

            int rowcount = zsize / _zcolcount;
            for (int i = 0; i < zsize; i++){
                out->printf("%g %g %g\n",
                            x[i], y[i], z[i]);
                if ((i + 1) % rowcount == 0){
                    out->printf("\n");
                }
            }
        } else {
            int col = y.size();
            for (int xsub = 0; xsub < x.size(); xsub++){
                for (int ysub = 0; ysub < col; ysub++){
                    out->printf("%g %g %g\n",
                                x[xsub], y[ysub], z[xsub * col + ysub]);
                }
                out->printf("\n");
            }
        }
    }

    string Mesh::getLineStyle(int line){
        return "set pal maxcolor 0\n"
            "unset pm3d\n"
            "set surface\n"
            "set isosamples 40, 40; set samples 40, 40\n"
            "set view 60, 315,1,1\n";
    }

    string Line3::getLineStyle(int line){
        return "set surface\n";
    }

    string Line3::getUsingClause(){
        return "'-' using 1:2:3";
    }

    string Line3::getWithClause(int lineId){
        stringstream ss;
        ss << "with line linestyle " << lineId;
        return ss.str();
    }

    string Contour3::getLineStyle(int line){
        return "set contour surface\n"
            "unset surface\n"
            "set cntrparam levels 10\n";
            "set view 60, 315,1,1\n";
            //"set dgrid3d\n";            
    }

    void Contour3::cleanup(GnuplotOutput *out) {
        out->printf("unset contour\n");
        out->printf("set surface\n");
        //out->printf("unset dgrid3d\n");
    }

    string Contour3::getUsingClause(){
        return "'-' using 1:2:3 with lines nosurf";
    }

    string Contour::getLineStyle(int line){
        return "set contour base\n"
            "unset surface\n"
            "set view map\n"
            //"set dgrid3d\n"
            "set cntrparam levels 10\n";
    }

    void Contour::cleanup(GnuplotOutput *out) {
        out->printf("unset contour\n");
        out->printf("unset view\n");
        //out->printf("unset dgrid3d\n");
    }

    string Stem::getWithClause(int lineId){
        stringstream ss;
        ss << "with impulses linestyle " << lineId;
        return ss.str();
    }

    string Loglog::getLineStyle(int line) {
        string style = Line::getLineStyle(line);
        stringstream ss;
        ss << style << "set logscale xy\n";
        return ss.str();
    }

    void Loglog::cleanup(GnuplotOutput *out){
        out->printf("unset logscale xy\n");
    }

    string Semilogx::getLineStyle(int line) {
        string style = Line::getLineStyle(line);
        stringstream ss;
        ss << style << "set logscale x\n";
        return ss.str();
    }

    void Semilogx::cleanup(GnuplotOutput *out){
        out->printf("unset logscale x\n");
    }

    string Semilogy::getLineStyle(int line) {
        string style = Line::getLineStyle(line);
        stringstream ss;
        ss << style << "set logscale y\n";
        return ss.str();
    }

    void Semilogy::cleanup(GnuplotOutput *out){
        out->printf("unset logscale y\n");
    }

    Text::Text(){
        double handle = m_handlePool->allocHandle();
        m_ps.push_back(Property("type", string("text"), PropertyType::STRING) );
        m_ps.push_back(Property("handle", handle, PropertyType::DOUBLE));
        m_ps.push_back(Property("string", string(""), PropertyType::STRING) );
        m_ps.push_back(Property("color", Color(string("black")), PropertyType::COLOR));
        m_ps.push_back(Property("fontname", string("arial"), PropertyType::STRING));
        m_ps.push_back(Property("fontsize", double(12), PropertyType::DOUBLE));
        m_ps.push_back(Property("fontweight", string("normal"), PropertyType::STRING));
        m_ps.push_back(Property("fontangle", string("regular"), PropertyType::STRING));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("x", double(0), PropertyType::DOUBLE));
        m_ps.push_back(Property("y", double(0), PropertyType::DOUBLE));
        m_ps.push_back(Property("z", double(0), PropertyType::DOUBLE));

        m_objectMap[handle] = this;
    }

    void Text::update(GnuplotOutput *out){
        double x = getPropertyValue("x").Scalar();
        double y = getPropertyValue("y").Scalar();
        double z = getPropertyValue("z").Scalar();
        string text = getPropertyValue("string").StringVal();

        string color = getPropertyValue("color").ColorString();
        string fontname = getPropertyValue("fontname").StringVal();
        double fontsize = getPropertyValue("fontsize").Scalar();
        string fontweight = getPropertyValue("fontweight").StringVal();
        string fontangle = getPropertyValue("fontangle").StringVal();

        stringstream ss;
        if (fontangle == "italic"){
            ss << "{/:Italic ";
        }
        if (fontweight == "bold"){
            ss << "{/:Bold ";
        }
        ss << text;
        if (fontweight == "bold"){
            ss << "}";
        }
        if (fontangle == "italic"){
            ss << "}";
        }

        out->printf("set label \"%s\" at %g,%g,%g font \"%s,%g\" tc rgb \"%s\" enhanced\n",
                    ss.str().c_str(), x, y, z,
                    fontname.c_str(), fontsize,
                    color.c_str());
    }

    void Text::update(GnuplotOutput *out, const string type){
        string text = getPropertyValue("string").StringVal();
        string color = getPropertyValue("color").ColorString();
        string fontname = getPropertyValue("fontname").StringVal();
        double fontsize = getPropertyValue("fontsize").Scalar();
        string fontweight = getPropertyValue("fontweight").StringVal();
        string fontangle = getPropertyValue("fontangle").StringVal();

        stringstream ss;
        if (fontangle == "italic"){
            ss << "{/:Italic ";
        }
        if (fontweight == "bold"){
            ss << "{/:Bold ";
        }
        ss << text;
        if (fontweight == "bold"){
            ss << "}";
        }
        if (fontangle == "italic"){
            ss << "}";
        }

        out->printf("set %s \"%s\" font \"%s,%g\" tc rgb \"%s\" enhanced\n",
                    type.c_str(), ss.str().c_str(),
                    fontname.c_str(), fontsize,
                    color.c_str());
    }

    HandlePool::HandlePool()
        :m_gen(m_rd()),
         m_dis(std::uniform_real_distribution<>(10, 20))
    {
    }

    double HandlePool::allocHandle(){
        double rnd = m_dis(m_gen);
        while (m_pool.find(rnd) != m_pool.end() ){
            rnd = m_dis(m_gen);
        }
        m_pool.insert(rnd);
        return rnd;
    }

    void HandlePool::releaseHandle(double h){
        set<double>::iterator it = m_pool.find(h);
        m_pool.erase(it);
    }

    bool HandlePool::has(double h){
        return m_pool.count(h) > 0;
    }

    FigureHandlePool::FigureHandlePool()
        :m_nextHandle(1){
    }

    int FigureHandlePool::allocHandle(){
        int h = m_nextHandle;
        setHandleUsed(h);
        return h;
    }

    bool FigureHandlePool::isFreeHandle(int h){
        return (m_pool.count(h) == 0);
    }

    void FigureHandlePool::setHandleUsed(int h){
        if (! isFreeHandle(h) ){
            assert(0);          // should never happen
        }
        m_pool.insert(h);
        if (m_nextHandle == h){
            ++m_nextHandle;
        }
        while (! isFreeHandle(m_nextHandle) ){
            ++m_nextHandle;
        }
    }

    void FigureHandlePool::releaseHandle(int h){
        m_pool.erase(h);
        if (m_nextHandle > h){
            m_nextHandle = h;
        }
    }

    CurrencyAndColor::CurrencyAndColor(int c){
        m_currency = Currency(double(c));
        m_type = CURRENCY;
    };
        
    CurrencyAndColor::CurrencyAndColor(double c){
        m_currency = Currency(c);
        m_type = CURRENCY;
    };
        
    CurrencyAndColor::CurrencyAndColor(const char* c){
        m_currency = Currency(string(c));
        m_type = CURRENCY;
    };

    CurrencyAndColor::CurrencyAndColor(std::string c){
        m_currency = Currency(c);
        m_type = CURRENCY;
    };

    CurrencyAndColor::CurrencyAndColor(std::vector<double> c){
        m_currency = Currency(c);
        m_type = CURRENCY;
    };

    CurrencyAndColor::CurrencyAndColor(Currency c){
        m_currency = c;
        m_type = CURRENCY;
    };
    CurrencyAndColor::CurrencyAndColor(Color c){
        m_color = c;
        m_type = COLOR;
    }

    double CurrencyAndColor::Scalar(){
        return m_currency.Scalar();
    }

    std::string CurrencyAndColor::StringVal(){
        return m_currency.StringVal();
    }

    std::vector<double> CurrencyAndColor::Vector(){
        return m_currency.Vector();
    }

    void *CurrencyAndColor::BoundObject() {
        return m_currency.BoundObject();
    }

    std::string CurrencyAndColor::ColorString(){
        return m_color.getString();
    }

    Currency CurrencyAndColor::getCurrency(){
        return m_currency;
    }

    vector<double> CurrencyAndColor::getColor(){
        return m_color.getComponent();
    }

    bool CurrencyAndColor::isCurrency(){
        return m_type == CURRENCY;
    }

    bool CurrencyAndColor::isColor(){
        return m_type == COLOR;
    }

    
}
