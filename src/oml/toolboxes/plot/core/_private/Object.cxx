/**
* @file Object.cxx
* @date May 2018
* Copyright (C) 2018-2022 Altair Engineering, Inc.  
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
        :m_type(ObjectType::NO_TYPE)
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

    bool Object::setPropertyValue(const string& name, VALUETYPE value){
        Property &p = getProperty(name);
		if (p.getType() == UNSUPPORTED)
			return false;
        if (name == "color" ||
            name == "xcolor" ||
            name == "ycolor" ||
            name == "zcolor" ||
            name == "facecolor" || 
			name == "markerfacecolor"){
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
		return true;
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
        Object* o = getObject(h);
        return o != nullptr;
    }

    bool Object::isFigure(){
        return false;
    }

    bool Object::isAxes(){
        return false;
    }

    bool Object::isDrawable()
    {
        return false;
    }

    bool Object::isText()
    {
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

    ObjectType Object::getObjectType() const
    {
        return m_type;
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

    Object* Object::getParentObject()
    {
        double parent = getPropertyValue("parent").Scalar();
        return getObject(parent);
    }

    Root::Root()
        :Object()
    {
        m_type = ObjectType::ROOT;

        m_ps.push_back(Property("type", string("root"), PropertyType::STRING) );
        m_ps.push_back(Property("handle", double(0), PropertyType::DOUBLE) );
        m_ps.push_back(Property("currentfigure", double(0), PropertyType::DOUBLE) );

        m_ps.push_back(Property("units", string("pixels"), PropertyType::STRING));

        m_objectMap[double(0)] = this;
    }

    void Root::repaint(){
        update(nullptr);
    }

    void Root::update(GnuplotOutput *out){
        vector<double> children = getPropertyValue("children").Vector();
        vector<double>::iterator it = children.begin();
        for ( ; it != children.end(); ++it){
            Object* o = getObject(*it);
            if (o && o->isFigure())
                static_cast<Figure*>(o)->update();
        }
    }

    Figure::Figure()
        :Object(), m_out(new GnuplotOutput), m_tLabel(new Text), m_rLabel(new Text), 
        m_lLabel(new Text), m_bLabel(new Text), m_rows(1), m_cols(1), m_isGridLayout(false)
    {
        m_type = ObjectType::FIGURE;

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
        m_ps.push_back(Property("toplabel", Currency(m_tLabel.get(), "text"), PropertyType::POINTER));
        m_ps.push_back(Property("bottomlabel", Currency(m_bLabel.get(), "text"), PropertyType::POINTER));
        m_ps.push_back(Property("leftlabel", Currency(m_lLabel.get(), "text"), PropertyType::POINTER));
        m_ps.push_back(Property("rightlabel", Currency(m_rLabel.get(), "text"), PropertyType::POINTER));
        m_tLabel->setPropertyValue("visible", "off");
        m_bLabel->setPropertyValue("visible", "off");
        m_rLabel->setPropertyValue("visible", "off");
        m_lLabel->setPropertyValue("visible", "off");

        m_objectMap[handle] = this;
        m_out->setWindowId(handle);
        m_out->printf("clear\n");
        m_out->flush();

		// not yet supported properties		
		m_ps.push_back(Property("colormap", string("colormap"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("name", string("name"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("numbertitle", string("numbertitle"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("uicontextmenu", string("uicontextmenu"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("windowstyle", string("windowstyle"), PropertyType::UNSUPPORTED));
    }

    Figure::~Figure(){
        int fh = int(getHandle());
        m_figureHandlePool->releaseHandle(fh);
    }

    bool Figure::isFigure(){
        return true;
    }

    void Figure::repaint(){
        update(nullptr);
    }

    void Figure::update(GnuplotOutput *out){
        m_out->printf("clear\n");
        m_out->printf("unset label\n");
        m_out->printf("set multiplot\n");
        m_out->flush();

        m_tLabel->update(m_out.get(), "tlabel");
        m_rLabel->update(m_out.get(), "rlabel");
        m_lLabel->update(m_out.get(), "llabel");
        m_bLabel->update(m_out.get(), "blabel");

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

    bool Figure::setPropertyValue(const string& name, VALUETYPE value)
    {
        bool ret = false;
        if (name == "toplabel") {
            ret = m_tLabel->setPropertyValue("string", value);
            m_tLabel->setPropertyValue("visible", "on");
        }
        else if (name == "leftlabel") {
            ret = m_lLabel->setPropertyValue("string", value);
            m_lLabel->setPropertyValue("visible", "on");
        }
        else if (name == "rightlabel") {
            ret = m_rLabel->setPropertyValue("string", value);
            m_rLabel->setPropertyValue("visible", "on");
        }
        else if (name == "bottomlabel") {
            ret = m_bLabel->setPropertyValue("string", value);
            m_bLabel->setPropertyValue("visible", "on");
        }
        else {
            ret = Object::setPropertyValue(name, value);
        }
        return ret;
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
        
        m_tLabel->setPropertyValue("string", "");
        m_bLabel->setPropertyValue("string", "");
        m_rLabel->setPropertyValue("string", "");
        m_lLabel->setPropertyValue("string", "");
        m_tLabel->setPropertyValue("visible", "off");
        m_bLabel->setPropertyValue("visible", "off");
        m_rLabel->setPropertyValue("visible", "off");
        m_lLabel->setPropertyValue("visible", "off");
       
        m_rows = 1;
        m_cols = 1;
        m_isGridLayout = false;

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

    void Figure::setGridLayout(int rows, int cols)
    {
        m_rows = rows;
        m_cols = cols;
        m_isGridLayout = true;
    }

    bool Figure::isGridLayout()
    {
        return m_isGridLayout;
    }

    void Figure::getSubplotPosition(int subplotIdx, double& x, double& y, double& w, double& h)
    {
        double wMult = 1.0;
        double xOffset = 0.0;
        if (m_lLabel->getVisible()) {
            wMult -= 0.05;
            xOffset = 0.05;
        }
        if (m_rLabel->getVisible()) {
            wMult -= 0.05;
        }
        double hMult = 1.0;
        double yOffset = 0.0;
        if (m_bLabel->getVisible()) {
            hMult -= 0.05;
            yOffset = 0.05;
        }
        if (m_tLabel->getVisible()) {
            hMult -= 0.05;
        }

        int r = (subplotIdx - 1) / m_cols;
        int c = (subplotIdx - 1) % m_cols;
        int rr = m_rows - 1 - r;
        
        w = wMult / m_cols;
        h = hMult / m_rows;
        x = xOffset + c * w;
        y = yOffset + rr * h;
    }

    void Figure::out(string s){
        m_out->printf("%s\n", s.c_str());
        m_out->flush();
    }

    Axes::Axes()
        :Object(), m_title(new Text), m_xlabel(new Text), m_ylabel(new Text), m_zlabel(new Text),
        _borderOn(true), _colorbarVisible(true), m_barlayout("grouped"), m_barWidth(0.8), 
        m_secYAxisVisible(false), m_updateXAxisRange(false), m_updateYAxisRange(false),
        m_updateY2AxisRange(false), m_updateZAxisRange(false), m_gridIndex(1), m_tailX(0.0),
        m_tailY(0.0), m_tailTheta(0.0), m_tailR(0.0)
    {
        m_type = ObjectType::AXES;

        double handle = m_handlePool->allocHandle();
        m_secYAxis.reset(new SecondaryYAxis(handle));
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
        m_ps.push_back(Property("fontangle", string("regular"), PropertyType::STRING));
        m_ps.push_back(Property("fontname", string("arial"), PropertyType::STRING));
        m_ps.push_back(Property("fontsize", 7.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("fontweight", string("normal"), PropertyType::STRING));
        m_ps.push_back(Property("xnumericformat", string("auto"), PropertyType::STRING));
        m_ps.push_back(Property("xnumericprecision", 5.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("ynumericformat", string("auto"), PropertyType::STRING));
        m_ps.push_back(Property("ynumericprecision", 5.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("znumericformat", string("auto"), PropertyType::STRING));
        m_ps.push_back(Property("znumericprecision", 5.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("barlabels", string("off"), PropertyType::STRING));
        m_ps.push_back(Property("barlabelsfontname", string("arial"), PropertyType::STRING));
        m_ps.push_back(Property("barlabelsfontsize", 7.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("barlabelsfontweight", string("normal"), PropertyType::STRING));
        m_ps.push_back(Property("barlabelsfontangle", string("regular"), PropertyType::STRING));
        m_ps.push_back(Property("bargap", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("xcategories", string(""), PropertyType::CELL));
        m_ps.push_back(Property("colororder", Currency(), PropertyType::MATRIX));
        m_ps.push_back(Property("colorbarscale", string("linear"), PropertyType::STRING));
        // default colormap 
        hwMatrix* map = new hwMatrix(4, 3, hwMatrix::DataType::REAL);
        (*map)(0, 0) = 0; (*map)(0, 1) = 0; (*map)(0, 2) = 255;
        (*map)(1, 0) = 0; (*map)(1, 1) = 255; (*map)(1, 2) = 0;
        (*map)(2, 0) = 255; (*map)(2, 1) = 255; (*map)(2, 2) = 0;
        (*map)(3, 0) = 255; (*map)(3, 1) = 0; (*map)(3, 2) = 0;
        m_ps.push_back(Property("colormap", Currency(map), PropertyType::MATRIX));
        //vector<double> cl;
        //pos.push_back(256);
        m_ps.push_back(Property("colorlevels", Currency(), PropertyType::VEC_DOUBLE));
        m_ps.push_back(Property("secondaryyaxis", Currency(m_secYAxis.get(),"yaxis"), PropertyType::POINTER));
        m_ps.push_back(Property("polarmethod", string("phasevsmag"), PropertyType::STRING));
        m_ps.push_back(Property("polartiptotail", string("off"), PropertyType::STRING));

        m_objectMap[handle] = this;
        m_title->setParent(this);

		// not yet supported properties
		m_ps.push_back(Property("xtick", string("xtick"), PropertyType::UNSUPPORTED));
		m_ps.push_back(Property("ytick", string("ytick"), PropertyType::UNSUPPORTED));
		m_ps.push_back(Property("ztick", string("ztick"), PropertyType::UNSUPPORTED));
		m_ps.push_back(Property("mouseclickcallback", string("mouseclickcallback"), PropertyType::UNSUPPORTED));
		m_ps.push_back(Property("contourtype", string("contourtype"), PropertyType::UNSUPPORTED));
		m_ps.push_back(Property("ycategories", string("ycategories"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("barorientation", string("barorientation"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xtickmethod", string("xtickmethod"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("ytickmethod", string("ytickmethod"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("ztickmethod", string("ztickmethod"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xaxislocation", string("xaxislocation"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("yaxislocation", string("yaxislocation"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("colorbar", string("colorbar"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("legend", string("legend"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("zerolinecolor", string("zerolinecolor"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("framecolor", string("framecolor"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("gridcolor", string("gridcolor"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xminortick", string("xminortick"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("yminortick", string("yminortick"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("zminortick", string("zminortick"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xdb10reference", string("xdb10reference"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("ydb10reference", string("ydb10reference"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xdb20reference", string("xdb20reference"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("ydb20reference", string("ydb20reference"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xweighting", string("xweighting"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("yweighting", string("yweighting"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xticklabel", string("xticklabel"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("yticklabel", string("yticklabel"), PropertyType::UNSUPPORTED));
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

        _borderOn = true;
        _colorbarVisible = true;
        _colorbarRange.clear();

        m_updateXAxisRange = false;
        m_updateYAxisRange = false;
        m_updateY2AxisRange = false;
        m_updateZAxisRange = false;
        m_secYAxisVisible = false;
        repaint();
    }

    void Axes::repaint(){
        vector<double> parent = getPropertyValue("parent").Vector();
        Object* o = getObject(parent[0]);
        if (o && o->isFigure())
            static_cast<Figure*>(o)->repaint();
    }        

    void Axes::update(GnuplotOutput *out){
        string v = this->getPropertyValue("visible").StringVal();
        if (v == "off"){
            return;
        }

        // reset grids and tics
        out->printf("set border 31\nset xtics mirror\n");
        out->printf("unset polar\nset size nosquare\nset grid nopolar\n");
        out->printf("set grid noxtics\nset grid noytics\nset grid noztics\n");
        out->printf("unset mxtics\nunset mytics\nunset mztics\n");

        // reset all ranges - if this is not done here then previous range will remain even if
        // the plot is cleared
        out->printf("unset xrange\n");
        out->printf("unset yrange\n");
        out->printf("unset y2range\n");
        out->printf("unset zrange\n");
        out->printf("unset y2tics\nunset y2label\n");

        // set the range only for the axes that a range has been specified
        if (m_updateXAxisRange) {
            vector<double> xrange = getPropertyValue("xlim").Vector();
            out->printf("set xrange [%g:%g]\n", xrange[0], xrange[1]);
        }
        if (m_updateYAxisRange) {
            vector<double> yrange = getPropertyValue("ylim").Vector();
            out->printf("set yrange [%g:%g]\n", yrange[0], yrange[1]);
        }
        if (m_updateZAxisRange) {
            vector<double> zrange = getPropertyValue("zlim").Vector();
            out->printf("set zrange [%g:%g]\n", zrange[0], zrange[1]);
        }
        if (m_updateY2AxisRange) {
            vector<double> y2range = m_secYAxis->getPropertyValue("ylim").Vector();
            out->printf("set y2range [%g:%g]\n", y2range[0], y2range[1]);
        }

        string bgcolor = getPropertyValue("color").ColorString();
        out->printf("set obj 1 rectangle from graph 0,0 to graph 1,1 behind fc rgb '%s' fs border rgb '%s'\n", bgcolor.c_str(), bgcolor.c_str());
        
        string tickFontName = getPropertyValue("fontname").StringVal();
        double tickFontSize = getPropertyValue("fontsize").Scalar();
        string tickFontWeight = getPropertyValue("fontweight").StringVal();
        string tickFontAngle = getPropertyValue("fontangle").StringVal();

        string boldString = tickFontWeight == "bold" ? "/:Bold" : "";
        string italicString = tickFontAngle == "italic" ? "/:Italic" : "";

        string xcolor = getPropertyValue("xcolor").ColorString();
        string xnf = getPropertyValue("xnumericformat").StringVal();
        int xnp = int(getPropertyValue("xnumericprecision").Scalar());
        string f = "g";
        if (xnf == "scientific")
            f = "e";
        else if (xnf == "fixed")
            f = "f";
        char numF[24];
        sprintf(numF, "%s.%d%s", "%", xnp, f.c_str());
        out->printf("set xtics format \"{%s {%s %s}}\" textcolor rgb '%s' font \"%s,%g\" \n", boldString.c_str(), 
            italicString.c_str(), numF, xcolor.c_str(), tickFontName.c_str(), tickFontSize);
        
        string ycolor = getPropertyValue("ycolor").ColorString();
        string ynf = getPropertyValue("ynumericformat").StringVal();
        int ynp = int(getPropertyValue("ynumericprecision").Scalar());
        f = "g";
        if (ynf == "scientific")
            f = "e";
        else if (ynf == "fixed")
            f = "f";
        sprintf(numF, "%s.%d%s", "%", ynp, f.c_str());
        out->printf("set ytics format \"{%s {%s %s}}\" textcolor rgb '%s' font \"%s,%g\" \n", boldString.c_str(),
            italicString.c_str(), numF, ycolor.c_str(), tickFontName.c_str(), tickFontSize);
        
        
        string zcolor = getPropertyValue("zcolor").ColorString();
        string znf = getPropertyValue("znumericformat").StringVal();
        int znp = int(getPropertyValue("znumericprecision").Scalar());
        f = "g";
        if (znf == "scientific")
            f = "e";
        else if (znf == "fixed")
            f = "f";
        sprintf(numF, "%s.%d%s", "%", znp, f.c_str());
        out->printf("set ztics format \"{%s {%s %s}}\" textcolor rgb '%s' font \"%s,%g\" \n", boldString.c_str(),
            italicString.c_str(), numF, zcolor.c_str(), tickFontName.c_str(), tickFontSize);

        string xscale = getPropertyValue("xscale").StringVal();
        if (xscale == "log"){
            out->printf("set logscale x\n");
        }
        else {
            out->printf("unset logscale x\n");
        }
        string yscale = getPropertyValue("yscale").StringVal();
        if (yscale == "log"){
            out->printf("set logscale y\n");
        }
        else {
            out->printf("unset logscale y\n");
        }
        string zscale = getPropertyValue("zscale").StringVal();
        if (zscale == "log"){
            out->printf("set logscale z\n");
        }
        else {
            out->printf("unset logscale z\n");
        }

        m_title->update(out, "title");
        m_xlabel->update(out, "xlabel");
        m_ylabel->update(out, "ylabel");
        m_zlabel->update(out, "zlabel");

        if (m_secYAxisVisible)
            m_secYAxis->update(out);

        vector<double> pos = getPropertyValue("position").Vector();
        double x = pos[0];
        double y = pos[1];
        double w = pos[2];
        double h = pos[3];
        Object* pobj = getParentObject();
        if (pobj && pobj->isFigure()) {
            Figure* fig = static_cast<Figure*>(pobj);
            if (fig->isGridLayout() || (x == 0 && y == 0 && w == 1 && h == 1)) {
                fig->getSubplotPosition(m_gridIndex, x, y, w, h);
            }
        }
        
        out->printf("set origin %f,%f\n", x, y);
        out->printf("set size %f,%f\n", w, h);
        {
            out->printf("unset key\nplot '-' with points ps 0\n0 0\n1 1\ne\n");
        }

        string xgrid = getPropertyValue("xgrid").StringVal();
        if (xgrid == "on") {
            out->printf("set grid xtics\n");
        }

        string xmgrid = getPropertyValue("xminorgrid").StringVal();
        if (xmgrid == "on") {
            out->printf("set mxtics 5\nset grid mxtics\n");
        }
        
        string ygrid = getPropertyValue("ygrid").StringVal();
        if (ygrid == "on") {
            out->printf("set grid ytics\n");
        }
        
        string ymgrid = getPropertyValue("yminorgrid").StringVal();
        if (ymgrid == "on") {
            out->printf("set mytics 5\nset grid mytics\n");
        }
        
        string zgrid = getPropertyValue("zgrid").StringVal();
        if (zgrid == "on") {
            out->printf("set grid ztics\n");
        }
        
        string zmgrid = getPropertyValue("zminorgrid").StringVal();
        if (zmgrid == "on") {
            out->printf("set mztics 5\nset grid mztics\n");
        }
        
        string lState = getPropertyValue("legend").StringVal();
        if (lState == "on") {
            out->printf("set key on\n");
        }
        else {
            out->printf("set key off\n");
        }

        if (_borderOn)
            out->printf("set border 31\nset xtics mirror\n");
        else if (m_secYAxisVisible)
            out->printf("set border 43\nset xtics nomirror\n");
        else
            out->printf("set border 3\nset xtics nomirror\n");

        if (m_secYAxisVisible || !_borderOn)
            out->printf("set ytics nomirror\n");
        else
            out->printf("set ytics mirror\n");
        
        vector<double> children = getPropertyValue("children").Vector();
        Drawable* pChild = nullptr;
        vector<string> linestyle;
        vector<string> usingclause;
        vector<string> legend;
        vector<string> withclause;
        vector<vector<double> > zdata;

        bool isBarPlot = false;
        bool isPolarPlot = false;
        int lineId = 1;
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it) {
            Object* obj = getObject(*it);
            if (obj->isDrawable()) {
                
                pChild = static_cast<Drawable*>(obj);
                if (!pChild->getVisible())
                    continue;

                linestyle.push_back(pChild->getLineStyle(lineId));
                usingclause.push_back(pChild->getUsingClause());
                legend.push_back(pChild->getLegend());
                withclause.push_back(pChild->getWithClause(lineId));
                zdata.push_back(pChild->getZdata());

                if (pChild->getObjectType() == ObjectType::BAR) {
                    isBarPlot = true;
                }
                else if (pChild->getObjectType() == ObjectType::POLAR_LINE) {
                    isPolarPlot = true;
                }
            }
            else if (obj->isText()) {
                Text* pChild = static_cast<Text*>(obj);
                if (pChild->getVisible()) {
                    pChild->update(out);
                }
            }
            ++lineId;
        }

        bool is3D = zdata.size() > 0 && zdata[0].size() > 0;

        if (is3D && _colorbarVisible){
            setupColorbar(out);
        }
        else
        {
            out->printf("unset colorbox\n");
        }
        // common style for bar
        if (isBarPlot) {
            // get bar width from the first BAR child
            double barwidth = 0.8;
            for (it = children.begin(); it != children.end(); ++it) {
                Object* obj = getObject(*it);
                if (obj->getObjectType() == ObjectType::BAR){
                    barwidth = obj->getPropertyValue("barwidth").Scalar();
                    break;
                }
            }
            string layout = "cluster";
            if (m_barlayout == "stacked")
                layout = "rowstacked";
            int bargap = int(getPropertyValue("bargap").Scalar());

            out->printf("set yrange [*<0:]\n set style data histogram\n");
            out->printf("set style histogram %s gap %g \n", layout.c_str(), bargap);
            out->printf("set boxwidth %g \n", barwidth);
        }

        // common style for polar
        if (isPolarPlot) {
            out->printf("set polar\n set ttics 45\nset grid polar 45\nset rtics\n");
            out->printf("set size square\nunset border\n");
            out->printf("unset xtics\nunset ytics\n");
        }

        size_t lineCount = usingclause.size();
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
            string cmd = ss.str();
            out->printf(cmd.c_str());

            string labelsData;
            if (isBarPlot && getPropertyValue("barlabels").StringVal() == "on" &&
                m_barlayout == "stacked") {
                labelsData = drawStackedBarPlotLabels(out);
            }

            out->printf("\n");

            vector<double>::iterator it = children.begin();
            for (; it != children.end(); ++it){
                Object* obj = getObject(*it);
                if (obj && obj->isDrawable()) {
                    pChild = static_cast<Drawable*>(obj);
                    if (pChild->getVisible())
                        pChild->putData(out);
                }
            }

            if (!labelsData.empty())
                out->printf(labelsData);
        }
        // reset tip-to-tail helper coordinates
        m_tailX = 0;
        m_tailY = 0;
        m_tailR = 0;
        m_tailTheta = 0;

        out->printf("set title ''\n");
        out->printf("set xlabel ''\n");
        out->printf("set ylabel ''\n");
        out->printf("set zlabel ''\n");
        out->printf("set key off\n");
        out->printf("\n\n");
    }

    bool Axes::isAxes(){
        return true;
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

    void Axes::setAxisNeedRepaint(int axisID){
        switch (axisID) {
        case 0:
            m_updateXAxisRange = true;
            break;
        case 1:
            m_updateYAxisRange = true;
            break;
        case 2:
            m_updateY2AxisRange = true;
            break;
        case 3:
            m_updateZAxisRange = true;
            break;
        }
    }

    void Axes::setBorder(bool state)
    {
        _borderOn = state;
    }

    bool Axes::getBorder() const
    {
        return _borderOn;
    }

    void Axes::setColorbarVisible(bool state)
    {
        _colorbarVisible = state;
    }

    bool Axes::getColorbarVisible() const
    {
        return _colorbarVisible;
    }

    void Axes::setColorbarRange(const std::vector<double>& range)
    {
        _colorbarRange.clear();
        _colorbarRange = range;
    }

    std::vector<double> Axes::getColorbarRange()
    {
        if (_colorbarRange.size() == 2)
            return _colorbarRange;
        // else get the values from the curves
        std::vector<double> ret;
        vector<double> children = getPropertyValue("children").Vector();
        for (vector<double>::iterator it = children.begin(); it != children.end(); ++it) 
        {
            Object* obj = getObject(*it);
            ObjectType t = obj->getObjectType();
            if (t == ObjectType::SURFACE || t == ObjectType::MESH ||
                t == ObjectType::CONTOUR_3D || t == ObjectType::CONTOUR) {
                
                Surface* pChild = static_cast<Surface*>(obj);
                double min, max;
                pChild->getMinMaxZ(min, max);
                if (ret.empty())
                {
                    ret.push_back(min);
                    ret.push_back(max);
                }
                else
                {
                    if (min < ret[0])
                        ret[0] = min;
                    if (max > ret[1])
                        ret[1] = max;
                }
            }
        }
        if (ret.empty())
        {
            ret.push_back(0);
            ret.push_back(10);
        }
        return ret;
    }

    void Axes::getBarNumberAndIdx(double barHandle, int& barNumber, int& barIndex)
    {
        barNumber = 0;
        vector<double> children = getPropertyValue("children").Vector();
        for (vector<double>::const_iterator it = children.cbegin(); it != children.cend(); ++it) {
            Object* obj = getObject(*it);
            if (obj->getObjectType() == ObjectType::BAR &&
                static_cast<HggroupBar*>(obj)->getVisible()) {
                barNumber++;
                if (barHandle == (*it))
                    barIndex = barNumber - 1;
            }
        }
    }

    void Axes::setBarLayout(const string& barLayout)
    {
        if (m_barlayout == barLayout)
            return;
        m_barlayout = barLayout;
        vector<double> children = getPropertyValue("children").Vector();
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it) {
            Object* obj = getObject(*it);
            if (obj->getObjectType() == ObjectType::BAR) {
                obj->setPropertyValue("barlayout", barLayout);
            }
        }
    }

    void Axes::setBarWidth(double barWidth)
    {
        if (m_barWidth == barWidth)
            return;

        m_barWidth = barWidth;
        vector<double> children = getPropertyValue("children").Vector();
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it) {
            Object* obj = getObject(*it);
            if (obj->getObjectType() == ObjectType::BAR) {
                obj->setPropertyValue("barwidth", barWidth);
            }
        }
    }

    void Axes::setSecondaryYAxisVisible(bool visible)
    {
        m_secYAxisVisible = visible;
    }

    bool Axes::setPropertyValue(const string& name, VALUETYPE value)
    {
        bool ret = Object::setPropertyValue(name, value);
        if (name == "fontname" || name == "fontsize" || name == "fontweight" || name == "fontangle") {
            m_secYAxis->setPropertyValue(name, value);
        }
        else if (name == "polartiptotail") {
            vector<double> children = getPropertyValue("children").Vector();
            vector<double>::iterator it = children.begin();
            for (; it != children.end(); ++it) {
                Object* obj = getObject(*it);
                if (obj->getObjectType() == ObjectType::POLAR_LINE) {
                    obj->setPropertyValue("tiptotail", value);
                }
            }
        }
        if (name == "xlim")
            setAxisNeedRepaint(0);
        else if (name == "ylim")
            setAxisNeedRepaint(1);
        else if (name == "zlim")
            setAxisNeedRepaint(3);
        return ret;
    }

    double Axes::getSecondaryYAxisHandle() const
    {
        return m_secYAxis->getHandle();
    }

    void Axes::setGridIndex(int idx)
    {
        m_gridIndex = idx;
    }

    void Axes::getTipToTailCoordinates(double& x, double& y, double& t, double& r) const
    {
        x = m_tailX;
        y = m_tailY;
        r = m_tailR;
        t = m_tailTheta;
    }

    void Axes::setTipToTailCoordinates(double x, double y, double t, double r)
    {
        m_tailX = x;
        m_tailY = y;
        m_tailR = r;
        m_tailTheta = t;
    }

    string Axes::drawStackedBarPlotLabels(GnuplotOutput* out)
    {
        stringstream ss;
        bool sumInitialized = false;
        vector<double> dataX;
        vector<double> sumY;

        vector<double> children = getPropertyValue("children").Vector();
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it) {
            Object* obj = getObject(*it);
            if (obj->getObjectType() == ObjectType::BAR) {
                HggroupBar* pChild = static_cast<HggroupBar*>(obj);
                if (pChild && pChild->getVisible()) {
                    vector<double> x = pChild->getXdata();
                    vector<double> y = pChild->getYdata();
                    if (!sumInitialized)
                    {
                        dataX = x;
                        sumY = y;
                        sumInitialized = true;
                    }
                    else
                    {
                        for (int i = 0; i < min(sumY.size(), y.size()); ++i)
                            sumY[i] += y[i];
                    }
                }
            }
        }

        string boldString = getPropertyValue("barlabelsfontweight").StringVal() == "bold" ? "/:Bold" : "";
        string italicString = getPropertyValue("barlabelsfontangle").StringVal() == "italic" ? "/:Italic" : "";

        double maxVal = 0;
        size_t count = min(dataX.size(), sumY.size());
        for (int j = 0; j < count; j++) {

            ss << dataX[j] << " " << sumY[j] 
               << " \"{" << boldString << "{" << italicString << "{" << sumY[j] << "}}}\"\n";

            if (maxVal < sumY[j] * 0.1)
                maxVal = sumY[j] * 0.1;
        }
        ss << "e\n";
        out->printf(",'' using ($0):($2 + %g):($3) with labels title \"\" font \"%s,%g\"\n", maxVal,
            getPropertyValue("barlabelsfontname").StringVal().c_str(),
            getPropertyValue("barlabelsfontsize").Scalar());

        return ss.str();
    }

    void Axes::setupColorbar(GnuplotOutput* out) {
        out->printf("set colorbox\n");

        string colorbarScale = getPropertyValue("colorbarscale").StringVal();
        if (colorbarScale == "log")
            out->printf("set log cb\nset cbtics 2\n");
        else
            out->printf("unset log cb\n");


        vector<double> colorlevels;
        Currency c = getPropertyValue("colorlevels").getCurrency();
        if (c.IsVector())
            colorlevels = c.Vector();

        vector<double> paletteTics;
        if (colorlevels.size() == 1) {
            out->printf("set palette maxcolors %g\n", colorlevels.front());
        }
        else if (colorlevels.size() > 1) {
            double range = colorlevels.back() - colorlevels.front();
            paletteTics.push_back(0);
            for (int i = 1; i < colorlevels.size() - 1; ++i) {
                double val = paletteTics[i-1] + (colorlevels[i] - colorlevels[i - 1]) / range;
                paletteTics.push_back(val);
            }
            paletteTics.push_back(1);
        }
        
        Currency colormap = getPropertyValue("colormap").getCurrency();
        vector<string> colorStrs;
        if (colormap.IsMatrix()) {
            const hwMatrix* colorMat = colormap.Matrix();
            int numColors = colorMat->M();
            int numColComponents = colorMat->N();

            if (numColComponents == 3) {
                for (int i = 0; i < numColors; ++i) {
                    double r = (*colorMat)(i, 0);
                    double g = (*colorMat)(i, 1);
                    double b = (*colorMat)(i, 2);
                    if (r <= 1 && g <= 1 && b <= 1) {
                        r *= 255; g *= 255; b *= 255;
                    }
                    vector<int> c = { (int)r,(int)g,(int)b };
                    Color cl(c);
                    colorStrs.push_back("\""+cl.getString()+"\"");
                }
            }
        }


        if (colorStrs.empty()) {
            colorStrs.push_back(string("\"blue\""));
            colorStrs.push_back(string("\"green\""));
            colorStrs.push_back(string("\"yellow\""));
            colorStrs.push_back(string("\"red\""));
        }
        if (paletteTics.empty()) {
            for (int i = 0; i < colorStrs.size(); ++i)
                paletteTics.push_back(i);
        }

        out->printf("set palette defined (");
        int sizeCM = int(colorStrs.size());
        for (int i = 0; i < paletteTics.size(); ++i) {
            out->printf("%g %s", paletteTics[i], colorStrs[i % sizeCM].c_str());
            if (i < (paletteTics.size() - 1))
                out->printf(",");
        }
        out->printf(")\n");

        if (_colorbarRange.size() == 2)
            out->printf("set cbrange [%g:%g]\n", _colorbarRange[0], _colorbarRange[1]);
        else
            out->printf("unset cbrange\n");


    }

    Drawable::Drawable()
        :Object(), _xcolcount(0), _ycolcount(0), _zcolcount(0), 
        m_xaxisRef("x1"), m_yaxisRef("y1")
    {
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

    bool Drawable::isDrawable()
    {
        return true;
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

    bool Drawable::getVisible()
    {
        return getPropertyValue("visible").StringVal() == "on";
    }

    void Drawable::putData(GnuplotOutput *out){
        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> z = getZdata();

        double xOffset = getPropertyValue("dataxoffset").Scalar();
        double yOffset = getPropertyValue("datayoffset").Scalar();
        double xScale = getPropertyValue("dataxscale").Scalar();
        double yScale = getPropertyValue("datayscale").Scalar();

        if (z.size() == 0){     // 2d data
            size_t count = min(x.size(), y.size());
            for (int j = 0; j < count; j++){
                out->printf("%g %g\n", x[j] * xScale + xOffset, y[j] * yScale + yOffset);
            }
        } else {                // 3d data
            double zOffset = getPropertyValue("datazoffset").Scalar();
            double zScale = getPropertyValue("datazscale").Scalar();

            size_t count = min(x.size(), y.size());
            count = min((int)count, (int)z.size());
            for (int j = 0; j < count; j++){
                out->printf("%g %g %g\n", x[j] * xScale + xOffset, y[j] * yScale + yOffset, z[j] * zScale +zOffset);
            }
        }
        // end of data
        out->printf("e\n");
    }

    void Drawable::plotOnSecondaryYAxis(bool plot)
    {
        m_yaxisRef = plot ? string("y2") : string("y1");
    }

    Line::Line()
        :Drawable()
    {
        m_type = ObjectType::LINE;

        m_ps.push_back(Property("type", string("line"), PropertyType::STRING) );
        m_ps.push_back(Property("linestyle", string("-"), PropertyType::STRING));
        m_ps.push_back(Property("color", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("linewidth", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("marker", string(""), PropertyType::STRING));
        m_ps.push_back(Property("markerfacecolor", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("markersize", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("dataxoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("dataxscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datazoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datazscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("markerevery", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("tiptotail", string("off"), PropertyType::STRING));
		
        // not yet supported properties
        m_ps.push_back(Property("units", string("units"), PropertyType::UNSUPPORTED));
    }

    void Line::init(const LineData &ld){
        Drawable::init(ld);

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        // update line-specific properties
        LineStyle ls = LineStyle(ld, colorOrder);
        setPropertyValue("linestyle", ls.m_lineStyle);
        setPropertyValue("color", ls.m_lineColor);
        setPropertyValue("linewidth", ls.m_lineWidth);
        setPropertyValue("marker", ls.m_markerStyle);
        setPropertyValue("markerfacecolor", ls.m_markerColor);
        setPropertyValue("markersize", ls.m_markerSize);
        setPropertyValue("displayname", ls.m_legend);

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++){
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string Line::getUsingClause(){
        stringstream ss;
        ss << "'-' using 1:2 axis " << m_xaxisRef << m_yaxisRef;
        return ss.str();
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
        int markerEvery = int(getPropertyValue("markerevery").Scalar());

        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line << " dt " << dt
           << " lc rgb '" << lineColor << "' lw " << lineWidth
           << " pt " << pt << " ps " << markerSize
           << " pi " << markerEvery
            //<< " tc rgb '" << markerColor << "'"
           << "\n";
        return ss.str();
    }

    void Line::cleanup(GnuplotOutput *out){
    }

    Polar::Polar()
        :Line()
    {
        m_type = ObjectType::POLAR_LINE;
    }

    void Polar::putData(GnuplotOutput* out)
    {
        vector<double> x = getXdata();
        vector<double> y = getYdata();

        if (x.empty() || y.empty())
        {
            out->printf("e\n");
            return;
        }

        double xOffset = getPropertyValue("dataxoffset").Scalar();
        double yOffset = getPropertyValue("datayoffset").Scalar();
        double xScale = getPropertyValue("dataxscale").Scalar();
        double yScale = getPropertyValue("datayscale").Scalar();

        string polarMethod;
        bool tipToTail = false;
        Axes* axes = nullptr;
        Object* obj = getParentObject();
        if (obj && obj->getObjectType() == ObjectType::AXES) {
            polarMethod = obj->getPropertyValue("polarmethod").StringVal();
            tipToTail = obj->getPropertyValue("polartiptotail").StringVal() == "on";
            axes = static_cast<Axes*>(obj);
        }
        if (polarMethod == "vector") {
            
            if (tipToTail &&  getPropertyValue("tiptotail").StringVal() == "on") {

                double xp, yp, rp, tp;
                axes->getTipToTailCoordinates(xp, yp, rp, tp);
                
                double r = y[0] * yScale + yOffset;
                double t = x[0] * xScale + xOffset;
                double x = r * cos(t);
                double y = r * sin(t);

                double x1 = xp + x;
                double y1 = yp + y;
                double r1 = sqrt(x1 * x1 + y1 * y1);
                double t1 = 0.0;

                if (y1 == 0)
                    t1 = 0;
                else if (x1 == 0 && y1 < 0)
                    t1 = PI / 2;
                else if (x1 == 0 && y1 > 0)
                    t1 = -PI / 2;
                else
                    t1 = atan(y1 / x1);

                if (x1 < 0)
                    t1 += PI;

                // start of vector
                out->printf("%g %g\n", rp, tp);
                // end of vector
                out->printf("%g %g\n", t1, r1);
                out->printf("e\n");

                axes->setTipToTailCoordinates(x1, y1, t1, r1);
            }
            else {
                size_t count = min(x.size(), y.size());
                out->printf("0 0\n");
                out->printf("%g %g\n", x[0] * xScale + xOffset, y[0] * yScale + yOffset);
                out->printf("e\n");
            }
        }
        else if (polarMethod == "radar") {

            double maxTheta = x.front();
            for (vector<double>::const_iterator it = x.cbegin(); it != x.cend(); ++it) {
                if (maxTheta < (*it))
                    maxTheta = (*it);
            }
            if (maxTheta == 0)
                maxTheta = 1;
            vector<double>::const_iterator yit = y.cbegin();
            for (vector<double>::const_iterator it = x.cbegin(); it != x.cend() && yit != y.cend(); ++it, ++yit) {
                out->printf("%g %g\n", 2* PI * (*it)/maxTheta, (*yit));
            }
            out->printf("e\n");
        }
        else {
            Drawable::putData(out);
        }
    }

    Fill::Fill()
        :Drawable()
    {
        m_type = ObjectType::FILL;

        m_ps.push_back(Property("facecolor", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("type", string("patch"), PropertyType::STRING) );
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("dataxoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("dataxscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datazoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datazscale", 1.0, PropertyType::DOUBLE));

        // not yet supported properties
        m_ps.push_back(Property("units", string("units"), PropertyType::UNSUPPORTED));
    }

    void Fill::init(const LineData &ld){
        Drawable::init(ld);

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        LineStyle ls = LineStyle(ld, colorOrder);
        setPropertyValue("facecolor", ls.m_lineColor);

        size_t propSize = ld.properties.size();
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
        string lineColor = getPropertyValue("facecolor").ColorString();

        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line << " linetype " << dt
           << " lc rgb '" << lineColor
           << "\n";
        return ss.str();
    }

    void Fill::cleanup(GnuplotOutput *out){
    }

    Area::Area()
        :Drawable()
    {
        m_type = ObjectType::AREA;

        m_ps.push_back(Property("type", string("hggroup"), PropertyType::STRING) );
        m_ps.push_back(Property("basevalue", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("facecolor", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("dataxoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("dataxscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("areagroup", -1.0, PropertyType::DOUBLE));

        // not yet supported properties
        m_ps.push_back(Property("units", string("units"), PropertyType::UNSUPPORTED));
    }

    void Area::init(const LineData& ld){
        Drawable::init(ld);

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        LineStyle ls = LineStyle(ld, colorOrder);
        setPropertyValue("facecolor", ls.m_lineColor);
        setPropertyValue("areagroup", getPropertyValue("handle").Scalar());

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
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

    HggroupBar::HggroupBar()
        :Drawable()
    {
        m_type = ObjectType::BAR;

        m_ps.push_back(Property("type", string("hggroup"), PropertyType::STRING) );
        m_ps.push_back(Property("barwidth", double(0.8), PropertyType::DOUBLE));
        m_ps.push_back(Property("barlayout", string("grouped"), PropertyType::STRING));
        m_ps.push_back(Property("facecolor", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("barstyle", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("bargroup", -1.0, PropertyType::DOUBLE));

        // not yet supported properties
        m_ps.push_back(Property("units", string("units"), PropertyType::UNSUPPORTED));
    }

    void HggroupBar::init(const LineData &ld){
        Drawable::init(ld);

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        // update line-specific properties
        LineStyle ls = LineStyle(ld, colorOrder);
        setPropertyValue("facecolor", ls.m_lineColor);
        setPropertyValue("displayname", ls.m_legend);
        setPropertyValue("bargroup", getPropertyValue("handle").Scalar());

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++){
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string HggroupBar::getUsingClause(){
        string color = getPropertyValue("facecolor").ColorString();        
        int barstyle = int(getPropertyValue("barstyle").Scalar());
        int f = 3;
        if (barstyle == 1)
            f = 6;
        else if (barstyle == 2)
            f = 7;
        else if (barstyle == 3)
            f = 2;
        else if (barstyle == 4)
            f = 1;
        string legend = getPropertyValue("displayname").StringVal();

        stringstream ss;
        ss << "'-' using 2:xtic(1) lc rgb '" << color << "' fill pattern " << f
            << " title '" << legend << "'";

        Object* p = getParentObject();
        if (p && p->isAxes() && p->getPropertyValue("barlabels").StringVal() == "on" && 
            getPropertyValue("barlayout").StringVal() == "grouped") {

            int numBars = 0;
            int barIdx = -1;
            static_cast<Axes*>(p)->getBarNumberAndIdx(getPropertyValue("handle").Scalar(), numBars, barIdx);

            if (barIdx >= 0) {
                int bargap = int(p->getPropertyValue("bargap").Scalar());
                double u = 1 / double(numBars + bargap);

                int lidx = -numBars / 2 + barIdx;
                double xOffset = 0;
                if (numBars % 2 == 0)
                    xOffset = (lidx + 0.5) * u;
                else
                    xOffset = lidx * u;

                vector<double> y = getYdata();
                double maxVal = 0;
                for (vector<double>::const_iterator it = y.cbegin();
                    it != y.cend(); ++it) {
                    if (maxVal < (*it) * 0.1)
                        maxVal = (*it) * 0.1;
                }

                ss << ",'' using ($0 +" << xOffset << "):($2 + " << maxVal << "):($3) with labels title \"\" font \""
                    << p->getPropertyValue("barlabelsfontname").StringVal() << ","
                    << p->getPropertyValue("barlabelsfontsize").Scalar() << "\"";
            }
        }
        return ss.str();
    }

    string HggroupBar::getLegend(){
        return "";
    }

    string HggroupBar::getWithClause(int lineId){
        return "";
    }

    string HggroupBar::getLineStyle(int line){
        return "";
    }

    void HggroupBar::cleanup(GnuplotOutput *out) {
    }

    void HggroupBar::putData(GnuplotOutput* out)
    {
        Object* p = getParentObject();
        if (!p)
            return;
        // check if x axis has string categories
        Currency catCurrency = p->getProperty("xcategories").getValue().getCurrency();
        bool hasStringCats = catCurrency.IsCellArray();
        HML_CELLARRAY* catCA = nullptr;
        if (catCurrency.IsCellArray() && !catCurrency.CellArray()->IsEmpty())
            catCA = catCurrency.CellArray();

        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> z = getZdata();

        if (z.size() == 0) {     // 2d data

            if (catCA) {
                size_t count = min(catCA->Size(), int(y.size()));
                for (int j = 0; j < count; j++) {
                    Currency cat = (*catCA)(j);
                    if (cat.IsString())
                        out->printf("%s %g\n", cat.StringVal().c_str(), y[j]);
                    else if (cat.IsScalar())
                        out->printf("%g %g\n", cat.Scalar(), y[j]);
                    else
                        out->printf("%g %g\n", x[j], y[j]);
                }
            }
            else {
                size_t count = min(x.size(), y.size());
                for (int j = 0; j < count; j++) {
                    out->printf("%g %g\n", x[j], y[j]);
                }
            }
        }
        else {                // 3d data
            size_t count = min(x.size(), y.size());
            count = min((int)count, (int)z.size());
            for (int j = 0; j < count; j++) {
                out->printf("%g %g %g\n", x[j], y[j], z[j]);
            }
        }
        // end of data
        out->printf("e\n");

        // data for bar labels
        if (p->getPropertyValue("barlabels").StringVal() == "on" &&
            getPropertyValue("barlayout").StringVal() == "grouped") {

            string boldString = p->getPropertyValue("barlabelsfontweight").StringVal() == "bold" ? "/:Bold" : "";
            string italicString = p->getPropertyValue("barlabelsfontangle").StringVal() == "italic" ? "/:Italic" : "";

            size_t count = min(x.size(), y.size());
            for (int j = 0; j < count; j++) {
                out->printf("%g %g \"{%s{%s{%g}}}\"\n", x[j], y[j], boldString.c_str(), italicString.c_str(), y[j]);
            }
            out->printf("e\n");
        }
    }

    bool HggroupBar::setPropertyValue(const string& name, VALUETYPE value)
    {
        bool ret = Object::setPropertyValue(name, value);
        if (name == "barlayout") {
            Object* p = getParentObject();
            if (p && p->isAxes()) {
                static_cast<Axes*>(p)->setBarLayout(value.StringVal());
            }
        }
        else if (name == "barwidth") {
            Object* p = getParentObject();
            if (p && p->isAxes()) {
                static_cast<Axes*>(p)->setBarWidth(value.Scalar());
            }
        }
        return ret;
    }

    Hist::Hist()
        :HggroupBar()
    {
        m_type = ObjectType::HIST;
        setPropertyValue("barwidth", 0.4);
    }

    void Hist::init(const LineData &ld){
        HggroupBar::init(ld);
    }

    HggroupScatter::HggroupScatter()
        :Drawable()
    {
        m_type = ObjectType::SCATTER;

        m_ps.push_back(Property("type", string("hggroup"), PropertyType::STRING) );
        m_ps.push_back(Property("marker", string("o"), PropertyType::STRING));
        m_ps.push_back(Property("markerfacecolor", string("blue"), PropertyType::COLOR));
        m_ps.push_back(Property("markersize", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("dataxoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("dataxscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datazoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datazscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("markerevery", 1.0, PropertyType::DOUBLE));

        // not yet supported properties
        m_ps.push_back(Property("units", string("units"), PropertyType::UNSUPPORTED));
    }

    void HggroupScatter::init(const LineData& ld)
    {
        Drawable::init(ld);

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        // update line-specific properties
        LineStyle ls = LineStyle(ld, colorOrder);
        if (!ls.m_markerStyle.empty())
            setPropertyValue("marker", ls.m_markerStyle);
        setPropertyValue("markerfacecolor", ls.m_markerColor);
        setPropertyValue("markersize", ls.m_markerSize);
        setPropertyValue("displayname", ls.m_legend);

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string HggroupScatter::getUsingClause(){
        int markerEvery = int(getPropertyValue("markerevery").Scalar());
        stringstream ss;
        ss << "'-' using 1:2 every " \
            << markerEvery << " ";
        return ss.str();
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
        int pt = 7;
        if (marker == "s"){
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
        }
            
        string markerColor = getPropertyValue("markerfacecolor").ColorString();
        int markerSize = int(getPropertyValue("markersize").Scalar());

        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line << " dt " << 0
           << " lc rgb '" << markerColor << "'"
           << " pt " << pt << " ps " << markerSize
           << "\n";
        return ss.str();
    }

    void HggroupScatter::cleanup(GnuplotOutput *out) {
    }

    HggroupScatter3::HggroupScatter3()
        :HggroupScatter()
    {
        m_type = ObjectType::SCATTER_3D;
    }
    
    string HggroupScatter3::getUsingClause(){
        int markerEvery = int(getPropertyValue("markerevery").Scalar());
        stringstream ss;
        ss << "'-' using 1:2:3 every " \
            << markerEvery << " ";
        return ss.str();
    }

    Surface::Surface()
        :Drawable(), _minZ(0), _maxZ(1)
    {
        m_type = ObjectType::SURFACE;

        m_ps.push_back(Property("type", string("surface"), PropertyType::STRING) );
        m_ps.push_back(Property("meshlines", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("color", string("blue"), PropertyType::COLOR));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));

		// not yet supported properties
        m_ps.push_back(Property("units", string("units"), PropertyType::UNSUPPORTED));
    }

    void Surface::init(const LineData& ld){
        Drawable::init(ld);

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        // update line-specific properties
        LineStyle ls = LineStyle(ld, colorOrder);
        setPropertyValue("color", ls.m_lineColor);

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string Surface::getUsingClause(){
        return "'-' using 1:2:3 ";
    }

    string Surface::getLegend(){
        return "title ''";
    }

    string Surface::getWithClause(int lineId){
        return "with pm3d ";
    }

    string Surface::getLineStyle(int line){
        stringstream ss;
                
        if (getPropertyValue("meshlines").StringVal() == "on") {
            ss << "unset style line " << line << "\n";
            ss << "set style line " << line << " lc \"black\" \n";
            ss << "set pm3d border linestyle " << line << "\n";
        }
        else {
            ss << "unset style line " << line << "\n";
            ss << "set pm3d\n";
        }
        ss << "set view 60, 315,1,1\n"
           << "unset contour\n"
           << "unset view\n";

        return ss.str();
    }

    void Surface::cleanup(GnuplotOutput *out) {
    }

    void Surface::getMinMaxZ(double& min, double& max)
    {
        min = _minZ;
        max = _maxZ;
    }

    void Surface::putData(GnuplotOutput *out){
        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> z = getZdata();


        _minZ = 0;
        _maxZ = 1;
        if (!z.empty())
        {
            _minZ = z[0];
            _maxZ = z[0];
        }

        if (_xcolcount > 1 && _ycolcount > 1){
            size_t xsize = x.size();
            size_t ysize = y.size();
            size_t zsize = z.size();
            if (! ((xsize == ysize) &&
                   (xsize == zsize)) ){
                throw;
            }

            int rowcount = (int)zsize / _zcolcount;
            for (int i = 0; i < zsize; i++){
                double zVal = z[i];
                out->printf("%g %g %g\n",
                            x[i], y[i], z[i]);
                if ((i + 1) % rowcount == 0){
                    out->printf("\n");
                }
                // find min/max for the colorbar
                if (zVal < _minZ)
                    _minZ = zVal;
                if (zVal > _maxZ)
                    _maxZ = zVal;
            }
        } else {
            size_t col = y.size();
            for (int xsub = 0; xsub < x.size(); xsub++){
                for (int ysub = 0; ysub < col; ysub++){
                    double zVal = z[xsub * col + ysub];
                    out->printf("%g %g %g\n",
                                x[xsub], y[ysub], zVal);
                    // find min/max for the colorbar
                    if (zVal < _minZ)
                        _minZ = zVal;
                    if (zVal > _maxZ)
                        _maxZ = zVal;
                }
                out->printf("\n");
            }
        }
        // end of data
        out->printf("e\n");
    }

    Mesh::Mesh()
        :Surface()
    {
        m_type = ObjectType::MESH;
    }

    string Mesh::getUsingClause(){
        return "'-' using 1:2:3 ";
    }

    string Mesh::getWithClause(int lineId){
        stringstream ss;
        ss << "with lines linestyle " << lineId;
        return ss.str();
    }
    
    string Mesh::getLineStyle(int line){
        string color = getPropertyValue("color").ColorString();
        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line << " lc rgb '" << color << "'\n";
        ss << "set pal maxcolor 0\n"
           << "unset pm3d\n"
           << "set surface\n"
           << "set isosamples 40, 40; set samples 40, 40\n"
           << "set view 60, 315,1,1\n"
           << "unset contour\n"
           << "unset view\n";
        return ss.str();
    }

    Line3::Line3()
        :Line()
    {
        m_type = ObjectType::LINE_3D;
    }
    
    string Line3::getUsingClause(){
        return "'-' using 1:2:3";
    }

    Contour3::Contour3()
        :Surface()
    {
        m_type = ObjectType::CONTOUR_3D;
    }

    string Contour3::getLineStyle(int line){
        return "set contour surface\n"
            "unset surface\n"
            "unset pm3d\n"
            "set cntrparam levels 10\n"
            "set view 60, 315,1,1\n"
            "unset view\n";
    }

    string Contour3::getWithClause(int)
    {
        return string();
    }

    void Contour3::cleanup(GnuplotOutput *out) {
        out->printf("unset contour\n");
        out->printf("set surface\n");
        //out->printf("unset dgrid3d\n");
    }

    string Contour3::getUsingClause(){
        return "'-' using 1:2:3 with lines nosurf";
    }

    Contour::Contour()
        :Contour3()
    {
        m_type = ObjectType::CONTOUR;
    }
    
    string Contour::getLineStyle(int line){
        return "set contour base\n"
            "unset surface\n"
            "unset pm3d\n"
            "set view map\n"
            "set cntrparam levels 10\n";
    }

    void Contour::cleanup(GnuplotOutput *out) {
        out->printf("unset contour\n");
        out->printf("unset view\n");
    }

    Stem::Stem()
        :Line()
    {
        m_type = ObjectType::STEM;
    }
    
    string Stem::getWithClause(int lineId){
        stringstream ss;
        ss << "with impulses linestyle " << lineId;
        return ss.str();
    }

    Loglog::Loglog()
        :Line()
    {
        m_type = ObjectType::LOGLOG;
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

    Semilogx::Semilogx()
        :Line()
    {
        m_type = ObjectType::SEMILOGX;
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

    Semilogy::Semilogy()
        :Line()
    {
        m_type = ObjectType::SEMILOGY;
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

    Text::Text()
        :Object()
    {
        m_type = ObjectType::TEXT;
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
        m_ps.push_back(Property("horizontalalignment", string("center"), PropertyType::STRING));
        m_ps.push_back(Property("verticalalignment", string("middle"), PropertyType::STRING));
        m_ps.push_back(Property("offset", 1.0, PropertyType::DOUBLE));

        m_objectMap[handle] = this;

		// not yet supported properties
		m_ps.push_back(Property("borderwidth", string("borderwidth"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("units", string("units"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("interpreter", string("interpreter"), PropertyType::UNSUPPORTED));
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
        ss << getModifiedString(text);
        if (fontweight == "bold"){
            ss << "}";
        }
        if (fontangle == "italic"){
            ss << "}";
        }


        string ha = getPropertyValue("horizontalalignment").StringVal();
        string va = getPropertyValue("verticalalignment").StringVal();
        double off = getPropertyValue("offset").Scalar();
        string align = "center";
        double offsetX = 0.0;
        double offsetY = 0.0;
        int pointSize = 0;
        if (ha == "left") {
            offsetX = off;
            align = ha;
            pointSize = 2;
        }
        else if (ha == "right") {
            offsetX = -off;
            align = ha;
            pointSize = 2;
        }
        if (va == "top") {
            offsetY = off;
            pointSize = 2;
        }
        else if (va == "bottom") {
            offsetY = -off;
            pointSize = 2;
        }

        out->printf("set label \"%s\" at %g,%g,%g %s offset %g,%g point ps %d font \"%s,%g\" tc rgb \"%s\" enhanced\n",
                    ss.str().c_str(), x, y, z, align.c_str(), offsetX, offsetY, pointSize,
                    fontname.c_str(), fontsize, color.c_str());
    }

    void Text::update(GnuplotOutput *out, const string& type){
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
        ss << getModifiedString(text);
        if (fontweight == "bold"){
            ss << "}";
        }
        if (fontangle == "italic"){
            ss << "}";
        }

        string lt = type;
        string loc = "";
        if (type == "tlabel") {
            lt = "label";
            loc = "at screen 0.5, 0.97 center";
        }
        else if (type == "rlabel") {
            lt = "label";
            loc = "at screen 0.97, 0.5 center rotate by 90";
        }
        else if (type == "llabel") {
            lt = "label";
            loc = "at screen 0.03, 0.5 center rotate by 90";
        }
        else if (type == "blabel") {
            lt = "label";
            loc = "at screen 0.5, 0.03 center";
        }
        
        out->printf("set %s \"%s\" %s font \"%s,%g\" tc rgb \"%s\" enhanced\n",
                    lt.c_str(), ss.str().c_str(),loc.c_str(),
                    fontname.c_str(), fontsize,
                    color.c_str());
    }

    bool Text::isText()
    {
        return true;
    }

    bool Text::getVisible()
    {
        return getPropertyValue("visible").StringVal() == "on";
    }

    string Text::getModifiedString(const string& text)
    {
        stringstream modifiedStr;
        for (int i = 0; i < text.size(); ++i)
        {
            if (text[i] == '{' || text[i] == '}')
                modifiedStr << "\\\\" << text[i];
            else
                modifiedStr << text[i];
        }
        return modifiedStr.str();
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

    
    SecondaryYAxis::SecondaryYAxis(double parentHandle)
        :Object(), m_ylabel(new Text), m_parentAxesHandle(parentHandle)
    {
        m_type = ObjectType::SECONDARYYAXIS;

        double handle = m_handlePool->allocHandle();
        // all of the properties of axes are added here
        m_ps.push_back(Property("type", string("axes"), PropertyType::STRING));
        m_ps.push_back(Property("handle", handle, PropertyType::DOUBLE));
        vector<double> lim;
        lim.push_back(0); lim.push_back(1);
        m_ps.push_back(Property("ylim", lim, PropertyType::VEC_DOUBLE));
        m_ps.push_back(Property("ygrid", string("off"), PropertyType::STRING));
        m_ps.push_back(Property("ylabel", Currency(m_ylabel.get(), "text"), PropertyType::POINTER));
        m_ps.push_back(Property("yminorgrid", string("off"), PropertyType::STRING));
        m_ps.push_back(Property("yscale", string("linear"), PropertyType::STRING));
        m_ps.push_back(Property("ycolor", Color(string("black")), PropertyType::COLOR));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("fontangle", string("regular"), PropertyType::STRING));
        m_ps.push_back(Property("fontname", string("arial"), PropertyType::STRING));
        m_ps.push_back(Property("fontsize", 7.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("fontweight", string("normal"), PropertyType::STRING));
        m_ps.push_back(Property("ynumericformat", string("auto"), PropertyType::STRING));
        m_ps.push_back(Property("ynumericprecision", 5.0, PropertyType::DOUBLE));
        
        // not yet supported properties
        m_ps.push_back(Property("ytick", string("ytick"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("ytickmethod", string("ytickmethod"), PropertyType::UNSUPPORTED));

        m_objectMap[handle] = this;
    }

    void SecondaryYAxis::update(GnuplotOutput* out)
    {
        if (getPropertyValue("visible").StringVal() == "off")
            return;

        string tickFontName = getPropertyValue("fontname").StringVal();
        double tickFontSize = getPropertyValue("fontsize").Scalar();
        string tickFontWeight = getPropertyValue("fontweight").StringVal();
        string tickFontAngle = getPropertyValue("fontangle").StringVal();
        string color = getPropertyValue("ycolor").ColorString();
        string nf = getPropertyValue("ynumericformat").StringVal();
        int np = int(getPropertyValue("ynumericprecision").Scalar());

        string boldString = tickFontWeight == "bold" ? "/:Bold" : "";
        string italicString = tickFontAngle == "italic" ? "/:Italic" : "";
        string f = "g";
        if (nf == "scientific")
            f = "e";
        else if (nf == "fixed")
            f = "f";
        
        char numF[24];
        sprintf(numF, "%s.%d%s", "%", np, f.c_str());
        out->printf("set y2tics format \"{%s {%s %s}}\" textcolor rgb '%s' font \"%s,%g\" \n", boldString.c_str(),
            italicString.c_str(), numF, color.c_str(), tickFontName.c_str(), tickFontSize);

        string yscale = getPropertyValue("yscale").StringVal();
        if (yscale == "log") {
            out->printf("set logscale y2\n");
        }
        else {
            out->printf("unset logscale y2\n");
        }

        string ygrid = getPropertyValue("ygrid").StringVal();
        if (ygrid == "on") {
            out->printf("set grid y2tics\n");
        }
        else {
            out->printf("set grid noy2tics\n");
        }
        string ymgrid = getPropertyValue("yminorgrid").StringVal();
        if (ymgrid == "on") {
            out->printf("set my2tics 5\nset grid my2tics\n");
        }
        else {
            out->printf("set my2tics 1\n");
        }

        m_ylabel->update(out, "y2label");
    }

    bool SecondaryYAxis::isAxes()
    {
        return true;
    }

    double SecondaryYAxis::getParentAxesHandle() const
    {
        return m_parentAxesHandle;
    }
}
