/**
* @file Object.cxx
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

#include "Object.h"
#include "LineStyle.h"
#include "OML_Error.h"
#include <cctype>
#include <algorithm>
#include <cfloat>
#include <iostream>
#include <cassert>
#include <limits>
#include "CoreMain.h"
#include "GnuplotOutput.h"
#include "hwMatrixN_NMKL.h"

using namespace std;

namespace omlplot{

    typedef Property::ValueType VALUETYPE;

    bool isSameDouble(const double& a, const double& b){
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
        m_ps.push_back(Property("handlevisibility", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("tag", string(""), PropertyType::STRING));
    }

    Object::~Object(){
        double parent = getPropertyValue("parent").Scalar();
        Object *p = getObject(parent);
        p->removeChild(this);

        vector<double> children = getAllChildren();
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
        if (name == "children") {
            Object* root = getObject(0.0);
            bool showHidden = root->getPropertyValue("showhiddenhandles").StringVal() == "on";
            if (showHidden)
                return p.getValue();

            std::vector<double> ch = p.getValue().Vector();
            std::vector<double> nonhidden;
            std::vector<double>::const_iterator it = ch.cbegin();
            for (; it != ch.cend(); ++it) {
                Object* obj = getObject(*it);
                if (obj && obj->getPropertyValue("handlevisibility").StringVal() == "on") {
                    nonhidden.push_back(*it);
                }
            }
            return VALUETYPE(nonhidden);
        }
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
        if (p.getType() == PropertyType::POINTER)
            return true;

        if (!validatePropertyType(name, value))
            return true;

        if (name == "color"     || name == "xcolor"         ||
            name == "ycolor"    || name == "zcolor"         ||
            name == "facecolor" || name == "markerfacecolor"||
            name == "edgecolor" || name == "bordercolor"    ||
            name == "gridcolor" || name == "framecolor"     ||
            name == "zerolinecolor" || name == "backgroundcolor"){
            if (value.isCurrency()){
                Currency c = value.getCurrency();
                if (getObjectType() == ObjectType::RECTANGLE && name == "facecolor"){
                    if (c.IsString() && c.StringVal() == "none") {
                        static_cast<Rectangle*>(this)->SetFillNone(true);
                        return true;
                    }
                    else {
                        static_cast<Rectangle*>(this)->SetFillNone(false);
                    }
                }
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

    bool Object::isPropertySupported(const string& name) {
        Property& p = getProperty(name);
        return p.getType() != PropertyType::UNSUPPORTED;
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
        Property& pp = getProperty("parent");
        pp.setValue(parent);
    }

    void Object::addChild(Object *o){
        vector<double> children = getAllChildren();
        children.push_back(o->getHandle());
        Property& ch = getProperty("children");
        ch.setValue(children);

        Property& pp = o->getProperty("parent");
        pp.setValue(this->getHandle());
    }

    void Object::removeChild(Object *o){
        vector<double> children = getAllChildren();
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
        Property& chP = getProperty("children");
        chP.setValue(children);
    }

    vector<double> Object::getAllChildren()
    {
        Property& p = getProperty("children");
        return p.getValue().Vector();
    }

    void Object::repaint(){
        Currency parent = getPropertyValue("parent").getCurrency(); 
        if (parent.IsVector() && !parent.Vector().empty())
        {
            vector<double> p = parent.Vector();
            Object* o = getObject(p[0]);
            o->repaint();
        }
    }

    ObjectType Object::getObjectType() const
    {
        return m_type;
    }

    bool Object::validatePropertyType(const string& name, VALUETYPE value)
    {
        // colors
        if (name == "color"     || name == "xcolor"         || name == "ycolor"         || 
            name == "zcolor"    || name == "facecolor"      || name == "markerfacecolor"||
            name == "edgecolor" || name == "bordercolor"    || name == "gridcolor"      || 
            name == "framecolor"|| name == "zerolinecolor"  || name == "backgroundcolor")
        {
            //  colors are also validated in the Color constructor
            if (!value.isCurrency() && !value.isColor())
                throw OML_Error(OML_ERR_PLOT_INVALID_COLOR);
            return true;
        }

        // read-only properties
        if (name == "children"      || name == "handle" || name == "parent"     ||
            name == "type"          || name == "legend" || name == "colorbar"   ||
            name == "secondaryyaxis"|| name == "areagroup" || name == "bargroup")
        {
            CoreMain::getInstance()->AddWarningString("cannot edit [" + name + "]");
            return false;
        }

        // all other must be Currency
        if (!value.isCurrency())
            throw OML_Error(OML_ERR_OPTIONVAL);

        Currency cur = value.getCurrency();
        if (name == "handlevisibility"  || name == "showhiddenhandles"  || name == "visible"    ||
            name == "xgrid"             || name == "ygrid"              || name == "zgrid"      ||
            name == "xminorgrid"        || name == "yminorgrid"         || name == "zminorgrid" ||
            name == "barlabels"         || name == "polartiptotail"     || name == "tiptotail"  ||
            name == "meshlines"         || name == "showarrowhead"      || name == "autoscale")
        {
            // 'on'/'off' options
            if (cur.IsString())
            {
                string val = cur.StringVal();
                if (val != "on" && val != "off")
                {
                    throw OML_Error(OML_Error(OML_ERR_FUNCSWITCH).GetErrorMessage() + " to set [" + name + "]");
                }
            }
            else
            {
                throw OML_Error(OML_Error(OML_ERR_FUNCSWITCH).GetErrorMessage() + " to set [" + name + "]");
            }
        }
        else if (name == "tag"                  || name == "units"              || name == "toplabel"           || 
                 name == "bottomlabel"          || name == "rightlabel"         || name == "leftlabel"          ||
                 name == "title"                || name == "xlabel"             || name == "ylabel"             ||
                 name == "zlabel"               || name == "xscale"             || name == "yscale"             ||
                 name == "zscale"               || name == "nextplot"           || name == "fontangle"          ||
                 name == "fontname"             || name == "fontweight"         || name == "xnumericformat"     ||
                 name == "ynumericformat"       || name == "znumericformat"     || name == "barlabelsfontname"  ||
                 name == "barlabelsfontweight"  || name == "barlabelsfontangle" || name == "colorbarscale"      ||
                 name == "polarmethod"          || name == "xaxislocation"      || name == "yaxislocation"      ||
                 name == "displayname"          || name == "linestyle"          || name == "marker"             ||
                 name == "barlayout"            || name == "string"             || name == "horizontalalignment"||
                 name == "verticalalignment")
        {
            // string options
            if (!cur.IsString())
                throw OML_Error(OML_Error(OML_ERR_STRING).GetErrorMessage() + " to set [" + name + "]");
            
            validateStringPropertyValue(name, cur);
        }
        else if (name == "position" || name == "xlim" || name == "ylim" || name == "zlim" || name == "clim")
        {
            // vector
            if (!cur.IsVector())
                throw OML_Error(OML_Error(OML_ERR_VECTOR).GetErrorMessage() + " to set [" + name + "]");

        }
        else if (name == "colormap" || name == "colororder")
        {
            // matrix
            if (!cur.IsMatrix() || cur.Matrix()->N() != 3)
                throw OML_Error("must be an Mx3 matrix to set [" + name + "]");

        }
        else if (name == "fontsize"         || name == "xnumericprecision"  || name == "ynumericprecision"  ||
                 name == "znumericprecision"|| name == "barlabelsfontsize"  || name == "bargap"             ||
                 name == "xminortick"       || name == "yminortick"         || name == "zminortick"         ||
                 name == "linewidth"        || name == "markersize"         || name == "dataxoffset"        ||
                 name == "datayoffset"      || name == "datazoffset"        || name == "dataxscale"         ||
                 name == "datayscale"       || name == "datazscale"         || name == "markerevery"        ||
                 name == "basevalue"        || name == "barwidth"           || name == "barstyle"           ||
                 name == "x"                || name == "y"                  || name == "z"                  ||
                 name == "offset"           || name == "borderwidth"        || name == "autoscalefactor"    ||
                 name == "maxheadsize"      || name == "currentaxes"        || name == "currentfigure"      ||
                 name == "markerfacealpha")
        {
            // scalar
            if (!cur.IsScalar())
                throw OML_Error(OML_Error(OML_ERR_SCALAR).GetErrorMessage() + " to set [" + name + "]");
        }
        else if (name == "xcategories" || name == "ycategories" || name == "xticklabel" || name == "yticklabel")
        {
            // cell
            if (!cur.IsCellArray())
                throw OML_Error(OML_Error(OML_ERR_CELL).GetErrorMessage() + " to set [" + name + "]");
        }
        else if (name == "colorlevels" || name == "xtick" || name == "ytick" || name == "ztick")
        {
            // scalar or vector
            if (!cur.IsScalar() && !cur.IsVector())
                throw OML_Error(OML_Error(OML_ERR_SCALARVECTOR).GetErrorMessage() + " to set [" + name + "]");
        }
        else if (name == "xdata" || name == "ydata" || name == "zdata" || name == "udata" || name == "vdata")
        {
            // scalar or matrix
            if (!cur.IsScalar() && !cur.IsVector())
                throw OML_Error(OML_Error(OML_ERR_SCALARMATRIX).GetErrorMessage() + " to set [" + name + "]");
        }
        else if (name == "location")
        {
            // legend location, string or vector
            // scalar or matrix
            if (!cur.IsString() && !cur.IsVector())
                throw OML_Error("must be a string or vector to set [" + name + "]");
        }
        return true;
    }

    void Object::validateStringPropertyValue(const string& name, const Currency& value)
    {
        if (!value.IsString())
            throw OML_Error(OML_Error(OML_ERR_STRING).GetErrorMessage() + " to set [" + name + "]");

        string val = value.StringVal();
        string commonError = OML_Error(OML_ERR_OPTION).GetErrorMessage() + " for [" + name + "]";
        if (name == "xscale" || name == "yscale" || name == "zscale" || name == "colorbarscale")
        {
            if (val != "linear" && val != "log")
                throw OML_Error(commonError+"; must be 'linear' or 'log'");
        }
        else if (name == "fontangle" || name == "barlabelsfontangle")
        {
            if (val != "regular" && val != "italic")
                throw OML_Error(commonError+"; must be 'italic' or 'regular'");
        }
        else if (name == "fontweight" || name == "barlabelsfontweight")
        {
            if (val != "normal" && val != "bold")
                throw OML_Error(commonError+"; must be 'bold' or 'normal'");
        }
        else if (name == "xnumericformat" || name == "ynumericformat" || name == "znumericformat")
        {
            if (val != "auto" && val != "scientific" && val != "fixed")
                throw OML_Error(commonError+"; must be 'auto', 'fixed' or 'scientific'");
        }
        else if (name == "polarmethod")
        {
            if (val != "phasevsmag" && val != "vector" && val != "radar")
                throw OML_Error(commonError+"; must be 'phasevsmag', 'radar' or 'vector'");
        }
        else if (name == "xaxislocation")
        {
            if (val != "bottom" && val != "top")
                throw OML_Error(commonError+"; must be 'bottom' or 'top'");
        }
        else if (name == "yaxislocation")
        {
            if (val != "left" && val != "right")
                throw OML_Error(commonError+"; must be 'left' or 'right'");
        }
        else if (name == "linestyle")
        {
            if (!LineStyle::isAbbrLine(*val.c_str()))
                throw OML_Error(commonError + "; must be ' ', '-', '--', ':', '-.' or '-:'");
        }
        else if (name == "marker")
        {
            if (!val.empty() && !LineStyle::isAbbrMarker(*val.c_str()))
                throw OML_Error(commonError + "; must be ' ', 's', '^', 'v', 'x', 'o', 'd', '+', '*' or '.'");
        }
        else if (name == "barlayout")
        {
            if (val != "grouped" && val != "stacked")
                throw OML_Error(commonError+"; must be 'grouped' or 'stacked'");
        }
        else if (name == "horizontalalignment")
        {
            if (val != "left" && val != "right" && val != "center")
                throw OML_Error(commonError+"; must be 'left', 'center' or 'right'");
        }
        else if (name == "verticalalignment")
        {
            if (val != "top" && val != "bottom" && val != "middle")
                throw OML_Error(commonError+"; must be 'top', 'middle' or 'bottom'");
        }
        else if (name == "nextplot")
        {
            if (val != "replace" && val != "add")
                throw OML_Error(commonError+"; must be 'add' or 'replace'");
        }
        else if (name == "units")
        {
            if (val != "pixels" && val != "normalized")
                throw OML_Error(commonError+"; must be 'normalized' or 'pixels'");
        }
    }

    void Object::dump(Object *o, int level){
        for (int i = 0; i < level; i++){
            cout << '\t';
        }
        cout << o->getType() << " ( " << o->getHandle() << " )" << endl;

        vector<double> children = o->getAllChildren();
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
        m_ps.push_back(Property("showhiddenhandles", string("off"), PropertyType::STRING));
        m_ps.push_back(Property("units", string("pixels"), PropertyType::STRING));
        m_objectMap[double(0)] = this;
    }

    void Root::repaint(){
        update(nullptr);
    }

    void Root::update(GnuplotOutput *out){
        vector<double> children = getAllChildren();
        vector<double>::iterator it = children.begin();
        for ( ; it != children.end(); ++it){
            Object* o = getObject(*it);
            if (o && o->isFigure())
                static_cast<Figure*>(o)->update();
        }
    }

    Figure::Figure()
        :Object(), m_out(new GnuplotOutput), m_tLabel(new Text(false)), 
        m_rLabel(new Text(false)), m_lLabel(new Text(false)), 
        m_bLabel(new Text(false)), m_rows(1), m_cols(1), m_isGridLayout(false)
    {
        m_type = ObjectType::FIGURE;

        double handle = m_figureHandlePool->allocHandle();
        init(handle);
    }

    Figure::Figure(int fh)
        :Object(), m_out(new GnuplotOutput), m_tLabel(new Text(false)), 
        m_rLabel(new Text(false)), m_lLabel(new Text(false)), 
        m_bLabel(new Text(false)), m_rows(1), m_cols(1), m_isGridLayout(false)
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
        m_ps.push_back(Property("color", Color(string("white")), PropertyType::COLOR));
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
        // figure background color
        string color = getPropertyValue("color").ColorString();
        m_out->printf("set obj 1000 rectangle from screen 0,0 to screen 1,1 behind fc rgb '%s' fs noborder\n", color.c_str());

        m_out->printf("set multiplot\n");
        
        m_tLabel->update(m_out.get(), "tlabel");
        m_rLabel->update(m_out.get(), "rlabel");
        m_lLabel->update(m_out.get(), "llabel");
        m_bLabel->update(m_out.get(), "blabel");

        vector<double> children = getAllChildren();
        Object *o = nullptr;
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it){
            double child = *it;
            o = getObject(child);
            o->update(m_out.get());
            m_out->printf("unset object 1000\n");
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
        vector<double> children = getAllChildren();
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
        Property& chP = getProperty("children");
        chP.setValue(children);
        Property& caP = getProperty("currentaxes");
        caP.setValue(0);
        
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
    }

    void Figure::saveas(string _filename, string fmt, int width, int height){
        string filename = _filename;
        string extension;
        size_t extindex = filename.rfind(".");
        if (extindex != std::string::npos && (extindex+1) !=std::string::npos)
            extension = filename.substr(extindex+1, filename.length());
        
        if (fmt.empty())
        {
            if (extension.empty()) {
                fmt = "png";
                filename += "." + fmt;
            }
            else {
                fmt = extension;
            }
        }
        else {
            filename += "." + fmt;
        }

        std::transform(fmt.begin(), fmt.end(), fmt.begin(), ::tolower);

        if (fmt == "png" || fmt == "jpeg" || fmt == "svg") {
            m_out->printf("set terminal push\n");
            if (width > 0 && height > 0)
                m_out->printf("set terminal %s size %d, %d\n", fmt.c_str(), width, height);
            else
                m_out->printf("set terminal %s\n", fmt.c_str());
            m_out->printf("set output \'%s\'\n", filename.c_str());

            m_out->printf("replot\n");
            m_out->printf("set terminal pop\n");
            m_out->printf("set output\n");
            m_out->flush();
        }
        else {
            throw OML_Error(OML_ERR_PLOT_UNSUPPORTED_FORMAT);
        }
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
        :Object(), m_title(new Text(false)), m_xlabel(new Text(false)),
        m_ylabel(new Text(false)), m_zlabel(new Text(false)), _borderOn(true), 
        m_barlayout("grouped"), m_barWidth(0.8), m_secYAxisVisible(false), 
        m_updateXAxisRange(false), m_updateYAxisRange(false), m_updateY2AxisRange(false), 
        m_updateZAxisRange(false), m_gridIndex(1), m_tailX(0.0), m_tailY(0.0), 
        m_tailTheta(0.0), m_tailR(0.0), m_rotX(60), m_rotZ(315), m_legend(new Legend()), 
        m_colorbar(new Colorbar()), _axesOn(true), m_xdateTicks(false), m_ydateTicks(false)
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
        m_ps.push_back(Property("legend", Currency(m_legend.get(), "legend"), PropertyType::POINTER));
        m_ps.push_back(Property("colorbar", Currency(m_colorbar.get(), "colorbar"), PropertyType::POINTER));
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
        m_ps.push_back(Property("ycategories", string(""), PropertyType::CELL));
        m_ps.push_back(Property("colororder", Currency(), PropertyType::MATRIX));
        m_ps.push_back(Property("colorbarscale", string("linear"), PropertyType::STRING));
        // default colormap 
        hwMatrix* map = new hwMatrix(4, 3, hwMatrix::DataType::REAL);
        (*map)(0, 0) = 0; (*map)(0, 1) = 0; (*map)(0, 2) = 255;
        (*map)(1, 0) = 0; (*map)(1, 1) = 255; (*map)(1, 2) = 0;
        (*map)(2, 0) = 255; (*map)(2, 1) = 255; (*map)(2, 2) = 0;
        (*map)(3, 0) = 255; (*map)(3, 1) = 0; (*map)(3, 2) = 0;
        m_ps.push_back(Property("colormap", Currency(map), PropertyType::MATRIX));
        m_ps.push_back(Property("colorlevels", Currency(), PropertyType::VEC_DOUBLE));
        m_ps.push_back(Property("secondaryyaxis", Currency(m_secYAxis.get(),"yaxis"), PropertyType::POINTER));
        m_ps.push_back(Property("polarmethod", string("phasevsmag"), PropertyType::STRING));
        m_ps.push_back(Property("polartiptotail", string("off"), PropertyType::STRING));
        m_ps.push_back(Property("xminortick", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("yminortick", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("zminortick", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("xtick", 0.0, PropertyType::VEC_DOUBLE));
        m_ps.push_back(Property("ytick", 0.0, PropertyType::VEC_DOUBLE));
        m_ps.push_back(Property("ztick", 0.0, PropertyType::VEC_DOUBLE));
        m_ps.push_back(Property("xticklabel", string(""), PropertyType::CELL));
        m_ps.push_back(Property("yticklabel", string(""), PropertyType::CELL));
        m_ps.push_back(Property("xaxislocation", string("bottom"), PropertyType::STRING));
        m_ps.push_back(Property("yaxislocation", string("left"), PropertyType::STRING));
        vector<int> gridColor = { 185, 185, 185 };
        m_ps.push_back(Property("gridcolor", Color(gridColor), PropertyType::COLOR));
        m_ps.push_back(Property("framecolor", Color(string("white")), PropertyType::COLOR));
        m_ps.push_back(Property("zerolinecolor", Color(string("white")), PropertyType::COLOR));

        m_objectMap[handle] = this;
        m_title->setParent(this);
        m_legend->setParent(this);
        m_colorbar->setParent(this);
        m_secYAxis->setParent(this);
        m_xlabel->setParent(this);
        m_ylabel->setParent(this);
        m_zlabel->setParent(this);

		// not yet supported properties
		m_ps.push_back(Property("mouseclickcallback", string("mouseclickcallback"), PropertyType::UNSUPPORTED));
		m_ps.push_back(Property("contourtype", string("contourtype"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("barorientation", string("barorientation"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xtickmethod", string("xtickmethod"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("ytickmethod", string("ytickmethod"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("ztickmethod", string("ztickmethod"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xdb10reference", string("xdb10reference"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("ydb10reference", string("ydb10reference"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xdb20reference", string("xdb20reference"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("ydb20reference", string("ydb20reference"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("xweighting", string("xweighting"), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("yweighting", string("yweighting"), PropertyType::UNSUPPORTED));
    }

    Axes::~Axes()
    {
        Object* parent = getParentObject();
        double h = parent->getPropertyValue("currentaxes").Scalar();
        if (isSameDouble(h, getHandle()))
        {
            Property& p = parent->getProperty("currentaxes");
            p.setValue(0);
        }
    }

    void Axes::clear(){
        vector<double> children = getAllChildren();
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
        Property& p = getProperty("children");
        p.setValue(children);

        m_title->setPropertyValue("string", string(""));

        _borderOn = true;
        _axesOn = true;
        m_colorbar->setVisible(true);

        m_updateXAxisRange = false;
        m_updateYAxisRange = false;
        m_updateY2AxisRange = false;
        m_updateZAxisRange = false;
        m_secYAxisVisible = false;

        m_axisOption = "";

        m_xdateTicks = false;
        m_ydateTicks = false;
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

        double x = 0, y = 0, w = 1, h = 1;
        getAxesPosition(x, y, w, h);

        bool isPolarPlot = false;
        bool isBarPlot = false;
        bool isContourPlot = false;
        bool isPColorPlot = false;
        bool isBar3Plot = false;
        vector<double> children = getAllChildren();
        vector<double>::const_iterator it = children.cbegin();
        for (; it != children.cend(); ++it) {
            Object* obj = getObject(*it);
            if (!obj)
                continue;

            if (obj->getObjectType() == ObjectType::POLAR_LINE) {
                isPolarPlot = true;
                break;
            }
            else if (obj->getObjectType() == ObjectType::BAR) {
                isBarPlot = true;
                break;
            }
            else if (obj->getObjectType() == ObjectType::CONTOUR) {
                isContourPlot = true;
                break;
            }
            else if (obj->getObjectType() == ObjectType::PCOLOR) {
                isPColorPlot = true;
                break;
            }
            else if (obj->getObjectType() == ObjectType::HGGROUPBAR3) {
                isBar3Plot = true;
                break;
            }
        }

        // unset any arrows (xline/yline)
        out->printf("unset arrow\n");
        // unset axes rectangle object
        out->printf("unset object 1001\n");

        // reset grids and tics
        out->printf("set border 31\n");
        out->printf("unset polar\nset size nosquare\nset grid nopolar\n");
        out->printf("set grid noxtics\nset grid noytics\nset grid noztics\n");
        out->printf("unset xtics\nunset mxtics\n");
        out->printf("unset x2tics\nunset mx2tics\n");
        out->printf("unset ytics\nunset mytics\n");
        out->printf("unset y2tics\nunset my2tics\n");
        out->printf("unset ztics\nunset mztics\n");
        out->printf("unset xlabel\nunset ylabel\n");
        out->printf("unset x2label\n");
        out->printf("unset y2label\nunset zlabel\n");
        out->printf("unset xzeroaxis\nunset x2zeroaxis\n");
        out->printf("unset yzeroaxis\nunset y2zeroaxis\n");
        out->printf("unset zzeroaxis\n");
        out->printf("unset colorbox\n");


        // reset all ranges - if this is not done here then previous range will remain even if
        // the plot is cleared
        out->printf("unset xrange\n");
        out->printf("unset yrange\n");
        out->printf("unset y2range\n");
        out->printf("unset zrange\n");
        out->printf("unset logscale x\n");
        out->printf("unset logscale y\n");
        out->printf("unset logscale y2\n");
        out->printf("unset logscale z\n");

        string framecolor = getPropertyValue("framecolor").ColorString();
        out->printf("set obj 1001 rectangle from screen %g,%g to screen %g,%g behind fc rgb '%s' fs border rgb '%s'\n",
            x, y, x + w, y + h, framecolor.c_str(), framecolor.c_str());
        string bgcolor = getPropertyValue("color").ColorString();
        out->printf("set obj rectangle from graph 0,0 to graph 1,1 behind fc rgb '%s' fs border rgb '%s'\n", bgcolor.c_str(), bgcolor.c_str());

        bool topXAxis = !is3DPlot() && !isPolarPlot && !isContourPlot && 
            !isPColorPlot && getPropertyValue("xaxislocation").StringVal() == "top";
        bool rightYAxis = !m_secYAxisVisible && !is3DPlot() && !isPolarPlot && !isContourPlot &&
            !isPColorPlot && getPropertyValue("yaxislocation").StringVal() == "right";

        // set the range only for the axes that a range has been specified
        if (m_updateXAxisRange && !isPolarPlot) {
            vector<double> xrange = getPropertyValue("xlim").Vector();
            out->printf("set %s [%g:%g]\n", topXAxis ? "x2range" : "xrange", xrange[0], xrange[1]);
        }
        if (m_updateYAxisRange && !isPolarPlot) {
            vector<double> yrange = getPropertyValue("ylim").Vector();
            out->printf("set %s [%g:%g]\n", rightYAxis ? "y2range" : "yrange", yrange[0], yrange[1]);
        }
        if (m_updateZAxisRange) {
            vector<double> zrange = getPropertyValue("zlim").Vector();
            out->printf("set zrange [%g:%g]\n", zrange[0], zrange[1]);
        }
        if (m_updateY2AxisRange) {
            vector<double> y2range = m_secYAxis->getPropertyValue("ylim").Vector();
            out->printf("set y2range [%g:%g]\n", y2range[0], y2range[1]);
        }

        if (_axesOn) {
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

            string xtag(topXAxis ? "x2tics" : "xtics");

            out->printf("set %s format \"{%s {%s %s}}\" textcolor rgb '%s' font \"%s,%g\" \n", 
                xtag.c_str(), boldString.c_str(), italicString.c_str(), numF, xcolor.c_str(), 
                tickFontName.c_str(), tickFontSize);
            int xtick = (int)getPropertyValue("xminortick").Scalar();
            if (xtick > 0)
                out->printf("set m%s %d\n", xtag.c_str(), xtick);

        // override xticks if it is a 3d bar plot and 'xcategories' are set
        Currency xcat = getPropertyValue("xcategories").getCurrency();
        if (isBar3Plot && xcat.IsCellArray() && xcat.CellArray()->Size() > 0)
        {
            // find all x values from the hggroupbar3 objects
            std::vector<double> xValues;
            it = children.cbegin();
            for (; it != children.cend(); ++it) {
                Object* obj = getObject(*it);
                if (!obj)
                    continue;

                if (obj->getObjectType() == ObjectType::HGGROUPBAR3) {
                    Currency tmpX = obj->getPropertyValue("xdata").getCurrency();
                    if (!tmpX.IsVector())
                        continue;
                 
                    std::vector<double> tmp = tmpX.Vector();
                    std::vector<double>::const_iterator xit = tmp.cbegin();
                    for (; xit != tmp.cend(); ++xit)
                    {
                        if (std::find(xValues.begin(), xValues.end(), (*xit)) == xValues.end())
                            xValues.push_back(*xit);
                    }
                }
            }

            std::sort(xValues.begin(), xValues.end());

            HML_CELLARRAY* ca = xcat.CellArray();
            out->printf("set %s (", xtag.c_str());
            for (int i = 0; i < ca->Size() && i<static_cast<int>(xValues.size()); ++i) {
                string lbl = "";
                if ((*ca)(i).IsString())
                    lbl = (*ca)(i).StringVal();
                out->printf("\"%s\" %g, ", lbl.c_str(), xValues[i]);
            }
            out->printf(")\n");
        }
        else
        {
            Currency tmp = getPropertyValue("xtick").getCurrency();
            if (tmp.IsScalar()) {
                if (!IsZero(tmp.Scalar()))
                    out->printf("set %s %g\n", xtag.c_str(), tmp.Scalar());
            }
            else {
                vector<double> xticks = tmp.Vector();
                Currency tickLbl = getPropertyValue("xticklabel").getCurrency();
                if (tickLbl.IsCellArray() && tickLbl.CellArray()->Size() > 0) {
                    HML_CELLARRAY* ca = tickLbl.CellArray();
                    int S = min((int)xticks.size(), ca->Size());
                    out->printf("set %s (", xtag.c_str());
                    for (int i = 0; i < S; ++i) {
                        string lbl = "";
                        if ((*ca)(i).IsString())
                            lbl = (*ca)(i).StringVal();
                        out->printf("\"%s\" %g, ", lbl.c_str(), xticks[i]);
                    }
                    out->printf(")\n");
                }
                else {
                    vector<double>::const_iterator it = xticks.cbegin();
                    out->printf("set %s (", xtag.c_str());
                    for (; it != xticks.cend(); ++it) {
                        out->printf("%g,", *it);
                    }
                    out->printf(")\n");
                }
            }
        }

        string ycolor = getPropertyValue("ycolor").ColorString();
        string ynf = getPropertyValue("ynumericformat").StringVal();
        int ynp = int(getPropertyValue("ynumericprecision").Scalar());
        f = "g";
        if (ynf == "scientific")
            f = "e";
        else if (ynf == "fixed")
            f = "f";
        sprintf(numF, "%s.%d%s", "%", ynp, f.c_str());

            string ytag(rightYAxis ? "y2tics" : "ytics");
            out->printf("set %s format \"{%s {%s %s}}\" textcolor rgb '%s' font \"%s,%g\" \n", ytag.c_str(), boldString.c_str(),
            italicString.c_str(), numF, ycolor.c_str(), tickFontName.c_str(), tickFontSize);
            int ytick = (int)getPropertyValue("yminortick").Scalar();
            if (ytick > 0)
                out->printf("set m%s %d\n", ytag.c_str(), ytick);

        // override yticks if it is a 3d bar plot and 'ycategories' are set
        Currency ycat = getPropertyValue("ycategories").getCurrency();
        if (isBar3Plot && ycat.IsCellArray() && ycat.CellArray()->Size() > 0)
        {
            // find all x values from the hggroupbar3 objects
            std::vector<double> yValues;
            it = children.cbegin();
            for (; it != children.cend(); ++it) {
                Object* obj = getObject(*it);
                if (!obj)
                    continue;

                if (obj->getObjectType() == ObjectType::HGGROUPBAR3) {
                    Currency tmpY = obj->getPropertyValue("ydata").getCurrency();
                    if (!tmpY.IsVector())
                        continue;

                    std::vector<double> tmp = tmpY.Vector();
                    std::vector<double>::const_iterator yit = tmp.cbegin();
                    for (; yit != tmp.cend(); ++yit)
                    {
                        if (std::find(yValues.begin(), yValues.end(), (*yit)) == yValues.end())
                            yValues.push_back(*yit);
                    }
                }
            }

            std::sort(yValues.begin(), yValues.end());

            HML_CELLARRAY* ca = ycat.CellArray();
            out->printf("set %s (", ytag.c_str());
            for (int i = 0; i < ca->Size() && i < static_cast<int>(yValues.size()); ++i) {
                string lbl = "";
                if ((*ca)(i).IsString())
                    lbl = (*ca)(i).StringVal();
                out->printf("\"%s\" %g, ", lbl.c_str(), yValues[i]);
            }
            out->printf(")\n");
        }
        else
        {
            Currency tmp = getPropertyValue("ytick").getCurrency();
            if (tmp.IsScalar()) {
                if (!IsZero(tmp.Scalar()))
                    out->printf("set %s %g\n", ytag.c_str(), tmp.Scalar());
            }
            else {
                vector<double> yticks = tmp.Vector();
                Currency tickLbl = getPropertyValue("yticklabel").getCurrency();
                if (tickLbl.IsCellArray() && tickLbl.CellArray()->Size() > 0) {
                    HML_CELLARRAY* ca = tickLbl.CellArray();
                    int S = min((int)yticks.size(), ca->Size());
                    out->printf("set %s (", ytag.c_str());
                    for (int i = 0; i < S; ++i) {
                        string lbl = "";
                        if ((*ca)(i).IsString())
                            lbl = (*ca)(i).StringVal();
                        out->printf("\"%s\" %g, ", lbl.c_str(), yticks[i]);
                    }
                    out->printf(")\n");
                }
                else {
                    vector<double>::const_iterator it = yticks.cbegin();
                    out->printf("set %s (", ytag.c_str());
                    for (; it != yticks.cend(); ++it) {
                        out->printf("%g,", *it);
                    }
                    out->printf(")\n");
                }
            }
        }

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
            int ztick = (int)getPropertyValue("zminortick").Scalar();
            if (ztick > 0)
                out->printf("set mztics %d\n", ztick);

        Currency tmp = getPropertyValue("ztick").getCurrency();
            if (tmp.IsScalar()) {
                if (!IsZero(tmp.Scalar()))
                    out->printf("set ztics %g\n", tmp.Scalar());
        }
        else {
                vector<double> zticks = tmp.Vector();
                vector<double>::const_iterator it = zticks.cbegin();
                out->printf("set ztics (");
                for (; it != zticks.cend(); ++it) {
                    out->printf("%g,", *it);
        }
                out->printf(")\n");
            }
        }
        else {
            out->printf("set xtics format \"\" scale 0\nset ytics format \"\" scale 0\n");
            out->printf("set y2tics format \"\" scale 0\nset ztics format \"\" scale 0\n");
            out->printf("set x2tics format \"\" scale 0\n");
        }

        string xscale = getPropertyValue("xscale").StringVal();
        if (xscale == "log")
            out->printf("set logscale %s\n", topXAxis ? "x2" : "x");
        string yscale = getPropertyValue("yscale").StringVal();
        if (yscale == "log")    
            out->printf("set logscale y\n");
        string zscale = getPropertyValue("zscale").StringVal();
        if (zscale == "log")
            out->printf("set logscale z\n");

        m_title->update(out, "title");
        if (_axesOn) {
            m_xlabel->update(out, topXAxis ? "x2label" : "xlabel");
            m_ylabel->update(out, rightYAxis ? "y2label" : "ylabel");
        m_zlabel->update(out, "zlabel");
        }

        if (m_secYAxisVisible)
            m_secYAxis->update(out);

        out->printf("set origin %f,%f\n", x, y);
        out->printf("set size %f,%f\n", w, h);
            out->printf("unset key\nplot '-' with points ps 0\n0 0\n1 1\ne\n");

        string gridcolor = getPropertyValue("gridcolor").ColorString();
        string gridStyle = "lt 1 lc rgb '" + gridcolor + "', lt 0 lc 0";

        string xgrid = getPropertyValue("xgrid").StringVal();
        if (xgrid == "on")
            out->printf("set grid %s %s\n", topXAxis ? "x2tics" : "xtics", gridStyle.c_str());

        string xmgrid = getPropertyValue("xminorgrid").StringVal();
        if (xmgrid == "on") 
            out->printf("set grid %s\n", topXAxis ? "mx2tics" : "mxtics");

        string ygrid = getPropertyValue("ygrid").StringVal();
        if (ygrid == "on")
            out->printf("set grid %s %s\n", rightYAxis ? "y2tics" : "ytics", gridStyle.c_str());

        string ymgrid = getPropertyValue("yminorgrid").StringVal();
        if (ymgrid == "on") 
            out->printf("set grid %s\n", rightYAxis ? "my2tics" : "mytics");

        string zgrid = getPropertyValue("zgrid").StringVal();
        if (zgrid == "on") 
            out->printf("set grid ztics %s\n", gridStyle.c_str());

        string zmgrid = getPropertyValue("zminorgrid").StringVal();
        if (zmgrid == "on") 
            out->printf("set grid mztics\n");

        if (!is3DPlot() && !isPolarPlot && !isContourPlot && !isPColorPlot) {
            string zerocolor = getPropertyValue("zerolinecolor").ColorString();
            out->printf("set %s lt 1 lc rgb '%s'\n", topXAxis ? "x2zeroaxis" : "xzeroaxis", zerocolor.c_str());
            out->printf("set %s lt 1 lc rgb '%s'\n", rightYAxis ? "y2zeroaxis" : "yzeroaxis", zerocolor.c_str());
            if (m_secYAxisVisible)
                out->printf("set y2zeroaxis lt 1 lc rgb '%s'\n", zerocolor.c_str());
        }
        m_legend->update(out);

        if (_borderOn) {
            out->printf("set border 31\n");
            if (m_secYAxisVisible)
                out->printf("set ytics nomirror\n");
            else if (rightYAxis)
                out->printf("set y2tics mirror\n");
            if (topXAxis)
                out->printf("set x2tics mirror\n");

        }
        else {
            if (topXAxis) {
                if (m_secYAxisVisible)
                    out->printf("set border 14\n");
                else if (rightYAxis)
                    out->printf("set border 12\n");
        else
                    out->printf("set border 6\n");
                out->printf("set x2tics nomirror\n");
                out->printf("set %s nomirror\n", rightYAxis ? "y2tics" : "ytics");
            }
            else {
                if (m_secYAxisVisible)
                    out->printf("set border 11\n");
                else if (rightYAxis)
                    out->printf("set border 9\n");
        else
                    out->printf("set border 3\n");
                out->printf("set xtics nomirror\n");
                out->printf("set %s nomirror\n", rightYAxis ? "y2tics" : "ytics");
            }
        }

        bool legendBold = m_legend->getPropertyValue("fontweight").StringVal() == "bold";
        bool legendItalic = m_legend->getPropertyValue("fontangle").StringVal() == "italic";

        Drawable* pChild = nullptr;
        vector<string> linestyle;
        vector<string> usingclause;
        vector<string> legend;
        vector<string> withclause;
        vector<vector<double> > zdata;
        vector<string> arrowusingclause;
        vector<string> arrowwithclause;
        vector<string> shapeclause;

        int lineId = 1;
        int textBoxIndex = 1;
        it = children.begin();
        for (; it != children.end(); ++it) {
            Object* obj = getObject(*it);
            if (obj->isDrawable()) {

                pChild = static_cast<Drawable*>(obj);
                if (!pChild->getVisible())
                    continue;

                pChild->setXAxisReference(topXAxis ? "x2" : "x1");
                if (!m_secYAxisVisible)
                    pChild->setYAxisReference(rightYAxis ? "y2" : "y1");

                linestyle.push_back(pChild->getLineStyle(lineId));

                if (pChild->getObjectType() == ObjectType::XLINE ||
                    pChild->getObjectType() == ObjectType::YLINE) {
                    arrowusingclause.push_back(pChild->getUsingClause());
                    arrowwithclause.push_back(pChild->getWithClause(lineId));
                }
                else if (pChild->getObjectType() == ObjectType::ELLIPSE ||
                    pChild->getObjectType() == ObjectType::RECTANGLE) {
                    shapeclause.push_back(pChild->getUsingClause());
                }
                else {
                    usingclause.push_back(pChild->getUsingClause());
                    legend.push_back(pChild->getLegend(legendBold, legendItalic));
                    withclause.push_back(pChild->getWithClause(lineId));
                    zdata.push_back(pChild->getZdata());
                }

               
            }
            else if (obj->isText()) {
                Text* pChild = static_cast<Text*>(obj);
                if (pChild->getVisible()) {
                    pChild->setIndex(textBoxIndex++);
                    pChild->update(out);
                }
            }
            ++lineId;
        }

        bool is3D = zdata.size() > 0 && zdata[0].size() > 0;
        // set view
        if (is3D) {
            if (isContourPlot || isPColorPlot)
                out->printf("set view map\n");
            else
                out->printf("set view %g, %g, 1, 1\n", m_rotX, m_rotZ);
        }
        if (!m_axisOption.empty())
            out->printf(m_axisOption);

        if (is3D)
            m_colorbar->update(out);

        if (isBar3Plot)
        {
            out->printf("set pm3d depthorder base\n");
            it = children.cbegin();
            for (; it != children.cend(); ++it) {
                Object* obj = getObject(*it);
                if (!obj)
                    continue;

                if (obj->getObjectType() == ObjectType::HGGROUPBAR3) {
                    Currency tmpY = obj->getPropertyValue("ydata").getCurrency();
                    if (tmpY.IsVector() && tmpY.Vector().size() > 2)
                    {
                        std::vector<double> tmp = tmpY.Vector();
                        out->printf("set boxdepth %g\n", (tmp[1] - tmp[0]) / 2); 
                    }
                    else
                    {
                        out->printf("set boxdepth 0.4\n");
                    }
                    break;
                }
            }
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
        vector<string>::const_iterator sit = linestyle.cbegin();
        for (; sit != linestyle.cend(); ++sit) {
            out->printf(*sit);
        }

        // draw shapes
        stringstream ss;
        for (size_t i = 0; i < shapeclause.size(); ++i) {
            ss << shapeclause[i] << " ";
        }

        // draw xline, yline arrows
        for (size_t i = 0; i < arrowusingclause.size(); ++i) {
            ss << arrowusingclause[i] << " ";
            ss << arrowwithclause[i] << " ";
        }
        out->printf(ss.str().c_str());

        if (lineCount > 0){
            out->printf("clear\n");
            stringstream ss;
            if (is3D){
                ss << "splot ";
            } else {
                ss << "plot ";
            }
            for (size_t i = 0; i < usingclause.size(); ++i){
                ss << usingclause[i] << " ";
                ss << legend[i] << " ";
                ss << withclause[i] << " ";
                if (i < (usingclause.size() - 1) )
                    ss << ", ";
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
                    ObjectType ot = obj->getObjectType();
                    if (ot == ObjectType::XLINE || ot == ObjectType::YLINE ||
                        ot == ObjectType::RECTANGLE || ot == ObjectType::ELLIPSE)
                        continue;
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

        if (!shapeclause.empty() && usingclause.empty()) {
            vector<double>::iterator it = children.begin();
            for (; it != children.end(); ++it) {
                Object* obj = getObject(*it);
                if (obj->isDrawable()) {

                    pChild = static_cast<Drawable*>(obj);
                    if (!pChild->getVisible())
                        continue;

                    if (pChild->getObjectType() == ObjectType::ELLIPSE) {
                        out->printf(static_cast<Ellipse*>(pChild)->getDummyLine());
                    }
                    else if (pChild->getObjectType() == ObjectType::RECTANGLE) {
                        out->printf(static_cast<Rectangle*>(pChild)->getDummyLine());
                    }
                }
            }
        }
        if (!is3D && m_legend->getVisible())
            out->printf("set key opaque\n");

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
        m_colorbar->setVisible(state);
    }

    bool Axes::getColorbarVisible() const
    {
        return m_colorbar->getVisible();
    }

    std::vector<double> Axes::getColorbarRange()
    {
        // else get the values from the curves
        std::vector<double> ret;
        vector<double> children = getAllChildren();
        for (vector<double>::iterator it = children.begin(); it != children.end(); ++it) 
        {
            Object* obj = getObject(*it);
            ObjectType t = obj->getObjectType();
            if (t == ObjectType::SURFACE || t == ObjectType::MESH ||
                t == ObjectType::CONTOUR_3D || t == ObjectType::CONTOUR ||
                t == ObjectType::TRIMESH || t == ObjectType::TRISURF) {
                
                Drawable* pChild = static_cast<Drawable*>(obj);
                double min, max;
                if (pChild->getMinMaxZ(min, max))
                {
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
        vector<double> children = getAllChildren();
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
        vector<double> children = getAllChildren();
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
        vector<double> children = getAllChildren();
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

    Property& Axes::getProperty(const string& name)
    {
        if (name == "colorbarscale" || name == "colormap" ||
            name == "colorlevels") {
            return m_colorbar->getProperty(name);
        }
        return Object::getProperty(name);
    }

    bool Axes::setPropertyValue(const string& name, VALUETYPE value)
    {
        if (name == "colorbarscale" || name == "colormap" ||
            name == "colorlevels") {
            return  m_colorbar->setPropertyValue(name, value);
        }

        bool ret = Object::setPropertyValue(name, value);
        if (name == "fontname" || name == "fontsize" || name == "fontweight" || name == "fontangle") {
            m_secYAxis->setPropertyValue(name, value);
        }
        else if (name == "polartiptotail") {
            vector<double> children = getAllChildren();
            vector<double>::iterator it = children.begin();
            for (; it != children.end(); ++it) {
                Object* obj = getObject(*it);
                if (obj->getObjectType() == ObjectType::POLAR_LINE) {
                    obj->setPropertyValue("tiptotail", value);
                }
            }
        }
        else  if (name == "xlim") {
            setAxisNeedRepaint(0);
        }
        else if (name == "ylim") {
            setAxisNeedRepaint(1);
        }
        else if (name == "zlim") {
            setAxisNeedRepaint(3);
        }
        else if (name == "xtick") {
            if (value.isCurrency() && value.getCurrency().IsScalar()) {
                m_xdateTicks = false;
            }
            else {
                if (m_xdateTicks)
                    CoreMain::getInstance()->setUpdateDatetickFlag(getHandle(), "x");
            }
        }
        else if (name == "ytick") {
            if (value.isCurrency() && value.getCurrency().IsScalar()) {
                m_ydateTicks = false;
            }
            else {
                if (m_ydateTicks)
                    CoreMain::getInstance()->setUpdateDatetickFlag(getHandle(), "y");
            }
        }
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

    void Axes::setView(const std::vector<double>& view)
    {
        int N = (int)view.size();
        if (N == 1)
        {
            if (view[0] == 2)
            {
                m_rotX = 0;
                m_rotZ = 0;
            }
            else if (view[0] == 3)
            {
                m_rotX = 60;
                m_rotZ = 315;
            }
            else
            {
                throw OML_Error(OML_ERR_OPTIONVAL);
            }
        }
        else if (N == 2)
        {
            // azimuth - elevation
            m_rotX = std::fmod(360.0 + 90.0 - view[1], 360.0);
            m_rotZ = std::fmod(360 + view[0], 360);
        }
        else if (N == 3)
        {
            // x, y, z -> calculate azimuth - elevation
            double x = view[0], y = view[1], z = view[2];
            const double PI = acos(double(-1));
            double az = 0, el = 0;
            if (isSameDouble(x, 0) && isSameDouble(y, 0))
            {
                az = 0;
                el = 0;
                if (z > 0)
                    el = 90;
                else if (z < 0)
                    el = -90;
            }
            else
            {
                // one of x and y is not zero
                // calculate azimuth first
                if (isSameDouble(x, 0))
                {
                    az = y > 0 ? 180 : 0;
                }
                else if (isSameDouble(y, 0))
                {
                    az = x > 0 ? 90 : -90;
                }
                else
                {
                    az = atan(x / (y * -1)) * 180 / PI;
                    if (y > 0)
                        az += 180;
                }
                // calculate elevation
                double zz = sqrt(x * x + y * y);
                el = atan(z / zz) * 180 / PI;
            }
            std::vector<double> d;
            d.push_back(az);
            d.push_back(el);
            setView(d);
        }
    }

    void Axes::setLegendVisible(bool visible)
    {
        m_legend->setVisible(visible);
    }

    bool Axes::getLegendVisible() const
    {
        return m_legend->getVisible();
    }

    double Axes::getLegendHandle() const
    {
        return m_legend->getHandle();
    }

    double Axes::getColorbarHandle() const
    {
        return m_colorbar->getHandle();
    }

    string Axes::drawStackedBarPlotLabels(GnuplotOutput* out)
    {
        stringstream ss;
        bool sumInitialized = false;
        vector<double> dataX;
        vector<double> sumY;

        vector<double> children = getAllChildren();
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

    void Axes::setAxisRangeFromData()
    {
        vector<double> range;
        vector<double> children = getAllChildren();
        for (vector<double>::iterator it = children.begin(); it != children.end(); ++it)
        {
            Object* obj = getObject(*it);
            if (!obj)
                continue;

            int rS = (int)range.size();
            vector<double> objRange = obj->getMinMaxData();
            if (rS == 0) {
                range = objRange;
                continue;
            }
            
            if ((int)objRange.size() >= 4 && rS >=4) {
                if (objRange[0] < range[0])
                    range[0] = objRange[0];
                if (objRange[1] > range[1])
                    range[1] = objRange[1];
                if (objRange[2] < range[2])
                    range[2] = objRange[2];
                if (objRange[3] > range[3])
                    range[3] = objRange[3];
            }

            if ((int)objRange.size() == 6) {
                if (rS == 4) {
                    range.push_back(objRange[4]);
                    range.push_back(objRange[5]);
                }
                else if (rS == 6) {
                    if (objRange[4] < range[4])
                        range[4] = objRange[4];
                    if (objRange[5] > range[5])
                        range[5] = objRange[5];
                }
            }
        }
    
        if ((int)range.size() >= 4) {
            vector<double> xr = { range[0], range[1] };
            vector<double> yr = { range[2], range[3] };
            setPropertyValue("xlim", xr);
            setPropertyValue("ylim", yr);
        }
        if ((int)range.size() >= 6) {
            vector<double> zr = { range[4], range[5] };
            setPropertyValue("zlim", zr);
        }
    }

    bool Axes::is3DPlot()
    {
        bool is3D = false;
        vector<double> children = getAllChildren();
        vector<double>::iterator it = children.begin();
        for (; it != children.end(); ++it) {
            Object* obj = getObject(*it);
            if (!obj)
                continue;
            ObjectType ot = obj->getObjectType();
            if (ot == ObjectType::CONTOUR_3D || ot == ObjectType::LINE_3D ||
                ot == ObjectType::MESH || ot == ObjectType::SCATTER_3D ||
                ot == ObjectType::PCOLOR || ot == ObjectType::SURFACE ||
                ot == ObjectType::WATERFALL) {
                
                return true;
            }
        }
        return false;
    }

    void Axes::initColorbarRange()
    {
        std::vector<double> clim = getColorbarRange();
        m_colorbar->setPropertyValue("clim", clim);
    }

    void Axes::setAxisOption(const string& option)
    {
        if (option == "on") {
            _axesOn = true;
        }
        else if (option == "off") {
            _axesOn = false;
        }
        else if (option == "tight") {
            setAxisRangeFromData();
        }
        else if (option == "equal") {
            m_axisOption = "set size ratio -1\n";
        }
        else if (option == "normal") {
            m_axisOption = "set size noratio\n";
        }
        else if (option == "square") {
            m_axisOption = "set size square\n";
        }
        else if (option == "cubical") {
            m_axisOption = "set view noequal_axes\n";
        }
        else if (option == "unscaled") {
            m_axisOption = "set view equal xyz\n";
        }
        else {
            m_axisOption = "";
        }
    }

    void Axes::getAxesPosition(double& x, double& y, double& w, double& h)
    {
        vector<double> pos = getPropertyValue("position").Vector();
        x = pos[0];
        y = pos[1];
        w = pos[2];
        h = pos[3];
        Object* pobj = getParentObject();
        if (pobj && pobj->isFigure()) {
            Figure* fig = static_cast<Figure*>(pobj);
            if (fig->isGridLayout() || (x == 0 && y == 0 && w == 1 && h == 1)) {
                fig->getSubplotPosition(m_gridIndex, x, y, w, h);
            }
        }
    }

    void Axes::setAxisDatetick(const string& axis, const std::string& datefmt, int datefmtIdx)
    {
        if (axis == "x") {
            m_xdateTickFmt = datefmt;
            m_xdateTickFmtIdx = datefmtIdx;
            m_xdateTicks = true;
        }
        else if (axis == "y") {
            m_ydateTickFmt = datefmt;
            m_ydateTickFmtIdx = datefmtIdx;
            m_ydateTicks = true;
        }
    }

    void Axes::getAxisDatetickOptions(const string& axis, bool& enabled, std::string& fmt, int& fmtIdx)
    {
        if (axis == "x") {
            enabled = m_xdateTicks;
            fmt = m_xdateTickFmt;
            fmtIdx = m_xdateTickFmtIdx;
        }
        else if (axis == "y") {
            enabled = m_ydateTicks;
            fmt = m_ydateTickFmt;
            fmtIdx = m_ydateTickFmtIdx;
        }
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

        _xcolcount = ld.xcolcount;
        _ycolcount = ld.ycolcount;
        _zcolcount = ld.zcolcount;

        for (int i = 0; i < ld.properties.size(); i++){
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    void Drawable::update(GnuplotOutput *){
    }

    bool Drawable::isDrawable()
    {
        return true;
    }

    vector<double> Drawable::getXdata(){
        Currency val = getPropertyValue("xdata").getCurrency();
        if (val.IsVector()) {
            return val.Vector();
        }
        else if (val.IsScalar()) {
            vector<double> ret;
            ret.push_back(val.Scalar());
            return ret;
        }
        return vector<double>();
    }

    vector<double> Drawable::getYdata(){
        Currency val = getPropertyValue("ydata").getCurrency();
        if (val.IsVector()) {
            return val.Vector();
        }
        else if (val.IsScalar()) {
            vector<double> ret;
            ret.push_back(val.Scalar());
            return ret;
        }
        return vector<double>();
    }

    vector<double> Drawable::getZdata(){
        std::vector<Property>::iterator fit = std::find_if(m_ps.begin(),
            m_ps.end(), [](Property p) {return p.getName() == "zdata"; });
        if (fit == m_ps.cend())
            return vector<double>();

        Currency val = getPropertyValue("zdata").getCurrency();
        if (val.IsVector()) {
            return val.Vector();
        }
        else if (val.IsScalar()) {
            vector<double> ret;
            ret.push_back(val.Scalar());
            return ret;
        }
        return vector<double>();
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

    vector<double> Drawable::getMinMaxData()
    {
        vector<double> ret;
        vector<double> x = getXdata();
        vector<double> y = getYdata();
        double xOffset = 0, xScale = 1;
        double yOffset = 0, yScale = 1;

        vector<string> propNames = getPropertyNames();
        if (std::find(propNames.cbegin(), propNames.cend(), "dataxoffset") != propNames.cend())
            xOffset = getPropertyValue("dataxoffset").Scalar();
        if (std::find(propNames.cbegin(), propNames.cend(), "datayoffset") != propNames.cend())
            yOffset = getPropertyValue("datayoffset").Scalar();
        if (std::find(propNames.cbegin(), propNames.cend(), "dataxscale") != propNames.cend())
            xScale = getPropertyValue("dataxscale").Scalar();
        if (std::find(propNames.cbegin(), propNames.cend(), "datayscale") != propNames.cend())
            yScale = getPropertyValue("datayscale").Scalar();

        vector<double> z;
        if (std::find(propNames.cbegin(), propNames.cend(), "zdata") != propNames.cend()) {
             z = getZdata();
        }

        if (z.empty()) {
            size_t count = min(x.size(), y.size());
            double minX = 0, maxX = 1, minY = 0, maxY = 1;
            if (count > 0) {
                minX = maxX = x[0];
                minY = maxY = y[0];
            }
            double tmpX, tmpY;
            for (int j = 0; j < count; j++) {
                tmpX = x[j] * xScale + xOffset;
                tmpY = y[j] * yScale + yOffset;
                if (tmpX < minX)
                    minX = tmpX;
                if (tmpX > maxX)
                    maxX = tmpX;
                if (tmpY < minY)
                    minY = tmpY;
                if (tmpY > maxY)
                    maxY = tmpY;
            }
            ret.push_back(minX);
            ret.push_back(maxX);
            ret.push_back(minY);
            ret.push_back(maxY);
        }
        else {
            double zOffset = 0, zScale = 1;
            if (std::find(propNames.cbegin(), propNames.cend(), "datazoffset") != propNames.cend())
                zOffset = getPropertyValue("datazoffset").Scalar();
            if (std::find(propNames.cbegin(), propNames.cend(), "datazscale") != propNames.cend())
                zScale = getPropertyValue("datazscale").Scalar();

            size_t count = min(x.size(), y.size());
            count = min(y.size(), z.size());
            double minX = 0, maxX = 1, minY = 0, maxY = 1, minZ = 0, maxZ = 1;
            if (count > 0) {
                minX = maxX = x[0];
                minY = maxY = y[0];
                minZ = maxZ = z[0];
            }
            double tmpX, tmpY, tmpZ;
            for (int j = 0; j < count; j++) {
                tmpX = x[j] * xScale + xOffset;
                tmpY = y[j] * yScale + yOffset;
                tmpZ = z[j] * zScale + zOffset;
                if (tmpX < minX)
                    minX = tmpX;
                if (tmpX > maxX)
                    maxX = tmpX;
                if (tmpY < minY)
                    minY = tmpY;
                if (tmpY > maxY)
                    maxY = tmpY;
                if (tmpZ < minZ)
                    minZ = tmpZ;
                if (tmpZ > maxZ)
                    maxZ = tmpZ;
            }
            ret.push_back(minX);
            ret.push_back(maxX);
            ret.push_back(minY);
            ret.push_back(maxY);
            ret.push_back(minZ);
            ret.push_back(maxZ);
        }
        return ret;
    }

    void Drawable::setXAxisReference(const string& xref)
    {
        m_xaxisRef = xref;
    }

    void Drawable::setYAxisReference(const string& yref)
    {
        m_yaxisRef = yref;
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

    string Line::getLegend(bool bold, bool italic){

        string legend = getPropertyValue("displayname").StringVal();
        stringstream ss;
        ss << "title \"";
        if (italic) {
            ss << "{/:Italic ";
        }
        if (bold) {
            ss << "{/:Bold ";
        }
        ss << legend;
        if (bold) {
            ss << "}";
        }
        if (italic) {
            ss << "}";
        }
        ss << "\"";
        return ss.str();
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
        }
        else if (linestyle == "-:") {
            dt = 5;
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

    string Fill::getLegend(bool bold, bool italic){
        string legend = getPropertyValue("displayname").StringVal();
        stringstream ss;
        ss << "title \"";
        if (italic) {
            ss << "{/:Italic ";
        }
        if (bold) {
            ss << "{/:Bold ";
        }
        ss << legend;
        if (bold) {
            ss << "}";
        }
        if (italic) {
            ss << "}";
        }
        ss << "\"";
        return ss.str();
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
        Property& ap = getProperty("areagroup");
        ap.setValue(getPropertyValue("handle").Scalar());

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string Area::getUsingClause(){
        return "'-' using 1:2";
    }

    string Area::getLegend(bool bold, bool italic){
        
        string legend = getPropertyValue("displayname").StringVal();
        stringstream ss;
        ss << "title \"";
        if (italic) {
            ss << "{/:Italic ";
        }
        if (bold) {
            ss << "{/:Bold ";
        }
        ss << legend;
        if (bold) {
            ss << "}";
        }
        if (italic) {
            ss << "}";
        }
        ss << "\"";
        return ss.str();
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
        Property& bp = getProperty("bargroup");
        bp.setValue(getPropertyValue("handle").Scalar());

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
        bool legendBold = false;
        bool legendItalic = false;
        Object* p = getParentObject();
        Axes* pAxes = p && p->getObjectType() == ObjectType::AXES ? static_cast<Axes*>(p) : nullptr;
        if (pAxes) {
            Object* lg = getObject(pAxes->getLegendHandle());
            if (lg && lg->getObjectType() == ObjectType::LEGEND) {
                legendBold = lg->getPropertyValue("fontweight").StringVal() == "bold";
                legendItalic = lg->getPropertyValue("fontangle").StringVal() == "italic";
            }
        }

        stringstream ss;
        ss << "'-' using 2:xtic(1) lc rgb '" << color << "' fill pattern " << f <<" title \"";
        if (legendItalic)   ss << "{/:Italic ";
        if (legendBold)     ss << "{/:Bold ";
        ss << legend;
        if (legendBold)     ss << "}";
        if (legendItalic)   ss << "}";
        ss << "\"";

        if (pAxes && pAxes->getPropertyValue("barlabels").StringVal() == "on" && 
            getPropertyValue("barlayout").StringVal() == "grouped") {

            int numBars = 0;
            int barIdx = -1;
            pAxes->getBarNumberAndIdx(getPropertyValue("handle").Scalar(), numBars, barIdx);

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

    string HggroupBar::getLegend(bool bold, bool italic){
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
        m_ps.push_back(Property("markerfacealpha", 1.0, PropertyType::DOUBLE));

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

    string HggroupScatter::getLegend(bool bold, bool italic){
        string legend = getPropertyValue("displayname").StringVal();
        stringstream ss;
        ss << "title \"";
        if (italic) {
            ss << "{/:Italic ";
        }
        if (bold) {
            ss << "{/:Bold ";
        }
        ss << legend;
        if (bold) {
            ss << "}";
        }
        if (italic) {
            ss << "}";
        }
        ss << "\"";
        return ss.str();
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
        double markerAlpha = getPropertyValue("markerfacealpha").Scalar();
        if (markerAlpha >= 0 && markerAlpha <= 1)
            markerColor = getPropertyValue("markerfacecolor").ColorStringWithAlpha(markerAlpha);

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
        :Drawable(), _minZ(0), _maxZ(1), m_recalcMinMax(true)
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
        if (p && p->isAxes()) {
            colorOrder = p->getPropertyValue("colororder").getCurrency();
            static_cast<Axes*>(p)->initColorbarRange();
        }

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

    string Surface::getLegend(bool bold, bool italic){
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
        ss << "unset contour\n";
        return ss.str();
    }

    void Surface::cleanup(GnuplotOutput *out) {
    }

    bool Surface::getMinMaxZ(double& min, double& max)
    {
        if (m_recalcMinMax) {
            vector<double> z = getZdata();
            _minZ = 0;
            _maxZ = 1;
            if (!z.empty())
            {
                _minZ = z[0];
                _maxZ = z[0];
            }
            for (int i = 0; i < (int)z.size(); i++) {
                double zVal = z[i];
                if (zVal < _minZ)   
                    _minZ = zVal;
                if (zVal > _maxZ)   
                    _maxZ = zVal;
            }
            m_recalcMinMax = false;
        }
        min = _minZ;
        max = _maxZ;
        return true;
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
           << "unset contour\n";
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

    Text::Text(bool canBeDeleted)
        :Object(), m_canBeDeleted(canBeDeleted)
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
        m_ps.push_back(Property("backgroundcolor", Color(string("white")), PropertyType::COLOR));

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
        string bgcolor = getPropertyValue("backgroundcolor").ColorString();
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

        out->printf("set style textbox %d opaque margins 6.0, 2 fc rgb \"%s\"\n", m_index, bgcolor);
        out->printf("set label \"%s\" at %g,%g,%g %s offset %g,%g point ps %d font \"%s,%g\" tc rgb \"%s\" enhanced front boxed bs %d\n",
                    ss.str().c_str(), x, y, z, align.c_str(), offsetX, offsetY, pointSize,
                    fontname.c_str(), fontsize, color.c_str(), m_index);
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

    std::string CurrencyAndColor::ColorStringWithAlpha(double alpha)
    {
        return m_color.getStringWithAlpha(alpha);
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
        :Object(), m_ylabel(new Text(false)), m_parentAxesHandle(parentHandle),
        m_ydateTicks(false)
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
        m_ps.push_back(Property("yminortick", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("ytick", 0.0, PropertyType::VEC_DOUBLE));
        m_ps.push_back(Property("yticklabel", string(""), PropertyType::CELL));
        
        // not yet supported properties
        m_ps.push_back(Property("ytickmethod", string("ytickmethod"), PropertyType::UNSUPPORTED));

        m_objectMap[handle] = this;
    }

    void SecondaryYAxis::update(GnuplotOutput* out)
    {
        if (getPropertyValue("visible").StringVal() == "off")
            return;

        Axes* ax = getParentObject() ? dynamic_cast<Axes*>(getParentObject()):nullptr;
        bool axisVisible = ax ? ax->getAxesOn() : true;
        if (axisVisible) {
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
            int ytick = (int)getPropertyValue("yminortick").Scalar();
            if (ytick > 0)
                out->printf("set my2tics %d\n", ytick);

            Currency tmp = getPropertyValue("ytick").getCurrency();
            if (tmp.IsScalar()) {
                if (!IsZero(tmp.Scalar()))
                    out->printf("set y2tics %g\n", tmp.Scalar());
        }
        else {
                vector<double> yticks = tmp.Vector();
                Currency tickLbl = getPropertyValue("yticklabel").getCurrency();
                if (tickLbl.IsCellArray() && tickLbl.CellArray()->Size()>0) {
                    HML_CELLARRAY* ca = tickLbl.CellArray();
                    int S = min((int)yticks.size(), ca->Size());
                    out->printf("set y2tics (");
                    for (int i = 0; i < S; ++i) {
                        string lbl = "";
                        if ((*ca)(i).IsString())
                            lbl = (*ca)(i).StringVal();
                        out->printf("\"%s\" %g, ", lbl.c_str(), yticks[i]);
        }
                    out->printf(")\n");
        }
        else {
                    vector<double>::const_iterator it = yticks.cbegin();
                    out->printf("set y2tics (");
                    for (; it != yticks.cend(); ++it) {
                        out->printf("%g,", *it);
        }
                    out->printf(")\n");
        }
        }
        }

        if (getPropertyValue("yscale").StringVal() == "log")
            out->printf("set logscale y2\n");

        string ygrid = getPropertyValue("ygrid").StringVal();
        if (ygrid == "on") 
            out->printf("set grid y2tics\n");
        
        string ymgrid = getPropertyValue("yminorgrid").StringVal();
        if (ymgrid == "on")
            out->printf("set grid my2tics\n");
        
        if (axisVisible)
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
    
    void SecondaryYAxis::setAxisDatetick(const std::string& datefmt, int datefmtIdx)
    {
        m_ydateTickFmt = datefmt;
        m_ydateTickFmtIdx = datefmtIdx;
        m_ydateTicks = true;
    }

    void SecondaryYAxis::getAxisDatetickOptions(bool& enabled, std::string& fmt, int& fmtIdx)
    {
        enabled = m_ydateTicks;
        fmt = m_ydateTickFmt;
        fmtIdx = m_ydateTickFmtIdx;
    }

    bool SecondaryYAxis::setPropertyValue(const string& name, VALUETYPE value)
    {
        Object::setPropertyValue(name, value);
        if (name == "ytick") {
            if (value.isCurrency() && value.getCurrency().IsScalar()) {
                m_ydateTicks = false;
            }
            else {
                if (m_ydateTicks)
                    CoreMain::getInstance()->setUpdateDatetickFlag(getHandle(), "x");
            }
        }
        return true;
    }
    
    XLine::XLine()
        :Line()
    {
        m_type = ObjectType::XLINE;
    }

    void XLine::init(const LineData& ld)
    {
        Line::init(ld);

        setPropertyValue("xdata", ld.y);
        setPropertyValue("ydata", ld.x);
        setPropertyValue("zdata", ld.z);

        for (int i = 0; i < ld.properties.size(); i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }

        _xcolcount = ld.ycolcount;
        _ycolcount = ld.xcolcount;
        _zcolcount = ld.zcolcount;
    }
    
    string XLine::getUsingClause(){
        vector<double> x = getXdata();
        double xVal = x.empty() ? 0 : x.front();
        //set arrow 10 from 3, graph 0 to 3, graph 1 nohead
        stringstream ss;
        ss << "set arrow from " << xVal << ", graph 0 to "
            << xVal << ", graph 1 nohead";
        return ss.str();
    }

    string XLine::getWithClause(int lineId)
    {
        stringstream ss;
        ss << " linestyle " << lineId << "\n";
        return ss.str();
    }
    
    YLine::YLine()
        :Line()
    {
        m_type = ObjectType::YLINE;
    }

    string YLine::getUsingClause(){
        vector<double> y = getYdata();
        double yVal = y.empty() ? 0 : y.front();
        //set arrow from graph 0,first 5.6 to graph 1, first 5.6
        stringstream ss;
        ss << "set arrow from graph 0, first " << yVal 
            << " to graph 1, first "<< yVal << " nohead";
        return ss.str();
    }
    string YLine::getWithClause(int lineId)
    {
        stringstream ss;
        ss << " linestyle " << lineId << "\n";
        return ss.str();
    }

    Waterfall::Waterfall()
        :Surface()
    {
        m_type = ObjectType::WATERFALL;
        m_ps.push_back(Property("linewidth", 1.0, PropertyType::DOUBLE));
    }
    string Waterfall::getUsingClause()
    {
        vector<double> y = getYdata();
        vector<double> z = getZdata();

        double _minZ = 0;
        vector<double> yVals;

        int numOfSlices = (int)z.size() / _zcolcount;
        int numOfPoints = _zcolcount;

        int k = 0;
        for (int i = 0; i < numOfSlices; ++i) {
            double yPoint = y[i];
            yVals.push_back(i < y.size() ? y[i] : 0.0);
            for (int j = 0; j < numOfPoints; ++j, ++k)
            {
                if (k < z.size()) {
                    double zVal = z[k];
                    // find min z
                    if (zVal < _minZ)
                        _minZ = zVal;
                }
            }
        }

        //std::string style("fs solid 0.2 fc 'grey' notitle linestyle");
        std::string style("fc '#eeeeee' linestyle ");
        
        std::stringstream ss;
        for (int i = ((int)yVals.size()-1); i > 0; --i) {
            ss << "'-' using 1:(" << yVals[i] << "):2:(" << _minZ << "):2 with zerrorfill "
                << style << m_lineId << " notitle ,";
        }
        ss << "'-' using 1:(" << yVals[0] << "):2:(" << _minZ << "):2 with zerrorfill "
            <<style << m_lineId;
        return ss.str();
    }
    string Waterfall::getWithClause(int line)
    {
        m_lineId = line;
        return string();
    }
    string Waterfall::getLineStyle(int line)
    {
        m_lineId = line;

        string lineColor = getPropertyValue("color").ColorString();
        int lineWidth = int(getPropertyValue("linewidth").Scalar());

        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line << " lc rgb '" 
            << lineColor << "' lw " << lineWidth << "\n";
        return ss.str();
    }
    void Waterfall::putData(GnuplotOutput* out)
    {
        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> z = getZdata();

        if (z.empty())
        {
            out->printf("e\n");
            return;
        }

        int numOfSlices = (int)z.size() / _zcolcount;
        int numOfPoints = _zcolcount;
        bool matInput = _xcolcount > 1 && _ycolcount > 1;
        for (int i = (numOfSlices - 1); i >= 0; --i) {
            for (int j = 0; j < numOfPoints; ++j)
            {
                int idx = j * numOfSlices + i;
                double xVal = 0;
                if (matInput) {
                    if (idx < x.size())
                        xVal = x[idx];
                }
                else {
                    if (j < x.size())
                        xVal = x[j];
                }
                
                double zVal = idx < z.size() ? z[idx] : 0.0;
                out->printf("%g %g\n", xVal, zVal);
            }
            out->printf("e\n");
        }
    }
    Ellipse::Ellipse()
        :Drawable()
    {
        m_type = ObjectType::ELLIPSE;

        // remove "Drawable" properties
        std::vector<string> propsToRemove = { "displayname", "xdata", "ydata" , "zdata" };
        std::vector<string>::const_iterator it = propsToRemove.cbegin();
        for (; it != propsToRemove.cend(); ++it)
        {
            string tmp = *it;
            std::vector<Property>::iterator fit = std::find_if(m_ps.begin(),
                m_ps.end(), [&tmp](Property p) {return p.getName() == tmp; });
            if (fit != m_ps.end())
                m_ps.erase(fit);
        }

        m_ps.push_back(Property("type", string("hggroup"), PropertyType::STRING));
        m_ps.push_back(Property("facecolor", Color(string("blue")), PropertyType::UNSUPPORTED));
        m_ps.push_back(Property("edgecolor", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("linestyle", string("-"), PropertyType::STRING));
        m_ps.push_back(Property("linewidth", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        vector<double> pos = { 0, 0, 1, 1 };
        m_ps.push_back(Property("position", pos, PropertyType::STRING));
    }

    void Ellipse::init(const LineData& ld)
    {
        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        // update line-specific properties
        LineStyle ls = LineStyle(ld, colorOrder);
        setPropertyValue("linestyle", ls.m_lineStyle);
        setPropertyValue("edgecolor", ls.m_lineColor);
        setPropertyValue("linewidth", ls.m_lineWidth);
        setPropertyValue("facecolor", ls.m_lineColor);

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string Ellipse::getUsingClause()
    {
        // style 
        string linestyle = getPropertyValue("linestyle").StringVal();
        int dt = 1;
        if (linestyle == "-")
            dt = 1;
        else if (linestyle == "--")
            dt = 2;
        else if (linestyle == ":")
            dt = 3;
        else if (linestyle == "-.")
            dt = 4;
        else if (linestyle == "")
            dt = 0;

        string edgeColor = getPropertyValue("edgecolor").ColorString();
        int lineWidth = int(getPropertyValue("linewidth").Scalar());

        vector<double> pos = getPropertyValue("position").Vector();
        if (pos.size() < 4)
        {
            for (size_t i = pos.size(); i < 4; ++i)
                pos.push_back(0);
        }
        std::stringstream ss;

        double lX = pos[2] / 2;
        double lY = pos[3] / 2;

        // draw the border
        ss << "set object ellipse center " << pos[0] + lX << "," << pos[1] + lY << " size "
            << lX << "," << lY;
        ss << " dt " << dt << " lw " << lineWidth
            << " fs empty border lc rgb '" << edgeColor << "'\n";

        return ss.str();
    }

    string Ellipse::getWithClause(int lineId)
    {
        return string();
    }

    string Ellipse::getLineStyle(int line)
    {
        return string();
    }

    void Ellipse::cleanup(GnuplotOutput*)
    {
    }

    string Ellipse::getDummyLine()
    {
        vector<double> pos = getPropertyValue("position").Vector();
        if (pos.size() < 4)
        {
            for (size_t i = pos.size(); i < 4; ++i)
                pos.push_back(0);
        }
        std::stringstream ss;
        ss << "plot '-' using 1:2 notitle ps 0\n"
            << pos[0] << " " << pos[1] << "\n"
            << pos[0] + pos[2] << " " << pos[1] + pos[3] << "\n"
            << "e\n";

        return ss.str();
    }

    Rectangle::Rectangle()
        :Drawable(), m_fillNone(false)
    {
        m_type = ObjectType::RECTANGLE;

        // remove "Drawable" properties
        std::vector<string> propsToRemove = { "displayname", "xdata", "ydata" , "zdata" };
        std::vector<string>::const_iterator it = propsToRemove.cbegin();
        for (; it != propsToRemove.cend(); ++it)
        {
            string tmp = *it;
            std::vector<Property>::iterator fit = std::find_if(m_ps.begin(),
                m_ps.end(), [&tmp](Property p) {return p.getName() == tmp; });
            if (fit != m_ps.end())
                m_ps.erase(fit);
        }

        m_ps.push_back(Property("type", string("hggroup"), PropertyType::STRING));
        m_ps.push_back(Property("facecolor", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("edgecolor", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("linestyle", string("-"), PropertyType::STRING));
        m_ps.push_back(Property("linewidth", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        vector<double> pos = { 0, 0, 1, 1 };
        m_ps.push_back(Property("position", pos, PropertyType::STRING));
        m_ps.push_back(Property("curvature", 1, PropertyType::UNSUPPORTED));
    }

    void Rectangle::init(const LineData& ld)
    {
        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        // update line-specific properties
        LineStyle ls = LineStyle(ld, colorOrder);
        setPropertyValue("linestyle", ls.m_lineStyle);
        setPropertyValue("edgecolor", ls.m_lineColor);
        setPropertyValue("linewidth", ls.m_lineWidth);
        setPropertyValue("facecolor", ls.m_lineColor);

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }
    
    string Rectangle::getUsingClause()
    {
        // style 
        string linestyle = getPropertyValue("linestyle").StringVal();
        int dt = 1;
        if (linestyle == "-")
            dt = 1;
        else if (linestyle == "--")
            dt = 2;
        else if (linestyle == ":")
            dt = 3;
        else if (linestyle == "-.")
            dt = 4;
        else if (linestyle == "")
            dt = 0;
        
        string edgeColor = getPropertyValue("edgecolor").ColorString();
        int lineWidth = int(getPropertyValue("linewidth").Scalar());
        string faceColor = getPropertyValue("facecolor").ColorString();

        vector<double> pos = getPropertyValue("position").Vector();
        if (pos.size() < 4)
        {
            for (size_t i = pos.size(); i < 4; ++i)
                pos.push_back(0);
        }
        std::stringstream ss;

        if (!m_fillNone) {
            // draw the inner
            ss << "set object rect from " << pos[0] << "," << pos[1] << " to "
                << pos[0] + pos[2] << "," << pos[1] + pos[3];
            ss << " fc '" << faceColor << "' fs noborder\n";
        }

        // draw the border
        ss << "set object rect from " << pos[0] << "," << pos[1] << " to "
            << pos[0] + pos[2] << "," << pos[1] + pos[3];
        ss << " dt " << dt << " lw " << lineWidth
            << " fs empty border lc rgb '" << edgeColor << "'\n";
        
        return ss.str();
    }
    string Rectangle::getWithClause(int lineId)
    {
        return string();
    }
    string Rectangle::getLineStyle(int line)
    {
        return string();
    }

    void Rectangle::cleanup(GnuplotOutput*)
    {
    }

    void Rectangle::SetFillNone(bool set)
    {
        m_fillNone = set;
    }

    string Rectangle::getDummyLine()
    {
        vector<double> pos = getPropertyValue("position").Vector();
        if (pos.size() < 4)
        {
            for (size_t i = pos.size(); i < 4; ++i)
                pos.push_back(0);
        }
        std::stringstream ss;
        ss << "plot '-' using 1:2 notitle ps 0 \n"
            << pos[0] << " " << pos[1] << "\n"
            << pos[0] + pos[2] << " " << pos[1] + pos[3] << "\n"
            << "e\n";
        return ss.str();
    }

    PColor::PColor()
        :Surface()
    {
        m_type = ObjectType::PCOLOR;

        // remove "Drawable" properties
        string tmp("meshlines");
        std::vector<Property>::iterator fit = std::find_if(m_ps.begin(),
            m_ps.end(), [&tmp](Property p) {return p.getName() == tmp; });
        if (fit != m_ps.end())
            m_ps.erase(fit);
        m_ps.push_back(Property("meshlines", string("on"), PropertyType::UNSUPPORTED));
    }
    string PColor::getWithClause(int)
    {
        return "with image";
    }
    
    void PColor::putData(GnuplotOutput* out)
    {
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

        if (_xcolcount > 1 && _ycolcount > 1) {
            size_t xsize = x.size();
            size_t ysize = y.size();
            size_t zsize = z.size();
            if (!((xsize == ysize) &&
                (xsize == zsize))) {
                throw;
            }

            int M = (int)z.size() / _zcolcount;
            int N = _zcolcount;
            for (int i = 0; i < M - 1; ++i) {
                for (int j = 0; j < N - 1; ++j)
                {
                    int idx = i * N + j;
                    double tmp = x[idx];
                    double xVal = tmp + (x[(i + 1) * N + j] - tmp) / 2.0;
                    tmp = y[idx];
                    double yVal = tmp + (y[i * N + j + 1] - tmp) / 2.0;
                    double zVal = z[idx];
                    out->printf("%g %g %g\n", xVal, yVal, zVal);

                    // find min/max for the colorbar
                    if (zVal < _minZ)
                        _minZ = zVal;
                    if (zVal > _maxZ)
                        _maxZ = zVal;
                }
            }
        }
        else {
            size_t col = y.size();
            for (int xsub = 0; xsub < x.size() - 1; xsub++) {

                for (int ysub = 0; ysub < col - 1; ysub++) {
                    double xVal = x[xsub] + (x[xsub + 1] - x[xsub]) / 2.0;
                    double yVal = y[ysub] + (y[ysub + 1] - y[ysub]) / 2.0;
                    double zVal = z[xsub * col + ysub];
                    out->printf("%g %g %g\n", xVal, yVal, zVal);

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

    Patch::Patch()
        :Drawable(), m_is3D(false)
    {
        m_type = ObjectType::PATCH;

        m_ps.push_back(Property("facecolor", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("type", string("patch"), PropertyType::STRING));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("dataxoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("dataxscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datazoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datazscale", 1.0, PropertyType::DOUBLE));
    }

    void Patch::init(const LineData& ld)
    {
        setPropertyValue("xdata", ld.x);
        setPropertyValue("ydata", ld.y);
        setPropertyValue("zdata", ld.z);

        _xcolcount = ld.xcolcount;
        _ycolcount = ld.ycolcount;
        _zcolcount = ld.zcolcount;

        m_is3D = !ld.z.empty();

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        LineStyle ls = LineStyle(ld, colorOrder);
        setPropertyValue("facecolor", ls.m_lineColor);
        setPropertyValue("displayname", ls.m_legend);

        int vertIdx = -1;
        int facesIdx = -1;
        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++) {
            if (ld.properties[i] == "vertices")
            {
                vertIdx = i;
                continue;
            }
            if (ld.properties[i] == "faces")
            {
                facesIdx = i;
                continue;
            }
            setPropertyValue(ld.properties[i], ld.values[i]);
        }

        if (vertIdx >= 0 && facesIdx >= 0 && 
            ld.values[vertIdx].IsMatrix() && ld.values[facesIdx].IsMatrix())
        {
            const hwMatrix* vertices = ld.values[vertIdx].Matrix();
            const hwMatrix* faces = ld.values[facesIdx].Matrix();

            vector<double> xdata, ydata, zdata;
            xdata.reserve(faces->Size());
            ydata.reserve(faces->Size());
            if (vertices->N() == 3)
            {
                m_is3D = true;
                zdata.reserve(faces->Size());
            }

            int numFaces = faces->M();
            int numPoints = faces->N();
            // for each face create the matrices x and y
            for (int i = 0; i < faces->M(); ++i)
            {
                hwMatrix fi;
                faces->ReadRow(i, fi);
                int numPoints = 0;
                for (int j = 0; j < fi.Size(); ++j)
                {
                    int vidx = (int)fi(j);
                    if ((vidx - 1) >= 0 && (vidx - 1) < vertices->M())
                    {
                        numPoints++;
                        xdata.push_back((*vertices)(vidx - 1, 0));
                        ydata.push_back((*vertices)(vidx - 1, 1));
                        if (m_is3D)
                            zdata.push_back((*vertices)(vidx - 1, 2));
                    }
                    else
                    {
                        xdata.push_back(std::numeric_limits<double>::quiet_NaN());
                        ydata.push_back(std::numeric_limits<double>::quiet_NaN());
                        if (m_is3D)
                            zdata.push_back(std::numeric_limits<double>::quiet_NaN());
                    }
                }
            }

            setPropertyValue("xdata", xdata);
            setPropertyValue("ydata", ydata);
            setPropertyValue("zdata", zdata);

            _xcolcount = numFaces;
        }

        if (m_is3D)
        {
            m_ps.push_back(Property("cdata", Currency(), PropertyType::MATRIX));
        }

    }

    string Patch::getUsingClause()
    {
        std::stringstream ss;
        string lineColor = getPropertyValue("facecolor").ColorString();
        if (m_is3D)
        {
            ss << "'-' using 1:2:3 with polygons fc rgb \""<<lineColor<<"\"";
            }
        else
        {
           ss << "'-' using 1:2 with polygons fs solid lc rgb '" << lineColor << "'";
        }
        return ss.str();
    }

    string Patch::getWithClause(int lineId)
    {
        return string();
    }

    string Patch::getLineStyle(int line)
    {
        if (m_is3D)
            return string("set pm3d depthorder base\nset pm3d border lc \"black\"\n");
        return string();
        }

    void Patch::putData(GnuplotOutput* out)
    {
        if (m_is3D)
        {
            put3DData(out);
        }
        else
        {
            put2DData(out);
        }
    }

    string Patch::getLegend(bool bold, bool italic)
    {
        string legend = getPropertyValue("displayname").StringVal();
        stringstream ss;
        ss << "title \"";
        if (italic) {
            ss << "{/:Italic ";
        }
        if (bold) {
            ss << "{/:Bold ";
        }
        ss << legend;
        if (bold) {
            ss << "}";
        }
        if (italic) {
            ss << "}";
        }
        ss << "\"";
        return ss.str();
    }

    void Patch::cleanup(GnuplotOutput*)
    {
    }
    
    void Patch::put2DData(GnuplotOutput* out)
    {
        vector<double> x = getXdata();
        vector<double> y = getYdata();

        double xOffset = getPropertyValue("dataxoffset").Scalar();
        double yOffset = getPropertyValue("datayoffset").Scalar();
        double xScale = getPropertyValue("dataxscale").Scalar();
        double yScale = getPropertyValue("datayscale").Scalar();

        int numFaces = _xcolcount;
        int numPoints = (int)x.size() / numFaces;
        int k = 0;
        for (int i = 0; i < numFaces; ++i)
        {
            for (int j = 0; j < numPoints; ++j, ++k)
            {
                if (k < x.size() && k < y.size() && !std::isnan(x[k]) && !std::isnan(y[k]))
                {
                    out->printf("%g %g\n", x[k] * xScale + xOffset, y[k] * yScale + yOffset);
                }
            }
            out->printf("\n");
        }
            out->printf("e\n");

    }

    void Patch::put3DData(GnuplotOutput* out)
    {
        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> z = getZdata();

        double xOffset = getPropertyValue("dataxoffset").Scalar();
        double yOffset = getPropertyValue("datayoffset").Scalar();
        double zOffset = getPropertyValue("datazoffset").Scalar();
        double xScale = getPropertyValue("dataxscale").Scalar();
        double yScale = getPropertyValue("datayscale").Scalar();
        double zScale = getPropertyValue("datazscale").Scalar();

        int numFaces = _xcolcount;
        int numPoints = (int)x.size() / numFaces;
        int k = 0;
        for (int i = 0; i < numFaces; ++i)
        {
            for (int j = 0; j < numPoints; ++j, ++k)
            {
                if (k < x.size() && k < y.size() && k < z.size() 
                    && !std::isnan(x[k]) && !std::isnan(y[k]) && !std::isnan(z[k]))
                {
                    out->printf("%g %g %g\n", x[k] * xScale + xOffset, y[k] * yScale + yOffset, z[k] * zScale + zOffset);
                }
            }
            out->printf("\n");
        }
            out->printf("e\n");
        }

    Triplot::Triplot()
        :Patch()
    {
        m_type = ObjectType::TRIPLOT;
        m_is3D = false;

        std::vector<Property>::iterator fit = std::find_if(m_ps.begin(),
            m_ps.end(), [](Property p) {return p.getName() == "facecolor"; });
        if (fit != m_ps.end())
            m_ps.erase(fit);

        m_ps.push_back(Property("color", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("linestyle", string("-"), PropertyType::STRING));
        m_ps.push_back(Property("linewidth", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("marker", string(""), PropertyType::STRING));
        m_ps.push_back(Property("markerfacecolor", Color(string("blue")), PropertyType::COLOR));
        m_ps.push_back(Property("markersize", 1.0, PropertyType::DOUBLE));
    }

    void Triplot::init(const LineData& ld)
    {
        std::vector<double> tri = ld.tri;
        std::vector<double> x = ld.x;
        std::vector<double> y = ld.y;
        std::string style = ld.style;

        // all matrices must have values
        if (x.empty() || y.empty() || (x.size() != y.size()))
            throw OML_Error(OML_ERR_PLOT_DIM_NOT_MATCH);

        // check the tri indices
        int M = static_cast<int>(x.size());
        for (size_t i = 0; i < tri.size(); ++i)
        {
            int t = (int)tri[i];
            t -= 1;
            if (t < 0 || t >= M)
            {
                throw OML_Error(OML_ERR_INVALIDINDEX);
            }
        }

        int nTri = ld.triCount;
        // Create the x and y matrices
        std::vector<double> xdata;
        std::vector<double> ydata;
        xdata.reserve(3 * nTri);
        ydata.reserve(3 * nTri);

        for (size_t i = 0; i < tri.size(); ++i)
        {
            int vidx = (int)(tri[i]) - 1;
            xdata.push_back(x[vidx]);
            ydata.push_back(y[vidx]);
        }

        setPropertyValue("xdata", xdata);
        setPropertyValue("ydata", ydata);
        setPropertyValue("zdata", ld.z);

        _xcolcount = nTri;
        _ycolcount = ld.ycolcount;
        _zcolcount = ld.zcolcount;

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        LineStyle ls = LineStyle(ld, colorOrder);
        setPropertyValue("linestyle", ls.m_lineStyle);
        setPropertyValue("color", ls.m_lineColor);
        setPropertyValue("displayname", ls.m_legend);
        setPropertyValue("linewidth", ls.m_lineWidth);
        setPropertyValue("marker", ls.m_markerStyle);
        setPropertyValue("markerfacecolor", ls.m_markerColor);
        setPropertyValue("markersize", ls.m_markerSize);

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string Triplot::getUsingClause()
    {
        std::stringstream ss;
        string lineColor = getPropertyValue("color").ColorString();
        ss << "'-' using 1:2 axis " << m_xaxisRef << m_yaxisRef;
        return ss.str();
    }

    void Triplot::putData(GnuplotOutput* out)
    {
        vector<double> x = getXdata();
        vector<double> y = getYdata();

        double xOffset = getPropertyValue("dataxoffset").Scalar();
        double yOffset = getPropertyValue("datayoffset").Scalar();
        double xScale = getPropertyValue("dataxscale").Scalar();
        double yScale = getPropertyValue("datayscale").Scalar();

        int numFaces = _xcolcount;
        int numPoints = static_cast<int>(std::min(x.size(), y.size()));
        std::vector<double> xP, yP;
        xP.resize(4);
        yP.resize(4);
        for (int i = 0; i < numFaces; ++i)
        {
            bool nextFace = false;
            bool breakLoop = false;
            // check there are no NaNs in all three points
            for (int j = 0; j < 3; ++j)
            {
                int idx = 3 * i + j;
                breakLoop = idx >= numPoints;
                if (std::isnan(x[idx]) || std::isnan(y[idx]))
                {
                    nextFace = true;
                    break;
                }

                xP[j] = x[idx] * xScale + xOffset;
                yP[j] = y[idx] * yScale + yOffset;
            }
            if (breakLoop)
                break;
            if (nextFace)
                continue;
            xP[3] = xP[0];
            yP[3] = yP[0];

            for (int j = 0; j < 4; ++j)
            {
                out->printf("%g %g\n", xP[j], yP[j]);
            }
            out->printf("\n");
        }
        out->printf("e\n");
    }

    string Triplot::getWithClause(int lineId) {
        string linestyle = getPropertyValue("linestyle").StringVal();
        stringstream ss;
        if (linestyle == "") {
            ss << "with points linestyle " << lineId;
        }
        else {
            ss << "with lp linestyle " << lineId;
        }
        return ss.str();
    }

    string Triplot::getLineStyle(int line) {
        string linestyle = getPropertyValue("linestyle").StringVal();
        int dt = 1;
        if (linestyle == "-") {
            dt = 1;
        }
        else if (linestyle == "--") {
            dt = 2;
        }
        else if (linestyle == ":") {
            dt = 3;
        }
        else if (linestyle == "-.") {
            dt = 4;
        }
        else if (linestyle == "-:") {
            dt = 5;
        }
        else if (linestyle == "") {
            dt = 0;
        }
        string lineColor = getPropertyValue("color").ColorString();
        int lineWidth = int(getPropertyValue("linewidth").Scalar());
        string marker = getPropertyValue("marker").StringVal();
        int pt = 0;
        if (marker == "") {
            pt = 0;
        }
        else if (marker == "s") {
            pt = 5;
        }
        else if (marker == "^") {
            pt = 9;
        }
        else if (marker == "v") {
            pt = 11;
        }
        else if (marker == "x") {
            pt = 2;
        }
        else if (marker == "o") {
            pt = 7;
        }
        else if (marker == "d") {
            pt = 13;
        }
        else if (marker == "+") {
            pt = 1;
        }
        else if (marker == "*") {
            pt = 3;
        }
        else if (marker == ".") {
            pt = 0;
        }
        else {
            pt = 6;
        }

        string markerColor = getPropertyValue("markerfacecolor").ColorString();
        int markerSize = int(getPropertyValue("markersize").Scalar());
        
        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line << " dt " << dt
            << " lc rgb '" << lineColor << "' lw " << lineWidth
            << " pt " << pt << " ps " << markerSize
            //<< " tc rgb '" << markerColor << "'"
            << "\n";
        return ss.str();
    }

    Trimesh::Trimesh()
        :Patch(), m_flatShading(true)
    {
        m_type = ObjectType::TRIMESH;
        m_is3D = true;

        std::vector<Property>::iterator fit = std::find_if(m_ps.begin(),
            m_ps.end(), [](Property p) {return p.getName() == "facecolor"; });
        if (fit != m_ps.end())
            m_ps.erase(fit);

        m_ps.push_back(Property("edgecolor", Color("blue"), PropertyType::COLOR));
        m_ps.push_back(Property("cdata", Currency(), PropertyType::MATRIX));
    }

    void Trimesh::init(const LineData& ld)
    {
        // validate inputs
        std::vector<double> tri = ld.tri;
        std::vector<double> x = ld.x;
        std::vector<double> y = ld.y;
        std::vector<double> z = ld.z;
        Currency tmp = ld.cData;
        const hwMatrix* cdata = tmp.Matrix();
        std::string style = ld.style;

        int faceCount = ld.triCount;
        int vertexCount = static_cast<int>(x.size());
        // helper variable to identify the type of cdata
        // 0: no cdata provided
        // 1: scalar value per vertex
        // 2: scalar value per face
        // 3: rgb color per face 
        int cdataType = 0;
        // all matrices must have values
        if (!tri.empty() && !x.empty() && !y.empty() && !z.empty())
        {
            if (x.size() != y.size() || x.size() != z.size())
                throw OML_Error(OML_ERR_PLOT_DIM_NOT_MATCH);
            
            if (!cdata->IsEmpty())
            {
                // cdata M must be either z->M() (color per vertex) or 
                // tri->M() (color per face)
                if (cdata->IsVector())
                {
                    if (cdata->Size() != vertexCount && cdata->Size() != faceCount)
                        throw OML_Error("Error: the number of colors must be equal to the number of triangles or the number of vertices");

                    cdataType = (cdata->Size() == vertexCount) ? 1 : 2;
                }
                else
                {
                    if (cdata->M() != faceCount)
                        throw OML_Error("Error: the number of colors must be equal to the number of triangles or the number of vertices");

                    // cdata N must be 3 (rgb color)
                    if (cdata->N() != 3)
                        throw OML_Error("Error: must be an Mx3 matrix");

                    cdataType = 3;
                }
            }
        }
        else
        {
            throw OML_Error(OML_ERR_PLOT_UNKNOWN_ERROR);
        }

        // check the tri indices
        for (int i = 0; i < static_cast<int>(tri.size()); ++i)
        {
            int t = (int)tri[i];
            t -= 1;
            if (t < 0 || t >= vertexCount)
            {
                throw OML_Error(OML_ERR_INVALIDINDEX);
            }
        }

        std::vector<double> xdata, ydata, zdata;
        xdata.reserve(3 * faceCount);
        ydata.reserve(3 * faceCount);
        zdata.reserve(3 * faceCount);
        
        Currency dataC;
        if (cdataType == 1)
        {
            dataC = Currency(new hwMatrix(3, faceCount, hwMatrix::REAL));
        }
        else if (cdataType == 2)
        {
            dataC = Currency(new hwMatrix(1, faceCount, hwMatrix::REAL));
        }
        else if (cdataType == 3)
        {
            std::vector<int> dims;
            dims.push_back(1);
            dims.push_back(faceCount);
            dims.push_back(3);
            dataC = Currency(new hwMatrixN(dims, hwMatrixN::REAL));
        }

        for (int i = 0; i < faceCount; ++i)
        {
            if (cdataType == 3)
            {
                double r, g, b;
                r = cdata->IsReal() ? (*cdata)(i, 0) : cdata->z(i, 0).Real();
                g = cdata->IsReal() ? (*cdata)(i, 1) : cdata->z(i, 1).Real();
                b = cdata->IsReal() ? (*cdata)(i, 2) : cdata->z(i, 2).Real();
                hwMatrixN* tmp = dataC.GetWritableMatrixN();
                (*tmp)(i) = r;
                (*tmp)(faceCount + i) = g;
                (*tmp)(2 * faceCount + i) = b;
            }
            else if (cdataType == 2)
            {
                double r = cdata->IsReal() ? (*cdata)(i) : cdata->z(i).Real();
                hwMatrix* tmp = dataC.GetWritableMatrix();
                (*tmp)(i) = r;
            }

            // Get the vertices of this face
            for (int j = 0; j < 3; ++j)
            {
                int idx = i * 3 + j;
                int t = (int)tri[idx];
                t -= 1;
                xdata.push_back(x[t]);
                ydata.push_back(y[t]);
                zdata.push_back(z[t]);
                if (cdataType == 1)
                {
                    hwMatrix* tmp = dataC.GetWritableMatrix();
                    (*tmp)(j, i) = cdata->IsReal() ? (*cdata)(t) : cdata->z(t).Real();
                }
            }
        }

        setPropertyValue("xdata", xdata);
        setPropertyValue("ydata", ydata);
        setPropertyValue("zdata", zdata);
        if (cdataType == 0)
        {
            Currency tmp(zdata);
            hwMatrix* mat = tmp.GetWritableMatrix();
            mat->Reshape(3, faceCount);
            setPropertyValue("cdata", tmp);

        }
        else
        {
            setPropertyValue("cdata", dataC);
        }
        
        _xcolcount = faceCount;
        _ycolcount = ld.ycolcount;
        _zcolcount = ld.zcolcount;

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        LineStyle ls = LineStyle(ld, colorOrder);
        setPropertyValue(getColorPropertyName(), ls.m_lineColor);
        setPropertyValue("displayname", ls.m_legend);

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string Trimesh::getUsingClause()
    {
        Currency cDataCur = getPropertyValue("cdata").getCurrency();
        std::stringstream ss;
        if (cDataCur.IsNDMatrix() && m_flatShading)
        {
            ss << "'-' using 1:2:3:(rgb($4,$5,$6)) with lines lc rgb variable";
        }
        else if (m_flatShading)
        {
            ss << "'-' using 1:2:3:4 with lines lc palette";
        }
        else
        {
            string edgeColor = getPropertyValue(getColorPropertyName()).ColorString();
            ss << "'-' using 1:2:3 with lines lc rgb \"" << edgeColor << "\"";
        }

        return ss.str();
    }

    void Trimesh::putData(GnuplotOutput* out)
    {
        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> z = getZdata();

        double xOffset = getPropertyValue("dataxoffset").Scalar();
        double yOffset = getPropertyValue("datayoffset").Scalar();
        double zOffset = getPropertyValue("datazoffset").Scalar();
        double xScale = getPropertyValue("dataxscale").Scalar();
        double yScale = getPropertyValue("datayscale").Scalar();
        double zScale = getPropertyValue("datazscale").Scalar();

        Currency cDataCur = getPropertyValue("cdata").getCurrency();
        int numFaces = _xcolcount;
        int numPoints = static_cast<int>(std::min(x.size(), y.size()));
        numPoints = std::min(numPoints, static_cast<int>(z.size()));
        bool isRGBPerFace = cDataCur.IsNDMatrix() && m_flatShading;

        vector<string> colorStringPerFace = getColorStringPerFace();
        std::vector<double> xP, yP, zP;
        xP.resize(4);
        yP.resize(4);
        zP.resize(4);
        stringstream ss;
        for (int i = 0; i < numFaces; ++i)
        {
            bool nextFace = false;
            bool breakLoop = false;
            // check there are no NaNs in all three points
            for (int j = 0; j < 3; ++j)
            {
                int idx = 3 * i + j;
                breakLoop = idx >= numPoints;
                if (std::isnan(x[idx]) || std::isnan(y[idx]) || std::isnan(z[idx]))
                {
                    nextFace = true;
                    break;
                }

                xP[j] = x[idx] * xScale + xOffset;
                yP[j] = y[idx] * yScale + yOffset;
                zP[j] = z[idx] * zScale + zOffset;
            }
            if (breakLoop)
                break;
            if (nextFace)
                continue;
            xP[3] = xP[0];
            yP[3] = yP[0];
            zP[3] = zP[0];

            if (isRGBPerFace || m_flatShading)
            {
                string color = colorStringPerFace[i];
                for (int j = 0; j < 4; ++j)
                {
                    ss << xP[j] << " " << yP[j] << " " << zP[j] << " " << color.c_str() << "\n";
                }
            }
            else
            {
                for (int j = 0; j < 4; ++j)
                {
                    ss << xP[j] << " " << yP[j] << " " << zP[j] << "\n";
                }
            }
            ss << "\n\n";
        }
        out->printf(ss.str() + "e\n");
    }

    string Trimesh::getLineStyle(int)
    {
        return std::string("rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b)\n");
    }

    VALUETYPE Trimesh::getPropertyValue(const string& name)
    {
        if (name == getColorPropertyName() && m_flatShading)
            return VALUETYPE("flat");
        return Patch::getPropertyValue(name);
    }

    bool Trimesh::setPropertyValue(const string& name, VALUETYPE value)
    {
        if (name == getColorPropertyName())
        {
            if (value.isCurrency())
            {
                Currency c = value.getCurrency();
                if (c.IsString())
                {
                    string val = c.StringVal();
                    if (val == "flat")
                    {
                        m_flatShading = true;
                        return true;
                    }
                    else if (val == "interp")
                    {
                        CoreMain::getInstance()->AddWarningString("value 'interp' for [" + name + "] is not supported in Open Matrix");
                        return true;
                    }
                    else
                    {
                        Patch::setPropertyValue(name, value);
                        m_flatShading = false;
                        return true;
                    }
                }
            }
        }
        return Patch::setPropertyValue(name, value);
    }

    bool Trimesh::getMinMaxZ(double& min, double& max)
    {
        bool ret = false;
        if (m_flatShading)
        {
            Currency cDataCur = getPropertyValue("cdata").getCurrency();
            if (cDataCur.IsMatrix() && !cDataCur.IsEmpty())
            {
                const hwMatrix* mat = cDataCur.Matrix();
                bool isReal = mat->IsReal();
                min = isReal ? (*mat)(0) : mat->z(0).Real();
                max = min;
                if ((_xcolcount != 1 && mat->IsVector()) || (_xcolcount == 1 && mat->Size() == 1))
                {
                    // cdata is a scalar value per face
                    for (int i = 0; i < mat->Size(); i++) {
                        double zVal = isReal ? (*mat)(i) : mat->z(i).Real();
                        if (zVal < min)
                            min = zVal;
                        if (zVal > max)
                            max = zVal;
                    }
                    ret = true;
                }
                else
                {
                    int cN = mat->N();
                    for (int i = 0; i < _xcolcount; ++i)
                    {
                        if (i >= cN)
                            break;

                        double zVal = isReal? (*mat)(0, i) : mat->z(0, i).Real();
                        if (zVal < min)
                            min = zVal;
                        if (zVal > max)
                            max = zVal;
                    }
                }
                ret = true;
            }
        }
        return ret;
    }

    vector<string> Trimesh::getColorStringPerFace()
    {
        Currency cDataCur = getPropertyValue("cdata").getCurrency();
        int numFaces = _xcolcount;
        bool isRGBPerFace = cDataCur.IsNDMatrix() && m_flatShading;
        bool singleFace = numFaces == 1;

        std::vector<std::string> colorStringPerFace;
        colorStringPerFace.reserve(numFaces);
        if (isRGBPerFace)
        {
            // rgb per face
            const hwMatrixN* cdata = cDataCur.MatrixN();
            double r, g, b;
            std::vector<int> dims = cdata->Dimensions();
            int numFacesCData = 0;
            if (dims.size() > 1)
                numFacesCData = dims[1];

            bool cDataIsReal = cdata->IsReal();
            std::vector<std::string> colorsPerFace;
            // set the same color value to all vertices of a face
            for (int i = 0; i < numFaces; ++i)
            {
                r = g = b = 0.0;
                if ((2 * numFacesCData + i) < cdata->Size())
                {
                    r = cDataIsReal ? (*cdata)(i) : cdata->z(i).Real();
                    g = cDataIsReal ? (*cdata)(numFacesCData + i) : cdata->z(numFacesCData + i).Real();
                    b = cDataIsReal ? (*cdata)(2 * numFacesCData + i) : cdata->z(2 * numFacesCData + i).Real();
                }
                if (r <= 1 && g <= 1 && b <= 1)
                {
                    r *= 255;
                    g *= 255;
                    b *= 255;
                }

                char buffer[100];
                sprintf(buffer, "%d %d %d", (int)r, (int)g, (int)b);
                colorStringPerFace.push_back(std::string(buffer));
            }
        }
        else if (m_flatShading)
        {
            if (cDataCur.IsMatrix() && cDataCur.Matrix()->Size() > 0)
            {
                const hwMatrix* cdata = cDataCur.Matrix();
                bool cDataIsReal = cdata->IsReal();

                if ((!singleFace && cdata->IsVector()) || (singleFace && cdata->Size() == 1))
                {
                    int S = cdata->Size();
                    // cdata is a scalar value per face
                    for (int i = 0; i < numFaces; ++i)
                    {
                        double val = 0;
                        if (i < S)
                            val = cDataIsReal ? (*cdata)(i) : cdata->z(i).Real();

                        char buffer[100];
                        sprintf(buffer, "%g", val);
                        colorStringPerFace.push_back(std::string(buffer));
                    }
                }
                else
                {
                    int cN = cdata->N();
                    for (int i = 0; i < numFaces; ++i)
                    {
                        double val = 0;
                        if (i < cN)
                            val = cDataIsReal ? (*cdata)(0, i) : cdata->z(0, i).Real();
                        char buffer[100];
                        sprintf(buffer, "%g", val);
                        colorStringPerFace.push_back(std::string(buffer));
                    }
                }
            }
            else
            {
                colorStringPerFace.insert(colorStringPerFace.end(), numFaces, "0");
            }
        }
        return colorStringPerFace;
    }

    Trisurf::Trisurf()
        :Trimesh()
    {
        m_type = ObjectType::TRISURF;

        std::vector<Property>::iterator fit = std::find_if(m_ps.begin(),
            m_ps.end(), [](Property p) {return p.getName() == "edgecolor"; });
        if (fit != m_ps.end())
            m_ps.erase(fit);

        m_ps.push_back(Property("facecolor", Color("blue"), PropertyType::COLOR));
    }

    string Trisurf::getUsingClause()
    {
        Currency cDataCur = getPropertyValue("cdata").getCurrency();
        std::stringstream ss;
        if (cDataCur.IsNDMatrix() && m_flatShading)
        {
            ss << "'-' using 1:2:3:(rgb($4,$5,$6)) with pm3d fc rgb variable";
        }
        else if (m_flatShading)
        {
            ss << "'-' using 1:2:3:4 with pm3d lc palette";
        }
        else
        {
            string color = getPropertyValue(getColorPropertyName()).ColorString();
            ss << "'-' using 1:2:3 with polygons fc rgb \"" << color << "\"";
        }

        return ss.str();
    }

    void Trisurf::putData(GnuplotOutput* out)
    {
        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> z = getZdata();

        double xOffset = getPropertyValue("dataxoffset").Scalar();
        double yOffset = getPropertyValue("datayoffset").Scalar();
        double zOffset = getPropertyValue("datazoffset").Scalar();
        double xScale = getPropertyValue("dataxscale").Scalar();
        double yScale = getPropertyValue("datayscale").Scalar();
        double zScale = getPropertyValue("datazscale").Scalar();

        Currency cDataCur = getPropertyValue("cdata").getCurrency();
        int numFaces = _xcolcount;
        int numPoints = static_cast<int>(std::min(x.size(), y.size()));
        numPoints = std::min(numPoints, static_cast<int>(z.size()));
        bool isRGBPerFace = cDataCur.IsNDMatrix() && m_flatShading;
        
        vector<string> colorStringPerFace = getColorStringPerFace();
        std::vector<double> xP, yP, zP;
        xP.resize(3);
        yP.resize(3);
        zP.resize(3);
        stringstream ss;
        for (int i = 0; i < numFaces; ++i)
        {
            bool nextFace = false;
            bool breakLoop = false;
            // check there are no NaNs in all three points
            for (int j = 0; j < 3; ++j)
            {
                int idx = 3 * i + j;
                breakLoop = idx >= numPoints;
                if (std::isnan(x[idx]) || std::isnan(y[idx]) || std::isnan(z[idx]))
                {
                    nextFace = true;
                    break;
                }

                xP[j] = x[idx] * xScale + xOffset;
                yP[j] = y[idx] * yScale + yOffset;
                zP[j] = z[idx] * zScale + zOffset;
            }
            if (breakLoop)
                break;
            if (nextFace)
                continue;

            if (isRGBPerFace || m_flatShading)
            {
                string color = colorStringPerFace[i];
                for (int j = 0; j < 2; ++j)
                {
                    ss << xP[j] << " " << yP[j] << " " << zP[j] << " " << color.c_str() << "\n";
                }
                ss << "\n";
                ss << xP[2] << " " << yP[2] << " " << zP[2] << " " << color.c_str() << "\n";
                ss << xP[2] << " " << yP[2] << " " << zP[2] << " " << color.c_str() << "\n";
                ss << "\n\n";
            }
            else
            {
                for (int j = 0; j < 3; ++j)
                {
                    ss << xP[j] << " " << yP[j] << " " << zP[j] << "\n";
                }
                ss << "\n";
            }
        }
        out->printf(ss.str() + "e\n");
    }

    string Trisurf::getLineStyle(int line)
    {
        return Trimesh::getLineStyle(line) + Patch::getLineStyle(line);
    }
    
    Legend::Legend()
        :Object(), m_legendLocation("top right")
    {
        m_type = ObjectType::LEGEND;
        double handle = m_handlePool->allocHandle();
        m_ps.push_back(Property("type", string("legend"), PropertyType::STRING));
        m_ps.push_back(Property("handle", handle, PropertyType::DOUBLE));
        m_ps.push_back(Property("bordercolor", Color(string("black")), PropertyType::COLOR));
        m_ps.push_back(Property("borderwidth", 1, PropertyType::DOUBLE));
        m_ps.push_back(Property("fontname", string("arial"), PropertyType::STRING));
        m_ps.push_back(Property("fontsize", double(8), PropertyType::DOUBLE));
        m_ps.push_back(Property("fontweight", string("normal"), PropertyType::STRING));
        m_ps.push_back(Property("fontangle", string("regular"), PropertyType::STRING));
        m_ps.push_back(Property("location",string("northeast"), PropertyType::STRING));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));

        m_objectMap[handle] = this;

        // not yet supported properties
        m_ps.push_back(Property("units", string("units"), PropertyType::UNSUPPORTED));
    }

    void Legend::update(GnuplotOutput* out)
    {
        if (getVisible())
        {
            string bc = getPropertyValue("bordercolor").ColorString();
            double bw = getPropertyValue("borderwidth").Scalar();
            string fn = getPropertyValue("fontname").StringVal();
            double fs = getPropertyValue("fontsize").Scalar();
            string fw = getPropertyValue("fontweight").StringVal();
            string fa = getPropertyValue("fontangle").StringVal();
            
            stringstream ss;
            ss << "set key on " << m_legendLocation;
            ss << " box linestyle 1000 lc rgb '" << bc << "'";
            ss << " linewidth " << bw;
            ss << " font \"" << fn << "," << fs << "\"\n";
            out->printf(ss.str());
        }
        else
        {
            out->printf("set key off\n");
        }
    }

    void Legend::setLegendLocation(const std::string& legendLoc)
    {
        Object* p = getParentObject();
        if (p->getObjectType() != ObjectType::AXES)
            return;

        if (static_cast<Axes*>(p)->is3DPlot()) {
            if (legendLoc == "northeast") {
                m_legendLocation = "outside top right vertical";
            }
            else if (legendLoc == "northwest") {
                m_legendLocation = "outside top left vertical";
            }
            else if (legendLoc == "southeast") {
                m_legendLocation = "outside bottom right vertical";
            }
            else if (legendLoc == "southwest") {
                m_legendLocation = "outside bottom left vertical";
            }
            else {
                throw OML_Error(OML_Error(OML_ERR_BAD_STRING).GetErrorMessage() + " for [location]");
            }
        }
        else {
            if (legendLoc == "north") {
                m_legendLocation = "inside top center vertical";
            }
            else if (legendLoc == "south") {
                m_legendLocation = "inside bottom center vertical";
            }
            else if (legendLoc == "east") {
                m_legendLocation = "inside center right vertical";
            }
            else if (legendLoc == "west") {
                m_legendLocation = "inside center left vertical";
            }
            else if (legendLoc == "northeast") {
                m_legendLocation = "inside top right vertical";
            }
            else if (legendLoc == "northwest") {
                m_legendLocation = "inside top left vertical";
            }
            else if (legendLoc == "southeast") {
                m_legendLocation = "inside bottom right vertical";
            }
            else if (legendLoc == "southwest") {
                m_legendLocation = "inside bottom left vertical";
            }
            else if (legendLoc == "northoutside") {
                m_legendLocation = "outside top center horizontal";
            }
            else if (legendLoc == "southoutside") {
                m_legendLocation = "outside bottom center horizontal";
            }
            else if (legendLoc == "eastoutside") {
                m_legendLocation = "outside center right vertical";
            }
            else if (legendLoc == "westoutside") {
                m_legendLocation = "outside center left vertical";
            }
            else {
                throw OML_Error(OML_Error(OML_ERR_BAD_STRING).GetErrorMessage() + " for [location]");
            }
        }
    }

    void Legend::setLegendLocation(const std::vector<double>& legendLoc)
    {
        Object* p = getParentObject();
        if (p->getObjectType() != ObjectType::AXES)
            return;

        if (static_cast<Axes*>(p)->is3DPlot())
            throw OML_Error(OML_Error(OML_ERR_BAD_STRING).GetErrorMessage() + " for [location]");

        if (legendLoc.size() != 2)
            throw OML_Error(OML_ERR_VECTOR2);

        double xpos = legendLoc[0];
        double ypos = legendLoc[1];
        if (xpos >= 0 && xpos <= 1 && ypos >= 0 && ypos <= 1) {
            stringstream ss;
            ss << "at graph " << xpos << "," << (1 - ypos) << " top left vertical";
            m_legendLocation = ss.str();
        }
        else {
            throw OML_Error(OML_ERR_PLOT_NEED_NORM_DATA);
        }
    }

    std::string Legend::getLegendLocation() const
    {
        return m_legendLocation;
    }

    void Legend::setVisible(bool visible)
    {
        setPropertyValue("visible", visible ? string("on") : string("off"));
    }

    bool Legend::getVisible()
    {
        return getPropertyValue("visible").StringVal() == "on";
    }

    bool Legend::setPropertyValue(const string& name, VALUETYPE value)
    {
        if (name == "location") {
            if (!value.isCurrency())
                throw OML_Error(OML_ERR_OPTIONVAL);
            Currency c = value.getCurrency();
            if (c.IsString())
                setLegendLocation(c.StringVal());
            else if (c.IsVector())
                setLegendLocation(c.Vector());
            else
                throw OML_Error(OML_ERR_OPTIONVAL);
            Object::setPropertyValue(name, value);
            return true;
        }
        return Object::setPropertyValue(name, value);
    }

    Colorbar::Colorbar()
        :Object(), m_colorbarLocation("default")
    {
        m_type = ObjectType::COLORBAR;
        double handle = m_handlePool->allocHandle();
        m_ps.push_back(Property("type", string("colorbar"), PropertyType::STRING));
        m_ps.push_back(Property("handle", handle, PropertyType::DOUBLE));
        m_ps.push_back(Property("clim", Currency(), PropertyType::VEC_DOUBLE));
        m_ps.push_back(Property("colorbarscale", string("linear"), PropertyType::STRING));
        vector<double> in;
        in.push_back(9);
        m_ps.push_back(Property("colorlevels", in, PropertyType::VEC_DOUBLE));
        hwMatrix* map = new hwMatrix(4, 3, hwMatrix::DataType::REAL);
        (*map)(0, 0) = 0; (*map)(0, 1) = 0; (*map)(0, 2) = 255;
        (*map)(1, 0) = 0; (*map)(1, 1) = 255; (*map)(1, 2) = 0;
        (*map)(2, 0) = 255; (*map)(2, 1) = 255; (*map)(2, 2) = 0;
        (*map)(3, 0) = 255; (*map)(3, 1) = 0; (*map)(3, 2) = 0;
        m_ps.push_back(Property("colormap", Currency(map), PropertyType::MATRIX));
        m_ps.push_back(Property("fontname", string("arial"), PropertyType::STRING));
        m_ps.push_back(Property("fontsize", double(8), PropertyType::DOUBLE));
        m_ps.push_back(Property("location", "northwest", PropertyType::STRING));
        m_ps.push_back(Property("numericformat", "fixed", PropertyType::STRING));
        m_ps.push_back(Property("numericprecision", 2, PropertyType::DOUBLE));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));

        m_objectMap[handle] = this;

        // not yet supported properties
        m_ps.push_back(Property("units", string("units"), PropertyType::UNSUPPORTED));

        setPropertyValue("location", "northwest");
    }

    void Colorbar::update(GnuplotOutput* out)
    {
        if (getVisible()) {
            double x = 0, y = 0, w = 1, h = 1;
            Object* obj = getParentObject();
            if (obj->isAxes()) {
                static_cast<Axes*>(obj)->getAxesPosition(x, y, w, h);
        }
            string loc;
            out->printf("set colorbox vertical user origin ");
            if (m_colorbarLocation == "northeast") {
                out->printf(" %g, %g size %g, %g\n", (x + w) - 0.12 * w, (y + h) - 0.5 * h, 0.02 * w, 0.48 * h);
            }
            else if (m_colorbarLocation == "northwest") {
                out->printf(" %g, %g size %g, %g\n", x + 0.01*w, (y + h) - 0.5 * h, 0.02 * w, 0.48 * h);
            }
            else if (m_colorbarLocation == "southeast") {
                out->printf(" %g, %g size %g, %g\n", (x + w) - 0.12 * w, y + 0.01 * h, 0.02 * w, 0.48 * h);
            }
            else if (m_colorbarLocation == "southwest") {
                out->printf(" %g, %g size %g, %g\n", x + 0.01 * w, y + 0.01 * h, 0.02 * w, 0.48 * h);
            }
        }
        else {
            out->printf("unset colorbox\n");
        }
        
        // palette
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
                    colorStrs.push_back("\"" + cl.getString() + "\"");
                }
            }
        }

        if (colorStrs.empty()) {
            colorStrs.push_back(string("\"blue\""));
            colorStrs.push_back(string("\"green\""));
            colorStrs.push_back(string("\"yellow\""));
            colorStrs.push_back(string("\"red\""));
        }

        out->printf("set palette defined (");
        int sizeCM = int(colorStrs.size());
        for (int i = 0; i < sizeCM; ++i) {
            out->printf("%d %s", i, colorStrs[i].c_str());
            if (i < (sizeCM - 1))
                out->printf(",");
        }
        out->printf(")\n");

        // tics
        vector<double> colorlevels;
        Currency c = getPropertyValue("colorlevels").getCurrency();
        if (c.IsVector() && c.Vector().size()>1) {
            vector<double> cbl = c.Vector();
            out->printf("set cbrange [%g:%g]\n", cbl.front(), cbl.back());
            out->printf("set cbtics (");
            for (int i = 0; i < (int)cbl.size(); ++i) {
                out->printf("%g", cbl[i]);
                if (i < ((int)cbl.size() - 1))
                    out->printf(",");
            }
            out->printf(")\n");
        }
        else {
            int numLevels = 10;
            if (c.IsScalar())
                numLevels = (int)c.Scalar();
            else if (c.IsVector() && !c.Vector().empty())
                numLevels = (int)c.Vector().front();

            Object* p = getParentObject();

            vector<double> cr;
            Currency clim = getPropertyValue("clim").getCurrency();
            if (clim.IsVector() && clim.Vector().size() == 2) {
                cr = clim.Vector();
            }
            else if (p && p->isAxes()) {
                cr = static_cast<Axes*>(p)->getColorbarRange();
            }
            else {
                cr.push_back(0);
                cr.push_back(1);
            }
            if (cr.size() == 2) {
                out->printf("set cbrange [%g:%g]\n", cr[0], cr[1]);
                out->printf("set cbtics %g\n", (cr[1] - cr[0]) / numLevels);
            }
            else {
                out->printf("unset cbrange\n");
            }
        }

        // tics format
        string cbf = getPropertyValue("numericformat").StringVal();
        int cbnp = int(getPropertyValue("numericprecision").Scalar());
        string f = "g";
        if (cbf == "scientific")
            f = "e";
        else if (cbf == "fixed")
            f = "f";
        char numF[24];
        sprintf(numF, "%s.%d%s", "%", cbnp, f.c_str());
        string fn = getPropertyValue("fontname").StringVal();
        double fs = getPropertyValue("fontsize").Scalar();
        out->printf("set cbtics format \"%s\" font \"%s,%g\" \n", numF, fn.c_str(), fs);

        string colorbarScale = getPropertyValue("colorbarscale").StringVal();
        if (colorbarScale == "log")
            out->printf("set log cb\nset cbtics 2\n");
        else
            out->printf("unset log cb\n");

    }

    void Colorbar::setVisible(bool visible)
    {
        setPropertyValue("visible", visible ? string("on") : string("off"));
    }

    bool Colorbar::getVisible()
    {
        return getPropertyValue("visible").StringVal() == "on";
    }

    bool Colorbar::setPropertyValue(const string& name, VALUETYPE value)
    {
        if(name == "location") {
            if (!value.isCurrency() || !value.getCurrency().IsString())  
                throw OML_Error(OML_ERR_OPTIONVAL);
            
            string loc = value.getCurrency().StringVal();
            if (loc == "northeast"|| loc == "northwest" ||
                loc == "southeast" || loc == "southwest") {
                m_colorbarLocation = loc;
            }
            else {
                throw OML_Error(OML_Error(OML_ERR_BAD_STRING).GetErrorMessage() + " for [location]");
            }
            Object::setPropertyValue(name, value);
            return true;
        }
        return Object::setPropertyValue(name, value);
    }

    Stem3::Stem3()
        :Line3()
    {
        m_type = ObjectType::STEM3;
    }

    string Stem3::getWithClause(int lineId)
    {
        stringstream ss;
        ss << "with impulses linestyle " << lineId;
        return ss.str();
    }

    HggroupVector::HggroupVector()
        :Drawable()
    {
        m_type = ObjectType::HGGROUPVECTOR;

        std::vector<Property>::iterator fit = std::find_if(m_ps.begin(),
            m_ps.end(), [](Property p) {return p.getName() == "zdata"; });
        if (fit != m_ps.end())
            m_ps.erase(fit);

        m_ps.push_back(Property("type", string("hggroup"), PropertyType::STRING));
        m_ps.push_back(Property("color", string("blue"), PropertyType::COLOR));
        m_ps.push_back(Property("dataxoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("dataxscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayoffset", 0.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("datayscale", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("linestyle", string("-"), PropertyType::STRING));
        m_ps.push_back(Property("linewidth", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("marker", string("o"), PropertyType::STRING));
        m_ps.push_back(Property("markerfacecolor", string("blue"), PropertyType::COLOR));
        m_ps.push_back(Property("markersize", 1.0, PropertyType::DOUBLE));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
        vector<double> d;
        m_ps.push_back(Property("udata", d, PropertyType::VEC_DOUBLE));
        m_ps.push_back(Property("vdata", d, PropertyType::VEC_DOUBLE));
        m_ps.push_back(Property("showarrowhead", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("autoscale", string("on"), PropertyType::STRING));
        m_ps.push_back(Property("autoscalefactor", 0.9, PropertyType::DOUBLE));
        m_ps.push_back(Property("maxheadsize", 0.3, PropertyType::DOUBLE));

        // not yet supported properties
        m_ps.push_back(Property("units", string("units"), PropertyType::UNSUPPORTED));
        
    }
    void HggroupVector::init(const LineData& ld)
    {
        setPropertyValue("xdata", ld.x);
        setPropertyValue("ydata", ld.y);
        setPropertyValue("udata", ld.u);
        setPropertyValue("vdata", ld.v);

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        // update line-specific properties
        LineStyle ls = LineStyle(ld, colorOrder);
        if (ls.m_lineStyle.empty())
            ls.m_lineStyle = "-";
        setPropertyValue("linestyle", ls.m_lineStyle);
        setPropertyValue("color", ls.m_lineColor);
        setPropertyValue("linewidth", ls.m_lineWidth);
        setPropertyValue("marker", ls.m_markerStyle);
        setPropertyValue("markerfacecolor", ls.m_markerColor);
        setPropertyValue("markersize", ls.m_markerSize);
        setPropertyValue("displayname", ls.m_legend);
        if (!ls.m_markerStyle.empty())
            setPropertyValue("showarrowhead", "off");
        if (ld.x.size()==1)
            setPropertyValue("autoscalefactor", 1.0);

        for (int i = 0; i < ld.properties.size(); i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string HggroupVector::getUsingClause()
    {
        return "'-' using 1:2:3:4 ";
    }

    string HggroupVector::getWithClause(int line)
    {
        stringstream ss;

        double headSize = getPropertyValue("maxheadsize").Scalar();
        string sa = getPropertyValue("showarrowhead").StringVal();
        if (sa == "on") {
            ss << "with vectors filled head linestyle " << line;
            ss << " size graph " << headSize/10 << ",20,60";
        }
        else {
            ss << "with vectors nohead linestyle " << line;
        }

        string marker = getPropertyValue("marker").StringVal();
        if (!marker.empty()) {
            int pt = 0;
            if (marker == "") {
                pt = 0;
            }
            else if (marker == "s") {
                pt = 5;
            }
            else if (marker == "^") {
                pt = 9;
            }
            else if (marker == "v") {
                pt = 11;
            }
            else if (marker == "x") {
                pt = 2;
            }
            else if (marker == "o") {
                pt = 7;
            }
            else if (marker == "d") {
                pt = 13;
            }
            else if (marker == "+") {
                pt = 1;
            }
            else if (marker == "*") {
                pt = 3;
            }
            else if (marker == ".") {
                pt = 0;
            }
            else {
                pt = 6;
            }
            int markerSize = int(getPropertyValue("markersize").Scalar());
            string markerColor = getPropertyValue("markerfacecolor").ColorString();

            ss << ", '-' using 1:2 with points lc rgb '" << markerColor << "' ";
            ss<< " pt "<<pt<<" ps "<<markerSize<<" notitle";

        }
        return ss.str();
    }

    void HggroupVector::putData(GnuplotOutput* out)
    {
        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> u = getUdata();
        vector<double> v = getVdata();

        double xOffset = getPropertyValue("dataxoffset").Scalar();
        double yOffset = getPropertyValue("datayoffset").Scalar();
        double xScale = getPropertyValue("dataxscale").Scalar();
        double yScale = getPropertyValue("datayscale").Scalar();

        string autoscale = getPropertyValue("autoscale").StringVal();
        double autoscaleFactor = getPropertyValue("autoscalefactor").Scalar();
        double arrowScale = 1;
        if (autoscale == "on") {
            arrowScale = autoscaleFactor == 0 ? 1.0 : autoscaleFactor;
        }

        size_t count = min(x.size(), y.size());
        count = min(count, u.size());
        count = min(count, v.size());

        double maxMag = 0;
        // note - uPoints and vPoints matrix dimensions have been validated earlier
        for (int i = 0; i < count; ++i)
        {
            double mag = std::sqrt(u[i] * u[i] + v[i] * v[i]);
            if (maxMag < mag)
                maxMag = mag;
        }

        // first non-zero x-y vector
        double xyVectorMag = 0;
        for (int i = 0; i < (count - 1); ++i)
        {
            double xVal0 = x[i];
            double xVal1 = x[i + 1];
            double yVal0 = y[i];
            double yVal1 = y[i + 1];
            double mag = std::sqrt((xVal0 - xVal1) * (xVal0 - xVal1) + (yVal0 - yVal1) * (yVal0 - yVal1));
            if (xyVectorMag < mag)
            {
                xyVectorMag = mag;
                break;
            }
        }
        // arrow scale ratio
        if (xyVectorMag != 0 && maxMag != 0)
            arrowScale = xyVectorMag == 0 ? arrowScale : arrowScale * xyVectorMag / maxMag;

        stringstream ss;
        for (int j = 0; j < count; j++) {
            double xd = x[j] * xScale + xOffset;
            double yd = y[j] * yScale + yOffset;
            out->printf("%g %g %g %g\n", xd, yd, u[j] * arrowScale, v[j] * arrowScale);
            ss << xd << " " << yd << "\n";
        }
        // end of data
        out->printf("e\n");

        string marker = getPropertyValue("marker").StringVal();
        if (!marker.empty()) {
            out->printf(ss.str());
            out->printf("e\n");
        }
    }

    string HggroupVector::getLegend(bool bold, bool italic)
    {
        string legend = getPropertyValue("displayname").StringVal();
        stringstream ss;
        ss << "title \"";
        if (italic) {
            ss << "{/:Italic ";
        }
        if (bold) {
            ss << "{/:Bold ";
        }
        ss << legend;
        if (bold) {
            ss << "}";
        }
        if (italic) {
            ss << "}";
        }
        ss << "\"";
        return ss.str();
    }
    string HggroupVector::getLineStyle(int line)
    {
        string linestyle = getPropertyValue("linestyle").StringVal();
        int dt = 1;
        if (linestyle == "-") {
            dt = 1;
        }
        else if (linestyle == "--") {
            dt = 2;
        }
        else if (linestyle == ":") {
            dt = 3;
        }
        else if (linestyle == "-.") {
            dt = 4;
        }
        else if (linestyle == "-:") {
            dt = 5;
        }
        else if (linestyle == "") {
            dt = 0;
        }

        string lineColor = getPropertyValue("color").ColorString();
        int lineWidth = int(getPropertyValue("linewidth").Scalar());

        stringstream ss;
        ss << "unset style line " << line << "\n";
        ss << "set style line " << line << " dt " << dt
            << " lc rgb '" << lineColor << "' lw " << lineWidth<<"\n";
        return ss.str();
    }
    void HggroupVector::cleanup(GnuplotOutput*)
    {
    }
    vector<double> HggroupVector::getUdata()
    {
        Currency val = getPropertyValue("udata").getCurrency();
        if (val.IsVector()) {
            return val.Vector();
        }
        else if (val.IsScalar()) {
            vector<double> ret;
            ret.push_back(val.Scalar());
            return ret;
        }
        return vector<double>();
    }
    vector<double> HggroupVector::getVdata()
    {
        Currency val = getPropertyValue("vdata").getCurrency();
        if (val.IsVector()) {
            return val.Vector();
        }
        else if (val.IsScalar()) {
            vector<double> ret;
            ret.push_back(val.Scalar());
            return ret;
        }
        return vector<double>();
    }

    HggroupBar3D::HggroupBar3D()
        :m_numRows(1)
    {
        m_type = ObjectType::HGGROUPBAR3;

        m_ps.push_back(Property("type", string("hggroup"), PropertyType::STRING));
        m_ps.push_back(Property("color", Currency(), PropertyType::MATRIX));
        m_ps.push_back(Property("visible", string("on"), PropertyType::STRING));
    }

    void HggroupBar3D::init(const LineData& ld)
    {
        Drawable::init(ld);

        Currency colorOrder;
        Object* p = getParentObject();
        if (p && p->isAxes())
            colorOrder = p->getPropertyValue("colororder").getCurrency();

        // Calculate number of rows
        int dataSize = static_cast<int>(ld.z.size());
        if (_zcolcount == 0)
            throw OML_Error(OML_ERR_PLOT_UNKNOWN_ERROR);
        m_numRows = dataSize / _zcolcount;


        const hwMatrix* co = colorOrder.Matrix();
        int numColors = co->M();
        hwMatrix* map = new hwMatrix(m_numRows, 3, hwMatrix::DataType::REAL);
        for (int i = 0; i < m_numRows; ++i)
        {
            if (co->Size() > 0 && co->N() == 3)
            {
                (*map)(i, 0) = (*co)(i % numColors, 0);
                (*map)(i, 1) = (*co)(i % numColors, 0);
                (*map)(i, 2) = (*co)(i % numColors, 0);
            }
            else
            {
                std::vector<double> col =LineStyle::ColorAtIndex(ld.index + i);
                (*map)(i, 0) = col[0];
                (*map)(i, 1) = col[1];
                (*map)(i, 2) = col[2];
            }
        }
        // update line-specific properties
        LineStyle ls = LineStyle(ld, colorOrder);
        if (!ld.style.empty())
        {
            Color cl =  Color(ld.style);
            setPropertyValue("color", cl);
        }
        else
        {
            setPropertyValue("color", Currency(map));
        }
        setPropertyValue("displayname", ls.m_legend);

        size_t propSize = ld.properties.size();
        for (int i = 0; i < propSize; i++) {
            setPropertyValue(ld.properties[i], ld.values[i]);
        }
    }

    string HggroupBar3D::getUsingClause()
    {
        string legend = getPropertyValue("displayname").StringVal();
        bool legendBold = false;
        bool legendItalic = false;
        Object* p = getParentObject();
        Axes* pAxes = p && p->getObjectType() == ObjectType::AXES ? static_cast<Axes*>(p) : nullptr;
        if (pAxes) {
            Object* lg = getObject(pAxes->getLegendHandle());
            if (lg && lg->getObjectType() == ObjectType::LEGEND) {
                legendBold = lg->getPropertyValue("fontweight").StringVal() == "bold";
                legendItalic = lg->getPropertyValue("fontangle").StringVal() == "italic";
            }
        }

        stringstream ss;
        ss << "'-' using 1:2:3:4:5 with boxes fc rgb variable title \"";
        if (legendItalic)   ss << "{/:Italic ";
        if (legendBold)     ss << "{/:Bold ";
        ss << legend;
        if (legendBold)     ss << "}";
        if (legendItalic)   ss << "}";
        ss << "\"";
        return ss.str();
    }
    string HggroupBar3D::getWithClause(int)
    {
        return string();
    }
    string HggroupBar3D::getLegend(bool bold, bool italic)
    {
        return string();
    }
    string HggroupBar3D::getLineStyle(int)
    {
        return string();
    }
    void HggroupBar3D::cleanup(GnuplotOutput*)
    {
    }
    void HggroupBar3D::putData(GnuplotOutput* out)
    {
        Object* p = getParentObject();
        if (!p)
            return;

        vector<double> x = getXdata();
        vector<double> y = getYdata();
        vector<double> z = getZdata();
        int dataSize = static_cast<int>(z.size());
        if (_zcolcount == 0)
            throw OML_Error(OML_ERR_PLOT_UNKNOWN_ERROR);
        // set numrows in case z data changed
        m_numRows = dataSize / _zcolcount;

        Currency colorCur = getPropertyValue("color").getCurrency();
        const hwMatrix* colors = colorCur.IsMatrix() ? colorCur.Matrix() : nullptr;
        
        int zi = 0;
        int xSize = std::min(static_cast<int>(x.size()), _zcolcount);
        int ySize = std::min(static_cast<int>(y.size()), m_numRows);
        double width = 0.4;
        if (xSize> 1)
            width = (x[1] - x[0]) / 2;
        for (int i = 0; i < ySize; ++i)
        {
            int rgb = 255;
            if (colors && i < colors->M())
            {
                double r = (*colors)(i, 0);
                double g = (*colors)(i, 1);
                double b = (*colors)(i, 2);
                if (r <= 1 && g <= 1 && b <= 1)
                {
                    r *= 255;
                    g *= 255;
                    b *= 255;
                }
                rgb = (int)r;
                rgb = (rgb << 8) + (int)g;
                rgb = (rgb << 8) + (int)b;
            }
            double tmpY = y[i];
            for (int j = 0; j < xSize; ++j)
            {
                int zIdx = i * _zcolcount + j;
                out->printf("%g %g %g %g %d\n", x[j], tmpY, z[zIdx], width, rgb);
                ++zi;
            }
        }
        // end of data
        out->printf("e\n");
    }

    bool HggroupBar3D::setPropertyValue(const string& name, VALUETYPE value)
    {
        if (name == "color")
        {
            Property& p = getProperty(name);
            if (value.isCurrency() && value.getCurrency().IsMatrix())
            {
                // check matrix dimensions
                const hwMatrix* matIn = value.getCurrency().Matrix();
                if (matIn->M() == m_numRows)
                {
                    p.setValue(value);
                }
                else
                {
                    Currency colorCur = getPropertyValue("color").getCurrency();
                    hwMatrix* colors = nullptr;
                    if (colorCur.IsMatrix())
                    {
                        colors = colorCur.GetWritableMatrix();
                    }
                    else
                    {
                        colorCur = Currency(new hwMatrix(m_numRows, 3, hwMatrix::REAL));
                        colors = colorCur.GetWritableMatrix();
                    }
                    for (int i = 0; i < m_numRows && i < matIn->M(); ++i)
                    {
                        (*colors)(i, 0) = (*matIn)(i, 0);
                        (*colors)(i, 1) = (*matIn)(i, 1);
                        (*colors)(i, 2) = (*matIn)(i, 2);
                    }
                    p.setValue(colorCur);
                }
                return true;
            }
            std::vector<double> singleColor;
            if (value.isCurrency() && value.getCurrency().IsString())
            {
                Color clr(value.getCurrency().StringVal());
                singleColor = clr.getComponent();
            }
            else if (value.isColor())
            {
                singleColor = value.getColor();
            }

            if (static_cast<int>(singleColor.size()) >= 3)
            {
                vector<double> z = getZdata();
                int dataSize = static_cast<int>(z.size());
                if (_zcolcount == 0)
                    throw OML_Error(OML_ERR_PLOT_UNKNOWN_ERROR);

                int numRows = dataSize / _zcolcount;
                hwMatrix* cmap = new hwMatrix(numRows, 3, hwMatrix::DataType::REAL);
                for (int i = 0; i < numRows; ++i)
                {
                    (*cmap)(i, 0) = singleColor[0];
                    (*cmap)(i, 1) = singleColor[1];
                    (*cmap)(i, 2) = singleColor[2];
                }
                Currency tmp(cmap);
                p.setValue(tmp);
            }
            return true;
        }
        return Drawable::setPropertyValue(name, value);
    }
   
}
