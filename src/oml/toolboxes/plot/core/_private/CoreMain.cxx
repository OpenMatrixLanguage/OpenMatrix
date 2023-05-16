/**
* @file CoreMain.cxx
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

#include "CoreMain.h"
#include <algorithm>
#include "OML_Error.h"
#include <thread>
#include <future>
#include <functional>

namespace omlplot{

    class RepaintTimer {
    public:
        RepaintTimer(const std::chrono::milliseconds& time, const std::function<void()>& func)
            :m_time(time), m_func(func), m_doCall(true)
        {
            m_thread = std::thread([this] {
                std::this_thread::sleep_for(m_time);
                if (m_doCall)
                    m_func();
                });
        }

        void StopTimer()
        {
            m_doCall = false;
        }

        ~RepaintTimer()
        {
            m_doCall = false;
            m_thread.join();
        }

    private:
        std::chrono::milliseconds m_time;
        std::function<void()> m_func;
        std::atomic<bool> m_doCall;
        std::thread m_thread;

    };

    CoreMain *CoreMain::m_instance = nullptr;

    CoreMain::CoreMain()
        :root(new Root), m_timer(nullptr){
    }

    CoreMain::~CoreMain(){
        delete m_timer;
        root.reset(nullptr);
    }

    void CoreMain::datetick(double handle, const std::string& datefmt, int datefmtIdx, const std::string& axis){
        // check that the element is an Axes
        Object* obj = getObject(handle);
        if (!obj)
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);

        if (!obj->isAxes() && obj->getObjectType() != ObjectType::SECONDARYYAXIS)
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);

        if (obj->getObjectType() == ObjectType::SECONDARYYAXIS) {
            static_cast<SecondaryYAxis*>(obj)->setAxisDatetick(datefmt, datefmtIdx);
        }
        else {
            static_cast<Axes*>(obj)->setAxisDatetick(axis, datefmt, datefmtIdx);
        }
    }

    void CoreMain::setUpdateDatetickFlag(double axesHandle, const string& axis)
    {
        m_updateDatetick.push_back(std::make_pair(axesHandle, axis));
    }

    std::vector<std::pair<double, string>> CoreMain::getUpdateDatetickFlag() {
        std::vector<std::pair<double, string>> ret = m_updateDatetick;
        m_updateDatetick.clear();
        return ret;
    }

    void CoreMain::getAxisDatetickOptions(double axesHandle, const string& axis, bool& enabled, std::string& fmt, int& fmtIdx)
    {
        Object* obj = getObject(axesHandle);
        if (!obj)
            return;

        if (obj->getObjectType() == ObjectType::SECONDARYYAXIS) {
            static_cast<SecondaryYAxis*>(obj)->getAxisDatetickOptions(enabled, fmt, fmtIdx);
        }
        else {
            static_cast<Axes*>(obj)->getAxisDatetickOptions(axis, enabled, fmt, fmtIdx);
        }
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
        Object* o = getObject(h);
        return o && o->isFigure();
    }

    bool CoreMain::isAxes(double h){
        Object* o = getObject(h);
        return o && o->isAxes();
    }


    Object *CoreMain::getObject(double h) const {
        return root->getObject(h);
    }

    Object::VALUETYPE CoreMain::getObjectPropertyValue(double h, string name) const {
        Object* o = getObject(h);
        return o->getPropertyValue(name);
    }

    vector<string> CoreMain::getObjectPropertyNames(double h){
        Object *o = getObject(h);
        return o->getPropertyNames();
    }

    bool CoreMain::isPropertySupported(double h, const string& name) {
        Object* o = getObject(h);
        return o->isPropertySupported(name);
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
        repaintLater(h);
        return h;
    }

    // represents a rectange, (x0, y0) is the bottom-left corner, (x1, y1) is the top-right corner
    class SubplotRect
    {
    public:
        SubplotRect(double x0, double y0, double width, double height) 
        :m_x0(x0), m_y0(y0), m_width(width), m_height(height) {

        }
        ~SubplotRect() {}

        double X0() const { return m_x0; }
        double Y0() const { return m_y0; }
        double X1() const { return m_x0 + m_width; }
        double Y1() const { return m_y0 + m_height; }
        double Width() const { return m_width; }
        double Height() const { return m_height; }
        
        //!
        //! Returns true if the rectangle intersects with the other rectangle
        //!
        bool Intersects(const SubplotRect& other) const {
            bool intersects = !(other.X0() > X1() || other.X1() < X0() ||
                other.Y0() > Y1() || other.Y1() < Y0());
            return intersects;
        }
        //!
        //! Returns the rectangle that is produced by the intersection of this 
        //! rectangle and the other rectangle
        //!
        SubplotRect Intersected(const SubplotRect& other) const {
            double x0 = std::max(X0(), other.X0());
            double x1 = std::min(X1(), other.X1());
            double y0 = std::max(Y0(), other.Y0());
            double y1 = std::min(Y1(), other.Y1());
            if (x0 > x1 || y0 > y1)
                return SubplotRect(0, 0, 0, 0);

            return SubplotRect(x0, y0, x1 - x0, y1 - y0);
        }
        //!
        //! Returns true if the rectangle is the same as the other rectangle
        //! 
        bool operator ==(const SubplotRect& other) const {
            bool ret = isSameDouble(m_x0, other.X0()) && isSameDouble(m_y0, other.Y0()) &&
                isSameDouble(m_width, other.Width()) && isSameDouble(m_height, other.Height());
            return ret;
        }
    private:
        double m_x0, m_y0, m_width, m_height;
    };

    double CoreMain::subplot(int row, int col, const vector<int>& active){
        Object *figure = getObject(gcf());

        // calculate position
        std::vector<int> cp = active;
        std::sort(cp.begin(), cp.end());
        int maxIdx = row * col;
        if (cp.front() < 1 || cp.back() > maxIdx)
        {
            throw OML_Error("Cannot create subplot; invalid index");
        }

        std::vector<int> gridLimits = { row, 0, col, 0 };
        std::vector<int>::const_iterator it = active.cbegin();
        for (; it != active.cend(); ++it)
        {
            int idx = *it;
            int r = (idx - 1) / col + 1;
            int c = (idx - 1) % col + 1;
            if (r < gridLimits[0]) gridLimits[0] = r;
            if (r > gridLimits[1]) gridLimits[1] = r;
            if (c < gridLimits[2]) gridLimits[2] = c;
            if (c > gridLimits[3]) gridLimits[3] = c;
        }

        int r = gridLimits[1] - 1;
        int c = gridLimits[2] - 1;

        int W = gridLimits[3] - gridLimits[2] + 1;
        int H = gridLimits[1] - gridLimits[0] + 1;

        double x = double(c) / col;
        double y = double(row - 1 - r) / row;
        double width = 1.0 / col * W;
        double height = 1.0 / row * H;

        SubplotRect newSubplot(x, y, width, height);
        // check if any other plots in that position
        std::vector<double> foundSubplots;
        std::vector<double> ch = figure->getAllChildren();
        std::vector<double>::const_iterator ait = ch.cbegin();
        for (; ait != ch.cend(); ++ait) {
            Object* ax = getObject(*ait);
            if (ax && ax->isAxes()){
                std::vector<double> posVec = ax->getPropertyValue("pos").Vector();
                if ((int)posVec.size() != 4)
                    continue;
                
                SubplotRect tmp(posVec[0], posVec[1], posVec[2] , posVec[3]);
                if (newSubplot == tmp) {
                    // exactly same position - set the other subplot active and return
                    figure->setPropertyValue("currentaxes", *ait);
                    return *ait;
                }
                if (newSubplot.Intersects(tmp)) {
                    SubplotRect intersected = newSubplot.Intersected(tmp);
                    // if both dimensions are > 0.001 long then it is an intersection.
                    if (intersected.Width() > 0.001 && intersected.Height() > 0.001)
                        foundSubplots.push_back(*ait);
                }
                
            }
        }

        // remove any plots in the same position as the new
        std::vector<double>::const_iterator dit = foundSubplots.cbegin();
        for (; dit != foundSubplots.cend(); ++dit) {
            Object* ax = getObject(*dit);
            if (ax && ax->objectCanBeDeleted())
                delete ax;
        }

        // create new 
        Axes *a = allocObject<Axes>();
        vector<double> pos;
        pos.push_back(x); pos.push_back(y);
        pos.push_back(width); pos.push_back(height);
        a->setPropertyValue("pos", pos);

        figure->addChild(a);
        double h = a->getHandle();
        figure->setPropertyValue("currentaxes", h);

        repaintLater(gcf());
        return h;
    }

    vector<double> CoreMain::plot(vector<LineData> &ldVec){
        return _T_2D_PLOT<Line>(ldVec);
    }

    vector<double> CoreMain::bar(vector<LineData> &ldVec){
        std::vector<double> ret = _T_2D_PLOT<HggroupBar>(ldVec);        
        if (!ldVec.empty() && ldVec.front().xCategories.IsCellArray())
            getObject(gca())->setPropertyValue("xcategories", ldVec.front().xCategories);
        return ret;
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
        repaintLater(ah);
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

    vector<double> CoreMain::waterfall(vector<LineData>& ldVec) {
        return _T_3D_PLOT<Waterfall>(ldVec);
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

    vector<Currency> CoreMain::plotyy(vector<LineData>& ldVec)
    {
        vector<Currency> ret;
        std::vector<double> axesVec;
        std::vector<double> lh = _T_2D_PLOT<Line>(ldVec);
        if (lh.size() == 2) {
            Drawable* lineObject = dynamic_cast<Drawable*>(getObject(lh[1]));
            if (lineObject)
                lineObject->plotOnSecondaryYAxis(true);
            Axes* axes = dynamic_cast<Axes*>(getObject(gca()));
            if (axes) {
                axes->setSecondaryYAxisVisible(true);
                axesVec.push_back(axes->getHandle());
                axesVec.push_back(axes->getSecondaryYAxisHandle());
            }
        }
        ret.push_back(axesVec);
        ret.insert(ret.end(), lh.begin(), lh.end());
        repaintLater(gca());
        return ret;
    }

    vector<double> CoreMain::xline(vector<LineData>& ldVec){
        vector<LineData>::const_iterator it = ldVec.cbegin();
        vector<LineData> newLdVec;
        for (; it != ldVec.cend(); ++it) {
            for (size_t i = 0; i < it->y.size(); ++i)
            {
                LineData dt;
                dt.parent = it->parent;
                dt.y.push_back(it->y[i]);
                dt.properties = it->properties;
                dt.values = it->values;
                newLdVec.push_back(dt);
            }
        }
        return _T_2D_PLOT<XLine>(newLdVec);
    }

    vector<double> CoreMain::yline(vector<LineData>& ldVec){
        vector<LineData>::const_iterator it = ldVec.cbegin();
        vector<LineData> newLdVec;
        for (; it != ldVec.cend(); ++it) {
            for (size_t i = 0; i < it->y.size(); ++i)
            {
                LineData dt;
                dt.parent = it->parent;
                dt.y.push_back(it->y[i]);
                dt.properties = it->properties;
                dt.values = it->values;
                newLdVec.push_back(dt);
            }
        }
        return _T_2D_PLOT<YLine>(newLdVec);
    }

    vector<Currency> CoreMain::fanplot(vector<LineData>& ldVec) {

        if (ldVec.size() != 1)
            return vector<Currency>();

        vector<Currency> ret;
        
        LineData data = ldVec[0];
        LineData surfData(data);
        surfData.properties.clear();
        surfData.values.clear();

        // create the surface
        vector<LineData> vld;
        vld.push_back(surfData);
        vector<double> s = surf(vld);
        ret.push_back(s.front());

        // hold on
        hold(gca(), true);
        
        int nLines = (int)data.x.size() / data.xcolcount;
        int nPoints = data.xcolcount;

        // extract curve indices and curve names
        std::vector<double> curveIndices;
        std::vector<string> curveNames;
        if (!data.properties.empty())
        {
            for (int i = 0; i < data.properties.size(); ++i)
            {
                if (data.properties[i] == "curves")
                {
                    if (data.values[i].IsRealVector())
                        curveIndices = data.values[i].Vector();
                }
                else if (data.properties[i] == "text")
                {
                    if (data.values[i].IsCellArray())
                    {
                        HML_CELLARRAY* strCell = data.values[i].CellArray();
                        for (int j = 0; j < strCell->Size(); ++j)
                        {
                            if (!(*strCell)(j).IsString())
                                curveNames.push_back(string(""));
                            else
                                curveNames.push_back((*strCell)(j).StringVal());
                        }
                    }
                }
            }
        }

        bool plotAll = curveIndices.empty();
        std::vector<double> lineHandles, textHandles;
        lineHandles.reserve(nLines);
        for (int i = 0; i < nLines; i++)
        {
            // Create a line, indexing in the 'curves' vector starts from 1!
            std::vector<double>::const_iterator it = std::find(curveIndices.cbegin(), curveIndices.cend(), i + 1);
            if (it != curveIndices.cend() || plotAll)
            {
                int orderIndex = it != curveIndices.cend() ? (int)(it - curveIndices.cbegin()) : -1;
                LineData tmp;
                // set the data
                std::vector<double> tmpX, tmpY, tmpZ;
                for (int j = 0; j < nPoints; ++j)
                {
                    int idx = j * nLines + i;
                    tmpX.push_back(data.x[idx]);
                    tmpY.push_back(data.y[idx]);
                    tmpZ.push_back(data.z[idx]);
                }
                tmp.x = tmpX;
                tmp.y = tmpY;
                tmp.z = tmpZ;

                std::vector<LineData> vl;
                vl.push_back(tmp);
                std::vector<double>  lh = _T_3D_PLOT<Line3>(vl);
                lineHandles.push_back(lh.front());
                
                // increase the line width and set color
                Object* obj = getObject(lineHandles.back());
                obj->setPropertyValue("linewidth", 3);

                std::string cname;
                if (orderIndex < curveNames.size())
                {
                    cname = curveNames[orderIndex];
                    obj->setPropertyValue("displayname", cname);
                }
            }
        }

        Axes* pAxes = dynamic_cast<Axes*>(getObject(gca()));
        pAxes->setLegendVisible(true);
        getObject(pAxes->getLegendHandle())->setPropertyValue("location", "northwest");

        ret.push_back(lineHandles);
        // for compatibility
        ret.push_back(Currency());
        return ret;
    }

    vector<double> CoreMain::ellipse(vector<LineData>& ldVec)
    {
        double parentAxes = -1;
        if (!ldVec.empty())
        {
            double parent = ldVec[0].parent;
            if (isFigure(parent))
            {
                // find the axes
                Object* obj = getObject(parent);
                double h = obj->getPropertyValue("currentaxes").Scalar();
                if (isAxes(h))
                    parentAxes = h;
                else
                    parentAxes = gca();
                ldVec[0].parent = parentAxes;
            }
            else
            {
                parentAxes = gca();
            }
        }
        else {
            parentAxes = gca();
        }
        
        bool h = ishold(parentAxes);
        hold(parentAxes, true);

        vector<double> ret = _T_2D_PLOT<Ellipse>(ldVec);
        hold(parentAxes, h);

        return ret;
    }

    vector<double> CoreMain::rectangle(vector<LineData>& ldVec)
    {
        double parentAxes = -1;
        if (!ldVec.empty())
        {
            double parent = ldVec[0].parent;
            if (isFigure(parent))
            {
                // find the axes
                Object* obj = getObject(parent);
                double h = obj->getPropertyValue("currentaxes").Scalar();
                if (isAxes(h))
                    parentAxes = h;
                else
                    parentAxes = gca();
                ldVec[0].parent = parentAxes;
            }
            else
            {
                parentAxes = gca();
            }
        }
        else
        {
            parentAxes = gca();
        }
        bool h = ishold(parentAxes);
        hold(parentAxes, true);

        vector<double> ret = _T_2D_PLOT<Rectangle>(ldVec);
        hold(parentAxes, h);

        return ret;
    }

    vector<double> CoreMain::pcolor(vector<LineData>& ldVec)
    {
        return _T_3D_PLOT<PColor>(ldVec);
    }

    vector<double> CoreMain::patch(vector<LineData>& ldVec)
    {
        // get current plot
        double ah = 0;
        if (ldVec.size() > 0) {
            ah = ldVec[0].parent;
        }
        if (!isAxes(ah)) {
            ah = gca();
        }

        Axes* axes = dynamic_cast<Axes*>(getObject(ah));
        if (!axes->ishold()) {
            axes->clear();
        }
        vector<double> res;
        Drawable* pLine = nullptr;
        vector<LineData>::iterator it = ldVec.begin();
        for (; it != ldVec.end(); ++it) {
            LineData ld = *it;
            ld.index = (int)axes->getAllChildren().size();
            try {
                pLine = allocObject<Patch>();
                axes->addChild(pLine);
                pLine->init(ld);
                for (int i = 0; i < ld.properties.size(); i++) {
                    if (ld.properties[i] == "vertices" || ld.properties[i] == "faces")
                        continue;
                    pLine->setPropertyValue(ld.properties[i], ld.values[i]);
                }
                res.push_back(pLine->getHandle());
            }
            catch (OML_Error& e) {
                delete pLine;
                throw e;
            }
        }
        repaintLater(ah);
        return res;
    }

    vector<double> CoreMain::stem3(vector<LineData>& ldVec)
    {
        return _T_3D_PLOT<Stem3>(ldVec);
    }

    vector<double> CoreMain::quiver(vector<LineData>& ldVec)
    {
        return _T_3D_PLOT<HggroupVector>(ldVec);
    }

    vector<double> CoreMain::bar3(vector<LineData>& ldVec)
    {
        vector<double> ret;
        // should be only one LineDta in the vector
        if (ldVec.empty())
            return ret;

        LineData ld = ldVec.front();
        if (ld.bar3PerRow)
        {
            vector<LineData> newData;
            // create a new LineData for each z row
            int dataSize = static_cast<int>(ld.z.size());
            if (ld.zcolcount == 0)
                throw OML_Error(OML_ERR_PLOT_UNKNOWN_ERROR);
            int numRows = dataSize / ld.zcolcount;
            for (int i = 0; i < numRows; ++i)
            {
                LineData tmp;
                tmp.style = ld.style;
                int start = i * ld.zcolcount;
                int end = (i + 1) * ld.zcolcount;
                tmp.x = ld.x;
                tmp.y.push_back(i);
                if ((end - 1) < ld.z.size())
                    tmp.z.insert(tmp.z.end(), ld.z.begin() + start, ld.z.begin() + end);
                tmp.zcolcount = ld.zcolcount;
                newData.push_back(tmp);
            }

            ret = _T_3D_PLOT<HggroupBar3D>(newData);
        }
        else
        {
            ret = _T_3D_PLOT<HggroupBar3D>(ldVec);
        }
        return ret;
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
            if (handle > 0)
            repaintLater(handle);
        }
    }

    double CoreMain::gcf(){
        double h = 0;
        h = root->getPropertyValue("currentfigure").Scalar();
        if (h > 0) {
            Object* root = getObject(0.0);
            if (root->getPropertyValue("showhiddenhandles").StringVal() == "on")
                return h;

            Object* fig = getObject(h);
            if (fig && fig->getPropertyValue("handlevisibility").StringVal() == "off")
                h = 0;
        }
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
        if (h > 0) {
            Object* root = getObject(0.0);
            if (root->getPropertyValue("showhiddenhandles").StringVal() == "on")
                return h;

            Object* ax = getObject(h);
            if (ax && ax->getPropertyValue("handlevisibility").StringVal() == "off")
                h = 0;
        }
        
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
        repaintLater(f);
        return (int)fig->getHandle();
    }

    double CoreMain::cla(double a){
        if (! isAxes(a)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *axes = dynamic_cast<Axes *>(getObject(a));
        axes->clear();
        repaintLater(a);
        return axes->getHandle();
    }

    void CoreMain::hold(double axes, bool hold){
        Object* obj = getObject(axes);
        if (!obj)
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
            
        if (obj->getObjectType() == ObjectType::SECONDARYYAXIS) {
            double parenthandle = static_cast<SecondaryYAxis*>(obj)->getParentAxesHandle();
            Object* parent = getObject(parenthandle);
            if (parent->getObjectType() == ObjectType::AXES)
                static_cast<Axes*>(parent)->hold(hold);
            else
                throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        else {

            if (!isAxes(axes)) {
                throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
            }
            Axes* a = dynamic_cast<Axes*>(getObject(axes));
            a->hold(hold);
        }
    }

    bool CoreMain::ishold(double axes){
        Object* obj = getObject(axes);
        if (!obj)
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);

        if (obj->getObjectType() == ObjectType::SECONDARYYAXIS) {
            double parenthandle = static_cast<SecondaryYAxis*>(obj)->getParentAxesHandle();
            Object* parent = getObject(parenthandle);
            if (parent->getObjectType() == ObjectType::AXES)
                return static_cast<Axes*>(parent)->ishold();
            return false;
        }

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
        repaintLater(axes);
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
        repaintLater(axes);
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
        repaintLater(axes);
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
        repaintLater(axes);
    }

    double CoreMain::title(double axes, string str){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        Axes *pAxes = dynamic_cast<Axes *>(getObject(axes));
        VALUETYPE c = pAxes->getPropertyValue("title");
        Object *text = (Object *)c.BoundObject();
        text->setPropertyValue("string", str);
        repaintLater(axes);
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
        repaintLater(axes);
        return text->getHandle();
    }

    double CoreMain::ylabel(double axes, string str){
        if (! isAxes(axes)){
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }

        Text* text = nullptr;
        Object* obj = getObject(axes);
        Axes *pAxes = dynamic_cast<Axes *>(obj);
        if (pAxes) {
        VALUETYPE c = pAxes->getPropertyValue("ylabel");
            text = (Text*)c.BoundObject();
        }
        else {
            SecondaryYAxis* secAxis = dynamic_cast<SecondaryYAxis*>(obj);
            if (!secAxis)
                throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
            Object* pObj = getObject(secAxis->getParentAxesHandle());
            pAxes = dynamic_cast<Axes*>(pObj);
            if (!axes)
                throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);

            VALUETYPE c = secAxis->getPropertyValue("ylabel");
            text = (Text*)c.BoundObject();
        }

        text->setPropertyValue("string", str);
        repaintLater(axes);
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
        repaintLater(axes);
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
        pAxes->setAxisNeedRepaint(0);
        if (data.size() == 4){
            limits[0] = data[2]; limits[1] = data[3];
            pAxes->setPropertyValue("ylim", limits);
            pAxes->setAxisNeedRepaint(1);
        }
        if (data.size() == 6){
            limits[0] = data[4]; limits[1] = data[5];
            pAxes->setPropertyValue("zlim", limits);
            pAxes->setAxisNeedRepaint(3);
        }
        repaintLater(axes);

    }

    void CoreMain::axis(const string& option) {
        Axes* pAxes = dynamic_cast<Axes*>(getObject(gca()));
        if (!pAxes)
            return;

        if (pAxes->is3DPlot()) {
            if (option == "on" || option == "off" || 
                option == "cubical" || option == "tight" ||
                option == "unscaled") {
                pAxes->setAxisOption(option);
                repaintLater(pAxes->getHandle());
                return;
            }
        }
        else {
            if (option == "on" || option == "off" ||
                option == "equal" || option == "normal" ||
                option == "square" || option == "tight") {
                pAxes->setAxisOption(option);
                repaintLater(pAxes->getHandle());
                return;
            }
        }
        AddWarningString("option [" + option + "] is not supported");
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
        pAxes->setAxisNeedRepaint(0);
        repaintLater(axes);
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
        Object* obj = getObject(axes);
        Axes *pAxes = dynamic_cast<Axes *>(obj);
        if (pAxes) {
        pAxes->setPropertyValue("ylim", limits);
            pAxes->setAxisNeedRepaint(1);
        }
        else {
            SecondaryYAxis * secAxis = dynamic_cast<SecondaryYAxis*>(obj);
            if (!secAxis)
                throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
            Object* pObj = getObject(secAxis->getParentAxesHandle());
            pAxes = dynamic_cast<Axes*>(pObj);
            if (!axes)
                throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);

            secAxis->setPropertyValue("ylim", limits);
            pAxes->setAxisNeedRepaint(2);
        }
        repaintLater(axes);
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
        pAxes->setAxisNeedRepaint(3);
        repaintLater(axes);
    }

    double CoreMain::legend(unique_ptr<LegendData> &ld){
        Axes *pAxes = dynamic_cast<Axes *>(getObject(gca()));

        bool locProp = !ld->properties.empty() && ld->properties[0] == "location" && !ld->values.empty();
        // validate legend location
        if (locProp) {
            Currency locCur = ld->values[0];
            if (!(locCur.IsString() || (locCur.IsVector() && locCur.Vector().size() == 2)))
                locProp = false;
        }
        vector<string> ls = ld->legends;
        if (ls.size() == 0){
            if (locProp || !ld->properties.empty()) {
                pAxes->setLegendVisible(true);
            }
            else {
                pAxes->setLegendVisible(!pAxes->getLegendVisible());
            }
        } else if ((ls.size() == 1) && (ls[0] == "on" || ls[0] == "off")){
            pAxes->setLegendVisible(ls[0] == "on");
        } else {
            // hide all
            vector<double> children = pAxes->getAllChildren();
            vector<double>::const_iterator cit = children.cbegin();
            for (; cit != children.cend(); ++cit) {
                Object* obj = getObject(*cit);
                if (obj && dynamic_cast<Drawable*>(obj)) {
                    dynamic_cast<Drawable*>(obj)->setPropertyValue("displayname", "");
                }
            }

            if (ld->handles.empty()) {
                vector<double> ch = pAxes->getPropertyValue("children").Vector();
                int count = min((int)ls.size(), (int)ch.size());
                Drawable* pLine = nullptr;
                for (int i = 0; i < count; i++) {
                    pLine = dynamic_cast<Drawable*>(getObject(ch[i]));
                    if (pLine)
                        pLine->setPropertyValue("displayname", ls[i]);
                }
            }
            else {
                vector<double>::const_iterator hit = ld->handles.cbegin();
                vector<string>::const_iterator lit = ld->legends.cbegin();
                for (; hit != ld->handles.cend() && lit != ld->legends.cend(); ++hit, ++lit) {
                    Object* obj = getObject(*hit);
                    if (obj && dynamic_cast<Drawable*>(obj)) {
                        dynamic_cast<Drawable*>(obj)->setPropertyValue("displayname", *lit);
                    }
                }
            }
            pAxes->setLegendVisible(true);
        }

        Object* lg = getObject(pAxes->getLegendHandle());
        if (lg) {
            for (int i = 0; i < (int)ld->properties.size(); ++i) {
                bool ret = lg->setPropertyValue(ld->properties[i], ld->values[i]);
                if (!ret) {
                    AddWarningString("Property [" + ld->properties[i] + "] is not supported in OpenMatrix");
                }
            }
        }
        
        repaintLater(gca());
        return pAxes->getLegendHandle();
    }

    vector<double> CoreMain::text(unique_ptr<TextData> &td){
        Axes *pAxes = dynamic_cast<Axes *>(getObject(gca()));
        vector<double> res;
        Text *pText = nullptr;

        bool hasZ = td->zpos.size() == td->xpos.size();
        size_t size = td->xpos.size();
        for (int i = 0; i < size; i++){
            try {
                pText = allocObject<Text>();
                pText->setPropertyValue("x", td->xpos[i]);
                pText->setPropertyValue("y", td->ypos[i]);
                if (hasZ)
                    pText->setPropertyValue("z", td->zpos[i]);
                pText->setPropertyValue("string", td->text[i]);
                pAxes->addChild(pText);
                res.push_back(pText->getHandle());
            } catch (OML_Error &e){
                delete pText;
                throw e;
            }
        }
        // set any properties
        vector<double>::const_iterator it = res.cbegin();
        for (; it != res.cend(); ++it) {
            Object* po = getObject(*it);
            for (int i = 0; i < (int)td->properties.size(); ++i) {
                bool ret = po->setPropertyValue(td->properties[i], td->values[i]);
                // show warning only once per property
                if (!ret && it == res.cbegin()) {
                    AddWarningString("Property [" + td->properties[i] + "] is not supported in OpenMatrix");
                }
            }
        }
        repaintLater(gca());
        return res;
    }

    void CoreMain::saveas(double h, string filename, string fmt, int width, int height){
        if (! isFigure(h)){
            throw OML_Error(OML_ERR_PLOT_INVALID_FIGURE_HANDLE);
        }
        if (m_timer)
            drawnow();

        Figure *figure = dynamic_cast<Figure*>(getObject(h));
        figure->saveas(filename, fmt, width, height);
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
        repaintLater(axes);
    }

    void CoreMain::box(double axes, string state)
    {
        if (!isAxes(axes)) {
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }

        Axes* pAxes = dynamic_cast<Axes*>(getObject(axes));
        pAxes->setBorder(state == "on");
        repaintLater(axes);
    }

    double CoreMain::colorbar(double axes, unique_ptr<ColorbarData>& cd)
    {
        if (!isAxes(axes)) {
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }
        bool doRepaint = false;
        Axes* pAxes = static_cast<Axes*>(getObject(axes));
        if (cd->toggleVisibility) {
            pAxes->setColorbarVisible(!pAxes->getColorbarVisible());
            doRepaint = true;
        }
        if (cd->visible) {
            doRepaint = !pAxes->getColorbarVisible();
            pAxes->setColorbarVisible(true);
        }

        Object* cb = getObject(pAxes->getColorbarHandle());
        if (cb && !cd->properties.empty()) {
            for (int i = 0; i < (int)cd->properties.size(); ++i) {
                bool ret = cb->setPropertyValue(cd->properties[i], cd->values[i]);
                if (!ret) {
                    AddWarningString("Property [" + cd->properties[i] + "] is not supported in OpenMatrix");
                }
            }
            doRepaint = true;
        }
        if (doRepaint)
            repaintLater(gca());
        return pAxes->getColorbarHandle();
    }

    std::vector<double> CoreMain::colorbarRange(double axes)
    {
        if (!isAxes(axes)) {
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);
        }

        Axes* pAxes = dynamic_cast<Axes*>(getObject(axes));
        return pAxes->getColorbarRange();
    }
    
    Currency CoreMain::colormap(double h) const
    {
        Object::VALUETYPE cmap = getObjectPropertyValue(h, "colormap");
        return cmap.getCurrency();
    }
    
    void CoreMain::colormap(double h, const Currency& cmap)
    {
        Object* po = getObject(h);
        if (!po)
            throw OML_Error(OML_ERR_PLOT_INVALID_AXES_HANDLE);

        bool ret = po->setPropertyValue("colormap", cmap);
        repaintLater(h);
    }

    void CoreMain::drawnow()
    {
        if (m_timer)
        {
            m_timer->StopTimer();
            delete m_timer;
            m_timer = nullptr;
        }

        repaint();
    }

    void CoreMain::getHandlesForSearch(double h, std::vector<double>& searchObj, int currentDepth, int depth, bool searchAll)
    {
        Object* obj = getObject(h);
        if (!obj)
            return;

        // don't add if already contained
        std::vector<double>::const_iterator it = std::find_if(searchObj.cbegin(),
            searchObj.cend(), [&h](double sh) {return isSameDouble(h, sh); });
        if (it == searchObj.cend())
            searchObj.push_back(h);

        if (currentDepth < depth || depth == -1)
        {
            std::vector<double> ch;
            if (searchAll)
                ch = obj->getAllChildren();
            else
                ch = obj->getPropertyValue("children").Vector();
            std::vector<double>::const_iterator it = ch.cbegin();
            for (; it != ch.cend(); ++it)
                getHandlesForSearch(*it, searchObj, currentDepth + 1, depth, searchAll);
        }
    }

    static bool identicalCurrency(const Currency& c1, const Currency& c2)
    {
        if (c1.IsScalar() && c2.IsScalar())
        {
            return isSameDouble(c1.Scalar(), c2.Scalar());
        }
        else if (c1.IsString() && c2.IsString())
        {
            return c1.StringVal() == c2.StringVal();
        }
        else if (c1.IsMatrix() && c2.IsMatrix())
        {
            const hwMatrix* m1 = c1.Matrix();
            const hwMatrix* m2 = c2.Matrix();
            return (*m1) == (*m2);
        }
        else
        {
            return false;
        }
    }

    std::vector<double> CoreMain::findobj(std::unique_ptr<QueryData>& data, bool searchAll)
    {
        std::vector<double> res;
        std::vector<double> handlesToSearch;
        if (data->handles.empty())
        {
            getHandlesForSearch(0, handlesToSearch, 0, data->m_depth, searchAll);
        }
        else
        {
            std::vector<double>::const_iterator it = data->handles.cbegin();
            for (; it != data->handles.cend(); ++it)
            {
                // throw error if first level item is not a PlotElement
                if (!isHandle(*it))
                    throw OML_Error(OML_ERR_PLOT_INVALID_OBJECT_HANDLE);

                getHandlesForSearch(*it, handlesToSearch, 0, data->m_depth, searchAll);
            }
        }

        std::vector<double>::const_iterator it = handlesToSearch.cbegin();
        for (; it != handlesToSearch.cend(); ++it)
        {
            Object* obj = getObject(*it);
            if (!obj)
                continue;

            std::vector<std::string> elementProps = obj->getPropertyNames();

            int numprops = static_cast<int>(data->properties.size());
            bool hasOperators = !data->ops.empty();
            std::vector<bool> propFound;
            propFound.resize(data->properties.size());
            for (int i = 0; i < numprops; ++i)
            {
                propFound[i] = false;
                // prop exists ?
                bool propExists = std::find(elementProps.cbegin(), elementProps.cend(), data->properties[i]) != elementProps.cend();
                if (propExists)
                {
                    if (data->m_onlyProperty)
                    {
                        propFound[i] = true;
                    }
                    else
                    {
                        CurrencyAndColor p = obj->getPropertyValue(data->properties[i]);
                        Currency pVal;
                        if (p.isColor())
                            pVal = p.getColor();
                        else
                            pVal = p.getCurrency();

                        if (identicalCurrency(pVal, data->values[i]))
                            propFound[i] = true;
                    }
                }

                if (!propFound[i])
                {
                    if (hasOperators)
                        continue;
                    else
                        break;
                }
            }

            if (hasOperators)
            {
                bool P = propFound[0];
                int j = 0;
                for (int i = 1; i < propFound.size() && j < data->ops.size(); ++i, ++j)
                {
                    int OP = data->ops[j];
                    switch (OP)
                    {
                    case QueryData::OP_AND:
                        P = P && propFound[i];
                        break;
                    case QueryData::OP_OR:
                        P = P || propFound[i];
                        break;
                    case QueryData::OP_XOR:
                        P = P != propFound[i];
                        break;
                    case QueryData::OP_NOT:
                        P = P && !propFound[i];
                        break;
                    default:
                        break;
                    }
                }
                if (P)
                    res.push_back(*it);
            }
            else
            {
                auto fit = std::find(propFound.cbegin(), propFound.cend(), false);
                if (fit == propFound.cend())
                    res.push_back(*it);
            }
        }
        return res;
    }

    bool CoreMain::deleteHandle(double h)
    {
        Object* obj = getObject(h);
        if (!obj)
            return false;

        Object* parent = obj->getParentObject();
        if (!parent)
            return false;
        if (obj->objectCanBeDeleted())
            delete obj;
        else
            throw OML_Error(std::string("Element cannot be deleted"));
        repaintLater(parent->getHandle());
        return true;
    }

    void CoreMain::view(double h, const std::vector<double> viewVal)
    {
        if (h != -1 && isAxes(h))
        {
            Axes* a = dynamic_cast<Axes*>(getObject(h));
            a->getParentObject()->setPropertyValue("currentaxes", h);
            root->setPropertyValue("currentfigure", a->getParentObject()->getHandle());
        }

        Axes* ca = dynamic_cast<Axes*>(getObject(gca()));
        ca->setView(viewVal);
        repaintLater(ca->getHandle());
    }

    std::string CoreMain::GetWarningString() {
        std::string ret;
        if (!m_warningStr.empty()) {
            ret = std::string("Warning(s):\n" + m_warningStr);
            m_warningStr = std::string();
        }
        return ret;
    }

    void CoreMain::AddWarningString(const std::string& wrn) {
        if (!m_warningStr.empty())
            m_warningStr += "\n";
        m_warningStr += wrn;
    }

    void CoreMain::repaintLater(double handle)
    {
        if (IsZero(handle))
        {
            // root, add all figures to m_handlesToRepaint
            Object* obj = getObject(handle);
            if (!obj)
                return;
            std::vector<double> fs = obj->getAllChildren();
            std::vector<double>::const_iterator fit = fs.cbegin();
            for (; fit != fs.cend(); ++fit)
            {
                auto it = std::find(m_handlesToRepaint.cbegin(), m_handlesToRepaint.cend(), *fit);
                if (it == m_handlesToRepaint.cend())
                    m_handlesToRepaint.push_back(*fit);
            }
        }
        else
        {
            double tmpHandle = handle;
            while (!isFigure(tmpHandle))
            {
                // find parent figure
                Object* obj = getObject(tmpHandle);
                if (!obj)
                    return;
                Object* parentObj = obj->getParentObject();
                if (!parentObj || parentObj->getObjectType() == ObjectType::ROOT)
                {
                    // no parent figure found
                    tmpHandle = -1;
                    break;
                }
                tmpHandle = parentObj->getHandle();
            }
            if (tmpHandle > 0)
            {
                auto it = std::find(m_handlesToRepaint.cbegin(), m_handlesToRepaint.cend(), tmpHandle);
                if (it == m_handlesToRepaint.cend())
                    m_handlesToRepaint.push_back(tmpHandle);
            }
            else
            {
                return;
            }
        }

        if (m_timer)
        {
            m_timer->StopTimer();
            delete m_timer;
            m_timer = nullptr;
        }
        if (!m_timer)
        {
            std::function<void()> func = std::bind(&CoreMain::repaint, this);
            m_timer = new RepaintTimer(std::chrono::milliseconds(30), func);
        }
    }

    void CoreMain::repaint()
    {
        if (m_handlesToRepaint.empty())
        {
            double h = 0;
            h = root->getPropertyValue("currentfigure").Scalar();
            if (h == 0)
                return;
            Object* obj = getObject(h);
            if (!obj)
                return;
            Figure* fig = dynamic_cast<Figure*>(obj);
            if (!fig)
                return;
            fig->repaint();
        }
        else
        {
            std::vector<double>::const_iterator it = m_handlesToRepaint.cbegin();
            for (; it != m_handlesToRepaint.cend(); ++it)
            {
                Object* obj = getObject(*it);
                if (!obj)
                    continue;
                obj->repaint();
            }
            m_handlesToRepaint.clear();
        }
    }
}
