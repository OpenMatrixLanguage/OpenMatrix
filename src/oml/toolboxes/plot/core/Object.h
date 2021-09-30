/**
* @file Object.h
* @date May 2017
* Copyright (C) 2017-2021 Altair Engineering, Inc.  
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

#ifndef _AXES_H_
#define _AXES_H_

#include "OmlPlotExport.h"
#include <string>
#include <vector>
#include <set>
#include <map>
#include <random>
#include <memory>
#include <iostream>
#include "boost/any.hpp"
#include "DataType.h"
#include "Currency.h"
#include "OML_Error.h"
#include "GnuplotOutput.h"

using namespace std;
using namespace boost;

namespace omlplot{

    enum PropertyType{
        NIL = 0x1,
        BOOL = 0x2,
        STRING = 0x4,
        DOUBLE = 0x8,
        MATRIX = 0x10,
        POINTER = 0x20,
        VEC_DOUBLE = 0x40,
        COLOR = 0x80,
        CELL = 0x100,
		UNSUPPORTED = 0x200 //Helper for identifying properties that are not supported in OM yet
    };

    enum class ObjectType {
        NO_TYPE, ROOT, FIGURE, AXES, SECONDARYYAXIS, TEXT, 
        LINE, POLAR_LINE, FILL, AREA, LINE_3D, BAR, HIST, 
        SCATTER, SCATTER_3D, SURFACE, MESH, CONTOUR_3D, 
        CONTOUR, STEM, LOGLOG, SEMILOGX, SEMILOGY
    };

    template <typename T1 = string, typename T2 = Currency,
              typename T3 = PropertyType>
    class PropertyTemplate{
    public:
        typedef T2 ValueType;

        explicit PropertyTemplate(T1 x, T2 y, T3 z)
            :m_name(x), m_value(y), m_type(z){
        }
        T1 getName(){return m_name;}
        T2 getValue() {
            return m_value;
        }
        T3 getType() {return m_type;}
        void setValue(T2 v) {
            m_value = v;
        }

    private:
        T1 m_name;
        T2 m_value;
        T3 m_type;
    };

    class OMLPLOT_EXPORT CurrencyAndColor{
    public:
        enum CCType{
            CURRENCY, COLOR
        };

        CurrencyAndColor(int c);        
        CurrencyAndColor(double c);        
        CurrencyAndColor(const char* c);
        CurrencyAndColor(std::string c);
        CurrencyAndColor(std::vector<double> c);
        CurrencyAndColor(Currency c);
        CurrencyAndColor(Color c);

        double Scalar();
        std::string StringVal();
        std::vector<double> Vector();
        void *BoundObject();
        std::string ColorString();

        Currency getCurrency();
        vector<double> getColor();

        bool isCurrency();
        bool isColor();

    private:
        Currency m_currency;
        Color m_color;
        CCType m_type;
    };

    typedef PropertyTemplate<string, CurrencyAndColor, PropertyType> Property;

    class HandlePool;
    class FigureHandlePool;

    class OMLPLOT_EXPORT Object{
    public:
        typedef Property::ValueType VALUETYPE;

        Object();
        virtual ~Object();

        vector<string> getPropertyNames();
        Property &getProperty(const string&);
        VALUETYPE getPropertyValue(const string&);
        virtual bool setPropertyValue(const string& , VALUETYPE);

        double getHandle();
        string getType();

        bool isHandle(double);
        virtual bool isFigure();
        virtual bool isAxes();
        virtual bool isDrawable();
        virtual bool isText();
        Object *getObject(double);
        void dump(Object *, int);
        Object* getParentObject();

        void setParent(Object *);
        void addChild(Object *);
        void removeChild(Object *);
        virtual void repaint();
        virtual void update(GnuplotOutput *) = 0;

        ObjectType getObjectType() const;

    protected:
        vector<Property> m_ps;         // properties
        static HandlePool *m_handlePool;
        static FigureHandlePool * m_figureHandlePool;
        static map<double, Object *> m_objectMap;
        ObjectType m_type;
    };

    class OMLPLOT_EXPORT Root: public Object{
    public:
        Root();
        void repaint() override;
        void update(GnuplotOutput *out = nullptr) override;
    };

    class Text;

    class OMLPLOT_EXPORT Figure: public Object{
    public:
        Figure();
        Figure(int);
        ~Figure();
        bool isFigure() override;
        void repaint() override;
        void update(GnuplotOutput *out = nullptr) override;
        virtual bool setPropertyValue(const string&, VALUETYPE) override;

        void clear();
        void out(string);
        void saveas(string , string);
        void setGridLayout(int rows, int cols);
        bool isGridLayout();
        void getSubplotPosition(int subplotIdx, double& x, double& y, double& w, double& h);
    private:
        void init(double h);
    private:
        unique_ptr<GnuplotOutput> m_out;
        unique_ptr<Text> m_tLabel, m_rLabel, m_lLabel, m_bLabel;       
        int m_rows, m_cols;
        bool m_isGridLayout;
    };

    class SecondaryYAxis;

    class OMLPLOT_EXPORT Axes: public Object{
    public:
        Axes();
        bool isAxes() override;
        void update(GnuplotOutput *) override;
        void repaint() override;

        void clear();
        void hold(bool);
        bool ishold();

        void setAxisNeedRepaint(int axisID);
        void setBorder(bool state);
        bool getBorder() const;
        void setColorbarVisible(bool state);
        bool getColorbarVisible() const;
        void setColorbarRange(const std::vector<double>& range);
        std::vector<double> getColorbarRange();
        void getBarNumberAndIdx(double barHandle, int& barNumber, int& barIndex);
        void setBarLayout(const string& barLayout);
        void setBarWidth(double barWidth);
        void setSecondaryYAxisVisible(bool visible);
        bool setPropertyValue(const string& name, VALUETYPE value) override;
        double getSecondaryYAxisHandle() const;
        void setGridIndex(int idx);
        void getTipToTailCoordinates(double& x, double& y, double& t, double& r) const;
        void setTipToTailCoordinates(double x, double y, double t, double r);
    private:
        string drawStackedBarPlotLabels(GnuplotOutput* out);
        void setupColorbar(GnuplotOutput* out);

        unique_ptr<Text> m_title;
        unique_ptr<Text> m_xlabel;
        unique_ptr<Text> m_ylabel;
        unique_ptr<Text> m_zlabel;
        unique_ptr<SecondaryYAxis> m_secYAxis;
        bool _borderOn;
        bool _colorbarVisible;
        std::vector<double> _colorbarRange;
        string m_barlayout;
        double m_barWidth;
        bool m_secYAxisVisible;
        bool m_updateXAxisRange, m_updateYAxisRange, m_updateY2AxisRange, m_updateZAxisRange;
        int m_gridIndex;
        double m_tailX, m_tailY, m_tailTheta, m_tailR;
    };

    class OMLPLOT_EXPORT SecondaryYAxis : public Object {
    public:
        SecondaryYAxis(double parentHandle);
        void update(GnuplotOutput*) override;
        bool isAxes() override;
        double getParentAxesHandle() const;

    private:
        unique_ptr<Text> m_ylabel;
        double m_parentAxesHandle;
    };

    class OMLPLOT_EXPORT Drawable : public Object{
    public:
        Drawable();
        virtual void init(const LineData &);
        vector<double> getXdata();
        vector<double> getYdata();
        vector<double> getZdata();
        bool getVisible();

        void update(GnuplotOutput *) override;
        bool isDrawable() override;

        virtual string getUsingClause() = 0;
        virtual string getWithClause(int) = 0;
        virtual string getLegend() = 0;
        virtual string getLineStyle(int) = 0;
        virtual void putData(GnuplotOutput *);
        virtual void cleanup(GnuplotOutput *) = 0;

        void plotOnSecondaryYAxis(bool plot);

    protected:
        int _xcolcount;
        int _ycolcount;
        int _zcolcount;
        std::string m_xaxisRef, m_yaxisRef;
    };

    class OMLPLOT_EXPORT Line: public Drawable{
    public:
        Line();
        void init(const LineData &ld) override;
        string getUsingClause() override;
        string getWithClause(int) override;
        string getLegend() override;
        string getLineStyle(int) override;
        void cleanup(GnuplotOutput *) override;
    };

    class OMLPLOT_EXPORT Polar : public Line{
    public:
        Polar();
        virtual void putData(GnuplotOutput* out) override;
    };

    class OMLPLOT_EXPORT Fill : public Drawable{
    public:
        Fill();
        void init(const LineData &ld) override;
        string getUsingClause() override;
        string getWithClause(int) override;
        string getLegend() override;
        string getLineStyle(int) override;
        void cleanup(GnuplotOutput *) override;
    };

    class OMLPLOT_EXPORT Area : public Drawable{
    public:
        Area();
        void init(const LineData& ld) override;
        string getUsingClause() override;
        string getWithClause(int) override;
        string getLegend() override;
        string getLineStyle(int) override;
        void cleanup(GnuplotOutput *) override;
    };

    class OMLPLOT_EXPORT Line3 : public Line{
    public:
        Line3();
        string getUsingClause() override;
    };

    class OMLPLOT_EXPORT HggroupBar : public Drawable{
    public:
        HggroupBar();

        void init(const LineData &ld) override;
        string getUsingClause() override;
        string getWithClause(int) override;
        string getLegend() override;
        string getLineStyle(int) override;
        void cleanup(GnuplotOutput *) override;
        void putData(GnuplotOutput*) override;
        bool setPropertyValue(const string& name, VALUETYPE value) override;
    };

    class OMLPLOT_EXPORT Hist : public HggroupBar{
    public:
        Hist();

        void init(const LineData &) override;
    };

    class OMLPLOT_EXPORT HggroupScatter : public Drawable{
    public:
        HggroupScatter();

        void init(const LineData& ld) override;
        string getUsingClause() override;
        string getWithClause(int) override;
        string getLegend() override;
        string getLineStyle(int) override;
        void cleanup(GnuplotOutput *) override;
    };

    class OMLPLOT_EXPORT HggroupScatter3 : public HggroupScatter{
    public:
        HggroupScatter3();
        string getUsingClause() override;
    };

    class OMLPLOT_EXPORT Surface : public Drawable{
    public:
        Surface();

        void init(const LineData& ld) override;
        string getUsingClause() override;
        string getWithClause(int) override;
        string getLegend() override;
        string getLineStyle(int) override;
        void putData(GnuplotOutput *) override;
        void cleanup(GnuplotOutput *) override;
        void getMinMaxZ(double& min, double& max);
    private:
        double _minZ, _maxZ;
    };

    class OMLPLOT_EXPORT Mesh : public Surface{
    public:
        Mesh();
        string getUsingClause() override;
        string getWithClause(int) override;
        string getLineStyle(int) override;
    };

    class OMLPLOT_EXPORT Contour3 : public Surface{
    public:
        Contour3();
        string getLineStyle(int) override;
        string getWithClause(int) override;
        string getUsingClause() override;
        void cleanup(GnuplotOutput *) override;
    };

    class OMLPLOT_EXPORT Contour : public Contour3{
    public:
        Contour();
        string getLineStyle(int) override;
        void cleanup(GnuplotOutput *) override;
    };

    class OMLPLOT_EXPORT Stem : public Line{
    public:
        Stem();
        string getWithClause(int lineId) override;        
    };

    class OMLPLOT_EXPORT Loglog : public Line{
    public:
        Loglog();
        string getLineStyle(int) override;
        void cleanup(GnuplotOutput *) override;
    };

    class OMLPLOT_EXPORT Semilogx : public Line{
    public:
        Semilogx();
        string getLineStyle(int) override;
        void cleanup(GnuplotOutput *) override;
    };

    class OMLPLOT_EXPORT Semilogy : public Line{
    public:
        Semilogy();
        string getLineStyle(int) override;
        void cleanup(GnuplotOutput *) override;
    };

    class OMLPLOT_EXPORT Text: public Object{
    public:
        Text();
        void update(GnuplotOutput *) override;
        void update(GnuplotOutput *, const string& type);
        bool isText() override;
        
        bool getVisible();
    private:
        string getModifiedString(const string& text);
    };

    class OMLPLOT_EXPORT HandlePool{
    public:
        HandlePool();

        double allocHandle();
        void releaseHandle(double);
        bool has(double);
    private:
        set<double> m_pool;
        random_device m_rd;
        mt19937 m_gen;
        uniform_real_distribution<> m_dis;
    };

    class OMLPLOT_EXPORT FigureHandlePool{
    public:
        FigureHandlePool();
        int allocHandle();
        bool isFreeHandle(int);
        void setHandleUsed(int);
        void releaseHandle(int);
    private:
        int m_nextHandle;
        set<int> m_pool;
    };

}

#endif
