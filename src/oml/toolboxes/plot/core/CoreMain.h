/**
* @file CoreMain.h
* @date May 2017
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

#ifndef _COREMAIN_H_
#define _COREMAIN_H_


#include "Object.h"
#include "DataType.h"
#include <memory>

using namespace std;

namespace omlplot{

    class OMLPLOT_EXPORT CoreMain{
        typedef Object::VALUETYPE VALUETYPE;
    public:
        static CoreMain *getInstance();
        static void releaseInstance();

        bool isHandle(double);
        bool isFigure(double);
        bool isAxes(double);
        template <typename T>
        bool isObjectPropertyName(string name){
            unique_ptr<T> obj(new T);
            vector<string> names = obj->getPropertyNames();

            transform(name.begin(), name.end(),
                      name.begin(), ::tolower); // ignore the case
            int count = 0;
			vector<string>::iterator it = std::find(names.begin(), names.end(), name);
            return it != names.cend();
        }        

        Object *getObject(double);
        Property getObjectProperty(double, string);
        vector<string> getObjectPropertyNames(double);

        template <typename T>
        vector<double> _T_2D_PLOT(vector<LineData> &ldVec){
            double ah = 0;
            if (ldVec.size() > 0){
                ah = ldVec[0].parent;
            }
            if (! isAxes(ah)){
                ah = gca();
            }
            Axes *axes = dynamic_cast<Axes *>(getObject(ah));
            if (! axes->ishold()){
                axes->clear();
            }
            vector<double> res;
            Drawable *pLine = nullptr;
            vector<LineData>::iterator it = ldVec.begin();
            for (; it != ldVec.end(); ++it){
                LineData ld = *it;
                ld.index = axes->getPropertyValue("children").Vector().size();
                try {
                    pLine = allocObject<T>();
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
        
        template <typename T>
        vector<double> _T_3D_PLOT(vector<LineData> &ldVec){
            double ah = 0;
            if (ldVec.size() > 0){
                ah = ldVec[0].parent;
            }
            if (! isAxes(ah)){
                ah = gca();
            }
            Axes *axes = dynamic_cast<Axes *>(getObject(ah));
            if (! axes->ishold()){
                axes->clear();
            }
            vector<double> res;
            Drawable *pLine = nullptr;
            vector<LineData>::iterator it = ldVec.begin();
            for (; it != ldVec.end(); ++it){
                LineData ld = *it;
                try {
                    pLine = allocObject<T>();
                    pLine->init(ld);
                    for (int i = 0; i < ld.properties.size(); i++){
                        pLine->setPropertyValue(ld.properties[i], ld.values[i]);
                    }
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
        
        vector<double> plot(vector<LineData> &);
        vector<double> bar(vector<LineData> &);
        vector<double> fill(vector<LineData> &);
        vector<double> hist(vector<LineData> &);
        vector<double> area(vector<LineData> &);
        vector<double> polar(vector<LineData> &);
        vector<double> line(vector<LineData> &);
        vector<double> scatter(vector<LineData> &);
        vector<double> plot3(vector<LineData> &);
        vector<double> scatter3(vector<LineData> &);
        vector<double> surf(vector<LineData> &);
        vector<double> mesh(vector<LineData> &);
        vector<double> waterfall(vector<LineData> &);
        vector<double> contour3(vector<LineData> &);
        vector<double> contour(vector<LineData> &);
        vector<double> stem(vector<LineData> &);
        vector<double> loglog(vector<LineData> &);
        vector<double> semilogx(vector<LineData> &);
        vector<double> semilogy(vector<LineData> &);

        double figure(unique_ptr<FigureData> &);
        double axes(unique_ptr<AxesData> &);
        double subplot(int, int, int);
        void set(unique_ptr<SetData> &, vector<string>& );
        double gcf();
        double gca();
        void close(int);
        int clf(int f);
        double cla(double a);
        void hold(double axes, bool hold);
        bool ishold(double axes);
        void grid(double);
        void grid(double, string);
        void minorgrid(double);
        void minorgrid(double, string);
        double title(double, string);
        double xlabel(double, string);
        double ylabel(double, string);
        double zlabel(double, string);
        vector<double> axis(double);
        void axis(double, vector<double>);
        vector<double> xlim(double);
        void xlim(double, vector<double>);
        vector<double> ylim(double);
        void ylim(double, vector<double>);
        vector<double> zlim(double);
        void zlim(double, vector<double>);
        void legend(unique_ptr<LegendData> &);
        vector<double> text(unique_ptr<TextData> &);
        void saveas(double, string, string);
        void dump();
        void out(string);

    private:
        CoreMain();
        ~CoreMain();
        template <typename T>
        T *allocObject(){
            return new T;
        }
        template <typename T, typename U>
        T *allocObject(U id){
            return new T(id);
        }
        static CoreMain *m_instance;

        unique_ptr<Root> root;
    };
    
}

#endif
