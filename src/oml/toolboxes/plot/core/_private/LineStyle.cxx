/**
* @file LineStyle.cxx
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

#include "LineStyle.h"
#include "OML_Error.h"

namespace omlplot{

    std::vector<int> LineStyle::c1;
    std::vector<int> LineStyle::c2;
    std::vector<int> LineStyle::c3;
    std::vector<int> LineStyle::c4;
    std::vector<int> LineStyle::c5;
    bool LineStyle::helper = LineStyle::initHelper();

    bool LineStyle::initHelper(){
        c1.push_back(0);   c1.push_back(128); c1.push_back(0);
        c2.push_back(0);   c2.push_back(192); c2.push_back(192);
        c3.push_back(192); c3.push_back(0);   c3.push_back(192);
        c4.push_back(192); c4.push_back(192); c4.push_back(0);
        c5.push_back(64);  c5.push_back(64);  c5.push_back(64);
        return true;
    }

    static Color InitColor[] = {Color(string("blue")),
                                Color(LineStyle::c1),
                                Color(string("red")),
                                Color(LineStyle::c2),
                                Color(LineStyle::c3),
                                Color(LineStyle::c4),
                                Color(LineStyle::c5)
    };
	static const int COLORCOUNT = sizeof(InitColor) / sizeof(InitColor[0]);

	LineStyle::LineStyle(const LineData& ld) {
		try {
			const string& fmt = ld.style;
			string style, legend;
			SplitFormat(fmt, style, legend);
			m_lineStyle = "";
			m_lineColor = InitColor[ld.index % COLORCOUNT];
			m_lineWidth = 1;
			m_markerStyle = "";
			m_markerColor = InitColor[ld.index % COLORCOUNT];
			m_markerSize = 1;
			m_legend = legend;

			if (style == "square") { // full name marker or color
				m_markerStyle = "s";
			}
			else if (style == "diamond") {
				m_markerStyle = "d";
			}
			else if (style == "red" || style == "green" ||
				style == "blue" || style == "cyan" ||
				style == "magenta" || style == "yellow" ||
				style == "white" || style == "black") {
				m_lineColor = Color(style);
			}
			else {                // abbreviated style
				size_t pos = 0;
				size_t size = style.size();
				char cc[2];
				cc[1] = 0;
				char& c = cc[0];
				while (pos < size) {
					c = style[pos];
					if (isAbbrColor(c)) {
						m_lineColor = Color(string(&c));
						pos++;
					}
					else if (isAbbrLine(c)) {
						if ('-' == style[pos]) {
							if ((pos + 1) < size) {
								if ('-' == style[pos + 1]) {
									m_lineStyle = "--";
									pos += 2;
								}
								else if ('.' == style[pos + 1]) {
									m_lineStyle = "-.";
									pos += 2;
								}
								else {
									m_lineStyle = string(&c);
									pos++;
								}
							}
							else {
								m_lineStyle = string(&c);
								pos++;
							}
						}
						else {
							m_lineStyle = string(&c);
							pos++;
						}
					}
					else if (isAbbrMarker(c)) {
						m_markerStyle = string(&c);
						pos++;
					}
					else {
						throw OML_Error(OML_ERR_OPTION);
					}
				}
			}
			if ((m_markerStyle == "") &&
				(m_lineStyle == "")) {
				m_lineStyle = "-";
			}
			if (m_legend == "") {
				char buf[10];
				sprintf(buf, "Curve %d", ld.index + 1);
				m_legend = string(buf);
			}
		}
		catch (...) {
			throw OML_Error(OML_ERR_OPTION);
		}
	}

    void LineStyle::SplitFormat(const string fmt, string &style, string &legend){
        string delim = ";";
        string::size_type pos1 = fmt.find_first_of(delim);
        if (pos1 == string::npos){   // no ";"
            style = fmt;
            legend = "";
        } else {
            string::size_type pos2 = fmt.find_first_of(delim, pos1 + 1);
            if (pos2 == string::npos){
                throw OML_Error(OML_ERR_OPTION);
            }
            style = fmt.substr(0, pos1);
            legend = fmt.substr(pos1 + 1, pos2 - 1 - pos1);
        }
    }

    static char ColorSymbols[] = {'r', 'g', 'b',
                                  'c', 'm', 'y',
                                  'k', 'w'};
    
    bool LineStyle::isAbbrColor(const char &c){
        int count = sizeof(ColorSymbols) / sizeof(ColorSymbols[0]);
        for (int i = 0; i < count; i++){
            if (c == ColorSymbols[i]){
                return true;
            }
        }
        return false;
    }

    static string LineStyles[] = {" ", "-", "--", ":", "-."};
    bool LineStyle::isAbbrLine(const char &l){
        int count = sizeof(LineStyles) / sizeof(LineStyles[0]);
        for (int i = 0; i < count; i++){
            if (l == LineStyles[i][0]){
                return true;
            }
        }
        return false;
    }

    static string MarkerStyles[] = {" ", 
                                    "s", "^", "v",
                                    "x", "o", "d",
                                    "+", "*", "."};

    bool LineStyle::isAbbrMarker(const char &m){
        int count = sizeof(MarkerStyles) / sizeof(MarkerStyles[0]);
        for (int i = 0; i < count; i++){
            if (m == MarkerStyles[i][0]){
                return true;
            }
        }
        return false;
    }
}
