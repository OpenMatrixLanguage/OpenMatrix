/**
* @file GnuplotOutput.h
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

#ifndef _GNUPLOT_OUTPUT_H_
#define _GNUPLOT_OUTPUT_H_

#include <cstdio>
#include <string>
#include <cstdarg>
#include "OmlPlotExport.h"

namespace omlplot{
    class OMLPLOT_EXPORT GnuplotOutput{
    public:
        GnuplotOutput();
        ~GnuplotOutput();

        static bool isGnuplotReady();
        void setWindowId(double);

#ifdef OS_WIN
        void printf(const std::string fmt, ...){
            va_list args;
            va_start(args, fmt);
            vfprintf(m_handle, fmt.c_str(), args);
#ifdef _DEBUG
            vfprintf(m_debug, fmt.c_str(), args);
#endif
            va_end(args);
        }
#else
        template <typename ... Args>
        void printf(const std::string fmt, Args... args){
            fprintf(m_handle, fmt.c_str(), args...);
#ifdef _DEBUG
            fprintf(m_debug, fmt.c_str(), args...);
#endif
        }
#endif

        void flush();

    private:
        FILE *m_handle;
#ifdef _DEBUG
        FILE *m_debug;
#endif
    };
}

#endif
