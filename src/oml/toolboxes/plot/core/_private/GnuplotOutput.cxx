/**
* @file GnuplotOutput.cxx
* @date May 2020
* Copyright (C) 2020 Altair Engineering, Inc.  
* This file is part of the OpenMatrix Language ("OpenMatrix") software.
* Open Source License Information:
* OpenMatrix is free software. You can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
* OpenMatrix is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for more details.
* You should have received a copy of the GNU Affero General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.
* 
* Commercial License Information: 
* For a copy of the commercial license terms and conditions, contact the Altair Legal Department at Legal@altair.com and in the subject line, use the following wording: Request for Commercial License Terms for OpenMatrix.
* Altair's dual-license business model allows companies, individuals, and organizations to create proprietary derivative works of OpenMatrix and distribute them - whether embedded or bundled with other software - under a commercial license agreement.
* Use of Altair's trademarks and logos is subject to Altair's trademark licensing policies.  To request a copy, email Legal@altair.com and in the subject line, enter: Request copy of trademark and logo usage policy.
*/

#include "GnuplotOutput.h"
#include <iostream>
#ifdef OS_WIN
#include "process.h"
#else
#include "unistd.h"
#include "sys/types.h"
#include "sys/wait.h"
#endif

namespace omlplot{
    GnuplotOutput::GnuplotOutput(){
#ifdef OS_WIN
        m_handle = _popen("gnuplot", "w");
#else
        m_handle = popen("gnuplot", "w");
#endif  // OS_WIN

#ifdef _DEBUG
        m_debug = fopen("plot.log", "w+");
#endif  // _DEBUG

    }

    GnuplotOutput::~GnuplotOutput(){
#ifdef OS_WIN
        _pclose(m_handle);
#else
        pclose(m_handle);
#endif

#ifdef _DEBUG
        fclose(m_debug);
#endif

    }

    void GnuplotOutput::setWindowId(double h){
        int id = int(h);
#ifdef OS_WIN
        fprintf(m_handle, "set terminal wxt title 'figure %d'\n", id);
#else
        fprintf(m_handle, "set terminal x11 title 'figure %d'\n", id);
#endif
    }

    void GnuplotOutput::flush(){
        fflush(m_handle);
        //#ifndef OS_WIN
#ifdef _DEBUG
        fflush(m_debug);
#endif
        //#endif
    }

    bool GnuplotOutput::isGnuplotReady(){
#ifdef OS_WIN
        intptr_t res = _spawnlp(_P_WAIT, "gnuplot", "gnuplot", "-V", NULL);
        if (res == -1){
            return false;
        }
        return true;
#else
        int pid = fork();

        if (pid < 0){
            exit(EXIT_FAILURE);
        } else if (pid == 0){
            freopen("/dev/null", "w", stdout);
            int res = execlp("gnuplot", "gnuplot", "-V", NULL);
            if (res == -1){
                exit(11);
            }
            exit(EXIT_SUCCESS);
        } else if (pid >= 1){
            int wait_status = -1;
            int options = 0;
            waitpid(pid, &wait_status, options);
            if (wait_status == 11){
                return false;
            } else if (wait_status == 0){
                return true;
            } else {
                return false;
            }
        }
        return false;
#endif  // OS_WIN
    }
}
