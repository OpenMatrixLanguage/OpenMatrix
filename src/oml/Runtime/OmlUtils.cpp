/**
* @file OmlUtils.cpp
* @date October 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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

// Begin defines/includes
#include "OmlUtils.h"

#ifdef OS_WIN
#    include <Windows.h>
#else
#    include <unistd.h>
#    include <glob.h>
#    include <libgen.h>
#    include <dirent.h>
#    include <sys/resource.h>
#    include <sys/times.h>
#    include <time.h>
#    include <dlfcn.h>
#endif

#include "ErrorInfo.h"
#include "OML_Error.h"
// End defines/includes

//------------------------------------------------------------------------------
//! Gets current working directory
//------------------------------------------------------------------------------
std::string OmlUtils::GetCurrentWorkingDir()
{
    std::string workingdir;

#ifdef OS_WIN
    int   expectedSize = GetCurrentDirectory(0, nullptr);
    char* cwd          = new char[expectedSize+1];
    memset(cwd, 0, sizeof(cwd));
    if (!GetCurrentDirectory(expectedSize, cwd))
        throw OML_Error(HW_ERROR_NOTFINDCURWORKDIR);

    workingdir = cwd;
    delete [] cwd;
    cwd = NULL;
#else
    char *cwd = getcwd(nullptr, 0);
    if (!cwd)
        throw OML_Error(HW_ERROR_NOTFINDCURWORKDIR);
    workingdir = cwd;
    free(cwd);
#endif

    return workingdir;
}
//------------------------------------------------------------------------------
//! Gets file names matching pattern in current directory
//! \param[in] pattern Pattern to match
//------------------------------------------------------------------------------
std::vector<std::string> OmlUtils::GetMatchingFiles(const std::string& pattern)
{
    std::vector<std::string> files;
#ifdef OS_WIN
    WIN32_FIND_DATA data;
    HANDLE          handle = FindFirstFile(pattern.c_str(), &data);

    while (handle != INVALID_HANDLE_VALUE)
    {
        files.push_back(std::string(data.cFileName));

        if (!FindNextFile(handle, &data))
            break;
    }
#else
    glob_t data;
    int retcode = glob(pattern.c_str(), GLOB_TILDE, nullptr, &data);

    if (retcode)
    {
        globfree(&data);
        if (retcode == GLOB_NOSPACE)
            throw OML_Error(HW_ERROR_OUTMEM);
        else if (retcode == GLOB_ABORTED)
            throw OML_Error(HW_ERROR_READ);
    }
    else
    {
        for (size_t i = 0; i < data.gl_pathc; ++i)
        {
            files.push_back(basename(data.gl_pathv[i]));
        }

        globfree(&data);
    }
#endif

    return files;
}
