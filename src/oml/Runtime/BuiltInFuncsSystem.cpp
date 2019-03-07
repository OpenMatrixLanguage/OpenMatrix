/**
* @file BuiltInFuncsSystem.cpp
* @date October 2016
* Copyright (C) 2016-2019 Altair Engineering, Inc.  
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

// Begin defines/includes
#include "BuiltInFuncsSystem.h"

#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#ifdef OS_WIN
#    include "windows.h"
#else
#    include <unistd.h>
#    include <dirent.h>
#    include <stdio.h>
#    include <sys/resource.h>
#    include <sys/times.h>
#endif

#include "BuiltInFuncsUtils.h"
#include "CurrencyDisplay.h"
#include "Evaluator.h"
#include "OML_Error.h"
#include "StructData.h"

// End defines/includes

//------------------------------------------------------------------------------
// Returns true after listing contents of a directory [ls command]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Ls(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.empty() ? 0 : inputs.size();

    std::string options;
    std::string filenames;
    std::string out;

    if (nargin > 0)
    {
        Currency cur = inputs[0];
        if (!cur.IsString())
            throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
        options = cur.StringVal();
    }

    if (nargin > 1)
    {
        Currency cur = inputs[1];
        if (!cur.IsString())
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
        filenames = cur.StringVal();
    }

    std::string strcmd;
#ifdef OS_WIN
    strcmd = "dir";
#else
    strcmd =  "ls";
#endif

    if (!options.empty())
        strcmd += " " + options;
    if (!filenames.empty())
        strcmd += " " + filenames;
   
    FILE* cmdoutput = nullptr;

#ifdef OS_WIN
    cmdoutput = _popen(strcmd.c_str(), "r");
#else
    cmdoutput = popen(strcmd.c_str(), "r");
#endif

    while (1) 
    {
        char buf[256];
        memset(buf, 0, sizeof(buf));
        if (fgets(buf, sizeof(buf), cmdoutput) <= 0)
            break;
        out += buf;
    }

#ifdef OS_WIN
    _pclose(cmdoutput);
#else
    pclose(cmdoutput);
#endif

    
    if (eval.GetNargoutValue() <= 0)
    {
         Currency output(out);
         output.DispOutput();
        eval.PrintResult(output);
    }
    else
        outputs.push_back(out);

    return true;
}
//------------------------------------------------------------------------------
// Gets info whether the given path is an existing directory or file
//------------------------------------------------------------------------------
void BuiltInFuncsSystem::IsDirectoryOrFile(const std::string& path,
                                           bool&              isdir,
                                           bool&              isfile) const
{
    if (path.empty())
        return;

    std::string in(path);
    if (BuiltInFuncsUtils::IsRootDir(path))
        in += "/";

    struct stat st;
    int result = stat(in.c_str(), &st);

    if (result != 0)
        return;

    if (st.st_mode & S_IFDIR)
        isdir = true;
#ifdef OS_WIN
    else if (st.st_mode & _S_IFREG)
        isfile = true;
#else
    // Find an alternate option for Linux - these are not defined
    //else if ((st.st_mode & S_ISREG) || (st.st_mode & S_ISFIFO) ||
    //         (st.st_mode & S_ISCHR) || (st.st_mode & S_ISBLK))
    //    isfile = true;
#endif
}
//------------------------------------------------------------------------------
// Returns true after listing directory contents [dir command]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Dir(EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs)
{
    size_t numargs = inputs.empty() ? 0 : inputs.size();
    if (numargs > 1) 
        throw OML_Error(OML_ERR_NUMARGIN);

    std::string dir;
    if (numargs == 0)
        dir = BuiltInFuncsUtils::GetCurrentWorkingDir(); // Use current working dir
    else
    {
        const Currency& cur = inputs[0];
        if (!cur.IsString()) 
            throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
        dir = cur.StringVal();
    }

    bool isdir  = false;
    bool isfile = false;

    BuiltInFuncsUtils::StripTrailingSlash(dir);
    if (BuiltInFuncsUtils::IsRootDir(dir))
        dir += "\\"; 

    BuiltInFuncsSystem funcs;
    funcs.IsDirectoryOrFile(dir, isdir, isfile);

    std::string pattern(dir);
    if (isdir)
        pattern = dir + DIRECTORY_DELIM + std::string("*");

    if (eval.GetNargoutValue() == 0)  // Just print the results
    {
        std::vector<std::string> files (BuiltInFuncsUtils::GetMatchingFiles(pattern));
        if (files.empty())
        {
            std::string warn("Warning: Cannot find directory/files matching '"
                             + pattern + "'");
            BuiltInFuncsUtils::SetWarning(eval, warn);
            return true;
        }
        for (std::vector<std::string>::const_iterator itr = files.begin();
             itr != files.end(); ++itr)
        {
            Currency name  (*itr);
            name.DispOutput();
            eval.PrintResult(name);
        }
        return true;
    }

    // Create an array of structs
    std::vector<std::string> files (BuiltInFuncsUtils::GetMatchingFiles(pattern));

    if (files.empty())
    {
        std::string warn("Warning: Cannot find directory/files matching '"
                            + pattern + "'");
        BuiltInFuncsUtils::SetWarning(eval, warn);
        outputs.push_back(new StructData());
        return true;
    }

    StructData *sd       = new StructData();
    int         numfiles = static_cast<int>(files.size());
    sd->DimensionNew(numfiles, 1);

    // Check if there are wild cards in the given dir
    std::string startupPath;
    if (isdir) 
    {
        startupPath =  dir + DIRECTORY_DELIM;
    }
    else if (isfile || (!isfile && !isdir && dir.find("*") != std::string::npos))
    {
        startupPath = BuiltInFuncsUtils::GetBaseDir(dir) + DIRECTORY_DELIM;
    }

    int i = 0;     
    for (std::vector<std::string>::const_iterator itr = files.begin(); 
         itr != files.end(); ++itr, ++i)
    {
        std::string name (*itr);
        sd->SetValue(i, 0, "name",  name);

        std::string fullname (startupPath + name);
        struct stat st;
        std::string timestr;
        int         numbytes = 0;
        int         isdir    = 0;

        int result = stat(fullname.c_str(), &st);
        if (result == 0)
        {
            timestr  = funcs.GetTimeString(st.st_mtime);
            numbytes = static_cast<int>(st.st_size);
            if (st.st_mode & S_IFDIR)
                isdir = 1;
        }
        sd->SetValue(i, 0, "date",  timestr);
        sd->SetValue(i, 0, "bytes", numbytes);
        sd->SetValue(i, 0, "isdir", isdir);
    }
    outputs.push_back(sd);
    return true;
}
//------------------------------------------------------------------------------
// Gets timestamp as string
//------------------------------------------------------------------------------
std::string BuiltInFuncsSystem::GetTimeString(time_t rawtime)
{
    struct tm * tinfo = localtime (&rawtime);
    if (!tinfo)
        return "";

    char tmp[128];
    sprintf(tmp, "%02d", tinfo->tm_mday);

    std::string out(tmp);

    switch(tinfo->tm_mon)
    {
        case 0:  out += "-Jan-"; break;
        case 1:  out += "-Feb-"; break;
        case 2:  out += "-Mar-"; break;
        case 3:  out += "-Apr-"; break;
        case 4:  out += "-May-"; break;
        case 5:  out += "-Jun-"; break;
        case 6:  out += "-Jul-"; break;
        case 7:  out += "-Aug-"; break;
        case 8 : out += "-Sep-"; break;
        case 9:  out += "-Oct-"; break;
        case 10: out += "-Nov-"; break;
        case 11: out += "-Dec-"; break;
        default: break;
    }

    out += std::to_string(static_cast<long long>(tinfo->tm_year + 1900));
    out += " ";

    sprintf(tmp, "%02d:%02d:%02d", tinfo->tm_hour, tinfo->tm_min, tinfo->tm_sec);
    out += tmp;

    return out;
}
//------------------------------------------------------------------------------
// Returns true after running a system command [system command]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::System(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs, 
                                std::vector<Currency>&       outputs)
{
    size_t nargin = inputs.size();
    if (nargin < 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }
    int  nargout    = eval.GetNargoutValue();
    bool saveoutput = (nargout >= 2);
    bool echo       = (nargout == 0);
    
    FILE* pipe = nullptr;

#ifdef OS_WIN
    std::string command = "\"" + inputs[0].StringVal() + "\" 2>&1";
	pipe = _popen(command.c_str(), "r");
#else
    std::string in(inputs[0].StringVal());

    // Strip leading spaces
    if (!in.empty())
    {
        size_t pos = in.find_first_not_of(" \t");
        if (pos != std::string::npos && pos != 0)
        {
            in = in.substr(pos);
        }
    }
    
    // Extra quotes needed on linux with spaces
    if (!in.empty() && in[0] != '\"' && in.find(" ") != std::string::npos)
    {
        in = "\"" + in + "\"";  
    }
    std::string command = "\"" + in + "\" 2>&1";
    pipe = popen(command.c_str(), "r");
#endif

    if (!pipe)
    {
        throw OML_Error(HW_ERROR_PROBOPENPIPE);
    }

    std::string output;

    // Set the chained display id so that if one of the currencies in this 
    // group is deleted, all of them will be deleted
    time_t ltime;
    time(&ltime);
    long id = static_cast<long>(ltime);

    while (1)
    {
        char buf[256];
        if (fgets(buf, sizeof(buf), pipe) <= 0)
        {
            break;
        }
        output += buf;
        if (echo) // Echo it instead of saving output
        {
            Currency tmp = buf;
            tmp.DispOutput();
            CurrencyDisplay* display = tmp.GetDisplay();
            if (display)
            {
                display->SetChainedDisplay(id);
            }
            eval.PrintResult(tmp);
        }
    }

    int returncode = 0;
#ifdef OS_WIN
    returncode = _pclose(pipe);
#else
    returncode = pclose(pipe);
#endif

    outputs.push_back(returncode);
    if (saveoutput)
    {
        outputs.push_back(output);
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns true after executing the system command on Unix [unix]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Unix(EvaluatorInterface           eval,
                              const std::vector<Currency>& inputs,
                              std::vector<Currency>&       outputs)
{
#ifdef OS_WIN
    outputs.push_back(0);
    int nargout = eval.GetNargoutValue();
    for (int i = 1; i < nargout; ++i)
    {
        outputs.push_back("");
    }
    return true;
#endif

    return BuiltInFuncsSystem::System(eval, inputs, outputs);
}
//------------------------------------------------------------------------------
// Returns true after deleting given file(s) [delete command]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Delete(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
    int nargin = (inputs.empty()) ? 0 : static_cast<int>(inputs.size());
    if (nargin < 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    std::string warn;
    for (int i = 0; i < nargin; ++i)
    {
        Currency cur = inputs[i];
        if (!cur.IsString())
        {
            throw OML_Error(OML_ERR_STRING, i + 1, OML_VAR_TYPE);
        }
        std::string name (cur.StringVal());
        if (name.empty())
        {
            continue;
        }
        name = BuiltInFuncsUtils::Normpath(name);
        int result = true;
        if (name.find("*") == std::string::npos)  // No wildcards
        {
            if (remove(name.c_str()))
            {
                result = false;
            }
        }
        else
        {
            // Use the system command to as the name could contain wildcards
            // System command will not return info about missing files unless
            // the output is sent to a pipe, which is not done here
            std::string cmd;
#ifdef OS_WIN
            cmd = "del \"" + name + "\" /F /Q";
#else
            cmd = "rm -rf " + name + " > /dev/null";
#endif

            if (system(cmd.c_str()) != 0)
            {
                result = false;
            }
        }
        if (!result)
        {
            warn += "File [" + name + "]";
            std::string err = strerror(errno);
            if (!err.empty() && err != "No error")
            {
                warn += ", error [" + err + "]";
            } 
            warn += "\n";
        }
    }

    if (!warn.empty())
    {
        if (warn[warn.size() - 1] == '\n')
        {
            warn.pop_back();
        }
        if (!warn.empty())
        {
            BuiltInFuncsUtils::SetWarning(eval, 
                "Warning: failed to delete the following file(s):\n" + warn);
        }
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns true and gets current environment handle [getcurrentenv]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::GetCurrentEnv(EvaluatorInterface           eval, 
                                       const std::vector<Currency>& inputs, 
                                       std::vector<Currency>&       outputs)
{
	outputs.push_back(eval.GetCurrentEnvHandle());
	return true;
}
//------------------------------------------------------------------------------
// Returns true and gets base environment handle [getbaseenv]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::GetBaseEnv(EvaluatorInterface           eval, 
                                    const std::vector<Currency>& inputs, 
                                    std::vector<Currency>&       outputs)
{
	outputs.push_back(eval.GetBaseEnvHandle());
	return true;
}
//------------------------------------------------------------------------------
// Returns true and gets new environment handle [getnewenv]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::GetNewEnv(EvaluatorInterface           eval, 
                                   const std::vector<Currency>& inputs, 
                                   std::vector<Currency>&       outputs)
{
	outputs.push_back(eval.GetNewEnvHandle());
	return true;
}
//------------------------------------------------------------------------------
// Returns true and gets environment variable value [getenvvalue]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::GetEnvVal(EvaluatorInterface           eval, 
                                   const std::vector<Currency>& inputs, 
                                   std::vector<Currency>&       outputs)
{
    size_t numargs = inputs.empty() ? 0 : inputs.size();
	if (numargs < 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

	if (!inputs[0].IsInteger())
    {
        throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_TYPE);
    }
	int handle = static_cast<int>(inputs[0].Scalar());

	if (!inputs[1].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
    }

	outputs.push_back(eval.GetEnvValue(handle, inputs[1].StringVal()));
	return true;
}
//------------------------------------------------------------------------------
// Returns true and sets environment variable value [setenvvalue]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::SetEnvVal(EvaluatorInterface           eval, 
                                   const std::vector<Currency>& inputs, 
                                   std::vector<Currency>&       outputs)
{
    size_t numargs = inputs.empty() ? 0 : inputs.size();
	if (numargs < 3)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

	if (!inputs[0].IsInteger())
    {
        throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_TYPE);
    }
	int handle = static_cast<int>(inputs[0].Scalar());

	if (!inputs[1].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
    }

	eval.SetEnvValue(handle, inputs[1].StringVal(), inputs[2]);
	return true;
}
//------------------------------------------------------------------------------
// Returns true and clears environment variable value [clearenvvalue]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::ClearEnvVal(EvaluatorInterface           eval, 
                                     const std::vector<Currency>& inputs, 
                                     std::vector<Currency>&       outputs)
{
    size_t numargs = inputs.empty() ? 0 : inputs.size();
	if (numargs < 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

	if (!inputs[0].IsInteger())
    {
        throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_TYPE);
    }
	int handle = static_cast<int>(inputs[0].Scalar());

	if (!inputs[1].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
    }

	eval.RemoveEnvValue(handle, inputs[1].StringVal());
	return true;
}
//------------------------------------------------------------------------------
// Returns true and imports environment handle [importenv]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::ImportEnv(EvaluatorInterface           eval, 
                                   const std::vector<Currency>& inputs, 
                                   std::vector<Currency>&       outputs)
{
	if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

	if (!inputs[0].IsInteger())
    {
        throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_TYPE);
    }
	int env1 = static_cast<int>(inputs[0].Scalar());

	eval.ImportEnv(env1, eval.GetCurrentEnvHandle());
	return true;
}
//------------------------------------------------------------------------------
// Returns true and imports given environment handle [importenvin]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::ImportEnvIn(EvaluatorInterface           eval, 
                                     const std::vector<Currency>& inputs, 
                                     std::vector<Currency>&       outputs)
{
	if (inputs.size() < 2)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

	if (!inputs[0].IsInteger())
    {
        throw OML_Error(OML_ERR_INTEGER, 1, OML_VAR_TYPE);
    }
	int env1 = static_cast<int>(inputs[0].Scalar());

	if (!inputs[1].IsInteger())
    {
        throw OML_Error(OML_ERR_INTEGER, 2, OML_VAR_TYPE);
    }
	int env2 = static_cast<int>(inputs[1].Scalar());

	eval.ImportEnv(env1, env2);
	return true;
}
//------------------------------------------------------------------------------
// Returns true and clones environment [cloneenv]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::CloneEnv(EvaluatorInterface           eval,
                                  const std::vector<Currency>& inputs, 
                                  std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

	// Don't check for input type for clone env
	int env1 = static_cast<int>(inputs[0].Scalar());

    int newenv = eval.GetNewEnvHandle();
    eval.ImportEnv(env1, newenv);
    outputs.push_back(newenv);
	return true;
}
//------------------------------------------------------------------------------
// Returns true after changing current directory [cd]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Cd(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    size_t nargin = (!inputs.empty()) ? inputs.size() : 0;
    if (nargin > 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    
    BuiltInFuncsUtils utils;
    if (nargin == 1)
    {
        if (!inputs[0].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
        }
        std::string val(inputs[0].StringVal());
        std::string dir;

#ifdef OS_WIN
        std::wstring wstr    = utils.StdString2WString(val);
        std::wstring abspath = utils.GetAbsolutePathW(wstr);
        abspath = utils.StripTrailingSlashW(abspath);

        BuiltInFuncsSystem funcs;
        if (!funcs.IsDir(abspath))
        {
            throw OML_Error("Error: invalid directory name in argument 1");
        }

        // Root directories need a trailing slash
        if (abspath.length() == 2 && iswalpha(abspath[0]) && abspath[1] == ':')
        {
            abspath += L"\\";
        }
        if (!SetCurrentDirectoryW(abspath.c_str()))
        {
            throw OML_Error("Error: cannot set current working directory in argument 1");
        }
        std::wstring cwd = utils.GetCurrentWorkingDirW();
        dir = utils.WString2StdString(cwd);
#else
        dir = utils.GetAbsolutePath(val);
        utils.StripTrailingSlash(dir);
        if (!utils.IsDir(dir))
        {
            throw OML_Error("Error: invalid directory name in argument 1");
        }
        if (chdir(dir.c_str()))
        {
            throw OML_Error("Error: cannot set current working directory in argument 1");
        }
#endif
        eval.ResetFuncSearchCache();
        eval.OnChangeDir(dir);       // broadcast change in current dir
    }
        
        
    if (nargin == 0 || eval.GetNargoutValue() > 0)
    {
#ifdef OS_WIN
        std::wstring cwd = utils.GetCurrentWorkingDirW();
        outputs.push_back(utils.WString2StdString(cwd));
#else
        outputs.push_back(utils.GetCurrentWorkingDir());
#endif
    }

    return true;
}
#ifdef OS_WIN
//------------------------------------------------------------------------------
// Returns true if given path is a directory, supports unicode on Windows
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::IsDir(std::wstring& path)
{
    if (path.empty())
    {
        return false;
    }

    // Root directories need a trailing slash
    std::wstring tmppath = path;
    if (tmppath.length() == 2 && iswalpha(tmppath[0]) && tmppath[1] == ':')
    {
        tmppath += L"\\";
    }
    struct _stat64i32 filestat;
    if (_wstat(tmppath.c_str(), &filestat) == -1)
    {
        return false;
    }

    int returncode = (filestat.st_mode & _S_IFDIR);
    if (returncode == 0)
    {
        return false;
    }
    return true;
}
#endif
//------------------------------------------------------------------------------
// Returns true after deleing given directory [rmdir]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Rmdir(EvaluatorInterface           eval,
                               const std::vector<Currency>& inputs,
                               std::vector<Currency>&       outputs)
{
    if (inputs.size() < 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
    }

    bool removeSubdir = false;
    if (inputs.size() > 1)
    {
        if (!inputs[1].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
        }
        std::string val(inputs[1].StringVal());
        if (val != "s" && val != "S")
        {
            throw OML_Error("Error: invalid option; must be 's' in argument 2");
        }
        removeSubdir = true;
    }

    BuiltInFuncsUtils utils;
    std::string dir(inputs[0].StringVal());
    dir = utils.GetAbsolutePath(dir);

    std::string error;
    std::string msgid;
    bool        result = true;

    // Just removes an empty directory
    if (!removeSubdir)
    {
#ifdef OS_WIN
        BOOL returncode = RemoveDirectory((LPCSTR)dir.c_str());
        if (!returncode)
        {
            result = false;
            msgid  = "rmdir";
            DWORD errid = GetLastError();
            if (errid == 0)
            {
                error = "Error removing directory";
            }
            else
            {
                LPSTR buff = nullptr;
                size_t len = FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER |
                    FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                    NULL, errid, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                    (LPSTR)&buff, 0, NULL);
                error = std::string(buff, len);
                utils.StripTrailingNewline(error);
            }
        }
#else
        int returncode = rmdir(dir.c_str());
        if (returncode)
        {
            result = false;
            msgid = "rmdir";
            error = std::string(strerror(errno));
        }
#endif
        outputs.push_back(Currency(result));
        outputs.push_back(error);
        outputs.push_back(msgid);
        return true;
    }

    // With deleting sub folders
    // Use the system command to delete recursively.
    // System command will not return info about missing files unless
    // the output is sent to a pipe, which is not done here
#ifdef OS_WIN
    std::wstring wstr = utils.StdString2WString(dir);
    std::wstring cmd  = L"rmdir \"" + wstr + L"\" /S /Q";
    if (_wsystem(cmd.c_str()) != 0)
    {
        result = false;
        msgid = "rmdir";
        std::string err = strerror(errno);
        error = "Error deleting [" + utils.WString2StdString(wstr) +
            +"]";
        if (!err.empty() && err != "No error")
        {
            error += "; " + err;
        }
    }
#else
    std::string cmd;
    cmd = "rm -rf " + dir + " > /dev/null";
    if (system(cmd.c_str()) != 0)
    {
        result = false;
        msgid = "rmdir";
        std::string err = strerror(errno);
        error = "Error deleting [" + dir + "]";
        if (!err.empty())
        {
            error += "; " + err;
        }
    }
#endif

    outputs.push_back(Currency(result));
    outputs.push_back(error);
    outputs.push_back(msgid);

    return true;
}
