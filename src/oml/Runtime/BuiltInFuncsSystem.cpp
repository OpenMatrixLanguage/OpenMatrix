/**
* @file BuiltInFuncsSystem.cpp
* @date October 2016
* Copyright (C) 2016-2023 Altair Engineering, Inc.  
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
#include <cassert>
#include <limits.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <thread>
#include <time.h>
#include <utility>

#ifdef OS_WIN
#    include "windows.h"
#else
#    include <unistd.h>
#    include <dirent.h>
#    include <pthread.h>
#    include <stdio.h>
#    include <sys/resource.h>
#    include <sys/times.h>
#endif

#include "BuiltInFuncsUtils.h"
#include "CurrencyDisplay.h"
#include "Evaluator.h"
#include "OML_Error.h"
#include "SignalHandlerBase.h"
#include "StructData.h"

// Helper method to run a system command asynchronously
void RunSystemCommand(const std::string& command)
{
	system(command.c_str());
}

// End defines/includes

//------------------------------------------------------------------------------
// Returns true after listing contents of a directory [ls command]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Ls(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>& outputs)
{
    size_t nargin = inputs.empty() ? 0 : inputs.size();

    std::string options;
    std::string filenames;
    std::string out;

    BuiltInFuncsUtils utils;
    if (nargin > 0)
    {
        if (!inputs[0].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
        }
        options = inputs[0].StringVal();

        // Check whether this is an option or a filename
        if (!options.empty() && nargin < 2)
        {
#ifdef OS_WIN
            if (options[0] != '/' && options[0] != '-') // This is a filename
            {
                filenames = utils.Normpath(options);
            }
#else
            if (options[0] != '-')                      // This is a filename
            {
                filenames = utils.Normpath(options);
                options = "";
            }
#endif
        }
    }

    if (nargin > 1)
    {
        if (!inputs[1].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 2, OML_VAR_TYPE);
        }
        filenames = utils.Normpath(inputs[1].StringVal());
    }

    bool echo = (eval.GetNargoutValue() == 0);

#ifdef OS_WIN
    std::string dir;
    if (filenames.empty())
    {
        std::wstring wstr = utils.GetCurrentWorkingDirW();
        dir = utils.WString2StdString(wstr);
    }
    else
    {
        dir = filenames;
    }
    bool isdir  = false;
    bool isfile = false;

    utils.StripTrailingSlash(dir);
    if (utils.IsRootDir(dir))
    {
        dir += "\\";
    }
    BuiltInFuncsSystem funcs;
    funcs.IsDirectoryOrFile(dir, isdir, isfile);

    std::string pattern(dir);
    if (isdir)
    {
        pattern = dir + DIRECTORY_DELIM + std::string("*");
    }
    std::vector<std::string> files(utils.GetMatchingFiles(pattern));
    if (files.empty())
    {
        std::string warn("Warning: Cannot find directory/files matching '"
            + pattern + "'");
        utils.SetWarning(eval, warn);
        return true;
    }
    std::string newline;
    for (std::vector<std::string>::const_iterator itr = files.begin();
        itr != files.end(); ++itr)
    {
        std::string val(*itr);
        isdir = (val == "." || val == "..");
        if (!isdir)
        {
            funcs.IsDirectoryOrFile(val, isdir, isfile);
        }
        if (isdir)
        {
            val = "[" + val + "]";
        }
        if (echo)
        {
            Currency name(val);
            name.DispOutput();
            eval.PrintResult(name);
        }
        else
        {
            out += newline + val;
            newline = "\n";
        }
    }
    if (!echo)
    {
        outputs.push_back(out);
    }
    return true;
#else

    std::string strcmd("ls");

    if (!options.empty())
    {
        strcmd += " " + options;
    }
    if (!filenames.empty())
    {
        strcmd += " " + filenames;
    }

    strcmd += " 2>&1"; // Redirects 2 (standard error) into 1 (standard output)
    std::cout << std::flush;
    FILE* cmdoutput = popen(strcmd.c_str(), "r");
    if (!cmdoutput)
    {
        std::string err(strerror(errno));
        std::string msg = (!err.empty() && err != "No error") ? err : 
            OML_Error(HW_ERROR_PROBOPENPIPE).GetErrorMessage();

        throw OML_Error(msg);
    }

    while (1)
    {
        char buf[256];
        memset(buf, 0, sizeof(buf));
        if (fgets(buf, sizeof(buf), cmdoutput) <= 0)
        {
            std::cout << std::flush;
            break;
        }
        std::cout << std::flush;
        out += buf;

        if (echo)
        {
            std::string line(buf);
            utils.StripTrailingNewline(line);
            if (line.empty())
            {
                continue;
            }
            Currency tmp = line;
            tmp.DispOutput();
            eval.PrintResult(tmp);
        }
    }

    pclose(cmdoutput);
#endif

    if (!echo)
    {
        outputs.push_back(out);
    }
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
    in = BuiltInFuncsUtils::Normpath(in);
    BuiltInFuncsUtils utils;
    if (utils.IsRootDir(path))
        in += "/";

    int result = 0;
#ifdef OS_WIN
    std::wstring wstr(utils.StdString2WString(in));
    struct _stat64i32 st;
    result = _wstat(wstr.c_str(), &st);
#else
    struct stat st;
    result = stat(in.c_str(), &st);
#endif

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

    BuiltInFuncsUtils utils;
    std::string dir;
    if (numargs == 0)
    {
#ifdef OS_WIN
        std::wstring wstr = utils.GetCurrentWorkingDirW();
        dir = utils.WString2StdString(wstr);
#else
        dir = utils.GetCurrentWorkingDir(); // Use current working dir
#endif
    }
    else
    {
        const Currency& cur = inputs[0];
        if (!cur.IsString()) 
            throw OML_Error(OML_ERR_STRING, 1, OML_VAR_TYPE);
        dir = cur.StringVal();
    }

    dir = utils.Normpath(dir);

    bool isdir  = false;
    bool isfile = false;
    utils.StripTrailingSlash(dir);
    if (utils.IsRootDir(dir))
    {
        dir += "\\";
        dir = utils.Normpath(dir);
    }

    BuiltInFuncsSystem funcs;
    funcs.IsDirectoryOrFile(dir, isdir, isfile);

    std::string pattern(dir);
    if (isdir)
        pattern = dir + DIRECTORY_DELIM + std::string("*");

    bool echo = (eval.GetNargoutValue() == 0);  // Just print the results
    std::vector<std::string> files(utils.GetMatchingFiles(pattern));

    std::unique_ptr<StructData> sd(EvaluatorInterface::allocateStruct());
    if (files.empty())
    {
        std::string warn("Warning: Cannot find directory/files matching '"
            + pattern + "'");
        utils.SetWarning(eval, warn);
        if (!echo)
        {
            outputs.push_back(sd.release());
        }
        return true;
    }

    if (echo)  // Just print the results
    {
        for (std::vector<std::string>::const_iterator itr = files.begin();
             itr != files.end(); ++itr)
        {
            std::string val(*itr);
            isdir = (val == "." || val == "..");
            if (!isdir)
            {
                funcs.IsDirectoryOrFile(val, isdir, isfile);
            }
            if (isdir)
            {
                val = "[" + val + "]";
            }
            Currency name(val);
            name.DispOutput();
            eval.PrintResult(name);
        }
        return true;
    }

    // Create an array of structs
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
        long long      numbytes = 0;
        int         isdir    = 0;

        int result = stat(fullname.c_str(), &st);
        if (result == 0)
        {
            timestr  = funcs.GetTimeString(st.st_mtime);
            numbytes = static_cast<long long>(st.st_size);
            isdir = (st.st_mode & S_IFDIR) ? 1 : 0;
        }
        else // (errno == 132) // File maybe too big
        {
            // Try and open the file
            std::ifstream ifs(fullname, std::ios::binary);
            if (ifs.is_open())
            {
                ifs.seekg(0, std::ios::end);
                numbytes = static_cast<long long>(ifs.tellg());
                ifs.close();
            }
        }
        sd->SetValue(i, 0, "date",  timestr);
        sd->SetValue(i, 0, "bytes", static_cast<double>(numbytes));
        sd->SetValue(i, 0, "isdir", isdir);
    }
    outputs.push_back(sd.release());
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
    std::vector<Currency>& outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    else if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }
    std::string in(inputs[0].StringVal());
    in = BuiltInFuncsUtils::LTrim(in);
    in = BuiltInFuncsUtils::RTrim(in);
    if (in.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }

    BuiltInFuncsUtils utils;

    bool returnoutput      = false; // True if output from system command is returned
    bool async             = false; // True if system command runs asynchronously
    bool maxThreadAffinity = false; // True if async command needs to run on all cores on Linux
    int  nargin       = static_cast<int>(inputs.size());
    for (int i = 1; i < nargin; ++i)
    {
        const Currency& cur = inputs[i];
        if (cur.IsLogical())
        {
            int val = static_cast<int>(cur.Scalar());
            returnoutput = (val == 1) ? true : false;
        }
        else if (cur.IsInteger())
        {
            int val = static_cast<int>(cur.Scalar());
            if (val != 0 && val != 1)
            {
                throw OML_Error(OML_ERR_LOGICAL, i + 1);
            }
            returnoutput = (val == 1) ? true : false;
        }
        else if (cur.IsString())
        {
            std::string opt(cur.StringVal());
            if (opt.empty())
            {
                throw OML_Error(OML_ERR_NONEMPTY_STR, i + 1);
            }
            std::transform(opt.begin(), opt.end(), opt.begin(), ::tolower);
            if (opt == "async")
            {
                async = true;
            }
            else if (opt == "sync")
            {
                async = false;
                returnoutput = true;
            }
            else if (opt == "threadaffinity")
            {
                // Get the value
                ++i;
                if (i >= nargin)
                {
                    throw OML_Error(OML_ERR_NUMARGIN);
                }
                const Currency& val = inputs[i];
                if (!val.IsString())
                {
                    throw OML_Error(OML_ERR_STRING, i + 1);
                }
                std::string threadopt(val.StringVal());
                if (threadopt.empty())
                {
                    throw OML_Error(OML_ERR_NONEMPTY_STR, i + 1);
                }
                std::transform(threadopt.begin(), threadopt.end(), threadopt.begin(), ::tolower);
                if (threadopt == "max")
                {
                    maxThreadAffinity = true;
                }
                else if (threadopt == "default")
                {
                    maxThreadAffinity = false;
                }
                else
                {
                    throw OML_Error(OML_ERR_BAD_STRING, i + 1);
                }
            }
            else
            {
                throw OML_Error(OML_ERR_BAD_STRING, i + 1);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_BAD_STRING, i + 1);
        }
    }


    int  nargout = eval.GetNargoutValue();
    bool saveoutput = false;

    if (nargout >= 2 || returnoutput)
    {
        if (async)
        {
            throw OML_Error("Error: invalid option 'async'; must be 'sync' to return output");
        }
        saveoutput = true;
    }

    if (maxThreadAffinity && !async)
    {
        throw OML_Error("Error: invalid option 'threadaffinity'; can be used only in 'async' mode");
    }

    // Construct the command
    std::string command(GetInputForSystemCommand(in));

    if (async)
    {

#ifndef OS_WIN
        cpu_set_t cpuset;
        if (maxThreadAffinity)
        {
            const auto ncores = std::thread::hardware_concurrency();
            // Set affinity mask to include all cores
            CPU_ZERO(&cpuset);
            for (int j = 0; j < ncores; ++j)
            {
                CPU_SET(j, &cpuset);
            }
        }
        std::thread t(RunSystemCommand, command);

        if (maxThreadAffinity)
        {
            // Set thread affinity to use all cores
            int s = pthread_setaffinity_np(t.native_handle(), sizeof(cpuset), &cpuset);
        }

		std::stringstream ss;
		ss << t.get_id();
		outputs.push_back(static_cast<int>(std::stoull(ss.str())));
        t.detach();    
#else
        std::thread t(RunSystemCommand, command);

        std::stringstream ss;
        ss << t.get_id();
        outputs.push_back(static_cast<int>(std::stoull(ss.str())));
        t.detach();
#endif
		return true;
	}

    bool echo = (nargout < 2);
    FILE* pipe = nullptr;

    int returncode = 0;
    std::string       output;
    // Set the chained display id so that if one of the currencies in this 
    // group is deleted, all of them will be deleted
    long id = 0;
    if (saveoutput || echo)
    {
        id = utils.CreateChainedDisplayId();
    }

#ifdef OS_WIN
    char* oldlocale = setlocale(LC_ALL, nullptr);
    std::string cachelocale("C");
    if (!oldlocale)
    {
        utils.SetWarning(eval, "Warning: Unable to determine locale");
    }
    else
    {
        cachelocale = oldlocale;
    }
    UINT oldoutputcp = GetConsoleOutputCP();

	pipe = _wpopen(utils.StdString2WString(command).c_str(), L"r");
#else
    pipe = popen(command.c_str(), "r");
#endif

    if (!pipe)
    {
        throw OML_Error(HW_ERROR_PROBOPENPIPE);
    }

#ifdef OS_WIN
    BuiltInFuncsSystem funcs;
    std::wstring line;
    if (saveoutput || echo)
    {
        if (oldlocale)
        {
            funcs.SetToUtf8Locale(); // Set to utf locale before getting char
            fflush(stdout);
        }

        while (!feof(pipe))
        {
            wint_t c = fgetwc(pipe);
            if (c == L'\r' || c == L'\n' || c == WEOF)
            {
                if (c == WEOF && line.empty())
                    continue;
                std::string buf(utils.WString2StdString(line));
                output += buf + "\n";
                line = L""; // Reset line
                if (echo)
                {
                    Currency tmp(buf);
                    tmp.DispOutput();
                    CurrencyDisplay* display = tmp.GetDisplay();
                    if (display)
                    {
                        display->SetChainedDisplay(id);
                    }
                    if (oldlocale)               //  Reset to original locale
                    {
                        setlocale(LC_ALL, cachelocale.c_str());
                        SetConsoleOutputCP(oldoutputcp);
                        fflush(stdout);
                    }
                    eval.PrintResult(tmp);
                    if (oldlocale)
                    {
                        funcs.SetToUtf8Locale(); // Set to utf locale before getting char
                        fflush(stdout);
                    }
                }
            }
            else if (c != WEOF)
            {
                line += c;
            }
        }
        if (oldlocale)               //  Reset to original locale
        {
            setlocale(LC_ALL, cachelocale.c_str());
            SetConsoleOutputCP(oldoutputcp);
            fflush(stdout);
        }

    }
    returncode = _pclose(pipe);
#else
	if (saveoutput || echo)
	{
		while (1)
		{
			char buf[256];
			if (fgets(buf, sizeof(buf), pipe) <= 0)
			{
				break;
			}
			if (echo) // Echo it instead of saving output
			{
                std::string val(buf);
                utils.StripTrailingNewline(val);
				Currency tmp = val;
				tmp.DispOutput();
				CurrencyDisplay* display = tmp.GetDisplay();
				if (display)
				{
					display->SetChainedDisplay(id);
				}
				eval.PrintResult(tmp);
			}
            else
            {
                output += buf;
            }
		}
	}
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
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
#ifdef OS_WIN
    outputs.push_back(1);
    int nargout = eval.GetNargoutValue();
    for (int i = 1; i < nargout; ++i)
    {
        outputs.push_back("");
    }
    return true;
#else
    return BuiltInFuncsSystem::System(eval, inputs, outputs);
#endif
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

    BuiltInFuncsUtils utils;
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

#ifdef OS_WIN
        std::wstring wname(utils.StdString2WString(name));
        if (wname.find(L"*") == std::wstring::npos)  // No wildcards
        {
            if (_wremove(wname.c_str()))
            {
                result = false;
            }
        }
        else
        {
            // Use the system command to as the name could contain wildcards
            // System command will not return info about missing files unless
            // the output is sent to a pipe, which is not done here
            std::wstring cmd = L"del \"" + wname + L"\" /F /Q";
            if (_wsystem(cmd.c_str()) != 0)
            {
                result = false;
            }
        }
        if (!result)
        {
            warn += "File [" + utils.WString2StdString(wname) + "]";
            std::string err = strerror(errno);
            if (!err.empty() && err != "No error")
            {
                warn += ", error [" + err + "]";
            }
            warn += "\n";
        }

#else
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
            std::string cmd = "rm -rf " + name + " > /dev/null";
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

#endif
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
// Returns true and gets new environment handle [getnewenv]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::DeleteEnv(EvaluatorInterface eval,
	const std::vector<Currency>& inputs,
	std::vector<Currency>& outputs)
{
	if (inputs.empty())
	{
		throw OML_Error(OML_ERR_NUMARGIN);
	}

	int env1 = static_cast<int>(inputs[0].Scalar());
	eval.DeleteEnv(env1);
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
        val = utils.Normpath(val);
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
        if (dir.length() > 1)
        {
            utils.StripTrailingSlash(dir);
        }
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
    dir = utils.Normpath(dir);

#ifndef OS_WIN
    if (utils.GetCurrentWorkingDir() == dir)
    {
        // Undefined behavior with rmdir when deleting current working dir
        std::string newcwd(utils.GetBaseDir(dir));
        std::vector<Currency> inputs2;
        inputs2.push_back(newcwd);
        std::vector<Currency> outputs2;
        BuiltInFuncsSystem::Cd(eval, inputs2, outputs2);
    }
#endif

    std::string error;
    std::string msgid;
    bool        result = true;

    // Just removes an empty directory
    if (!removeSubdir)
    {
#ifdef OS_WIN
        std::wstring wdir(utils.StdString2WString(inputs[0].StringVal()));
        wdir = utils.GetAbsolutePathW(wdir);

        BOOL returncode = RemoveDirectoryW((LPCWSTR)wdir.c_str());
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
    std::wstring wdir(utils.StdString2WString(inputs[0].StringVal()));
    wdir = utils.GetAbsolutePathW(wdir);
    std::wstring cmd  = L"rmdir \"" + wdir + L"\" /S /Q";
    if (_wsystem(cmd.c_str()) != 0)
    {
        result = false;
        msgid = "rmdir";
        std::string err = strerror(errno);
        error = "Error deleting [" + utils.WString2StdString(wdir) +
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
//------------------------------------------------------------------------------
// Gives the current time in seconds since Jan 1, 1970 [time]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Time(EvaluatorInterface           eval,
                              const std::vector<Currency>& inputs,
                              std::vector<Currency>&       outputs)
{
    if (!inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    time_t seconds = time(nullptr);
    outputs.push_back(static_cast<double>(seconds));
    return true;
}
//------------------------------------------------------------------------------
// Gives the current time as a string [ctime]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::CTime(EvaluatorInterface           eval,
                               const std::vector<Currency>& inputs,
                               std::vector<Currency>&       outputs)
{
    if (inputs.size() != 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    if (!inputs[0].IsScalar())
    {
        throw OML_Error(OML_ERR_NONNEGATIVE_SCALAR, 1);
    }

    double val = inputs[0].Scalar();
    if (val < 0)
    {
        throw OML_Error(OML_ERR_NONNEGATIVE_SCALAR, 1);
    }
    else if (IsInf_T(val) || IsNaN_T(val) || IsNegInf_T(val))
    {
        throw OML_Error(OML_ERR_INVALIDTIMEVAL, 1);
    }
    time_t rawtime = static_cast<time_t>(val);

    struct tm* mytime = localtime(&rawtime);
    if (!mytime)
    {
        throw OML_Error(OML_ERR_INVALIDTIMEVAL, 1);
    }
    char buff[256];
    if (strftime(buff, sizeof(buff), "%a %b %d %H:%M:%S %Y", mytime))
    {
        outputs.push_back(std::string(buff));
    }
    else
    {
        throw OML_Error(OML_ERR_INVALIDTIMEVAL, 1);
    }
    return true;
}
//------------------------------------------------------------------------------
// Gets a path constructed from the given dir/sub-directories [genpath]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Genpath(EvaluatorInterface           eval,
                                 const std::vector<Currency>& inputs,
                                 std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    BuiltInFuncsUtils utils;
    std::string path (inputs[0].StringVal());
    path = utils.Normpath(path);
    if (!utils.DoesPathExist(path))
    {
        outputs.push_back("");
        return true;
    }

    int nargin = static_cast<int>(inputs.size());

#ifdef OS_WIN
    // Directories that need to be excluded
    std::vector< std::pair<std::wstring, bool> > exclude;
    exclude.reserve(nargin);
    for (int i = 1; i < nargin; ++i)
    {
        if (!inputs[i].IsString())
        {
            throw OML_Error(OML_ERR_STRING, i + 1);
        }

        std::string tmp(inputs[i].StringVal());
        if (tmp.empty())
        {
            continue;
        }
        std::wstring wtmp = utils.StdString2WString(tmp);
        std::transform(wtmp.begin(), wtmp.end(), wtmp.begin(), ::tolower);

        bool isabsolute = utils.IsAbsolutePath(tmp);
        if (isabsolute)
        {
            wtmp = utils.GetNormpathW(wtmp);
        }
        exclude.push_back(std::pair<std::wstring, bool>(wtmp, isabsolute));
    }

    std::wstring wpath = utils.StdString2WString(path);
    wpath = utils.GetAbsolutePathW(wpath);
    
    std::wstring paths = wpath;  // Combined paths

    BuiltInFuncsSystem funcs;
    funcs.ListDirW(wpath, exclude, paths);
    outputs.push_back(utils.WString2StdString(paths));

#else
    std::vector< std::pair<std::string, bool> > exclude;
    exclude.reserve(nargin);
    for (int i = 1; i < nargin; ++i)
    {
        if (!inputs[i].IsString())
        {
            throw OML_Error(OML_ERR_STRING, i + 1);
        }

        std::string tmp(inputs[i].StringVal());
        if (tmp.empty())
        {
            continue;
        }
        std::transform(tmp.begin(), tmp.end(), tmp.begin(), ::tolower);

        bool isabsolute = utils.IsAbsolutePath(tmp);
        if (isabsolute)
        {
            tmp = utils.Normpath(tmp);
        }
        exclude.push_back(std::pair<std::string, bool>(tmp, isabsolute));
}

    path = utils.GetAbsolutePath(path);

    std::string paths (path);
    BuiltInFuncsSystem funcs;
    funcs.ListDir(path, exclude, paths);
    outputs.push_back(paths);

#endif
    return true;
}
#ifdef OS_WIN
//------------------------------------------------------------------------------
// Lists directories, recursively; helper for genpath. Supports unicode
//------------------------------------------------------------------------------
void BuiltInFuncsSystem::ListDirW(const std::wstring&   parent,
    const std::vector< std::pair<std::wstring, bool> >& exclude,
    std::wstring&                                       paths)
{
    BuiltInFuncsUtils utils;
    WIN32_FIND_DATAW  find_data;

    std::wstring base = parent + L"\\*";

    HANDLE handle = FindFirstFileW((LPCWSTR)base.c_str(), &find_data);

    while (handle != INVALID_HANDLE_VALUE)
    {
        std::wstring name    = find_data.cFileName;
        bool         addpath = true;
        if (name == L"." || name == L"..")
        {
            addpath = false; // Ignore this
        }
        else if (find_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
        {
            std::wstring subdir = parent;
            utils.AddTrailingSlashW(subdir);
            subdir += name;

            std::wstring lowername   = name;
            std::wstring lowersubdir = subdir;
            std::transform(lowersubdir.begin(), lowersubdir.end(), 
                           lowersubdir.begin(), ::tolower);
            for (std::vector< std::pair<std::wstring, bool> > ::const_iterator itr = exclude.begin();
                 itr != exclude.end(); ++itr)
            {
                std::pair<std::wstring, bool> val = *itr;

                std::wstring lower = lowername;
                std::wstring ignorepath = (*itr).first;

                if ((*itr).second)  // If absolute use the subdirectory
                {
                    lower = lowersubdir;
                }
                if (lower == ignorepath)
                {
                    addpath = false;
                    break;
                }
            } 
            if (addpath)
            {
                std::wstring quotes;
                if (subdir.find(L" ") != std::wstring::npos)
                {
                    quotes = L"\"";
                }

                paths += L";" + quotes + subdir + quotes;
                ListDirW(subdir, exclude, paths);
            }
        }
        if (!FindNextFileW(handle, &find_data))
        {
            break;
        }
    }    
    FindClose(handle);
}
#else
//------------------------------------------------------------------------------
// Lists directories on Linux, recursively; helper for genpath.
//------------------------------------------------------------------------------
void BuiltInFuncsSystem::ListDir(const std::string& parent,
    const std::vector< std::pair<std::string, bool> >& exclude,
    std::string& paths)
{

    BuiltInFuncsUtils utils;

    std::string base(parent);
    utils.AddTrailingSlash(base);

    DIR* dir = opendir(base.c_str());
    struct dirent* epdf;
    if (dir)
    {
        while (1)
        {
            struct dirent* contents = readdir(dir);
            if (!contents)
            {
                break;
            }
            std::string name(contents->d_name);
            bool        addpath = true;
            if (name == "." || name == "..")
            {
                addpath = false; // Ignore this
            }
            else if (contents->d_type == DT_DIR)
            {
                std::string subdir = base + name;

                std::string lowername = name;
                std::string lowersubdir = subdir;
                for (std::vector< std::pair<std::string, bool> > ::const_iterator itr = exclude.begin();
                    itr != exclude.end(); ++itr)
                {
                    std::pair<std::string, bool> val = *itr;

                    std::string lower = lowername;
                    std::string ignorepath = (*itr).first;

                    if ((*itr).second)  // If absolute use the subdirectory
                    {
                        lower = lowersubdir;
                    }
                    if (lower == ignorepath)
                    {
                        addpath = false;
                        break;
                    }
                }
                if (addpath)
                {
                    if (!paths.empty())
                    {
                        paths += ":";
                    }
                    paths += subdir;
                    ListDir(subdir, exclude, paths);
                }
            }
        }
    }
    closedir(dir);
}
#endif
//------------------------------------------------------------------------------
// Gets process id [getpid]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::GetPid(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs)
{
#ifdef OS_WIN
    DWORD pid = GetCurrentProcessId();
    outputs.push_back(static_cast<int>(pid));
#else
    pid_t pid = getpid();
    outputs.push_back(static_cast<int>(pid));
#endif

    return true;
}
//------------------------------------------------------------------------------
// Creates an absolute path name [make_absolute_filename]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::MakeAbsFilename(EvaluatorInterface           eval,
    const std::vector<Currency>& inputs,
    std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    std::string in(inputs[0].StringVal());
    if (in.empty())
    {
        outputs.push_back(std::string());
        return true;
    }

    BuiltInFuncsUtils utils;
    std::string path;

#ifdef OS_WIN
    std::wstring wstr = utils.StdString2WString(in);
    std::wstring abspath = utils.GetAbsolutePathW(wstr); // Does not need path to exist

    path = utils.WString2StdString(abspath);
    path = utils.Normpath(path);
#else
    path = utils.GetAbsolutePath(in); // Needs path to exist
    if (path.empty())
    {
        in = utils.Normpath(in);
        path = utils.GetCurrentWorkingDir();
        if (in.find("./") == 0)
        {
            path += "/" + in.substr(2);
        }
        else if (in[0] == '/')
        {
            path += in;
        }       
        else 
        {
            path += "/" + in;
        }
    }
#endif

    outputs.push_back(path);
    return true;

}
#ifdef OS_WIN
//------------------------------------------------------------------------------
// Helper method to set locale and output mode to support Unicode
//------------------------------------------------------------------------------
void BuiltInFuncsSystem::SetToUtf8Locale()
{
    setlocale(LC_ALL, "en_US.UTF-8");
    SetConsoleOutputCP(65001);
}
#endif
//------------------------------------------------------------------------------
// Helper method to correctly quote command that is sent to system
//------------------------------------------------------------------------------
std::string BuiltInFuncsSystem::GetInputForSystemCommand(const std::string& str)
{
#ifdef OS_WIN
    return std::string('"' + str + "\" 2>&1");
#else
    //system('"ls /tmp"')
    //system('ls /tmp')
    //system('"ls" "/tmp"')
    //system('""ls" "/tmp""')
    //system('/users/user1/test.sh')
    //system('/users/user1/test.sh 1')
    //system('"/users/user1/test.sh" 1')
    //system('"<exe> -f getcmdinput_test.oml -input [1 3 4] string "test""');
    //system('"<exe> -input [1 3 4] string "test" -f getcmdinput_test.oml"');

    std::string in(str);
    std::string cmd;
    size_t pos = in.find(" ");
    if (pos != std::string::npos)
    {
        cmd = in.substr(0, pos);
        in = in.substr(pos);
    }
    else
    {
        // No space, so strip all extra quotes, if any
        cmd = in;
        in = "";
    }
    // This is the actual system command that needs to be run, so strip any
    // quotes from that
    
    size_t numLeadingQuotes = cmd.find_first_not_of("\"");
    size_t numTrailingQuotes = 0;
    if (numLeadingQuotes == std::string::npos)
    {
        numLeadingQuotes = 0;
    }
    else
    {
        size_t len1 = cmd.length();
        size_t numquotesremoved = 0;
        std::string tmp;
        for (size_t i = 0; i < len1; ++i)
        {
            char ch = cmd[i];
            if (ch != '"')
            {
                tmp += ch;
            }
            else
            {
                numquotesremoved++;
            }
        }
        cmd = tmp;

        if (numquotesremoved == numLeadingQuotes)
        {
            numTrailingQuotes = numLeadingQuotes;
        }
        else if (numquotesremoved == numLeadingQuotes * 2)
        {
            numTrailingQuotes = 0;
        }
        else if (numquotesremoved > numLeadingQuotes)
        {
            numTrailingQuotes = numquotesremoved - numLeadingQuotes;
        }
        if (numTrailingQuotes > 0 && !in.empty())
        {
            in = in.substr(0, in.length() - numTrailingQuotes);
        }
    }
    
    cmd += in;
    return std::string(cmd + " 2>&1");
#endif
}
//------------------------------------------------------------------------------
// Returns system dependant path separator [filesep]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Filesep(EvaluatorInterface           eval,
    const std::vector<Currency>& inputs,
    std::vector<Currency>& outputs)
{
#ifdef OS_WIN
    std::string sep;
    if (!inputs.empty())
    {
        if (!inputs[0].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 1);
        }
        std::string val(inputs[0].StringVal());
        if (!val.empty())
        {
            std::transform(val.begin(), val.end(), val.begin(), ::tolower);
            if (val != "all")
            {
                throw OML_Error(OML_ERR_BAD_STRING, 1);
            }
            sep += '/';
        }
    }
    sep += '\\';
    outputs.push_back(sep);
#else
    outputs.push_back("/");
#endif
    return true;
}
//------------------------------------------------------------------------------
// Returns true after starting/stopping writing to diary [diary]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::Diary(EvaluatorInterface           eval,
                               const std::vector<Currency>& inputs,
                               std::vector<Currency>&       outputs)
{
    std::ofstream& ofs = eval.GetDiary();
    std::string    name("diary.txt");
    bool           close = false;

    if (inputs.empty())
    {
        close = ofs.is_open();
    }
    else
    {
        const Currency& cur = inputs[0];
        if (cur.IsLogical())
        {
            int val = static_cast<int>(cur.Scalar());
            if (val == 0)
            {
                close = true;
            }
            else if (val != 1)
            {
                throw(OML_ERR_BAD_STRING, 1);
            }
        }
        else if (cur.IsString())
        {
            name = cur.StringVal();
            if (name.empty())
            {
                throw (OML_ERR_NONEMPTY_STR, 1);
            }
            std::string lower(name);
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            if (lower == "off")
            {
                close = true;
            }
            else if (lower == "on")
            {
                name = ""; // Open diary with default name
            }
        }
    }

    if (ofs.is_open())
    {
        ofs.close(); // Force close if diary is open
    }

    if (!close)
    {
        BuiltInFuncsUtils utils;
        if (name.empty())
        {
            name = "diary.txt";
        }
        name = utils.StripMultipleSlashesAndNormalize(name);
        ofs.open(name, std::ios_base::out | std::ios_base::trunc);
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns true after starting/stopping writing to outputlog [outputlog]
//------------------------------------------------------------------------------
bool BuiltInFuncsSystem::OutputLog(EvaluatorInterface           eval,
                                   const std::vector<Currency>& inputs,
                                   std::vector<Currency>& outputs)
{
    std::ofstream& ofs = CurrencyDisplay::GetOutputLog();
    std::wstring    name(L"omloutputlog.txt");
    bool           close = false;

    if (inputs.empty())
    {
        close = ofs.is_open();
    }
    else
    {
        const Currency& cur = inputs[0];
        if (cur.IsLogical())
        {
            int val = static_cast<int>(cur.Scalar());
            if (val == 0)
            {
                close = true;
            }
            else if (val != 1)
            {
                throw(OML_ERR_BAD_STRING, 1);
            }
        }
        else if (cur.IsString())
        {
            std::string val (cur.StringVal());
            if (val.empty())
            {
                throw (OML_ERR_NONEMPTY_STR, 1);
            }
            std::string lower(val);
            std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
            if (lower == "off")
            {
                close = true;
            }
            else if (lower != "on")
            {
                name = BuiltInFuncsUtils::StdString2WString(val);
            }
        }
    }

    if (ofs.is_open())
    {
        ofs.close(); // Force close if outputlog is open
    }

    if (!close)
    {
        if (name.empty())
        {
            name = L"omloutputlog.txt";
        }
        name = BuiltInFuncsUtils::GetNormpathW(name);
        CurrencyDisplay::SetOutputLogName(name);
#ifdef OS_WIN
        ofs.open(name, std::ios_base::out | std::ios_base::trunc);
#else
        std::string logname(BuiltInFuncsUtils::WString2StdString(name));
        ofs.open(logname, std::ios_base::out | std::ios_base::trunc);
#endif
    }
    return true;
}
