/**
* @file BuiltInFuncsSystem.h
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

#ifndef __BUILTINFUNCSSYSTEM__
#define __BUILTINFUNCSSYSTEM__

#include "EvaluatorInt.h"
//------------------------------------------------------------------------------
//!
//! \brief Class for built-in functions implementing system/time commands
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS BuiltInFuncsSystem
{
public:
    //!
    //! Destructor
    //!
    ~BuiltInFuncsSystem() {}
    //!
    //! Returns true after listing directory contents [dir command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! \todo implement a return value on *nix -- need to return a struct
    //!
    static bool Dir( EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns true after listing contents of a directory [ls command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Ls( EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs, 
                    std::vector<Currency>&       outputs);    
    //!
    //! Returns true after running a system command [system command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool System(EvaluatorInterface           eval,
                       const std::vector<Currency>& inputs, 
                       std::vector<Currency>&       outputs);   
    //!
    //! Returns true after running a unix system command [unix command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Unix(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs, 
                     std::vector<Currency>&       outputs);    
    //!
    //! Returns true after deleting given file(s) [delete command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Delete(EvaluatorInterface           eval,
                       const std::vector<Currency>& inputs, 
                       std::vector<Currency>&       outputs);    

    // Client specific environment functions
    //!
    //! Returns true and clones given environment handle [cloneenv]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool CloneEnv(EvaluatorInterface           eval, 
                         const std::vector<Currency>& inputs, 
                         std::vector<Currency>&       outputs);
    //!
    //! Returns true and gets base environment handle [getbaseenv]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool GetBaseEnv(EvaluatorInterface           eval, 
                           const std::vector<Currency>& inputs, 
                           std::vector<Currency>&       outputs);
    //!
    //! Returns true and gets current environment handle [getcurrentenv]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool GetCurrentEnv(EvaluatorInterface           eval, 
                              const std::vector<Currency>& inputs, 
                              std::vector<Currency>&       outputs);
    //!
    //! Returns true and gets new environment handle [getnewenv]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool GetNewEnv(EvaluatorInterface           eval, 
                          const std::vector<Currency>& inputs, 
                          std::vector<Currency>&       outputs);
    //!
    //! Returns true and imports environment handle [importenv]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool ImportEnv(EvaluatorInterface           eval, 
                          const std::vector<Currency>& inputs, 
                          std::vector<Currency>&       outputs);
    //!
    //! Returns true and imports given environment handle [importenvin]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool ImportEnvIn(EvaluatorInterface           eval, 
                            const std::vector<Currency>& inputs, 
                            std::vector<Currency>&       outputs);

    //!
    //! Returns true and gets environment variable value [getenvvalue]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool GetEnvVal(EvaluatorInterface           eval, 
                          const std::vector<Currency>& inputs, 
                          std::vector<Currency>&       outputs);
    //!
    //! Returns true and sets environment variable value [setenvvalue]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool SetEnvVal(EvaluatorInterface           eval, 
                          const std::vector<Currency>& inputs, 
                          std::vector<Currency>&       outputs);
    //!
    //! Returns true and clears environment variable value [clearenvvalue]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool ClearEnvVal(EvaluatorInterface           eval, 
                            const std::vector<Currency>& inputs, 
                            std::vector<Currency>&       outputs);
    //!
    //! Returns true after changing current directory [cd]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! \todo implement a return value on *nix -- need to return a struct
    //!
    static bool Cd(EvaluatorInterface           eval,
                   const std::vector<Currency>& inputs,
                   std::vector<Currency>&       outputs);
    //!
    //! Returns true after deleing given directory [rmdir]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //! \todo implement a return value on *nix -- need to return a struct
    //!
    static bool Rmdir(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Gives the current time in seconds since Jan 1, 1970 [time]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Time(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Gives the current time as a string [ctime]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool CTime(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Gets a path constructed from the given dir/sub-directories [genpath]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Genpath(EvaluatorInterface           eval,
                        const std::vector<Currency>& inputs,
                        std::vector<Currency>&       outputs);
    //!
    //! Gets process id [getpid]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool GetPid(EvaluatorInterface           eval,
                       const std::vector<Currency>& inputs,
                       std::vector<Currency>&       outputs);
    //!
    //! Creates an absolute path name [make_absolute_filename]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool MakeAbsFilename(EvaluatorInterface           eval,
                                const std::vector<Currency>& inputs,
                                std::vector<Currency>&       outputs);
private:
    //!
    //! Constructor
    //!
    BuiltInFuncsSystem() {}
    //!
    //! Gets info whether the given path is an existing directory or file
    //! \param path Given path
    //! \param isdir  True if given path refers to a directory
    //! \param isfile True if given path refers to a file
    //!
    void IsDirectoryOrFile( const std::string& path,
                            bool&              isdir,
                            bool&              isfile) const;
    //!
    //! Gets timestamp as string
    //! \param rawtime Time
    //!
    std::string GetTimeString( time_t rawtime);

#ifdef OS_WIN
    //!
    //! Returns true if path is a directory, supports unicode on Windows
    //! \param path Given path
    //!
    bool IsDir(std::wstring& path);
    //!
    //! Lists directories, recursively; helper for genpath. Supports unicode
    //! \param parent     Parent dir
    //! \param exclude    Dirs to ignore, normalized and in lower case
    //! \param paths      String with all the valid paths, separated by ';'
    //!
    void ListDirW(const std::wstring&                                 parent,
                  const std::vector< std::pair<std::wstring, bool> >& exclude,
                  std::wstring&                                       paths);
    //!
    //! Helper method to set locale and output mode to support Unicode
    //!
    void SetToUtf8Locale();
#else
    //!
    //! Lists directories, recursively; helper for genpath.
    //! \param parent     Parent dir
    //! \param ignorepath Dirs to ignore, normalized and in lower case
    //! \param paths      String with all the valid paths, separated by ':'
    //!
    void ListDir(const std::string&                                 parent,
                 const std::vector< std::pair<std::string, bool> >& exclude,
                 std::string&                                       paths);

#endif
};
#endif // __BUILTINFUNCSSYSTEM__
