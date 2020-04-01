/**
* @file BuiltInFuncsCore.h
* @date February 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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

#ifndef __BUILTINFUNCSCORE__
#define __BUILTINFUNCSCORE__

#include "EvaluatorInt.h"

//------------------------------------------------------------------------------
//!
//! \brief Utility class for built-in functions implementing core commands
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS BuiltInFuncsCore
{
public:
    //!
    //! Destructor
    //!
    ~BuiltInFuncsCore() {}
    //!
    //! Returns true after adding a toolbox [addtoolbox command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool AddToolbox( EvaluatorInterface           eval, 
                            const std::vector<Currency>& inputs, 
                            std::vector<Currency>&       outputs);
    //!
    //! Returns true after getting class information [class command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Class( EvaluatorInterface           eval,
                       const std::vector<Currency>& inputs,
                       std::vector<Currency>&       outputs);

	//!
    //! Returns true after determining if the input is of the specified class/type [class command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool IsA( EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);

    //!
    //! Returns true after specifying format options [format command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Format( EvaluatorInterface           eval,
                        const std::vector<Currency>& inputs,
                        std::vector<Currency>&       outputs);
    //!
    //! Returns true after getting user input. [input command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs - string to be displayed to user
    //! \param outputs Vector of outputs - string returned by user
    //!
    static bool Input( EvaluatorInterface           eval, 
	                   const std::vector<Currency>& inputs, 
			           std::vector<Currency>&       outputs);
    //!
    //! Returns true after pausing a script [pause command]
    //! \param[in] eval     Evaluator interface
    //! \param[in] inputs   Vector of inputs - string to be displayed to user
    //! \param[out] outputs Vector of outputs - string returned by user
    //!
    static bool Pause( EvaluatorInterface           eval, 
	                   const std::vector<Currency>& inputs, 
			           std::vector<Currency>&       outputs);
    //!
    //! Returns true, gets information of variables defined in current scope [whos command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Whos( EvaluatorInterface           eval, 
	                   const std::vector<Currency>& inputs, 
			           std::vector<Currency>&       outputs);
    //!
    //! Returns true after displaying a count of functions [funccount command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Funccount( EvaluatorInterface           eval, 
                           const std::vector<Currency>& inputs, 
                           std::vector<Currency>&       outputs);
    //!
    //! Returns true after displaying a list of functions along with total count [funclist command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Funclist( EvaluatorInterface           eval, 
                          const std::vector<Currency>& inputs, 
                          std::vector<Currency>&       outputs);
    //!
    //! Returns true after displaying a list of variables in the interpreter [varlist command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Varlist( EvaluatorInterface           eval, 
                         const std::vector<Currency>& inputs, 
                         std::vector<Currency>&       outputs);

	//!
	//! Returns true after displaying a list of keywords in the interpreter [varlist command]
	//! \param eval    Evaluator interface
	//! \param inputs  Vector of inputs
	//! \param outputs Vector of outputs
	//!
	static bool Keywordlist(EvaluatorInterface           eval,
		const std::vector<Currency>& inputs,
		std::vector<Currency>&       outputs);
    //!
    //! Returns true after getting current file name being processed [omlfilename command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Omlfilename( EvaluatorInterface           eval, 
                             const std::vector<Currency>& inputs, 
                             std::vector<Currency>&       outputs);
    //!
    //! Returns true after checking if given variable/file/path/function exists [exist command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Exist( EvaluatorInterface           eval, 
                       const std::vector<Currency>& inputs, 
                       std::vector<Currency>&       outputs);
    //!
    //! Returns true and returns a description of type of the given input [type command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Type( EvaluatorInterface           eval, 
                      const std::vector<Currency>& inputs, 
                      std::vector<Currency>&       outputs);
    //!
    //! Gets build from the given file
    //! \param versionfile Path to the version file
    //!
    static std::string GetBuildNumber( const std::string& versionfile);
    //!
    //! Enables/disables pagination (omlpaginate)
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool OmlPaginate(EvaluatorInterface           eval, 
                            const std::vector<Currency>& inputs, 
                            std::vector<Currency>&       outputs);
    //!
    //! Sleeps for the given time period (sleep)
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool OmlSleep(EvaluatorInterface           eval, 
                         const std::vector<Currency>& inputs, 
                         std::vector<Currency>&       outputs);
    //!
    //! Dynamically loads a library and returns handle to it
    //! \param name Library name
    //!
    static void* DyLoadLibrary(const std::string& name);
    //! Returns function pointer from a dynamically loaded library
    //! \param handle Handle to dynamically loaded library
    //! \param name   Function name
    //!
    static void* DyGetFunction(void*              handle,
                               const std::string& name);  
    //!
    //! Releases dynamically loaded library
    //! \param handle Handle to dynamically loaded library
    //!
    static void DyFreeLibrary(void* handle);

	//!
	//! Gets the path of a loaded library
	//! \param handle Handle to dynamically loaded library
	//!
	static std::string DyGetLibraryPath(void* handle);

private:
    //!
    //! Constructor
    //!
    BuiltInFuncsCore() {}
    //!
    //! Sleeps for the given amount of seconds
    //! \param seconds Seconds to sleep
    //!
    void SleepFunc( double seconds) const;
    //!
    //! Converts dimensions to string - whos helper method
    //! \param dims Dimensions
    //!
    std::string DimensionToString( const std::vector<int>& dims) const;
    //!
    //! Gets class for given currency
    //! \param cur Given currency
    //!
    std::string GetCurrencyClass( const Currency& cur) const;
    //!
    //! Gets number of bytes used by Currency
    //! \param cur Given currency
    //!
    int GetsBytesUsed( const Currency& cur) const;
};

#endif // __BUILTINFUNCSCORE__
