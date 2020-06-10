/**
* @file OmlPythonBridge.h
* @date December, 2017
* Copyright (C) 2015-2020 Altair Engineering, Inc.  
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
#ifndef OML_PYTHON_BRIDGE_H__
#define OML_PYTHON_BRIDGE_H__

// Begin defines/includes

#include <string>
#include "Currency.h"
#include "EvaluatorInt.h" 
#ifndef Py_PYTHON_H
#define PyObject void
#endif
// End defines/includes

//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//! Class for providing Python initialization and OML Python interaction
//-------------------------------------------------------------------------------------
class OmlPythonBridge
{
public:
    //!
    //! Gets instance
    //! \param eval Evaluator interface
    //!
    static OmlPythonBridge* GetInstance(EvaluatorInterface eval);
    static void ReleaseInstance();
    //! Creates/Updates Python variable in the __main__ scope and sets the value
    //! \param name is name of python variable
    //! \param value is value to be set to Python variable
    //! \return true if setting value to Python variable value did not produce errors, false otherwise
    bool SetPythonVariable(const std::string& name, const Currency& value);

    //! Retrieve Python variable value from __main__ scope
    //! \param name is name of python variable value to be retrieve
    //! \param outputs is value of Python variable
    //! \return true if retrieve Python variable value did not produce errors, false otherwise
    bool GetPythonVariable(const std::string& name, std::vector<Currency>& outputs);

    std::string GetLastErrorMessage();
    void SetLastErrorMessage(std::string error);

    //! Evaluate python file in the __main__ scope
    //! \param python_file is python file to evaluate
    //! \return true if evaluation attempt did not produce errors, false otherwise
    bool EvalPythonFile(const std::string &python_file);

    //! Evaluate script text in the __main__ scope
    //! \param script python code to evaluate
    //! \return true if script evaluation attempt did not produce errors, false otherwise
    bool EvalPythonScript(const std::string &script);
private:
    //!
    //! Constructor
    //! \param eval Evaluator interface
    //!
    OmlPythonBridge(EvaluatorInterface eval);
    //! Copy Constructor
    OmlPythonBridge( const OmlPythonBridge& ) = delete;
    //! Assignment Operator
    OmlPythonBridge& operator=( const OmlPythonBridge& ) = delete;
    //! Destructor
    ~OmlPythonBridge();

    //! Evaluate python file
    //! \param python_file is python file to evaluate
    //! \param globals is dictionay of global variables
    //! \param locals is dictionay of local variables
    //! \return true if evaluation attempt did not produce errors, false otherwise
    bool RunFile(const std::string& python_file, PyObject* globals, PyObject* locals);
    //! Initializes numpy
    void* InitNumpy();
    //!
    //! Initializes command line arguments
    //! \param eval Evaluator interface
    //!
    void InitializeCommandLineArgs(EvaluatorInterface eval);

    static OmlPythonBridge* _instance;
    void* _threadState;
};

#endif  OML_PYTHON_BRIDGE_H__