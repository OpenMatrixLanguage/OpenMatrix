/**
* @file OmlPythonBridge.cxx
* @date December, 2017
* Copyright (C) 2015-2018 Altair Engineering, Inc.  
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

/* Remove it if debug version of python is available. */
#if defined(_DEBUG)
#   undef _DEBUG
#   define PYTHON_DEBUG
#endif
/* End change */

#include "Python.h"
#include "marshal.h"
#include "frameobject.h"
#include "arrayobject.h"

/* Remove it if debug version of python is available. */
#if defined(PYTHON_DEBUG)
#   define _DEBUG 1
#endif
/* End change */
#include <vector>
#include "StructData.h"
#include "OmlPythonBridge.h"
#include "EvaluatorInt.h" 
#include <iostream>
#include "BuiltInFuncsUtils.h"
#include "OmlPythonBridgeCore.h"

#ifdef OS_WIN
#include <fcntl.h>
#include <io.h>
#endif
// End defines/includes

OmlPythonBridge* OmlPythonBridge::_instance = NULL;

OmlPythonBridge* OmlPythonBridge::GetInstance()
{
    if (!_instance)
    {
        _instance = new OmlPythonBridge();
    }
    return _instance;
}

void OmlPythonBridge::ReleaseInstance()
{
    delete _instance;
    _instance = NULL;
}

//! Constructor
//! Initialize Python
OmlPythonBridge::OmlPythonBridge(): _threadState(NULL)
{
#ifdef OS_WIN
    int result_stdin, result_stdout, result_stderr;
    result_stdin = _setmode(_fileno(stdin), O_TEXT);
    result_stdout = _setmode(_fileno(stdout), O_TEXT);
    result_stderr = _setmode(_fileno(stderr), O_TEXT);
#endif

    char* python_home = getenv("OML_PYTHONHOME");
    if (python_home)
    {
        std::string python_home_str = python_home;
        static std::wstring widePythonHome = std::wstring(python_home_str.begin(), python_home_str.end());;
        Py_SetPythonHome  (const_cast <wchar_t*>(widePythonHome.c_str()));
    }

    Py_Initialize();
#ifdef OS_WIN
    if (result_stdin != -1)
        _setmode(_fileno(stdin), result_stdin);
    if (result_stdout != -1)
        _setmode(_fileno(stdout), result_stdout);
    if (result_stderr != -1)
        _setmode(_fileno(stderr), result_stderr);
#endif

    if (Py_IsInitialized())
    {
        PyEval_InitThreads();
        InitNumpy();
        _threadState = PyEval_SaveThread();
    }
}

void* OmlPythonBridge::InitNumpy()
{
    if (PyArray_API == NULL)
    {
        import_array();
    }
}

//! Destructor
//! UnInitialize Python
OmlPythonBridge::~OmlPythonBridge()
{
    if (Py_IsInitialized())
    {
        PyEval_RestoreThread((PyThreadState*)_threadState);
        Py_Finalize();
    }
}

bool OmlPythonBridge::SetPythonVariable(const std::string& name, const Currency& value)
{
    bool success = false;
    PyGILState_STATE state = PyGILState_Ensure();
    PyObject* m = PyImport_AddModule("__main__");
    PyObject* valueobj = NULL;
    PyObject* nameobj = NULL;
    int status = -1;

    OmlPythonBridgeCore::GetInstance()->ConvertCurrencyToPyObject(valueobj, value);
    if (valueobj != NULL)
    {
        OmlPythonBridgeCore::GetInstance()->ConvertCurrencyToPyObject(nameobj, Currency(name));
                
        if (nameobj != NULL)
        {
            status = PyObject_SetAttr(m, nameobj, valueobj);
            Py_DECREF(nameobj);
        }

        Py_DECREF(valueobj);
    }
            
    if (status == 0)
    {
        success = true;
    }
    
    OmlPythonBridgeCore::GetInstance()->HandleException();
    PyGILState_Release(state);
    return success;
}

bool OmlPythonBridge::GetPythonVariable(const std::string& name, std::vector<Currency>& outputs)
{
    bool success = false;
    PyGILState_STATE state = PyGILState_Ensure();
    PyObject *m = PyImport_AddModule("__main__");
    PyObject *nameobj = NULL;
    OmlPythonBridgeCore::GetInstance()->ConvertCurrencyToPyObject(nameobj, Currency(name));
            
    if (PyObject_HasAttr(m,nameobj))
    {
        PyObject* valueobj = PyObject_GetAttr(m,nameobj);			
        bool is_type_supported = false;		

        if (valueobj)
        {
            //check for limitations:  1. list to cell 2. dict to struct
            if (PyList_Check(valueobj) || PyDict_Check(valueobj))
            {
                is_type_supported = OmlPythonBridgeCore::GetInstance()->IsTypeSupported(valueobj);
            } 
            else
            {
                is_type_supported = true;
            }
                
            if (is_type_supported)
                success = OmlPythonBridgeCore::GetInstance()->ConvertPyObjectToCurrency(outputs, valueobj);
                    
            Py_DECREF(valueobj);
        }
            
        if (!is_type_supported)
        {
            SetLastErrorMessage("Type not supported.");
        }
    } 
    else
    {
        SetLastErrorMessage("NameError: name '"+name+"' is not defined");
    }

    Py_XDECREF(nameobj);
    PyGILState_Release(state);

    return success;
}

std::string OmlPythonBridge::GetLastErrorMessage()
{
    return OmlPythonBridgeCore::GetInstance()->GetErrorMessage();
}

void OmlPythonBridge::SetLastErrorMessage(std::string error)
{
    OmlPythonBridgeCore::GetInstance()->SetErrorMessage(error);
}

bool OmlPythonBridge::EvalPythonFile(const std::string &python_file)
{
    PyGILState_STATE state = PyGILState_Ensure();

    PyObject* main = PyImport_AddModule("__main__");
    PyObject* globals = PyModule_GetDict(main);
    bool success = RunFile(python_file, globals, globals);

    PyGILState_Release(state);
    return success;
}

bool OmlPythonBridge::EvalPythonScript(const std::string &script)
{
    bool success = false;
    PyGILState_STATE state = PyGILState_Ensure();

    PyObject* globals = PyModule_GetDict(PyImport_AddModule("__main__"));
    if(NULL != PyRun_StringFlags(script.c_str(), Py_file_input, globals, globals, 0))
    {
        success = true;
    }
    OmlPythonBridgeCore::GetInstance()->HandleException();
    
    PyGILState_Release(state);

    return success;
}

bool OmlPythonBridge::RunFile(const std::string& python_file, PyObject* globals, PyObject* locals)
{
    // Caller is responsible for PyGILState_Ensure()
    bool success = false;
    size_t dot_index = python_file.rfind('.');

    if (dot_index == std::string::npos)
    {
        std::string msg("file extension is missing.");
        PyErr_SetString(PyExc_IOError, msg.c_str());
        OmlPythonBridgeCore::GetInstance()->HandleException();
        return success;
    }
    std::string file = BuiltInFuncsUtils::Normpath(python_file);
    std::string ext(BuiltInFuncsUtils::GetFileExtension(file));
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    bool is_pyc = false;

    if ( (0 == ext.compare("py")) && !BuiltInFuncsUtils::FileExists(file) &&
        BuiltInFuncsUtils::FileExists(file.substr(0, dot_index)+".pyc"))
    {
        ext = "pyc";
        is_pyc = true;
        file = file.substr(0, dot_index) + ".pyc";
    }
    else if (0 == ext.compare("pyc"))
    {
        is_pyc = true;
    }

    const char* read_mode = is_pyc ? "rb" : "r";

    // Update the value of __file__. It is restored after file execution.
    PyObject* old__file__ = PyDict_GetItemString(globals, "__file__");
    PyObject* new__file__ = PyUnicode_FromString(file.c_str());
    PyDict_SetItemString(globals, "__file__", new__file__);

    // Python must open FILE, otherwise PyRun_SimpleFile will crash in debug
    // mode on windows (because python uses release mode regardless)
#ifdef IS_PY3
    PyObject* fileobj = PyUnicode_FromString(file.c_str());
    FILE* fp = _Py_fopen_obj(fileobj, read_mode);
    bool closeit = true;
#else
    PyObject* fileobj = PyFile_FromString((char *)file.c_str(), (char *) read_mode);
    FILE* fp = PyFile_AsFile (fileobj);
    bool closeit = false;
#endif

    if(NULL == fp)
    {
        std::string msg("could not open file - "); msg += file;
        PyErr_SetString(PyExc_IOError, msg.c_str());
    }
    else if(is_pyc)
    {
        if (PyMarshal_ReadLongFromFile(fp) != PyImport_GetMagicNumber())
        {
            std::string msg("pyc file is compiled with a different python version - "); msg += file;
            PyErr_SetString(PyExc_IOError, msg.c_str());
        }
        else
        {
            PyMarshal_ReadLongFromFile(fp);
#ifdef IS_PY3
            PyMarshal_ReadLongFromFile(fp);
#endif
            PyObject* filedata = PyMarshal_ReadLastObjectFromFile(fp);
            if((NULL == filedata) || !PyCode_Check(filedata))
            {
                std::string msg ("pyc file does not contain valid code objects - "); msg += file;
                PyErr_SetString(PyExc_IOError, msg.c_str());
            }
            else
            {
#ifdef IS_PY3
                PyObject* code = PyEval_EvalCode(filedata, globals, locals);
#else
                PyObject* code = PyEval_EvalCode((PyCodeObject*) filedata, globals, locals);
#endif
                success = (NULL != code)?(true):(false);
                Py_XDECREF (code);
            }
            Py_XDECREF (filedata);
        }
    }
    else
    {
        PyObject* code = PyRun_FileExFlags(fp, file.c_str(), Py_file_input, globals, locals, closeit, NULL);
        success = (code != NULL)?(true):(false);
        Py_XDECREF(code);
    }

    OmlPythonBridgeCore::GetInstance()->HandleException();
    
    // Restore __file__ to its state prior to this function call.
    if(old__file__ == NULL)
    {
        PyDict_DelItemString(globals, "__file__");
    }
    else
    {
        PyDict_SetItemString(globals, "__file__", old__file__);
    }
    Py_XDECREF (new__file__);

    return success;
}