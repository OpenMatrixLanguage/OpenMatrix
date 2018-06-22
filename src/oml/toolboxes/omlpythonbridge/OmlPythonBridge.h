/////////////////////////////////////////////////////////////////////////
 // File : OmlPythonBridge.h
 // Copyright (c) 2017 solidThinking Inc.  All Rights Reserved
 // Contains trade secrets of solidThinking, Inc.  Copyright notice
 // does not imply publication.  Decompilation or disassembly of this
 // software is strictly prohibited.
 ////////////////////////////////////////////////////////////////////////
#ifndef OML_PYTHON_BRIDGE_H__
#define OML_PYTHON_BRIDGE_H__

// Begin defines/includes

#include <string>
#include "Currency.h"
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

    static OmlPythonBridge* GetInstance();
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

    void ClearLastErrorMessage();
    std::string GetLastErrorMessage();

    //! Evaluate python file in the __main__ scope
    //! \param python_file is python file to evaluate
    //! \return true if evaluation attempt did not produce errors, false otherwise
    bool EvalPythonFile(const std::string &python_file);

    //! Evaluate script text in the __main__ scope
    //! \param script python code to evaluate
    //! \return true if script evaluation attempt did not produce errors, false otherwise
    bool EvalPythonScript(const std::string &script);
private:
    
    //! Constructor
    OmlPythonBridge();
    //! Copy Constructor
    OmlPythonBridge( const OmlPythonBridge& ) = delete;
    //! Assignment Operator
    OmlPythonBridge& operator=( const OmlPythonBridge& ) = delete;
    //! Destructor
    ~OmlPythonBridge();

    //! Retrieve error message from python interpreter and stores it in internal variable
    void HandleException(void);

    //! Verifies if the given python variable value can be represented by types in oml
    //! \param obj is the python variable name
    //! \return true if python variable value can be represented in oml type, false otherwise
    bool IsTypeSupported(PyObject* const& obj);

    //! Convert python variable value in to string representation
    //! \param obj is the python variable name
    //! \return string representation of python variable value
    std::string GetPyObjectAsString(PyObject* const& obj);
    //! Convert OML variable value into data type in Python
    //! OMLDataType     PythonDataType Limitation
    //! Logical         Bool
    //! Number          Float
    //! Complex         Complex
    //! String          String
    //! Cell(1, n)      List           Cell with one dimension only supported.
    //! n:number of elements in list
    //! Struct        dict             Struct with one dimension only supported.
    //! Matrix        Numpy - matrix   Data types supported in OML : matrix Bool, Int, long, Float, Complex.
    //! ND Matrix     Numpy - Ndarray  Data types supported in OML : matrix Bool, Int, long, Float, Complex.
    //! \param[in] obj is python variable name
    //! \param[in] var is oml variable name
    //! \return true if oml variable value can be represented in python
    bool ConvertCurrencyToPyObject(PyObject*& obj, const Currency& var);

    //! Convert Python variable value into data type in OML
    //! PythonDataType  OMLDataType      Limitation
    //! Bool            Logical
    //! long, Float     Number
    //! Complex         Complex
    //! List            Cell(1, n)    Does not support if list contains Dict(with limitation), Tupple, Set
    //! n:number of elements in list
    //! dict            Struct        Supports only if keys in dict are string or char.
    //! Numpy - array   Matrix        Data types supported in OML : matrix, Bool, Ing, long, Float, Complex.
    //! Numpy - matrix  Matrix        Data types supported in OML : matrix, Bool, Ing, long, Float, Complex.
    //! Numpy Ndarray   ND Matrix     Data types supported in OML : matrix Bool, Int, long, Float, Complex.
    //! \param[out] outputs python variable value in OML type
    //! \param[in] obj is python variable name
    //! \return true if python variable value can be represented in OML
    bool ConvertPyObjectToCurrency (std::vector<Currency>& outputs, PyObject* const& obj);

    //! Evaluate python file
    //! \param python_file is python file to evaluate
    //! \param globals is dictionay of global variables
    //! \param locals is dictionay of local variables
    //! \return true if evaluation attempt did not produce errors, false otherwise
    bool RunFile(const std::string& python_file, PyObject* globals, PyObject* locals);
    //! Initializes numpy
    void* InitNumpy();

    static OmlPythonBridge* _instance;
    void* _threadState;
    std::string _errorMessage;
};

#endif  OML_PYTHON_BRIDGE_H__