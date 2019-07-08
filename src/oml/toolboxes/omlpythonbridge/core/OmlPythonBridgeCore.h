/**
* @file OmlPythonBridgeCore.h
* @date February, 2019
* Copyright (C) 2015-2019 Altair Engineering, Inc.
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
#ifndef OML_PYTHON_BRIDGE_CORE_H__
#define OML_PYTHON_BRIDGE_CORE_H__

#include "OmlPythonBridgeCoreDefs.h"
#include "Currency.h"

#ifndef PyObject_HEAD
struct _object;
typedef _object PyObject;
#endif

#if PY_MAJOR_VERSION >= 3
#  define IS_PY3
#  define PyString_FromString PyUnicode_FromString
#  define PyInt_FromLong PyLong_FromLong
#  define PyInt_Check PyLong_Check
#  define PyInt_AsLong PyLong_AsLong
#  define PyString_Check PyBytes_Check
#endif
//-------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------
//! Class for providing data exchange between Python and OML
//-------------------------------------------------------------------------------------
class OMLPYTHONBRIDGECORE_DECLS OmlPythonBridgeCore
{

public:
    static OmlPythonBridgeCore* GetInstance();
    static void ReleaseInstance();
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
    bool ConvertPyObjectToCurrency(std::vector<Currency>& outputs, PyObject* const& obj);
    //! Convert python variable value in to string representation
    //! \param obj is the python variable name
    //! \return string representation of python variable value
    std::string GetPyObjectAsString(PyObject* const& obj);
    //! Verifies if the given python variable value can be represented by types in oml
    //! \param obj is the python variable name
    //! \return true if python variable value can be represented in oml type, false otherwise
    bool IsTypeSupported(PyObject* const& obj);

    void HandleException(void);
    void SetErrorMessage(const std::string &error);
    std::string GetErrorMessage();

private:
    //! Constructor
    OmlPythonBridgeCore();
    //! Copy Constructor
    OmlPythonBridgeCore(const OmlPythonBridgeCore&) = delete;
    //! Assignment Operator
    OmlPythonBridgeCore& operator=(const OmlPythonBridgeCore&) = delete;
    //! Destructor
    ~OmlPythonBridgeCore();

    static OmlPythonBridgeCore* _instance;
    std::string m_errorMessagePython;
};

#endif