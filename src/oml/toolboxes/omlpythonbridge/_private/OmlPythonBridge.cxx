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

#if PY_MAJOR_VERSION >= 3
#   define IS_PY3
#   define PyString_FromString PyUnicode_FromString
#   define PyInt_FromLong PyLong_FromLong
#   define PyInt_Check PyLong_Check
#   define PyInt_AsLong PyLong_AsLong
#   define PyString_Check PyBytes_Check
#endif
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
OmlPythonBridge::OmlPythonBridge(): _threadState(NULL),_errorMessage("")
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

    ConvertCurrencyToPyObject(valueobj, value);
    if (valueobj != NULL)
    {
        ConvertCurrencyToPyObject(nameobj, Currency(name));
                
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
    
    HandleException();
    PyGILState_Release(state);
    return success;
}

bool OmlPythonBridge::GetPythonVariable(const std::string& name, std::vector<Currency>& outputs)
{
    bool success = false;
    PyGILState_STATE state = PyGILState_Ensure();
    PyObject *m = PyImport_AddModule("__main__");
    PyObject *nameobj = NULL;
    ConvertCurrencyToPyObject(nameobj, Currency(name));
            
    if (PyObject_HasAttr(m,nameobj))
    {
        PyObject* valueobj = PyObject_GetAttr(m,nameobj);			
        bool is_type_supported = false;		

        if (valueobj)
        {
            //check for limitations:  1. list to cell 2. dict to struct
            if (PyList_Check(valueobj) || PyDict_Check(valueobj))
            {
                is_type_supported = IsTypeSupported(valueobj);
            } 
            else
            {
                is_type_supported = true;
            }
                
            if (is_type_supported)
                success = ConvertPyObjectToCurrency(outputs, valueobj);
                    
            Py_DECREF(valueobj);
        }
            
        if (!is_type_supported)
        {
            _errorMessage = "Type not supported.";
        }
    } 
    else
    {
        _errorMessage = "NameError: name '"+name+"' is not defined";
    }

    Py_XDECREF(nameobj);
    PyGILState_Release(state);

    return success;
}

void OmlPythonBridge::ClearLastErrorMessage()
{
    _errorMessage = "";
}

std::string OmlPythonBridge::GetLastErrorMessage()
{
    return _errorMessage;
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
    HandleException();
    
    PyGILState_Release(state);

    return success;
}

void OmlPythonBridge::HandleException(void)
{
    PyObject * expn = NULL, * xval = NULL, * xtrc = NULL;
    PyErr_Fetch(&expn, &xval, &xtrc);

    if (NULL != expn)
    {
        _errorMessage = "Unknown Exception";

        if (NULL != xval)
        { 
            _errorMessage = GetPyObjectAsString(xval);
        }

        PyErr_Restore(expn, xval, xtrc);

        if (PyErr_Occurred())
        {
            PyErr_Clear();
        }
    }
}

bool OmlPythonBridge::IsTypeSupported(PyObject* const& obj)
{
    bool is_supported = false;

    if (PyString_Check(obj) || PyFloat_Check(obj) || PyLong_Check(obj) || PyBool_Check(obj) || PyComplex_Check(obj))
    {
        is_supported = true;
    }
    else if (PyDict_Check(obj))
    {
        PyObject* keys = PyDict_Keys(obj);
        bool areKeysValid = true;
        
        if (keys)
        {
            Py_ssize_t size = PyList_GET_SIZE(keys);

            for (Py_ssize_t index = 0; index < size; index++)
            {
                PyObject* key = PyList_GetItem(keys, index);

                if (!(PyUnicode_Check(key) || PyString_Check(key)))
                {
                    areKeysValid = false;
                    break;
                }
            }
        }

        Py_XDECREF (keys);

        if (areKeysValid)
        {
            PyObject* values = PyDict_Values(obj);
            is_supported = true;

            if (values)
            {
                Py_ssize_t size = PyList_GET_SIZE(values);

                for (Py_ssize_t index = 0; index < size; index++)
                {
                    PyObject* value = PyList_GetItem(values, index);
                
                    if (!IsTypeSupported(value))
                    {
                        is_supported = false;
                        break;
                    }
                }
            }
            Py_XDECREF (values);
        }
    }
    else if (PyList_Check(obj))
    {
        Py_ssize_t size = PyList_GET_SIZE(obj);
        is_supported = true;

        for (Py_ssize_t index = 0; index < size; index++)
        {
            PyObject* value = PyList_GetItem(obj, index);
        
            if (!IsTypeSupported(value))
            {
                is_supported = false;
                break;
            }
        }
    }
    else if (PyUnicode_Check(obj))
    {
        PyObject *bytesObj = PyUnicode_AsASCIIString(obj);

        if (bytesObj)
        {
            is_supported = true;
        }
        Py_XDECREF (bytesObj);
    }
#if PY_MAJOR_VERSION < 3
    else if (PyInt_Check(obj))
    {
        is_supported = true;
    }
#endif
    else {
        if (PyArray_Check(obj))
        {
            is_supported = true;
        }
    }

    return is_supported;
}

std::string OmlPythonBridge::GetPyObjectAsString(PyObject* const& obj)
{
    std::string value = "";
    PyObject* strobj = PyObject_Str(obj);

    if (NULL != strobj)
    {
        #ifdef IS_PY3
            PyObject* bytesObj = PyUnicode_AsUTF8String(strobj);
            if (bytesObj)
            {
                value = PyBytes_AS_STRING(bytesObj);
            }
            Py_XDECREF (bytesObj); 
        #else
            value = PyString_AsString(strobj);
        #endif
    }
        
    Py_XDECREF(strobj);

    return value;
}

bool OmlPythonBridge::ConvertCurrencyToPyObject(PyObject*& obj, const Currency& var)
{
    if (var.IsScalar())
    {
        if (var.IsLogical())
        {
            obj = PyBool_FromLong(var.Scalar());
        }
        //ToDo: uncomment following code once oml interpreter supports true integer
        //else if (var.IsInteger())
        //{
        //	obj = PyInt_FromLong(var.Scalar());
        //}
        else
        {
            obj = PyFloat_FromDouble(var.Scalar());
        }
    }
    else if (var.IsString())
    {
        obj = PyString_FromString((var.StringVal()).c_str());
    }
    else if (var.IsComplex())
    {
        double real = var.Real();
        double imag = var.Imag();
        obj = PyComplex_FromDoubles(real,imag);
    }
    else if (var.IsMatrix())
    {
        const hwMatrix* data = var.Matrix();
        npy_intp dims[2] = {data->M(),data->N()};
        int nd = 2;
        
        if (data->IsReal())
        {
            obj = PyArray_ZEROS (nd, dims, NPY_DOUBLE, 1);
            double* dptr = (npy_double *) PyArray_DATA(obj);

            for (int i = 0; i < data->Size(); i++)
            {
                dptr[i] = (*data)(i);
            }
        }
        else
        {
            const hwTComplex<double>* cmpxdata = data->GetComplexData();
            obj = PyArray_ZEROS (nd, dims, NPY_COMPLEX128, 1);
            npy_complex128* outdata = (npy_complex128 *) PyArray_DATA(obj);

            for (int i = 0; i < data->Size(); i++) 
            {
                outdata[i].real = (cmpxdata[i]).Real();
                outdata[i].imag = (cmpxdata[i]).Imag();
            }
        }
    }
    else if (var.IsFunctionHandle())
    {
        //ToDo:implement when in need
    }
    else if (var.IsStruct())
    {
        const StructData* data = var.Struct();
        if ((1 == data->M()) && (1 == data->N()))
        {
            PyObject* dict = PyDict_New();
            bool status = true;
            
            if (NULL != dict)
            {
                const std::map<std::string, int>& fieldnames =  data->GetFieldNames();
                
                for (std::map<std::string, int>::const_iterator it=fieldnames.begin(); it!=fieldnames.end(); ++it)
                {
                    const Currency& val = data->GetValue(0, 0, it->first);
                    PyObject* value_obj = NULL;
                    status = ConvertCurrencyToPyObject (value_obj,val);

                    if (!status)
                    {
                        Py_XDECREF(dict);
                        break;
                    }

                    if (-1 == PyDict_SetItemString(dict, (it->first).c_str(), value_obj))
                    {
                        Py_XDECREF(dict);
                        status = false;
                        break;
                    }
                }
            }

            if (status)
            {
                obj = dict;
            }
        }
        else
        {
            _errorMessage = "Multi dimension struct export is not supported.";
        }
        //ToDo:implement struct array exporting when in need
    }
    else if (var.IsCellArray())
    {
        const HML_CELLARRAY* cell = var.CellArray();
        //cell of size 1,n only exported to list
        if (1 == cell->M() )
        {
            bool status = true;
            int size = cell->N();
            PyObject* lst = PyList_New(size);
            
            if (NULL != lst)
            {
                for (int index = 0; index < size; index++)
                {
                    PyObject* valueObj = NULL;
                    const Currency& val = (*cell)(0, index);
                    status = ConvertCurrencyToPyObject (valueObj,val);
                
                    if (!status)
                    {
                        Py_XDECREF(lst);
                        break;
                    }
                    
                    if ( -1 == PyList_SetItem(lst, index, valueObj))
                    {
                        Py_XDECREF(lst);
                        status = false;
                        break;
                    }
                }
            
                if (status)
                {
                    obj = lst;
                }
            }
        }
        else if (0 == cell->M() )
        {
            PyObject* lst = PyList_New(0);
            
            if (NULL != lst)
            {
                obj = lst;
            }
        }
        else
        {
            _errorMessage = "Multi dimension cell export is not supported.";
        }

        //ToDo:implement when in need
    }
    else if (var.IsNDMatrix()) 
    {
        const hwMatrixN* data = var.MatrixN();		
        std::vector<int> dimensions = data->Dimensions();
        npy_intp* dims = new npy_intp[dimensions.size()];		
        size_t nd = dimensions.size();

        for (int i = 0; i < nd; i++)
        {
            dims[i] = dimensions[i];
        }

        if (data->IsReal())
        {
            obj = PyArray_ZEROS (nd, dims, NPY_DOUBLE, 1);
            double* dptr = (npy_double *) PyArray_DATA(obj);

            for (int i = 0; i < data->Size(); i++)
            {
                dptr[i] = (*data)(i);
            }
        }
        else
        {
            obj = PyArray_ZEROS (nd, dims, NPY_COMPLEX128, 1);
            npy_complex128* outdata = (npy_complex128 *) PyArray_DATA(obj);

            for (int i = 0; i < data->Size(); i++) 
            {
                outdata[i].real = (data->z(i)).Real();
                outdata[i].imag = (data->z(i)).Imag();
            }
        }
    }
    else if (var.IsBoundObject())
    {
        //ToDo:implement when in need
    }

    return (NULL != obj)?(true):(false);
}

bool OmlPythonBridge::ConvertPyObjectToCurrency (std::vector<Currency>& outputs, PyObject* const& obj)
{
    bool success = true;	
    
    if (PyDict_Check(obj))
    {
        bool status = true;
        StructData*  dict = EvaluatorInterface::allocateStruct();
        PyObject* keys = PyDict_Keys(obj);
        dict->Dimension(0,0);

        if (keys)
        {
            Py_ssize_t size = PyList_GET_SIZE(keys);
            
            for (Py_ssize_t index = 0; index < size; index++)
            {
                PyObject* key = PyList_GetItem(keys, index);
                std::vector<Currency> field;
                std::vector<Currency> value;
                status = ConvertPyObjectToCurrency(field,key);
                
                if (!status)
                {
                    status = false;
                    break;
                }

                PyObject* value_obj = PyDict_GetItem(obj, key);
                status = ConvertPyObjectToCurrency(value, value_obj);

                if (!status)
                {
                    status = false;
                    break;
                }

                if(field[0].IsCharacter())
                {
                    dict->SetValue(0,0,field[0].StringVal(),value[0]);
                } 
                else if (field[0].IsString())
                {
                    dict->SetValue(0,0,field[0].StringVal(),value[0]);
                }
                else
                {
                    status = false;
                    break;
                }
            }
        }

        Py_XDECREF (keys);

        if(status)
        {
            outputs.push_back(dict);
        }
        else
        {
            success = false;
            outputs.push_back("");
            delete dict;
        }
    }
    else if (PyList_Check(obj))
    {
        Py_ssize_t size = PyList_GET_SIZE(obj);
        HML_CELLARRAY* list = EvaluatorInterface::allocateCellArray(1,size);
        bool status = true;
        for (Py_ssize_t index = 0; index < size; index++)
        {
            std::vector<Currency> elements;
            PyObject* value = PyList_GetItem(obj, index);
            status = ConvertPyObjectToCurrency(elements,value);

            if (!status)
            {
                break;
            }
            (*list)(0,   index) = elements[0];
        }
        
        if (status)
        {
            outputs.push_back(list);
        }
        else
        {
            success = false;
            outputs.push_back("");
            delete list;
        }
    }
    else if (PyTuple_Check(obj) || PyAnySet_Check(obj))
    {
        //ToDo:Implement once find matching data type available in OML
        success = false;
        outputs.push_back("");
    }
    else if (PyBool_Check(obj))
    {
        outputs.push_back(PyObject_IsTrue(obj) == 1);
    }
    else if (PyUnicode_Check(obj))
    {
        PyObject* bytesObj = PyUnicode_AsASCIIString(obj);

        if (!bytesObj)
        {
            success = false;
            outputs.push_back("");
        }
        else
        {
            outputs.push_back(PyBytes_AsString(bytesObj));
        }

        Py_XDECREF (bytesObj);
    }
    else if (PyFloat_Check(obj))
    {
        outputs.push_back((double)PyFloat_AsDouble(obj));
    }	
    else if (PyLong_Check(obj))
    {
        outputs.push_back((double)PyLong_AsLong(obj));
    }
#if PY_MAJOR_VERSION < 3
    else if (PyInt_Check(obj))
    {
        outputs.push_back((int)PyInt_AsLong(obj));
    }
#endif
    else if (PyComplex_Check(obj))
    {
        hwTComplex<double> cplx;
        cplx.Set(PyComplex_RealAsDouble(obj),PyComplex_ImagAsDouble(obj));
        outputs.push_back(cplx);
    }
    else if (PyString_Check(obj))
    {
        outputs.push_back(PyBytes_AsString(obj));
    }
    else 
    { 
        if (PyArray_Check(obj))
        {
            int nd = PyArray_NDIM((PyArrayObject *)obj);
            npy_intp* dimenstions = PyArray_DIMS((PyArrayObject *)obj);
            int size = PyArray_SIZE(obj);
            std::vector<int> dims;
            bool hasZeroDimSize = false;

            if (nd > 1)
            {
                for (int i = 0; i < nd; i++)
                {
                    dims.push_back(dimenstions[i]);
                    if (0 == dimenstions[i])
                        hasZeroDimSize = true;
                }
            }
            else if (1 == nd)
            {
                dims.push_back(1);
                dims.push_back(dimenstions[0]);

                if (0 == dimenstions[0])
                    hasZeroDimSize = true;
            }

            int type = PyArray_TYPE((PyArrayObject *)obj);
            hwMatrixN* matN =  EvaluatorInterface::allocateMatrixN();
            NpyIter* iter = NULL;
            
            if (!hasZeroDimSize)
                iter = NpyIter_New ((PyArrayObject *)obj,NPY_ITER_READONLY | NPY_ITER_ALIGNED, NPY_FORTRANORDER,NPY_NO_CASTING, NULL);
            
            if (NULL != iter)
            {
                NpyIter_IterNextFunc* iternext = NpyIter_GetIterNext(iter, NULL);
                
                if (NULL != iternext)
                {
                    int i = 0;
                    char** dataptr =  NpyIter_GetDataPtrArray(iter);

                    if (PyArray_ISCOMPLEX((PyArrayObject *)obj))
                    {
                        matN->Dimension(dims,hwMatrixN::COMPLEX);
                        PyObject* obj_itr = PyArray_IterNew(obj);
                        
                        do
                        {
                            char* data = *dataptr;

                            switch (type)
                            {
                                case NPY_CFLOAT:
                                {
                                    npy_cfloat* outdata = (npy_cfloat*) data;
                                    matN->z(i) = hwTComplex<double> ((*outdata).real,(*outdata).imag);
                                    break;
                                }

                                case NPY_CDOUBLE:
                                {
                                    npy_cdouble* outdata = (npy_cdouble*) data;
                                    matN->z(i) = hwTComplex<double> ((*outdata).real,(*outdata).imag);
                                    break;
                                }

                                case NPY_CLONGDOUBLE:
                                {
                                    npy_clongdouble*outdata = (npy_clongdouble*) data;
                                    matN->z(i) = hwTComplex<double> ((*outdata).real,(*outdata).imag);
                                    break;
                                }
                                default:
                                {
                                    //ToDo:implement when in need
                                }
                            }
                            
                            ++i;

                        } while (iternext(iter));

                        outputs.push_back(matN);
                    }
                    else if (PyArray_ISNUMBER((PyArrayObject *)obj))
                    {
                        matN->Dimension(dims,hwMatrixN::REAL);

                        do
                        {
                            char* data = *dataptr;

                            switch (type)
                            {
                                case NPY_BOOL:
                                {
                                    npy_bool* outdata = (npy_bool*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_BYTE:
                                {
                                    npy_byte* outdata = (npy_byte*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_UBYTE:
                                {
                                    npy_ubyte* outdata = (npy_ubyte*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_SHORT:
                                {
                                    npy_short* outdata = (npy_short*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_USHORT:
                                {
                                    npy_ushort* outdata = (npy_ushort*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_INT:
                                {
                                    npy_int* outdata = (npy_int*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_UINT:
                                {
                                    npy_uint* outdata = (npy_uint*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_LONG:
                                {
                                    npy_long* outdata = (npy_long*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_ULONG:
                                {
                                    npy_ulong* outdata = (npy_ulong*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_LONGLONG:
                                {
                                    npy_longlong* outdata = (npy_longlong*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_ULONGLONG:
                                {
                                    npy_ulonglong* outdata = (npy_ulonglong*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_FLOAT:
                                {
                                    npy_float* outdata = (npy_float*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_DOUBLE:
                                {
                                    npy_double* outdata = (npy_double*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                case NPY_LONGDOUBLE:
                                {
                                    npy_longdouble* outdata = (npy_longdouble*) data;
                                    (*matN)(i) = (*outdata);
                                    break;
                                }
                                default:
                                {
                                    //ToDo:implement when in need
                                }
                            }

                            ++i;

                        } while (iternext(iter));

                        outputs.push_back(matN);
                    }
                    else
                    {
                        success = false; 
                        outputs.push_back("");
                    }
                }
                else
                {
                    success = false; 
                    outputs.push_back("");
                }

                NpyIter_Deallocate(iter);
            }
            else
            {
                success = true;
                matN->Dimension(dims,hwMatrixN::REAL);
                outputs.push_back(matN);
            }
        }
        else
        {
            success = false; 
            outputs.push_back("");
        }
    }

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
        HandleException();
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

    HandleException();
    
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