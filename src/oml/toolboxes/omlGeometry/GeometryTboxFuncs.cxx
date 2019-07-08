/**
* @file GeometryTboxFuncs.cxx
* @date November, 2018
* Copyright (C) 2018-2019 Altair Engineering, Inc.  
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
#include <memory>       // for std::unique_ptr
#include <assert.h>

#ifdef OS_WIN
#include <windows.h>    // for makeTempFile2()
#endif

#include "BuiltInFuncs.h"
#include "BuiltInFuncsUtils.h"
#include "SignalHandlerBase.h"
#include "OML_Error.h"
#include "GeometryTboxFuncs.h"
#include "GeometryFuncs.h"

#define GEOM "Geometry"
#define TBOXVERSION 2019.1

//------------------------------------------------------------------------------
// Entry point which registers oml Geometry functions with oml
//------------------------------------------------------------------------------
int InitDll(EvaluatorInterface eval)
{
    eval.RegisterBuiltInFunction("convhull", &OmlConvHull,
                                 FunctionMetaData(-3, -2, GEOM));
    eval.RegisterBuiltInFunction("convhulln", &OmlConvHulln,
                                 FunctionMetaData(-2, -2, GEOM));
    eval.RegisterBuiltInFunction("delaunay", &OmlDelaunay,
                                 FunctionMetaData(-3, 1, GEOM));
    eval.RegisterBuiltInFunction("delaunayn", &OmlDelaunayn,
                                 FunctionMetaData(-2, 1, GEOM));
    return 1;
}
//------------------------------------------------------------------------------
// Creates temporary file and returns file pointer
//------------------------------------------------------------------------------
static FILE* makeTempFile2()
{
    // copied from BuiltInFuncs.cpp; discard if modifications not retained
#ifdef OS_WIN
    DWORD len = GetTempPath(0, nullptr);

    if (len != 0)
    {
        LPSTR dirpath = new CHAR[len + 1];
        len = GetTempPath(len + 1, dirpath);
        if (len != 0)
        {
            LPSTR path = new CHAR[MAX_PATH + 1];
            if (GetTempFileName(dirpath, "HW_TEMP_FILE_BUFFER", 0, path))
            {
                FILE* f = fopen(path, "w");
                // FILE* f = fopen(path, "w", stderr); // redirect stderr
                delete[] path;
                if (f)
                {
                    delete[] dirpath;
                    return f;
                }
            }
        }
        delete[] dirpath;
    }
#else
    FILE* f = tmpfile();
    if (f)
        return f;
#endif
    throw OML_Error(HW_ERROR_PROBCREATTEMPF);
}
//------------------------------------------------------------------------------
// Computes the 2D convex hull [convhull command]
//------------------------------------------------------------------------------
bool OmlConvHull(EvaluatorInterface           eval,
                 const std::vector<Currency>& inputs,
                 std::vector<Currency>&       outputs)
{
    // unpack inputs
    int nargin = eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin < 2 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix())
    {
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_TYPE);
    }

    if (!inputs[1].IsMatrix())
    {
        throw OML_Error(OML_ERR_VECTOR, 2, OML_VAR_TYPE);
    }

    const hwMatrix* x = inputs[0].Matrix();
    const hwMatrix* y = inputs[1].Matrix();
    std::string options = "Qt";

    if (nargin > 2)
    {
        if (inputs[2].IsString())
        {
            options = inputs[2].StringVal();
        }
        else if (inputs[2].IsCellArray())
        {
            HML_CELLARRAY* cell = inputs[2].CellArray();

            for (int i = 0; i < cell->Size(); ++i)
            {
                const Currency& idx = (*cell)(i);

                if (idx.IsString())
                {
                    if (i == 0)
                    {
                        options = idx.StringVal();
                    }
                    else
                    {
                        options += " " + idx.StringVal();
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_STRING_STRINGCELL, 3, OML_VAR_TYPE);
                }
            }
        }
        else if (inputs[2].IsMatrix())
        {
            const hwMatrix* mat = inputs[2].Matrix();

            if (mat->M() != 0 || mat->N() != 0)
            {
                throw OML_Error(OML_ERR_STRING_STRINGCELL, 3, OML_VAR_TYPE);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_STRING_STRINGCELL, 3, OML_VAR_TYPE);
        }
    }

    // check for batch mode and set error file
    bool batchmode = false;
    FILE* errfile = nullptr;

    SignalHandlerBase* handler = eval.GetSignalHandler();
    if (handler && (handler->IsInGuiMode() || handler->IsInConsoleMode()))
    {
        errfile = stderr;
    }
    else // batch mode
    {
        batchmode = true;
        errfile = makeTempFile2();
    }

    // call QHull
    hwMatrixI hull;
    double    area;

    if (nargout == 2)
        area = 1.0;
    else
        area = -1.0;

    hwMathStatus status = ConvexHull(*x, *y, options, hull, area, errfile);

    if (batchmode)
    {
        if (status == HW_MATH_ERR_QHULL)
        {
            // do something with the contents
        }

        int success = fclose(errfile);
    }

    BuiltInFuncsUtils::CheckMathStatus(eval, status);

    hwMatrix* hullidx = EvaluatorInterface::allocateMatrix(hull.Size(), 1, hwMatrix::REAL);

    for (int i = 0; i < hullidx->M(); ++i)
    {
        (*hullidx)(i) = static_cast<int> (hull(i));
    }

    outputs.push_back(hullidx);

    if (nargout == 2)
        outputs.push_back(area);

    return true;
}
//------------------------------------------------------------------------------
// Computes the ND convex hull [convhull command]
//------------------------------------------------------------------------------
bool OmlConvHulln(EvaluatorInterface           eval,
                  const std::vector<Currency>& inputs,
                  std::vector<Currency>&       outputs)
{
    // unpack inputs
    int nargin = eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix())
    {
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_TYPE);
    }

    const hwMatrix* P = inputs[0].Matrix();
    std::string options;

    if (P->N() < 5)
        options = "Qt";
    else
        options = "Qt Qx";

    if (nargin > 1)
    {
        if (inputs[1].IsString())
        {
            options = inputs[1].StringVal();
        }
        else if (inputs[1].IsCellArray())
        {
            HML_CELLARRAY* cell = inputs[1].CellArray();

            for (int i = 0; i < cell->Size(); ++i)
            {
                const Currency& idx = (*cell)(i);

                if (idx.IsString())
                {
                    if (i == 0)
                    {
                        options = idx.StringVal();
                    }
                    else
                    {
                        options += " " + idx.StringVal();
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_STRING_STRINGCELL, 2, OML_VAR_TYPE);
                }
            }
        }
        else if (inputs[1].IsMatrix())
        {
            const hwMatrix* mat = inputs[1].Matrix();

            if (mat->M() != 0 || mat->N() != 0)
            {
                throw OML_Error(OML_ERR_STRING_STRINGCELL, 2, OML_VAR_TYPE);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_STRING_STRINGCELL, 2, OML_VAR_TYPE);
        }
    }

    // check for batch mode and set error file
    bool batchmode;
    FILE* errfile = nullptr;

    SignalHandlerBase* handler = eval.GetSignalHandler();
    if (handler && (handler->IsInGuiMode() || handler->IsInConsoleMode()))
    {
        errfile = stderr;
    }
    else // batch mode
    {
        batchmode = true;
        errfile = makeTempFile2();
    }

    // call QHull
    double volume;

    if (nargout == 2)
        volume = 1.0;
    else
        volume = -1.0;

    hwMatrixI    hull;
    hwMatrix     Ptrans;
    hwMathStatus status = Ptrans.Transpose(*P);

    status = ConvexHulln(Ptrans, options, hull, volume, errfile);

    if (batchmode)
    {
        if (status == HW_MATH_ERR_QHULL)
        {
            // do something with the contents
        }

        int success = fclose(errfile);
    }

    BuiltInFuncsUtils::CheckMathStatus(eval, status);

    hwMatrix* hullidx = EvaluatorInterface::allocateMatrix(hull.M(), hull.N(), hwMatrix::REAL);

    for (int i = 0; i < hullidx->Size(); ++i)
    {
        (*hullidx)(i) = static_cast<int> (hull(i));
    }

    outputs.push_back(hullidx);

    if (nargout == 2)
        outputs.push_back(volume);

    return true;
}
//------------------------------------------------------------------------------
// Computes the 2D or 3D Delaunay triangulation [delaunay command]
//------------------------------------------------------------------------------
bool OmlDelaunay(EvaluatorInterface           eval,
	             const std::vector<Currency>& inputs,
	             std::vector<Currency>&       outputs)
{
    // unpack inputs
    int nargin = eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin == 1)
    {
        return OmlDelaunayn(eval, inputs, outputs);
    }

    if (nargin == 2 && (inputs[1].IsString() || inputs[1].IsCellArray()))
    {
        return OmlDelaunayn(eval, inputs, outputs);
    }

    if (nargin < 2 || nargin > 4)
        throw OML_Error(OML_ERR_NUMARGIN);

    int numDim;

    if (nargin > 2 && (inputs[nargin-1].IsString() || inputs[nargin-1].IsCellArray()))
    {
        numDim = nargin - 1;
    }
    else
    {
        numDim = nargin;
    }

    if (numDim > 3)
        throw OML_Error(hwMathStatus(HW_MATH_ERR_QHULL_DIMS23));

    if (!inputs[0].IsMatrix())
    {
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_TYPE);
    }

    const hwMatrix* v = inputs[0].Matrix();

    if (!v->IsReal())
    {
        throw OML_Error(OML_ERR_REAL, 1, OML_VAR_TYPE);
    }

    if (!v->IsVector())
    {
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_TYPE);
    }

    int numPts = v->Size();
    hwMatrix* P = EvaluatorInterface::allocateMatrix(numPts, numDim, hwMatrix::REAL);

    for (int j = 0; j < numPts; ++j)
    {
        (*P)(j, 0) = (*v)(j);
    }

    for (int i = 1; i < numDim; ++i)
    {
        if (!inputs[i].IsMatrix())
        {
            throw OML_Error(OML_ERR_VECTOR, i+1, OML_VAR_TYPE);
        }

        const hwMatrix* v = inputs[i].Matrix();

        if (!v->IsReal())
        {
            throw OML_Error(OML_ERR_REAL, i + 1, OML_VAR_TYPE);
        }

        if (!v->IsVector())
        {
            throw OML_Error(OML_ERR_VECTOR, i + 1, OML_VAR_TYPE);
        }

        if (v->Size() != numPts)
        {
            throw OML_Error(OML_ERR_ARRAYSIZE, 1, i + 1);
        }

        for (int j = 0; j < numPts; ++j)
        {
            (*P)(j, i) = (*v)(j);
        }
    }

    std::vector<Currency> inputs2;
    inputs2.push_back(P);

    if (numDim == nargin - 1)
        inputs2.push_back(inputs[nargin - 1]);

    return OmlDelaunayn(eval, inputs2, outputs);
}
//------------------------------------------------------------------------------
// Computes the ND Delaunayn triangulation [delaunayn command]
//------------------------------------------------------------------------------
bool OmlDelaunayn(EvaluatorInterface           eval,
                  const std::vector<Currency>& inputs,
                  std::vector<Currency>&       outputs)
{
    // unpack inputs
    int nargin = static_cast<int> (inputs.size()); // could be different from eval.GetNarginValue();
    int nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsMatrix() && !inputs[0].IsScalar())
    {
        throw OML_Error(OML_ERR_VECTOR, 1, OML_VAR_TYPE);
    }

    const hwMatrix* P = inputs[0].Matrix();
    std::string options;

    if (P->N() < 4)
        options = "Qt Qbb Qc Qz";
    else
        options = "Qt Qbb Qc Qx";

    if (nargin > 1)
    {
        if (inputs[1].IsString())
        {
            options = inputs[1].StringVal();
        }
        else if (inputs[1].IsCellArray())
        {
            HML_CELLARRAY* cell = inputs[1].CellArray();

            for (int i = 0; i < cell->Size(); ++i)
            {
                const Currency& idx = (*cell)(i);

                if (idx.IsString())
                {
                    if (i == 0)
                    {
                        options = idx.StringVal();
                    }
                    else
                    {
                        options += " " + idx.StringVal();
                    }
                }
                else
                {
                    throw OML_Error(OML_ERR_STRING_STRINGCELL, 2, OML_VAR_TYPE);
                }
            }
        }
        else if (inputs[1].IsMatrix())
        {
            const hwMatrix* mat = inputs[1].Matrix();

            if (mat->M() != 0 || mat->N() != 0)
            {
                throw OML_Error(OML_ERR_STRING_STRINGCELL, 2, OML_VAR_TYPE);
            }
        }
        else
        {
            throw OML_Error(OML_ERR_STRING_STRINGCELL, 2, OML_VAR_TYPE);
        }
    }

    // check for batch mode and set error file
    bool batchmode = false;
    FILE* errfile = nullptr;

    SignalHandlerBase* handler = eval.GetSignalHandler();
    if (handler && (handler->IsInGuiMode() || handler->IsInConsoleMode()))
    {
        errfile = stderr;
    }
    else // batch mode
    {
        batchmode = true;
        errfile = makeTempFile2();
    }
    // call QHull
    hwMatrixI    triang;
    hwMatrix     Ptrans;
    hwMathStatus status = Ptrans.Transpose(*P);

    status = Delaunayn(Ptrans, options, triang, errfile);

    if (batchmode)
    {
        if (status == HW_MATH_ERR_QHULL)
        {
            // do something with the contents
        }

        int success = fclose(errfile);
    }

    BuiltInFuncsUtils::CheckMathStatus(eval, status);

    hwMatrix* triangidx = EvaluatorInterface::allocateMatrix(triang.M(), triang.N(), hwMatrix::REAL);

    for (int i = 0; i < triangidx->Size(); ++i)
    {
        (*triangidx)(i) = static_cast<int> (triang(i));
    }

    outputs.push_back(triangidx);

    return true;
}
//------------------------------------------------------------------------------
// Returns toolbox version
//------------------------------------------------------------------------------
double GetToolboxVersion(EvaluatorInterface eval)
{
	return TBOXVERSION;
}
