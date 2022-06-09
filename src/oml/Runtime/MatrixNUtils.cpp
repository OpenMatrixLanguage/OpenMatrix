/**
* @file MatrixNUtils.cpp
* @date January 2018
* Copyright (C) 2018 Altair Engineering, Inc.  
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
#include "hwMatrix.h"
#include "hwMatrixN.h"
#include "OML_Error.h"
#include "Evaluator.h"
#include "MatrixNUtils.h"
#include <memory>

// Apply a function to each element of an ND matrix. The syntax is oml_func(ND, ...).
bool oml_MatrixNUtil1(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                      const OML_func1 oml_func)
{
    if (!inputs[0].IsNDMatrix())
    {
        // throw OML_Error(OML_ERR_NDMATRIX, 1, OML_STR_VARIABLE);
    }

    const hwMatrixN* matrix = inputs[0].MatrixN();
    const std::vector<int>& dims = matrix->Dimensions();

    // convert from ND to vector
    hwMatrix* slice1D = new hwMatrix;
    matrix->ConvertNDto1D(*slice1D);

    // call computation function
    std::vector<Currency> inputs2;
    std::vector<Currency> outputs2;

    inputs2.push_back(slice1D);

    for (int i = 1; i < inputs.size(); ++i)
    {
        inputs2.push_back(inputs[i]);
    }

    oml_func(eval, inputs2, outputs2);

    // check type of outputs2; not used for assert(ND)
    if (outputs2.size() && outputs2[0].IsMatrix())
    {
        // convert from vector to ND
        hwMatrixN* outMatrix = new hwMatrixN;
        outMatrix->Convert1DtoND(*outputs2[0].GetWritableMatrix(), dims);

   		Currency out(outMatrix);
		out.SetMask(outputs2[0].GetMask());
        outputs.push_back(out);
    }

    if (outputs2.size() == 2 && outputs2[1].IsMatrix())     // [F, E] = log2(ND)
    {
        // convert from vector to ND
        hwMatrixN* outMatrix = new hwMatrixN;
        outMatrix->Convert1DtoND(*outputs2[1].GetWritableMatrix(), dims);

   		Currency out(outMatrix);
		out.SetMask(outputs2[1].GetMask());
        outputs.push_back(out);
    }

    return true;
}

// Apply a function to each element pair from the arguments, with at least one ND matrix. The syntax is oml_func(ND1, ND2).
bool oml_MatrixNUtil2(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                      const OML_func1 oml_func)
{
    size_t nargin = inputs.size();

    if (nargin != 1 && nargin != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    const Currency& op1 = inputs[0];
    const Currency& op2 = inputs[1];
    std::vector<hwSliceArg> sliceArgs;
    sliceArgs.push_back(hwSliceArg());
    hwMatrixN* outMatrix = new hwMatrixN;
    std::vector<Currency> inputs2;
    std::vector<Currency> outputs2;

    if (op1.IsNDMatrix())
    {
        const hwMatrixN* matrix_1 = op1.MatrixN();
        const std::vector<int>& dims_1 = matrix_1->Dimensions();

        // convert from ND to vector
        hwMatrix* slice1D_1 = new hwMatrix;
        matrix_1->ConvertNDto1D(*slice1D_1);
        Currency op1_new(slice1D_1);

        if (op2.IsNDMatrix())
        {
            const hwMatrixN* matrix_2 = op2.MatrixN();
            const std::vector<int>& dims_2 = matrix_2->Dimensions();

            // dims check here
            if (dims_1 != dims_2)
            {
                throw OML_Error(OML_ERR_ARRAYSIZE);
            }

            // convert from ND to vector
            hwMatrix* slice1D_2 = new hwMatrix;
            matrix_2->ConvertNDto1D(*slice1D_2);
            Currency op2_new(slice1D_2);

            // call computation function
            inputs2.push_back(op1_new);
            inputs2.push_back(op2_new);
            oml_func(eval, inputs2, outputs2);

            // convert from vector to ND
            outMatrix->Convert1DtoND(*outputs2[0].GetWritableMatrix(), dims_1);
        }
        else if (!op2.IsMatrix())
        {
            // call computation function
            inputs2.push_back(op1_new);
            inputs2.push_back(op2);
            oml_func(eval, inputs2, outputs2);

            // convert from vector to ND
            outMatrix->Convert1DtoND(*outputs2[0].GetWritableMatrix(), dims_1);
        }
        else
        {
            throw OML_Error(HW_ERROR_UNSUPOP);
        }
    }
    else if (!op1.IsMatrix() && op2.IsNDMatrix())
    {
        const hwMatrixN* matrix_2 = op2.MatrixN();
        const std::vector<int>& dims_2 = matrix_2->Dimensions();

        // convert from ND to vector
        hwMatrix* slice1D_2 = new hwMatrix;
        matrix_2->ConvertNDto1D(*slice1D_2);
        Currency op2_new(slice1D_2);

        // call computation function
        inputs2.push_back(op1);
        inputs2.push_back(op2_new);
        oml_func(eval, inputs2, outputs2);

        // convert from vector to ND
        outMatrix->Convert1DtoND(*outputs2[0].GetWritableMatrix(), dims_2);
    }
    else
    {
        throw OML_Error(HW_ERROR_UNSUPOP);
    }

    outputs.push_back(outMatrix);

    return true;
}

// Apply a function to each vector of an ND matrix in the specified dimension. The
// syntax is oml_func(ND, dim). This utility applies to functions for which the
// output for each vector input is a scalar.
bool oml_MatrixNUtil3(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                      const OML_func1 oml_func, int dimArg)
{
    size_t nargin = inputs.size();
    size_t nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsNDMatrix())
    {
        // throw OML_Error(OML_ERR_NDMATRIX, 1, OML_STR_VARIABLE);
    }

    const hwMatrixN* matrix = inputs[0].MatrixN();
    const std::vector<int>& dims = matrix->Dimensions();
    int numDim = static_cast<int>(dims.size());
    int dim = 0;

    bool passArg2 = false;

    if (dimArg == 23)
    {
        // this is a kludge to tell the function to pass the 2nd argument
        passArg2 = true;
        dimArg -= 20;
    }

    if (!dimArg)
    {
        for (int i = 0; i < numDim; ++i)
        {
            if (dims[i] != 1)
            {
                dim = i;
                break;
            }
        }
    }
    else if (dimArg > 0)
    {
        if (inputs[dimArg-1].IsPositiveInteger())
        {
            dim = static_cast<int>(inputs[dimArg-1].Scalar()) - 1;
        }
        else if (inputs[dimArg-1].IsMatrix())
        {
            if (inputs[dimArg-1].Matrix()->M() || inputs[dimArg-1].Matrix()->N())
            {
                throw OML_Error(OML_ERR_POSINTEGER, dimArg, OML_VAR_VARIABLE);
            }

            for (int i = 0; i < numDim; ++i)
            {
                if (dims[i] != 1)
                {
                    dim = i;
                    break;
                }
            }
        }
        else
        {
            throw OML_Error(OML_ERR_POSINTEGER, dimArg, OML_VAR_VARIABLE);
        }
    }

    // dimension output
    std::vector<int> sliceDims(dims);

    if (dims[dim] != 0)
        sliceDims[dim] = 1;

    hwMatrixN* outMatrix1 = new hwMatrixN(sliceDims, matrix->Type());
    hwMatrixN* outMatrix2 = nullptr;

    if (nargout == 2)
    {
        outMatrix2 = new hwMatrixN(sliceDims, hwMatrixN::REAL);
    }

    // initialize slices
    std::vector<hwSliceArg> sliceArgs;

    for (int i = 0; i < numDim; ++i)
    {
        if (i == dim)
            sliceArgs.push_back(hwSliceArg());
        else
            sliceArgs.push_back(0);
    }

    // apply function to each vector in specified dimension
    std::vector<int> sliceDims2D(2);

    sliceDims2D[0] = -1;
    sliceDims2D[1] = 1;
    int numVecs;
    
    if (dim < dims.size())
    {
        if (dims[dim] == 0)
            numVecs = 0;
        else
            numVecs = matrix->Size() / dims[dim];
    }
    else
    {
        numVecs = matrix->Size();
    }

    for (int i = 0; i < numVecs; ++i)
    {
        // slice matrix to retrieve vector
        hwMatrixN slice;
        hwMatrix* slice2D = new hwMatrix;
        matrix->SliceRHS(sliceArgs, slice);
        slice.Reshape(sliceDims2D);     // reshape to a column
        slice.ConvertNDto2D(*slice2D, false);

        // call computation function
        std::vector<Currency> inputs2;
        std::vector<Currency> outputs2;

        inputs2.push_back(slice2D);

        if (nargin == 3 && inputs[2].IsString())     // geometric mean
        {
            inputs2.push_back(1);
            inputs2.push_back(inputs[2]);
        }
        else if (nargin > 1 && passArg2)
        {
            inputs2.push_back(inputs[1]);
        }

        oml_func(eval, inputs2, outputs2);

        if (outMatrix1->IsReal())
        {
            (*outMatrix1)(i) = outputs2[0].Scalar();
        }
        else
        {
            if (outputs2[0].IsScalar())
                outMatrix1->z(i) = outputs2[0].Scalar();
            else
                outMatrix1->z(i) = outputs2[0].Complex();
        }

        if (nargout == 2)
        {
            (*outMatrix2)(i) = outputs2[1].Scalar();
        }

        // advance slice indices
        for (int j = 0; j < numDim; ++j)
        {
            if (j == dim)
            {
                continue;
            }

            // increment index j if possible
            if (sliceArgs[j].Scalar() < (int) dims[j]-1)
            {
                ++sliceArgs[j].Scalar();
                break;
            }

            // index j is maxed out, so reset and continue to j+1
            sliceArgs[j].Scalar() = 0;
        }
    }

    outputs.push_back(outMatrix1);

    if (nargout == 2)
    {
        outputs.push_back(outMatrix2);
    }

    return true;
}

// Apply a function to each vector of an ND matrix in the specified dimension. The
// syntax is oml_func(ND, ..., dim). This utility applies to functions for which the
// output for each vector input is a vector.
bool oml_MatrixNUtil4(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                      const OML_func1 oml_func, int dimArg, int ndArg)
{
    size_t nargin = inputs.size();
    size_t nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 5)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[ndArg-1].IsNDMatrix())
    {
        // throw OML_Error(OML_ERR_NDMATRIX, 1, OML_STR_VARIABLE);
    }

    const hwMatrixN* matrix = inputs[ndArg-1].MatrixN();
    const std::vector<int>& dims = matrix->Dimensions();
    int numDim = static_cast<int>(dims.size());
    int dim = 0;

    if (!dimArg)
    {
        // use first non-singleton dimension
        for (int i = 0; i < numDim; ++i)
        {
            if (dims[i] != 1)
            {
                dim = i;
                break;
            }
        }
    }
    else if (dimArg > 0)
    {
        if (!inputs[dimArg-1].IsPositiveInteger())     // sort(ND, dim, mode), cumsum(ND, dim)
            throw OML_Error(OML_ERR_POSINTEGER, 1, OML_VAR_TYPE);
        
        dim = static_cast<int>(inputs[dimArg-1].Scalar()) - 1;
    }
    else
    {
        dim = -dimArg - 1;  // dimArg is the dimension when < 0
        dimArg = 0;
    }

    // initialize slices
    std::vector<hwSliceArg> sliceArgs;

    for (int i = 0; i < numDim; ++i)
    {
        if (i == dim)
            sliceArgs.push_back(hwSliceArg());
        else
            sliceArgs.push_back(0);
    }

    // apply function to each vector in specified dimension
    hwMatrixN* outMatrix1 = nullptr;
    hwMatrixN* outMatrix2 = nullptr;
    std::vector<int> sliceDims2D(2);

    sliceDims2D[0] = -1;
    sliceDims2D[1] = 1;
    int numVecs;

    if (dim < dims.size())
    {
        if (dims[dim] == 0)
            numVecs = 0;
        else
            numVecs = matrix->Size() / dims[dim];
    }
    else
    {
        numVecs = matrix->Size();
    }

    for (int i = 0; i < numVecs; ++i)
    {
        // slice matrix to retrieve vector
        hwMatrixN slice;
        hwMatrix* slice2D = new hwMatrix;
        matrix->SliceRHS(sliceArgs, slice);
        slice.Reshape(sliceDims2D);     // reshape to a column
        slice.ConvertNDto2D(*slice2D, false);

        // call computation function
        std::vector<Currency> inputs2;
        std::vector<Currency> outputs2;

        switch (dimArg)
        {
        case 0:
            for (int j = 0; j < ndArg - 1; ++j)
                inputs2.push_back(inputs[j]);           // interp1(x, ND, ...), filtfilt(b, a, ND)

            inputs2.push_back(slice2D);

            for (int j = ndArg; j < nargin; ++j)
                inputs2.push_back(inputs[j]);           // sort(ND, mode), circshift(ND, n)
            break;
        case 2:
            inputs2.push_back(slice2D);

            if (nargin == 3 && inputs[2].IsString())    // sort(ND, dim, mode)
                inputs2.push_back(inputs[2]);
            break;
        case 3:
            inputs2.push_back(slice2D);

            if (nargin == 3 && inputs[1].IsInteger())   // circshift(ND, n, dim)
                inputs2.push_back(inputs[1]);
            break;
        case 4:
            inputs2.push_back(slice2D);                 // resample(ND, p, q)
            inputs2.push_back(inputs[1]);               // implicitly, dim = 0
            inputs2.push_back(inputs[2]);
            break;
        case 5:
            for (int j = 0; j < ndArg - 1; ++j)
                inputs2.push_back(inputs[j]);           // filter(b, a, ND, [], dim)

            inputs2.push_back(slice2D);

            for (int j = ndArg; j < dimArg - 2; ++j)
                inputs2.push_back(inputs[j]);
            break;
        }

        oml_func(eval, inputs2, outputs2);

        // convert from vector to ND
        slice = hwMatrixN();
        slice.Convert2DtoND(*outputs2[0].GetWritableMatrix(), false);

        if (i == 0)
        {
            // dimension output
            std::vector<int> outdims(dims);

            if (dim < outdims.size())
                outdims[dim] = slice.Size();

            outMatrix1 = new hwMatrixN(outdims, slice.Type());

            if (nargout == 2)
                outMatrix2 = new hwMatrixN(outdims, hwMatrixN::REAL);
        }

        outMatrix1->SliceLHS(sliceArgs, slice);

        if (nargout == 2)
        {
            slice = hwMatrixN();
            slice.Convert2DtoND(*outputs2[1].GetWritableMatrix(), false);
            outMatrix2->SliceLHS(sliceArgs, slice);
        }

        // advance slice indices
        for (int j = 0; j < numDim; ++j)
        {
            if (j == dim)
            {
                continue;
            }

            // increment index j if possible
            if (sliceArgs[j].Scalar() < (int) dims[j]-1)
            {
                ++sliceArgs[j].Scalar();
                break;
            }

            // index j is maxed out, so reset and continue to j+1
            sliceArgs[j].Scalar() = 0;
        }
    }

    if (!numVecs)
    {
        // dimension output for the empty matrix case
        // this needs to be enhanced if ever dims[dim] != slice.Size()
        outMatrix1 = new hwMatrixN(dims, hwMatrixN::REAL);

        if (nargout == 2)
            outMatrix2 = new hwMatrixN(dims, hwMatrixN::REAL);
    }

    outputs.push_back(outMatrix1);

    if (nargout == 2)
        outputs.push_back(outMatrix2);

    return true;
}

// Apply a function to an argument list that contains vectors stored as ND matrices. The
// syntax is oml_func(..., ND, ...). Output vectors will be in the same dimension as the first
// argument if it is a vector.
// Examples: not yet in use
bool oml_MatrixNVecs(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                     const OML_func1 oml_func)
{
    size_t nargin = inputs.size();
    std::vector<Currency> inputs2;
    std::vector<Currency> outputs2;
    size_t dim = -1;

    for (size_t i = 0; i < nargin; ++i)
    {
        if (inputs[i].IsNDMatrix())
        {
            const hwMatrixN* matrixN = inputs[i].MatrixN();

            if (matrixN->IsVector())
            {
                if (i == 0)
                {
                    const std::vector<int>& dims = matrixN->Dimensions();

                    for (size_t j = dims.size()-1; j > 0; --j)
                    {
                        if (dims[j] != 1)
                        {
                            dim = j;
                            break;
                        }
                    }
                }

                std::vector<hwSliceArg> sliceArgs;
                sliceArgs.push_back(hwSliceArg());
                hwMatrixN slice;
                matrixN->SliceRHS(sliceArgs, slice);
                hwMatrix* vec = new hwMatrix;
                slice.ConvertNDto2D(*vec, false);
                inputs2.push_back(vec);
            }
            else
            {
                inputs2.push_back(inputs[i]);
            }
        }
        else
        {
            inputs2.push_back(inputs[i]);
        }
    }

    oml_func(eval, inputs2, outputs2);

    size_t nargout = outputs2.size();

    for (size_t i = 0; i < nargin; ++i)
    {
        if (outputs2[i].IsMatrix())
        {
            hwMatrix* matrix = outputs2[i].GetWritableMatrix();

            if (matrix->IsVector())
            {
                std::vector<int> dims(dim+1);

                for (int j = 0; j < dim; ++j)
                    dims[j] = 1;

                dims[dim] = matrix->Size();

                hwMatrixN* vec = new hwMatrixN;
                vec->Convert2DtoND(*matrix, false);
                vec->Reshape(dims);
                outputs.push_back(vec);
            }
            else
            {
                outputs.push_back(outputs2[i]);
            }
        }
        else
        {
            outputs.push_back(outputs2[i]);
        }
    }

    return true;
}

// ND support for dot and cross functions
bool oml_MatrixN_VecProd(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                         const OML_func1 oml_func, int vecLength)
{
    size_t nargin = inputs.size();
    size_t nargout = eval.GetNargoutValue();

    if (nargin != 2 && nargin != 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsNDMatrix() || !inputs[1].IsNDMatrix())
    {
        // throw OML_Error(OML_ERR_NDMATRIX, 1, OML_STR_VARIABLE);
    }

    const hwMatrixN* matrix1 = inputs[0].MatrixN();
    const hwMatrixN* matrix2 = inputs[1].MatrixN();
    const std::vector<int>& dims1 = matrix1->Dimensions();
    const std::vector<int>& dims2 = matrix2->Dimensions();
    int numDim = static_cast<int>(dims1.size());
    int dim = 0;

    if (dims1 != dims2)
    {
        throw OML_Error(OML_ERR_ARRAYSIZE, 1, 2);
    }

    if (nargin == 2)
    {
        for (int i = 0; i < numDim; ++i)
        {
            if ((vecLength == 1 && dims1[i] != 1) ||    // dot product
                (vecLength == 3 && dims1[i] == 3))      // cross product
            {
                dim = i;
                break;
            }
        }
    }
    else if (inputs[2].IsPositiveInteger())
    {
        dim = static_cast<int>(inputs[2].Scalar()) - 1;
    }
    else
    {
        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_DIM);
    }

    // dimension output
    std::vector<int> sliceDims(dims1);

    if (dim < dims1.size())
        sliceDims[dim] = vecLength;

    hwMatrixN::DataType type;

    if (matrix1->Type() == hwMatrixN::REAL && matrix2->Type() == hwMatrixN::REAL)
        type = hwMatrixN::REAL;
    else
        type = hwMatrixN::COMPLEX;

    hwMatrixN* outMatrix = new hwMatrixN(sliceDims, type);

    // initialize slices
    std::vector<hwSliceArg> sliceArgs;

    for (int i = 0; i < numDim; ++i)
    {
        if (i == dim)
            sliceArgs.push_back(hwSliceArg());
        else
            sliceArgs.push_back(0);
    }

    // apply function to each vector in specified dimension
    std::vector<int> sliceDims2D(2);

    sliceDims2D[0] = -1;
    sliceDims2D[1] = 1;
    int numVecs;
    
    if (dim < dims1.size())
        numVecs = matrix1->Size() / dims1[dim];
    else
        numVecs = matrix1->Size();

    for (int i = 0; i < numVecs; ++i)
    {
        // slice matrix to retrieve vector
        hwMatrixN slice1;
        hwMatrixN slice2;
        hwMatrix* slice2D_1 = new hwMatrix;
        matrix1->SliceRHS(sliceArgs, slice1);
        slice1.Reshape(sliceDims2D);     // reshape to a column
        slice1.ConvertNDto2D(*slice2D_1, false);

        hwMatrix* slice2D_2 = new hwMatrix;
        matrix2->SliceRHS(sliceArgs, slice2);
        slice2.Reshape(sliceDims2D);     // reshape to a column
        slice2.ConvertNDto2D(*slice2D_2, false);

        // call computation function
        std::vector<Currency> inputs2;
        std::vector<Currency> outputs2;

        inputs2.push_back(slice2D_1);
        inputs2.push_back(slice2D_2);
        oml_func(eval, inputs2, outputs2);

        if (outMatrix->IsReal())
        {
            if (vecLength == 1)     // dot product
            {
                (*outMatrix)(i) = outputs2[0].Scalar();
            }
            else                    // cross product
            {
                // convert from vector to ND
                slice1.Convert2DtoND(*outputs2[0].GetWritableMatrix(), false);
                outMatrix->SliceLHS(sliceArgs, slice1);
            }
        }
        else     // dot product
        {
            if (outputs2[0].IsScalar())
                outMatrix->z(i) = outputs2[0].Scalar();
            else
                outMatrix->z(i) = outputs2[0].Complex();
        }

        // advance slice indices
        for (int j = 0; j < numDim; ++j)
        {
            if (j == dim)
            {
                continue;
            }

            // increment index j if possible
            if (sliceArgs[j].Scalar() < (int) dims1[j]-1)
            {
                ++sliceArgs[j].Scalar();
                break;
            }

            // index j is maxed out, so reset and continue to j+1
            sliceArgs[j].Scalar() = 0;
        }
    }

    outputs.push_back(outMatrix);

    return true;
}

// ND support for sum function
bool oml_MatrixN_sum(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                     int dimArg)
{
    size_t nargin = inputs.size();
    size_t nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsNDMatrix())
    {
        // throw OML_Error(OML_ERR_NDMATRIX, 1, OML_STR_VARIABLE);
    }

    const hwMatrixN* matrix = inputs[0].MatrixN();
    const std::vector<int>& dims = matrix->Dimensions();
    int numDim = static_cast<int>(dims.size());
    int dim = 0;

    if (!dimArg)
    {
        for (int i = 0; i < numDim; ++i)
        {
            if (dims[i] != 1)
            {
                dim = i;
                break;
            }
        }
    }
    else if (dimArg > 0)
    {
        if (inputs[dimArg - 1].IsPositiveInteger())
        {
            dim = static_cast<int>(inputs[dimArg - 1].Scalar()) - 1;
        }
        else if (inputs[dimArg - 1].IsMatrix())
        {
            if (inputs[dimArg - 1].Matrix()->M() || inputs[dimArg - 1].Matrix()->N())
            {
                throw OML_Error(OML_ERR_POSINTEGER, dimArg, OML_VAR_VARIABLE);
            }

            for (int i = 0; i < numDim; ++i)
            {
                if (dims[i] != 1)
                {
                    dim = i;
                    break;
                }
            }
        }
        else
        {
            throw OML_Error(OML_ERR_POSINTEGER, dimArg, OML_VAR_VARIABLE);
        }
    }

    if (dim >= dims.size())
    {
        outputs.push_back(inputs[0]);
        return true;
    }

    // dimension output
    std::vector<int> sliceDims(dims);
    sliceDims[dim] = 1;

    hwMatrixN* sum = new hwMatrixN(sliceDims, matrix->Type());

    if (dims[dim] == 0)
    {
        sum->SetElements(0.0);
        outputs.push_back(sum);
        return true;
    }

    if (dim == 0)
    {
        // sum each vector in dimension 1
        int numVecs = matrix->Size() / dims[0];
        std::vector<int> rhsMatrixIndex(numDim);

        for (int i = 0; i < numVecs; ++i)
        {
            // set the rhsMatrix indices to the first index in each slice
            int start = matrix->Index(rhsMatrixIndex);

            if (matrix->IsReal())
            {
                const double* real = matrix->GetRealData() + start;
                double sumVec = 0.0;

                for (int k = 0; k < dims[dim]; ++k)
                    sumVec += *real++;

                (*sum)(i) = sumVec;
            }
            else
            {
                const hwComplex* cmplx = matrix->GetComplexData() + start;
                hwComplex sumVec(0.0, 0.0);

                for (int k = 0; k < dims[dim]; ++k)
                    sumVec += *cmplx++;

                sum->z(i) = sumVec;
            }

            // advance slice indices
            for (int j = 1; j < numDim; ++j)
            {
                // increment index j if possible
                if (rhsMatrixIndex[j] < static_cast<int> (dims[j]) - 1)
                {
                    ++rhsMatrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                rhsMatrixIndex[j] = 0;
            }
        }
    }
    else
    {
        sum->SetElements(0.0);
        int numVecs = matrix->Size() / dims[0];
        std::vector<int> rhsMatrixIndex(numDim);

        if (dims[0] > numVecs || dims[0] > 16)
        {
            // sum each dimension 1 vector into the specified dimension
            for (int i = 0; i < numVecs; ++i)
            {
                // set the rhsMatrix indices to the first index in each slice
                int start1 = matrix->Index(rhsMatrixIndex);
                int temp = rhsMatrixIndex[dim];
                rhsMatrixIndex[dim] = 0;
                int start2 = sum->Index(rhsMatrixIndex);
                rhsMatrixIndex[dim] = temp;

                if (matrix->IsReal())
                {
                    const double* realD = matrix->GetRealData() + start1;
                    const double* realS = sum->GetRealData() + start2;
                    hwMatrix dataVec(dims[0], 1, (void*)realD, hwMatrix::REAL);
                    hwMatrix sumVec(dims[0], 1, (void*)realS, hwMatrix::REAL);

                    sumVec += dataVec;
                }
                else
                {
                    const hwComplex* cmplxD = matrix->GetComplexData() + start1;
                    const hwComplex* cmplxS = sum->GetComplexData() + start2;
                    hwMatrix dataVec(dims[0], 1, (void*)cmplxD, hwMatrix::COMPLEX);
                    hwMatrix sumVec(dims[0], 1, (void*)cmplxS, hwMatrix::COMPLEX);

                    sumVec += dataVec;
                }

                // advance slice indices
                for (int j = 1; j < numDim; ++j)
                {
                    // increment index j if possible
                    if (rhsMatrixIndex[j] < static_cast<int> (dims[j]) - 1)
                    {
                        ++rhsMatrixIndex[j];
                        break;
                    }

                    // index j is maxed out, so reset and continue to j+1
                    rhsMatrixIndex[j] = 0;
                }
            }
        }
        else
        {
            // sum each vector along the dimension of interest
            int numVecs = matrix->Size() / dims[dim];
            int stride = matrix->Stride(dim);

            for (int i = 0; i < numVecs; ++i)
            {
                // set the rhsMatrix indices to the first index in each slice
                int start = matrix->Index(rhsMatrixIndex);

                if (matrix->IsReal())
                {
                    const double* realD = matrix->GetRealData() + start;
                    hwMatrix dataVec(stride, dims[dim], (void*)realD, hwMatrix::REAL);
                    (*sum)(i) = dataVec.RowSumReal(0);
                }
                else
                {
                    const hwComplex* cmplxD = matrix->GetComplexData() + start;
                    hwMatrix dataVec(stride, dims[dim], (void*)cmplxD, hwMatrix::COMPLEX);
                    sum->z(i) = dataVec.RowSumComplex(0);
                }

                // advance slice indices
                for (int j = 0; j < numDim; ++j)
                {
                    if (j == dim)
                        continue;

                    // increment index j if possible
                    if (rhsMatrixIndex[j] < static_cast<int> (dims[j]) - 1)
                    {
                        ++rhsMatrixIndex[j];
                        break;
                    }

                    // index j is maxed out, so reset and continue to j+1
                    rhsMatrixIndex[j] = 0;
                }
            }
        }
    }

    outputs.push_back(sum);

    return true;
}

// ND support for cumsum function
bool oml_MatrixN_cumsum(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                        int dimArg)
{
    size_t nargin = inputs.size();
    size_t nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsNDMatrix())
    {
        // throw OML_Error(OML_ERR_NDMATRIX, 1, OML_STR_VARIABLE);
    }

    const hwMatrixN* matrix = inputs[0].MatrixN();
    const std::vector<int>& dims = matrix->Dimensions();
    int numDim = static_cast<int>(dims.size());
    int dim = 0;

    if (!dimArg)
    {
        for (int i = 0; i < numDim; ++i)
        {
            if (dims[i] != 1)
            {
                dim = i;
                break;
            }
        }
    }
    else if (dimArg > 0)
    {
        if (inputs[dimArg - 1].IsPositiveInteger())
        {
            dim = static_cast<int>(inputs[dimArg - 1].Scalar()) - 1;
        }
        else if (inputs[dimArg - 1].IsMatrix())
        {
            if (inputs[dimArg - 1].Matrix()->M() || inputs[dimArg - 1].Matrix()->N())
            {
                throw OML_Error(OML_ERR_POSINTEGER, dimArg, OML_VAR_VARIABLE);
            }

            for (int i = 0; i < numDim; ++i)
            {
                if (dims[i] != 1)
                {
                    dim = i;
                    break;
                }
            }
        }
        else
        {
            throw OML_Error(OML_ERR_POSINTEGER, dimArg, OML_VAR_VARIABLE);
        }
    }

    if (dim >= dims.size())
    {
        outputs.push_back(inputs[0]);
        return true;
    }

    // dimension output
    std::vector<int> sliceDims(dims);
    hwMatrixN* cumsum = nullptr;

    if (dim == 0)
    {
        cumsum = new hwMatrixN(sliceDims, matrix->Type());

        if (dims[dim] == 0)
        {
            outputs.push_back(cumsum);
            return true;
        }

        // cumsum each vector in dimension 1
        int numVecs = matrix->Size() / dims[0];
        std::vector<int> rhsMatrixIndex(numDim);

        for (int i = 0; i < numVecs; ++i)
        {
            // set the rhsMatrix indices to the first index in each slice
            int start = matrix->Index(rhsMatrixIndex);

            if (matrix->IsReal())
            {
                const double* real = matrix->GetRealData() + start;
                double* result = cumsum->GetRealData() + start;
                double sum = 0.0;

                for (int k = 0; k < dims[0]; ++k)
                {
                    sum += *(real++);
                    *(result++) = sum;
                }
            }
            else
            {
                const hwComplex* cmplx = matrix->GetComplexData() + start;
                hwComplex* result = cumsum->GetComplexData() + start;
                hwComplex sum;

                for (int k = 0; k < dims[0]; ++k)
                {
                    sum += *(cmplx++);
                    *(result++) = sum;
                }
            }

            // advance slice indices
            for (int j = 1; j < numDim; ++j)
            {
                // increment index j if possible
                if (rhsMatrixIndex[j] < static_cast<int> (dims[j]) - 1)
                {
                    ++rhsMatrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                rhsMatrixIndex[j] = 0;
            }
        }
    }
    else
    {
        cumsum = new hwMatrixN(*matrix);

        if (dims[dim] == 0)
        {
            outputs.push_back(cumsum);
            return true;
        }

        // cumsum each dimension 1 vector into in specified dimension
        int numVecs = matrix->Size() / dims[0];
        std::vector<int> rhsMatrixIndex(numDim);
        const double* realR = cumsum->GetRealData();
        const hwComplex* cmplxR = cumsum->GetComplexData();

        for (int i = 0; i < numVecs; ++i)
        {
            // set the rhsMatrix indices to the first index in each slice
            int start1 = matrix->Index(rhsMatrixIndex);
            int index = rhsMatrixIndex[dim];

            if (index != 0)
            {
                rhsMatrixIndex[dim] = index - 1;
                int start2 = cumsum->Index(rhsMatrixIndex);
                rhsMatrixIndex[dim] = index;

                if (matrix->IsReal())
                {
                    hwMatrix cumsumVec1(dims[0], 1, (void*)(realR + start1), hwMatrix::REAL);
                    hwMatrix cumsumVec2(dims[0], 1, (void*)(realR + start2), hwMatrix::REAL);

                    cumsumVec1.AddEquals(cumsumVec2);
                }
                else
                {
                    hwMatrix cumsumVec1(dims[0], 1, (void*)(cmplxR + start1), hwMatrix::COMPLEX);
                    hwMatrix cumsumVec2(dims[0], 1, (void*)(cmplxR + start2), hwMatrix::COMPLEX);

                    cumsumVec1.AddEquals(cumsumVec2);
                }
            }

            // advance slice indices
            for (int j = 1; j < numDim; ++j)
            {
                // increment index j if possible
                if (rhsMatrixIndex[j] < static_cast<int> (dims[j]) - 1)
                {
                    ++rhsMatrixIndex[j];
                    break;
                }

                // index j is maxed out, so reset and continue to j+1
                rhsMatrixIndex[j] = 0;
            }
        }
    }

    outputs.push_back(cumsum);

    return true;
}

// ND support for diff function
bool oml_MatrixN_diff(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                      const OML_func1 oml_func)
{
    size_t nargin = inputs.size();
    size_t nargout = eval.GetNargoutValue();

    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsNDMatrix())
    {
        // throw OML_Error(OML_ERR_NDMATRIX, 1, OML_STR_VARIABLE);
    }

    const hwMatrixN* matrix = inputs[0].MatrixN();
    const std::vector<int>& dims = matrix->Dimensions();
    int numDim = static_cast<int>(dims.size());
    int dim = 0;

    if (nargin < 3)
    {
        for (int i = 0; i < numDim; ++i)
        {
            if (dims[i] != 1)
            {
                dim = i;
                break;
            }
        }
    }
    else if (nargin == 3 && inputs[2].IsPositiveInteger())
    {
        dim = static_cast<int>(inputs[2].Scalar()) - 1;
    }
    else
    {
        throw OML_Error(OML_ERR_POSINTEGER, 3, OML_VAR_VARIABLE);
    }

    // dimension output
    int cycles = 1;
    int extraCycles = 0;
    std::vector<int> sliceDims(dims);

    if (nargin > 1)
    {
        if (!inputs[1].IsPositiveInteger())
            throw OML_Error(OML_ERR_POSINTEGER, 2, OML_VAR_VARIABLE);

        cycles = static_cast<int>(inputs[1].Scalar());
    }

    if (dim < dims.size())
    {
        if (nargin < 3 && cycles >= dims[dim])
        {
            extraCycles = cycles - (dims[dim]-1);
            cycles -= extraCycles;
            sliceDims[dim] = 1;
        }

        if (dims[dim] > cycles)
        {
            sliceDims[dim] = dims[dim] - cycles;
        }
        else
        {
            sliceDims[dim] = 0;
        }
    }
    else
    {
        for (int i = (int) dims.size(); i < dim; ++i)
            sliceDims.push_back(0);
    }

    hwMatrixN* outMatrix1 = new hwMatrixN(sliceDims, matrix->Type());

    // initialize slices
    std::vector<hwSliceArg> sliceArgs;

    for (int i = 0; i < numDim; ++i)
    {
        if (i == dim)
            sliceArgs.push_back(hwSliceArg());
        else
            sliceArgs.push_back(0);
    }

    // apply function to each vector in specified dimension
    std::vector<int> sliceDims2D(2);

    sliceDims2D[0] = -1;
    sliceDims2D[1] = 1;
    int numVecs;
    
    if (dim < dims.size())
        numVecs = matrix->Size() / dims[dim];
    else
        numVecs = 0;    // actually matrix->Size(), but each yields []

    for (int i = 0; i < numVecs; ++i)
    {
        // slice matrix to retrieve vector
        hwMatrixN slice;
        hwMatrix* slice2D = new hwMatrix;
        matrix->SliceRHS(sliceArgs, slice);
        slice.Reshape(sliceDims2D);     // reshape to a column
        slice.ConvertNDto2D(*slice2D, false);

        // call computation function
        std::vector<Currency> inputs2;
        std::vector<Currency> outputs2;

        inputs2.push_back(slice2D);

        if (nargin > 1)
            inputs2.push_back(cycles);

        oml_func(eval, inputs2, outputs2);

        // convert from vector to ND
        slice.Convert2DtoND(*outputs2[0].GetWritableMatrix(), false);
        outMatrix1->SliceLHS(sliceArgs, slice);

        // advance slice indices
        for (int j = 0; j < numDim; ++j)
        {
            if (j == dim)
            {
                continue;
            }

            // increment index j if possible
            if (sliceArgs[j].Scalar() < (int) dims[j]-1)
            {
                ++sliceArgs[j].Scalar();
                break;
            }

            // index j is maxed out, so reset and continue to j+1
            sliceArgs[j].Scalar() = 0;
        }
    }

    if (extraCycles)
    {
        std::vector<Currency> inputs2;
        inputs2.push_back(outMatrix1);
        inputs2.push_back(extraCycles);

        return oml_func(eval, inputs2, outputs);
    }

    outputs.push_back(outMatrix1);

    return true;
}

// ND support for circshift function
bool oml_MatrixN_circshift(EvaluatorInterface eval, const std::vector<Currency>& inputs, std::vector<Currency>& outputs,
                           const OML_func1 oml_func)
{
    size_t nargin = inputs.size();

    if (nargin != 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    if (!inputs[0].IsNDMatrix())
    {
        // throw OML_Error(OML_ERR_NDMATRIX, 1, OML_STR_VARIABLE);
    }

    if (!inputs[1].IsInteger() && !inputs[1].IsVector())
        throw OML_Error(OML_ERR_SCALARVECTOR, 2);

    // dimension output
   const hwMatrixN* matrix = inputs[0].MatrixN();
   const std::vector<int>& dims = matrix->Dimensions();
   hwMatrixN* outMatrix = new hwMatrixN(dims, matrix->Type());

    // shift each element
    const hwMatrix* shift = inputs[1].Matrix();
    int numShifts = shift->Size();
    int numElems = matrix->Size();
    int numDim = static_cast<int> (dims.size());

    std::vector<int> offset(numDim);
    std::vector<int> lhsMatrixIndex(numDim);

    for (int i = 0; i < numShifts; ++i)
    {
        offset[i] = static_cast<int>((*shift)(i)) % dims[i];

        if (offset[i] < 0)
            offset[i] += dims[i];

        lhsMatrixIndex[i] = offset[i];
    }

    for (int i = 0; i < numElems; ++i)
    {
        if (matrix->IsReal())
            (*outMatrix)(lhsMatrixIndex) = (*matrix)(i);
        else
            outMatrix->z(lhsMatrixIndex) = matrix->z(i);

        // advance matrix indices
        for (int j = 0; j < numDim; ++j)
        {
            // increment j or wrap around
            if (lhsMatrixIndex[j] < dims[j] - 1)
                ++lhsMatrixIndex[j];
            else
                lhsMatrixIndex[j] = 0;

            if (lhsMatrixIndex[j] != offset[j])
                break;
        }
    }

    outputs.push_back(outMatrix);

    return true;
}

// Apply an operator function to each element of an ND matrix.
Currency oml_MatrixNUtil6(const Currency& op, OML_func2 oml_func)
{
    if (op.IsNDMatrix())
    {
        // throw OML_Error(OML_ERR_NDMATRIX, 1, OML_STR_VARIABLE);
    }

    const hwMatrixN* matrix = op.MatrixN();
    const std::vector<int>& dims = matrix->Dimensions();

    // convert from ND to vector
    hwMatrix* slice1D = new hwMatrix;
    matrix->ConvertNDto1D(*slice1D);

    // call computation function
    ExprTreeEvaluator dummy;
    Currency op_new(slice1D);
    Currency result;

    result = (dummy.*oml_func)(op_new);

    // convert from vector to ND
    hwMatrixN* outMatrix = new hwMatrixN;
    outMatrix->Convert1DtoND(*result.GetWritableMatrix(), dims);

    return outMatrix;
}

// Apply an arithmetic operator function to each element pair from the arguments, with at least one ND matrix.
Currency oml_MatrixNUtil7(const Currency& op1, const Currency& op2, const OML_func3 oml_func)
{
    static ExprTreeEvaluator dummy;
    hwMatrixN* outMatrix = new hwMatrixN;
    Currency result;

    if (op1.IsNDMatrix())
    {
        const hwMatrixN* matrix_1 = op1.MatrixN();
        const std::vector<int>& dims_1 = matrix_1->Dimensions();

        // convert from ND to vector
        hwMatrix* slice1D_1 = new hwMatrix;
        matrix_1->ConvertNDto1D(*slice1D_1);
        Currency op1_new(slice1D_1);

        if (op2.IsNDMatrix())
        {
            const hwMatrixN* matrix_2 = op2.MatrixN();
            const std::vector<int>& dims_2 = matrix_2->Dimensions();

            // dims check here
            if (dims_1 != dims_2)
            {
                throw OML_Error(OML_ERR_ARRAYSIZE);
            }

            // convert from ND to vector
            hwMatrix* slice1D_2 = new hwMatrix;
            matrix_2->ConvertNDto1D(*slice1D_2);
            Currency op2_new(slice1D_2);
            result = (dummy.*oml_func)(op1_new, op2_new);
        }
        else if (op2.IsScalar() || op2.IsComplex())
        {
            result = (dummy.*oml_func)(op1_new, op2);
        }
        else
        {
            throw OML_Error(HW_ERROR_UNSUPOP);
        }

        // convert from vector to ND
        outMatrix->Convert1DtoND(*result.GetWritableMatrix(), dims_1);
    }
    else if ((op1.IsScalar() || op1.IsComplex()) && op2.IsNDMatrix())
    {
        const hwMatrixN* matrix_2 = op2.MatrixN();
        const std::vector<int>& dims_2 = matrix_2->Dimensions();

        // convert from ND to vector
        hwMatrix* slice1D_2 = new hwMatrix;
        matrix_2->ConvertNDto1D(*slice1D_2);
        Currency op2_new(slice1D_2);

        result = (dummy.*oml_func)(op1, op2_new);

        // convert from vector to ND
        outMatrix->Convert1DtoND(*result.GetWritableMatrix(), dims_2);
    }
    else
    {
        throw OML_Error(HW_ERROR_UNSUPOP);
    }

    return outMatrix;
}

// Apply a comparison operator function to each element pair from the arguments, with at least one ND matrix.
Currency oml_MatrixNUtil8(const Currency& op1, const Currency& op2, int op, const OML_func4 oml_func)
{
    ExprTreeEvaluator dummy;
    hwMatrixN* outMatrix = new hwMatrixN;
    Currency result;

    if (op1.IsNDMatrix())
    {
        const hwMatrixN* matrix_1 = op1.MatrixN();
        const std::vector<int>& dims_1 = matrix_1->Dimensions();

        // convert from ND to vector
        hwMatrix* slice1D_1 = new hwMatrix;
        matrix_1->ConvertNDto1D(*slice1D_1);
        Currency op1_new(slice1D_1);

        if (op2.IsNDMatrix())
        {
            const hwMatrixN* matrix_2 = op2.MatrixN();
            const std::vector<int>& dims_2 = matrix_2->Dimensions();

            // dims check here
            if (dims_1 != dims_2)
            {
                throw OML_Error(OML_ERR_ARRAYSIZE);
            }

            // convert from ND to vector
            hwMatrix* slice1D_2 = new hwMatrix;
            matrix_2->ConvertNDto1D(*slice1D_2);
            Currency op2_new(slice1D_2);
            result = (dummy.*oml_func)(op1_new, op2_new, op);
        }
        else if (op2.IsScalar() || op2.IsComplex())
        {
            result = (dummy.*oml_func)(op1_new, op2, op);
        }
        else
        {
            throw OML_Error(HW_ERROR_UNSUPOP);
        }

        // convert from vector to ND
        outMatrix->Convert1DtoND(*result.GetWritableMatrix(), dims_1);
    }
    else if ((op1.IsScalar() || op1.IsComplex()) && op2.IsNDMatrix())
    {
        const hwMatrixN* matrix_2 = op2.MatrixN();
        const std::vector<int>& dims_2 = matrix_2->Dimensions();

        // convert from ND to vector
        hwMatrix* slice1D_2 = new hwMatrix;
        matrix_2->ConvertNDto1D(*slice1D_2);
        Currency op2_new(slice1D_2);

        result = (dummy.*oml_func)(op1, op2_new, op);

        // convert from vector to ND
        outMatrix->Convert1DtoND(*result.GetWritableMatrix(), dims_2);
    }
    else
    {
        throw OML_Error(HW_ERROR_UNSUPOP);
    }

    return outMatrix;
}

// Apply a bool function to a pair of arguments, with at least one ND matrix.
bool oml_MatrixNUtil9(const Currency& op1, const Currency& op2, OML_func5 oml_func)
{
    bool result;

    if (op1.IsNDMatrix())
    {
        const hwMatrixN* matrix_1 = op1.MatrixN();
        const std::vector<int>& dims_1 = matrix_1->Dimensions();

        // convert from ND to vector
        hwMatrix* slice1D_1 = new hwMatrix;
        matrix_1->ConvertNDto1D(*slice1D_1);
        Currency op1_new(slice1D_1);

        if (op2.IsNDMatrix())
        {
            const hwMatrixN* matrix_2 = op2.MatrixN();
            const std::vector<int>& dims_2 = matrix_2->Dimensions();

            // dims check here
            if (dims_1 != dims_2)
            {
                throw OML_Error(OML_ERR_ARRAYSIZE);
            }

            // convert from ND to vector
            hwMatrix* slice1D_2 = new hwMatrix;
            matrix_2->ConvertNDto1D(*slice1D_2);
            Currency op2_new(slice1D_2);
            result = oml_func(op1_new, op2_new);
        }
        else if (op2.IsScalar() || op2.IsComplex())
        {
            result = oml_func(op1_new, op2);
        }
        else
        {
            throw OML_Error(HW_ERROR_UNSUPOP);
        }
    }
    else if ((op1.IsScalar() || op1.IsComplex()) && op2.IsNDMatrix())
    {
        const hwMatrixN* matrix_2 = op2.MatrixN();

        // convert from ND to vector
        hwMatrix* slice1D_2 = new hwMatrix;
        matrix_2->ConvertNDto1D(*slice1D_2);
        Currency op2_new(slice1D_2);

        result = oml_func(op1, op2_new);
    }
    else
    {
        throw OML_Error(HW_ERROR_UNSUPOP);
    }

    return result;
}

// Apply a bool function to a pair of arguments, with at least one ND matrix and a tolerance.
bool oml_MatrixNUtil10(const Currency& op1, const Currency& op2, const Currency& tol, OML_func6 oml_func)
{
    bool result;

    if (op1.IsNDMatrix())
    {
        const hwMatrixN* matrix_1 = op1.MatrixN();
        const std::vector<int>& dims_1 = matrix_1->Dimensions();

        // convert from ND to vector
        hwMatrix* slice1D_1 = new hwMatrix;
        matrix_1->ConvertNDto1D(*slice1D_1);
        Currency op1_new(slice1D_1);

        if (op2.IsNDMatrix())
        {
            const hwMatrixN* matrix_2 = op2.MatrixN();
            const std::vector<int>& dims_2 = matrix_2->Dimensions();

            // dims check here
            if (dims_1 != dims_2)
            {
                throw OML_Error(OML_ERR_ARRAYSIZE);
            }

            // convert from ND to vector
            hwMatrix* slice1D_2 = new hwMatrix;
            matrix_2->ConvertNDto1D(*slice1D_2);
            Currency op2_new(slice1D_2);
            result = oml_func(op1_new, op2_new, tol);
        }
        else if (!op2.IsMatrix())
        {
            result = oml_func(op1_new, op2, tol);
        }
        else
        {
            throw OML_Error(HW_ERROR_UNSUPOP);
        }
    }
    else if (!op1.IsMatrix() && op2.IsNDMatrix())
    {
        const hwMatrixN* matrix_2 = op2.MatrixN();

        // convert from ND to vector
        hwMatrix* slice1D_2 = new hwMatrix;
        matrix_2->ConvertNDto1D(*slice1D_2);
        Currency op2_new(slice1D_2);

        result = oml_func(op1, op2_new, tol);
    }
    else
    {
        throw OML_Error(HW_ERROR_UNSUPOP);
    }

    return result;
}

// Apply a function to each vector of a 2D matrix in the specified dimension. The
// syntax is oml_func(2D, dim). This utility applies to functions for which the
// output for each vector is a scalar.
Currency oml_MatrixUtil(EvaluatorInterface& eval, const hwMatrix* mtx, int dim, double (*func)(EvaluatorInterface&, const hwMatrix*))
{
    // todo: evolve function to be more parallel to oml_MatrixNUtil3, at least the signature
    hwMatrix* result;
    Currency rescur;
    switch (dim)
    {
        case 1:
        {
            result = EvaluatorInterface::allocateMatrix(1, mtx->N(), hwMatrix::REAL);
            rescur = result;
            for (int i = 0; i < mtx->N(); ++i)
            {
                std::unique_ptr<const hwMatrix> holder(EvaluatorInterface::allocateColumn(mtx, i));
                (*result)(i) = (*func)(eval, holder.get());
            }
            return rescur;
        }
        case 2:
        {
            result = EvaluatorInterface::allocateMatrix(mtx->M(), 1, hwMatrix::REAL);
            std::unique_ptr<hwMatrix> holder(EvaluatorInterface::allocateMatrix());
            rescur = result;
            for (int i = 0; i < mtx->M(); ++i)
            {
                BuiltInFuncsUtils::CheckMathStatus(eval, mtx->ReadRow(i, *holder));
                (*result)(i) = (*func)(eval, holder.get());
            }
            return rescur;
        }
        default:
        {
            result = EvaluatorInterface::allocateMatrix(mtx->M(), mtx->N(), hwMatrix::REAL);
            std::unique_ptr<hwMatrix> holder(EvaluatorInterface::allocateMatrix(1, 1, mtx->Type()));
            rescur = result;
            for (int i = 0; i < mtx->Size(); ++i)
            {
                if (holder->IsReal())
                    (*holder)(0) = (*mtx)(i);
                else
                    holder->z(0) = mtx->z(i);

                (*result)(i) = (*func)(eval, holder.get());
            }
            return rescur;
        }
    }
}
