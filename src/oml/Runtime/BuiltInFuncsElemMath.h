/**
* @file BuiltInFuncsElemMath.h
* @date July 2016
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

#ifndef __BUILTINFUNCSELEMMATH__
#define __BUILTINFUNCSELEMMATH__

// Begin defines/includes
#include <deque>

#include "EvaluatorInt.h"

#include "hwComplex.h"
// End defines/includes

//------------------------------------------------------------------------------
//!
//! \brief Class for built-in functions implementing elementary math commands
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS BuiltInFuncsElemMath
{
public:
    //!
    //! Destructor
    //!
    ~BuiltInFuncsElemMath() {}
    //!
    //! Returns true after flipping matrix [flip command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Flip( EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns true after flipping matrix from left to right [fliplr command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Fliplr( EvaluatorInterface           eval,
                        const std::vector<Currency>& inputs,
                        std::vector<Currency>&       outputs);
    //!
    //! Returns true after flipping matrix from up to down [flipud command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Flipud( EvaluatorInterface           eval,
                        const std::vector<Currency>& inputs,
                        std::vector<Currency>&       outputs);
    //!
    //! Helper function for unique command for 2D matrices
    //! \param eval      Evaluator interface
    //! \param x         Input currency
    //! \param cmprows   True if whole rows are compared
    //! \param forward   True if forward search
    //! \param outputIdx True if idx(i) of each element of output(y) in given 
    //!                  curr(x) is returned such that y = x(i)
    //! \param inputIdx  True if idx(j) of each element of input(x) in 
    //!                  output(y) such that x = y(j)
    //! \param outputs   Vector of outputs
    //!
    static void UniqueHelperFuncMtx( EvaluatorInterface&    eval,
                                     const Currency&        x,
                                     bool                   forward,
                                     bool                   outputIdx,
                                     bool                   inputIdx,
                                     std::vector<Currency>& outputs);
    //!
    //! Helper function for unique command for ND matrices
    //! \param eval      Evaluator interface
    //! \param x         Input currency
    //! \param cmprows   True if whole rows are compared
    //! \param forward   True if forward search
    //! \param outputIdx True if idx(i) of each element of output(y) in given 
    //!                  curr(x) is returned such that y = x(i)
    //! \param inputIdx  True if idx(j) of each element of input(x) in 
    //!                  output(y) such that x = y(j)
    //! \param outputs   Vector of outputs
    //!
    static void UniqueHelperFuncMtxN( EvaluatorInterface&    eval,
                                      const Currency&        x,
                                      bool                   forward,
                                      bool                   outputIdx,
                                      bool                   inputIdx,
                                      std::vector<Currency>& outputs);

private:
    //!
    //! Constructor
    //!
    BuiltInFuncsElemMath() {}
    //!
    //! Helper function for Unique with real Matrix
    //! \param mtx Input matrix
    //!
    std::deque<double> UniqueHelperRealMtx( const hwMatrix* mtx);
    //!
    //! Helper function for Unique with complex Matrix
    //! \param mtx Input matrix
    //!
    std::deque<hwComplex> UniqueHelperComplexMtx( const hwMatrix* mtx);
    //!
    //! Helper function for Unique with real MatrixN
    //! \param mtx Input matrix
    //!
    std::deque<double> UniqueHelperRealMtxN( const hwMatrixN* mtx);
    //!
    //! Helper function for Unique with complex MatrixN
    //! \param mtx Input matrix
    //!
    std::deque<hwComplex> UniqueHelperComplexMtxN( const hwMatrixN* mtx);
    //!
    //! Gets indices of occurences of matrix elements(x) in values => y = x(i)
    //! \param x       Matrix to search in
    //! \param y       Given values
    //! \param forward True if forward search
    Currency GetMatrixIndices( const hwMatrix*           x, 
                               const std::deque<double>& y, 
                               bool                      forward);
    //!
    //! Gets indices of occurences of matrix elements(x) in values => y = x(i)
    //! \param x       Matrix to search in
    //! \param y       Given values
    //! \param forward True if forward search
    //!
    Currency GetMatrixIndices( const hwMatrix*              x, 
                               const std::deque<hwComplex>& y,
                               bool                         forward);
    //!
    //! Gets indices of occurences of unsorted values in a real matrix
    //! \param vals     Given values
    //! \param searchin Matrix to search in
    //! \param forward  True if searching forward
    //!
    Currency GetMatrixNIndices( const std::deque<double>& vals, 
                                const hwMatrixN*          in, 
                                bool                      forward);
    //!
    //! Gets indices of occurences of matrix elements(x) in values => y = x(i)
    //! \param x       Matrix to search in
    //! \param y       Given values
    //! \param forward True if forward search
    //!
    Currency GetMatrixNIndices( const hwMatrixN*             x,
                                const std::deque<hwComplex>& y,
                                bool                         forward);
    //!
    //! Gets indices of occurences of values in a complex matrix => x = y(j)
    //! \param x Matrix to search in
    //! \param y Given values
    //!
    Currency GetValueIndices( const hwMatrix*           x,
                              const std::deque<double>& y);
    //!
    //! Gets indices of occurences of values in a complex matrix => x = y(j)
    //! \param x Matrix to search in
    //! \param y Given values
    //!
    Currency GetValueIndices( const hwMatrix*              x,
                              const std::deque<hwComplex>& y);
    //!
    //! Gets indices of occurences of unsorted values in a real matrix
    //! \param vals Given values
    //! \param in   Matrix to search in
    //!
    Currency GetValueIndicesND( const std::deque<double>& vals, 
                                const hwMatrixN*          in);
    //!
    //! Gets indices of occurences of values in a complex matrix => x = y(j)
    //! \param x Matrix to search in
    //! \param y Given values
    //!
    Currency GetValueIndicesND( const hwMatrixN*             x,
                                const std::deque<hwComplex>& y);
    //!
    //! True if given ND matrix has a single row 
    //! \param mtx Given matrix
    //!
    bool IsSingleRowND( const hwMatrixN* mtx) const;
    //!
    //! Helper function for flip methods which returns flipped matrix 
    //! \param in  Input matrix
    //! \param dim Dimension to flip
    //!
    hwMatrix* FlipHelper( const hwMatrix* in, 
                          int             dim);
};

#endif // __BUILTINFUNCSELEMMATH__

