/**
* @file BuiltInFuncsMKL.h
* @date September 2021
* Copyright (C) 2016-2021 Altair Engineering, Inc.
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

#ifndef __BUILTINFUNCSMKL__
#define __BUILTINFUNCSMKL__

// Begin defines/includes
#include "EvaluatorInt.h"

// End defines/includes

//------------------------------------------------------------------------------
//!
//! \brief Class for built-in functions implementing MKL based functions
//!
//------------------------------------------------------------------------------
class OMLDLL_DECLS BuiltInFuncsMKL
{
public:
    //!
    //! Destructor
    //!
    ~BuiltInFuncsMKL() {}
    //!
    //! Returns cosine of input in radians [cos command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Cos(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns sine of input in radians [sin command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Sin(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns tangent of input in radians [tan command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Tan(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns secant of input in radians [sec command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Sec(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns cosecant of input in radians [csc command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Csc(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns cotangent of input in radians [cot command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Cot(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns cosine of input in degrees [cosd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool CosD(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns sine of input in degrees [sind command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool SinD(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns tangent of input in degrees [tand command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool TanD(EvaluatorInterface          eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns secant of input in degrees [secd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool SecD(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns cosecant of input in degrees [cscd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool CscD(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns cotangent of input in degrees [cotd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool CotD(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cosine of input in radians [acos command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCos(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse sine of input in radians [asin command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aSin(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse tangent of input in radians [atan command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aTan(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse secant of input in radians [asec command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aSec(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cosecant of input in radians [acsc command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCsc(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cotangent of input in radians [acot command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCot(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cosine of input in degrees [acosd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCosD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse sine of input in degrees [asind command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aSinD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse tangent of input in degrees [atand command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aTanD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse secant of input in degrees [asecd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aSecD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cosecant of input in degrees [acscd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCscD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse cotangent of input in degrees [acotd command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCotD(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse tangent of input in radians [atan2 command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aTan2(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse tangent of input in degrees [atan2d command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aTan2D(EvaluatorInterface           eval,
                       const std::vector<Currency>& inputs,
                       std::vector<Currency>&       outputs);
    //!
    //! Returns hyperbolic cosine of input [cosh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Cosh(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns hyperbolic sine of input [sinh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Sinh(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns hyperbolic tangent of input [tanh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Tanh(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);
    //!
    //! Returns inverse hyperbolic sine of input [asinh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aSinh(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns exponential of input [exp command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Exp(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! Returns inverse hyperbolic cosine of input [acosh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aCosh(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);
    //!
    //! Returns inverse hyperbolic tangent of input [atanh command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool aTanh(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);

    //!
    //! Returns natural logarithm of input [log command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Log(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Returns logarithm base 2 of input [log2 command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Log2(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);

    //!
    //! Returns logarithm base 10 of input [log10 command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Log10(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);

    //!
    //! Returns square root of input [sqrt command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Sqrt(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);

    //!
    //! Returns magnitude of input [abs command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Abs(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Returns angle of input [abs command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Arg(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Rounds input to the nearest integer [round command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Round(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);

    //!
    //! Rounds input toward inf [ceil command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Ceil(EvaluatorInterface           eval,
                     const std::vector<Currency>& inputs,
                     std::vector<Currency>&       outputs);

    //!
    //! Rounds input input toward -inf [floor command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Floor(EvaluatorInterface           eval,
                      const std::vector<Currency>& inputs,
                      std::vector<Currency>&       outputs);

    //!
    //! Rounds input toward zero [fix command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Fix(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Computes division remainder [rem command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Rem(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Computes modulo function [mod command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Mod(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Returns element-wise maximum of two inputs [max command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Max(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);

    //!
    //! Returns element-wise minimum of two inputs [min command]
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    static bool Min(EvaluatorInterface           eval,
                    const std::vector<Currency>& inputs,
                    std::vector<Currency>&       outputs);
    //!
    //! \param eval    Evaluator interface
    //! \param inputs  Vector of inputs
    //! \param outputs Vector of outputs
    //!
    // Set pivot threshold for MKL sparse matrix division
    static bool oml_mkl_sdpt(EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs);

private:
    //!
    //! Constructor
    //!
    BuiltInFuncsMKL() {}
};

#endif // __BUILTINFUNCSMKL__
