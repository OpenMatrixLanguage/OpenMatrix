/**
* @file MKL_Tutilities.h
* @date March 2022
* Copyright (C) 2022 Altair Engineering, Inc.
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

#ifndef __MKL_UTILITY_FUNCS__
#define __MKL_UTILITY_FUNCS__

// forward declarations
class hwMathStatus;
template <typename T> class hwTComplex;
template <typename T1, typename T2> class hwTMatrix;
template <typename T1, typename T2> class hwTMatrixS;
typedef hwTComplex<double> hwComplex;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;
typedef hwTMatrixS<double, hwTComplex<double> > hwMatrixS;

//------------------------------------------------------------------------------
//!
//! \brief Class to wrap MKL based utility functions
//!
//------------------------------------------------------------------------------
template < typename T1, typename T2 = hwTComplex<T1> > class MKL_Tutilities
{
public:
    //! Destructor
    ~MKL_Tutilities() {}

    // Return by reference
    static void Cos(const hwTMatrix<T1, T2>& X,
                    hwTMatrix<T1, T2>&       C);
    static void Sin(const hwTMatrix<T1, T2>& X,
                    hwTMatrix<T1, T2>&       S);
    static void Exp(const hwTMatrix<T1, T2>& X,
                    hwTMatrix<T1, T2>&       E);
    static void Add(const hwTMatrix<T1, T2>& A,
                    const hwTMatrix<T1, T2>& B,
                    hwTMatrix<T1, T2>&       sum);
    static void Subtr(const hwTMatrix<T1, T2>& A,
                      const hwTMatrix<T1, T2>& B,
                      hwTMatrix<T1, T2>&       diff);
    static void PowerByElems(const hwTMatrix<T1, T2>& A,
                             T1                       P,
                             hwTMatrix<T1, T2>&       B);
    static void PowerByElems(const hwTMatrix<T1, T2>&  A,
                             const T2&                 P,
                             hwTMatrix<T1, T2>&        B);
    static void PowerByElems(const hwTMatrix<T1, T2>& A,
                             const hwTMatrix<T1, T2>& P,
                             hwTMatrix<T1, T2>&       B);
    static void MultByElems(const hwTMatrix<T1, T2>& A,
                            const hwTMatrix<T1, T2>& B,
                            hwTMatrix<T1, T2>&       prod);
    static void DivideByElems(const hwTMatrix<T1, T2>& A,
                              const hwTMatrix<T1, T2>& B,
                              hwTMatrix<T1, T2>&       quot);
    static void Mod(const hwTMatrix<T1, T2>& A,
                    T1                       b,
                    hwTMatrix<T1, T2>&       mod);
    static void Mod(const hwTMatrix<T1, T2>& A,
                    const hwTMatrix<T1, T2>& B,
                    hwTMatrix<T1, T2>&       mod);
    static void Floor(const hwTMatrix<T1, T2>& A,
                      hwTMatrix<T1, T2>&       F);
    static void Conj(const hwTMatrix<T1, T2>& A,
                     hwTMatrix<T1, T2>& conj);
    static void Hypot(const hwTMatrix<T1, T2>& A,
                      const hwTMatrix<T1, T2>& B,
                      hwTMatrix<T1, T2>& hyp);
    static void SparseAdd(const hwTMatrixS<T1, T2>& A,
                          const hwTMatrixS<T1, T2>& B,
                          hwTMatrixS<T1, T2>&       sum);
    static void SparseMult(const hwTMatrixS<T1, T2>& A,
                           const hwTMatrix<T1, T2>&  B,
                           hwTMatrix<T1, T2>&        prod);
    static void SparseMult(const hwTMatrix<T1, T2>&  A,
                           const hwTMatrixS<T1, T2>& B,
                           hwTMatrix<T1, T2>&        prod);
    static void SparseMult(const hwTMatrixS<T1, T2>& A,
                           const hwTMatrixS<T1, T2>& B,
                           hwTMatrixS<T1, T2>&       prod);
    static void SparseDivideLeft(const hwTMatrixS<T1, T2>& A,
                                 const hwTMatrix<T1, T2>&  B,
                                 hwTMatrix<T1, T2>&        Q);
    static void SparseDivideRight(const hwTMatrix<T1, T2>&  A,
                                  const hwTMatrixS<T1, T2>& B,
                                  hwTMatrix<T1, T2>&        Q);
    static void SetIparam9(int _iparm9);

    // Return by value
    static hwTMatrix<T1, T2> Cos(const hwTMatrix<T1, T2>& X);
    static hwTMatrix<T1, T2> Sin(const hwTMatrix<T1, T2>& X);
    static hwTMatrix<T1, T2> Exp(const hwTMatrix<T1, T2>& E);
    static hwTMatrix<T1, T2> Add(const hwTMatrix<T1, T2>& A,
                                 const hwTMatrix<T1, T2>& B);
    static hwTMatrix<T1, T2> Subtr(const hwTMatrix<T1, T2>& A,
                                   const hwTMatrix<T1, T2>& B);
    static hwTMatrix<T1, T2> PowerByElems(const hwTMatrix<T1, T2>& A,
                                          T1                       P);
    static hwTMatrix<T1, T2> PowerByElems(const hwTMatrix<T1, T2>& A,
                                          const T2&                P);
    static hwTMatrix<T1, T2> PowerByElems(const hwTMatrix<T1, T2>& A,
                                          const hwTMatrix<T1, T2>& P);
    static hwTMatrix<T1, T2> MultByElems(const hwTMatrix<T1, T2>& A,
                                         const hwTMatrix<T1, T2>& B);
    static hwTMatrix<T1, T2> DivideByElems(const hwTMatrix<T1, T2>& A,
                                           const hwTMatrix<T1, T2>& B);
    static hwTMatrix<T1, T2> Mod(const hwTMatrix<T1, T2>& A,
                                 T1                       b);
    static hwTMatrix<T1, T2> Mod(const hwTMatrix<T1, T2>& A,
                                 const hwTMatrix<T1, T2>& B);
    static hwTMatrix<T1, T2> Floor(const hwTMatrix<T1, T2>& A);

private:
    //!
    //! Constructor
    //!
    MKL_Tutilities() {}
    //!
    //! Stubbed out copy constructor
    //!
    MKL_Tutilities(const MKL_Tutilities&);
    //!
    //! Stubbed out assignment operator
    //!
    MKL_Tutilities& operator=(const MKL_Tutilities&);
};

//! template implementation file
#include "MKL_Tutilities.cc"

#endif // __MKL_UTILITY_FUNCS__
