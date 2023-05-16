/**
* @file omlmatio.h
* @date August 2020
* Copyright (C) 2020-2023 Altair Engineering, Inc.
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
#ifndef __OMLMATIO_H__
#define __OMLMATIO_H__

#include "MatioTboxDefs.h"

#include <map>
#include <string>

#include "Currency.h"
#include "EvaluatorInt.h"

struct matvar_t;

//------------------------------------------------------------------------------
//! \class OmlMatio
//! \brief Interfaces with matio library for loading/saving variables
//------------------------------------------------------------------------------
class MATIOOMLTBOX_DECLS OmlMatio
{
public:
    //!
    //! \enum: Mat file version
    //!
    enum MATFILEVERSION
    {
        MATFILEVERSION_73,   //!< Version 7.3
        MATFILEVERSION_5     //!< Version 5
    };
    //!
    //! Constructor
    //! \param Filename
    //! \param Verbose level
    //! \param Version
    //!
    OmlMatio(const std::string&, int, MATFILEVERSION);
    //!
    //! Destructor
    //!
    ~OmlMatio() {}

    //!
    //! Adds warning
    //! \param Warning
    //!
    void AddWarning(const std::string&);
    //!
    //! Gets warning
    //!
    std::string GetWarning();

    //!
    //! Gets variable name, creating one if needed
    //! \param Mat variable
    //! \param True if warning about auto generated name needs to be created
    //!
    std::string GetName(matvar_t*, bool = true);

    //!
    //! Converts mat variable to currency
    //! \param Mat variable
    //! \param Evaluator interface
    //!
    Currency MatVarToCurrency(matvar_t*,
                              EvaluatorInterface);
    //!
    //! Converts currency to matvar
    //! \param Currency
    //!
    matvar_t* CurrencyToMatVar(const char* name, const Currency&);

    //!
    //! Sets the version
    //! \param Version
    //!
    void SetVersion(MATFILEVERSION v) { _version = v; }
    //!
    //! Gets type description
    //
    std::string GetTypeString(matvar_t*);  
    //!
    //! Resets data
    //! 
    void Reset();
    //!
    //! Sets the name prefix
    //! 
    void SetNamePrefix(const std::string&);
    //!
    //! Helper method to get children
    //! \param Evaluator
    //! \param Parent
    //! 
    std::vector<matvar_t*> GetChildren(EvaluatorInterface, matvar_t*);

private:
    std::string    _warn;       //!< Warning
    std::string    _name;       //!< Used in name if matvar_t->name is null
    int            _idx;        //!< Used in name if matvar_t->name is null
    int            _verbose;    //!< Verbose level
    MATFILEVERSION _version;    //!< Matio version

    OmlMatio();                             // Stubbed out 
    OmlMatio(const OmlMatio&);              // Stubbed out 
    OmlMatio& operator=(const OmlMatio&);   // Stubbed out


    // Conversion methods
    //!
    //! Converts 2D/ND cell array to currency
    //! \param Variable
    //! \param Evaluator Interface
    //!
    Currency CellArray2Currency(matvar_t*,
                                EvaluatorInterface);
    //!
    //! Converts char array to currency
    //! \param Variable
    //!
    Currency Char2Currency(matvar_t*);
    //!
    //! Converts double to currency
    //! \param Variable
    //!
    Currency Double2Currency(matvar_t*);
    //!
    //! Converts real integer data to currency
    //! \param Variable
    //!
    Currency IntComplex2Currency(matvar_t*);
    //!
    //! Converts real integer data to currency
    //! \param Variable
    //!
    Currency IntReal2Currency(matvar_t*);
    //!
    //! Converts 2D matrix to currency
    //! \param Variable
    //! \param Name
    //!
    Currency Matrix2Currency(matvar_t*,
                             const std::string&);
    //!
    //! Converts number to currency
    //! \param Variable
    //! \param Name
    //!
    Currency Number2Currency(matvar_t*,
                             const std::string&);
    //!
    //! Converts structure array class to currency
    //! \param Variable
    //!
    Currency Sparse2Currency(matvar_t*);
    //!
    //! Converts structure array class to currency
    //! \param Variable
    //! \param Evaluator interface
    //!
    Currency Struct2Currency(matvar_t*,
                             EvaluatorInterface);

    //!
    //! Gets dummy variable of value NaN
    //! \param Name
    //!
    matvar_t* GetDummyMatVar(const std::string& name);
    //!
    //! Converts currency to matvar of type cell
    //! \param Name,
    //! \param Currency
    //!
    matvar_t* Currency2MatCell(const char*,
                               const Currency&);
    //!
    //! Converts currency to matvar of type ND cell
    //! \param Name,
    //! \param Currency
    //!
    matvar_t* Currency2MatCellND(const char*,
                                 const Currency&);
    //!
    //! Converts currency to matvar of type matrix
    //! \param Name,
    //! \param Currency
    //!
    matvar_t* Currency2MatMatrix(const char*, 
                                 const Currency&);
    //!
    //! Converts currency to matvar of type matrix ND
    //! \param Name,
    //! \param Currency
    //!
    matvar_t* Currency2MatMatrixND(const char*,
                                   const Currency&);
    //!
    //! Converts currency to matvar_t of type sparse matrix
    //! \param Name
    //! \param Currency
    //!
    matvar_t* Currency2MatSparse(const char*     name, 
                                 const Currency& cur);
    //!
    //! Converts currency to matvar of type struct
    //! \param Name
    //! \param Currency
    //!
    matvar_t* Currency2MatStruct(const char*     name,
                                 const Currency& cur);

    //!
    //! Converts currency to matvar of type object
    //! \param Name
    //! \param Currency
    //!
    matvar_t* Currency2MatObj(const char*     name,
                              const Currency& cur);

    // Error handling
    //!
    //! Returns empty currency and adds invalid dimension warning
    //! \param Variable
    //! \param Name
    //!
    Currency HandleInvalidDims(matvar_t*, 
                               const std::string&);
    //!
    //! Returns empty currency and adds warning about unsupported type
    //! \param Variable
    //! \param Name
    //!
    Currency HandleInvalidType(matvar_t*,
                               const std::string&);
    //!
    //! Returns empty currency and adds warning about invalid rank
    //! \param Variable
    //! \param Name
    //!
    Currency HandleInvalidRank(matvar_t*,
                               const std::string&);

    //!
    //! Adds warning about invalid dimensions
    //! \param Variable
    //! \param Name
    //!
    void AddInvalidDimsWarning(matvar_t*, const std::string&);
    //!
    //! Adds warning about invalid rank
    //! \param Variable
    //! \param Name
    //!
    void AddInvalidRankWarning(matvar_t*, const std::string&);

};

#endif


