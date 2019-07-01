/**
* @file MatrixNDisplay.h
* @date July 2016
* Copyright (C) 2016-2019 Altair Engineering, Inc.  
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

#ifndef __MATRIXNDISPLAY__
#define __MATRIXNDISPLAY__

// Begin defines/includes
#include "Hml2Dll.h"
#include "Currency.h"
#include "CurrencyDisplay.h"
#include "hwMatrixN.h"

// End defines/includes

//------------------------------------------------------------------------------
//!
//! \class MatrixNDisplay
//! \brief Displays ND matrices
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS MatrixNDisplay : public CurrencyDisplay
{
public:
    friend class Currency; //! Only currency is allowed to construct

    //!
    //! Destructor
    //!
    virtual ~MatrixNDisplay() {}

    //!
    //! Gets output - called from Currency::GetOutputString
    //! \param fmt Output format
    //! \param os  Output stream
    //!
    std::string GetOutput(const OutputFormat* fmt,
                          std::ostringstream& os) const;
    //!
    //! Gets number of rows and cols in given currency
    //! \param rows Number of rows
    //! \param cols Number of columns
    //!
    virtual void GetCurrencySize(int& rows, 
                                 int& cols) const;
    //!
    //! Helper to slice a given ND matrix to 2D matrices
    //! \param mtx    Given matrix
    //! \param curs   Vector of currencies of 2D component matrices
    //! \param labels Vector of slice labels
    //!
    static void GetSlices(const hwMatrixN*          mtx,
                          std::vector<Currency>&    curs,
                          std::vector<std::string>& labels);
    //!
    //! Returns true if parent is ND matrix
    //!
    virtual bool IsNDMatrixDisplay() const { return true; }
protected:
    //!
    //! Sets data for forward pagination
    //!
    virtual void SetForwardDisplayData();
    //!
    //! Sets data for back pagination
    //!
    virtual void SetBackDisplayData();

    //!
    //! Gets number of rows that can be fit
    //!
    virtual int GetNumRowsToFit() const;

    //!
    //! Gets values as a string
    //! \param fmt Output format
    //!
    virtual std::string GetValues( const OutputFormat* fmt) const;

private:
    //!
    //! Constructor - Only currency is allowed to construct
    //! \param cur Currency associated with this display
    //!
    MatrixNDisplay(const Currency& cur);

    MatrixNDisplay();                                       // Stubbed out 
    MatrixNDisplay(            const MatrixNDisplay& src) ; // Stubbed out 
    MatrixNDisplay& operator=( const MatrixNDisplay& src);  // Stubbed out

    //!
    //! Gets output with no pagination
    //! \param fmt Format
    //!
    std::string GetOutputNoPagination(const OutputFormat* fmt) const;
    //!
    //! Gets outputfor forward/down pagination
    //! \param fmt Format
    //!
    std::string GetOutputForwardPagination(const OutputFormat* fmt) const;
    //!
    //! Gets outputfor back/up pagination
    //! \param fmt Format
    //!
    std::string GetOutputBackPagination(const OutputFormat* fmt) const;

    //!
    //! Helper to slice a given ND matrix to 2D matrices
    //! \param mtx    Given matrix
    //! \param slices Slices
    //! \param curs   Vector of currencies of 2D component matrices
    //! \param labels Vector of slice labels
    static void GetSlicesHelper(const hwMatrixN*          mtx,
                                std::vector<hwSliceArg>   slices,
                                std::vector<Currency>&    curs,
                                std::vector<std::string>& sliceLabels);
    //!
    //! Updates the number of rows to fit
    //!
    void UpdateNumLinesPrinted();
    //!
    //! True if this ND matrix has only empty slices
    //! \param slices Slices
    //!
    bool HasOnlyEmptySlices(const std::vector<Currency>& slices) const;
    //!
    //! Gets display string if matrix ND has only empty slices
    //! 
    std::string GetOutputEmpty() const;
};
#endif

// End of file:

