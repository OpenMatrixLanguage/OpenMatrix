/**
* @file MatrixDisplay.h
* @date November, 2015
* Copyright (C) 2015-2024 Altair Engineering, Inc.  
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

#ifndef __MATRIXDISPLAY_H__
#define __MATRIXDISPLAY_H__

// Begin defines/includes
#include "OMLDll.h"
#include "Currency.h"
#include "CurrencyDisplay.h"

#include <string>
#include <vector>

#include "hwComplex.h"

template <typename T1, typename T2> class hwTMatrix;
typedef hwTMatrix<double, hwTComplex<double> > hwMatrix;

class OutputFormat;

// End defines/includes
//------------------------------------------------------------------------------
//!
//! \class MatrixDisplay
//! \brief Displays formatted matrix output
//!
//------------------------------------------------------------------------------
class OMLDLL_DECLS MatrixDisplay : public CurrencyDisplay
{
public:
    friend class Currency; //! Only currency is allowed to construct

    //!
    //! Destructor 
    //!
    virtual ~MatrixDisplay() {}

    //!
    //! Initialize
    //! \param Output format   
    //! \param Interpreter
    //! \param Parent display
    virtual void Initialize(const OutputFormat*, Interpreter*, CurrencyDisplay* = 0) override;
    //!
    //! True if columns can be paginated
    //!
    virtual bool CanPaginateColumns() const override { return true; }
    //!
    //! Gets number of rows and cols in matrix
    //! \param Number of rows
    //! \param Number of columns
    //!
    virtual void GetCurrencySize(int&, int&) const override;
    //!
    //! Sets indices for back pagination
    //!
	virtual void SetBackDisplayData() override;
    //!
	//! Sets indices for forward pagination
    //!
	virtual void SetForwardDisplayData() override;
    //!
	//! Sets indices for right pagination
    //!
	virtual void SetRightDisplayData() override;
    //!
	//! Sets indices for left pagination
    //!
	virtual void SetLeftDisplayData() override;
    //!
	//! Sets indices for down pagination
    //!
	virtual void SetDownDisplayData() override;
    //!
	//! Sets indices for up pagination
    //!
	virtual void SetUpDisplayData() override;

    //!
    //! Utility which returns matrix values as string
    //! \param Input currency
    //! \param Format
    //! \param Delimiter for rows
    //! \param Delimiter for columns
    //! \param Custom precision for real part
    //! \param Custom precision for imaginary part
    //! \param Number of columns which need cdelim prepended
    //!   
    static std::string GetOutputValues(const Currency&, const OutputFormat*, const std::string&, const std::string&, const std::string&, const std::string&, int);
    //!
    //! Utility which returns matrix values as string, without formatting output
    //! \param Input currency
    //! \param Row delimiter
    //! \param Column delimiter
    //! \param Number of columns which need col delim prepended - col offset
    //!
    static void WriteNonFormattedOutputValues(const Currency&, const std::string&, const std::string&, int, std::FILE*);
    //!
    //! Gets values as a string
    //! \param Output format
    //!
    virtual std::string GetValues(const OutputFormat*) const override;
    //!
    //! Returns true if end of pagination message needs to be printed
    //! \param Additional message that needs to be printed
    virtual bool GetPaginationEndMsg(std::string&) const override;

private:
    mutable bool             _haslargeint;     //!< True if matrix has large ints
    mutable bool             _uppercase;       //!< True if scientific uppercase
    mutable int              _formatinteger;   //!< Integer part for format, if applicable
    mutable int              _formatdecimal;   //!< Decimal part for format, if applicable
    mutable DisplayFormat    _displayformat;   //!< Display format for matrix
    mutable std::streamsize  _precision;       //!< Precision
    mutable std::vector<int> _realwidth;       //!< Widest real value width/column
    mutable std::vector<int> _imagwidth;       //!< Widest imag value width/column
    mutable std::string      _delimiter;       //!< Delimiter between cols
    mutable int              _skipformatwidth; //!< Width for skip format option

    mutable long long   _maxdigits;          //!< Format - max digits
    mutable long double _maxfloat;           //!< Format - max possible float
    mutable long long   _maxint;             //!< Format - max possible int
    mutable long double _minfloat;           //!< Format - min possible float
    mutable size_t      _totaldigits;        //!< Format - total digits

    //!
    //! Constructor - Only currency is allowed to construct
    //! \param cur Currency associated with this display
    //!
    MatrixDisplay(const Currency& cur);
    
    MatrixDisplay();                                      // Stubbed out 
    MatrixDisplay(            const MatrixDisplay& src) ; // Stubbed out 
    MatrixDisplay& operator=( const MatrixDisplay& src);  // Stubbed out
    
    //!
	//! Gets pagination info for printing
    //! \param Number of rows
    //! \param Number of columns
    //!
	std::string GetPaginationHeader(int, int) const;
    //!
    //! Gets output
    //! \param Output format
    //! \param Output stream
    //!
    std::string GetOutput(const OutputFormat*, std::ostringstream&) const override;
    //!
	//! Gets matrix data with no pagination - using defaults
	//! \param Format
    //!
	std::string GetOutputNoPagination(const OutputFormat*) const;
    //!
	//! Gets matrix data with no pagination
	//! \param Format
    //! \param Delimiter for rows
    //! \param Appends a new line to the start of a row
    //! \param Number of columns which need cdelim prepended
    //!
	std::string GetOutputNoPagination(const OutputFormat*, const std::string&, bool, int) const;
    //!
	//! Gets matrix data with forward pagination
	//! \param Format
    //!
	std::string GetOutputForwardPagination(const OutputFormat*) const;
    //!
	//! Gets matrix data with back pagination
	//! \param Format
    //!
	std::string GetOutputBackPagination(const OutputFormat*) const; 
    //!
    //! Gets output
    //! \param Row
    //! \param Column
    //! \param Width of real field
    //! \param Width of imaginary field
    //! \param True if this is a real mtx
    //! \param True if this is a complex mtx with no imaginary parts
    //! \param Output
    //!
    void GetOutput(int, int, int, int, bool, bool, std::string&) const;
    //!
    //! Gets number of columns to fit
    //! \param Column to start from
    //! \param Total number of columns
    //! \param True if forward paginating
    //! \param True if there is only real data
    //!
    int GetNumColumnsToFit(int, int, bool, bool) const;
    //!
    //! Converts real value to formatted string
    //! \param Value
    //!
    std::string RealToString(double val) const { return RealToString(val, _displayformat); }
    //!
    //! Converts real value to formatted string
    //! \param Value
    //! \param Format
    //!
    std::string RealToString(double, DisplayFormat) const; 
    //!
    //! Utility which returns real value to string
    //! \param Value
    //! \param Custom precision
    //!
    static std::string RealToString(double, const std::string&);
    //!
    //! Gets format for given double
    //! \param Given value
    //!
    DisplayFormat GetFormat(double) const;
    //!
    //! Returns format for complex numbers
    //! \param Given value
    //!
    DisplayFormat GetFormat(const hwComplex&) const;
    //!
    //! Sets matrix format
    //! \param Format
    //! \param Interpreter
    //!
    void SetFormat(const OutputFormat*, Interpreter*);
    //!
    //! Resets format
    //!
    void ResetFormat();

    //!
    //! Scan matrix and set the width of different columns
    //! \param Interpreter
    //!
    void SetWidth(Interpreter*);
    //!
    //! Sets the delimiter
    //! \param Delimiter to set, if empty, default delim will be used
    //!
    void SetDelimiter(const std::string&);
    //!
    //! Returns true if matrix was paginating
    //!
    virtual bool WasPaginating()  const;
    //!
    //! Returns true if matrix is paginating
    //!
    bool IsPaginating()  const;
    //!
    //! Returns formatted string length
    //! \param value
    //! \param display format type
    //!
    size_t GetFormattedStringLength(double, DisplayFormat) const;

    //!
    //! Gets real component width
    //! \param index
    int GetRealWidth(int) const;
    //!
    //! Gets imaginary component width
    //! \param index
    //! \param True if this is real data
    //!
    int GetImagWidth(int, bool) const;
};

#endif
// End of file:
