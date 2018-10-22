/**
* @file MatrixDisplay.h
* @date November, 2015
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

#ifndef __MATRIXDISPLAY_H__
#define __MATRIXDISPLAY_H__

// Begin defines/includes
#include "Hml2Dll.h"
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
class HML2DLL_DECLS MatrixDisplay : public CurrencyDisplay
{
public:
    friend class Currency; //! Only currency is allowed to construct

    //!
    //! Destructor 
    //!
    virtual ~MatrixDisplay() {}

    //!
    //! Initialize
    //! \param fmt    Output format   
    //! \param interp Interpreter
    //! \param parent Parent display
    virtual void Initialize(const OutputFormat* fmt,
                            Interpreter*        interp,
                            CurrencyDisplay*    parent = 0);
    //!
    //! True if columns can be paginated
    //!
    virtual bool CanPaginateColumns() const { return true; }
    //!
    //! Gets number of rows and cols in given currency
    //! \param cur  Given cell array/matrix
    //! \param rows Number of rows
    //! \param cols Number of columns
    //!
    virtual void GetCurrencySize(int& rows, 
                                 int& cols) const;
    //!
    //! Sets indices for back pagination
    //!
	virtual void SetBackDisplayData();
    //!
	//! Sets indices for forward pagination
    //!
	virtual void SetForwardDisplayData();
    //!
	//! Sets indices for right pagination
    //!
	virtual void SetRightDisplayData();
    //!
	//! Sets indices for left pagination
    //!
	virtual void SetLeftDisplayData();
    //!
	//! Sets indices for down pagination
    //!
	virtual void SetDownDisplayData();
    //!
	//! Sets indices for up pagination
    //!
	virtual void SetUpDisplayData();


    //!
    //! Utility which returns matrix values as string
    //! \param in            Input currency
    //! \param fmt           Format
    //! \param rdelim        Delimiter for rows
    //! \param cdelim        Delimiter for columns
    //! \param precisionreal Custom precision for real part
    //! \param precisionimag Custom precision for imaginary part
    //! \param coffset       Number of columns which need cdelim prepended
    //!   
    static std::string GetOutputValues(const Currency&     in,
                                       const OutputFormat* fmt,
                                       const std::string&  rdelim,
                                       const std::string&  cdelim,
                                       const std::string&  precisionreal,
                                       const std::string&  precisionimag,
                                       int                 coffset);
    //!
    //! Gets values as a string
    //! \param fmt Output format
    //!
    virtual std::string GetValues(const OutputFormat* fmt) const;
    //!
    //! Returns true if end of pagination message needs to be printed
    //! \param msg Additional message that needs to be printed
    virtual bool GetPaginationEndMsg(std::string& msg) const;

private:
    mutable bool             _haslargeint;   //!< True if matrix has large ints
    mutable bool             _uppercase;     //!< True if scientific uppercase
    mutable int              _formatinteger; //!< Integer part for format, if applicable
    mutable int              _formatdecimal; //!< Decimal part for format, if applicable
    mutable DisplayFormat    _displayformat; //!< Display format for matrix
    mutable std::streamsize  _precision;     //!< Precision
    mutable std::vector<int> _realwidth;     //!< Widest real value width/column
    mutable std::vector<int> _imagwidth;     //!< Widest imag value width/column
    mutable std::string      _delimiter;     //!< Delimiter between cols

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
    //! \param rows Number of rows
    //! \param cols Number of columns
    //!
	std::string GetPaginationHeader(int rows,
                                    int cols) const;
    //!
    //! Gets output
    //! \param fmt Output format
    //! \param os  Output stream
    //!
    std::string GetOutput(const OutputFormat* fmt,
                          std::ostringstream& os) const;
    //!
	//! Gets matrix data with no pagination - using defaults
	//! \param fmt Format
    //!
	std::string GetOutputNoPagination(const OutputFormat* fmt) const;
    //!
	//! Gets matrix data with no pagination
	//! \param fmt              Format
    //! \param rowdelim         Delimiter for rows
    //! \param startwithnewline Appends a new line to the start of a row
    //! \param coffset          Number of columns which need cdelim prepended
    //!
	std::string GetOutputNoPagination(const OutputFormat* fmt,
                                      const std::string&  rowdelim,
                                      bool                startwithnewline,
                                      int                 coffset) const;
    //!
	//! Gets matrix data with forward pagination
	//! \param fmt Format
    //!
	std::string GetOutputForwardPagination(const OutputFormat* fmt) const;
    //!
	//! Gets matrix data with back pagination
	//! \param fmt Format
    //!
	std::string GetOutputBackPagination(const OutputFormat* fmt) const; 
    //!
    //! Gets output
    //! \param row        Row
    //! \param col        Column
    //! \param realwidth  Width of real field
    //! \param imagwidth  Width of imaginary field
    //! \param isreal     True if this is a real mtx
    //! \param isrealdata True if this is a complex mtx with no imaginary parts
    //! \param output     Output
    //!
    void GetOutput(int          row,
                   int          col,
                   int          realwidth,
                   int          imagwidth,
                   bool         isreal,
                   bool         isrealdata,
                   std::string& output) const;
    //!
    //! Gets number of columns to fit
    //! \param refcol  Column to start from
    //! \param numcols Total number of columns
    //! \param forward True if forward paginating
    //! \param isreal  True if there is only real data
    //!
    int GetNumColumnsToFit(int  refcol,
                           int  numcols,
                           bool forward,
                           bool isreal) const;
    //!
    //! Converts real value to formatted string
    //! \param fmt Format
    //! \param val Value
    //!
    std::string RealToString(double val) const { return RealToString(val, _displayformat); }
    //!
    //! Converts real value to formatted string
    //! \param val Value
    //! \param fmt Format
    //!
    std::string RealToString(double        val,
                             DisplayFormat fmt) const; 
    //!
    //! Utility which returns real value to string
    //! \param val       Value
    //! \param precision Custom precision
    //!
    static std::string RealToString(double             val,
                                    const std::string& precision);
    //!
    //! Gets format for given double
    //! \param val Given value
    //!
    DisplayFormat GetFormat(double val) const;
    //!
    //! Returns format for complex numbers
    //! \param val Given value
    //!
    DisplayFormat GetFormat(const hwComplex& val) const;
    //!
    //! Sets matrix format
    //! \param fmt    Format
    //! \param interp Interpreter
    //!
    void SetFormat(const OutputFormat* fmt, 
                   Interpreter*        interp);
    //!
    //! Resets format
    //!
    void ResetFormat();

    //!
    //! Scan matrix and set the width of different columns
    //! \param interp Interpreter
    //!
    void SetWidth(Interpreter* interp);
    //!
    //! Sets the delimiter
    //! \param delim Delimiter to set, if empty, default delim will be used
    //!
    void SetDelimiter(const std::string& delim);
    //!
    //! Returns true if matrix was paginating
    //!
    bool WasPaginating()  const;
    //!
    //! Returns true if matrix is paginating
    //!
    bool IsPaginating()  const;
    //!
    //! True if header needs to be printed
    //!
    bool IsHeaderPrinted() const;
};

#endif
// End of file:
