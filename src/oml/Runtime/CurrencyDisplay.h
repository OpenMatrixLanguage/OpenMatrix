/**
* @file CurrencyDisplay.h
* @date January 2016
* Copyright (C) 2016-2020 Altair Engineering, Inc.  
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

#ifndef __CURRENCYDISPLAY_H__
#define __CURRENCYDISPLAY_H__

// Begin defines/includes
#include "Hml2Dll.h"

#include "Currency.h"

class Interpreter;
class OutputFormat;
class SignalHandlerBase;
// End defines/includes

//------------------------------------------------------------------------------
//! \class CurrencyDisplay
//! \brief Displays currency
//------------------------------------------------------------------------------
class HML2DLL_DECLS CurrencyDisplay
{
public:
    friend class Currency; //!< Only currency is allowed to construct

    //! Destructor 
    virtual ~CurrencyDisplay() {}

    //!
    //! \enum PAGINATE
    //!
    enum PAGINATE
    {
        PAGINATE_OFF = 0,              //!< No pagination
        PAGINATE_ON,                   //!< Paginate
        PAGINATE_INTERACTIVE           //!< Interactive pagination
    };

    //!
    //! Sets maximum columns for display
    //! \param val Max columns, should be greater than or equal to 0
    //!
    static void SetMaxCols(int val);
    //!
    //! Sets maximum rows for display
    //! \param val Max rows, should be greater than or equal to 0
    //!
    static void SetMaxRows(int val);
    
    //!
    //! True if there is a valid visible display size
    //!
    static bool IsValidDisplaySize();
    //!
    //! Clears count of lines displayed
    //!
    static void ClearLineCount() { m_linesPrinted = 0; }
    //!
    //! Deletes display
    //! \param display Given display
    //!
    static void DeleteDisplay(CurrencyDisplay* display);
    //!
    //! Gets pagination mode
    //!
    static CurrencyDisplay::PAGINATE GetPaginate() { return m_paginate; }
    //!
    //! Sets pagination mode
    //! \param Pagination mode
    //!
    static void SetPaginate(CurrencyDisplay::PAGINATE val) { m_paginate = val; }
    //!
    //! Returns true if pagination is off
    //!
    static bool IsPaginateOff() { return (m_paginate == CurrencyDisplay::PAGINATE_OFF); }
    //!
    //! Returns true if pagination is on
    //!
    static bool IsPaginateOn() { return (m_paginate == CurrencyDisplay::PAGINATE_ON); }
    //!
    //! Returns true if pagination is interactive
    //!
    static bool IsPaginateInteractive() { return (m_paginate == CurrencyDisplay::PAGINATE_INTERACTIVE); }
    //!
    //! True if given currency can paginate
    //! \param cur Given currency
    //!
    static bool CanPaginate(const Currency& cur);
    //!
    //! True if rows are being processed during pagination
    //!
    bool IsPaginatingRows() const;
    //!
    //! True if columns are being processed during pagination
    //!
    bool IsPaginatingCols() const;
    //!
    //! True if currency can paginate columns
    //!
    virtual bool CanPaginateColumns() const { return false; }
    //!
    //! Returns true if end of pagination message needs to be printed
    //! \param msg Additional message that needs to be printed
    //!
    virtual bool GetPaginationEndMsg( std::string& msg) const;

    //!
    //! Gets output - needs to be handled in derived classes
    //! \param fmt Output format
    //! \param os  Output stream
    //!
    virtual std::string GetOutput( const OutputFormat* fmt,
                                   std::ostringstream& os) const = 0;

    //!
    //! \enum DISPLAYMODE
    //!
    enum DISPLAYMODE
    {
        DISPLAYMODE_QUIT,              //!< Quit display
        DISPLAYMODE_FORWARD,           //!< Display forward
        DISPLAYMODE_BACK,              //!< Display back
        DISPLAYMODE_UP,                //!< Display up
        DISPLAYMODE_DOWN,              //!< Display down
        DISPLAYMODE_LEFT,              //!< Display left
        DISPLAYMODE_RIGHT,             //!< Display right
        DISPLAYMODE_EXIT               //!< Exit all (nested) pagination
    };
    //!
    //! Gets mode
    //!
    DISPLAYMODE GetMode() const { return m_mode; }
    //!
    //! Sets mode
    //! \param mode Mode to set
    //!
    void SetMode(DISPLAYMODE mode) { m_mode = mode; }
    //!
    //! Sets indices and data for mode
    //!
    void SetModeData();

    //!
    //! Gets currency
    //!
    const Currency& GetCurrency() const { return m_currency; }
    //!
    //! Gets parent display
    //!
    CurrencyDisplay* GetParentDisplay() const { return m_parentDisplay; }
    //!
    //! Sets parent display
    //!
    void SetParentDisplay(CurrencyDisplay* parent) { m_parentDisplay = parent; }

    //!
    //! Initialize
    //! \param fmt    Output format   
    //! \param interp Interpreter
    //! \param parent Parent display
    //!
    virtual void Initialize(const OutputFormat* fmt,
                            Interpreter*        interp,
                           CurrencyDisplay*    parent = 0);
    //!
    //! Gets number of rows and cols in given currency
    //! \param cur   Given cell array/matrix
    //! \param rows Number of rows
    //! \param cols Number of columns
    //!
    virtual void GetCurrencySize(int& rows, 
                                 int& cols) const = 0;
    //!
    //! Returns true if display was paginating
    //!
    bool WasPaginating() const;

    //!
    //! Gets real, imaginary values and sign for complex number
    //! \param val       Given complex number
    //! \param realval  Real part of complex number
    //! \param imagval  Imaginary part of complex number
    //! \param imagsign Sign of complex number
    //!
    static void GetComplexNumberVals(const hwComplex& val,
                                     double&          realval,
                                     double&          imagval,
                                     std::string&     imagsign);
    //!
    //! True if there are rows that can be printed
    //!
    bool CanPrintRows() const;
    //!
    //! True if given string can paginate
    //! \param str Given string
    //!
    static bool CanPaginate(const std::string& str);

    //!
    //! Returns true if this is a chained display, which needs to be deleted
    //! along with given currency display
    //! \param display Given display
    //!
    virtual bool IsChainedDisplay(CurrencyDisplay* display) const { return false; }
    //!
    //! Gets chained display
    //!
    virtual long GetChainedDisplay() const { return -1; }
    //!
    //! Sets chained display
    //! \param id Chained display id
    //!
    virtual void SetChainedDisplay(long id) {}

    //!
    //! Gets indentation factor
    //!
    int GetIndent() const { return m_indent; }
    //!
    //! Sets indent value
    //! \param val Indent value
    //!
    void SetIndent(int val) { m_indent = val; }
    //!
    //! Returns true if parent is ND matrix
    //!
    virtual bool IsNDMatrixDisplay() const { return false; }

    //!
    //! Returns true if first line needs to be deleted in client before printing
    //!
    bool GetDeleteLine() const { return m_deleteLine; }
    //!
    //! Sets to true if an extra line needs to be deleted in line before print
    //! \param val True if first line in client needs to be deleted before printing
    //!
    void SetDeleteLine(bool val) { m_deleteLine = val;}

    //!
    //! Return estimated length of formatted string for value
    //! \param Format specification
    //! \param Value
    //!
    static size_t GetFormattedStringLength(const char*, 
                                           double);
    //!
    //! Return formatted string for value
    //! \param Format specification
    //! \param Value
    //!
    static std::string GetFormattedString(const char*,
                                          double);
    //!
    //! Sets mode data for interactive pagination
    //!
    void SetPaginateOnModeData();

protected:  
    //!
    //! \enum DisplayFormat
    //!
    enum DisplayFormat
    {
        DisplayFormatInt,              //!<  Integer format
        DisplayFormatFloat,            //!< Float format
        DisplayFormatScientific,       //!< Scientific format
    };

    static int m_maxCols;              //!< Max columns (chars) for display        
    static int m_maxRows;              //!< Max rows (lines) for display
    static int m_linesPrinted;         //!< of lines displayed
    static PAGINATE m_paginate;        //!< Pagination mode
    
    mutable int m_colBegin;            //!< 0-based index of first column displayed
    mutable int m_colEnd;              //!< 0-based index of last  column displayed
    mutable int m_rowBegin;            //!< 0-based index of first row    displayed
    mutable int m_rowEnd;              //!< 0-based index of last  row    displayed
    mutable int m_indent;              //!< Indentation factor

    mutable bool       m_deleteLine;     //!< True if an extra line needs to be deleted before printing
    bool               m_initialized;    //!< True if initialized
    DISPLAYMODE        m_mode;           //!< Display mode
    CurrencyDisplay*   m_parentDisplay;  //!< Parent display - for nested currencies
    SignalHandlerBase* m_signalHandler;  //!< Signal handler - for client communication
    Currency           m_currency;       //!< Currency associated with this display

    // Note: Currency is a copy and not a pointer because when results are 
    // displayed, the original currency could have already been destroyed
    
    //!
    //! Constructor - Only currency or derived classes can access constructor
    //! \param cur Currency associated with this display
    //!
    CurrencyDisplay(const Currency& cur);

    //!
    //! True if header needs to be printed
    //!
    bool IsHeaderPrinted() const;

    //!
    //! Sets data for forward pagination
    //!
    virtual void SetForwardDisplayData() {}
    //!
    //! Sets data for back pagination
    //!
    virtual void SetBackDisplayData() {}
    //!
    //! Sets data for right pagination
    //!
    virtual void SetRightDisplayData() {}
    //!
    //! Sets data for left pagination
    //!
    virtual void SetLeftDisplayData() {}
    //!
    //! Sets data for up pagination
    //!
    virtual void SetUpDisplayData() {}
    //!
    //! Sets data for down pagination
    //!
    virtual void SetDownDisplayData() {}

    //!
    //! Gets number of rows that can be fit
    //!
    virtual int GetNumRowsToFit() const;

    //!
    //! Get output string for complex numbers - experimental mode
    //! \param fmt Output format
    //! \param val Given number
    //! \param os  Output stream
    //!
    static std::string ComplexToString(const OutputFormat* fmt,
                                       const hwComplex&    val,
                                       std::ostringstream& os);
    //!
    //! Gets values as a string - needs to be handled in derived classes
    //! \param fmt Output format
    //!
    virtual std::string GetValues(const OutputFormat* fmt) const = 0;
    //!
    //! Strip trailing new line
    //! \param str Given string
    //!
    void StripEndline(std::string& str) const;
    //!
    //! Gets indent string
    //! \param indent value to use
    //!
    std::string GetIndentString(int val) const;
    
private:
    CurrencyDisplay();                                        //!< Stubbed out 
    CurrencyDisplay(const CurrencyDisplay& src) ;             //!< Stubbed out 
    CurrencyDisplay& operator=(const CurrencyDisplay& src);  //!< Stubbed out

    //!
    //! Get output string for integers
    //! \param fmt Output format
    //! \param val Given value
    //! \param os  Output stream
    //!
    static std::string IntToString(const OutputFormat* fmt,
                                   double              val,
                                   std::ostringstream& os);
    //!
    //! Get output string for scalars
    //! \param fmt Output format
    //! \param val Given real number
    //! \param os  Output stream
    //!
    static std::string ScalarToString(const OutputFormat* fmt,
                                      double              val,
                                      std::ostringstream& os);
    //!
    //! Get format type and string for output
    //! \param  fmt    Output format
    //! \param  val    Given real number
    //! \param  fmtstr Format string for output
    //!
    static DisplayFormat GetFormatInfo(const OutputFormat* fmt,
                                       double              val,
                                       std::string&        fmtstr,
                                       bool&               hasLargeInt);
    //!
    //! Get formatted output string for scalars
    //! \param fmt Output format
    //! \param val Given real number
    //! \param os  Output stream
    //!
    static std::string GetFormattedValue(double              val,
                                         DisplayFormat       type,
                                         const OutputFormat* fmt,
                                         const std::string&  fmtstr);
};
#endif
// End of file:
