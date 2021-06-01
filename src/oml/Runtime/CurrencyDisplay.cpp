/**
* @file CurrencyDisplay.cpp
* @date January 2016
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

#include "CurrencyDisplay.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <climits>

#ifndef OS_WIN
#    include <float.h>
# ifndef DECIMAL_DIG
#     define DECIMAL_DIG 21
# endif
#endif

#include "BuiltInFuncsUtils.h"
#include "Interpreter.h"
#include "OutputFormat.h"
#include "SignalHandlerBase.h"
#include "StructData.h"

#include "hwMatrix.h"
#include "hwMatrixS.h"
#include "math/kernel/GeneralFuncs.h"

int CurrencyDisplay::m_skipFormat   = 10000;
int CurrencyDisplay::m_maxRows      = 0;
int CurrencyDisplay::m_maxCols      = 0;
int CurrencyDisplay::m_linesPrinted = 0;
CurrencyDisplay::PAGINATE CurrencyDisplay::m_paginate = CurrencyDisplay::PAGINATE_ON;

std::ofstream CurrencyDisplay::_outputlog;
std::wstring CurrencyDisplay::_outputlogname = L"outputlog.txt";

//# define CURRENCYDISPLAY_DBG 1  // Uncomment to print debug info
#ifdef CURRENCYDISPLAY_DBG
#    define CURRENCYDISPLAY_PRINT(str,val) { std::cout << str << val << std::endl; }
#else
#    define CURRENCYDISPLAY_PRINT(str,val) 0
#endif
// End defines/includes

//------------------------------------------------------------------------------
// Sets maximum columns for display
//------------------------------------------------------------------------------
void CurrencyDisplay::SetMaxCols(int val)
{
    (val >= 0) ? m_maxCols = val : m_maxCols = 0;
    CURRENCYDISPLAY_PRINT("CurrencyDisplay max cols: ", m_maxCols);
}
//------------------------------------------------------------------------------
// Sets maximum rows for display
//------------------------------------------------------------------------------
void CurrencyDisplay::SetMaxRows(int val)
{
    (val >= 0) ? m_maxRows = val : m_maxRows = 0;
    CURRENCYDISPLAY_PRINT("CurrencyDisplay max rows: ", m_maxRows);
}
//------------------------------------------------------------------------------
// True if given currency can paginate
//------------------------------------------------------------------------------
bool CurrencyDisplay::CanPaginate(const Currency& cur)
{
    if (cur.IsScalar() || cur.IsComplex())
    {
        return false;
    }
    else if (cur.IsCellArray())
    {
        HML_CELLARRAY* cells = cur.CellArray();
        return (cells && cells->M() > 0 && cells->N() > 0);
    }
    else if (cur.IsMatrix())
    {
        return (cur.Matrix() && cur.Matrix()->Size() > 1);
    }
    else if (cur.IsStruct())
    {
        return (cur.Struct() && cur.Struct()->Size() == 1); // Paginates only if size is 1
    }
    else if (cur.IsNDMatrix())
    {
        if (!cur.MatrixN()) return false;
        return true;
    }
    else if (cur.IsString() && CanPaginate(cur.StringVal()))
    {
        return true;
    }
    else if (cur.IsSparse())
    {
        const hwMatrixS* mtx = cur.MatrixS();
        if (!mtx || mtx->Size() == 0 || mtx->IsEmpty())
        {
            return false;
        }
        return true;
    }
    return false;
}
//------------------------------------------------------------------------------
// Constructor - Only currency or derived classes can access constructor
//------------------------------------------------------------------------------
CurrencyDisplay::CurrencyDisplay(const Currency&  cur)
    : m_colBegin      (-1)
    , m_colEnd        (-1)
    , m_rowBegin      (-1)
    , m_rowEnd        (-1)
    , m_mode          (DISPLAYMODE_FORWARD)
    , m_currency      (cur)
    , m_parentDisplay (0)
    , m_initialized   (false)
    , m_signalHandler (NULL)
    , m_indent        (0)
    , m_deleteLine    (false)
{
    // This will be a copy of the currency, so set display explicitly
    m_currency.SetDisplay(this);
}
//------------------------------------------------------------------------------
// Sets indices and data for display mode
//------------------------------------------------------------------------------
void CurrencyDisplay::SetModeData()
{
    switch(m_mode)
    {
        case DISPLAYMODE_DOWN:
            if (CanPaginateColumns())
            {
                SetDownDisplayData();
            }
            else
            {
                SetForwardDisplayData();
            }
            break;

        case DISPLAYMODE_UP:
            if (CanPaginateColumns())
            {
                SetUpDisplayData();
            }
            else
            {
                SetBackDisplayData();
            }
            break;

        case DISPLAYMODE_FORWARD: SetForwardDisplayData(); break;

        case DISPLAYMODE_BACK:    SetBackDisplayData();    break;

        case DISPLAYMODE_RIGHT:   SetRightDisplayData();   break;
        case DISPLAYMODE_LEFT:    SetLeftDisplayData();    break;
        default: break;
    }
}
//------------------------------------------------------------------------------
// True if rows are being processed during pagination
//------------------------------------------------------------------------------
bool CurrencyDisplay::IsPaginatingRows() const
{	
    if (m_currency.IsString())
    {
        std::string str = m_currency.StringVal();
        if (m_rowEnd != -1 && m_rowEnd < static_cast<int>(str.length()))
        {
            return true;
        }
        return false;
    }

    int numrows = 0;
    int numcols = 0;
    GetCurrencySize(numrows, numcols);

    if (numrows <= 0) return false;

    int lastrow = numrows - 1;
	return (m_rowEnd >= 0 && m_rowEnd < lastrow);
}
//------------------------------------------------------------------------------
// True if columns are being processed during pagination
//------------------------------------------------------------------------------
bool CurrencyDisplay::IsPaginatingCols() const
{
    if (m_currency.IsString())
    {
        return false;  // Strings do not paginate columns
    }

    int numrows = 0;
    int numcols = 0;
    GetCurrencySize(numrows, numcols);

    if (numcols <= 0) return false;
    
    int lastcol = numcols - 1;

	return (m_colEnd >= 0 && m_colEnd < lastcol);
}
//------------------------------------------------------------------------------
// Initialize
//------------------------------------------------------------------------------
void CurrencyDisplay::Initialize(const OutputFormat* fmt, 
                                 Interpreter*        interp,
                                 CurrencyDisplay*    parent)
{
    if (parent)
    {
        m_parentDisplay = parent;
    }
    if (m_initialized)
    {
        return;
    }

    m_colBegin = 0;
    m_colEnd   = 0;
    m_rowBegin = 0;
    m_rowEnd   = 0;
    m_mode     = DISPLAYMODE_FORWARD;

    m_signalHandler = interp ? interp->GetSignalHandler() : 0;
    m_initialized   = true;
}
//------------------------------------------------------------------------------
// Gets number of rows that can be fit
//------------------------------------------------------------------------------
int CurrencyDisplay::GetNumRowsToFit() const
{
    if (m_parentDisplay && CurrencyDisplay::IsPaginateOn())
    {
        m_linesPrinted = 0;
       // return m_maxRows;
    }
    int numrows = m_maxRows - m_linesPrinted;
    if (numrows < 1)
    {
        numrows = 1;  // Force print at least one row
    }

    return numrows;
}
//------------------------------------------------------------------------------
// Returns true if end of pagination message needs to be printed
//------------------------------------------------------------------------------
bool CurrencyDisplay::GetPaginationEndMsg(std::string& msg) const
{
    int totalrows = 0;
    int totalcols = 0;
    GetCurrencySize(totalrows, totalcols);

    if (m_rowBegin == 0 && m_rowEnd >= totalrows - 1 &&
        m_colBegin == 0 && m_colEnd >= totalcols - 1)
        return false;  // No need to print end of pagination message

    // Don't print pagination message for empty elements
    if (totalrows == 0 || totalcols == 0) return false;

    return true;
}
//------------------------------------------------------------------------------
// Deletes display
//------------------------------------------------------------------------------
void CurrencyDisplay::DeleteDisplay(CurrencyDisplay* display)
{
    if (display)
    {
        display->m_currency.SetDisplay(nullptr);
        delete display;
        display = nullptr;
    }
}
//------------------------------------------------------------------------------
// Get formatted output string for scalars
//------------------------------------------------------------------------------
std::string CurrencyDisplay::ScalarToString(const OutputFormat* fmt,
                                            double              val,
                                            std::ostringstream& os)
{
    if (IsNegInf_T (val)) return std::string (os.str() + "-Inf");
    if (IsInf_T(val))     return std::string (os.str() + "Inf");
    if (IsNaN_T(val))     return std::string (os.str() + "NaN");

    bool precShort = (fmt && fmt->GetPrecision() == OutputFormat::PRECISION_SHORT);
    bool precLong  = (fmt && fmt->GetPrecision() == OutputFormat::PRECISION_LONG);

    int intpart = fmt ? fmt->GetIntegerPart() : 0;
    int decpart = fmt ? fmt->GetDecimalPart() : 0;

    if (BuiltInFuncsUtils::IsInt(val))
        return IntToString(fmt, val, os);

    char formatString [1024];
    memset(formatString, 0, sizeof(formatString));

    bool scientific      = (fmt && fmt->GetFlags() == std::ios::scientific);
    bool scientificupper = (fmt && fmt->GetFlags() == (std::ios::scientific | std::ios::uppercase));
    bool floatformat     = (!scientific && !scientificupper);

    if (floatformat)
    {
        double aval = fabs(val);
        if (precShort)
        {
            if (aval > 1e-6)
                sprintf(formatString, "%s", "%.5f");
            else                                     // Very small floats
                sprintf(formatString, "%s", "%.5e");
        }
        else if (precLong)
        {
            if (aval > 1e-9)
                sprintf(formatString, "%s", "%.8f");
            else
                sprintf(formatString, "%s", "%.8e");
        }
        else if (intpart > 0 && decpart > 0)
        {
            long long minfloat = (long long)std::pow(10.0, -(decpart+1));
            size_t len = intpart + decpart + 4; // decimal, e, buffer
            bool isScientific = false;
            if (aval > minfloat)
			{
				sprintf(formatString, "%%-%d.%df", intpart, decpart);
            }
            else
            {
                sprintf(formatString, "%%-%d.%de", intpart, decpart);
                isScientific = true;
            }

            if (len >= 1024)
            {
                std::ostringstream os;
                os << std::left;
                if (scientific)
                    os << std::scientific;
                else
                    os << std::fixed;
                os << std::setprecision(static_cast<std::streamsize>(decpart));
                os << val;
                return os.str();
            }
        }
        else if (!(aval > 1e-6))
            sprintf(formatString, "%s", "%.5e");
        else
        {
            if (fmt)
            {
                os.setf(fmt->GetFlags(), std::ios::floatfield);
                os.precision(static_cast<std::streamsize>(fmt->GetPrecision()));
            }
            else
            {
                os.setf(static_cast<std::ios_base::fmtflags>(0), std::ios::floatfield);
                os.precision(static_cast<std::streamsize>(9));
            }
            os << val;
            return std::string (os.str());
        }
    }
    else
    {
        // Scientific notation
        if (precShort)
        {
            scientificupper ? sprintf(formatString, "%s", "%.5E"):
                              sprintf(formatString, "%s", "%.5e");
        }
        else if (precLong)
        {
            scientificupper ? sprintf(formatString, "%s", "%.8E"):
                              sprintf(formatString, "%s", "%.8e");
        }
    }
    char tmp[1024];
    sprintf(tmp, formatString, val);
    return std::string (os.str() +  tmp);
}
//------------------------------------------------------------------------------
// Gets real, imaginary values and sign for complex number
//------------------------------------------------------------------------------
void CurrencyDisplay::GetComplexNumberVals(const hwComplex& val,
                                           double&          realval,
                                           double&          imagval,
                                           std::string&     imagsign)
{	
    realval = val.Real();
    imagval = val.Imag();

    // Add the sign
    imagsign = " + ";
    if (IsInf_T(imagval) || IsNaN_T(imagval)) return;

    if (IsNegInf_T(imagval) || std::signbit(static_cast<long double>(imagval)))
    {
         imagsign = " - ";
         imagval  = fabs(imagval);
    }
}
//------------------------------------------------------------------------------
// Get output string for integers
//------------------------------------------------------------------------------
std::string CurrencyDisplay::IntToString(const OutputFormat* fmt,
                                         double              val,
                                         std::ostringstream& os)
{
    if (IsNegInf_T (val)) return "-Inf";
    if (IsInf_T(val))     return "Inf";
    if (IsNaN_T(val))     return "NaN";

    bool precShort = (fmt && fmt->GetPrecision() == OutputFormat::PRECISION_SHORT);
    bool precLong  = (fmt && fmt->GetPrecision() == OutputFormat::PRECISION_LONG);

    int intpart = fmt ? fmt->GetIntegerPart() : 0;
    int decpart = fmt ? fmt->GetDecimalPart() : 0;

    double          aval        = fabs(val);
    std::streamsize totaldigits = static_cast<std::streamsize>(9);
    bool  customfmt = false;
    if (fmt)
    {
        totaldigits = fmt->GetPrecision();

        // Check for very small floats
        bool scientific = (IsZero(aval) && !(aval < 1e-16));  

        // Check for large integers
        if (aval > LLONG_MAX) 
            scientific = true;

        if (precShort)
        {
            long long maxint = (long long)1e+12;       // max possible int
            if (!(aval < maxint))
                scientific = true;

            if (scientific)
                totaldigits = static_cast<std::streamsize>(5);
            else
                totaldigits = static_cast<std::streamsize>(12);
        }

        else if (precLong)
        {
            long long maxdigits = (long long)1e+15;       // max digits
            if (!(aval < maxdigits))
                scientific = true;
            if (scientific)
                totaldigits = static_cast<std::streamsize>(15);
            else
                totaldigits = static_cast<std::streamsize>(24);
        }
        else if (intpart > 0 && decpart > 0)
        {
            customfmt = true;
            totaldigits = static_cast<std::streamsize>(intpart + decpart);
            long long maxdigits = static_cast<long long>(std::pow(10.0, static_cast<double>(totaldigits)));
            if (!(aval < maxdigits))
                scientific = true;
        }

        if (!scientific ||
            (scientific && fmt->GetFlags() == (std::ios::scientific | std::ios::uppercase)))
            os.setf(fmt->GetFlags(), std::ios::floatfield);        
        else if (scientific)
            os.setf(std::ios::scientific, std::ios::floatfield); 
        
    }
    if (customfmt)
        os.precision(decpart);
    else
        os.precision(totaldigits);
    os << val;
    return std::string(os.str());
}
//------------------------------------------------------------------------------
// Returns true if display was paginating
//------------------------------------------------------------------------------
bool CurrencyDisplay::WasPaginating() const
{
    return (!(m_colBegin == -1 && m_colEnd == -1 && m_rowBegin == -1 && m_rowEnd == -1));
}
//------------------------------------------------------------------------------
// Get formatted output string for scalars
//------------------------------------------------------------------------------
CurrencyDisplay::DisplayFormat CurrencyDisplay::GetFormatInfo(
    const OutputFormat* fmtIn,
    double              val,
    std::string&        fmtstr,
    bool&               hasLargeInt)
{
    if (IsNegInf_T(val) || IsInf_T(val) || IsNaN_T(val)) return DisplayFormatInt;

    OutputFormat* fmt = fmtIn ? new OutputFormat(*fmtIn) : nullptr;

    bool scientific   = (fmt && fmt->GetFlags() == std::ios::scientific);
    bool scientificup = (fmt && fmt->GetFlags() == (std::ios::scientific | 
                                                    std::ios::uppercase)); 
    bool precshort = (fmt && fmt->GetPrecision() == OutputFormat::PRECISION_SHORT);
    bool preclong  = (fmt && fmt->GetPrecision() == OutputFormat::PRECISION_LONG);
    bool precscalar = (!precshort && !preclong);

    // Define limits to determine format
    int intpart = fmt ? fmt->GetIntegerPart() : 0;
    int decpart = fmt ? fmt->GetDecimalPart() : 0;

    std::string     partialFmt;
    std::streamsize precision  = fmt ? fmt->GetPrecision() : 
                                 OutputFormat::PRECISION_SCALAR;
    int         totaldigits = 12;
    double      aval        = fabs(val);
    long double maxfloat    = 1e+6;
    long double minfloat    = 1e-6;

    long long maxdigits = static_cast<long long>(1e+15);
    long long maxint    = precscalar ? maxdigits : static_cast<long long>(1e+12);

    if (intpart > 0 && decpart > 0)
    {
        int len = intpart + decpart + 4; // decimal, e, buffer
        if (len >= 1024) // Very large format string
        {
            std::ostringstream os;
            os << std::left << std::fixed;
            os << std::setprecision(static_cast<std::streamsize>(decpart));
            os << val;
            std::string val (os.str());
            bool setdefault = true;
            if (!val.empty())
            {
                size_t pos = val.find(".");
                if (pos != std::string::npos)
                {
                    setdefault = false;
                    std::string intstr = val.substr(0, pos);
                    std::string decstr = val.substr(pos + 1);
                    
                    intpart = intstr.empty() ? 9 : static_cast<int>(intstr.size());
                    decpart = decstr.empty() ? 8 : static_cast<int>(decstr.size());
                }
            }
            if (setdefault)
            {
                intpart = std::min(15, intpart);
                decpart = std::min(9,  decpart);
            }
            fmt->SetCustomFormat(intpart, decpart);  
        }
        totaldigits = intpart + decpart;
        maxfloat    = std::pow(10.0, intpart);
        minfloat    = std::pow(10.0, -(decpart + 1));
        partialFmt  = "%" + std::to_string(static_cast<long long>(intpart)) +
                        "." + std::to_string(static_cast<long long>(decpart));
    }
    else if (preclong)
    {
        totaldigits = 24;
        minfloat    = 1e-9;
    }

    bool stdformat = false;
    if (partialFmt.empty())
    {
        if (fmt && !precscalar)
        {
            std::ostringstream os;
            os << "%." << static_cast<std::streamsize>(fmt->GetPrecision());
            partialFmt = os.str();
        }
        else
        {
            stdformat = true;
            partialFmt = "%.8";
        }
    }        

    delete fmt;
    fmt = nullptr;

    DisplayFormat thisformat = (scientific || scientificup) ? 
                               DisplayFormatScientific : DisplayFormatInt;

    if (thisformat == DisplayFormatInt && 
        BuiltInFuncsUtils::IsInt(val)) // Integer, compare by values
    {
        if ((IsZero(aval) && !(aval < 1e-16))     || // Very small floats
            (precshort    && !(aval < maxint))    || // Very large ints
            (precscalar   && !(aval < maxdigits)) || 
            (preclong     && !(aval < maxdigits)) ||
             aval > LLONG_MAX)
             thisformat = DisplayFormatScientific;
        else
        {
            std::string tmp (GetFormattedValue(val, DisplayFormatInt, fmt, ""));
            if (!tmp.empty())
            {
                if (tmp.find("e") != std::string::npos ||
                    tmp.find("E") != std::string::npos ||
                    tmp.size()    > totaldigits) 
                    thisformat = DisplayFormatScientific;
            }

            // Additional check if we have large integers
            if (!hasLargeInt && !(aval < maxfloat) && (precscalar || precshort))
                hasLargeInt = true;

            if (!hasLargeInt && thisformat == DisplayFormatInt) 
            {
                return DisplayFormatInt;
            }
        }
    }

    if (hasLargeInt || !(aval > minfloat)) // Very large ints/very small floats
        thisformat = DisplayFormatScientific;

    if (thisformat == DisplayFormatScientific)
    {
        fmtstr = partialFmt;
        if (scientificup)
            fmtstr += "E";
        else
            fmtstr += "e";

        return DisplayFormatScientific;
    }

    if (precshort)
    {
        if (aval > 1e-6) 
        {
            fmtstr = "%.5lf";
            return DisplayFormatFloat;  // sprintf(formatString, "%s", "%.5f");
        }
        // Very small floats
        fmtstr = "%.5e";
        return DisplayFormatScientific; // sprintf(formatString, "%s", "%.5e");
    }
    else if (preclong)
    {
        if (aval > 1e-9)
        {
            fmtstr = "%.8lf";
            return DisplayFormatFloat; //sprintf(formatString, "%s", "%.8f");
        }
        fmtstr = "%.8e";
        return DisplayFormatScientific; // sprintf(formatString, "%s", "%.8e");
    }
        
    if (intpart > 0 && decpart > 0)
    {
        if (aval > minfloat)
        {
            fmtstr = partialFmt + "lf";
            return DisplayFormatFloat;
        }
        fmtstr = partialFmt + "e";
        return DisplayFormatScientific; //
    }
    if (!(aval > 1e-6))
    {
        fmtstr = "%.9e"; //fmtstr = "%.5e";
        return DisplayFormatScientific; // sprintf(formatString, "%s", "%.5e");
    }

    std::string tmp (GetFormattedValue(val, DisplayFormatFloat, fmt, ""));
    if (!tmp.empty())
    {
        if (tmp.find("e") != std::string::npos ||
            tmp.find("E") != std::string::npos ||
            tmp.size()    > totaldigits) 
            thisformat = DisplayFormatScientific;
        else if (tmp.find(".") == std::string::npos)
        {
            fmtstr = "";
            return DisplayFormatInt;
        }
        if (!stdformat)
            fmtstr = partialFmt + "lf";

        return DisplayFormatFloat;
    }

    // Scientific notation
    if (precshort)
        fmtstr = scientificup ? "%.5E": "%.5e";
    else if (preclong)
        fmtstr = scientificup ? "%.8E": "%.8e";
    else
    {
        fmtstr = partialFmt;
        if (scientificup)
            fmtstr += "E";
        else
            fmtstr += "e";
    }

    return DisplayFormatScientific; 
}
//------------------------------------------------------------------------------
// Get output string for complex numbers
//------------------------------------------------------------------------------
std::string CurrencyDisplay::ComplexToString(const OutputFormat* fmt,
                                             const hwComplex&    val,
                                             std::ostringstream& os)
{
    // Split the complex number into real, imaginary and sign
    double rval = 0.0;
    double ival = 0.0;
    std::string isign;
    GetComplexNumberVals(val, rval, ival, isign);

    bool haslargeint = false;

    // Get the display format, format string for real value
    std::string   rfmtstr;
    DisplayFormat rfmt = GetFormatInfo(fmt, rval, rfmtstr, haslargeint);

    DisplayFormat commonfmt = rfmt;
    std::string   commonfmtstr (rfmtstr);
    if (rfmt != DisplayFormatScientific)
    {
        // Get the display format, format string for imaginary value
        std::string   ifmtstr;    
        DisplayFormat ifmt = GetFormatInfo(fmt, ival, ifmtstr, haslargeint);

        // The common format will be the max of the two formats
        commonfmt    = (ifmt > rfmt) ? ifmt    : rfmt;
        commonfmtstr = (ifmt > rfmt) ? ifmtstr : rfmtstr;
    }

    if (commonfmtstr.empty() && commonfmt == DisplayFormatFloat)
    {
        std::ostringstream tmp1;
        std::string rstr = ScalarToString(fmt, rval, tmp1);
        int precR = 0;
        if (!rstr.empty())
        {
            size_t pos = rstr.find(".");
            if (pos != std::string::npos)
                precR = int(rstr.length() - pos -1);
        }
        std::ostringstream tmp2;
        std::string istr = ScalarToString(fmt, ival, tmp2);
        int precI = 0;
        if (!istr.empty())
        {
            size_t pos = istr.find(".");
            if (pos != std::string::npos)
                precI = int(istr.length() - (int)pos-1);
        }
        int prec = std::max(precR, precI);
        if (prec > 8)
            prec = 8;
        commonfmtstr = "%0." + std::to_string(static_cast<long long>(prec)) + "f";
    }
 
    // Assemble the output
    std::string out (os.str());

    std::string tmp = GetFormattedValue(rval, commonfmt, fmt, commonfmtstr);
    if (!tmp.empty())
    {
        size_t pos = tmp.find_first_not_of(" ");
        if (pos != std::string::npos)
        {
            out += tmp.substr(pos);
        }
        else
        {
            out += tmp;
        }
    }
    out += isign; 
    tmp = GetFormattedValue(ival, commonfmt, fmt, commonfmtstr) + "i";
    if (!tmp.empty())
    {
        size_t pos = tmp.find_first_not_of(" ");
        if (pos != std::string::npos)
        {
            out += tmp.substr(pos);
        }
        else
        {
            out += tmp;
        }
    }

    return out;
}
//------------------------------------------------------------------------------
// Get formatted output string for scalars
//------------------------------------------------------------------------------
std::string CurrencyDisplay::GetFormattedValue(double val,
        DisplayFormat commonfmt,
        const OutputFormat* fmt,
        const std::string&  fmtstr)
{
    std::ostringstream os;
    if (commonfmt == DisplayFormatInt)
        return CurrencyDisplay::IntToString(fmt, val, os);

    if (fmtstr.empty())
    {
        std::ostringstream os;
        if (fmt)
        {
            os.setf(fmt->GetFlags(), std::ios::floatfield);
            os.precision(static_cast<std::streamsize>(fmt->GetPrecision()));
        }
        else
        {
            os.setf(static_cast<std::ios_base::fmtflags>(0), std::ios::floatfield);
            os.precision(static_cast<std::streamsize>(9));
        }
        os << val;
        return std::string (os.str());
    }
        
    char tmp[1024];
    sprintf(tmp, fmtstr.c_str(), val);
    return std::string (tmp);
}
//------------------------------------------------------------------------------
// True if there is a valid visible display size and pagination is enabled
//------------------------------------------------------------------------------
bool CurrencyDisplay::IsValidDisplaySize()
{
    if (m_paginate == CurrencyDisplay::PAGINATE_OFF)
    {
        return false;
    }
    return (m_maxRows > 0 && m_maxCols > 0);
}
//------------------------------------------------------------------------------
// True if there are rows that can be printed
//------------------------------------------------------------------------------
bool CurrencyDisplay::CanPrintRows() const
{
    if (GetNumRowsToFit() - m_linesPrinted >= 1)
    {
        return true;
    }
    return false;
}
//------------------------------------------------------------------------------
// True if given string can paginate
//------------------------------------------------------------------------------
bool CurrencyDisplay::CanPaginate(const std::string& str)
{
    if (!str.empty())
    {
        size_t pos = str.find("\n");
        if (pos != std::string::npos && pos < str.length() - 1)
        {
            // Paginate only if new line found is not the last character
            return true;
        }
    }
    return false;
}
//------------------------------------------------------------------------------
// Strip trailing new line
//------------------------------------------------------------------------------
void CurrencyDisplay::StripEndline(std::string& str) const
{
    if (!str.empty() && str[str.size() - 1] == '\n')
    {
        str.pop_back();
    }
    if (!str.empty() && str[str.size() - 1] == '\r')
    {
        str.pop_back();
    }
}
//------------------------------------------------------------------------------
// Gets indent string
//------------------------------------------------------------------------------
std::string CurrencyDisplay::GetIndentString(int val) const
{
    if (val > 0)
    {
        return std::string(2 * val, ' ');
    }
    return "";
}
//------------------------------------------------------------------------------
// True if header needs to be printed
//------------------------------------------------------------------------------
bool CurrencyDisplay::IsHeaderPrinted() const
{
    if (m_parentDisplay && m_parentDisplay->IsNDMatrixDisplay() &&
        WasPaginating())
    {
        return false;
    }
    return true;
}
//------------------------------------------------------------------------------
// Return estimated length of formatted string for given value
//------------------------------------------------------------------------------
size_t CurrencyDisplay::GetFormattedStringLength(const char* fmt, double val)
{
    size_t len = 0;
    if (!fmt)
    {
        return len;
    }
    try
    {
        len = snprintf(nullptr, 0, fmt, val) + 1; // Extra space for '\0'
        if (len > 0)
        {
            return (len - 1);  // Don't include '\0'
        }
    }
    catch (...)
    {
        return 0;
    }
    return len;
}
//------------------------------------------------------------------------------
// Return formatted string for value
//------------------------------------------------------------------------------
std::string CurrencyDisplay::GetFormattedString(const char* fmt, double val)
{
    size_t len = GetFormattedStringLength(fmt, val);
    if (len > 0)
    {
        try
        {
            std::unique_ptr<char[]> buf(new char[len + 1]);  // Include '\0'
            snprintf(buf.get(), len + 1, fmt, val);
            return std::string(buf.get(), buf.get() + len); 
        }
        catch (...)
        {
            return "";
        }
    }
    return "";
}
//------------------------------------------------------------------------------
// Sets mode data for interactive pagination
//------------------------------------------------------------------------------
void CurrencyDisplay::SetPaginateOnModeData()
{
    SetMode(DISPLAYMODE_FORWARD);
    SetModeData();
}
//------------------------------------------------------------------------------
// Sets size, after which printing will be done skipping format of each element
// Value > 1 will result in formatting being checked (traditional printing)
// Value < 1 will result in format of first element applied to matrix (fast printing)
//------------------------------------------------------------------------------
void CurrencyDisplay::SetSkipFormat(int val)
{
    if (val < 1 || IsNegInf_T(val) || IsInf_T(val) || IsNaN_T(val))
    {
        m_skipFormat = -1;
    }
    else
    {
        m_skipFormat = val;
    }
}
//------------------------------------------------------------------------------
// Sets name of outputlog
//------------------------------------------------------------------------------
void CurrencyDisplay::SetOutputLogName(const std::wstring& name)
{
    if (name.empty())
    {
        _outputlogname = L"omloutputlog.txt";
    }
    else
    {
        size_t pos = name.find_last_of(L"/\\");
        if (pos == std::wstring::npos)
        {
            _outputlogname = name;
        }
        _outputlogname = name.substr(pos + 1);
    }
}
//------------------------------------------------------------------------------
// Utility to convert double to string without precision loss
//------------------------------------------------------------------------------
std::string CurrencyDisplay::NonFormattedDoubleToString(double val)
{
    if (IsNegInf_T(val))
    {
        return "-Inf";
    }
    else if (IsInf_T(val))
    {
        return "Inf";
    }
    else if (IsNaN_T(val))
    {
        return "NaN";
    }
    else
    {
        char* tmp = new char[128];
#ifdef OS_WIN
        sprintf(tmp, "%.*g", DBL_DECIMAL_DIG, val);
#else
        sprintf(tmp, "%.*g", DECIMAL_DIG, val);
#endif

        std::string out(tmp);

        delete[] tmp;
        tmp = nullptr;

        return out;
    }
    return "";
}
//------------------------------------------------------------------------------
// Utility to convert double to string without precision loss
//------------------------------------------------------------------------------
std::string CurrencyDisplay::NonFormattedComplexToString(const hwComplex& val)
{
    std::string out (NonFormattedDoubleToString(val.Real()));

    double imagval = val.Imag();

    if (IsNegInf_T(imagval) || std::signbit(static_cast<long double>(imagval)))
    {
        out += "-" + NonFormattedDoubleToString(fabs(imagval)) + "i";
    }
    else
    {
        out += "+" + NonFormattedDoubleToString(imagval) + "i";
    }
    return out;
}
