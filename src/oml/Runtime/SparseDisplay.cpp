/**
* @file SparseDisplay.cpp
* @date June, 2019
* Copyright (C) 2019 Altair Engineering, Inc.
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
#include "SparseDisplay.h"

#include <cassert>
#include <iomanip>
#include <sstream>

#include "CurrencyDisplay.h"
#include "Interpreter.h"

template <typename T1, typename T2> class hwTMatrixS;
typedef hwTMatrixS<double, hwTComplex<double> > hwMatrixS;

#include "hwMatrixS.h"
// End defines/includes

//------------------------------------------------------------------------------
// Constructor - Only currency is allowed to construct
//------------------------------------------------------------------------------
SparseDisplay::SparseDisplay(const Currency& cur)
    : CurrencyDisplay(cur)
    , _realwidth   (0)
    , _imagwidth   (0)
{
}
//------------------------------------------------------------------------------
// Gets number of rows and cols in currency
//------------------------------------------------------------------------------
void SparseDisplay::GetCurrencySize(int& rows, int& cols) const
{
    const hwMatrixS* mtx = m_currency.MatrixS();
    if (mtx)  // Need to only display non-zero elements
    {
        rows = mtx->NNZ();
        cols = 1;
    }
}
//------------------------------------------------------------------------------
// Initializes
//------------------------------------------------------------------------------
void SparseDisplay::Initialize(const OutputFormat* fmt,
                               Interpreter*        interp,
                               CurrencyDisplay*    parent)
{
    if (!m_initialized)
    {
        const hwMatrixS* mtx = m_currency.MatrixS();
        if (!mtx || mtx->Size() == 0 || mtx->IsEmpty())
        {
            return;
        }
        _formatvars.Reset();          // Reset format and widths
        SetFormat(fmt, interp); // Sets the format
        SetWidth(interp);       // Scan matrix and set the width of columns
    }

    CurrencyDisplay::Initialize(fmt, interp, parent);
}
//------------------------------------------------------------------------------
// Gets output - called from Currency::GetOutputString
//------------------------------------------------------------------------------
std::string SparseDisplay::GetOutput(const OutputFormat* fmt,
                                     std::ostringstream& os) const
{
    if (m_mode == DISPLAYMODE_QUIT || m_mode == DISPLAYMODE_EXIT)
    {
        _formatvars.Reset();
        return "";   // Do nothing
    }

    // Print header
    os << GetIndentString(m_indent) << m_currency.GetTypeString();
    std::string header(os.str());

    const hwMatrixS* mtx = m_currency.MatrixS();
    if (!mtx || mtx->Size() == 0 || mtx->IsEmpty())
    {
        return header;
    }
    if (!IsHeaderPrinted())
    {
        header = "";
    }

    std::string output;
    if (!IsValidDisplaySize()) 
    {
        std::string output(GetOutputNoPagination(fmt));
        return (header + output);
    }

    if (!IsHeaderPrinted())
    {
        header = "";
    }
    if (!m_initialized)
    {
        _formatvars.Reset();
        const_cast<SparseDisplay*>(this)->Initialize(fmt, nullptr);
    }
    if (!m_initialized)
    {
        return "";
    }

    if (m_mode == DISPLAYMODE_FORWARD || m_mode == DISPLAYMODE_DOWN)
    {
        output = GetOutputForwardPagination(fmt);
    }
    else if (m_mode == DISPLAYMODE_BACK || m_mode == DISPLAYMODE_UP)
    {
        output = GetOutputBackPagination(fmt);
    }

    if (!output.empty())
    {
        return header + output;
    }
    else if (mtx->NNZ() == 0)
    {
        return header;
    }
    return "";
}
//------------------------------------------------------------------------------
// Gets data with no pagination
//------------------------------------------------------------------------------
std::string SparseDisplay::GetOutputNoPagination(const OutputFormat* fmt) const
{

    const hwMatrixS* mtx = m_currency.MatrixS();
    if (!mtx || mtx->Size() == 0 || mtx->IsEmpty())
    {
        return "";
    }

    if (!m_initialized)
    {
        _formatvars.Reset();
        const_cast<SparseDisplay*>(this)->Initialize(fmt, nullptr);
    }
    if (!m_initialized)
    {
        return "";
    }

    int nnz         = mtx->NNZ();
    bool isreal     = mtx->IsReal();
    bool isrealdata = mtx->IsRealData();
    std::string output;

    std::string myindent = GetIndentString(m_indent);

    for (int i = 0; i < nnz; ++i)
    {
        output += "\n" + myindent;
        GetOutput(i, isreal, isrealdata, output);
    }
    _formatvars.Reset();
    return output;
}
//------------------------------------------------------------------------------
// Gets output
//------------------------------------------------------------------------------
void SparseDisplay::GetOutput(int          index, 
                              bool         isreal,
                              bool         isrealdata,
                              std::string& output) const
{
    const hwMatrixS* mtx = m_currency.MatrixS();

    assert(mtx);
    assert(index >= 0 && index < mtx->NNZ());

    int row = 0;
    int col = 0;

    double      realval = 0.0;
    double      imagval = 0.0;
    std::string imagsign;
    if (isreal)
    {
        mtx->NZinfo(index, row, col, realval);
    }
    else
    {
        hwComplex val(0, 0);
        mtx->NZinfo(index, row, col, val);
        GetComplexNumberVals(val, realval, imagval, imagsign);
    }

    std::ostringstream os;
    os << "[" << row + 1 << "," << col + 1 << "] ";
    os << std::right << std::setfill(' ') << std::setw(_realwidth);
    os << _formatvars.RealToString(realval);

    output += os.str();
    if (!isrealdata)
    {
        output += imagsign;
        std::ostringstream ios;
        ios << std::right << std::setfill(' ') << std::setw(_imagwidth);
        ios << (_formatvars.RealToString(imagval) + "i");

        output += ios.str();
    }
}
//------------------------------------------------------------------------------
// Sets matrix format
//------------------------------------------------------------------------------
void SparseDisplay::SetFormat(const OutputFormat* fmt, Interpreter* interp)
{
    _formatvars.Initialize(fmt);

    const hwMatrixS* mtx = m_currency.MatrixS();
    if (!mtx || mtx->Size() == 0 || mtx->IsEmpty())
    {
        return;
    }

    int nnz = mtx->NNZ();

    // Scan matrix and set format
    bool isreal = mtx->IsReal();

    for (int i = 0; i < nnz; ++i)
    {
        if (_formatvars.GetFormatType() == DisplayFormatVars::TypeScientific)
        {
            break;  // This is the max format that can be applied
        }
        if (interp && interp->IsInterrupt())
        {
            _formatvars.Reset();
            return;
        }
        int row = 0;
        int col = 0;
        DisplayFormatVars::Type thisformat = DisplayFormatVars::TypeInt;
        if (isreal)
        {
            double val = 0;
            mtx->NZinfo(i, row, col, val);
            thisformat = _formatvars.GetFormatType(val);
        }
        else
        {
            hwComplex val(0, 0);
            mtx->NZinfo(i, row, col, val);
            thisformat = _formatvars.GetFormatType(val);
        }
        if (_formatvars.GetFormatType() < thisformat)
        {
            _formatvars.SetFormatType(thisformat);
        }
    }
}
//------------------------------------------------------------------------------
// Scan matrix and set the width of different columns
//------------------------------------------------------------------------------
void SparseDisplay::SetWidth(Interpreter* interp)
{
    const hwMatrixS* mtx = m_currency.MatrixS();
    if (!mtx || mtx->Size() == 0 || mtx->IsEmpty())
    {
        return;
    }

    int nnz = mtx->NNZ();
    bool isreal = mtx->IsReal();
    bool isrealdata = mtx->IsRealData();

    _realwidth = 0;
    _imagwidth = 0;
    for (int i = 0; i < nnz; ++i)
    {
        if (interp && interp->IsInterrupt())
        {
            return;
        }
        int row = 0;
        int col = 0;


        double      realval = 0;
        double      imagval = 0;
        std::string imagsign;          

        if (isreal)
        {
            mtx->NZinfo(i, row, col, realval);
        }
        else
        {
            hwComplex cval(0, 0);
            mtx->NZinfo(i, row, col, cval);
            GetComplexNumberVals(cval, realval, imagval, imagsign);
        }

        std::string realstr(_formatvars.RealToString(realval));
        assert(!realstr.empty());

        _realwidth = std::max(_realwidth, static_cast<int>(realstr.size()));

        if (isrealdata)
        {
            continue;
        }

        std::string imagstr(_formatvars.RealToString(imagval));
        assert(!imagstr.empty());
        _imagwidth = std::max(_imagwidth, static_cast<int>(imagstr.size()) + 1);
    }
}
//------------------------------------------------------------------------------
// Returns true if paginating
//------------------------------------------------------------------------------
bool SparseDisplay::IsPaginating() const
{
    return (IsPaginatingRows());
}
//------------------------------------------------------------------------------
// Gets values as a string
//------------------------------------------------------------------------------
std::string SparseDisplay::GetValues(const OutputFormat* fmt) const
{
    std::string out = GetOutputNoPagination(fmt);
    return out;
}
//------------------------------------------------------------------------------
// Gets matrix data with forward pagination
//------------------------------------------------------------------------------
std::string SparseDisplay::GetOutputForwardPagination(const OutputFormat* fmt) const
{
    const hwMatrixS* mtx = m_currency.MatrixS();
    if (!mtx || mtx->Size() == 0 || mtx->IsEmpty())
    {
        return "";
    }

    int maxlines = std::max(GetNumRowsToFit() - 1, 1); // Add pagination msg
    if (m_linesPrinted >= m_maxRows - 1)
    {
        m_linesPrinted = m_maxRows - 2;  // Force print 1 line
    }

    int nnz = mtx->NNZ();
    int startrow = (m_rowBegin >= 0) ? m_rowBegin : 0;
    int totalrows = std::min(nnz, startrow + maxlines);

    if (nnz - totalrows == 1)  // Print the last line too instead of paginating
    {
        totalrows += 1;
    }

    // If there is only one row left, print instead of paginating
    if (m_parentDisplay)
    {
        if (totalrows == nnz - 1)
        {
            maxlines += 1;
            totalrows += 1;
        }
        else if (totalrows == nnz - 2)
        {
            maxlines += 2;
            totalrows += 2;
        }
    }

    bool isreal = mtx->IsReal();
    bool isrealdata = mtx->IsRealData();

    std::string myindent = GetIndentString(m_indent);
    std::string output;

    bool waspaginating = IsPaginating();
    int  prevstartrow = m_rowBegin;

    for (int i = m_rowBegin; i < totalrows && i < nnz; ++i)
    {
        output += "\n" + myindent;
        GetOutput(i, isreal, isrealdata, output);
        m_linesPrinted++;
        m_rowEnd = i;
    }

    std::string header;
    if (IsPaginating() || (waspaginating && prevstartrow > 0))
    {
        header = GetPaginationHeader(nnz);
    }
    return (header + output);
}
//------------------------------------------------------------------------------
// Gets outputfor back/up pagination
//------------------------------------------------------------------------------
std::string SparseDisplay::GetOutputBackPagination(const OutputFormat* fmt) const
{
    const hwMatrixS* mtx = m_currency.MatrixS();
    if (!mtx || mtx->Size() == 0 || mtx->IsEmpty())
    {
        return "";
    }

    int nnz = mtx->NNZ();
    int endrow = (m_rowEnd >= 0 && m_rowEnd < nnz) ? m_rowEnd : nnz - 1;
    int linestofit = GetNumRowsToFit();

    bool isreal = mtx->IsReal();
    bool isrealdata = mtx->IsRealData();

    std::string myindent (GetIndentString(m_indent));
    std::string output;

    m_rowEnd = endrow;
    for (int i = endrow; i >= 0 && m_linesPrinted <= linestofit; --i)
    {
        std::string rowdata("\n" + myindent);
        GetOutput(i, isreal, isrealdata, rowdata);
        m_linesPrinted++;
        m_rowBegin = i;

        output.insert(0, rowdata);
    }
    std::string header(GetPaginationHeader(nnz));
    return (header + output);
}
//------------------------------------------------------------------------------
// Returns true if display was paginating
//------------------------------------------------------------------------------
bool SparseDisplay::WasPaginating() const
{
    bool wasPaginating = IsPaginating();
    if (wasPaginating && m_rowBegin == 0 && m_rowEnd == 0)
    {
        wasPaginating = false;
    }
    return wasPaginating;
}
//------------------------------------------------------------------------------
// Gets pagination info for printing
//------------------------------------------------------------------------------
std::string SparseDisplay::GetPaginationHeader(int rows) const
{
    std::string msg;
    if (m_rowBegin >= 0 || m_rowEnd + 1 < rows)
    {
        int lastrow = std::min(m_rowEnd + 1, rows);
        std::ostringstream os;
        os << ", nnz rows [";
        if (m_rowBegin + 1 != lastrow)
        {
            os << m_rowBegin + 1 << ":";
        }
        os << lastrow << "]";
        msg = os.str();
    }
    return msg;
}
//------------------------------------------------------------------------------
// Sets data for forward display
//------------------------------------------------------------------------------
void SparseDisplay::SetForwardDisplayData()
{
    const hwMatrixS* mtx = m_currency.MatrixS();
    if (!mtx || mtx->Size() == 0 || mtx->IsEmpty())
    {
        return;
    }

    if (IsPaginatingRows())               // Row pagination
    {
        m_rowBegin = m_rowEnd + 1;
    }
    else if (m_rowEnd >= mtx->NNZ())
    {
        m_rowEnd = -1;
        m_rowBegin = -1;
        return;
    }

    // Make sure indices are positive
    m_rowBegin = std::max(m_rowBegin, 0);

    // Figure out how many lines have been printed and set it
    if (!m_parentDisplay ||  // First currency being paginated
        m_rowBegin > 0)      // Only this currency is being paginated
    {
        m_linesPrinted = 0;
    }
}
//------------------------------------------------------------------------------
// Sets indices for back pagination
//------------------------------------------------------------------------------
void SparseDisplay::SetBackDisplayData()
{
    m_linesPrinted = 0;     // Reset number of lines printed

    // We have reached the beginning of the matrix, so we cannot go back further
    if (m_rowBegin == 0)
    {
        return;
    }
    const hwMatrixS* mtx = m_currency.MatrixS();
    if (!mtx || mtx->Size() == 0 || mtx->IsEmpty())
    {
        return;
    }

    m_rowEnd = std::max(m_rowBegin - 1, 0);
}

// End of file:
