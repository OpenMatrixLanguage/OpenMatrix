/**
* @file CellDisplay.cxx
* @date June 2016
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

// Begin defines/includes
#include "CellDisplay.h"

#include <cassert>

#include "Currency.h"
#include "Interpreter.h"
#include "OutputFormat.h"
#include "SignalHandlerBase.h"

// End defines/includes

//------------------------------------------------------------------------------
// Constructor - Only currency or derived classes can access constructor
//------------------------------------------------------------------------------
CellDisplay::CellDisplay(const Currency& cur)
    : CurrencyDisplay(cur)
{
}
//------------------------------------------------------------------------------
// Gets number of rows and cols in currency
//------------------------------------------------------------------------------
void CellDisplay::GetCurrencySize(int& rows, int& cols) const
{
    HML_CELLARRAY* cells = m_currency.CellArray();
    if (!cells) return;

    rows = cells->M();
    cols = cells->N();
}
//------------------------------------------------------------------------------
// Gets output - called from Currency::GetOutputString
//------------------------------------------------------------------------------
std::string CellDisplay::GetOutput(const OutputFormat* fmt, 
                                   std::ostringstream& os) const
{
    if (m_mode == DISPLAYMODE_QUIT || m_mode == DISPLAYMODE_EXIT) return ""; // Do nothing

    int numrows = 0;
    int numcols = 0;
    GetCurrencySize(numrows, numcols);

    if (numrows == 0 || numcols == 0 || !IsValidDisplaySize())
    {
		return GetOutputNoPagination(fmt, os);
    }
    const_cast<CellDisplay*>(this)->Initialize(fmt, NULL, NULL); 

    if (m_mode == DISPLAYMODE_FORWARD || m_mode == DISPLAYMODE_DOWN)
        return GetOutputForwardPagination(fmt, os);

    if (m_mode == DISPLAYMODE_BACK || m_mode == DISPLAYMODE_UP)
        return GetOutputBackPagination(fmt);

    return "";
}
//------------------------------------------------------------------------------
// Gets output with no pagination
//------------------------------------------------------------------------------
std::string CellDisplay::GetOutputNoPagination(const OutputFormat* fmt,
                                               std::ostringstream& os) const
{
    HML_CELLARRAY* cells = m_currency.CellArray();
    assert(cells);
    
    int numrows = cells->M();
    int numcols = cells->N();

    int         cindent     = m_indent + 1;
    std::string childindent = GetIndentString(cindent);
    std::string myindent    = GetIndentString(m_indent);

    os << std::endl << myindent << "{" << std::endl;

	for (int i = 0; i < numrows; ++i)
	{	
		for (int j = 0; j < numcols; ++j)
		{
            Currency         cur     = (*cells)(i,j);
            CurrencyDisplay* display = cur.GetDisplay();
            if (display)
            {
                display->SetIndent(cindent);
            }

            os << childindent << "[" << i+1  << "," << j+1 << "] " 
               << cur.GetOutputString(fmt) << std::endl;
                        
            CurrencyDisplay::DeleteDisplay(display);
            cur.SetDisplay(nullptr);
        }
	}
    os << myindent << "}";
    std::string output (os.str());
	return output;
}
//------------------------------------------------------------------------------
// Sets data for forward display
//------------------------------------------------------------------------------
void CellDisplay::SetForwardDisplayData()
{
    int numrows = 0;
    int numcols = 0;
    GetCurrencySize(numrows, numcols);

	if (IsPaginatingCols())                    // Col pagination	
    {
		m_colBegin = m_colEnd + 1;
        m_rowBegin = std::max(0, m_rowEnd);    // Start with row that ended
    }
	else if (IsPaginatingRows())               // Row pagination
	{
        m_colBegin = 0;                        // Start row pagination
		m_rowBegin = m_rowEnd + 1;
	}
    else if (m_colEnd >= numcols && m_rowEnd >= numrows)
    {
        m_colEnd = -1;
        m_rowEnd = -1;
        m_rowBegin = -1;
        m_colBegin = -1;
        return;
    }

    // Make sure indices are positive
    m_rowBegin = std::max(m_rowBegin, 0);
    m_colBegin = std::max(m_colBegin, 0);

    // Reset the end indices
    m_colEnd = -1;
    m_rowEnd = -1;

    // Figure out how many lines have been printed and set it
    if (!m_parentDisplay ||                // First currency being paginated
       (m_rowBegin > 0 || m_colBegin > 0)) // Only this currency is being paginated
       m_linesPrinted = 0;
}
//------------------------------------------------------------------------------
// Gets data with forward pagination
//------------------------------------------------------------------------------
std::string CellDisplay::GetOutputForwardPagination(const OutputFormat* fmt,
                                                    std::ostringstream& os) const
{
    HML_CELLARRAY* cells = m_currency.CellArray();
    assert(cells);
    
    const_cast<CellDisplay*>(this)->UpdateNumLinesPrinted();

    bool printclosebraces = false;

    int         cindent     = m_indent + 1;
    std::string childindent = GetIndentString(cindent);
    std::string myindent    = GetIndentString(m_indent);

    if (m_colBegin == 0 && m_rowBegin == 0) // Need header only for first row/col
    {
        printclosebraces = true;
        os << std::endl << myindent << "{" << std::endl;
    }
    else
    {
        os.str("");
    }

    int numrows = cells->M();
    int numcols = cells->N();

    int linestofit = GetNumRowsToFit();
    int totalrows  = m_linesPrinted + linestofit;
    int startcol   = m_colBegin;
	for (int i = m_rowBegin; i < numrows && m_linesPrinted < totalrows; ++i)
	{	
        m_rowEnd = i;
        if (i != m_rowBegin)
            os << std::endl;

		for (int j = startcol; j < numcols && m_linesPrinted < totalrows; ++j)
		{
            m_colEnd = j;

            if (j != startcol)
                os << std::endl;

            Currency cur = (*cells)(i,j);
            bool canpaginate = CanPaginate(cur);
            if (canpaginate)
            {
                CurrencyDisplay* disp = cur.GetDisplay();
                if (disp)
                {
                    disp->SetParentDisplay(const_cast<CellDisplay*>(this));
                    disp->SetIndent(cindent);
                }
            }

            std::string tmp (cur.GetOutputString(fmt));
            StripEndline(tmp);
           // if (!tmp.empty())
            {
                os << childindent << "["<< i+1 << ","<< j+1 << "] " << tmp;
            }

            if (!canpaginate)
                m_linesPrinted++;

            if (canpaginate)
            {
                CurrencyDisplay* nesteddisplay = cur.GetDisplay();
                if (nesteddisplay)
                {
                    if (nesteddisplay->IsPaginatingCols() ||
                        nesteddisplay->IsPaginatingRows())
                    {
                        if (m_signalHandler)
                            m_signalHandler->OnAddDisplayHandler(nesteddisplay);
                        std::string output (os.str());
                        return output;
                    }
                }
            }
            CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
            cur.SetDisplay(NULL);

            // Add trailing '}' for cells
            //if (cur.IsCellArray())
            //{
            //    os << std::endl;
            //    os << "}";
            //    m_linesPrinted++;
            //}
        }
        startcol = 0;  // Reset the columns
	}

    if (printclosebraces && !IsPaginatingRows() && !IsPaginatingCols())
    {
        os << std::endl;
        os << myindent << "}";
        m_linesPrinted++;
    }
    std::string output (os.str());

    return output;
}
//------------------------------------------------------------------------------
//! Sets data for back pagination
//------------------------------------------------------------------------------
void CellDisplay::SetBackDisplayData()
{
    m_linesPrinted = 0;

	// We have reached the beginning of the matrix, so we cannot go back further
	if (m_rowBegin == 0 && m_colBegin == 0) return;

    int numcols = 0;
    int numrows = 0;
    GetCurrencySize(numrows, numcols);

    if (m_colBegin >= 0 && m_colBegin < numcols-1) // Paginating cols
    {
        m_colEnd = std::max(m_colBegin - 1,  0);
        m_rowEnd = m_rowBegin;

        if (m_colEnd <= 0)
        {
            m_rowEnd = std::max(m_rowEnd - 1, 0);
            m_colEnd = numcols - 1; // Start row pagination
        }
    }
    else 
    {
        m_rowEnd = std::max(m_rowBegin - 1, 0);
        m_colEnd = numcols - 1; // Start row pagination
    }
}
//------------------------------------------------------------------------------
// Gets data with back pagination
//------------------------------------------------------------------------------
std::string CellDisplay::GetOutputBackPagination(const OutputFormat* fmt) const
{
    HML_CELLARRAY* cells = m_currency.CellArray();
    assert(cells);
    
    int numrows = cells->M();
    int numcols = cells->N();

	int endrow     = (m_rowEnd >= 0 && m_rowEnd < numrows) ? m_rowEnd : numrows - 1;
    int endcol     = m_colEnd;
    int linestofit = GetNumRowsToFit();

    int         cindent     = m_indent + 1;
    std::string childindent = GetIndentString(cindent);
    std::string myindent    = GetIndentString(m_indent);

    std::string output;
	for (int i = endrow; i >= 0 && m_linesPrinted <= linestofit; --i)
	{	
        m_rowBegin = i;
		for (int j = endcol; j >= 0 && m_linesPrinted <= linestofit; --j)
		{
            std::ostringstream os;
            m_colBegin = j;

            Currency cur = (*cells)(i,j);
            bool canpaginate = CanPaginate(cur);
            if (canpaginate)
            {
                CurrencyDisplay* disp = cur.GetDisplay();
                if (disp)
                {
                    disp->SetParentDisplay(const_cast<CellDisplay*>(this));
                    disp->SetIndent(cindent);
                }
            }
            std::string tmp(cur.GetOutputString(fmt));
            StripEndline(tmp);
            if (!tmp.empty())
            {
                os << childindent << "["<< i+1 << ","<< j+1 << "] " << tmp;
            }

            if (!(j == endcol && i == endrow))
                os << std::endl;

            if (output.empty())
                output = os.str();
            else
                output.insert(0, std::string(os.str()));

            if (canpaginate)
            {
                CurrencyDisplay* nesteddisplay = cur.GetDisplay();
                if (nesteddisplay)
                {
                    if (nesteddisplay->IsPaginatingCols() ||
                        nesteddisplay->IsPaginatingRows())
                    {
                        if (m_signalHandler)
                            m_signalHandler->OnAddDisplayHandler(nesteddisplay);
                        return output;
                    }
                }
            }
            else
                m_linesPrinted++;

            CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
            cur.SetDisplay(NULL);
        }
        endcol = numcols -1;
	}

    return output;
}
//------------------------------------------------------------------------------
// Returns true if end of pagination message needs to be printed
//------------------------------------------------------------------------------
bool CellDisplay::GetPaginationEndMsg(std::string& msg) const
{
    int totalrows = 0;
    int totalcols = 0;
    GetCurrencySize(totalrows, totalcols);

    bool printmsg = CurrencyDisplay::GetPaginationEndMsg(msg);

    std::string myindent (GetIndentString(m_indent));
    if (!printmsg)
    {
        // Add trailing '}' if the currency can paginate
        if (!(totalrows == 0 || totalcols == 0)) 
            msg = myindent + "}";
        return false; // No need to print end of pagination message
    }

    if (m_rowEnd >= totalrows - 1 && m_colEnd >= totalcols - 1)
        msg = myindent + "}";

    return true;
}
//------------------------------------------------------------------------------
// Updates the number of rows to fit
//------------------------------------------------------------------------------
void CellDisplay::UpdateNumLinesPrinted()
{
    if (m_rowBegin == 0 && m_colBegin == 0)
    {
        m_linesPrinted ++;       // '\n{\n'

        if (!m_currency.IsDispOutput())
            m_linesPrinted ++;  // Assignment statement
    }
}
//------------------------------------------------------------------------------
// Gets values as a string
//------------------------------------------------------------------------------
std::string CellDisplay::GetValues(const OutputFormat* fmt) const
{
    HML_CELLARRAY* cells = m_currency.CellArray();
    assert(cells);
    if (!cells) return "";
    
    std::ostringstream os;

    int numrows = cells->M();
    int numcols = cells->N();

    for (int j = 0; j < numcols; ++j)
	{	
		for (int i = 0; i < numrows; ++i)
		{
            Currency cur = (*cells)(i,j);

            os << "["<< i+1 << "," << j+1 << "] " << cur.GetValues(fmt);
            
            if (i < numrows - 1)
                os << std::endl;
                        
            CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
            cur.SetDisplay(NULL);
        }
	}

    std::string output (os.str());
	return output;
}
// End of file:

