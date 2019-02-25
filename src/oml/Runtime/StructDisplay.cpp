/**
* @file StructDisplay.cpp
* @date February 2016
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
#include "StructDisplay.h"

#include <cassert>
#include <map>

#include "Currency.h"
#include "SignalHandlerBase.h"
#include "StructData.h"
// End defines/includes

//------------------------------------------------------------------------------
// Constructor - Only currency or derived classes can access constructor
//------------------------------------------------------------------------------
StructDisplay::StructDisplay(const Currency& cur)
    : CurrencyDisplay(cur)
{
}
//------------------------------------------------------------------------------
// Gets number of rows and cols in currency
//------------------------------------------------------------------------------
void StructDisplay::GetCurrencySize(int& rows, int& cols) const
{
    StructData* structdata = m_currency.Struct();
    assert(structdata);

    if (structdata->Size() != 1)
    {
        rows = structdata->M();
        cols = structdata->N();
    }
    else
    {
        std::map<std::string, int> fields (structdata->GetFieldNames());
        cols =  (fields.empty() ? 0 : static_cast<int>(fields.size()));
        rows = 1; 
    }
}
//------------------------------------------------------------------------------
// Gets output - called from Currency::GetOutputString
//------------------------------------------------------------------------------
std::string StructDisplay::GetOutput(const OutputFormat* fmt, 
                                     std::ostringstream& os) const
{
    if (m_mode == DISPLAYMODE_QUIT || m_mode == DISPLAYMODE_EXIT) 
    {
        return ""; // Do nothing
    }

    if (!CanPaginate(m_currency) || !IsValidDisplaySize())
		return GetOutputNoPagination(fmt, os);

    const_cast<StructDisplay*>(this)->Initialize(fmt, NULL, NULL); 

    if (m_mode == DISPLAYMODE_FORWARD || m_mode == DISPLAYMODE_DOWN)
        return GetOutputForwardPagination(fmt, os);

    if (m_mode == DISPLAYMODE_BACK || m_mode == DISPLAYMODE_UP)
        return GetOutputBackPagination(fmt);

    return "";
}
//------------------------------------------------------------------------------
// Gets output with no pagination
//------------------------------------------------------------------------------
std::string StructDisplay::GetOutputNoPagination(const OutputFormat* fmt,
                                                 std::ostringstream& os) const
{
    StructData* sd = m_currency.Struct();
    assert(sd);

    std::map<std::string, int> fields (sd->GetFieldNames());

    std::string myindent   = GetIndentString(m_indent);
    bool        isindented = (!myindent.empty());

	if (m_currency.IsObject())
		os << "object [\n";
	else
		os << "struct [\n";

    int structsize = sd->Size();
    if (structsize != 1 && (!fields.empty() || sd->M() > 0 || sd->N() > 0))
    { 
        // Just print the name of the fields
        os << myindent << "Struct array of size " << sd->M() << " x " << sd->N();
		os << " with fields:" << "\n";
    }

    for (std::map<std::string, int>::const_iterator itr = fields.begin(); 
         itr != fields.end(); ++itr)
    {
		if (itr != fields.begin())
			os << "\n";

        std::string name (itr->first);
        os << myindent << name;

        if (structsize != 1) continue; // Show values only for size = 1

        Currency cur (sd->GetValue(0, 0, name));
        CurrencyDisplay* display = cur.GetDisplay();
        if (display && isindented)
        {
            display->SetIndent(m_indent);
        }

        os << ": " << cur.GetOutputString(fmt);

        CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
        cur.SetDisplay(NULL);
	}

    os << std::endl << myindent << "]";
	return std::string(os.str());
}
//------------------------------------------------------------------------------
// Returns true if end of pagination message needs to be printed
//------------------------------------------------------------------------------
bool StructDisplay::GetPaginationEndMsg(std::string& msg) const
{
    int totalrows = 0;
    int totalcols = 0;
    GetCurrencySize(totalrows, totalcols);

    std::string myindent = GetIndentString(m_indent);
    if (m_colBegin == 0 && m_colEnd >= totalcols - 1)
    {
        msg = myindent + "]";
        return false;  // No need to print end of pagination message
    }

    if (m_colEnd >= totalcols - 1)
    {
        msg = myindent + "]";
    }
    return true;
}
//------------------------------------------------------------------------------
// Sets data for forward display
//------------------------------------------------------------------------------
void StructDisplay::SetForwardDisplayData()
{
    int numrows = 0;
    int numcols = 0;
    GetCurrencySize(numrows, numcols);

    m_rowBegin = 0;  // Structs do not paginate rows
    m_rowEnd   = 0;

	if (IsPaginatingCols())                    // Col pagination	
		m_colBegin = m_colEnd + 1;    
    else if (m_colEnd >= numcols)
    {
        m_colEnd = -1;
        m_rowEnd = -1;
        m_rowBegin = -1;
        m_colBegin = -1;
        return;
    }

    m_colBegin = std::max(m_colBegin, 0);  // Make sure indices are positive
    m_colEnd   = -1;                       // Reset the end indices

    //// Figure out how many lines have been printed and set it
    //if (!m_parentDisplay ||   // First currency being paginated
    //    m_colBegin > 0)       // Only this currency is being paginated
       m_linesPrinted = 0;
}
//------------------------------------------------------------------------------
// Sets data for back pagination
//------------------------------------------------------------------------------
void StructDisplay::SetBackDisplayData()
{
    m_linesPrinted = 0;          // Reset lines printed
    m_rowBegin     = 0;          // Structs do not paginate rows
    m_rowEnd       = 0;

	if (m_colBegin == 0) return; // Reached the beginning, can't go back further

    int numcols = 0;
    int numrows = 0;
    GetCurrencySize(numrows, numcols);

    if (m_colBegin == m_colEnd)             // Just finished paginating last col
        m_colEnd = std::max(m_colBegin, 0);

    else if (m_colBegin >= 0 && m_colBegin < numcols-1) // Paginating cols
        m_colEnd = std::max(m_colBegin - 1,  0);
}
//------------------------------------------------------------------------------
// Gets data with back pagination
//------------------------------------------------------------------------------
std::string StructDisplay::GetOutputBackPagination(const OutputFormat* fmt) const
{
    StructData* sd = m_currency.Struct();
    assert(sd);

    int numcols = 0;
    int numrows = 0;
    GetCurrencySize(numrows, numcols);

    // Note: for struct paginations, rows are always 1
    m_rowBegin = 0;
    m_rowEnd   = 0;

    int endcol     = m_colEnd;
    int linestofit = GetNumRowsToFit();

    std::string output;
    std::map<std::string, int> fields (sd->GetFieldNames());
    std::map<std::string, int>::reverse_iterator ritr = fields.rbegin();

    std::string myindent = GetIndentString(m_indent);
    bool        isindented = (!myindent.empty());

    for (int j = numcols - 1; j >= 0 && m_linesPrinted <= linestofit && ritr != fields.rend();
         --j, ++ritr)
	{ 
        if (j > m_colEnd) continue;  // Already printed

        std::ostringstream os;
        m_colBegin = j;

        std::string field (ritr->first);

        Currency cur (sd->GetValue(0, 0, field));
        CurrencyDisplay* disp = cur.GetDisplay();
        if (disp && isindented)
        {
            disp->SetIndent(m_indent);
        }

        bool canpaginate = CanPaginate(cur);
        if (canpaginate)
        {
            if (disp)
            {
                disp->SetParentDisplay(const_cast<StructDisplay*>(this));
            }
        }
        std::string tmp (cur.GetOutputString(fmt));
        StripEndline(tmp);
        os << myindent << field << ": " << tmp;

        if (j != endcol)
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

    return output;
}
//------------------------------------------------------------------------------
//! Gets data with forward pagination
//! \param[in] fmt Format
//! \param[in] os  Output stream
//------------------------------------------------------------------------------
std::string StructDisplay::GetOutputForwardPagination(const OutputFormat* fmt,
                                                      std::ostringstream& os) const
{
    StructData* sd = m_currency.Struct();
    assert (sd);

    const_cast<StructDisplay*>(this)->UpdateNumLinesPrinted();

    std::string myindent   = GetIndentString(m_indent);
    bool        isindented = (!myindent.empty());

    bool printclosebraces = false;
    if (m_colBegin == 0) // Need header only for first col
    {
        printclosebraces = true;
        os << "struct [" << std::endl;
    }
    else
    {
        os.str("");
    }
    
    int numcols = 0;
    int numrows = 0;
    GetCurrencySize(numrows, numcols);

    int linestofit = GetNumRowsToFit();
    int totalrows  = m_linesPrinted + linestofit;
    int startcol   = m_colBegin;

    std::map<std::string, int> fields (sd->GetFieldNames());

    // Note: for struct paginations, rows are always 1
    m_rowBegin = 0;
    m_rowEnd   = 0;

    std::map<std::string, int>::const_iterator itr = fields.begin();
    for (int j = 0; j < numcols && m_linesPrinted < totalrows && itr != fields.end();
         ++j, ++itr)
    {
        if (j < startcol) continue;  // Already printed
            
        m_colEnd = j;

        if (j != startcol)
            os << std::endl;

        std::string field (itr->first);

        Currency cur (sd->GetValue(0, 0, field));
        CurrencyDisplay* disp = cur.GetDisplay();
        if (disp && isindented)
        {
            disp->SetIndent(m_indent);
        }

        bool canpaginate = CanPaginate(cur);
        if (disp && canpaginate)
        {
           disp->SetParentDisplay(const_cast<StructDisplay*>(this));
        }

        std::string tmp (cur.GetOutputString(fmt));
        StripEndline(tmp);
        os << myindent << field << ": " << tmp;

        if (!canpaginate)
            m_linesPrinted++;

        else // Nested pagination
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
    }

    if (printclosebraces && !IsPaginatingRows() && !IsPaginatingCols())
    {
        os << std::endl << myindent << "]";
        m_linesPrinted++;
    }

    std::string output (os.str());
    return output;
}
//------------------------------------------------------------------------------
//! Updates the number of rows to fit
//------------------------------------------------------------------------------
void StructDisplay::UpdateNumLinesPrinted()
{
    if (m_rowBegin == 0 && m_colBegin == 0)
        m_linesPrinted ++;       // Assignment
}
//------------------------------------------------------------------------------
//! Gets values as a string
//------------------------------------------------------------------------------
std::string StructDisplay::GetValues(const OutputFormat* fmt) const
{
    StructData* sd = m_currency.Struct();
    assert(sd);
    if (!sd) return "";

    int structsize = sd->Size();
    std::ostringstream os;
    std::map<std::string, int> fields (sd->GetFieldNames());

    for (std::map<std::string, int>::const_iterator itr = fields.begin(); 
         itr != fields.end(); ++itr)
    {
		if (itr != fields.begin())
			os << "\n";

        std::string name (itr->first);
        os << name;

        if (structsize != 1) continue; // Show values only for size = 1

        Currency cur (sd->GetValue(0, 0, name));
        os << ": " << cur.GetValues(fmt);

        CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
        cur.SetDisplay(NULL);
	}

	return std::string(os.str());

}
// End of file:

