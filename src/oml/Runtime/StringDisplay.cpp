/**
* @file StringDisplay.cpp
* @date February 2018
* Copyright (C) 2018 Altair Engineering, Inc.  
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
#include "StringDisplay.h"

#include <algorithm>
#include <cassert>

#include "Currency.h"
#include "Interpreter.h"
#include "OML_Error.h"
#include "OutputFormat.h"

// End defines/includes

//------------------------------------------------------------------------------
// Constructor - Only currency or derived classes can access constructor
//------------------------------------------------------------------------------
StringDisplay::StringDisplay(const Currency& cur)
    : CurrencyDisplay(cur)
    , _chaineddisplay (0)
{
    m_rowEnd = -1;
}
//------------------------------------------------------------------------------
// Gets output
//------------------------------------------------------------------------------
std::string StringDisplay::GetOutput(const OutputFormat* fmt,
                                     std::ostringstream& os) const
{
    if (m_mode == DISPLAYMODE_QUIT || m_mode == DISPLAYMODE_EXIT)
    {
		return "";   // Do nothing
    }

    if (!IsValidDisplaySize())
    {
		return GetOutputNoPagination(fmt, os);
    }

    const_cast<StringDisplay*>(this)->Initialize(fmt, NULL, NULL); 

    std::string output;
	if (m_mode == DISPLAYMODE_FORWARD || m_mode == DISPLAYMODE_DOWN)
    {
        output = GetOutputForwardPagination(fmt, os);
    }
	else if (m_mode == DISPLAYMODE_BACK || m_mode == DISPLAYMODE_UP)
    {
		output = GetOutputBackPagination(fmt, os);
    }

    return output;
}
//------------------------------------------------------------------------------
// Gets output with no pagination
//------------------------------------------------------------------------------
std::string StringDisplay::GetOutputNoPagination(const OutputFormat* fmt,
                                                 std::ostringstream& os) const
{
    os << m_currency.StringVal();
    std::string output (os.str());
	return output;
}
//------------------------------------------------------------------------------
//! Sets data for forward display
//------------------------------------------------------------------------------
void StringDisplay::SetForwardDisplayData()
{
    if (m_rowEnd == -1)
    {
        m_rowEnd = 0;
    }

    // Figure out how many lines have been printed and set it
    if (!m_parentDisplay ||   // First currency being paginated
        m_rowEnd > 0)         // Only this currency is being paginated
    {
       m_linesPrinted = 0;
    }
}
//------------------------------------------------------------------------------
//! Sets data for back pagination
//------------------------------------------------------------------------------
void StringDisplay::SetBackDisplayData()
{
    m_linesPrinted = 0;          // Reset lines printed
}
//------------------------------------------------------------------------------
// Gets data with back pagination
//------------------------------------------------------------------------------
std::string StringDisplay::GetOutputBackPagination(const OutputFormat* fmt,
                                                   std::ostringstream& os) const
{
    if (m_rowEnd == -1 || m_rowEnd == 0 || m_rowEnd <= m_maxCols)
    {
        m_linesPrinted = 0;
        m_rowEnd       = 0;    
        return GetOutputForwardPagination(fmt, os);
    }

    std::string str (m_currency.StringVal());
    std::string sub = str.substr(0, m_rowEnd);
    if (!sub.empty() && sub[sub.size() - 1] == '\n')
    {
        sub.pop_back();
        m_rowEnd --;
    }
    if (sub.empty())
    {
        return GetOutputForwardPagination(fmt, os);
    }

    std::string output;
    int         linestofit = GetNumRowsToFit();
    while (m_linesPrinted < linestofit && !sub.empty())
    {
        size_t pos = sub.rfind("\n");
        m_linesPrinted ++;
        if (pos == std::string::npos)
        {
            if (output.empty())
            {
                output = "\n" + sub;
            }
            else
            {
                output.insert(0, "\n" + sub);
            }
            break;
        }

        std::string line = sub.substr(pos);
        if (output.empty())
        {
            output = line;
        }
        else
        {
            output.insert(0, line);
        }
        if (pos > 0)
        {
            sub = sub.substr(0, pos);
        }
        else
        {
            break;
        }
    }
    if (m_linesPrinted > 0)
    {
        m_rowEnd -= static_cast<int>(output.length());
    }
    if (m_rowEnd < 0)
    {
        m_rowEnd = 0;
    }

    return os.str() + output;
}
//------------------------------------------------------------------------------
// Gets data with forward pagination
//------------------------------------------------------------------------------
std::string StringDisplay::GetOutputForwardPagination(const OutputFormat* fmt,
                                                      std::ostringstream& os) const
{  
    if (m_rowEnd == -1)
    {
        m_rowEnd = 0;
    }

    std::string str (m_currency.StringVal());
    std::string sub = str.substr(m_rowEnd);
    if (sub.empty())
    {
        m_linesPrinted++;
        return os.str();
    }

    std::string output;
    int         linestofit = GetNumRowsToFit();
    size_t      len        = sub.length();
    int         totalrows  = m_linesPrinted + linestofit;

    for (size_t i = 0; i < len && m_linesPrinted < totalrows; ++i, ++m_rowEnd)
    {
        char ch = sub[i];
        output += ch;

#ifdef OS_WIN
        if (ch == '\n')
#else
        if (ch == '\n' || ch == '\r')
#endif
        {
            ++m_linesPrinted;
        }
    }
    std::string out = os.str() + output;

    return out;
}
//------------------------------------------------------------------------------
// Gets values as a string
//------------------------------------------------------------------------------
std::string StringDisplay::GetValues(const OutputFormat* fmt) const
{
    return m_currency.StringVal();
}
//------------------------------------------------------------------------------
// Returns true if this is a chained display, which needs to be deleted
//------------------------------------------------------------------------------
bool StringDisplay::IsChainedDisplay(CurrencyDisplay* display) const
{
    return (display && display->GetChainedDisplay() == _chaineddisplay);
}