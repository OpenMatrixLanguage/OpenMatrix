/**
* @file CellNDisplay.cxx
* @date January 2019
* Copyright (C) 2019-2021 Altair Engineering, Inc.
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

#include "CellNDisplay.h"

#include <cassert>
#include <memory>        // For std::unique_ptr

#include "Interpreter.h"
#include "OutputFormat.h"
#include "SignalHandlerBase.h"


//------------------------------------------------------------------------------
// Constructor - Only currency or derived classes can access constructor
//------------------------------------------------------------------------------
CellNDisplay::CellNDisplay(const Currency& cur)
	: CurrencyDisplay(cur)
    , _childPaginating(false)
    , _lastRowPrinted (-1)
{
}
//------------------------------------------------------------------------------
// Gets number of rows and cols in currency
//------------------------------------------------------------------------------
void CellNDisplay::GetCurrencySize(int& rows, int& cols) const
{
	HML_ND_CELLARRAY* cell = m_currency.CellArrayND();
	if (!cell)
	{
		return;
	}

	std::vector<int> dims = cell->Dimensions();
	if (dims.empty())
	{
		return;
	}
	size_t numdims = dims.size();

	cols = 1;
	rows = 1;

	if (numdims == 3)
	{
		rows = dims[2];        // 3D cell array
	}
	else if (numdims > 3)      // ND cell array
	{
		for (std::vector<int>::const_iterator itr = dims.begin() + 2;
			itr != dims.end(); ++itr)
		{
			rows *= *itr;
		}
	}
}
//------------------------------------------------------------------------------
// Gets output - called from Currency::GetOutputString
//------------------------------------------------------------------------------
std::string CellNDisplay::GetOutput(const OutputFormat* fmt, 
	                                std::ostringstream& os) const
{
    _childPaginating = false;
	if (m_mode == DISPLAYMODE_QUIT || m_mode == DISPLAYMODE_EXIT)
	{
		return "";
	}

	int numrows = 0;
	int numcols = 0;
	GetCurrencySize(numrows, numcols);

	if (numrows == 0 || numcols == 0)
	{
		if (m_currency.IsDispOutput())
		{
			return "";
		}
		return std::string(os.str());
	}

    bool paginate = IsValidDisplaySize();
    if (!paginate)
    {
        m_rowBegin = 0;
        m_rowEnd = 0;
        m_colBegin = 0;
        m_colEnd = 0;
    }

	std::string header;
	if (m_rowBegin < 1 && !m_currency.IsDispOutput())
	{

        header = os.str();
        if (!header.empty() || m_parentDisplay)
        {
            header += '\n';
        }
	}

    std::string output;
    if (!paginate)
    {
        output = GetOutputNoPagination(fmt);
    }
    else
    {
        const_cast<CellNDisplay*>(this)->Initialize(fmt, NULL, NULL);

        if (m_mode == DISPLAYMODE_FORWARD || m_mode == DISPLAYMODE_DOWN)
        {
            output = GetOutputForwardPagination(fmt);
        }
        else if (m_mode == DISPLAYMODE_BACK || m_mode == DISPLAYMODE_UP)
        {
            output = GetOutputBackPagination(fmt);
        }
    }
    return (header + output);
}
//------------------------------------------------------------------------------
// Helper to slice a given ND cell array to 2D cells
//------------------------------------------------------------------------------
void CellNDisplay::GetSlices(HML_ND_CELLARRAY*         cell,
	                         std::vector<Currency>&    curs,
	                         std::vector<std::string>& labels)
{
	if (!cell)
	{
		return;
	}

	std::vector<hwSliceArg> base;
	base.push_back(hwSliceArg());
	base.push_back(hwSliceArg());

	GetSlicesHelper(cell, base, curs, labels);

	assert(curs.size() == labels.size());
}
//------------------------------------------------------------------------------
// Helper to slice a given ND cell array to 2D cells
//------------------------------------------------------------------------------
void CellNDisplay::GetSlicesHelper(HML_ND_CELLARRAY*         cell,
	                               std::vector<hwSliceArg>   slices,
	                               std::vector<Currency>&    curs,
	                               std::vector<std::string>& labels)
{
	if (!cell)
	{
		return;
	}

	std::vector<int> dims    = cell->Dimensions();
	size_t           numdims = dims.empty() ? 0 : dims.size();

	if (numdims < 3)
	{
		// Get the slice label
		std::string label("slice(");
		for (std::vector<hwSliceArg>::const_iterator itr = slices.begin();
			itr != slices.end(); ++itr)
		{
			hwSliceArg slice = *itr;

			if (itr != slices.begin())
			{
				label += ", ";
			}

			if (slice.IsColon())
			{
				label += ":";
			}
			else
			{
				label += std::to_string(static_cast<long long>(slice.Scalar() + 1));
			}
		}
		label += ")";

		labels.push_back(label);                        // Slice label for cell

		HML_CELLARRAY* temp = EvaluatorInterface::allocateCellArray();
		cell->ConvertNDto2D(*temp);
		curs.push_back(temp);  // Cell type currency

		return;
	}

	size_t curDimIndex = slices.size();
	if (curDimIndex == numdims)
	{
		std::unique_ptr<HML_ND_CELLARRAY> tmpslice(new HML_ND_CELLARRAY);
		cell->SliceRHS(slices, *(tmpslice.get()));
		CellNDisplay::GetSlicesHelper(tmpslice.release(), slices, curs, labels);
	}
	else if (curDimIndex < numdims)
	{
		int curdim = dims[curDimIndex];

		for (int j = 0; j < curdim; ++j)
		{
			slices.push_back(j);
			CellNDisplay::GetSlicesHelper(cell, slices, curs, labels);
			slices.pop_back();
		}
	}
}
//------------------------------------------------------------------------------
// Gets output with no pagination
//------------------------------------------------------------------------------
std::string CellNDisplay::GetOutputNoPagination(const OutputFormat* fmt) const
{
	HML_ND_CELLARRAY* cell = m_currency.CellArrayND();
	if (!cell)
	{
		return "";
	}

	std::vector<std::string> labels;
	std::vector<Currency>    curs;
	GetSlices(cell, curs, labels);

	std::ostringstream os;

	if (cell->Size() == 1)
	{
		Currency temp((*cell)(0));
		os << temp.GetOutputString(fmt) << '\n';
	}
	else
	{
		size_t numlabels = labels.empty() ? 0 : labels.size();
		std::vector<Currency>::const_iterator itr = curs.begin();

		for (size_t i = 0; itr != curs.end() && i < numlabels; ++itr, ++i)
		{
			Currency cur(*itr);

			os << labels[i] << " = ";

			std::string tmp(cur.GetOutputString(fmt));
			StripEndline(tmp);
			if (!tmp.empty())
			{
				os << tmp << '\n';
			}

			CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
			cur.SetDisplay(nullptr);
		}
	}

	std::string tmp(os.str());
	StripEndline(tmp);

	return tmp;
}
//------------------------------------------------------------------------------
// Gets values as a string
//------------------------------------------------------------------------------
std::string CellNDisplay::GetValues(const OutputFormat* fmt) const
{
	HML_ND_CELLARRAY* cell = m_currency.CellArrayND();
	if (!cell)
	{
		return "";
	}

	std::vector<std::string> labels;
	std::vector<Currency>    curs;
	GetSlices(cell, curs, labels);

	if (cell->Size() == 1)
	{
		Currency temp((*cell)(0));
		return temp.GetValues(fmt);
	}

	std::ostringstream os;
	size_t numlabels = labels.empty() ? 0 : labels.size();
	std::vector<Currency>::const_iterator itr = curs.begin();

	for (size_t i = 0; itr != curs.end() && i < numlabels; ++itr, ++i)
	{
		Currency cur(*itr);
		os << labels[i] << " = " << cur.GetValues(fmt) << std::endl;

		CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
		cur.SetDisplay(nullptr);
	}

	return os.str();
}
//------------------------------------------------------------------------------
// Sets data for forward display
//------------------------------------------------------------------------------
void CellNDisplay::SetForwardDisplayData()
{
    int numrows = 0;
    int numcols = 0;
    GetCurrencySize(numrows, numcols);

    m_colBegin = 0;

    // Figure out how many lines have been printed and set it
    if (!m_parentDisplay ||                // First currency being paginated
        (m_rowBegin > 0 || m_colBegin > 0)) // Only this currency is being paginated
    {
        m_linesPrinted = 0;
    }


    if (IsPaginatingRows())               // Row pagination
    {
        m_rowBegin = m_rowEnd;
    }
    else if (m_rowEnd >= numrows || m_rowBegin >= numrows)
    {
        m_colEnd = -1;
        m_rowEnd = -1;
        m_rowBegin = -1;
        m_colBegin = -1;
        _lastRowPrinted = -1;
        return;
    }
    else if (m_rowBegin < m_rowEnd)
    {
        m_rowBegin += 1;
    }

    // Make sure indices are positive
    m_rowBegin = std::max(m_rowBegin, 0);

    // Reset the end indices
    m_colEnd = -1;
    m_rowEnd = -1;
}
//------------------------------------------------------------------------------
// Gets data with forward pagination
//------------------------------------------------------------------------------
std::string CellNDisplay::GetOutputForwardPagination(const OutputFormat* fmt) const
{
    HML_ND_CELLARRAY* cell = m_currency.CellArrayND();
    if (!cell)
    {
        return "";
    }

    std::vector<std::string> labels;
    std::vector<Currency>    curs;
    GetSlices(cell, curs, labels);

    if (HasOnlyEmptySlices(curs))
    {
        return GetOutputEmpty();
    }

    int numrows = 0;
    int numcols = 1;
    CellNDisplay::GetCurrencySize(numrows, numcols);

    int linestofit = GetNumRowsToFit();
    int totalrows = m_linesPrinted + linestofit;
    m_colEnd = 0;


    std::ostringstream os;
    std::vector<Currency>::iterator itr = curs.begin() + m_rowBegin;
    int numlabels = labels.empty() ? 0 : static_cast<int>(labels.size());

    int prevrow = m_rowBegin;

    for (int i = m_rowBegin; i < numlabels && m_linesPrinted <= totalrows && 
         itr != curs.end(); ++i, ++itr)
    {

        Currency cur(*itr);
        bool canpaginate = CanPaginate(cur);
        m_linesPrinted++;

        CurrencyDisplay* nesteddisplay = cur.GetDisplay();
        int cachedLines = (m_paginate == PAGINATE_INTERACTIVE) ? m_linesPrinted : 0;

        if (nesteddisplay)
        {
            nesteddisplay->SetParentDisplay(const_cast<CellNDisplay*>(this));
        }
       
        std::string tmp (cur.GetOutputString(fmt));

        m_linesPrinted += cachedLines;

        _lastRowPrinted = i;
        StripEndline(tmp);

        if (canpaginate)
        {
            if (nesteddisplay)
            {
                if (nesteddisplay->IsPaginatingCols() ||
                    nesteddisplay->IsPaginatingRows())
                {
                    // Nested display is still paginating
                    if (m_signalHandler)
                        m_signalHandler->OnAddDisplayHandler(nesteddisplay);
                    std::string output(os.str() + labels[i] + tmp);
                    if (!output.empty() && output[output.size() - 1] == '\n')
                    {
                        output.pop_back();
                    }
                    _childPaginating = true;
                   
                    m_rowEnd = i + 1;
                    return output;
                }
            }
        }

        CurrencyDisplay::DeleteDisplay(nesteddisplay);
        cur.SetDisplay(nullptr);

        os << labels[i] << tmp << '\n';

        _childPaginating = false;

        m_linesPrinted++;
        prevrow++;
        m_rowEnd = prevrow;
    }

    _childPaginating = false;
    std::string output(os.str());
    if (!output.empty() && output[output.size() - 1] == '\n')
    {
        output.pop_back();
        if (m_linesPrinted >= 1)
        {
            m_linesPrinted--;
        }
    }
    return output;
}
//------------------------------------------------------------------------------
// Gets data with back pagination
//------------------------------------------------------------------------------
std::string CellNDisplay::GetOutputBackPagination(const OutputFormat* fmt) const
{
    HML_ND_CELLARRAY* cell = m_currency.CellArrayND();
    if (!cell)
    {
        return "";
    }

    std::vector<std::string> labels;
    std::vector<Currency>    curs;
    GetSlices(cell, curs, labels);

    if (HasOnlyEmptySlices(curs))
    {
        return GetOutputEmpty();
    }

    int numrows = 0;
    int numcols = 1;
    CellNDisplay::GetCurrencySize(numrows, numcols);

    int endrow = (m_rowEnd >= 0 && m_rowEnd < numrows) ? m_rowEnd : numrows - 1;
    int endcol = 0;
    int linestofit = GetNumRowsToFit();

    std::string output;
    m_colBegin = 0;

    std::vector<Currency>::reverse_iterator ritr = curs.rbegin() + endrow;
    std::vector<Currency>::reverse_iterator end = curs.rend();
    int numlabels = labels.empty() ? 0 : static_cast<int>(labels.size());
    for (int i = endrow; 
        i >= 0 && m_linesPrinted <= m_maxRows && ritr != end && i < numlabels;
        --i, ++ritr)
    {
        m_rowBegin = i;
        std::ostringstream os;

        Currency cur(*ritr);
        bool canpaginate = CanPaginate(cur);
        int cachedLines = (m_paginate == PAGINATE_INTERACTIVE) ? m_linesPrinted : 0;

        if (canpaginate)
        {
            CurrencyDisplay* disp = cur.GetDisplay();
            if (disp)
            {
                disp->SetParentDisplay(const_cast<CellNDisplay*>(this));
            }
        }
        os << '\n' << labels[i] << " = ";
        std::string tmp(cur.GetOutputString(fmt));
        StripEndline(tmp);
        os << tmp;
        m_linesPrinted += cachedLines;

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
                    m_rowEnd = i;
                    m_linesPrinted -= 1;
                    return output;
                }
            }
        }

        CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
        cur.SetDisplay(nullptr);
    }
    //m_colEnd = numcols - 1;
    StripEndline(output);
    if (!output.empty() && output[0] == '\n')
    {
        output.erase(output.begin());
    }
    return output;
}
//------------------------------------------------------------------------------
// Sets data for back pagination
//------------------------------------------------------------------------------
void CellNDisplay::SetBackDisplayData()
{
    m_linesPrinted = 0;

    // We have reached the beginning of the cell, so we cannot go back further
    if (m_rowBegin == 0 && m_colBegin == 0)
    {
        _lastRowPrinted = -1;
        return;
    }

    int numcols = 0;
    int numrows = 0;
    GetCurrencySize(numrows, numcols);

    if (m_colBegin >= 0 && m_colBegin < numcols - 1) // Paginating cols
    {
        m_colEnd = std::max(m_colBegin - 1, 0);
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
// Updates the number of rows to fit
//------------------------------------------------------------------------------
void CellNDisplay::UpdateNumLinesPrinted()
{
    if (m_rowBegin == 0 && m_colBegin == 0 && !m_currency.IsDispOutput())
    {
        m_linesPrinted++; 
    }
}
//------------------------------------------------------------------------------
// Gets number of rows that can be fit
//------------------------------------------------------------------------------
int CellNDisplay::GetNumRowsToFit() const
{
    int numrows = m_maxRows - m_linesPrinted;
    if (numrows < 1)
        numrows = 1;  // Force print at least one row

    return numrows;
}
//------------------------------------------------------------------------------
// True if this ND matrix has only empty slices
//------------------------------------------------------------------------------
bool CellNDisplay::HasOnlyEmptySlices(const std::vector<Currency>& slices) const
{
    if (slices.empty())
    {
        return false;
    }

    for (std::vector<Currency>::const_iterator itr = slices.begin();
        itr != slices.end(); ++itr)
    {
        if (!(*itr).IsNothing())
        {
            return false;
        }
    }
    return true;
}
//------------------------------------------------------------------------------
// Gets display string if matrix ND has only empty slices
//------------------------------------------------------------------------------
std::string CellNDisplay::GetOutputEmpty() const
{
    std::string out("[Cell] ");

    HML_ND_CELLARRAY* cell = m_currency.CellArrayND();
    if (!cell)
    {
        return "";
    }
    std::vector<int> dims(cell->Dimensions());

    size_t ndims = dims.size();

    std::vector<int>::const_iterator itr = dims.begin();
    for (size_t i = 0; itr != dims.end(); ++itr, ++i)
    {
        out += std::to_string(static_cast<long long>(*itr));
        if (i != ndims - 1)
        {
            out += " x ";
        }
    }
    m_linesPrinted++;

    m_rowBegin = -1;
    m_rowEnd = -1;
    return out;
}
//------------------------------------------------------------------------------
// True if rows are being processed during pagination
//------------------------------------------------------------------------------
bool CellNDisplay::IsPaginatingRows() const
{
    int m = 0;
    if (CurrencyDisplay::IsPaginatingRows())
    {
        return true;
    }
    if (_lastRowPrinted >= 0)
    {
        int n = 0;
        GetCurrencySize(m, n);
        if (_lastRowPrinted < m - 1)
        {
            return true;
        }
        return false;
    }
    if (_childPaginating)
    {
        if (m_rowBegin >= m)
        {
            _childPaginating = false;
            return false;
        }
        return true;
    }
        
    return false;
}
