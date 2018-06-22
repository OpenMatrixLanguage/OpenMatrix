/**
* @file MatrixNDisplay.cpp
* @date July 2016
* Copyright (C) 2016-2018 Altair Engineering, Inc.  
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

// Begin defines/includes
#include "MatrixNDisplay.h"

#include <cassert>

#include "Currency.h"
#include "Interpreter.h"
#include "OutputFormat.h"
#include "SignalHandlerBase.h"

// End defines/includes
//------------------------------------------------------------------------------
// Constructor - Only currency or derived classes can access constructor
//------------------------------------------------------------------------------
MatrixNDisplay::MatrixNDisplay(const Currency& cur)
    : CurrencyDisplay(cur)
{
}
//------------------------------------------------------------------------------
// Gets number of rows and cols in currency
//------------------------------------------------------------------------------
void MatrixNDisplay::GetCurrencySize(int& rows, int& cols) const
{
    const hwMatrixN* mtx = m_currency.MatrixN();
    if (!mtx) return;

    std::vector<int> dims = mtx->Dimensions();
    if (dims.empty()) return;
    size_t numdims = dims.size();      

    cols = 1;
    rows = 1;
      
    if (numdims == 3)
        rows = dims[2];        // 3D matrix
    else if (numdims > 3)      // ND matrix
    {
        for (std::vector<int>::const_iterator itr = dims.begin()+2; 
             itr != dims.end(); ++itr)
             rows *= *itr;
    }
}
//------------------------------------------------------------------------------
// Gets output - called from Currency::GetOutputString
//------------------------------------------------------------------------------
std::string MatrixNDisplay::GetOutput(const OutputFormat* fmt, 
                                      std::ostringstream& os) const
{
    if (m_mode == DISPLAYMODE_QUIT || m_mode == DISPLAYMODE_EXIT)    
		return "";

    int numrows = 0;
    int numcols = 0;
    GetCurrencySize(numrows, numcols);

    if (numrows == 0 || numcols == 0)
    {
        if (m_currency.IsDispOutput()) 
            return "";
        return std::string(os.str());
    }

	bool paginate = IsValidDisplaySize();
    if (!paginate)
    {
        m_rowBegin = 0;
        m_rowEnd   = 0;
        m_colBegin = 0;
        m_colEnd   = 0;
    }
    
    std::string header;
    if (m_rowBegin == 0 && !m_currency.IsDispOutput())
        header = os.str();

    std::string output;
	if (!paginate) 
        output = GetOutputNoPagination(fmt);
    else
    {
        const_cast<MatrixNDisplay*>(this)->Initialize(fmt, NULL, NULL); 
        //m_linesPrinted ++;       // Matrix header will be printed

        if (m_mode == DISPLAYMODE_FORWARD || m_mode == DISPLAYMODE_DOWN)
            output = GetOutputForwardPagination(fmt);

        else if (m_mode == DISPLAYMODE_BACK || m_mode == DISPLAYMODE_UP)
            output = GetOutputBackPagination(fmt);
    }

    return (header + output);
}
//------------------------------------------------------------------------------
// Gets output with no pagination
//------------------------------------------------------------------------------
std::string MatrixNDisplay::GetOutputNoPagination(const OutputFormat* fmt) const
{
    const hwMatrixN* mtx = m_currency.MatrixN();
    assert(mtx);
    
    std::vector<std::string> labels;
    std::vector<Currency>    curs;
    GetSlices(mtx, curs, labels);
    
    std::ostringstream os;

	if (mtx->Size() == 1)
	{
		Currency temp((*mtx)(0));
        os << temp.GetOutputString(fmt) << std::endl;
	}
	else
	{
		size_t numlabels = labels.empty() ? 0 : labels.size();   
		std::vector<Currency>::const_iterator itr = curs.begin();

		for (size_t i = 0; itr != curs.end() && i < numlabels; ++itr, ++i)
		{	
			Currency cur (*itr);

			os << std::endl << labels[i] << " = "  << std::endl;

            std::string tmp (cur.GetOutputString(fmt));
            StripEndline(tmp);
            if (!tmp.empty())
            {
			    os << tmp << std::endl;
            }

			CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
			cur.SetDisplay(NULL);
		}
	}

    std::string tmp (os.str());
    StripEndline(tmp);

	return tmp;
}
//------------------------------------------------------------------------------
// Helper to slice a given ND matrix to 2D matrices
//------------------------------------------------------------------------------
void MatrixNDisplay::GetSlicesHelper(const hwMatrixN*          mtx,
                                     std::vector<hwSliceArg>   slices,
                                     std::vector<Currency>&    curs,
                                     std::vector<std::string>& labels)
{
    if (!mtx) return;

    std::vector<int> dims    = mtx->Dimensions();
	size_t           numdims = dims.empty() ? 0 : dims.size();

    if (numdims < 3)
    {
        // Get the slice label
        std::string label ("slice(");
        for (std::vector<hwSliceArg>::const_iterator itr = slices.begin();
             itr != slices.end(); ++itr)
		{
			hwSliceArg slice = *itr;

            if (itr != slices.begin())
                label += ", ";

			if (slice.IsColon())
				label += ":";
			else
                label += std::to_string(static_cast<long long>(slice.Scalar() + 1));
		}
        label += ")";

        labels.push_back(label);                       // Slice label for matrix
        curs.push_back(Currency(ConvertNDto2D(mtx)));  // Matrix type currency

        return;
    }

	size_t curDimIndex = slices.size();
	if (curDimIndex == numdims)
	{
		hwMatrixN tmpslice;
		mtx->SliceRHS(slices, tmpslice);
		
        GetSlicesHelper(&tmpslice, slices, curs, labels);
	}
	else if (curDimIndex < numdims)
	{
		int curdim = dims[curDimIndex];

		for (int j = 0; j < curdim; ++j)
		{
			slices.push_back(j);
			MatrixNDisplay::GetSlicesHelper(mtx, slices, curs, labels);
			slices.pop_back();
		}
	}
}
//------------------------------------------------------------------------------
// Helper to slice a given ND matrix to 2D matrices
//------------------------------------------------------------------------------
void MatrixNDisplay::GetSlices(const hwMatrixN*          mtx,
                               std::vector<Currency>&    curs,
                               std::vector<std::string>& labels)
{
    if (!mtx) return;

    std::vector<hwSliceArg> base;
    base.push_back(hwSliceArg());
    base.push_back(hwSliceArg());

    GetSlicesHelper(mtx, base, curs, labels);

    assert(curs.size() == labels.size());
}
//------------------------------------------------------------------------------
// Sets data for forward display
//------------------------------------------------------------------------------
void MatrixNDisplay::SetForwardDisplayData()
{
    int numrows = 0;
    int numcols = 0;
    GetCurrencySize(numrows, numcols);

    m_colBegin = 0; 

	if (IsPaginatingRows())               // Row pagination
		m_rowBegin = m_rowEnd + 1;
	
    else if (m_rowEnd >= numrows)
    {
        m_colEnd   = -1;
        m_rowEnd   = -1;
        m_rowBegin = -1;
        m_colBegin = -1;
        return;
    }

    // Make sure indices are positive
    m_rowBegin = std::max(m_rowBegin, 0);

    // Reset the end indices
    m_colEnd = -1;
    m_rowEnd = -1;

    // Figure out how many lines have been printed and set it
    if (!m_parentDisplay ||                // First currency being paginated
       (m_rowBegin > 0 || m_colBegin > 0)) // Only this currency is being paginated
    {
        m_linesPrinted = 0;
    }
}
//------------------------------------------------------------------------------
// Gets data with forward pagination
//------------------------------------------------------------------------------
std::string MatrixNDisplay::GetOutputForwardPagination(const OutputFormat* fmt) const
{
    const hwMatrixN* mtx = m_currency.MatrixN();
    assert(mtx);
    
    int numrows = 0;
    int numcols = 1;
    MatrixNDisplay::GetCurrencySize(numrows, numcols);

    int linestofit = GetNumRowsToFit();
    int totalrows  = m_linesPrinted + linestofit;
    m_colEnd = 0;

    std::vector<std::string> labels;
    std::vector<Currency>    curs;
    GetSlices(mtx, curs, labels);

    std::ostringstream os;
    int                numlabels = labels.empty() ? 0 : static_cast<int>(labels.size());  

    std::vector<Currency>::iterator itr = curs.begin() + m_rowBegin;
	for (int i = m_rowBegin; 
         i < numrows && i < numlabels && m_linesPrinted < totalrows && itr != curs.end(); 
         ++i, ++itr)
	{	
        m_rowEnd = i;

        Currency cur (*itr);
        bool canpaginate = CanPaginate(cur);

        os << std::endl << labels[i] << " = " << std::endl;
        m_linesPrinted += 2;

        if (canpaginate)
        {
            CurrencyDisplay* disp = cur.GetDisplay();
            if (disp)
            {
                disp->SetParentDisplay(const_cast<MatrixNDisplay*>(this));
            }
        }

        std::string tmp (cur.GetOutputString(fmt));
        StripEndline(tmp);
        os << tmp;

        if (canpaginate)
        {
            CurrencyDisplay* nesteddisplay = cur.GetDisplay();
            if (nesteddisplay)
            {
                if (nesteddisplay->IsPaginatingCols() ||
                    nesteddisplay->IsPaginatingRows())
                {
                    // Nested display is still paginating
                    nesteddisplay->SetParentDisplay(const_cast<MatrixNDisplay*>(this));
                    if (m_signalHandler)
                        m_signalHandler->OnAddDisplayHandler(nesteddisplay);
                    std::string output (os.str());
                    m_rowBegin = i;
                    if (!output.empty() && output[output.size() - 1] == '\n')
                    {
                        output.pop_back();
                    }
                    return output;
                }
            }
        }
        CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
        cur.SetDisplay(NULL);

        os << std::endl;
        m_linesPrinted++;
    }

    std::string output (os.str());
    if (!output.empty() && output[output.size()-1] == '\n')
    {
        output.pop_back();
    }
    return output;
}
//------------------------------------------------------------------------------
// Gets data with back pagination
//------------------------------------------------------------------------------
std::string MatrixNDisplay::GetOutputBackPagination(const OutputFormat* fmt) const
{
    const hwMatrixN* mtx = m_currency.MatrixN();
    assert(mtx);
    
    int numrows = 0;
    int numcols = 1;
    MatrixNDisplay::GetCurrencySize(numrows, numcols);

	int endrow     = (m_rowEnd >= 0 && m_rowEnd < numrows) ? m_rowEnd : numrows - 1;
    int endcol     = 0;
    int linestofit = GetNumRowsToFit();

    std::string output;
    m_colBegin = 0;
    
    std::vector<Currency>    curs;
    std::vector<std::string> labels;
    GetSlices(mtx, curs, labels);

    std::vector<Currency>::reverse_iterator ritr = curs.rbegin() + endrow;
    std::vector<Currency>::reverse_iterator end  = curs.rend();
    int numlabels = labels.empty() ? 0 : static_cast<int>(labels.size()); 
	for (int i = endrow; i >= 0 && m_linesPrinted <= linestofit && ritr != end; 
         --i, ++ritr)
	{	
        m_rowBegin = i;
        std::ostringstream os;

        Currency cur (*ritr);
        bool canpaginate = CanPaginate(cur);
        if (canpaginate)
        {
            CurrencyDisplay* disp = cur.GetDisplay();
            if (disp)
            {
                disp->SetParentDisplay(const_cast<MatrixNDisplay*>(this));
            }
        }

        os << std::endl << labels[i] << " = " << std::endl;
        m_linesPrinted += 2;

        std::string tmp (cur.GetOutputString(fmt));
        StripEndline(tmp);
        os << tmp;

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
        else
            m_linesPrinted++;

        CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
        cur.SetDisplay(NULL);
	}

    StripEndline(output);
    return output;
}
//------------------------------------------------------------------------------
// Sets data for back pagination
//------------------------------------------------------------------------------
void MatrixNDisplay::SetBackDisplayData()
{
    m_linesPrinted = 0;

	// We have reached the beginning of the matrix, so we cannot go back further
	if (m_rowBegin == 0) return;

    int numcols = 1;
    int numrows = 0;
    GetCurrencySize(numrows, numcols);

    m_rowEnd = std::max(m_rowBegin - 1, 0);
    m_colEnd = 0;
}
//------------------------------------------------------------------------------
// Updates the number of rows to fit
//------------------------------------------------------------------------------
void MatrixNDisplay::UpdateNumLinesPrinted()
{
    if (m_rowBegin == 0 && m_colBegin == 0 && !m_currency.IsDispOutput())
        m_linesPrinted ++;       // Assignment
}
//------------------------------------------------------------------------------
// Gets number of rows that can be fit
//------------------------------------------------------------------------------
int MatrixNDisplay::GetNumRowsToFit() const
{
    int numrows = m_maxRows - m_linesPrinted;
    if (numrows < 1)
        numrows = 1;  // Force print at least one row

    return numrows;
}
//------------------------------------------------------------------------------
// Gets values as a string
//------------------------------------------------------------------------------
std::string MatrixNDisplay::GetValues(const OutputFormat* fmt) const
{
    const hwMatrixN* mtx = m_currency.MatrixN();
    assert(mtx);
    
    std::vector<std::string> labels;
    std::vector<Currency>    curs;
    GetSlices(mtx, curs, labels);
    
	if (mtx->Size() == 1)
	{
		Currency temp((*mtx)(0));
        return temp.GetValues(fmt);
	}

    std::ostringstream os;
	size_t numlabels = labels.empty() ? 0 : labels.size();   
	std::vector<Currency>::const_iterator itr = curs.begin();

	for (size_t i = 0; itr != curs.end() && i < numlabels; ++itr, ++i)
	{	
		Currency cur (*itr);

		os << std::endl << labels[i] << " = "  << std::endl;
		os << cur.GetValues(fmt) << std::endl;

		CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
		cur.SetDisplay(NULL);
	}

	return os.str();
}
// End of file:
