/**
* @file CellNDisplay.cxx
* @date January 2019
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
//#include "CellND.cc"

#include <cassert>
#include <memory>        // For std::unique_ptr

#include "CellNDisplay.h"
#include "EvaluatorInt.h"

//------------------------------------------------------------------------------
// Constructor - Only currency or derived classes can access constructor
//------------------------------------------------------------------------------
CellNDisplay::CellNDisplay(const Currency& cur)
	: CurrencyDisplay(cur)
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

	std::string header;
	if (m_rowBegin == 0 && !m_currency.IsDispOutput())
	{
		header = os.str();
	}

	// Disabled pagination for now
	std::string output = GetOutputNoPagination(fmt);
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
		os << temp.GetOutputString(fmt) << std::endl;
	}
	else
	{
		size_t numlabels = labels.empty() ? 0 : labels.size();
		std::vector<Currency>::const_iterator itr = curs.begin();

		for (size_t i = 0; itr != curs.end() && i < numlabels; ++itr, ++i)
		{
			Currency cur(*itr);

			os << std::endl << labels[i] << " = ";

			std::string tmp(cur.GetOutputString(fmt));
			StripEndline(tmp);
			if (!tmp.empty())
			{
				os << tmp << std::endl;
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

		os << std::endl << labels[i] << " = " << std::endl;
		os << cur.GetValues(fmt) << std::endl;

		CurrencyDisplay::DeleteDisplay(cur.GetDisplay());
		cur.SetDisplay(nullptr);
	}

	return os.str();
}
