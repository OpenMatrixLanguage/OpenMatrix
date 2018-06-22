/**
* @file StructData.cpp
* @date September 2013
* Copyright (C) 2013-2018 Altair Engineering, Inc.  
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

#include "StructData.h"
#include "ErrorInfo.h"
#include "OML_Error.h"

StructData::StructData(const StructData& in)
{
	field_names  = in.field_names;

	if (in.field_values)
		field_values = new HML_FIELDVALS(*in.field_values);	
	else
		field_values = NULL;

	next_key     = in.next_key;
}

const Currency& StructData::GetValue(int index_1, int index_2, std::string field) const
{
	static Currency _not_used = new hwMatrix();

	const Currency* ret = GetPointer(index_1, index_2, field);

	if (!ret)
		return _not_used;

	return *ret;
}

const Currency* StructData::GetPointer(int index_1, int index_2, std::string field) const
{
	std::map<std::string, int>::const_iterator temp;
	temp = field_names.find(field);

	if (temp == field_names.end())
		return NULL;

	int field_key = temp->second;

	if (index_2 == -1)
	{
		if (index_1 < 0)
			throw OML_Error(HW_ERROR_INDEXRANGE);
		else if (index_1 >= field_values->Size())
			throw OML_Error(HW_ERROR_INDEXRANGE);
	}
	else
	{
		if (index_1 < 0)
			throw OML_Error(HW_ERROR_INDEXRANGE);
		else if (index_1 >= field_values->M())
			throw OML_Error(HW_ERROR_INDEXRANGE);

		if (index_2 < 0)
			throw OML_Error(HW_ERROR_INDEXRANGE);
		else if (index_2 >= field_values->N())
			throw OML_Error(HW_ERROR_INDEXRANGE);
	}

	std::map<int, Currency>* field_vals;
	
	if (index_2 != -1)
		field_vals = &((*field_values)(index_1, index_2));
	else
		field_vals = &((*field_values)(index_1));

	std::map<int, Currency>::const_iterator temp2;
	temp2 = field_vals->find(field_key);

	if (temp2 == field_vals->end())
		return NULL;

	const Currency* ret_val = &(temp2->second);
	return ret_val;
}

void StructData::SetValue(int index_1, int index_2, std::string field, Currency value)
{
	value.ClearOutputName();

	if (!field_values)
	{
		if (index_2 == -1)
			DimensionNew(index_1+1);
		else
			DimensionNew(index_1+1, index_2+1);
	}
	else if (index_2 == -1)
	{
		if ((index_1 >= field_values->Size()))
			DimensionNew(index_1+1);
	}
	else if ((index_1 >= field_values->M()) || (index_2 >= field_values->N()))
	{
		int new_m = field_values->M()-1;
		if (index_1 > new_m)
			new_m = index_1;

		int new_n = field_values->N()-1;
		if (index_2 > new_n)
			new_n = index_2;

		DimensionNew(new_m+1, new_n+1);
	}

	if (field_names.find(field) == field_names.end())
	{
		field_names[field] = next_key;
		next_key++;
	}

	int field_key = field_names[field];
	
	if (index_2 == -1)
	{
		if (index_1 < 0)
			throw OML_Error(HW_ERROR_INDEXRANGE);

		std::map<int, Currency>& field_vals = (*field_values)(index_1);
		field_vals[field_key] = value;
	}
	else
	{
		if ((index_1 < 0) || (index_2 < 0))
			throw OML_Error(HW_ERROR_INDEXRANGE);

		std::map<int, Currency>& field_vals = (*field_values)(index_1, index_2);
		field_vals[field_key] = value;
	}
}

void StructData::Dimension(int num_rows, int num_cols, bool force)
{
	if (!field_values)
		field_values = new HML_FIELDVALS;

	if (field_values->IsEmpty())
	{
		if ((num_cols == -1) && !force)
			field_values->Dimension(num_rows+1, HML_FIELDVALS::REAL);
		else
			field_values->Dimension(num_rows+1, num_cols+1, HML_FIELDVALS::REAL);
	}
	else
	{
		if ((num_cols == -1) && (num_rows != -1))
			field_values->Resize(1, num_rows+1);
		else
			field_values->Resize(num_rows+1, num_cols+1);
	}
}

void StructData::DimensionNew(int num_rows, int num_cols)
{
	if (!field_values)
		field_values = new HML_FIELDVALS;

	if (field_values->IsEmpty())
		field_values->Dimension(num_rows, num_cols, HML_FIELDVALS::REAL);
	else
		field_values->Resize(num_rows, num_cols);
}

void StructData::DimensionNew(int num_cols)
{
	DimensionNew(1, num_cols);
}

StructData* StructData::GetElement(int index_1, int index_2)
{
	StructData* ret_val = new StructData;
	ret_val->DimensionNew(1, 1);
	ret_val->field_names = field_names;

	if (index_2 != -1)
	{
		if ((index_1 > field_values->M()) || (index_1 < 1))
			throw OML_Error(HW_ERROR_INDEXRANGE);

		if ((index_2 > field_values->N()) || (index_2 < 1))
			throw OML_Error(HW_ERROR_INDEXRANGE);

		(*ret_val->field_values)(0) = (*field_values)(index_1-1, index_2-1);
	}
	else
	{
		if ((index_1 > field_values->Size()) || (index_1 < 1))
			throw OML_Error(HW_ERROR_INDEXRANGE);

		(*ret_val->field_values)(0) = (*field_values)(index_1-1);
	}

	return ret_val;
}

void StructData::SetElement(int index_1, int index_2, StructData* sd)
{
	if (field_names.size() == 0)
		field_names = sd->field_names;

	if (index_2 != -1)
	{
		if (index_1 > field_values->M())
			throw OML_Error(HW_ERROR_INDEXRANGE);

		if (index_2 > field_values->N())
			throw OML_Error(HW_ERROR_INDEXRANGE);

		(*field_values)(index_1-1, index_2-1) = (*sd->field_values)(0);
	}
	else
	{
		if (index_1 > field_values->Size())
			throw OML_Error(HW_ERROR_INDEXRANGE);

		(*field_values)(index_1-1) = (*sd->field_values)(0);
	}
}

void StructData::addField(std::string& field)
{
	if (field_names.find(field) == field_names.end())
		field_names[field] = next_key++;
}

bool StructData::Contains(const std::string& field) const
{
	std::map<std::string, int>::const_iterator temp; temp = field_names.find(field);

	if (temp == field_names.end())
		return false;

	return true;
}

void StructData::removeField(const std::string& fieldname)
{
	if (Contains(fieldname))
	{
		int id = field_names[fieldname];
		field_names.erase(fieldname);
		for (int i = 0; i < field_values->Size(); ++i)
			(*field_values)(i).erase(id);
	}
}

bool StructData::IsEmpty() const
{
	// only considered empty if there are no elements OR
	// there is one element that has no data
	if (field_values->Size() == 0)
		return true;

	if (field_values->Size() == 1)
	{
		if ((*field_values)(0).empty())
			return true;
	}

	return false;
}

void StructData::IncrRefCount()
{
	if (!field_values)
		DimensionNew(0, 0);

	field_values->IncrRefCount();
}

void StructData::DecrRefCount()
{
	field_values->DecrRefCount();
}

int StructData::GetRefCount() const
{
	if (field_values)
		return field_values->GetRefCount();
	else
		return 1;
}