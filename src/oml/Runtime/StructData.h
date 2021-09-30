/**
* @file StructData.h
* @date September 2013
* Copyright (C) 2013-2018 Altair Engineering, Inc.  
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

#ifndef __StructData_h
#define __StructData_h

#include <string>
#include <vector>
#include <hwMatrix.h>
#include "Currency.h"
#include <map>

typedef hwTMatrix<std::map<int, Currency>, void*> HML_FIELDVALS;

class OMLDLL_DECLS StructData 
{
public:
	StructData():field_values(NULL), next_key(1) {}
	~StructData() { delete field_values; }
	StructData(const StructData&);

	const Currency& GetValue(int index_1, int index_2, const std::string& field) const;
	const Currency& GetValue(int index_1, int index_2, const std::string* field) const;
	void            SetValue(int index_1, int index_2, const std::string& field, Currency value);
	void            SetValue(int index_1, int index_2, const std::string* field, Currency value);
	const Currency* GetPointer(int index_1, int index_2, const std::string& field) const;
	const Currency* GetPointer(int index_1, int index_2, const std::string* field) const;
	StructData*     GetElement(int index_1, int index_2);
	void            SetElement(int index_1, int index_2, StructData* sd);
	void            Dimension(int index_1, int index_2, bool force = false);
	void            DimensionNew(int index_1);
	void            DimensionNew(int index_1, int index_2);
	bool            Contains(const std::string& field) const;
	bool            Contains(const std::string* field) const;
	bool            IsEmpty() const;
    inline int      M() const { return field_values == nullptr ? 0 : field_values->M(); }
    inline int      N() const { return field_values == nullptr ? 0 : field_values->N(); }
	void            addField(const std::string*);
    void            addField(const std::string&);
	void            removeField(const std::string*);
	void            removeField(const std::string&);
	void            IncrRefCount();
	void            DecrRefCount();
	int             GetRefCount() const;

	void            Transpose();

	const std::map<std::string, int> GetFieldNames() const;
    const std::map<const std::string*, int>& GetFieldNamePtrs() const { return field_names; }
    hwMathStatus Reshape (int m, int n) { return field_values->Reshape(m, n); }

    //! Returns the size
    inline int Size() const { return field_values == nullptr ? 0 : field_values->Size(); }
private:

	std::map<const std::string*, int> field_names;
	HML_FIELDVALS*             field_values;
	int                        next_key;
};

#endif