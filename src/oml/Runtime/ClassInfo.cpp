/**
* @file ClassInfo.cxx
* @date March 2016
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
#include "ClassInfo.h"

#include <cassert>

#include "FunctionInfo.h"
#include "MemoryScope.h"
#include "StructData.h"
// End defines/includes

//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
ClassInfo::~ClassInfo()
{
    for (std::vector<PropertyInfo*>::iterator itr = _properties.begin();
         itr != _properties.end();)
    {
        delete *itr;
        itr = _properties.erase(itr);
    }
}
//------------------------------------------------------------------------------
//! Gets the function info corresponding to the method name in the language
//! \param[in] name Name of the method in the language
//------------------------------------------------------------------------------
FunctionInfo* ClassInfo::GetFunctionInfo(const std::string& name) const
{
    if (name.empty()) 
		return NULL;

    std::map<std::string, FunctionInfo*>::const_iterator itr = _methods.find(name);
    if (itr != _methods.end()) 
        return itr->second;

	for (int j=0; j<_baseclass.size(); ++j)
	{
		FunctionInfo* fi = _baseclass[j]->GetFunctionInfo(name);
		if (fi)
			return fi;
	}

    return NULL;
}
//------------------------------------------------------------------------------
//! Adds a base class
//! \param[in] name Given base class name
//------------------------------------------------------------------------------
void ClassInfo::AddBaseClass(const ClassInfo* base_class)
{
    if (!base_class) return;

    if (_baseclass.empty() || 
        std::find(_baseclass.begin(), _baseclass.end(), base_class) == _baseclass.end())
        _baseclass.push_back(base_class);
}
//------------------------------------------------------------------------------
//! Returns true if given class name is a base class - used only in language
//! \param[in] baseclass Given base class name
//------------------------------------------------------------------------------
bool ClassInfo::IsSubclassOf(const std::string& baseclass) const
{
    if (baseclass.empty()) return false;

	for (int j=0; j<_baseclass.size(); j++)
	{
		const ClassInfo* base_class = _baseclass[j];
		if (base_class->_class_name == baseclass)
			return true;
	}

    return false;
}
//------------------------------------------------------------------------------
//! Adds class method
//! \param[in] name Name of the method in the language
//! \param[in] fi   Function info
//------------------------------------------------------------------------------
void ClassInfo::AddClassMethod(const std::string& name, FunctionInfo* fi)
{
    if (!name.empty())
        _methods[name] = fi;
}
//------------------------------------------------------------------------------
//! Returns true if given function pointer is a method - used only in language
//! \param[in] fi Given function pointer
//------------------------------------------------------------------------------
bool ClassInfo::IsClassMethod(FunctionInfo* fi) const
{
    if (!fi) return false;

	if (fi->IsConstructor())
	{
		if (fi->FunctionName() == _class_name)
			return true;
	}

	std::map<std::string, FunctionInfo*>::const_iterator iter = _methods.begin();
	for (; iter != _methods.end(); ++iter)
	{
		if (iter->second == fi)
			return true;
	}

	for (int j=0; j<_baseclass.size(); ++j)
	{
		if (_baseclass[j]->IsClassMethod(fi))
			return true;
	}

	return false;
}
//------------------------------------------------------------------------------
//! Registers a property
//! \param[in] name      Property name
//! \param[in] isPrivate True if property is private, defaults to false
//------------------------------------------------------------------------------
void ClassInfo::AddProperty(const std::string& name, bool isPrivate)
{
    if (name.empty() || HasProperty(name)) return;

    PropertyInfo* prop = new PropertyInfo(name, isPrivate);    
    _properties.push_back(prop);
}
//------------------------------------------------------------------------------
//! Returns true if a property with the given name exists
//! \param[in] name Name of the property
//------------------------------------------------------------------------------
bool ClassInfo::HasProperty(const std::string& name) const
{
    if (name.empty()) return false;

    for (std::vector<PropertyInfo*>::const_iterator itr = _properties.begin();
         itr != _properties.end(); ++itr)
    {
        if ((*itr)->Name() == name) 
            return true;
    }

	for (int j=0; j<_baseclass.size(); ++j)
	{
		if (_baseclass[j]->HasProperty(name))
			return true;
	}

    return false;
}
//------------------------------------------------------------------------------
//! Gets property with the given name
//! \param[in] name Name of the property
//------------------------------------------------------------------------------
PropertyInfo* ClassInfo::GetProperty(const std::string& name) const
{
    if (name.empty()) return NULL;

    for (std::vector<PropertyInfo*>::const_iterator itr = _properties.begin();
         itr != _properties.end(); ++itr)
    {
        if ((*itr)->Name() == name) 
            return (*itr);
    }
    return NULL;
}
//------------------------------------------------------------------------------
//! Returns true if property with the given name is private
//! \param[in] name Name of the property
//------------------------------------------------------------------------------
bool ClassInfo::IsPropertyPrivate(const std::string& name) const
{
    PropertyInfo* prop = GetProperty(name);
    return (prop && prop->IsPrivate());
}
// end of file:

Currency ClassInfo::CreateEmpty() const
{
	Currency cur;
	cur.MakeStruct();

	if (_baseclass.size())
	{
		// not allowing multiple inheritance for now
		cur = _baseclass[0]->CreateEmpty();
	}

	StructData* sd = cur.Struct();

	for (int j=0; j<_properties.size(); j++)
		sd->SetValue(0, 0, _properties[j]->Name(), new hwMatrix); 

	cur.SetClass(_class_name);

	return cur;
}