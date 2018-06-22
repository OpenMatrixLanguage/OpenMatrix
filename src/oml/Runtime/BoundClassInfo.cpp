/**
* @file BoundClassInfo.cpp
* @date October 2016
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
#include "BoundClassInfo.h"

#include <cassert>

#include "OML_Error.h"
// End defines/includes

//------------------------------------------------------------------------------
//! Constructor
//------------------------------------------------------------------------------
BoundClassInfo::BoundClassInfo() 
    : _methods (NULL)
    , _properties (NULL)
{
    _methods    = new std::map<std::string, FUNCPTR>; // To avoid warnings
    _properties = new std::vector<BoundProperty*>;
}
//------------------------------------------------------------------------------
//! Destructor
//------------------------------------------------------------------------------
BoundClassInfo::~BoundClassInfo() 
{
    for (std::vector<BoundProperty*>::iterator itr = _properties->begin();
         itr != _properties->end();)
    {
        delete *itr;
        itr = _properties->erase(itr);
    }
    delete _properties;
    delete _methods;
}
//------------------------------------------------------------------------------
//! Adds method
//------------------------------------------------------------------------------
void BoundClassInfo::AddMethod(const std::string& name, FUNCPTR funcptr)
{
    if (!name.empty())
        (*_methods)[name] = funcptr;
}
//------------------------------------------------------------------------------
//! Returns true if a method with the given name exists
//------------------------------------------------------------------------------
bool BoundClassInfo::HasMethod(const std::string& name) const
{
    return (_methods->find(name) != _methods->end());
}
//------------------------------------------------------------------------------
//! Searches methods, and property setters (if needed) and returns funcptr
//------------------------------------------------------------------------------
FUNCPTR BoundClassInfo::GetMethod(const std::string& name) const
{
    if (name.empty()) return nullptr;

    // Search through methods
    std::map<std::string, FUNCPTR>::const_iterator itr1 = _methods->find(name);
    if (itr1 != _methods->end()) 
        return itr1->second;

    // Search through properties
    BoundProperty* prop = GetProperty(name);
    if (prop)
        return prop->_getter;

    throw OML_Error("Error: Unknown method: " + name);
    return nullptr;
}
//------------------------------------------------------------------------------
//! Adds a property
//------------------------------------------------------------------------------
void BoundClassInfo::AddProperty(const std::string& name,
                                 FUNCPTR            getter,
                                 FUNCPTR            setter)
{
    if (name.empty() || HasProperty(name)) return;

    BoundProperty* prop = new BoundProperty(name, getter, setter);   
    _properties->push_back(prop);
}
//------------------------------------------------------------------------------
//! Returns true if a property with the given name exists
//------------------------------------------------------------------------------
bool BoundClassInfo::HasProperty(const std::string& name) const
{
    if (name.empty()) return false;

    for (std::vector<BoundProperty*>::const_iterator itr = _properties->begin();
         itr != _properties->end(); ++itr)
    {
        if ((*itr)->_name == name) 
            return true;
    }
    return false;
}
//------------------------------------------------------------------------------
//! Gets property setter
//------------------------------------------------------------------------------
FUNCPTR BoundClassInfo::GetPropertySetter(const std::string& name) const
{
    BoundProperty* prop = GetProperty(name);
    if (prop) 
        return prop->_setter;

    return nullptr;
}
//------------------------------------------------------------------------------
//! Gets property
//------------------------------------------------------------------------------
BoundProperty* BoundClassInfo::GetProperty(const std::string& name) const
{
    if (name.empty()) return NULL;

    for (std::vector<BoundProperty*>::const_iterator itr = _properties->begin();
        itr != _properties->end(); ++itr)
    {
        BoundProperty* prop = *itr;
        assert(prop);
        if (prop->_name == name) 
            return prop;
    }
    return NULL;
}
//------------------------------------------------------------------------------
//! Populates given vector with method names
//------------------------------------------------------------------------------
void BoundClassInfo::GetMethodNames( std::vector<std::string>& names) const
{
    for (std::map<std::string, FUNCPTR>::const_iterator itr = _methods->begin();
         itr != _methods->end(); ++itr)
    {
        std::string str (itr->first);
        if (std::find(names.begin(), names.end(), str) == names.end())
            names.push_back(str);
    }
}
