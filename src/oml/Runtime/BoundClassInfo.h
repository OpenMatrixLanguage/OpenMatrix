/**
* @file BoundClassInfo.h
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

#ifndef __BOUNDCLASSINFO_H__
#define __BOUNDCLASSINFO_H__

#include "EvaluatorInt.h"

#include <map>

class BoundProperty;

//------------------------------------------------------------------------------
//!
//! \brief Interface for internal Swig bound C++ classes
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS BoundClassInfo
{
public:
    //!
    //! Constructor
    //!
    BoundClassInfo();
    //!
    //! Destructor
    //!
    ~BoundClassInfo();

    //!
    //! Adds method
    //! \param name Method name
    //! \param fptr Function pointer to method
    void AddMethod( const std::string& name,
                    FUNCPTR            fptr);
    //!
    //! Gets method
    //! \param name Method name
    //!
    FUNCPTR GetMethod( const std::string& name) const;
    //!
    //! Searches methods, and property setters (if needed) and returns funcptr
    //! \param name Given name to match
    //!
    bool HasMethod( const std::string& name) const;
    //!
    //! Populates given vector with method names
    //! \param names Vector of method names
    //!
    void GetMethodNames( std::vector<std::string>& names) const;
    //!
    //! Adds a property
    //! \param name   Property name
    //! \param getter Function pointer to method that gets property
    //! \param setter Function pointer to method that sets property
    //!
    void AddProperty( const std::string& name,
                      FUNCPTR            getter,
                      FUNCPTR            setter);
    //!
    //! Returns true if a property with the given name exists
    //! \param name Name of the property
    //!
    bool HasProperty( const std::string& name) const;
    //!
    //! Gets property setter
    //! \param name Name of the property
    //!
    FUNCPTR GetPropertySetter( const std::string& name) const;

private:
    std::map<std::string, FUNCPTR>* _methods;    //!<  name-func: Methods
    std::vector<BoundProperty*>*    _properties; //!<  Properties

    //!
    //! Gets property
    //! \param name Property name
    //!
    BoundProperty* GetProperty( const std::string& name) const;
};
//------------------------------------------------------------------------------
//!
//! Interface for storing property information in (swig) bound classes
//!
//------------------------------------------------------------------------------
class HML2DLL_DECLS BoundProperty
{
    friend class BoundClassInfo;  //!< Special access to construct
public:
    //!
    //! Destructor
    //!
    ~BoundProperty() {}

private:
    std::string _name;            //!< Name
    FUNCPTR     _getter;          //!< Function pointer to getter
    FUNCPTR     _setter;          //!< Function pointer to setter

    //!
    //! Constructor
    //! \param name   Name of the property
    //! \param getter Function pointer to get value of property
    //! \param setter Function pointer to set value of property
    BoundProperty( const std::string& name,
                   FUNCPTR            getter,
                   FUNCPTR            setter)
                   : _name(name), _getter(getter), _setter(setter) {}
};
#endif // __BOUNDCLASSINFO_H__

