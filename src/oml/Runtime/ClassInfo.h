/**
* @file ClassInfo.h
* @date March 2016
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

#ifndef __ClassInfo_h
#define __ClassInfo_h

// Begin defines/includes
#include <string>
#include <vector>

#include "Evaluator.h" // indirectly gets HML2DLL_DECLS

class FunctionInfo;
class PropertyInfo;
// End defines/includes

//------------------------------------------------------------------------------
//! Registers external C++ classes defined in OML language with OML core
//------------------------------------------------------------------------------
class HML2DLL_DECLS ClassInfo 
{
public:
    //! Constructor
    //! \param[in] name Name of the external C++ class
    ClassInfo( const std::string& name) : _class_name(name) {}
    //! Destructor
    ~ClassInfo();

    // Methods used in OML language - not bound

    //! Gets the function info corresponding to the method name in the language
    //! \param[in] name Name of the method in the language
    FunctionInfo* GetFunctionInfo( const std::string& name) const;
    //! Adds class method
    //! \param[in] name Name of the method in the language
    //! \param[in] fi   Function info
    void AddClassMethod( const std::string& name,
                         FunctionInfo*      fi);
    //! Returns true if given function pointer is a method - used only in language
    //! \param[in] fi Given function pointer
	bool IsClassMethod( FunctionInfo* fi) const;

    //! Adds a base class
    //! \param[in] name Given base class name
    void AddBaseClass( const std::string& name);
    //! Returns true if given class name is a base class - used only in language
    //! \param[in] name Given base class name
	bool IsSubclassOf( const std::string& name) const;

    //! Registers a property
    //! \param[in] name      Property name
    //! \param[in] isPrivate True if property is private, defaults to false
    void AddProperty( const std::string& name,
                      bool               isPrivate = false);
    //! Gets property with the given name
    //! \param[in] name Name of the property
    PropertyInfo* GetProperty( const std::string& name) const;
    //! Returns true if a property with the given name exists
    //! \param[in] name Name of the property
    bool HasProperty( const std::string& name) const;
    //! Returns true if property with the given name is private
    //! \param[in] name Name of the property
    bool IsPropertyPrivate( const std::string& name) const;

private:
	std::string                          _class_name;   //! External class name
    std::vector<PropertyInfo*>           _properties;   //! Properties
    std::vector<std::string>             _baseclass;    //! Base classes
    std::map<std::string, FunctionInfo*> _methods;      //! Methods in language
}; 
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//! Property info class
//------------------------------------------------------------------------------
class HML2DLL_DECLS PropertyInfo 
{
public:
    //! Constructor
    //! \param[in] name      Name
    //! \param[in] isPrivate True if property is private, defaults to false
    PropertyInfo( const std::string& name,
                  bool               isPrivate = false) 
        : _name (name), _isPrivate (isPrivate) {}
    //! Destructor
    ~PropertyInfo() {}

    //! Returns property name
    std::string Name() const { return _name; }
    //! True if this is private
    bool IsPrivate() const { return _isPrivate; }

private:
    std::string    _name;       //! Name
    bool           _isPrivate;  //! True if property is private
};

#endif

// End of file: