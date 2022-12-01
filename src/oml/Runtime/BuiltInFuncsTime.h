/**
* @file BuiltInFuncsTime.h
* @date April 2019
* Copyright (C) 2019-2022 Altair Engineering, Inc.
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

#ifndef __BUILTINFUNCSTIME__
#define __BUILTINFUNCSTIME__

#include "OMLDll.h"

#include "EvaluatorInt.h"

//------------------------------------------------------------------------------
//!
//! \class BuiltInFuncsTime
//! \brief Utility class for built-in functions related to time commands
//!
//------------------------------------------------------------------------------
class OMLDLL_DECLS BuiltInFuncsTime
{
public:
    //!
    //! Destructor
    //!
    ~BuiltInFuncsTime() {}

    //!
    //! Returns true and gets the day in the given string [day]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Day(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Returns true and gets the hour in the given string [hour]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Hour(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Returns true and gets the minutes in the given string [minute]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Minute(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Returns true and gets the month from the given string [month]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Month(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Returns true and gets the seconds in the given string [second]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Second(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Returns true and gets the year from the given date string [year]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Year(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Returns true and gets the quarter from the given date string [quarter]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Quarter(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Sets the tic time [tic]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Tic(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Gets the time in seconds from the last tic call [toc]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Toc(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Returns true, gets the days, with Jan 1 0000 considered as day 1 [now]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Now(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
	//!
    //! Gets the day and time input as a number of days since January 1, 0000 [datenum]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Datenum(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
	//!
    //! Gets the date vector from a date number or a date string [datevec]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Datevec(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Gets the day of the week from a date number or a date string [weekday]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Weekday(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
    //!
    //! Gets the date string from a date number or a date string [datestr]
    //! \param Evaluator interface
    //! \param Vector of inputs
    //! \param Vector of outputs
    //! 
    static bool Datestr(EvaluatorInterface, const std::vector<Currency>&, std::vector<Currency>&);
private:
    static std::vector<std::string> _omlDateTimeFmt; //<! Valid date formats
    static std::vector<std::string> _omlMonth;       //<! Month names
    static bool                     g_ticstartSet;   //<! True if g_ticstart has been set


    //!
    //! Constructor
    //!
    BuiltInFuncsTime() {}

    //! 
    //! Gets the format index for the given date string
    //! \param in Given input string
    //!
    size_t GetFormat(const std::string& in);

    //! 
    //! Initializes all the date string formats
    static void InitializeFormats();
};

#endif  // __BUILTINFUNCSTIME__

