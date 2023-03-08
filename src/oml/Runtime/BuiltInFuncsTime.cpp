/**
* @file BuiltInFuncsTime.cpp
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

#include "BuiltInFuncsTime.h"

#include <cassert>
#include <chrono>
#include <iomanip>
#include <map>
#include <array>

#include "BuiltInFuncs.h"
#include "BuiltInFuncsUtils.h"
#include "hwMatrix.h"
#include "OML_Error.h"
#include "OutputFormat.h"

std::vector<std::string> BuiltInFuncsTime::_omlDateTimeFmt;  // Formats
std::vector<std::string> BuiltInFuncsTime::_omlMonth;        // Month list

static std::chrono::time_point<std::chrono::system_clock> g_ticstart;  // tic Start time
static std::map< double, std::chrono::time_point<std::chrono::system_clock> > g_ticmap; // Map of time stamps

bool BuiltInFuncsTime::g_ticstartSet = false; // True once time stamp is set

//------------------------------------------------------------------------------
// Initializes all the date string formats
//------------------------------------------------------------------------------
void BuiltInFuncsTime::InitializeFormats()
{
    if (!_omlDateTimeFmt.empty())
    {
        return;
    }
    _omlDateTimeFmt.reserve(32);  // Max available formats currently

    _omlDateTimeFmt.push_back("dd-mmm-yyyy HH:MM:SS");      // Index 0 
    _omlDateTimeFmt.push_back("dd-mmm-yyyy");               // Index 1
    _omlDateTimeFmt.push_back("mm/dd/yy");                  // Index 2
    _omlDateTimeFmt.push_back("mmm");                       // Index 3
    _omlDateTimeFmt.push_back("m");                         // Index 4
    _omlDateTimeFmt.push_back("mm");                        // Index 5
    _omlDateTimeFmt.push_back("mm/dd");                     // Index 6
    _omlDateTimeFmt.push_back("dd");                        // Index 7
    _omlDateTimeFmt.push_back("ddd");                       // Index 8
    _omlDateTimeFmt.push_back("d");                         // Index 9
    _omlDateTimeFmt.push_back("yyyy");                      // Index 10
    _omlDateTimeFmt.push_back("yy");                        // Index 11
    _omlDateTimeFmt.push_back("mmmyy");                     // Index 12
    _omlDateTimeFmt.push_back("HH:MM:SS");                  // Index 13
    _omlDateTimeFmt.push_back("HH:MM:SS PM");               // Index 14
    _omlDateTimeFmt.push_back("HH:MM");                     // Index 15
    _omlDateTimeFmt.push_back("HH:MM PM");                  // Index 16
    _omlDateTimeFmt.push_back("QQ-YY");                     // Index 17
    _omlDateTimeFmt.push_back("QQ");                        // Index 18
    _omlDateTimeFmt.push_back("dd/mm");                     // Index 19
    _omlDateTimeFmt.push_back("dd/mm/yy");                  // Index 20
    _omlDateTimeFmt.push_back("mmm.dd,yyyy HH:MM:SS");      // Index 21
    _omlDateTimeFmt.push_back("mmm.dd,yyyy");               // Index 22
    _omlDateTimeFmt.push_back("mm/dd/yyyy");                // Index 23
    _omlDateTimeFmt.push_back("dd/mm/yyyy");                // Index 24
    _omlDateTimeFmt.push_back("yy/mm/dd");                  // Index 25
    _omlDateTimeFmt.push_back("yyyy/mm/dd");                // Index 26
    _omlDateTimeFmt.push_back("QQ-YYYY");                   // Index 27
    _omlDateTimeFmt.push_back("mmmyyyy");                   // Index 28
    _omlDateTimeFmt.push_back("yyyy-mm-dd");                // Index 29
    _omlDateTimeFmt.push_back("yyyymmddTHHMMSS");           // Index 30
    _omlDateTimeFmt.push_back("yyyy-mm-dd HH:MM:SS");       // Index 31

    _omlMonth.reserve(12);
    _omlMonth.push_back("jan");
    _omlMonth.push_back("feb");
    _omlMonth.push_back("mar");
    _omlMonth.push_back("apr");
    _omlMonth.push_back("may");
    _omlMonth.push_back("jun");
    _omlMonth.push_back("jul");
    _omlMonth.push_back("aug");
    _omlMonth.push_back("sep");
    _omlMonth.push_back("oct");
    _omlMonth.push_back("nov");
    _omlMonth.push_back("dec");
}
//------------------------------------------------------------------------------
// Returns true and gets the year from the given date string [year]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Year(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    std::string in(inputs[0].StringVal());
    if (in.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }

    BuiltInFuncsTime funcs;
    funcs.InitializeFormats();
    assert(!funcs._omlDateTimeFmt.empty());
    if (funcs._omlDateTimeFmt.empty())
    {
        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
    }
    std::string fmt(funcs._omlDateTimeFmt[1]);  // Default format is dd-mmm-yyyy
    size_t      idx = 1;
    int         fmtidx = 1;
    if (inputs.size() > 1)
    {
        fmtidx = 2;
        if (!inputs[1].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 2);
        }
        idx = funcs.GetFormat(inputs[1].StringVal());
    }

    std::string str;
    switch (idx)
    {
        // dd-mmm-yyyy HH:MM:SS
        case 0: str = "%*d-%*3[^-]-%d %*d:%*d:%*d"; break;

        // dd-mmm-yyyy
        case 1: str = "%*d-%*3[^-]-%d"; break;
        
        // mm/dd/yy dd/mm/yy mm/dd/yyyy dd/mm/yyyy
        case 2:         
        case 20:       
        case 23:        
        case 24: str = "%*d/%*d/%d"; break;

        // yyyy yy
        case 10:        
        case 11: str = "%d"; break;

        // mmmyy mmmyyyy
        case 12:        
        case 28: str = (in.length() >= 5) ? "%*3s%d" : ""; break;

        // QQ-YY QQ-YYYY
        case 17:        
        case 27: str = (in.length() >= 5) ? "%*2[^-]-%d" : ""; break;

        // mmm.dd,yyyy HH:MM:SS
        case 21: str = "%*3[^.].%*d,%d %*d:%*d:%*d"; break;
            
        // mmm.dd,yyyy
        case 22: str = "%*3[^.].%*d,%d"; break;

        // yy/mm/dd yyyy/mm/dd
        case 25:        
        case 26: str = (in.length() >= 5) ? "%d/%*d/%*d" : ""; break;

        // yyyy-mm-dd
        case 29: str = (in.length() >= 5) ? "%d-%*d-%*d" : ""; break;

        // yyyymmddTHHMMSS
        case 30: str = (in.length() >= 15) ? "%4d%*d%*dT%*d%*d%*d" : ""; break;

        // yyyy-mm-dd HH:MM:SS
        case 31: str = "%d-%*d-%*d %*d:%*d:%*d"; break;

        default: break;
    }

    int year   = -1;
    int result = -1;
    if (!str.empty())
    {
        result = std::sscanf(in.c_str(), str.c_str(), &year);
    }
    if (result < 1 || year < 0)
    {
        throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, fmtidx);
    }

    outputs.push_back(year);
    return true;
}
//------------------------------------------------------------------------------
// Gets the format index for the given date string
//------------------------------------------------------------------------------
size_t BuiltInFuncsTime::GetFormat(const std::string& fmt)
{
    if (fmt.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 2);
    }
    std::vector<std::string>::iterator itr = std::find(
        _omlDateTimeFmt.begin(), _omlDateTimeFmt.end(), fmt);
    
    if (itr == _omlDateTimeFmt.end())
    {
        throw OML_Error(OML_ERR_INVALID_DATE_FMT, 2);
    }
    size_t idx = std::distance(_omlDateTimeFmt.begin(), itr);
    return idx;
}
//------------------------------------------------------------------------------
// Returns true and gets the day from the given date string [day]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Day(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    std::string in(inputs[0].StringVal());
    if (in.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }

    BuiltInFuncsTime funcs;
    funcs.InitializeFormats();
    assert(!funcs._omlDateTimeFmt.empty());
    if (funcs._omlDateTimeFmt.empty())
    {
        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
    }
    std::string fmt(funcs._omlDateTimeFmt[1]);  // Default format is dd-mmm-yyyy
    size_t      idx    = 1;
    int         fmtidx = 1;
    if (inputs.size() > 1)
    {
        fmtidx = 2;
        if (!inputs[1].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 2);
        }
        idx = funcs.GetFormat(inputs[1].StringVal());
    }

    std::string str;
    switch (idx)
    {
        // dd-mmm-yyyy HH:MM:SS
        case 0: str = "%d-%*3[^-]-%*d %*d:%*d:%*d"; break;   
            
        // dd-mmm-yyyy
        case 1: str = "%d-%*3[^-]-%*d"; break;

        // mm/dd/yy
        case 2: str = "%*d/%d/%*d"; break;
    
        // dd/mm
        case 19: str = "%d/%*d"; break;

        // dd/mm/yy	dd/mm/yyyy
        case 20:      	
        case 24: str = "%d/%*d/%*d"; break;

        // mmm.dd,yyyy HH:MM:SS
        case 21: str = "%*3[^.].%d,%*d %*d:%*d:%*d"; break;

        // mmm.dd,yyyy
        case 22: str = "%*3[^.].%d,%*d";  break;

        // mm/dd/yyyy
        case 23: str = "%*d/%d/%*d"; break;

        // yy/mm/dd yyyy/mm/dd
        case 25:        
        case 26: str = (in.length() >= 5) ? "%*d/%*d/%d" : ""; break;

        // yyyy-mm-dd
        case 29: str = (in.length() >= 5) ? "%*d-%*d-%d" : ""; break;

        // yyyymmddTHHMMSS
        case 30: str = (in.length() >= 15) ? "%*4d%*2d%dT%*d%*d%*d" : ""; break;

        // yyyy-mm-dd HH:MM:SS
        case 31: str = "%*d-%*d-%d %*d:%*d:%*d"; break;

        default: break;
    }

    int year   = -1;
    int result = -1;
    if (!str.empty())
    {
        result = std::sscanf(in.c_str(), str.c_str(), &year);
    }
    if (result < 1 || year < 0)
    {
        throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, fmtidx);
    }

    outputs.push_back(year);
    return true;
}
//------------------------------------------------------------------------------
// Returns true and gets the month from the given date string [month]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Month(EvaluatorInterface           eval,
                             const std::vector<Currency>& inputs,
                             std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    std::string in(inputs[0].StringVal());
    if (in.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }

    BuiltInFuncsTime funcs;
    funcs.InitializeFormats();
    assert(!funcs._omlDateTimeFmt.empty());
    if (funcs._omlDateTimeFmt.empty())
    {
        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
    }
    std::string fmt(funcs._omlDateTimeFmt[1]);  // Default format is dd-mmm-yyyy
    size_t      idx = 1;
    int         fmtidx = 1;
    if (inputs.size() > 1)
    {
        fmtidx = 2;
        if (!inputs[1].IsString())
        {
            throw OML_Error(OML_ERR_STRING, 2);
        }
        idx = funcs.GetFormat(inputs[1].StringVal());
    }

    std::string str;
    bool isstring = false;
    switch (idx)
    {
        // dd-mmm-yyyy HH:MM:SS
        case 0: str = "%*d-%3[^-]-%*d %*d:%*d:%*d"; isstring = true;  break;

        // dd-mmm-yyyy
        case 1: str = "%*d-%3[^-]-%*d"; isstring = true; break;

        // mm/dd/yy
        case 2: str = "%d/%*d/%*d"; break;

        // mmm
        case 3: str = "%3s"; isstring = true; break;

        // m mm
        case 4:
        case 5: str = "%d"; break;

        // mm/dd
        case 6: str = "%d/%*d"; break;

        // mmmyy mmmyyyy
        case 12:
        case 28: str = "%3s%*d"; isstring = true; break;

        // dd/mm
        case 19: str = "%*d/%d"; break;

        // dd/mm/yy	dd/mm/yyyy
        case 20:
        case 24: str = "%*d/%d/%*d"; break;

        // mmm.dd,yyyy HH:MM:SS
        case 21: str = "%3[^.].%*d,%*d %*d:%*d:%*d"; isstring = true;  break;

        // mmm.dd,yyyy
        case 22: str = "%3[^.].%*d,%*d"; isstring = true; break;

        // mm/dd/yyyy
        case 23: str = "%d/%*d/%*d"; break;

        // yy/mm/dd yyyy/mm/dd
        case 25:
        case 26: str = (in.length() >= 5) ? "%*d/%d/%*d" : ""; break;

        // yyyy-mm-dd
        case 29: str = (in.length() >= 5) ? "%*d-%d-%*d" : ""; break;

        // yyyymmddTHHMMSS
        case 30: str = (in.length() >= 15) ? "%*4d%2d%*2dT%*d%*d%*d" : ""; break;

        // yyyy-mm-dd HH:MM:SS
        case 31: str = "%*d-%d-%*d %*d:%*d:%*d"; break;

        default: break;
    }

    int month  = -1;
    int result = -1;
    std::string desc;
    if (isstring)
    {
        char buff[128];
        result = std::sscanf(in.c_str(), str.c_str(), buff);
        if (result < 1)
        {
            throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, fmtidx); 
        }
        desc = buff;
        std::transform(desc.begin(), desc.end(), desc.begin(), ::tolower);
        std::vector<std::string>::iterator itr  = std::find(_omlMonth.begin(),
            _omlMonth.end(), desc);
        if (itr == _omlMonth.end())
        {
            throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, fmtidx);
        }
        month = static_cast<int>(std::distance(_omlMonth.begin(), itr)) + 1;
    }
    else if (!str.empty())
    {
        result = std::sscanf(in.c_str(), str.c_str(), &month);
    }
    if (result < 1 || month < 1 || month > 12)
    {
        throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, fmtidx);
    }

    outputs.push_back(month);

    if (eval.GetNargoutValue() > 1)
    {
        if (desc.empty())
        {
            desc = _omlMonth[month - 1];
        }

        if (!desc.empty())
        {
            std::transform(desc.begin(), desc.begin() + 1, desc.begin(), ::toupper);
        }
        outputs.push_back(desc);
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns true and gets the minutes in the given date string [minute]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Minute(EvaluatorInterface           eval,
                              const std::vector<Currency>& inputs,
                              std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    std::string in(inputs[0].StringVal());
    if (in.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }

    if (inputs.size() == 1)
    {
        outputs.push_back(0);  // Minutes are 0 in default format
        return true;
    }

    BuiltInFuncsTime funcs;
    funcs.InitializeFormats();
    assert(!funcs._omlDateTimeFmt.empty());
    if (funcs._omlDateTimeFmt.empty())
    {
        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
    }

    std::string fmt(funcs._omlDateTimeFmt[1]);  // Default format is dd-mmm-yyyy
    if (!inputs[1].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 2);
    }
    size_t idx = funcs.GetFormat(inputs[1].StringVal());

    std::string str;
    switch (idx)
    {
        case 0:   // dd-mmm-yyyy HH:MM:SS
            str = "%*d-%*3[^-]-%*d %*d:%d:%*d"; 
            break;

        case 1:   // dd-mmm-yyyy
        case 2:   // mm/dd/yy
        case 3:   // mmm
        case 6:   // mm
        case 7:   // mm/dd
        case 8:   // dd
        case 10:  // ddd
        case 11:  // yyyy
        case 12:  // mmmyy
        case 19:  // dd/mm
        case 20:  // dd/mm/yy
        case 22:  // mmm.dd, yyyy
        case 23:  // mm/dd/yyyy
        case 24:  // dd/mm/yyyy
        case 28:  // mmmyyyy
            outputs.push_back(0);
            return true;

        case 13:  // HH:MM:SS
            str = "%*d:%d:%*d"; 
            break;

        case 14:  // HH:MM:SS PM
            str = "%*d:%d:%*d %*s";
            break;

        case 15:  // HH:MM
            str = "%*d:%d";
            break;

        case 16:  // HH:MM PM
            str = "%*d:%d %*s";
            break;

        case 21:  // mmm.dd,yyyy HH:MM:SS
            str = "%*3[^.].%*d,%*d %*d:%d:%*d";
            break;
        
        case 30:  // yyyymmddTHHMMSS
            str = (in.length() >= 15) ? "%*4d%*2d%*2dT%*2d%2d%*2d" : ""; 
            break;

        case 31:  // yyyy-mm-dd HH:MM:SS
            str = "%*d-%*d-%*d %*d:%d:%*d"; break;

        default: break;
    }

    int value = -1;
    int result = -1;
    if (!str.empty())
    {
        result = std::sscanf(in.c_str(), str.c_str(), &value);
    }
    if (result < 1 || value < 0 || value > 60)
    {
        throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, 2);
    }

    outputs.push_back(value);
    return true;
}
//------------------------------------------------------------------------------
// Returns true and gets the hour in the given date string [hour]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Hour(EvaluatorInterface           eval,
                            const std::vector<Currency>& inputs,
                            std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    std::string in(inputs[0].StringVal());
    if (in.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }

    if (inputs.size() == 1)
    {
        outputs.push_back(0);  // Hours are 0 in default format
        return true;
    }

    BuiltInFuncsTime funcs;
    funcs.InitializeFormats();
    assert(!funcs._omlDateTimeFmt.empty());
    if (funcs._omlDateTimeFmt.empty())
    {
        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
    }

    std::string fmt(funcs._omlDateTimeFmt[1]);  // Default format is dd-mmm-yyyy
    if (!inputs[1].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 2);
    }

    size_t idx = funcs.GetFormat(inputs[1].StringVal());

    std::string str;
    bool hasTimePeriod = false;
    switch (idx)
    {
        case 0:   // dd-mmm-yyyy HH:MM:SS
            str = "%*d-%*3[^-]-%*d %d:%*d:%*d";
            break;

        case 1:   // dd-mmm-yyyy
        case 2:   // mm/dd/yy
        case 3:   // mmm
        case 6:   // mm
        case 7:   // mm/dd
        case 8:   // dd
        case 10:  // ddd
        case 11:  // yyyy
        case 12:  // mmmyy
        case 19:  // dd/mm
        case 20:  // dd/mm/yy
        case 22:  // mmm.dd, yyyy
        case 23:  // mm/dd/yyyy
        case 24:  // dd/mm/yyyy
        case 28:  // mmmyyyy
            outputs.push_back(0);
            return true;

        case 13:  // HH:MM:SS
            str = "%d:%*d:%*d";
            break;

        case 14:  // HH:MM:SS PM
            str = "%d:%*d:%*d %s";
            hasTimePeriod = true;
            break;

        case 15:  // HH:MM
            str = "%d:%*d";
            break;

        case 16:  // HH:MM PM
            str = "%d:%*d %s";
            hasTimePeriod = true;
            break;

        case 21:  // mmm.dd,yyyy HH:MM:SS
            str = "%*3[^.].%*d,%*d %d:%*d:%*d";
            break;

        case 30:  // yyyymmddTHHMMSS
            str = (in.length() >= 15) ? "%*4d%*2d%*2dT%2d%*2d%*2d" : "";
            break;

        case 31:  // yyyy-mm-dd HH:MM:SS
            str = "%*d-%*d-%*d %d:%*d:%*d"; break;

        default: break;
    }

    int value = -1;
    int result = -1;
    if (!str.empty())
    {
        if (!hasTimePeriod)
        {
            result = std::sscanf(in.c_str(), str.c_str(), &value);
        }
        else
        {
            char tmp[128];
            result = std::sscanf(in.c_str(), str.c_str(), &value, tmp);
            std::string timePeriod(tmp);
            if (!timePeriod.empty())
            {
                std::transform(timePeriod.begin(), timePeriod.end(),
                    timePeriod.begin(), ::tolower);
                if (value < 12 && timePeriod == "pm")
                {
                    value += 12;
                }
            }
        }
    }
    if (result < 1 || value < 0 || value > 24)
    {
        value = 0; // Just return 0 if format is not applicable
    }

    outputs.push_back(value);
    return true;
}
//------------------------------------------------------------------------------
// Returns true and gets the seconds in the given string [second]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Second(EvaluatorInterface           eval,
                              const std::vector<Currency>& inputs,
                              std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    std::string in(inputs[0].StringVal());
    if (in.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }

    if (inputs.size() == 1)
    {
        outputs.push_back(0);  // Seconds are 0 in default format
        return true;
    }

    BuiltInFuncsTime funcs;
    funcs.InitializeFormats();
    assert(!funcs._omlDateTimeFmt.empty());
    if (funcs._omlDateTimeFmt.empty())
    {
        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
    }

    std::string fmt(funcs._omlDateTimeFmt[1]);  // Default format is dd-mmm-yyyy
    if (!inputs[1].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 2);
    }
    size_t idx = funcs.GetFormat(inputs[1].StringVal());

    std::string str;
    switch (idx)
    {
        case 0:   // dd-mmm-yyyy HH:MM:SS
            str = "%*d-%*3[^-]-%*d %*d:%*d:%d";
            break;

        case 1:   // dd-mmm-yyyy
        case 2:   // mm/dd/yy
        case 3:   // mmm
        case 6:   // mm
        case 7:   // mm/dd
        case 8:   // dd
        case 10:  // ddd
        case 11:  // yyyy
        case 12:  // mmmyy
        case 15:  // HH:MM
        case 16:  // HH:MM PM
        case 19:  // dd/mm
        case 20:  // dd/mm/yy
        case 22:  // mmm.dd, yyyy
        case 23:  // mm/dd/yyyy
        case 24:  // dd/mm/yyyy
        case 28:  // mmmyyyy
            outputs.push_back(0);
            return true;

        case 13:  // HH:MM:SS
            str = "%*d:%*d:%d";
            break;

        case 14:  // HH:MM:SS PM
            str = "%*d:%*d:%d %*s";
            break;

        case 21:  // mmm.dd,yyyy HH:MM:SS
            str = "%*3[^.].%*d,%*d %*d:%*d:%d";
            break;

        case 30:  // yyyymmddTHHMMSS
            str = (in.length() >= 15) ? "%*4d%*2d%*2dT%*2d%*2d%2d" : "";
            break;

        case 31:  // yyyy-mm-dd HH:MM:SS
            str = "%*d-%*d-%*d %*d:%*d:%d"; break;

        default: break;
    }

    int value = -1;
    int result = -1;
    if (!str.empty())
    {
        result = std::sscanf(in.c_str(), str.c_str(), &value);
    }
    if (result < 1 || value < 0 || value > 60)
    {
        throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, 2);
    }

    outputs.push_back(value);
    return true;
}
//------------------------------------------------------------------------------
// Returns true and gets the quarter from the given date string [quarter]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Quarter(EvaluatorInterface           eval,
                               const std::vector<Currency>& inputs,
                               std::vector<Currency>&       outputs)
{
    if (inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    if (!inputs[0].IsString())
    {
        throw OML_Error(OML_ERR_STRING, 1);
    }

    std::string in(inputs[0].StringVal());
    if (in.empty())
    {
        throw OML_Error(OML_ERR_NONEMPTY_STR, 1);
    }

    int quarterstart = 1;  // Month number when quarter starts
    if (inputs.size() > 1)
    {
        if (!inputs[1].IsPositiveInteger())
        {
            throw OML_Error(OML_ERR_POSINTEGER, 2);
        }
        int val = static_cast<int>(inputs[1].Scalar());
        if (val < 1 || val > 12)
        {
            throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, 2);
        }
        quarterstart = val;
    }

    BuiltInFuncsTime funcs;
    funcs.InitializeFormats();
    assert(!funcs._omlDateTimeFmt.empty());
    if (funcs._omlDateTimeFmt.empty())
    {
        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
    }

    int month = -1;
    std::string desc;

    char buff[128];
    int  result = std::sscanf(in.c_str(), "%*d-%3[^-]-%*d", buff);
    if (result < 1)
    {
        throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, 1);
    }
    desc = buff;
    std::transform(desc.begin(), desc.end(), desc.begin(), ::tolower);
    std::vector<std::string>::iterator itr = std::find(_omlMonth.begin(),
        _omlMonth.end(), desc);
    if (itr == _omlMonth.end())
    {
        throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, 1);
    }
    month = static_cast<int>(std::distance(_omlMonth.begin(), itr)) + 1;

    int monthDiff = month - quarterstart;
    if (monthDiff < 0)
    {
        monthDiff += 12;
    }
    int quarter = (monthDiff + 3) / 3;

    outputs.push_back(quarter);
    return true;
}
//------------------------------------------------------------------------------
// Sets the tic time [tic]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Tic(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    if (!inputs.empty())
    {
        BuiltInFuncsUtils::SetWarning(eval, "Warning: ignoring arguments in [tic]");
    }

    g_ticstartSet = true;

    if (eval.GetNargoutValue() <= 0)
    {
        g_ticstart = std::chrono::system_clock::now();
    }
    else
    {
        clock_t ct;

#       ifdef _SC_CLK_TCK
            ct = times(nullptr);
#       else
            ct = clock();
#       endif

        double key = static_cast<double>(ct);
        g_ticmap.insert(std::make_pair(key, std::chrono::system_clock::now()));
        outputs.emplace_back(key);
    }

    return true;
}
//------------------------------------------------------------------------------
// Gets the time in seconds from the last tic call [toc]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Toc(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    auto end = std::chrono::system_clock::now();

    if (inputs.size() > 1)
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }
    else if (!g_ticstartSet)
    {
        BuiltInFuncsUtils::SetWarning(eval, "Warning: tic must be called first");
        return true;
    }

    std::chrono::duration<double> diff = end - g_ticstart;

    if (inputs.size() == 1)
    {
        const Currency& in = inputs[0];
        if (!in.IsScalar() && !in.IsComplex())
        {
            throw OML_Error(OML_ERR_SCALAR_COMPLEX, 1, OML_VAR_TYPE);
        }

        double key = (in.IsScalar()) ? in.Scalar() : (in.Complex().Real());
        if (!g_ticmap.empty() && g_ticmap.find(key) != g_ticmap.end())
        {
            diff = end - g_ticmap[key];
        }
    }

    double timediff = diff.count();
    if (eval.GetNargoutValue() == 0)
    {
        std::ostringstream os;
        os << "Elapsed time is ";
        os << std::fixed << std::setprecision(OutputFormat::PRECISION_LONG);
        os << timediff << " seconds.";
        eval.PrintResult(os.str());
    }
    else
    {
        outputs.emplace_back(timediff);
    }

    return true;
}
//------------------------------------------------------------------------------
// assumes no leap year
//------------------------------------------------------------------------------
int getDaysInMonth(int month)
{
    switch (month)
    {
    case 1:
    case 3:
    case 5:
    case 7:
    case 8:
    case 10:
    case 12:
        return 31;
    case 2:
        return 28;
    case 4:
    case 6:
    case 9:
    case 11:
        return 30;
    default:
        throw OML_Error(GetHMathErrMsg(HW_MATH_ERR_INTERNALERROR));
    }
}
//------------------------------------------------------------------------------
// Returns true if leap year
//------------------------------------------------------------------------------
bool isLeapYear(int year)
{
    return !year || (isint(year / 4.0) && (isint(year / 400.0) || !isint(year / 100.0)));
}
//------------------------------------------------------------------------------
// assumes months is between 1 and 12 inclusively
//------------------------------------------------------------------------------
double daysInMonths(int monthsPassed, int year)
{
    double rv = 0.0;
    //int monthsPassed = (int) months;
    for (int i = 1; i <= monthsPassed; i++)
    {
        rv += getDaysInMonth(i);
    }
    // rv += getDaysInMonth(monthsPassed + 1) * (months - monthsPassed);
    if (monthsPassed >= 2 && isLeapYear(year))
        rv++;
    return rv;
}
//------------------------------------------------------------------------------
// if matrices are used as inputs, an array is returned instead of just one
// assumes order of year, month, day, hour, minute, second
// the latter three are optional
// no support yet for date vector
//------------------------------------------------------------------------------
std::vector<std::array<double, 6> > datesFromInput(const std::vector<Currency>& inputs, int* m, int* n)
{
    size_t nargin = inputs.size();

    if (nargin < 1 || nargin > 6 || nargin == 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    // get dimensions of date
    *m = *n = 1;

    for (size_t i = 0; i < nargin; i++)
    {
        const Currency& cur = inputs[i];
        if (cur.IsMatrix() || cur.IsString())
        {
            const hwMatrix* mtx = cur.Matrix();

            if (!mtx->IsRealData())
                throw OML_Error(HW_ERROR_DATENOTCOMP);
            if (mtx->Size() != 1)
            {
                if (*m == 1 && *n == 1)
                {
                    *m = mtx->M();
                    *n = mtx->N();
                }
                else if (!(*m == mtx->M() && *n == mtx->N()))
                    throw OML_Error(HW_ERROR_MATRIXDIM);
            }
        }
    }

    size_t size = (size_t)(*m * *n);
    std::vector<std::array<double, 6> > dates;

    for (size_t i = 0; i < size; i++)
    {
        std::array<double, 6> d = { 0, 0, 0, 0, 0, 0 };
        dates.push_back(d);
    }

    for (size_t i = 0; i < nargin; i++)
    {
        const Currency& cur = inputs[i];
        if (cur.IsScalar() || (cur.IsComplex() && iszero(cur.Complex().Imag())))
        {
            for (size_t j = 0; j < size; j++)
            {
                dates[j][i] = cur.IsScalar() ? cur.Scalar() : cur.Complex().Real();
            }
        }
        else if (cur.IsComplex())
        {
            throw OML_Error(HW_ERROR_DATENOTCOMP);
        }
        else if (cur.IsMatrix() || cur.IsString())
        {
            const hwMatrix* mtx = cur.Matrix();
            for (size_t j = 0; j < size; j++)
                dates[j][i] = (*mtx)((int)j);
        }
        else
            throw OML_Error(OML_ERR_REAL, 1, OML_VAR_VALUE);
    }

    // month-related corrections
    for (size_t i = 0; i < size; i++)
    {
        // don't allow months less than 1
        if (dates[i][1] < 1.0)
            dates[i][1] = 1.0;

        if (!isint(dates[i][1]))
            throw OML_Error(OML_ERR_POSINTEGER);
    }

    return dates;
}
//------------------------------------------------------------------------------
// Returns true, gets the days, with Jan 1 0000 considered as day 1 [now]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Now(EvaluatorInterface           eval,
                           const std::vector<Currency>& inputs,
                           std::vector<Currency>&       outputs)
{
    if (!inputs.empty())
    {
        throw OML_Error(OML_ERR_NUMARGIN);
    }

    time_t rawtime = time(nullptr);
    struct tm* curtime = localtime(&rawtime);

    double year = curtime->tm_year + 1900;
    double mon = curtime->tm_mon;
    double day = curtime->tm_mday;

    year += static_cast<int>(mon / 12.0);
    mon = rem(mon, 12.0);

    // deal with leap years
    double result = year-- * 365;
    result += floor(year / 4.0) - floor(year / 100.0) + floor(year / 400.0) + 1;

    //months
    result += daysInMonths((int)mon, (int)(year + 1.0));

    // days
    result += day;

    // Integer part of result will give the number of days and the fractional
    // part of results should give the number of seconds
    struct tm start;
    start.tm_year = curtime->tm_year;
    start.tm_mon = curtime->tm_mon;
    start.tm_mday = curtime->tm_mday;
    start.tm_hour = 0;
    start.tm_min = 0;
    start.tm_sec = 0;
    start.tm_isdst = -1;

    std::chrono::system_clock::time_point today =
        std::chrono::system_clock::from_time_t(std::mktime(&start));
    std::chrono::system_clock::duration d = std::chrono::system_clock::now() - today;

    std::chrono::seconds s = std::chrono::duration_cast<std::chrono::seconds>(d);
    long long secs = s.count();
    int numdigits = (secs < 1) ? 1 : static_cast<int>(log10(secs) + 1);

    result += secs / pow(10, numdigits);

    outputs.push_back(result);

    return true;
}
//------------------------------------------------------------------------------
// Returns the day/time input as a number of days since January 1, 0000 [datenum]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Datenum(EvaluatorInterface           eval, 
                               const std::vector<Currency>& inputs, 
                               std::vector<Currency>&       outputs)
{
    int nargout = eval.GetNargoutValue();
    int nargin = (int)inputs.size();
    int m, n;

    std::vector<std::array<double, 6> > dates = datesFromInput(inputs, &m, &n);

    hwMatrix* result = EvaluatorInterface::allocateMatrix(m, n, true);
    hwMatrix* resultSecs = EvaluatorInterface::allocateMatrix(m, n, true);

    double daynum;

    for (int i = 0; i < dates.size(); i++)
    {
        std::array<double, 6> date = dates[i];
        double years = date[0];
        double months = date[1] - 1;

        years += (int)(months / 12.0);
        months = rem(months, 12.0);

        // deal with leap years
        daynum = years-- * 365;
        daynum += floor(years / 4.0) - floor(years / 100.0) + floor(years / 400.0) + 1;

        //months
        daynum += daysInMonths((int)months, (int)(years + 1.0));

        // days
        daynum += date[2];

        // hours
        daynum += date[3] / 24.0;

        // minutes
        daynum += date[4] / 1440.0;

        // seconds
        daynum += date[5] / 86400.0;

        (*result)(i) = daynum;
        (*resultSecs)(i) = daynum * 60 * 60 * 24;
    }
    outputs.push_back(result);

    if (nargout > 1)
    {
        outputs.push_back(resultSecs);
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns the name of the month, month must be an integer between 1 and 12.
// If abbr is True then the abbreviated name (first 3 letters) is returned.
//------------------------------------------------------------------------------
std::string getMonthName(int month, bool abbr)
{
    std::string name;
    switch (month)
    {
    case 1:     name = std::string("January");      break;
    case 2:     name = std::string("February");     break;
    case 3:     name = std::string("March");        break;
    case 4:     name = std::string("April");        break;
    case 5:     name = std::string("May");          break;
    case 6:     name = std::string("June");         break;
    case 7:     name = std::string("July");         break;
    case 8:     name = std::string("August");       break;
    case 9:     name = std::string("September");    break;
    case 10:    name = std::string("October");      break;
    case 11:    name = std::string("November");     break;
    case 12:    name = std::string("December");     break;
    default:
        return std::string();
    }
    if (abbr)
        name.resize(3);
    return name;
}
//------------------------------------------------------------------------------
// Returns the name of day, day must be an integer between 1 and 7 (1 is Sunday)
// If abbr is True then the abbreviated name (first 3 letters) is returned.
//------------------------------------------------------------------------------
std::string getDayName(int day, bool abbr)
{
    std::string name;
    switch (day)
    {
    case 1: name = std::string("Sunday");       break;
    case 2: name = std::string("Monday");       break;
    case 3: name = std::string("Tuesday");      break;
    case 4: name = std::string("Wednesday");    break;
    case 5: name = std::string("Thursday");     break;
    case 6: name = std::string("Friday");       break;
    case 7: name = std::string("Saturday");     break;
    default:
        return std::string();
    }
    if (abbr)
        name.resize(3);
    return name;
}
//------------------------------------------------------------------------------
// Replaces the given pattern in the input date string.
//------------------------------------------------------------------------------
void findDatePatternAndReplace(std::string& date, const std::string& pattern, const std::string& replace)
{
    size_t s = 0;
    while ((s = date.find(pattern, s)) != std::string::npos)
    {
        date.replace(s, pattern.length(), replace);
        s += 1 + replace.length();
    }
}
//------------------------------------------------------------------------------
// Replaces the user date format pattern in userFMT with userReplace and the
// pattern in replaceMillisecFMT with millisecReplace.
// 
// Returns the number of times the pattern was found.
//------------------------------------------------------------------------------
int replacePatternInUserFMT(std::string& userFMT, std::string& millisecFMT,
    bool replaceMillisecFMT, const std::string& pattern,
    const std::string& userReplace,
    const std::string& millisecReplace)
{
    int N = 0;
    size_t s = 0;
    while ((s = userFMT.find(pattern, s)) != std::string::npos)
    {
        userFMT.replace(s, pattern.length(), userReplace);
        s += 1 + userReplace.length();
        ++N;
    }

    if (replaceMillisecFMT)
    {
        s = 0;
        while ((s = millisecFMT.find(pattern, s)) != std::string::npos)
        {
            millisecFMT.replace(s, pattern.length(), millisecReplace);
            s += 1 + millisecReplace.length();
        }
    }
    return N;
}
//------------------------------------------------------------------------------
// Converts the input date strings to date vectors. If userFmt is empty then the 
// predefined date formats (see help) will be used to determine the format.
// \param dates   Date strings
// \param userFmt Format of the date string
// \param refYear Reference year, for 2-digit year formats
// \param curYear Current year
// \param out     Nx6 matrix containing the dates as vector (output)
// 
// Returns 0 if successfull, 1 if the date could not be read, 
// 2 if the userFmt is wrong.
//------------------------------------------------------------------------------
int getDateVectorFromString(const std::vector<std::string>& dates, 
                            const std::string& userFmt, int refYear, 
                            int curYear, hwMatrix* out)
{
    assert(out && out->M() == dates.size() && out->N() == 6);
    if (!(out && out->M() == dates.size() && out->N() == 6))
        return 1;

    if (userFmt.empty())
    {
        // format id <-> fmt string mapping
        std::vector<std::string> dateFMTVec;
        dateFMTVec.push_back("%d-%b-%Y %H:%M:%S");
        dateFMTVec.push_back("%b.%d,%Y %H:%M:%S");
        dateFMTVec.push_back("%Y-%m-%d %H:%M:%S");
        dateFMTVec.push_back("%H:%M:%S %p");
        dateFMTVec.push_back("%H:%M:%S");
        dateFMTVec.push_back("%H:%M %p");
        dateFMTVec.push_back("%H:%M");
        dateFMTVec.push_back("%Y-%m-%d");
        dateFMTVec.push_back("%d-%b-%Y");
        dateFMTVec.push_back("%b%y");
        dateFMTVec.push_back("%b.%d,%Y");
        dateFMTVec.push_back("%b%Y");

        // try to find out the format of each date
        for (int i = 0; i < static_cast<int>(dates.size()); ++i)
        {
            std::string date(dates[i]);
            int year = 0, month = 1, day = 1, hour = 0, min = 0, sec = 0;

            // check format 'dd mmm yyyy HH:MM:SS' or 'dd mmmm yyyy HH:MM:SS' - these are not in the
            // predefined formats but others seem to support them
            bool containsDash = date.find("-") != std::string::npos;
            std::string tmpFmt = "%d %s %d %d:%d:%d";
            char monthBuf[32];
            int ret = std::sscanf(date.c_str(), tmpFmt.c_str(), &day, &monthBuf, &year, &hour, &min, &sec);
            if (!containsDash && (ret == 3 || ret == 6))
            {
                // find the month
                std::string mStr(monthBuf);
                bool abbr = mStr.length() == 3;
                month = -1;
                for (int i = 1; i <= 12; ++i)
                {
                    std::string mN = getMonthName(i, abbr);
                    if (mN == mStr)
                    {
                        month = i;
                        break;
                    }
                }
                if (month == -1 || day <= 0 || day > 31)
                    return 1;
            }
            // check formats containing '/'
            else if (date.find("/") != std::string::npos)
            {
                char buf1[12], buf2[12], buf3[12];
                std::string tmp = "%[^/]/%[^/]/%[^/]";
                int result = std::sscanf(date.c_str(), tmp.c_str(), &buf1, &buf2, &buf3);
                if (result != 3 && result != 2)
                    return 1;
                
                if (result == 3)
                {
                    std::string a(buf1), b(buf2), c(buf3);
                    if (a.length() == 4)
                    {
                        // yyyy/mm/dd
                        std::string tmp1 = "%d/%d/%d";
                        std::sscanf(date.c_str(), tmp1.c_str(), &year, &month, &day);
                    }
                    else if (c.length() == 4)
                    {
                        // mm/dd/yyyy
                        std::string tmp1 = "%d/%d/%d";
                        std::sscanf(date.c_str(), tmp1.c_str(), &month, &day, &year);
                    }
                    else
                    {
                        // assume it is mm/dd/yy
                        std::string tmp1 = "%d/%d/%d";
                        std::sscanf(date.c_str(), tmp1.c_str(), &month, &day, &year);
                        if (month > 12)
                        {
                            // it is yy/mm/dd
                            std::sscanf(date.c_str(), tmp1.c_str(), &year, &month, &day);
                        }
                        // adjust year based on the reference year
                        int d = (int)mod(refYear, 100);
                        if (year >= d)
                            year += (refYear / 100) * 100;
                        else
                            year += ((refYear + 100) / 100) * 100;
                    }
                }
                else if (result == 2)
                {
                    year = curYear;
                    // assume it is mm/dd
                    std::string tmp1 = "%d/%d";
                    std::sscanf(date.c_str(), tmp1.c_str(), &month, &day);
                    if (month > 12)
                    {
                        // it is dd/mm
                        int tmpd = day;
                        day = month;
                        month = tmpd;
                    }
                }
                // validate 1<= month <=12 and 1 <= day <= 31 
                if (month < 1 || month>12 || day < 1 || day > 31)
                    return 1;
            }
            else
            {
                // look into the other predefined formats               
                bool fmtFound = false;
                for (int i = 0; i < (int)dateFMTVec.size(); ++i)
                {
                    std::string tmpD(date);
                    std::string tmpFMT(dateFMTVec[i]);
                    bool hasPM = false;
                    // Linux is not reading the %p identifier - remove it
#ifndef OS_WIN
                    if (i == 3 || i == 5)
                    {
                        std::transform(tmpD.begin(), tmpD.end(), tmpD.begin(), ::tolower);
                        bool hasAM = tmpD.find("am") != std::string::npos;
                        hasPM = tmpD.find("pm") != std::string::npos;
                        if (!hasAM && !hasPM)
                            continue;
                        findDatePatternAndReplace(tmpD, "am", "");
                        findDatePatternAndReplace(tmpD, "pm", "");
                        findDatePatternAndReplace(tmpFMT, "%p", "");
                    }
#endif // !OS_WIN

                    tm t = {};
                    std::istringstream iss{ tmpD.c_str() };
                    iss >> std::get_time(&t, tmpFMT.c_str());
                    if (!iss.fail())
                    {
                        year = t.tm_year + 1900;
                        month = t.tm_mon + 1;
                        day = t.tm_mday;
                        hour = t.tm_hour;
                        min = t.tm_min;
                        sec = t.tm_sec;
                        if (hasPM)
                            hour += 12;
                        
                        // format has only the time
                        if (i >= 3 && i <= 6)
                        {
                            year = curYear;
                            month = 1;
                            day = 1;
                            if (i == 3 || i == 5)
                            {
                                if (hour == 24)
                                    hour = 12;
                                else if (hour == 12)
                                    hour = 0;
                            }
                        }
                        else if (i == 9 || i == 11)
                        {
                            day = 1;
                            if (i == 9)
                            {
                                // adjust year based on the reference year
                                int rd = (int)mod(refYear, 100);
                                int yd = (int)mod(year, 100);
                                if (yd >= rd)
                                    year = yd + (refYear / 100) * 100;
                                else
                                    year = yd + ((refYear + 100) / 100) * 100;
                            }
                        }

                        // date found, break;
                        fmtFound = true;
                        break;
                    }
                }
                if (!fmtFound)
                    return 1;
            }

            (*out)(i, 0) = year;
            (*out)(i, 1) = month;
            (*out)(i, 2) = day;
            (*out)(i, 3) = hour;
            (*out)(i, 4) = min;
            (*out)(i, 5) = sec;
        }
    }
    else
    {
        // validate the pattern first
        std::string fmtStr(userFmt);
        std::string millisecFmt(userFmt);

        // millisec
        size_t millisecIdx = userFmt.find("FFF");
        bool hasMilliseconds = millisecIdx != std::string::npos;
        int ms = replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "FFF", "", "%d");
        if (ms > 1) return 2;

        // year
        int yearOcc = replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "yyyy", "%Y", "%*d");
        yearOcc += replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "YYYY", "%Y", "%*d");
        if (yearOcc > 1)
            return 2;

        int yearShortOcc = replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "yy", "%y", "%*d");
        yearShortOcc += replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "YY", "%y", "%*d");
        if ((yearOcc + yearShortOcc) > 1)
            return 2;

        // month
        int monthOcc = replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "mmmm", "%B", "%*s");
        monthOcc += replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "mmm", "%h", "%*s");
        monthOcc += replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "mm", "%m", "%*d");
        if (monthOcc > 1)
            return 2;

        // day
        int dayOWOcc = replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "dddd", "%A", "%*s");
        dayOWOcc += replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "ddd", "%a", "%*s");
        if (dayOWOcc > 1)
            return 2;
        int dayNOWOcc = replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "dd", "%d", "%*d");
        if (dayNOWOcc > 1) return 2;

        // hour
        // am/pm
        int ampmOcc = replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "AM", "%p", "%*s");
        ampmOcc += replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "PM", "%p", "%*s");
        if (ampmOcc > 1)
            return 2;

        std::string tmp = ampmOcc == 1 ? "%I" : "%H";
        if (replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "HH", tmp, "%*d") > 1) return 2;
        // min
        if (replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "MM", "%M", "%*d") > 1) return 2;
        // sec
        if (replacePatternInUserFMT(fmtStr, millisecFmt, hasMilliseconds, "SS", "%S", "%*d") > 1) return 2;

        // check if the fmt contains '%p' - will need special handling on Linux
#ifndef OS_WIN
        if (ampmOcc == 1)
        {
            findDatePatternAndReplace(fmtStr, "%p", "");
        }
#endif

        // read each date
        for (int i = 0; i < static_cast<int>(dates.size()); ++i)
        {
            std::string date(dates[i]);

            // extract the milliseconds and replace the millisec number in the date
            int msec = 0;
            if (hasMilliseconds)
            {
                if (msec < 0 || msec>999)
                    return 1;
                int res = sscanf(date.c_str(), millisecFmt.c_str(), &msec);
                if (res == 1)
                {
                    std::string tmp = std::to_string(msec);
                    size_t s = 0;
                    while ((s = date.find(tmp, s)) != std::string::npos)
                    {
                        // replace only if it is in a position >= than the position
                        // in the original fmt
                        if (s >= millisecIdx)
                        {
                            date.replace(s, tmp.length(), " ");
                            break;
                        }
                        ++s;
                    }
                }
            }

            bool hasPM = false;
#ifndef OS_WIN
            if (ampmOcc == 1)
            {
                std::transform(date.begin(), date.end(), date.begin(), ::tolower);
                bool hasAM = date.find("am") != std::string::npos;
                hasPM = date.find("pm") != std::string::npos;
                findDatePatternAndReplace(date, "am", "");
                findDatePatternAndReplace(date, "pm", "");
            }
#endif // !OS_WIN

            tm t = {};
            std::istringstream iss{ date.c_str() };
            iss >> std::get_time(&t, fmtStr.c_str());
            if (iss.fail())
                return 1;

            (*out)(i, 0) = t.tm_year + 1900;
            (*out)(i, 1) = t.tm_mon + 1;
            (*out)(i, 2) = t.tm_mday;
            (*out)(i, 3) = t.tm_hour;
            (*out)(i, 4) = t.tm_min;
            (*out)(i, 5) = t.tm_sec + (double)msec / 1000;

            // might not have a date - use 1/1/currentYear
            if (yearOcc == 0 && yearShortOcc == 0)
                (*out)(i, 0) = curYear;
            if (monthOcc == 0)
                (*out)(i, 1) = 1;
            if (dayOWOcc == 0 && dayNOWOcc == 0)
                (*out)(i, 2) = 1;

            if (hasPM)
                (*out)(i, 3) += 12;

            if (yearShortOcc)
            {
                // adjust year based on the reference year
                int rd = (int)mod(refYear, 100);
                int yd = (int)mod((*out)(i, 0), 100);
                if (yd >= rd)
                    (*out)(i, 0) = yd + (refYear / 100) * 100;
                else
                    (*out)(i, 0) = yd + ((refYear + 100) / 100) * 100;
            }
        }
    }
    return 0;
}
//------------------------------------------------------------------------------
// Gets the date vector from a date number or a date string [datevec]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Datevec(EvaluatorInterface           eval,
                               const std::vector<Currency>& inputs,
                               std::vector<Currency>&       outputs)
{
    int nargin = (int)inputs.size();
    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency vecMat;
    if (inputs[0].IsScalar() || inputs[0].IsMatrix())
    {
        const hwMatrix* in = inputs[0].IsScalar() ? inputs[0].ConvertToMatrix() : inputs[0].Matrix();
        hwMatrix* out = EvaluatorInterface::allocateMatrix(in->Size(), 6, true);
        vecMat = Currency(out);
        for (int i = 0; i < in->Size(); ++i)
        {
            double year = 0, month = 1, day = 1;
            double days = in->IsReal() ? (*in)(i) : in->z(i).Real();
            double daysPassed = days < 0 ? ceil(days) : floor(days);
            double remTime = days - daysPassed;
            if (days < 0 && abs(remTime) > 0)
            {
                daysPassed -= 1;
                remTime = 1 - abs(remTime);
            }

            year = floor(daysPassed / 365);
            double leapYadj = floor((year - 1) / 4.0) - floor((year - 1) / 100.0) + floor((year - 1) / 400.0) + 1;
            double calcDays = year * 365 + leapYadj;
            double remainingDays = daysPassed - calcDays;
            bool ly = isLeapYear((int)year);
            while (remainingDays <= 0 || (remainingDays > 365 && !ly) || (remainingDays > 366 && ly))
            {
                if (remainingDays <= 0)
                    --year;
                else
                    ++year;
                leapYadj = floor((year - 1) / 4.0) - floor((year - 1) / 100.0) + floor((year - 1) / 400.0) + 1;
                calcDays = year * 365 + leapYadj;
                remainingDays = daysPassed - calcDays;
                ly = isLeapYear((int)year);

            }

            // calculate exact month
            if (remainingDays > 0)
            {
                double tmp = remainingDays;
                bool isLY = isLeapYear((int)year);
                for (int i = 1; i < 13; ++i)
                {
                    if (i == 12)
                    {
                        month = i;
                        day = tmp;
                        break;
                    }
                    int dm = getDaysInMonth(i);
                    if (i == 2 && isLY)
                        dm++;
                    if (tmp - dm > 0)
                    {
                        tmp -= dm;
                    }
                    else
                    {
                        month = i;
                        day = tmp;
                        break;
                    }
                }
            }

            // milliseconds
            double remTimeInMilliSecs = round(remTime * 24 * 360000);
            double hour = floor(remTimeInMilliSecs / 360000);
            double min = floor((remTimeInMilliSecs - hour * 360000) / 6000);
            double sec = (remTimeInMilliSecs - (hour * 360000) - min * 6000) / 100;

            (*out)(i, 0) = year;
            (*out)(i, 1) = month;
            (*out)(i, 2) = day;
            (*out)(i, 3) = hour;
            (*out)(i, 4) = min;
            (*out)(i, 5) = sec;
        }
    }
    else if (inputs[0].IsString() || inputs[0].IsCellArray())
    {
        std::vector<std::string> datesToProcess;
        if (inputs[0].IsString())
        {
            datesToProcess.push_back(inputs[0].StringVal());
        }
        else
        {
            HML_CELLARRAY* cl = inputs[0].CellArray();
            for (int i = 0; i < cl->Size(); ++i)
            {
                if ((*cl)(i).IsString())
                {
                    datesToProcess.push_back((*cl)(i).StringVal());
                }
                else
                {
                    throw OML_Error(OML_ERR_STRING_STRINGCELL);
                }
            }
        }

        // on Linux, get_time cannot read 1-digit (without zero padding) dates or times 
        // check each date, if a string has such numbers add the zero 
        // TODO: review when gcc version is updated
#ifndef OS_WIN
        std::vector<Currency> ins;
        ins.push_back("");
        ins.push_back("[^\\d]\\d(?!\\d)");
        std::vector<std::string>::iterator it = datesToProcess.begin();
        for (; it != datesToProcess.end(); ++it)
        {
            std::string cp = " " + std::string(*it);
            ins[0] = cp;
            Currency tmp = eval.CallFunction("regexp", ins);
            
            std::vector<double> fidx;
            if (tmp.IsScalar())
                fidx.push_back(tmp.Scalar());
            else if (tmp.IsVector())
                fidx = tmp.Vector();
            
            std::vector<double>::const_reverse_iterator rit = fidx.crbegin();
            for (; rit != fidx.crend(); ++rit)
            {
                size_t idx = (int)(*rit);
                std::string num = cp.substr(idx, 1);
                num = "0" + num;
                cp.replace(idx, 1, num.c_str());
            }
            (*it) = cp.substr(1, cp.length() - 1);
        }
#endif //OS_WIN

        time_t rawtime = time(nullptr);
        struct tm* curtime = localtime(&rawtime);
        int currentYear = (int)curtime->tm_year + 1900;
        int referenceYear = currentYear - 50;
        std::string userFmt;
        
        if (nargin == 2)
        {
            if (inputs[1].IsScalar())
                referenceYear = (int)inputs[1].Scalar();
            else if (inputs[1].IsString())
                userFmt = inputs[1].StringVal();
            else
                throw OML_Error(OML_ERR_STRING_INTEGER, 2);
        }
        else if (nargin == 3)
        {
            if (!inputs[1].IsString())
                throw OML_Error(OML_ERR_STRING, 2);
            if (!inputs[2].IsScalar())
                throw OML_Error(OML_ERR_INTEGER, 3);
            
            referenceYear = (int)inputs[2].Scalar();
            userFmt = inputs[1].StringVal();
        }

        hwMatrix* mat = EvaluatorInterface::allocateMatrix((int)datesToProcess.size(), 6, 0.0);
        vecMat = Currency(mat); 
        int res = getDateVectorFromString(datesToProcess, userFmt, referenceYear, currentYear, mat);
        if (res == 1)
            throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT);
        else if (res == 2)
            throw OML_Error(OML_ERR_INVALID_DATE_FMT);
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARSTRING);
    }

    int nargout = eval.GetNargoutValue();
    if (nargout > 1)
    {
        const hwMatrix* mat = vecMat.Matrix();

        // output matrices dimensions must the same as the input matrix (if input is matrix)
        int M = mat->M();
        int N = 1;
        if (inputs[0].IsMatrix())
        {
            M = inputs[0].Matrix()->M();
            N = inputs[0].Matrix()->N();
        }
        for (int i = 0; i < nargout && i < 6; ++i)
        {
            hwMatrix* col = EvaluatorInterface::allocateMatrix(M, N, true);
            mat->ReadColumn(i, *col);
            col->Reshape(M, N);
            outputs.push_back(col);
        }
    }
    else
    {
        outputs.push_back(vecMat);
    }
    return true;
}
//------------------------------------------------------------------------------
// Gets the day of the week from a date number or a date string [weekday]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Weekday(EvaluatorInterface           eval,
                               const std::vector<Currency>& inputs, 
                               std::vector<Currency>&       outputs)
{
    int nargout = getNumOutputs(eval);
    int nargin = (int)inputs.size();
    if (nargin < 1 || nargin > 2)
        throw OML_Error(OML_ERR_NUMARGIN);

    Currency out;
    if (inputs[0].IsScalar())
    {
        hwMatrix* mat = EvaluatorInterface::allocateMatrix(1, 1, true);
        (*mat)(0) = inputs[0].Scalar();
        out = Currency(mat);
    }
    else if (inputs[0].IsMatrix())
    {
        hwMatrix* mat = EvaluatorInterface::allocateMatrix(inputs[0].Matrix());
        out = Currency(mat);
    }
    else if (inputs[0].IsString() || inputs[0].IsCellArray())
    {
        // first convert to date vector
        std::vector<Currency> ins, outs;
        ins.push_back(inputs[0]);
        Currency tmp = eval.CallFunction("datevec", ins);
        
        if (!tmp.IsMatrix())
            throw OML_Error(OML_ERR_INVALID_DATE_FMT);

        const hwMatrix* datevecs = tmp.Matrix();
        if (datevecs->N() != 6)
            throw OML_Error(OML_ERR_INVALID_DATE_FMT);

        // initialize the output matrix, if the input is 
        // a cell of strings then the matrix must have the same dimensions
        hwMatrix* mat = inputs[0].IsString() ? EvaluatorInterface::allocateMatrix(1, 1, true):
            EvaluatorInterface::allocateMatrix(inputs[0].CellArray()->M(), inputs[0].CellArray()->N(), true);
        out = Currency(mat);

        for (int i=0;i<datevecs->M();++i)
        {
            ins.clear();
            ins.push_back((*datevecs)(i, 0));
            ins.push_back((*datevecs)(i, 1));
            ins.push_back((*datevecs)(i, 2));
            ins.push_back((*datevecs)(i, 3));
            ins.push_back((*datevecs)(i, 4));
            ins.push_back((*datevecs)(i, 5));
            
            std::vector<Currency> dnout;
            Datenum(eval, ins, dnout);

            if (!dnout.empty() && dnout[0].IsScalar())
                (*mat)(i) = dnout[0].Scalar();
        }
    }
    else
    {
        throw OML_Error(OML_ERR_SCALARSTRING, 1);
    }

    if (!out.IsMatrix())
        throw OML_Error(OML_ERR_INTERNAL);

    bool shortFormat = true;
    if (nargin == 2)
    {
        if (!inputs[1].IsString())
            throw OML_Error(OML_ERR_STRING, 2);

        std::string val = inputs[1].StringVal();
        std::transform(val.begin(), val.end(), val.begin(), ::tolower);
        if (val == "short" || val == "long")
        {
            shortFormat = val == "short";
        }
        else
        {
            throw OML_Error(OML_ERR_OPTIONVAL, 2);
        }
    }

    std::vector<std::string> outNames;
    hwMatrix* mat = out.GetWritableMatrix();
    for (int i = 0; i < mat->Size(); ++i)
    {
        // fixed point - Sunday Aug 7 2022
        double val = (*mat)(i);
        int diff = (int)val - 738740;
        int dayOfWeek = (int)mod(diff, 7) + 1;
        if (dayOfWeek == 0)
            dayOfWeek = 7;
        if (diff >= 0)
        {
            if (dayOfWeek == 0)
                dayOfWeek = 7;
        }
        (*mat)(i) = dayOfWeek;
        if (nargout == 2)
            outNames.push_back(shortFormat ? getDayName(dayOfWeek, true) : getDayName(dayOfWeek, false));
    }

    outputs.push_back(out);
    if (nargout == 2)
    {
        // convert from cell of strings to matrix
        std::vector<Currency> ins, outs;
        ins.push_back(outNames);
        Currency tmp = eval.CallFunction("str2mat", ins);
        outputs.push_back(tmp);
    }
    return true;
}
//------------------------------------------------------------------------------
// Returns the date matrix as a date string using the given dateFmt.
//------------------------------------------------------------------------------
std::vector<std::string> getDateStringFromDateVec(const hwMatrix* dateMat, const std::string& dateFmt,
    const Currency& dayOfWeek)
{
    // check if it's "Q" date format
    int qFmt = 0;
    if (dateFmt == "QQ-YY")
    {
        qFmt = 1;
    }
    else if (dateFmt == "QQ")
    {
        qFmt = 2;
    }
    else if (dateFmt == "QQ-YYYY")
    {
        qFmt = 3;
    }

    // replace format to symbols that strftime recognises
    std::string strftimeFmt(dateFmt);
    bool amTime = false, pmTime = false;
    bool hasMilliseconds = false, hasFirstLMonth = false, hasFirstLDay = false;
    bool has4DigitYear = false, has2DigitYear = false;
    if (qFmt == 0)
    {
        // year
        findDatePatternAndReplace(strftimeFmt, "yyyy", "%Y");
        findDatePatternAndReplace(strftimeFmt, "YYYY", "%Y");
        findDatePatternAndReplace(strftimeFmt, "yy", "%y");
        findDatePatternAndReplace(strftimeFmt, "YY", "%y");
        has4DigitYear = strftimeFmt.find("%Y") != std::string::npos;
        has2DigitYear = strftimeFmt.find("%y") != std::string::npos;


        // month
        findDatePatternAndReplace(strftimeFmt, "mmmm", "%B");
        findDatePatternAndReplace(strftimeFmt, "mmm", "%h");
        findDatePatternAndReplace(strftimeFmt, "mm", "%m");

        // day
        findDatePatternAndReplace(strftimeFmt, "dddd", "%A");
        findDatePatternAndReplace(strftimeFmt, "ddd", "%a");
        findDatePatternAndReplace(strftimeFmt, "dd", "%d");

        // hour
        amTime = strftimeFmt.find("AM") != std::string::npos;
        pmTime = strftimeFmt.find("PM") != std::string::npos;
        if (amTime || pmTime)
            findDatePatternAndReplace(strftimeFmt, "HH", "%I");
        else
            findDatePatternAndReplace(strftimeFmt, "HH", "%H");

        // min
        findDatePatternAndReplace(strftimeFmt, "MM", "%M");
        // sec
        findDatePatternAndReplace(strftimeFmt, "SS", "%S");

        // check for milliseconds
        hasMilliseconds = strftimeFmt.find("FFF") != std::string::npos;

        // check for first letter of month
        size_t s = 0;
        while ((s = strftimeFmt.find("m", s)) != std::string::npos)
        {
            if (s == 0 || !(strftimeFmt[s - 1] == '%'))
            {
                hasFirstLMonth = true;
                break;
            }
            ++s;
        }

        // check for first letter of day
        s = 0;
        while ((s = strftimeFmt.find("d", s)) != std::string::npos)
        {
            if (s == 0 || !(strftimeFmt[s - 1] == '%'))
            {
                hasFirstLDay = true;
                break;
            }
            ++s;
        }
    }

    std::vector<std::string> ret;
    for (int i = 0; i < dateMat->M(); ++i)
    {
        std::string date(strftimeFmt);
        int dW = -1;
        if (dayOfWeek.IsScalar())
            dW = (int)dayOfWeek.Scalar();
        else if (dayOfWeek.IsMatrix())
        {
            const hwMatrix* mt = dayOfWeek.Matrix();
            dW = (int)(*mt)(i);
        }

        double year = (*dateMat)(i, 0);
        int month = (int)(*dateMat)(i, 1);
        int day = (int)(*dateMat)(i, 2);
        int hour = (int)(*dateMat)(i, 3);
        int min = (int)(*dateMat)(i, 4);
        double sec = (*dateMat)(i, 5);

        // "Q" date formats
        if (qFmt == 1)
        {
            int q = (int)(floor((month - 1) / 3) + 1);
            char buf[10];
            sprintf(buf, "Q%d-%02.0f", q, mod(year, 100));
            ret.push_back(std::string(buf));
            continue;
        }
        else if (qFmt == 2)
        {
            int q = (int)(floor((month - 1) / 3) + 1);
            char buf[10];
            sprintf(buf, "Q%d", q);
            ret.push_back(std::string(buf));
            continue;
        }
        else if (qFmt == 3)
        {
            int q = (int)(floor((month - 1) / 3) + 1);
            char buf[10];
            sprintf(buf, "Q%d-%04.0f", q, year);
            ret.push_back(std::string(buf));
            continue;
        }

        // adjust for < 0 years
        if (year < 0)
        {
            findDatePatternAndReplace(date, "%Y", "-%Y");
            findDatePatternAndReplace(date, "%y", "-%y");
        }

        // AM <->PM
        if (amTime && hour >= 12)
            findDatePatternAndReplace(date, "AM", "PM");
        else if (pmTime && hour < 12)
            findDatePatternAndReplace(date, "PM", "AM");

        // milliseconds
        if (hasMilliseconds)
        {
            char buffer[10];
            double milli = round((sec - (int)sec)*1000);
            sprintf(buffer, "%03d", (int)milli);
            findDatePatternAndReplace(date, "FFF", std::string(buffer));
        }

        // replace 'm' with first letter of the month
        if (hasFirstLMonth)
        {
            size_t s = 0;
            std::string monthFirstLetter = getMonthName(month, true);
            monthFirstLetter.resize(1);
            while ((s = date.find("m", s)) != std::string::npos)
            {
                if (s == 0 || !(date[s - 1] == '%'))
                {
                    date.replace(s, 1, monthFirstLetter);
                }
                ++s;
            }
        }

        // replace 'd' with first letter of the day
        if (hasFirstLDay)
        {
            std::string dayFirstLetter = getDayName(dW, true);
            dayFirstLetter.resize(1);
            size_t s = 0;
            while ((s = date.find("d", s)) != std::string::npos)
            {
                if (s == 0 || !(date[s - 1] == '%'))
                {
                    date.replace(s, 1, dayFirstLetter);
                }
                ++s;
            }
        }

        // strftime cannot handle years with more than 4 digits and will crash
        if (year > 9999 || year < -9999)
        {
            if (has2DigitYear || has4DigitYear)
            {
                double y = has4DigitYear ? abs(year) : mod(abs(year),100);
                char tmp[128];
                sprintf(tmp, "%02.0f", y);
                std::string yearStr = std::string(tmp);
                std::string pat = has4DigitYear ? std::string("%Y") : std::string("%y");
                size_t s = 0;
                while ((s = date.find(pat, s)) != std::string::npos)
                {
                    date.replace(s, 2, yearStr);
                    ++s;
                }
            }
            else
            {
                // just in case
                year = 9999;
            }
        }

        tm t = {};
        t.tm_year = (int)abs(year) - 1900;
        t.tm_mon = month - 1;
        t.tm_mday = day;
        t.tm_hour = hour;
        t.tm_min = min;
        t.tm_sec = (int)sec;
        t.tm_wday = dW - 1;

        char buf[128];
        size_t res = std::strftime(buf, 128, date.c_str(), &t);
        if (res == 0)
            throw OML_Error(OML_ERR_INVALID_DATE_FMT);

        ret.push_back(std::string(buf));
    }
    return ret;
}
//------------------------------------------------------------------------------
// Gets the date string from a date number or a date string [datestr]
//------------------------------------------------------------------------------
bool BuiltInFuncsTime::Datestr(EvaluatorInterface           eval,
                               const std::vector<Currency>& inputs,
                               std::vector<Currency>& outputs)
{
    int nargin = (int)inputs.size();
    if (nargin < 1 || nargin > 3)
        throw OML_Error(OML_ERR_NUMARGIN);

    BuiltInFuncsTime funcs;
    funcs.InitializeFormats();
    if (!inputs[0].IsScalar() && !inputs[0].IsMatrix() && !inputs[0].IsString() && !inputs[0].IsCellArray())
        throw OML_Error(OML_ERR_SCALAR_VECTOR_STRING, 1);

    // validate 'fmt' input
    if (nargin > 1 && !inputs[1].IsString() && !inputs[1].IsScalar())
        throw OML_Error(OML_ERR_STRING, 2);
    // validate 'reference year' input
    if (nargin > 2 && !inputs[2].IsScalar())
        throw OML_Error(OML_ERR_INTEGER, 3);

    // first convert to datevec
    std::vector<Currency> ins, outs;
    ins.push_back(inputs[0]);

    if (inputs[0].IsString() || inputs[0].IsCellArray())
    {
        // check for reference year input
        if (nargin > 2)
        {
            ins.push_back(inputs[2].Scalar()); // already validated that it is a scalar
        }
    }
    Datevec(eval, ins, outs);

    if (outs.empty() || !outs[0].IsMatrix())
        throw OML_Error(OML_ERR_INTERNAL);

    const hwMatrix* dateMat = outs[0].Matrix();
    // get the date output format
    std::string dateFormat;
    int dateId = -1;
    if (nargin > 1)
    {
        if (inputs[1].IsString())
            dateFormat = inputs[1].StringVal();
        else if (inputs[1].IsInteger())
            dateId = (int)inputs[1].Scalar();
        else
            throw OML_Error(OML_ERR_INTSTRING, 2);
    }
    if (dateFormat.empty())
    {
        if (dateId >= 0)
        {
            if (dateId < 32)
                dateFormat = funcs._omlDateTimeFmt[dateId];
            else
                throw OML_Error(OML_ERR_INVALID_DATE_FMT);
        }
        else
        {
            // if all vecs have zeros in the time elements use format 1
            // if all vecs have zeros in the date elements use format 16 (only valid for string input, date cannot be 0/0/0)
            // if vecs have both date and time use format 0
            bool timeZero = true;
            bool daysZero = true;
            for (int i = 0; i < dateMat->M(); ++i)
            {
                if (!iszero((*dateMat)(i, 3)) || !iszero((*dateMat)(i, 4)) || !iszero((*dateMat)(i, 5)))
                {
                    timeZero = false;
                }
                if (!iszero((*dateMat)(i, 0)) || !iszero((*dateMat)(i, 1)) || !iszero((*dateMat)(i, 2)))
                {
                    daysZero = false;
                }
            }
            if (timeZero)
                dateFormat = funcs._omlDateTimeFmt[1];
            else if (daysZero)
                dateFormat = funcs._omlDateTimeFmt[16];
            else
                dateFormat = funcs._omlDateTimeFmt[0];
        }
    }

    // get day of the week in case it is needed in date format
    std::vector<Currency> weekdayIn, weekdayOut;
    weekdayIn.push_back(inputs[0]);
    Weekday(eval, weekdayIn, weekdayOut);
    if (weekdayOut.empty() || !weekdayOut[0].IsMatrix())
        throw OML_Error(OML_ERR_INTERNAL);

    std::vector<Currency> str2MatOut;
    ins.clear();
    // create the date strings and convert cell of strings to matrix
    ins.push_back(getDateStringFromDateVec(dateMat, dateFormat, weekdayOut[0]));
    Currency tmp = eval.CallFunction("str2mat", ins);
    outputs.push_back(tmp);
    return true;
}