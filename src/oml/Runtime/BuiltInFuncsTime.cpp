/**
* @file BuiltInFuncsTime.cpp
* @date April 2019
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

#include "BuiltInFuncsTime.h"

#include <cassert>

#include "OML_Error.h"

std::vector<std::string> BuiltInFuncsTime::_omlDateTimeFmt;  // Formats
std::vector<std::string> BuiltInFuncsTime::_omlMonth;        // Month list

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

        // mmmyy
        case 12: str = "%3s%*d"; isstring = true; break;

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
            str = "%d:%*d:%*d %*s";
            break;

        case 15:  // HH:MM
            str = "%d:%*d";
            break;

        case 16:  // HH:MM PM
            str = "%d:%*d %*s";
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
        result = std::sscanf(in.c_str(), str.c_str(), &value);
    }
    if (result < 1 || value < 0 || value > 24)
    {
        throw OML_Error(OML_ERR_CANNOTAPPLY_DATE_FMT, 2);
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