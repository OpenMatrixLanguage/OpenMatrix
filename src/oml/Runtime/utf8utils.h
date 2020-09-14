/**
* @file utf8utils.h
* @date August 2019
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

#ifndef __utf8utils_h
#define __utf8utils_h

#ifndef OS_WIN
	#include <cstddef>
#endif 

size_t utf8_get_char_size(unsigned char* ptr);
unsigned char* utf8_increment_pointer(unsigned char* ptr);
size_t utf8_strlen(unsigned char* ptr);
size_t utf8_byte_position_from_index(unsigned char* ptr, size_t index);

size_t utf8_get_char_size(const double* ptr);
const double* utf8_increment_pointer(const double* ptr);
size_t utf8_strlen(const double* ptr);
size_t utf8_byte_position_from_index(const double* ptr, size_t index);
#endif
