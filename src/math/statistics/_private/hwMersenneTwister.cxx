/**
* @file hwMersenneTwister.cxx
* @date June 2007
* Copyright (C) 2007-2018 Altair Engineering, Inc.  
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
#include "hwMersenneTwister.h"

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwMersenneTwisterState::hwMersenneTwisterState()
    : m_initialized (false)
{
}
//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
hwMersenneTwisterState::hwMersenneTwisterState(unsigned long seed)
{
    Initialize(seed);
}
//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
hwMersenneTwisterState::~hwMersenneTwisterState()
{
}

#define N 624
#define M 397
#define MATRIX_A 2567483615UL       // constant vector a, originally 0x9908b0dfUL
#define UPPER_MASK 2147483648UL     // most significant w-r bits, originally 0x80000000UL
#define LOWER_MASK 2147483647UL     // least significant r bits, originally 0x7fffffffUL

//------------------------------------------------------------------------------
// Initialize the random number generator
//------------------------------------------------------------------------------
void hwMersenneTwisterState::Initialize(unsigned long seed)
{
    mt[0]= seed & 4294967295UL;     // originally 0xffffffffUL

    for (mti=1; mti<N; mti++)
    {
        // See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. In the previous versions, MSBs of
        // the seed affect only MSBs of the array mt[]. 2002/01/09 modified by Makoto Matsumoto.
        mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 

        // for >32 bit machines
        mt[mti] &= 4294967295UL;
    }

    m_initialized = true;
}
//------------------------------------------------------------------------------
// Compute a random deviate from the distribution
//------------------------------------------------------------------------------
unsigned long hwMersenneTwisterState::NextValue()
{
    static unsigned long mag01[2] = {0UL, MATRIX_A};    // mag01[x] = x * MATRIX_A  for x=0,1
    unsigned long y;

    if (mti >= N)                           // generate N words at one time
    {
        int kk;

        // if (mti == N+1)                  // if init_genrand() has not been called,
        //     init_genrand(5489UL);        // a default initial seed is used

        for (kk = 0; kk < N-M; kk++)
        {
            y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 1UL];
        }

        for (; kk<N-1; kk++)
        {
            y = (mt[kk] & UPPER_MASK)|(mt[kk+1] & LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 1UL];
        }

        y = (mt[N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 1UL];

        mti = 0;
    }

    y = mt[mti++];

    // Tempering
    y ^= (y >> 11);
    y ^= (y << 7) & 2636928640UL;       // originally 0x9d2c5680UL
    y ^= (y << 15) & 4022730752UL;      // originally 0xefc60000UL
    y ^= (y >> 18);

    return y;
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK
