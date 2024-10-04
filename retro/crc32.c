/*
 * Calculations of standard CRC32
 * Copyright (C) 2010 Rat and Catcher Tech.
 *
 *  "This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details."
 *
 *   See the file COPYING for details of this license.
 *
 * Some code fragments based by the ALGLIB project (see comments in besseli0.* sources).
 * Both contributed files was slightly modified by the author.
 *
 * For each CRC computation, the code must provide
 * the temporary variable ('tcrc'), which initialized
 * by the crc32Init() and updated by crc32update().
 * The final CRC computed by value of tcrc by crc32final.
 */

#include "crc32.h"

// the  table for CRC32 computations
static unsigned crc32table_[256] =
{
 0, 0                                           // not initialized
};

/* helper -- initialization of CRC32 computation
*/
static unsigned ini_(unsigned temp)
{
 temp = temp >> 8 ^ crc32table_[(unsigned char)temp];
 temp = temp >> 8 ^ crc32table_[(unsigned char)temp];
 temp = temp >> 8 ^ crc32table_[(unsigned char)temp];
 temp = temp >> 8 ^ crc32table_[(unsigned char)temp];
 return temp;
}

/* initialization of computation
*/
void crc32init(TMP_CRC32 *t)
{
// the table
 if(!crc32table_[1])
 {
  int i, j;
  unsigned temp;

  for(i = 0; i < 256; ++i)
  {
   temp = i;
   for(j = 0; j < 8; ++j)
   {
    temp = temp >> 1 ^ -(int)(temp & 1) & POLYNOMIAL;
   }
   crc32table_[i] = temp;
  }
 }

// initial values
 t -> xOr = ~0;
 t -> temp = 0;
}

/* update CRC by block of data
*/
void crc32update(const void *data, unsigned len, TMP_CRC32 *t)
{
 unsigned char *pdata = (unsigned char *)data;

// 1-st 4 bytes take with inversion. If data portion less 4, it prolongated by zeros
 for ( ; t -> xOr && len; --len, ++pdata)
 {
  t -> temp = t -> temp >> 8 | ~*pdata << 24;
  t -> xOr >>= 8;
  if(0 == t -> xOr)
  {
   t -> temp  = ini_(t -> temp);
  }
 }
 
// "generic" CRC32 computation
 while(len)
 {
  t -> temp = crc32table_[(unsigned char)t -> temp ^ *pdata++] ^ (t -> temp >> 8);
  --len;
 }
}

/* return final CRC by last 'temp'
*/
unsigned crc32final(TMP_CRC32 *t)
{
 return ~(t -> xOr? t -> xOr ^ ini_(t -> temp) : t -> temp);
}

/* the end...
*/


