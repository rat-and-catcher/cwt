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
 * For each CRC computation, the code must provide
 * the temporary variable ('tcrc'), which initialized
 * by the crc32Init() and updated by crc32update().
 * The final CRC computed by value of tcrc by crc32final.
 */

#ifndef _crc32_h_
#define _crc32_h_

#ifdef __cplusplus
extern "C" {
#endif

/* the polynomial
*/
#define POLYNOMIAL	(0xEDB88320)

/* the struct for temporary data
*/
typedef struct tagTmpCrc32
{
 unsigned xOr;			// inversion mask
 unsigned temp;			// temporary value for CRC
} TMP_CRC32;

/* initialization of computation
*/
void crc32init(TMP_CRC32 *t);

/* update CRC by block of data
*/
void crc32update(const void *data, unsigned len, TMP_CRC32 *t);

/* return final CRC by last 'temp'
*/
unsigned crc32final(TMP_CRC32 *t);


#ifdef __cplusplus
}
#endif

#endif		// def _crc32_h_

/* the end...
*/

