/*
 * hilb_fir.h -- the declaration file for
 * the function to compute Hilbert transform FIR filter coeff.
 *
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
 */

/* the Kaiser's window approach based at:
 * A. Oppenheim, R. Shaffer
 * Discrete Time Signal Processing,
 * Russian translation
 * (c) 1999, published by Pearson Education, Inc, publishing
 * as Prentice Hall
 * ISBN 0-13-754920-2
*/

#ifndef _hilb_fir_h_
#define _hilb_fir_h_

#include <math.h>
#include "besseli0.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Front end function(s)
 */

/* compute n-th FIR Hilbert filter coefficient, n = [0..M]
*/
double FilterCoeff(int n, int M, double beta);

#ifdef __cplusplus
}
#endif

#endif          // def _hilb_fir_h_

/* the end...
*/

