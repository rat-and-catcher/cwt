/*
 * hilb_fir.c -- the imlementation file for
 * the function to compute Hilbert transform FIR filter coeff.
 * for the spectrum shifter program.
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
 * (c) 1999, published by Pearson Education, Inc, publishing
 * as Prentice Hall
 * ISBN 0-13-754920-2
 */

#include "hilb_fir.h"

#ifdef PI
#undef PI
#endif
#define PI      (3.141592653589793238462643383279502884197169399375105820974944)

/*
 * Static helpers
 */

/* Square of double
*/
static double sqrd(double x)
{
 return x * x;
}

/*
 * Front end functions
 */
/* compute n-th FIR Hilbert filter coefficient, n = [0..M]
*/
double FilterCoeff(int n, int M, double beta)
{
 // the Book's example 11.4
 int n_delay = M / 2;
 double n_nd;

 if(n > M || n < 0)             // unreachable
  return 0.0;

 n_nd = (double)(n - n_delay);

 return
        // Kaiser's window part
        (besselI0(beta * sqrt(1.0 - sqrd(n_nd / ((double)n_delay)))) /
        besselI0(beta)) *
        // the "true" pulse response part
#if 0
        // not so good
        (n != n_delay? 2.0 * sqrd(sin(PI * n_nd / 2.0)) / (PI * n_nd) : 0.0);
#else
        // very good
        (((n - n_delay) & 1)? 2.0 / (PI * n_nd) : 0.0);
#endif
}

/* the end...
*/

