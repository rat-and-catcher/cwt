/*
 * besseli0.c -- computation of modified Bessel function of order zero
 *
 * This file extracted from ALGLIB project. ALGLIB copyright notice
 * and conditions keep unchanged
 *
 * -- The source modified for the C-programming conventions
 * -- by Rat and Catcher Tech.; all rights reserved
 */

/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to
      pseudocode.

See subroutines comments for additional copyrights.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.

- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*************************************************************************/

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <math.h>
#include "besseli0.h"


/*************************************************************************
Internal subroutine;
-- The source modified for the C-programming conventions
-- by Rat and Catcher Tech.; all rights reserved

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
#define besselmfirstcheb(c, b0, b1, b2)         \
        do                                      \
        {                                       \
         (b0) = (c);                            \
         (b1) = 0.0;                            \
         (b2) = 0.0;                            \
        }                                       \
        while(0)


/*************************************************************************
Internal subroutine;
-- The source modified for the C-programming conventions
-- by Rat and Catcher Tech.; all rights reserved

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
#define besselmnextcheb(x, c, b0, b1, b2)       \
        do                                      \
        {                                       \
         (b2) = (b1);                           \
         (b1) = (b0);                           \
         (b0) = (x) * (b1) - (b2) + (c);        \
        }                                       \
        while(0)


/*************************************************************************
Modified Bessel function of order zero

Returns modified Bessel function of order zero of the
argument.

The function is defined as i0(x) = j0( ix ).

The range is partitioned into the two intervals [0,8] and
(8, infinity).  Chebyshev polynomial expansions are employed
in each interval.

ACCURACY:

                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE      0,30        30000       5.8e-16     1.4e-16

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 2000 by Stephen L. Moshier
*************************************************************************/
double besselI0(double x)
{
 double y;
 double v;
 double z;
 double b0;
 double b1;
 double b2;

 if(x < 0)
 {
  x = -x;
 }
 if(x <= 8.0)
 {
  y = x / 2.0 - 2.0;
  besselmfirstcheb(-4.41534164647933937950E-18, b0, b1, b2);
  besselmnextcheb(y, 3.33079451882223809783E-17, b0, b1, b2);
  besselmnextcheb(y, -2.43127984654795469359E-16, b0, b1, b2);
  besselmnextcheb(y, 1.71539128555513303061E-15, b0, b1, b2);
  besselmnextcheb(y, -1.16853328779934516808E-14, b0, b1, b2);
  besselmnextcheb(y, 7.67618549860493561688E-14, b0, b1, b2);
  besselmnextcheb(y, -4.85644678311192946090E-13, b0, b1, b2);
  besselmnextcheb(y, 2.95505266312963983461E-12, b0, b1, b2);
  besselmnextcheb(y, -1.72682629144155570723E-11, b0, b1, b2);
  besselmnextcheb(y, 9.67580903537323691224E-11, b0, b1, b2);
  besselmnextcheb(y, -5.18979560163526290666E-10, b0, b1, b2);
  besselmnextcheb(y, 2.65982372468238665035E-9, b0, b1, b2);
  besselmnextcheb(y, -1.30002500998624804212E-8, b0, b1, b2);
  besselmnextcheb(y, 6.04699502254191894932E-8, b0, b1, b2);
  besselmnextcheb(y, -2.67079385394061173391E-7, b0, b1, b2);
  besselmnextcheb(y, 1.11738753912010371815E-6, b0, b1, b2);
  besselmnextcheb(y, -4.41673835845875056359E-6, b0, b1, b2);
  besselmnextcheb(y, 1.64484480707288970893E-5, b0, b1, b2);
  besselmnextcheb(y, -5.75419501008210370398E-5, b0, b1, b2);
  besselmnextcheb(y, 1.88502885095841655729E-4, b0, b1, b2);
  besselmnextcheb(y, -5.76375574538582365885E-4, b0, b1, b2);
  besselmnextcheb(y, 1.63947561694133579842E-3, b0, b1, b2);
  besselmnextcheb(y, -4.32430999505057594430E-3, b0, b1, b2);
  besselmnextcheb(y, 1.05464603945949983183E-2, b0, b1, b2);
  besselmnextcheb(y, -2.37374148058994688156E-2, b0, b1, b2);
  besselmnextcheb(y, 4.93052842396707084878E-2, b0, b1, b2);
  besselmnextcheb(y, -9.49010970480476444210E-2, b0, b1, b2);
  besselmnextcheb(y, 1.71620901522208775349E-1, b0, b1, b2);
  besselmnextcheb(y, -3.04682672343198398683E-1, b0, b1, b2);
  besselmnextcheb(y, 6.76795274409476084995E-1, b0, b1, b2);
  v = 0.5 * (b0 - b2);
  return (exp(x) * v);
 }

 z = 32.0 / x - 2.0;
 besselmfirstcheb(-7.23318048787475395456E-18, b0, b1, b2);
 besselmnextcheb(z, -4.83050448594418207126E-18, b0, b1, b2);
 besselmnextcheb(z, 4.46562142029675999901E-17, b0, b1, b2);
 besselmnextcheb(z, 3.46122286769746109310E-17, b0, b1, b2);
 besselmnextcheb(z, -2.82762398051658348494E-16, b0, b1, b2);
 besselmnextcheb(z, -3.42548561967721913462E-16, b0, b1, b2);
 besselmnextcheb(z, 1.77256013305652638360E-15, b0, b1, b2);
 besselmnextcheb(z, 3.81168066935262242075E-15, b0, b1, b2);
 besselmnextcheb(z, -9.55484669882830764870E-15, b0, b1, b2);
 besselmnextcheb(z, -4.15056934728722208663E-14, b0, b1, b2);
 besselmnextcheb(z, 1.54008621752140982691E-14, b0, b1, b2);
 besselmnextcheb(z, 3.85277838274214270114E-13, b0, b1, b2);
 besselmnextcheb(z, 7.18012445138366623367E-13, b0, b1, b2);
 besselmnextcheb(z, -1.79417853150680611778E-12, b0, b1, b2);
 besselmnextcheb(z, -1.32158118404477131188E-11, b0, b1, b2);
 besselmnextcheb(z, -3.14991652796324136454E-11, b0, b1, b2);
 besselmnextcheb(z, 1.18891471078464383424E-11, b0, b1, b2);
 besselmnextcheb(z, 4.94060238822496958910E-10, b0, b1, b2);
 besselmnextcheb(z, 3.39623202570838634515E-9, b0, b1, b2);
 besselmnextcheb(z, 2.26666899049817806459E-8, b0, b1, b2);
 besselmnextcheb(z, 2.04891858946906374183E-7, b0, b1, b2);
 besselmnextcheb(z, 2.89137052083475648297E-6, b0, b1, b2);
 besselmnextcheb(z, 6.88975834691682398426E-5, b0, b1, b2);
 besselmnextcheb(z, 3.36911647825569408990E-3, b0, b1, b2);
 besselmnextcheb(z, 8.04490411014108831608E-1, b0, b1, b2);
 v = 0.5 * (b0 - b2);
 return (exp(x) * v / sqrt(x));
}

/* the end...
*/


