/*
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * Some has made by Rat and Catcher Tech. for internal needs only.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "kernel/ifftw.h"
#include "./ifftw-rc.h"
#include <stdio.h>
#include <stdlib.h>

/* +added by Rat and Catcher Tech. */
/* get / set assertion_failed() behavior -- internal
*/
int X(ratcat_fail_behavior_gs)(int mode)
{
    static int cur_mode = 0;

    if(mode >= 0)
        cur_mode = mode;

    return cur_mode;
}
/* -added by Rat and Catcher Tech. */

void X(assertion_failed)(const char *s, int line, const char *file)
{
/* +added by Rat and Catcher Tech. */
#if 1				        /* make an exception */
    volatile char *badptr__ = NULL;
    
    if(X(ratcat_fail_behavior_gs)(-1))
    {
     if(*badptr__)          /* exceptionally lazy solution */
      exit(EXIT_FAILURE);	/* unreachable */
    }
#endif
/* -added by Rat and Catcher Tech. */

     fflush(stdout);
     fprintf(stderr, "fftw: %s:%d: assertion failed: %s\n", file, line, s);
#ifdef HAVE_ABORT
     abort();
#else
     exit(EXIT_FAILURE);
#endif
}
