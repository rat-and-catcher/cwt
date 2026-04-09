/*
 * Copyright (c) 2026 Rat and Catcher Tech.
 * Some FFTW3's API addons
 *
 * Original FFTW3 infrastructure taken from fftw3.h:
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
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

#include "api/api.h"                /* fftw3.h defs -- here. */
#include "./fftw3-rc.h"             /* R&C API */
#include "../kernel/ifftw-rc.h"     /* R&C internals */

/* get / set assertion_failed() behavior -- frontend
*/
int X(rc_fail_behavior_gs)(int mode)
{
    return X(ratcat_fail_behavior_gs)(mode);
}

