/*
 * Copyright (c) 2026 Rat and Catcher Tech.
 *
 * Original FFTW3 infrastructure taken from fftw3.h:
 * Copyright (c) 2003, 2007-14 Matteo Frigo
 * Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
 *
 * The following statement of license applies *only* to this header file,
 * and *not* to the other files distributed with FFTW or derived therefrom:
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#if !defined(FFTW3_RC_H_)
#define FFTW3_RC_H_

#if defined(__cplusplus)
extern "C"
{
#endif /* __cplusplus */

/* The only docs for exported function(s)
 * --------------------------------------
 * -- get/set assert behavior mode:
 * int fftw*_rc_fail_behavior_gs)(int mode);
 * mode < 0 -- get current mode (=0/>0);
 * mode = 0 -- set fftw3 native mode (default, abort() when assert failed);
 * mode > 0 -- set memory exception when assert failed;
 */

/* we can't to include proper "fftw3.h" in the all cases of using this file
*/
#if !defined(FFTW_MANGLE_DOUBLE)
#error "You should to include original fftw3.h before fftw-rc.h (this file)"
#endif

#define FFTW_RC_DEFINE_API(X, R, C)                                     \
FFTW_EXTERN int                                                         \
FFTW_CDECL X(rc_fail_behavior_gs)(int mode);

FFTW_RC_DEFINE_API(FFTW_MANGLE_DOUBLE, double, fftw_complex)
FFTW_RC_DEFINE_API(FFTW_MANGLE_FLOAT, float, fftwf_complex)
FFTW_RC_DEFINE_API(FFTW_MANGLE_LONG_DOUBLE, long double, fftwl_complex)

/* __float128 (quad precision) is a gcc extension on i386, x86_64, and ia64
   for gcc >= 4.6 (compiled in FFTW with --enable-quad-precision) */
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 6)) \
 && !(defined(__ICC) || defined(__INTEL_COMPILER) || defined(__CUDACC__) || defined(__PGI)) \
 && (defined(__i386__) || defined(__x86_64__) || defined(__ia64__))
#  if !defined(FFTW_NO_Complex) && defined(_Complex_I) && defined(complex) && defined(I)
/* note: __float128 is a typedef, which is not supported with the _Complex
         keyword in gcc, so instead we use this ugly __attribute__ version.
         However, we can't simply pass the __attribute__ version to
         FFTW_DEFINE_API because the __attribute__ confuses gcc in pointer
         types.  Hence redefining FFTW_DEFINE_COMPLEX.  Ugh. */
#    undef FFTW_DEFINE_COMPLEX
#    define FFTW_DEFINE_COMPLEX(R, C) typedef _Complex float __attribute__((mode(TC))) C
#  endif
FFTW_RC_DEFINE_API(FFTW_MANGLE_QUAD, __float128, fftwq_complex)
#endif 

#if defined(__cplusplus)
}                       /* extern "C" */
#endif                  /* __cplusplus */

#endif                  /* FFTW3_RC_H_ */

