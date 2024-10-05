/*
 * cw.h -- the main declaration file windows-only
 * 2-ch wav to atalitic (complex) signal transformation;
 * (FIR based Hilbert transform + direct FFT-based transform)
 * This program can be distributed under GNU GPL
 *
 * Copyright (C) 2010-2012 Rat and Catcher Tech.
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
 * The code in THIS file also used functions from FFTW library, see www.fftw.org
 */
/*
 * HISTORY
 * Version V1.0.0 28-sep-2010
 * -- initial workable version; input 2ch 16bit only:
 * (-h, -M, -b, -g, -p, -t, -v, -c, -i[d|s|m], -1, -2)
 * Version V1.0.1 14-oct-2010
 * -- -f -- direct conversion via FFT(W), initial attempt
 * Version V1.0.2 21-oct-2010
 * -- -fs/"memory safe" FFT transform added; CWAVE format HCW_FMT_PCM_FLT32
 * added; some rearrange of TheProcess() and related code; fix truncate to round
 * for real parts of FFT results.
 * Version V1.0.3 5-nov-2010
 * -- -v[s][p] flags added. Some copyright notes corrected. Fix checking FFTW plans
 * before usage. -f[i]/-fs[i] added. Some mutithread FIR coding style improved.
 * MAX_SAMPLES_PER_TIME now 512 (256 was); mutltithread usage for FFTW;
 * fix writing temporary filel >= 4GB
 * Version V1.0.4 7-dec-2010
 * -- -rX -- FFT filtering added. Some FFT code rearrange; fix format for "Bad CRC"
 * message
 * Version V1.0.5 9-jan-2011
 * -- some fixes in FFT filtering as "multiband" approach
 * Version V1.0.6 6-jan-2012
 * -- move to VC2010 + FFTW 3.3. SIMD (SSE2) code. Native + MS (DLL +
 * static) == 12 exe's
 * Version V1.0.7 14-feb-2012   -- BUGGY - VOIDED!!
 * -- fix(es) (lots of) FFT processing (for even N to odd by default), -fe switch
 * Version V1.0.8 22-feb-2012   -- BUGGY - VOIDED!!
 * -- handle FFTW exceptions (actual fo modified FFTW only)
 * Version V1.0.9 12-apr-2012
 * -- major bug fix @processSpectrun().cw_fft.c (V1.0.7 and V1.0.8)
 * Version V1.1.0 2-oct-2012
 * fix some comments / messages / copyright notes only
 * Version V1.1.1 31-mar-2015
 * Move to VC2013; some x64 project fixes; some comments / messages fixes
 * Version V1.1.2 21-jun-2016
 * Some project fixes; FFTW3.3.4 added; -q; -r21; FFTW version info
 * Version V1.1.3 16-oct-2017
 * minor fixes in comments/messages
 * Version V1.1.4 05-0ct-2024
 * Fix some misudersood in NULL return value of FFTW planners
 */

#ifndef _cw_h_
#define _cw_h_

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <windows.h>
#include <process.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <io.h>

#include "hilb_fir.h"
#include "crc32.h"

#if !defined(FFTW334)
#include <fftw3.h>                      /* out of the project tree now from V1.0.6 */
#else
#include <fftw334.h>                    /* fftw3.h from FFTW 3.3.4 */
#endif

#include "cwave.h"

// version
#define VERSION         ("V1.1.4")      /* program version */

#if defined(_WIN64)                     /* platform name */
#define PLATFORM        ("x64")
#else
#define PLATFORM        ("x86")
#endif

#if defined(NDEBUG)                     /* debug/release */
#define TYPE_BUILD      ("Release")
#else
#define TYPE_BUILD      ("Debug")
#endif

#if defined(FFTW_DLL)                   /* DLL/static */
#define TYPE_LIB        ("DLL")
#else
#define TYPE_LIB        ("Static")
#endif


// defaults
#define DEF_M           (140002)        /* filter order */
#define DEF_BETA        (2.629)         /* Kaiser param */
#define DEF_GAIN        (1.0)           /* gain multiplier */

// the type(s)
// one channel data
typedef struct tagLRCH
{                               // (indexes are common)
 // values
 double s_in;                   // the input sample
 double s_real;                 // the current real ("I") sample
 double s_image;                // the current image ("Q") sample
 // buffer(s)
 double *s_buf;                 // sample buffer [k_M]
 // indexes
 int ixCur;                     // X[i] for current input/output data
 int ixOutRe;                   // X[i] for output Real ("I") (delayed) part of data
} LRCH;

// one channel thread context
#define MAX_SAMPLES_PER_TIME    (512)           /* size of samples buffers */
typedef struct tagTCH
{
 LRCH ch;                                       // channel descriptor
 double indata[MAX_SAMPLES_PER_TIME];           // input data
 double re_outdata[MAX_SAMPLES_PER_TIME];       // output data (real/I)
 double im_outdata[MAX_SAMPLES_PER_TIME];       // output data (image/Q)
 volatile int nsdata;                           // real number of samples
 volatile int isEnd;                            // != 0 -> EOF
 HANDLE waitIn;                                 // event for waiting input data
 HANDLE waitOut;                                // event for waiting the processing
} TCH;

/*
 * Compiler depended stuff
 * -------- -------- -----
 */
#define INLINE  __inline        /* MSVC */
#define ISATTY  _isatty
#define FILENO  _fileno

/*
 * The main global application object
 * --- ---- ------ ----------- ------
 */
typedef struct tagAPPLICATION_CW
{
// files
 char *nameif;                  // name of input file
 char *nameof;                  // name of output file
// FIR-specific conversion parameters
 int k_M;                       // Kaiser's window parameter, even
 double k_beta;                 // Kaiser's window parameter
 int n_delay;                   // group delay of Hilbert filter
 double *hpr;                   // Hilbert pulse response
// general program parameters
 double gain_mul;               // gain multiplier
 int isPrintH;                  // !=0 -> print Hilbert's FIR coefficients
 int isTestCRC;                 // !=0 -> test mode of input CWAVE
 int isVerbose;                 // !=0 -> progress output, and, probably, more
 int notCompHeadTail;           // !=0 -> do not compensate head/tail delay shift
 int isDirectScan;              // !=0 -> direct memory scan for cache optimization
 unsigned c_format;             // CWAVE format variation
 int nThr;                      // 1 - single thread, 2 - (two) threads, def == #CPU's in system
// FFT specific parameters
 int isFFT;                     // 1 - conversion via FFT(W)
 int isFFTeven;                 // 1 - even number of points in FFT, (0 == odd)
 int isFFTsafe;                 // 1 - safe but slow FFT
 int isFFTstat;                 // 1 - print FFTW plan statistics
 int isPlanOut;                 // 1 - write FFTW plans to stdout
 int isFFTnoSIMD;               // 1 - don't use SIMD (e.g. SSE2) instructions for FFTW
 unsigned nsFFT;                // _real_ number of FFT samples (strictly odd or strictly even)
 double lo_band;                // low frequency to pass
 double hi_band;                // high frequency to pass
 unsigned lo_rem;               // # of low spectral components to remove [0..lo_rem]
 unsigned hi_rem;               // # of high spectral components to remove [hi_rem, N_SAMPLES/2)
// misc.
 double progress;               // processing progress, %
 time_t tstart;                 // startup time
 FILE *fpif;                    // input file descriptor
 FILE *fpof;                    // output file descriptor
 TMP_CRC32 tcrc;                // CRC32 stuff
 long l_clips;                  // detected clips for the left channel
 long r_clips;                  // detected clips for the right channel
 void (*MakeHilbert)(LRCH *ch); // active method for the FIR-based Hilbert transform
 HCWAVE hcw;                    // complex wave header
} APPLICATION_CW;

/*
 * The globals
 */
// @cw.c
// -- variables
extern APPLICATION_CW app;
// -- functions

// @helpers.c
// -- variables
// -- functions
/* print an error message and terminate (no cleanup)
*/
void error(const char *fmt, ...);
/* open file with check
*/
FILE *cfopen(const char *name, const char *mode, const char *msg);
/* allocate memory with check
*/
void *cmalloc(size_t size, const char *msg);
/* make ftell() with check
*/
long cftell(FILE *fp);
/* make fseek() to the position with check
*/
void cfseek(FILE *fp, long pos);
/* create temporary file name (must be free())
*/
char *ctempfile(void);
/* check file extension
*/
int checkFileExt(const char *fname, const char *ext);
/* read and check WAV PCM header
*/
void readWavHeader(FILE *fp, unsigned *srate, unsigned *nsamples);
/* read and convert to double one sample
*/
void readWavSample(FILE *fp, double *ls, double *rs);
/* read complex wave (CWAWE) header
*/
void readCwaveHeader(FILE *fp, HCWAVE *hcw);
/* write complex wave (CWAWE) header
*/
void writeCwaveHeader(FILE *fp, HCWAVE *hcw);
/* read complex data (for CRC check only)
*/
void readComplex(FILE *fp, unsigned t_format, TMP_CRC32 *tcrc);
/* write complex data
*/
void writeComplex(FILE *fp, double l_re, double l_im,
          double r_re, double r_im, const HCWAVE *hcw,
          long *l_clips, long *r_clips, TMP_CRC32 *tcrc);
/* print CWAVE format information
*/
void PrintCwaveFormat(unsigned cw_format);
/* get number of CPU in the system
*/
int GetCpuNumber(void);
/* create auto non-signaling event with check
*/
HANDLE ccreate_event(void);
/* calculate processing time
*/
void CalcTime(time_t tstart);

// @cw_fft.c
// -- variables
// -- functions
/* calculate FFT pass-filter parameters (if any)
*/
void CalcFFTpass(void);
/* The FFT based processing function (unsafe version)
*/
void ProcessFFT(void);
/* The FFT based processing function (safe version)
*/
void ProcessFFT_Safe(void);

#endif                  // _cw_h_

/* the end...
*/

