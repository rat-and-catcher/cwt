/*
 * cwt.h -- the main declarationthe file posix-gcc-only tiny
 * 2-ch wav to atalitic (complex) signal transformation;
 * (direct FFT-based transform only)
 * This program can be distributed under GNU GPL
 *
 * Copyright (C) 2010-2012-2014 Rat and Catcher Tech.
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
 * The code in THIS uses functions from FFTW library, see www.fftw.org
 */
/*
 * HISTORY
 * Version V2.0.0 20-may-2014
 * -- rewrite Windows cw for posix, remove FIR-related code/switches/etc;
 * Version V2.0.1 26-may-2015
 * -- fix makefile to force static FFTW libs;
 */

#ifndef _cw_h_
#define _cw_h_

#ifndef	_CRT_SECURE_NO_WARNINGS		/* MS VC stuff, it's safe to leave this here */
#define	_CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#include "crc32.h"

#include <fftw3.h>			/* out of the project tree now from V1.0.6 */

#include "cwave.h"

// version
#define	VERSION		("V2.0.1")	/* program version */

// defaults
#define DEF_GAIN	(1.0)		/* gain multiplier */

/*
 * Compiler depended stuff
 * -------- -------- -----
 */
#define INLINE		inline		/* posix */
#define ISATTY		isatty
#define FILENO		fileno
#define	STRICMP		strcasecmp

#if 0			// "true" gcc C-exceptions (??)
#define	TRY		try
#define	CATCH		catch
#define	CATCH_ALL	catch(...)
#else			// fake exceptions -- old gcc
#define	TRY		if(1)
#define	CATCH		if(0)
#define	CATCH_ALL	if(0)
#endif

/*
 * The main global application object
 * --- ---- ------ ----------- ------
 */
typedef struct tagAPPLICATION_CW
{
// files
 char *nameif;			// name of input file
 char *nameof;			// name of output file
// general program parameters
 double gain_mul;		// gain multiplier
 int isTestCRC;			// !=0 -> test mode of input CWAVE
 int isVerbose;			// !=0 -> progress output, and, probably, more
 unsigned c_format;		// CWAVE format variation
 int nCPU;			// Number of CPU's (max threads for FFTW)
// FFT specific parameters
 int isFFTeven;			// 1 - even number of points in FFT, (0 == odd)
 int isFFTsafe;			// 1 - safe but slow FFT
 int isFFTstat;			// 1 - print FFTW plan statistics
 int isPlanOut;			// 1 - write FFTW plans to stdout
 int isFFTnoSIMD;		// 1 - don't use SIMD (e.g. SSE2) instructions for FFTW
 unsigned nsFFT;		// _real_ number of FFT samples (strictly odd or strictly even)
 double lo_band;		// low frequency to pass
 double hi_band;		// high frequency to pass
 unsigned lo_rem;		// # of low spectral components to remove [0..lo_rem]
 unsigned hi_rem;		// # of high spectral components to remove [hi_rem, N_SAMPLES/2)
// misc.
 double progress;		// processing progress, %
 time_t tstart;			// startup time
 FILE *fpif;			// input file descriptor
 FILE *fpof;			// output file descriptor
 TMP_CRC32 tcrc;		// CRC32 stuff
 long l_clips;			// detected clips for the left channel
 long r_clips;			// detected clips for the right channel
 HCWAVE hcw;			// complex wave header
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

#endif			// _cw_h_

/* the end...
*/

