/*
 * cw-fft.c -- the FFT implementation file for the
 * 2-ch wav to atalitic (complex) signal transformation (tiny / posix);
 * This program can be distributed under GNU GPL
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
 * The code in THIS file uses functions from FFTW library, see www.fftw.org
 */

#include "cwt.h"

/*
 * Forward declarations
 */
/* single channel processing
*/
static void ProcessSingleChFFT(fftw_complex *data,
	unsigned plflags, const char *nmCh);
/* make spectrum processing for one channel
*/
static void processSpectrum(fftw_complex *data);
/* check a plan for existence
*/
static void chkplan(const fftw_plan plan);
/* execute DFT with exception check (for modified FFTW only!)
*/
void myFFTWexecute(const fftw_plan plan);
/* print plan stuff
*/
static void printPlanStat(const char *plName, const fftw_plan plan);

/* calculate FFT pass-filter parameters (if any)
*/
void CalcFFTpass(void)
{
 if(app.lo_band >= 0.0 || app.hi_band > 0.0)
 {
  // (NOT n_samples) nsFFT =::= 2*pi =::= Fs
  double fres = (double)app.nsFFT /*.hcw.n_samples*/ /
	(double)app.hcw.sample_rate;

  if(app.lo_band >= (double)app.hcw.sample_rate / 2.0)
  {
   printf("Warning: low filter pass frequency too high\n");
   app.lo_band = -1.0;
  }
  if(app.hi_band >= (double)app.hcw.sample_rate / 2.0)
  {
   printf("Warning: high filter pass frequency too high\n");
   app.hi_band = -1.0;
  }
  if(app.lo_band < 0.0 && app.hi_band <= 0.0)
  {
   printf("Warning: FFT filter disabled\n");
  }
  else
  {
   app.lo_rem = (unsigned)(fres * app.lo_band /* + 0.5 */);
   app.hi_rem = (unsigned)(fres * app.hi_band /* + 0.5 */);
   printf("-- FFT pass filter [%.3f..%.3f] Hz\n",
	app.lo_band >= 0.0? (double)app.lo_rem / fres : 0.0,
	app.hi_band >  0.0? (double)app.hi_rem / fres : (double)app.hcw.sample_rate / 2.0);
  }
 }
}

/* The FFT based processing function (unsafe version)
*/
void ProcessFFT(void)
{
 fftw_complex *dataL = NULL, *dataR = NULL;
 fftw_plan pL, pR;
 unsigned i;
 unsigned plflags = FFTW_ESTIMATE;

 if(sizeof(char *) <= 4 &&			// limit for 2GB virtual memory
	app.hcw.n_samples > 50000000)		// for 32-bit system only
  error("Input too long for FFT processing @ 32- bits system");

 if(app.isFFTnoSIMD)				// silly / stupid
  plflags |= (FFTW_UNALIGNED | FFTW_NO_SIMD);

// FFTW direct plans
 if(app.isVerbose)
  printf("-- Creating Forward FFTW plans\n");
 dataL = (fftw_complex *)fftw_malloc(((size_t)app.nsFFT) * sizeof(fftw_complex));
 dataR = (fftw_complex *)fftw_malloc(((size_t)app.nsFFT) * sizeof(fftw_complex));
 if(NULL == dataL || NULL == dataR)
  error("Not enoght memory for FFTW buffers");

 chkplan(pL = fftw_plan_dft_r2c_1d(app.nsFFT, (double *)dataL, dataL, plflags));
 chkplan(pR = fftw_plan_dft_r2c_1d(app.nsFFT, (double *)dataR, dataR, plflags));
 printPlanStat("Forward", pL);

// prepare real data
 for(i = 0; i < app.hcw.n_samples; ++i)
 {
  readWavSample(app.fpif, &((double *)dataL)[i], &((double *)dataR)[i]);
  ((double *)dataL)[i] *= app.gain_mul;
  ((double *)dataR)[i] *= app.gain_mul;
 }
 if(app.hcw.n_samples < app.nsFFT)
 {
  ((double *)dataL)[app.hcw.n_samples] = 0.0;
  ((double *)dataR)[app.hcw.n_samples] = 0.0;
 }

// make forward FFT's
 if(app.isVerbose)
  printf("-- Make forward FFT's for the Left Channel\n");
 myFFTWexecute(pL);

 if(app.isVerbose)
  printf("-- Make forward FFT's for the Right Channel\n");
 myFFTWexecute(pR);

 fftw_destroy_plan(pR);
 fftw_destroy_plan(pL);

// processing complex spectrum
 processSpectrum(dataL);
 processSpectrum(dataR);

// FFTW backward plans
 if(app.isVerbose)
  printf("-- Creating Backward FFTW plans\n");

 chkplan(pL = fftw_plan_dft_1d(app.nsFFT, dataL, dataL, FFTW_BACKWARD, plflags));
 chkplan(pR = fftw_plan_dft_1d(app.nsFFT, dataR, dataR, FFTW_BACKWARD, plflags));
 printPlanStat("Backward", pL);

// make backward FFT's
 if(app.isVerbose)
  printf("-- Make backward FFT's for the Left Channel\n");
 myFFTWexecute(pL);

 if(app.isVerbose)
  printf("-- Make backward FFT's for the Right Channel\n");
 myFFTWexecute(pR);

 fftw_destroy_plan(pR);
 fftw_destroy_plan(pL);

// write the data
 for(i = 0; i < app.hcw.n_samples; ++i)
 {
  writeComplex(app.fpof, dataL[i][0], dataL[i][1],
	dataR[i][0], dataR[i][1], &app.hcw,
	&app.l_clips, &app.r_clips, &app.tcrc);
 }

 fftw_free(dataR);
 fftw_free(dataL);
 fftw_cleanup();
}

/* The FFT based processing function (safe version)
*/
void ProcessFFT_Safe(void)
{
 fftw_complex *data = NULL, dl;
 unsigned i;
 double dummy;
 long beginWavData;
 FILE *fpt;
 char *nmt;
 unsigned plflags = FFTW_ESTIMATE;

 if(sizeof(char *) < 8 &&			// limit for 2GB virtual memory
	app.hcw.n_samples > 65000000)		// for 32-bit system only
  error("Input too long for FFT processing");

 // prepare for processing
 if(app.isFFTnoSIMD)				// silly and stupid
  plflags |= (FFTW_UNALIGNED | FFTW_NO_SIMD);

 beginWavData = cftell(app.fpif);

 data = (fftw_complex *)fftw_malloc(((size_t)app.nsFFT) * sizeof(fftw_complex));
 if(NULL == data)
  error("Not enoght memory for FFTW buffer");


 // prepare real data -- left channel
 for(i = 0; i < app.hcw.n_samples; ++i)
 {
  readWavSample(app.fpif, &((double *)data)[i], &dummy);
  ((double *)data)[i] *= app.gain_mul;
 }
 if(app.hcw.n_samples < app.nsFFT)
  ((double *)data)[app.hcw.n_samples] = 0.0;

 // process left channel
 ProcessSingleChFFT(data, plflags, "Left");

 // don't need the second copy of plan statistics
 app.isFFTstat = 0;
 app.isPlanOut = 0;

 // keep the data
 if(app.isVerbose)
  printf("-- Saving temporary data for the Left channel\n");

 nmt = ctempfile();
 fpt = cfopen(nmt, "w+", "temporary data");
 // here we can't use singe fwrite -- some C-runtimes `fwrite` hang up
 // at files >= 4 GB
 for(i = 0; i < app.hcw.n_samples; ++i)
  if(fwrite(&data[i], sizeof(fftw_complex), 1, fpt) != 1)
  {
   fclose(fpt);
   error("Write Error for temporary data");
  }

 // prepare real data -- right channel
 cfseek(app.fpif, beginWavData);
 for(i = 0; i < app.hcw.n_samples; ++i)
 {
  readWavSample(app.fpif, &dummy, &((double *)data)[i]);
 ((double *)data)[i] *= app.gain_mul;
 }
 if(app.hcw.n_samples < app.nsFFT)
  ((double *)data)[app.hcw.n_samples] = 0.0;

 // process right channel
 ProcessSingleChFFT(data, plflags, "Right");

 // write target CWAVE file
 if(app.isVerbose)
  printf("-- Saving target CWAVE data\n");
 rewind(fpt);
 for(i = 0; i < app.hcw.n_samples; ++i)
 {
  if(fread(dl, sizeof(fftw_complex), 1, fpt) != 1)
   error("Read error for temporary data");
  writeComplex(app.fpof, dl[0], dl[1],
	data[i][0], data[i][1], &app.hcw,
	&app.l_clips, &app.r_clips, &app.tcrc);
 }

 fclose(fpt);
 remove(nmt);
 free(nmt);
 fftw_free(data);
}

/* single channel processing
*/
static void ProcessSingleChFFT(fftw_complex *data,
	unsigned plflags, const char *nmCh)
{
 fftw_plan p;

 // FFTW direct plan
 if(app.isVerbose)
  printf("-- Creating Forward FFTW plan for the %s channel\n", nmCh);
 chkplan(p = fftw_plan_dft_r2c_1d(app.nsFFT, (double *)data, data, plflags));
 printPlanStat("Forward", p);

 // make forward FFT's
 if(app.isVerbose)
  printf("-- Make forward FFT's for the %s channel\n", nmCh);
 myFFTWexecute(p);

 fftw_destroy_plan(p);
 fftw_cleanup();

 // processing complex spectrum
 processSpectrum(data);

 // FFTW backward plans
 if(app.isVerbose)
  printf("-- Creating Backward FFTW plan for the %s channel\n", nmCh);

 chkplan(p = fftw_plan_dft_1d(app.nsFFT, data, data, FFTW_BACKWARD, plflags));
 printPlanStat("Backward", p);

 // make backward FFT's
 if(app.isVerbose)
  printf("-- Make backward FFT's for the %s channel\n", nmCh);
 myFFTWexecute(p);

 fftw_destroy_plan(p);
 fftw_cleanup();
}

/* make spectrum processing for the one channel
*/
static void processSpectrum(fftw_complex *data)
{
 unsigned n = (app.nsFFT + 1U) / 2U, mid = app.nsFFT / 2U, i;
 double snorm = 2.0 / (double)app.nsFFT;

 for(i = 1; i < n; ++i)
 {
  data[i][0] *= snorm;
  data[i][1] *= snorm;
  data[i + mid][0] = 0.0;
  data[i + mid][1] = 0.0;
 }
 // adjust DC bin (must be zero ;-)
 data[0][0] /= (double)(app.nsFFT);
 data[0][1] /= (double)(app.nsFFT);
 // rest of negative frequencies
 if((app.nsFFT & 1) == 0)			// mid == n (!!)
  data[mid][0] = data[mid][1] = 0.0;
 // here we have "normal" spectrum to make complex/analitic signal
 // ..filtering..
 if(app.lo_band >= 0.0)
 {
  for(i = 0; i < app.lo_rem; ++i)
  {
   data[i][0] = data[i][1] = 0.0;
  }
 }
 if(app.hi_band > 0.0)
 {
  for(i = app.hi_rem; i < n; ++i)
  {
   data[i][0] = data[i][1] = 0.0;
  }
 }
}

/* check a plan for existence
*/
static void chkplan(const fftw_plan plan)
{
 if(NULL == plan)
  error("Not enough memory for FFTW plan creating");
}

/* execute DFT with exception check (for modified FFTW only!)
*/
void myFFTWexecute(const fftw_plan plan)
{
 volatile int isFail = 0;

 TRY
 {
  fftw_execute(plan);
 }
 CATCH_ALL
 {
  isFail = 1;
 }
 if(isFail)
  error("Sorry, out of memory during execute FFT; try different N points");
}

/* print plan stuff
*/
static void printPlanStat(const char *plName, const fftw_plan plan)
{
 double n_adds, n_muls, n_fmas;

 if((!app.isFFTstat) && (!app.isPlanOut))
  return;

 fftw_flops(plan, &n_adds, &n_muls, &n_fmas);
 if(app.isFFTstat)
 {
  printf("-- %s FFT[%u]: +:%.0f, *:%.0f, *+:%.0f\n",
	plName, app.nsFFT, n_adds, n_muls, n_fmas);
 }
 if(app.isPlanOut)
 {
  if(!app.isFFTstat)
  {
   printf("-- %s Plan[%u]: +:%.0f, *:%.0f, *+:%.0f\n",
	plName, app.nsFFT, n_adds, n_muls, n_fmas);
  }
  fftw_print_plan(plan);
  printf("\n");
 }
}

/* the end...
*/

