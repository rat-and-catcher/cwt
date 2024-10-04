/*
 * cw.c -- the main implementation file windows-only
 * 2-ch wav to atalitic (complex) signal transformation;
 * (FIR based Hilbert transform + direct FFT-based transform)
 * This program can be distributed under GNU GPL
 * Copyright (C) 2010-2013 Rat and Catcher Tech.
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

#include "cw.h"

/*
 * global variable(s)
 */
APPLICATION_CW app;

/*
 * the forward declarations
 */
/* set initial state of application
*/
static void VarInit(void);
/* command line parser
*/
static void parseCommandLine(int argc, char **argv);
/* show version info
*/
static void verinfo(void);
/* show help text
*/
static void help(void);
/* The main processing function
*/
static void TheProcess(void);
/* check CRC of CWAVE
*/
static void CheckCrcProcess(void);
/* the single thread processing routine
*/
static void SingleTreadProcess(void);
/* the multi tread processing routine
*/
static void MultiThreadProcess(void);
/* allocate and init all the need to LRCH
*/
static void PrepareChBufs(LRCH *ch, const char *comment);
/* free LRCH buffers
*/
static void FreeChBufs(LRCH *ch);
/* make Hilbert transform -- reverse memory scan order (as previous versions)
*/
static void rMakeHilbert(LRCH *ch);
/* make Hilbert transform -- direct memory scan order
*/
static void dMakeHilbert(LRCH *ch);
/* make Hilbert transform for buffer
*/
static unsigned WINAPI hilb_thr(TCH *ch);
/* allocate TCH and init it
*/
static TCH *PrepareThrChBufs(const char *comment);
/* free TCH
*/
static void FreeThrChBufs(TCH *ch);

/*
 * The code
 */
/* the main function
*/
int main(int argc, char **argv)
{
 printf("cw -- real(wav) to analitic(cwave) file converter, version %s\n", VERSION);

 VarInit();

 parseCommandLine(argc, argv);

 TheProcess();

 return 0;
}

/* set initial state of the application
*/
static void VarInit(void)
{
 app.nameif = app.nameof = NULL;

 app.k_M = DEF_M;
 app.k_beta = DEF_BETA;
 app.hpr = NULL;

 app.gain_mul = DEF_GAIN;
 app.isPrintH = 0;
 app.isTestCRC = 0;
 app.isVerbose = 0;
 app.notCompHeadTail = 0;
 app.isDirectScan = 0;
 app.c_format = HCW_FMT_BAD_FMT;
 app.nThr = GetCpuNumber();

 app.isFFT = 0;
 app.isFFTeven = 0;
 app.isFFTsafe = 0;
 app.isFFTstat = 0;
 app.isPlanOut = 0;
 app.isFFTnoSIMD = 0;
 app.nsFFT = 0;
 app.lo_band = app.hi_band = -1.0;
 app.lo_rem = app.hi_rem = 0;

 app.progress = -1.0;
 time(&app.tstart);
 app.fpif = app.fpof = NULL;
 crc32init(&app.tcrc);
 app.l_clips = app.r_clips = 0;
 app.MakeHilbert = NULL;
}

/* command line parser
*/
static void parseCommandLine(int argc, char **argv)
{
 static const char extWav[] = "wav", extCwave[] = "cwave";
 char *s, dmy;
 int k;

 while(--argc)
 {
  s = *++argv;
  if('-' == *s)         // key
  {
   switch(s[1])
   {
    case 'h':           // help message
    case '?':
     help();
     exit(0);
     break;

    case 'q':           // additional version info
     verinfo();
     exit(0);
     break;

    case 'f':           // conversion via FFT(W)
     app.isFFT = 1;
     for(k = 2; s[k]; ++k)
     {
      switch(s[k])
      {
       case 'e':        // make "even" FFT
        app.isFFTeven = 1;
        break;
       case 's':        // make "safe" FFT
        app.isFFTsafe = 1;
        break;
       case 'i':        // don't use SIMD (SSE2)
        app.isFFTnoSIMD = 1;
        break;
       default:         // BAD
        error("Bad FFT key modificator: %s", s);
        break;
      }
     }
     break;

    case 'r':           // remove frequencies (FFT only)
     switch(s[2])
     {
      case 'l':         // lower pass frequency
       if(--argc <= 0 || 1 != sscanf(*++argv, "%lf %c", &app.lo_band, &dmy))
       {
        error("Bad or absent lower frequency parameter");
       }
       if(app.lo_band < 0.0)
       {
        error("Lower frequency parameter must be non-negative");
       }
       break;
      case 'h':         // higher pass frequency
       if(--argc <= 0 || 1 != sscanf(*++argv, "%lf %c", &app.hi_band, &dmy))
       {
        error("Bad or absent higher frequency parameter");
       }
       if(app.hi_band <= 0.0)
       {
        error("Higher frequency parameter must be positive");
       }
       break;
      case '2':         // -r21: 21..21000 Hz
       if('1' == s[3])
       {
        app.lo_band = 21.0;
        app.hi_band = 21000.0;
       }
       else
       {
        error("Did you mean '-r21' (%s)?", s);
       }
       break;
      default:          // BAD
       error("Bad FFT filter modificator: %s", s);
       break;
     }
     break;

    case 'M':           // filter order
     if(--argc <= 0 || 1 != sscanf(*++argv, "%d %c", &app.k_M, &dmy))
     {
      error("Bad or absent FIR filter order value");
     }
     if(app.k_M <= 0)
     {
      error("FIR filter order must be positive");
     }
     break;

    case 'b':           // "beta" Kaiser's parameter
     if(--argc <= 0 || 1 != sscanf(*++argv, "%lf %c", &app.k_beta, &dmy))
     {
      error("Bad or absent Kaiser 'beta' parameter");
     }
     if(app.k_beta < 0.0)
     {
      error("Kaiser 'beta' parameter must be non-negative");
     }
     break;

    case 'g':           // gain multiplier
     if(--argc <= 0 || 1 != sscanf(*++argv, "%lf %c", &app.gain_mul, &dmy))
     {
      error("Bad or absent gain multiplier");
     }
     if(app.gain_mul <= 0.0)
     {
      error("Gain multiplier must be positive");
     }
     break;

    case 'p':           // print Hilbert's FIR filter coefficients
     app.isPrintH = 1;
     break;

    case 't':           // CRC check
     app.isTestCRC = 1;
     break;

    case 'v':           // verbose mode
     app.isVerbose = 1;
     for(k = 2; s[k]; ++k)
     {
      switch(s[k])
      {
       case 'p':        // print FFTW plans
        app.isPlanOut = 1;
        break;
       case 's':        // print FFTW statistics
        app.isFFTstat = 1;
        break;
       default:         // BAD
        error("Bad verbose key modificator: %s", s);
        break;
      }
     }
     break;

    case 'c':           // cancel time compensation for the group delay (FIR)
     app.notCompHeadTail = 1;
     break;

    case 'd':           // direct memory scan (FIR)
     app.isDirectScan = 1;
     break;

    case 'i':           // CWAVE format variation
     switch(s[2])
     {
      case 'd':         // double Re(L), Im(L), Re(R), Im(R), [-32768.0..32767.0]
       app.c_format = HCW_FMT_PCM_DBL64;
       break;
      case 's':         // signed short Re(L), Im(L), Re(R), Im(R), [-32768..32767]
       app.c_format = HCW_FMT_PCM_INT16;
       break;
      case 'm':         // signed short Re(.), float Im(.), [-32768..32767]/[-32768.0..32767.0]
       app.c_format = HCW_FMT_PCM_INT16_FLT32;
       break;
      case 'f':         // float Re(L), Im(L), Re(R), Im(R), [-32768.0..32767.0]
       app.c_format = HCW_FMT_PCM_FLT32;
       break;
      default:
       error("Illegal CWAVE format specification: '%s'", s);
      break;
     }
     break;

    case '1':           // force single-thread version
     app.nThr = 1;
     break;

    case '2':           // force multi-thread version (2 threads for FFTW)
     app.nThr = 2;
     break;

    default:
     error("Illegal key '%s'", s);
     break;
   }
  }
  else                  // file name
  {
   if(NULL == app.nameif)
   {
    app.nameif = s;
   }
   else
   {
    if(NULL == app.nameof)
    {
     app.nameof = s;
    }
    else
     error("More than two file names in a command line");
   }
  }
 }

 // checking parameters
 if(!app.isPrintH && !app.isTestCRC &&
        (NULL == app.nameif || NULL == app.nameof))
  error("You MUST specify input.WAV and output.CWAVE for processing");
 if(app.isPrintH && app.isTestCRC)
  error("-p and -t switches are incompatible");
 if(app.isPrintH && (NULL != app.nameif || NULL != app.nameof))
  error("-p switch incompatible with file processing");
 if(app.isTestCRC && (NULL == app.nameif || NULL != app.nameof))
  error("-t switch need only one input.CWAVE file");
 if(app.isFFT && (app.isPrintH || app.isTestCRC))
  error("-f incompatible with -p or -t switches");

 if(!app.isPrintH && !app.isTestCRC)
 {
  if(!checkFileExt(app.nameif, extWav) || !checkFileExt(app.nameof, extCwave))
   error("Input file must have .WAV, and output file - .CWAVE extensions");
 }

 if(app.isTestCRC)
 {
  if(!checkFileExt(app.nameif, extCwave))
   error("For CRC testing input file must have .CWAVE extension");
 }
 else
 {
  if(!app.isFFT)
  {                                     // FIR
   if(app.isFFTstat || app.isPlanOut)
   {
    printf("Warning: there is no FFTW statistics for FIR conversion\n");
    app.isFFTstat = 0;
    app.isPlanOut = 0;
   }
   if(app.lo_band >= 0.0 || app.hi_band > 0.0)
   {
    printf("Warning: Filtering possible for the FFT mode only\n");
    app.lo_band = app.hi_band = -1.0;
   }
   if(app.k_M & 1)
   {
    printf("Warning: Order of Hilbert FIR filter must be even\n");
    ++app.k_M;
   }
   if(((app.k_M / 2) & 1) == 0)
   {
    printf("Warning: Order / 2 of FIR Hilbert filter must be odd\n");
    app.k_M += 2;
   }
   if(app.k_M < 2)                      // unreachable
   {
    printf("Warning: Order of FIR Hilbert filter too low, set to min (2)\n");
    app.k_M = 2;
   }
   app.n_delay = app.k_M / 2;           // always >= 1 and ODD !!
   printf("-- Real to Analitic conversion: FIR filter order %d; Kaiser's beta is %g\n",
        app.k_M, app.k_beta);
  }
  else                                  // FFT
  {
   app.k_M = -1;
   app.k_beta = 0.0;
   printf("-- Real to Analitic conversion: using direct FFT algorithm\n");
   if(app.lo_band >= 0.0 && app.hi_band >= 0.0)
   {
    if(app.lo_band >= app.hi_band)
     error("Low frequency (%f Hz) must be lower than the high (%f Hz)",
        app.lo_band, app.hi_band);
   }
  }
  if(!app.isPrintH)
  {
   app.MakeHilbert = app.isDirectScan? dMakeHilbert : rMakeHilbert; 

   printf("-- Gain multiplier %g (%g%% == %g dB)\n",
        app.gain_mul, ((double)((int)(app.gain_mul * 1000.0 + 0.5))) / 10.0,
        20.0 * log10(app.gain_mul));
   if(HCW_FMT_BAD_FMT == app.c_format)
   {
    app.c_format = app.isFFT? HCW_FMT_PCM_FLT32 : HCW_FMT_PCM_INT16_FLT32;
   }
   PrintCwaveFormat(app.c_format);
  }
 }
}

/* show version info
*/
static void verinfo(void)
{
 printf("(Platform %s, %s, %s libs; FFTW: %s/%s/%s)\n",
        PLATFORM, TYPE_BUILD, TYPE_LIB,
        fftw_version, fftw_cc,
        fftw_codelet_optim[0]? fftw_codelet_optim : "std");
}

/* show help text
*/
static void help(void)
{
 verinfo();
 printf("Usage: cw [keys] input.wav output.cwave\n"
        "or cw -t [keys] input.cwave\n"
        "or cw -p [keys]\n"
        "or cw -f[e][s][i] [keys] input.wav output.cwave\n"
        "or cw -h\n"
        "or cw -q\n"
        "where *.wav - 16-bps stereo Windows PCM wav file\n"
        "*.cwave - target special wave format in complex form.\n"
        "The input and output files MUST have correct extensions.\n"
        "Possible keys:\n"
        "-h, -? - print this text\n"
        "-q - print advanced version info\n"
        "-M value - order of Hilbert FIR filter, (integer, even; M/2 must be odd)\n"
        "-b value - Kaiser's window parameter for Hilbert FIR filter (float)\n"
        "-g value - set gain multiplier for input samples\n"
        "-p - print only Hilbert FIR filter coefficients (w/o any files)\n"
        "-t - check input.CWAVE (w/o any output) for integrity\n"
        "-f[e][s][i] - make conversion via FFT, not Hilbert FIR filter\n"
        "-fe - same as -f, but make total number of FFT points strictly even\n"
        "      (default - strictly odd, regardless real number of samples)\n"
        "-fs - same as -f, but slowly and safely (recommended for big files)\n"
        "-fi - same as -f, but w/o using SIMD (SSE2) instructions (silly:);\n"
        "-rX freq - (FFT only) remove low/high frequences from the spectrum:\n"
        "-rl freq - from 0 to freq, Hz; -rh freq - form freq, Hz to maximum,\n"
        "-r21 equvalent for '-rl 21 -rh 21000' (standard for Audio CD tracks)\n"
        "-v[s][p] - verbose mode; [s] [p] modificators works for FFT only:\n"
        "-vs - print FFTW statistics; -vp - print FFTW plans\n"
        "-c - turn off compensate time shift for FIR filters group delay\n"
        "-d - strictly direct memory scan during FIR-based Hilbert transform\n"
        "   (probably give better performance for some (overloaded) systems)\n"
        "-iX - set sample format for the output CWAVE file\n"
        "(-id == double; -is == short(16-bit);\n"
        " -im == short+float; -if == float)\n"
        "-1 - force single-thread version\n"
        "-2 - force multi-thread version (both FIR/FFT)\n"
        "NOTE: most of options for -t and -p will be ignored\n"
        "NOTE: only -g, -v[s][p], -1, -2, -rX and -iX valid with -f option\n"
        "NOTE: if FFTW and cw using different C-runtime library,\n"
        "  -vp option may be very *DANGEROUS*, especially if you wish\n"
        "  redirect stdout (from wrong output to program crash).\n"
        "  So, please use -vp only for entertainment ;)\n"
        );

 printf("** The default values:\n"
        "-M %d; -b %f -g %f\n",
        DEF_M, DEF_BETA, DEF_GAIN);
 printf("** Default CWAVE type -- '-im' for FIR; '-if' for -f (FFT)\n");

 printf("\nCopyright (C) 2010-2016 Rat and Catcher Tech.\n"
        "\n"
        "This program is free software: you can redistribute it and/or modify\n"
        "it under the terms of the GNU General Public License as published by\n"
        "the Free Software Foundation, either version 3 of the License, or\n"
        "(at your option) any later version.\n"
        "\n"
        "This program is distributed in the hope that it will be useful,\n"
        "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
        "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the\n"
        "GNU General Public License for more details.\n"
        "\n"
        "Some code fragments based by the ALGLIB project (see comments in\n"
        "besseli0.* sources)\n"
        "This program also used FFTW library written by Matteo Frigo and\n"
        "Steven G. Johnson. Visit www.fftw.org for details.\n"
        );
}

/* The main processing function
*/
// inline helper(s)
static void INLINE ShowProgress(unsigned cur, unsigned tot)
{
 if(app.isVerbose)
 {
  double t = ((double)cur) * 100.0 / tot;

  if(app.progress < t - 1.0)
  {
   if(ISATTY(FILENO(stdout)))
   {
    app.progress = t;
    printf("-- Done %3d%%\r", (int)app.progress);
    fflush(stdout);
   }
  }
 }
}

// the function
static void TheProcess(void)
{
 int i;

 // very special case
 if(app.isTestCRC)
 {
  CheckCrcProcess();
  return;                               // and here the story end...
 }

 if(!app.isFFT)
 {
  // Hilbert coefficients
  app.hpr = (double *)cmalloc(sizeof(double) * (app.n_delay + 1), "Pulse response buffer");
  for(i = 0; i <= app.n_delay; ++i)
   app.hpr[i] = FilterCoeff(i * 2, app.k_M, app.k_beta);        // all odd values == 0

  if(app.isPrintH)
  {
   printf("-- Hilbert FIR filter coefficients (all odd values == 0):\n");
   for(i = 0; i <= app.n_delay; ++i)
    printf("h[%4d] = %23.14A (%24.16G)\n", i * 2, app.hpr[i], app.hpr[i]);
   if(app.hpr)
   {
    free(app.hpr);
    app.hpr = NULL;
   }
   return;                              // ..and here the end of story..
  }
 }

 // prepare for processing
 app.fpif = cfopen(app.nameif, "rb", "input WAV-data");
 app.fpof = cfopen(app.nameof, "wb", "output complex data");
 readWavHeader(app.fpif, &app.hcw.sample_rate, &app.hcw.n_samples);
 if(!app.isFFT)
  if(app.hcw.n_samples < (unsigned)(app.k_M * 2))
   error("Input too short");
 memcpy(&(app.hcw.magic[0]), HCW_MAGIC, sizeof(app.hcw.magic));
 app.hcw.hsize = sizeof(app.hcw);
 app.hcw.version = HCW_VERSION_BAD;
 app.hcw.format = app.c_format;
 app.hcw.n_channels = 2;
 app.hcw.k_M = app.k_M;
 app.hcw.n_CRC32 = 0;
 app.hcw.k_beta = app.k_beta;
 writeCwaveHeader(app.fpof, &app.hcw);

 // ..processing..
 if(app.isFFT)
 {
  if(!fftw_init_threads())
   error("Can't initialized multi-tread FFTW");
  fftw_plan_with_nthreads(app.nThr);

  app.nsFFT = app.isFFTeven?                            // STRICTLY even or STRICTLY odd
        (app.hcw.n_samples + 1U) & (~01U) : app.hcw.n_samples | 01U;

  CalcFFTpass();                // FFT filter calculations

  if(app.isFFTsafe)
   ProcessFFT_Safe();
  else
   ProcessFFT();

  fftw_cleanup_threads();
 }
 else
 {
  if(app.nThr < 2)
   SingleTreadProcess();
  else
   MultiThreadProcess();
 }

 // finish
 app.hcw.n_CRC32 = crc32final(&app.tcrc);
 app.hcw.version = HCW_VERSION_CUR;
 rewind(app.fpof);
 writeCwaveHeader(app.fpof, &app.hcw);
 if(app.isVerbose)
  CalcTime(app.tstart);
 if(HCW_FMT_PCM_INT16 == app.hcw.format ||
        (app.isFFT && HCW_FMT_PCM_INT16_FLT32 == app.hcw.format))
  printf("-- Clipping statistics: Left %ld, Right %ld clips\n",
        app.l_clips, app.r_clips);
 printf("-- Conversion OK, CRC is 0x%08X\n", app.hcw.n_CRC32);

 // cleanup
 if(app.fpof)
 {
  fclose(app.fpof);
  app.fpof = NULL;
 }
 if(app.fpif)
 {
  fclose(app.fpif);
  app.fpif = NULL;
 }
 if(app.hpr)
 {
  free(app.hpr);
  app.hpr = NULL;
 }
}

/* check CRC of CWAVE
*/
static void CheckCrcProcess(void)
{
 unsigned crc;
 unsigned i;

 app.fpif = cfopen(app.nameif, "rb", "input complex data");
 readCwaveHeader(app.fpif, &app.hcw);
 printf("-- Conv. parameters: Hilbert order %d; Kaiser's beta is %g\n",
        app.hcw.k_M, app.hcw.k_beta);
 PrintCwaveFormat(app.hcw.format);
 printf("-- Header version V%u; %u channels; %u samples; sample rate %u Hz\n",
        app.hcw.version, app.hcw.n_channels, app.hcw.n_samples, app.hcw.sample_rate);

 for(i = 0; i < app.hcw.n_samples; ++i)
 {
  readComplex(app.fpif, app.hcw.format, &app.tcrc);
  ShowProgress(i, app.hcw.n_samples);
 }

 crc = crc32final(&app.tcrc);

 if(app.isVerbose)
  CalcTime(app.tstart);

 if(app.hcw.version < HCW_VERSION_V2)
 {
  printf("** This file does not contain CRC. Calculated CRC is 0x%08X\n", crc);
 }
 else
 {
  if(crc == app.hcw.n_CRC32)
   printf("-- CRC32 (0x%08X) OK!\n", crc);
  else
   printf("** Real CRC (0x%08X) MISMATCH with header (0x%08X)\n"
        "** DATA SEEMS CORRUPTED!!!\n",
        crc, app.hcw.n_CRC32);
 }
 if(app.fpif)
 {
  fclose(app.fpif);
  app.fpif = NULL;
 }
}

/* the single thread processing routine
*/
static void SingleTreadProcess(void)
{
 LRCH lch;                              // left channel data
 LRCH rch;                              // right channel data 
 unsigned cnt, n;

 PrepareChBufs(&lch, "Data for the left channel");
 PrepareChBufs(&rch, "Data for the right channel");

 cnt = app.notCompHeadTail? app.hcw.n_samples : app.hcw.n_samples + app.n_delay;
 for(n = 0; n < cnt; ++n)
 {
  if(n < app.hcw.n_samples)
  {
   readWavSample(app.fpif, &(lch.s_in), &(rch.s_in));
   lch.s_in *= app.gain_mul;
   rch.s_in *= app.gain_mul;
  }
  else
  {
   lch.s_in = rch.s_in = 0.0;
  }
  (*app.MakeHilbert)(&lch);
  (*app.MakeHilbert)(&rch);

  // here we have a complex signal
  if(n >= (unsigned)app.n_delay || app.notCompHeadTail)
  {
   writeComplex(app.fpof, lch.s_real, lch.s_image,
        rch.s_real, rch.s_image, &app.hcw,
        &app.l_clips, &app.r_clips, &app.tcrc);
  }
  ShowProgress(n, cnt);
 }

 // cleanup
 FreeChBufs(&lch);
 FreeChBufs(&rch);
}

/* the multi tread processing routine
*/
static void MultiThreadProcess(void)
{
 HANDLE waitOut[2];
 HANDLE hthr_l, hthr_r;
 TCH *lch, *rch;
 unsigned cnt, full, tail, nr, nw, nb, i;

 lch = PrepareThrChBufs("left channel tread context");
 rch = PrepareThrChBufs("right channel thread context");

 waitOut[0] = lch -> waitOut;   // copy for WaitForMultiplyObjects()
 waitOut[1] = rch -> waitOut;

 hthr_l = (HANDLE)_beginthreadex(NULL, 0, hilb_thr, lch, 0, NULL);
 hthr_r = (HANDLE)_beginthreadex(NULL, 0, hilb_thr, rch, 0, NULL);
 if(NULL == hthr_l || NULL == hthr_r)
  error("Can't create working thread(s)");
 if(NORMAL_PRIORITY_CLASS == GetPriorityClass(GetCurrentProcess()))
 {
  SetThreadPriority(hthr_l, THREAD_PRIORITY_HIGHEST);
  SetThreadPriority(hthr_r, THREAD_PRIORITY_HIGHEST);
 }

 cnt = app.notCompHeadTail? app.hcw.n_samples : app.hcw.n_samples + app.n_delay;
 tail = cnt % MAX_SAMPLES_PER_TIME;
 full = cnt / MAX_SAMPLES_PER_TIME;

 // 1st stage -- full buffers
 lch -> nsdata = MAX_SAMPLES_PER_TIME;
 rch -> nsdata = MAX_SAMPLES_PER_TIME;
 for(nb = nr = nw = 0; nb < full; ++nb)
 {
  // fill buffers
  for(i = 0; i < MAX_SAMPLES_PER_TIME; ++i, ++nr)
  {
   if(nr < app.hcw.n_samples)
   {
    readWavSample(app.fpif, &(lch -> indata[i]), &(rch -> indata[i]));
   }
   else
   {
    lch -> indata[i] = rch -> indata[i] = 0.0;
   }
  }

  // run process
  SetEvent(lch -> waitIn);
  SetEvent(rch -> waitIn);
  WaitForMultipleObjects(2, waitOut, TRUE, INFINITE);

  // store results
  for(i = 0; i < MAX_SAMPLES_PER_TIME; ++i, ++nw)
  {
   if(nw >= (unsigned)app.n_delay || app.notCompHeadTail)
   {
    writeComplex(app.fpof, lch -> re_outdata[i], lch -> im_outdata[i],
        rch -> re_outdata[i], rch -> im_outdata[i], &app.hcw,
        &app.l_clips, &app.r_clips, &app.tcrc);
   }
  }
  ShowProgress(nw, cnt);
 }

 // 2nd stage -- the tail
 lch -> nsdata = tail;
 lch -> isEnd = 1;
 rch -> nsdata = tail;
 rch -> isEnd = 1;
 // fill buffers
 for(i = 0; i < tail; ++i, ++nr)
 {
  if(nr < app.hcw.n_samples)
  {
   readWavSample(app.fpif, &(lch -> indata[i]), &(rch -> indata[i]));
  }
  else
  {
   lch -> indata[i] = rch -> indata[i] = 0.0;
  }
 }

 // run process
 SetEvent(lch -> waitIn);
 SetEvent(rch -> waitIn);
 WaitForMultipleObjects(2, waitOut, TRUE, INFINITE);

 // store results
 for(i = 0; i < tail; ++i, ++nw)
 {
  if(nw >= (unsigned)app.n_delay || app.notCompHeadTail)
  {
   writeComplex(app.fpof, lch -> re_outdata[i], lch -> im_outdata[i],
        rch -> re_outdata[i], rch -> im_outdata[i], &app.hcw,
        &app.l_clips, &app.r_clips, &app.tcrc);
  }
  ShowProgress(nw, cnt);
 }

 // cleanup
 WaitForSingleObject(hthr_l, INFINITE);
 WaitForSingleObject(hthr_r, INFINITE);
 CloseHandle(hthr_l);
 CloseHandle(hthr_r);

 FreeThrChBufs(lch);
 FreeThrChBufs(rch);
}

/*
 * operations with LRCH
 */
/* allocate and init all the need to LRCH
*/
static void PrepareChBufs(LRCH *ch, const char *comment)
{
 int i;

 // just in case
 ch -> s_buf = NULL;

 // buffers
 ch -> s_buf = (double *)cmalloc(sizeof(double) * app.k_M, comment);
 for(i = 0; i < app.k_M; ++i)
 {
  ch -> s_buf[i] = 0.0;
 }

 // indexes
 if(app.isDirectScan)
 {
  ch -> ixCur = app.k_M - 1;        // X[i] for input/output data, decrease for each sample
  ch -> ixOutRe = app.n_delay - 1;  // X[i] for output of the Real part of data, decrease for each sample
 }
 else
 {
  ch -> ixCur = 0;                  // X[i] for input/output data, increase for each sample
  ch -> ixOutRe = app.n_delay;      // X[i] for output of the Real part of data, increase for each sample
 }
}

/* free LRCH buffers
*/
static void FreeChBufs(LRCH *ch)
{
 if(ch -> s_buf)
 {
  free(ch -> s_buf);
  ch -> s_buf = NULL;
 }
}

/* make Hilbert transform -- reverse memory scan order (as previous versions)
*/
static void rMakeHilbert(LRCH *ch)
{
 double hilb;
 int i, ix;

 hilb = ch -> s_in * app.hpr[0];
 for(i = 1, ix = ch -> ixCur; i <= app.n_delay; ++i)
 {
  ix -= 2;                                      // take the PREVIOUS value (step 2)
  if(ix < 0)
   ix += app.k_M;
  hilb += ch -> s_buf[ix] * app.hpr[i];
 }

 ch -> s_real  = ch -> s_buf[ch -> ixOutRe];
 ch -> s_image = hilb;

 ch -> s_buf[ch -> ixCur] = ch -> s_in;

 // indexes (we can use mod(%), because it make rarely)
 ch -> ixCur = (ch -> ixCur + 1) % app.k_M;
 ch -> ixOutRe = (ch -> ixOutRe + 1) % app.k_M;
}

/* make Hilbert transform -- direct memory scan order
*/
static void dMakeHilbert(LRCH *ch)
{
 double hilb;
 int i, ix;

 hilb = ch -> s_in * app.hpr[0];
 for(i = 1, ix = ch -> ixCur; i <= app.n_delay; ++i)
 {
  ix += 2;                              // take the PREVIOUS value (step 2)
  if(ix >= app.k_M)                     // we think, that the mod (%) operation is too slowly(!)
   ix -= app.k_M;                       // so ix = (ix + 2) % k_M; looks BAD
  hilb += ch -> s_buf[ix] * app.hpr[i];
 }

 ch -> s_real  = ch -> s_buf[ch -> ixOutRe];
 ch -> s_image = hilb;

 ch -> s_buf[ch -> ixCur] = ch -> s_in;

 // indexes
 --(ch -> ixCur);
 if(ch -> ixCur < 0)
  ch -> ixCur += app.k_M;

 --(ch -> ixOutRe);
 if(ch -> ixOutRe < 0)
  ch -> ixOutRe += app.k_M;
}

/*
 * Working thread function(s)
 */
/* make Hilbert transform for buffer
*/
static unsigned WINAPI hilb_thr(TCH *ch)
{
 int i, n, isEnd;

 do 
 {
  WaitForSingleObject(ch -> waitIn, INFINITE);
  n = ch -> nsdata;
  isEnd = ch -> isEnd;

  for(i = 0; i < n; ++i)
  {
   ch -> ch.s_in = ch -> indata[i] * app.gain_mul;
   (*app.MakeHilbert)(&(ch -> ch));
   ch -> re_outdata[i] = ch -> ch.s_real;
   ch -> im_outdata[i] = ch -> ch.s_image;
  }

  SetEvent(ch -> waitOut);      // portion complete!
// NOTE:: after SetEvent() ch->isEnd *not valid*, so we use its copy,
// which was received after WaitFor...
 }
 while(!isEnd);

 _endthreadex(0);
 return 0;
}

/*
 * operations with TCH
 */
/* allocate TCH and init it
*/
static TCH *PrepareThrChBufs(const char *comment)
{
 TCH *res = (TCH *)cmalloc(sizeof(TCH), comment);

 PrepareChBufs(&(res -> ch), comment);
 res -> nsdata = 0;
 res -> isEnd = 0;
 res -> waitIn = ccreate_event();
 res -> waitOut = ccreate_event();
 return res;
}

/* free TCH
*/
static void FreeThrChBufs(TCH *ch)
{
 if(ch)
 {
  CloseHandle(ch -> waitIn);
  CloseHandle(ch -> waitOut);
  FreeChBufs(&(ch -> ch));
  free(ch);
 }
}

/* the end...
*/
