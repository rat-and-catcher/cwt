/*
 * cwt.c -- the main implementation file posix-gcc-only tiny
 * 2-ch wav to atalitic (complex) signal transformation;
 * (direct FFT-based transform only)
 * This program can be distributed under GNU GPL
 *
 * Copyright (C) 2010-2013-2014 Rat and Catcher Tech.
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
/* show help text
*/
static void help(void);
/* The main processing function
*/
static void TheProcess(void);
/* check CRC of CWAVE
*/
static void CheckCrcProcess(void);

/*
 * The code
 */
/* the main function
*/
int main(int argc, char **argv)
{
 printf("cwt -- tiny real(wav) to analitic(cwave) file converter, version %s\n", VERSION);

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

 app.gain_mul = DEF_GAIN;
 app.isTestCRC = 0;
 app.isVerbose = 0;
 app.c_format = HCW_FMT_BAD_FMT;
 app.nCPU = 1;					// def. == single thread

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
  if('-' == *s)		// key
  {
   switch(s[1])
   {
    case 'h':		// help message
    case '?':
     help();
     exit(0);
     break;

    case 'f':		// FFT(W) options
     for(k = 2; s[k]; ++k)
     {
      switch(s[k])
      {
       case 'e':	// make "even" FFT
        app.isFFTeven = 1;
        break;
       case 's':	// make "safe" FFT
        app.isFFTsafe = 1;
        break;
       case 'i':	// don't use SIMD (SSE2)
        app.isFFTnoSIMD = 1;
        break;
       default:		// BAD
        error("Bad FFT key modifier: %s", s);
        break;
      }
     }
     break;

    case 'r':		// remove frequencies
     switch(s[2])
     {
      case 'l':		// lower pass frequency
       if(--argc <= 0 || 1 != sscanf(*++argv, "%lf %c", &app.lo_band, &dmy))
       {
        error("Bad or absent lower frequency parameter");
       }
       if(app.lo_band < 0.0)
       {
        error("Lower frequency parameter must be non-negative");
       }
       break;
      case 'h':		// higher pass frequency
       if(--argc <= 0 || 1 != sscanf(*++argv, "%lf %c", &app.hi_band, &dmy))
       {
        error("Bad or absent higher frequency parameter");
       }
       if(app.hi_band <= 0.0)
       {
        error("Higher frequency parameter must be positive");
       }
       break;
      case 's':		// "standard" 21 Hz .. 21 KHz band
      case 'c':		// same as "CD quality"
	app.lo_band = 21.0;
	app.hi_band = 21000.0;
	break;
      case '2':		// -r21 -- "standard" 21 Hz .. 21 KHz band
	if(s[3] != '1')
	{
	 printf("Warning: bad spec. case FFT filter modifier: %s; assume -r21\n", s);
	}
	app.lo_band = 21.0;
	app.hi_band = 21000.0;
	break;
      default:		// BAD
       error("Bad FFT filter modifier: %s", s);
       break;
     }
     break;

    case 'g':		// gain multiplier
     if(--argc <= 0 || 1 != sscanf(*++argv, "%lf %c", &app.gain_mul, &dmy))
     {
      error("Bad or absent gain multiplier");
     }
     if(app.gain_mul <= 0.0)
     {
      error("Gain multiplier must be positive");
     }
     break;

    case 't':		// CRC check
     app.isTestCRC = 1;
     break;

    case 'v':		// verbose mode
     app.isVerbose = 1;
     for(k = 2; s[k]; ++k)
     {
      switch(s[k])
      {
       case 'p':	// print FFTW plans
        app.isPlanOut = 1;
        break;
       case 's':	// print FFTW statistics
        app.isFFTstat = 1;
        break;
       default:		// BAD
        error("Bad verbose key modifier: %s", s);
        break;
      }
     }
     break;

    case 'i':		// CWAVE format variation
     switch(s[2])
     {
      case 'd':		// double Re(L), Im(L), Re(R), Im(R), [-32768.0..32767.0]
       app.c_format = HCW_FMT_PCM_DBL64;
       break;
      case 's':		// signed short Re(L), Im(L), Re(R), Im(R), [-32768..32767]
       app.c_format = HCW_FMT_PCM_INT16;
       break;
      case 'm':		// signed short Re(.), float Im(.), [-32768..32767]/[-32768.0..32767.0]
       app.c_format = HCW_FMT_PCM_INT16_FLT32;
       break;
      case 'f':		// float Re(L), Im(L), Re(R), Im(R), [-32768.0..32767.0]
       app.c_format = HCW_FMT_PCM_FLT32;
       break;
      default:
       error("Illegal CWAVE format specification: '%s'", s);
      break;
     }
     break;

    case 'j':		// set number of loaded CPU (maximum number of threads)
     if(--argc <= 0 || 1 != sscanf(*++argv, "%d %c", &app.nCPU, &dmy))
     {
      error("Bad or absent number of CPU's / threads");
     }
     if(app.nCPU < 1)
     {
      error("Bad number of CPU's / threads: %d; must be > 0", app.nCPU);
     }
     break;

    default:
     error("Illegal key '%s'", s);
     break;
   }
  }
  else			// file name
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
 if(!app.isTestCRC &&
	(NULL == app.nameif || NULL == app.nameof))
  error("You MUST specify input.WAV and output.CWAVE for processing");
 if(app.isTestCRC && (NULL == app.nameif || NULL != app.nameof))
  error("-t switch need only one input.CWAVE file");

 if(!app.isTestCRC)
 {
  if(!checkFileExt(app.nameif, extWav) || !checkFileExt(app.nameof, extCwave))
   error("Input file must have .WAV, and output file - .CWAVE extensions");
 }
 else
 {
  if(!checkFileExt(app.nameif, extCwave))
   error("For CRC testing input file must have .CWAVE extension");
 }

 if(app.lo_band >= 0.0 && app.hi_band >= 0.0)
 {
  if(app.lo_band >= app.hi_band)
   error("Low frequency (%f Hz) must be lower than the high (%f Hz)",
	app.lo_band, app.hi_band);
 }

 printf("-- Gain multiplier %g (%g%% == %g dB)\n",
	app.gain_mul, ((double)((int)(app.gain_mul * 1000.0 + 0.5))) / 10.0,
	20.0 * log10(app.gain_mul));

 if(HCW_FMT_BAD_FMT == app.c_format)
  app.c_format = HCW_FMT_PCM_FLT32;
 PrintCwaveFormat(app.c_format);
}

/* show help text
*/
static void help(void)
{
 printf("Usage: cwt [keys] input.wav output.cwave\n"
	"or cwt -t [keys] input.cwave\n"
	"or cwt -h\n"
	"where *.wav - 16-bit stereo Windows PCM wav file;\n"
	"      *.cwave - target special wave file contain complex samples.\n"
	"The input and output files MUST have correct extensions.\n"
	"Possible keys:\n"
	"-h, -? - print this text\n"
	"-g value - set gain multiplier for input samples\n"
        "-t - check input.CWAVE (w/o any output) for integrity\n"
	"-f[e][s][i] - set some FFT options (-f w/o letters take no effect):\n"
	"-fe - make total number of FFT points strictly even\n"
	"      (default - strictly odd, regardless real number of samples)\n"
	"-fs - make FFT procwssing slowly and safely (recommended for big files)\n"
	"-fi - prohibit using SIMD (SSE2) instructions (silly:);\n"
        "-rX freq - remove low/high frequences from the spectrum:\n"
        "(-rl freq - from 0 to freq, Hz; -rh freq - form freq, Hz to maximum)\n"
        "-rc or -rs or -r21 - same as '-rl 21 -rh 21000'\n"
	"-v[s][p] - verbose mode; [s] [p] modifiers set special FFTW output:\n"
	"-vs - print FFTW statistics; -vp - print FFTW plans\n"
	"-iX - set sample format for the output CWAVE file\n"
	"(-id == double; -is == short(16-bit);\n"
	" -im == short+float; -if == float)\n"
	"-j nCPU - set maximum number of loaded CPU (working threads)\n"
	"NOTE: most of options for -t will be ignored\n"
	);

 printf("\n"
	"-- The default values:\n"
	"-g %f, -j 1; CWAVE sample type -if\n",
	DEF_GAIN);

 printf("\n"
	"Copyright (C) 2010-2014 Rat and Catcher Tech.\n"
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
	"This program uses FFTW library written by Matteo Frigo and\n"
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
 // very special case
 if(app.isTestCRC)
 {
  CheckCrcProcess();
  return;				// and here the story finished...
 }

 // prepare for processing
 app.fpif = cfopen(app.nameif, "r", "input WAV-data");
 app.fpof = cfopen(app.nameof, "w", "output complex data");
 readWavHeader(app.fpif, &app.hcw.sample_rate, &app.hcw.n_samples);
 memcpy(&(app.hcw.magic[0]), HCW_MAGIC, sizeof(app.hcw.magic));
 app.hcw.hsize = sizeof(app.hcw);
 app.hcw.version = HCW_VERSION_BAD;
 app.hcw.format = app.c_format;
 app.hcw.n_channels = 2;
 app.hcw.k_M = -1;			// dumb FFT-based algorithm
 app.hcw.n_CRC32 = 0;
 app.hcw.k_beta = 0.0;			// dumb FFT-based algorithm
 writeCwaveHeader(app.fpof, &app.hcw);

 // ..processing..
 if(!fftw_init_threads())
  error("Can't initialized multi-tread FFTW");
 fftw_plan_with_nthreads(app.nCPU);

 app.nsFFT = app.isFFTeven?		// STRICTLY even or STRICTLY odd
	(app.hcw.n_samples + 1U) & (~01U) : app.hcw.n_samples | 01U;

 CalcFFTpass();				// FFT filter calculations

 if(app.isFFTsafe)
  ProcessFFT_Safe();
 else
  ProcessFFT();

 fftw_cleanup_threads();

 // finish
 app.hcw.n_CRC32 = crc32final(&app.tcrc);
 app.hcw.version = HCW_VERSION_CUR;
 rewind(app.fpof);
 writeCwaveHeader(app.fpof, &app.hcw);
 if(app.isVerbose)
  CalcTime(app.tstart);
 if(HCW_FMT_PCM_INT16 == app.hcw.format ||
	HCW_FMT_PCM_INT16_FLT32 == app.hcw.format)
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
}

/* check CRC of CWAVE
*/
static void CheckCrcProcess(void)
{
 unsigned crc;
 unsigned i;

 app.fpif = cfopen(app.nameif, "r", "input complex data");
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

/* the end...
*/
