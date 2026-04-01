/*
 * helpers.c -- the implementation file for
 * the helpers functions for the analitic (complex) signal
 * transform program.
 *
 * Copyright (C) 2010-2026 Rat and Catcher Tech.
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

#include "cw.h"

/* print an error message and terminate (no cleanup)
*/
void error(const char *fmt, ...)
{
 va_list ap;

 va_start(ap, fmt);
 // plain command line
 vfprintf(stderr, fmt, ap);
 fprintf(stderr, "\n");
 exit(1);
}

/* open file with check
*/
FILE *cfopen(const char *name, const char *mode, const char *msg)
{
 FILE *res = fopen(name, mode);

 if(NULL == res)
  error("Can't open %s for %s", name, msg);
 return res;
}

/* allocate memory with check
*/
void *cmalloc(size_t size, const char *msg)
{
 void *res = calloc(size, 1);

 if(NULL == res)
  error("Can't allocate %lu bytes for %s", (unsigned long)size, msg);
 return res;
}

/* make ftell() with check
*/
long cftell(FILE *fp)
{
 long res = ftell(fp);

 if(res < 0)
  error("Internal ftell() error");
 return res;
}

/* make fseek() to the position with check
*/
void cfseek(FILE *fp, long pos)
{
 long cpos;

 if(fseek(fp, pos, SEEK_SET))
  error("Internal fseek() error (pos = %ld)", pos);
 cpos = cftell(fp);
 if(cpos != pos)
  error("fseek() can't reach the requested pos (%ld != %ld)", cpos, pos);
}

/* create temporary file name (should be free())
*/

#if 0

char *ctempfile(void)
{
 char *t = getenv("TMP"), *res, uname[80];
 size_t l;
 static unsigned nn = 0;

 if(NULL == t)
 {
  t = getenv("TEMP");
  if(NULL == t)
   t = ".";
 }
 l = strlen(t);
 res = (char *)cmalloc(l + 80, "temp file name");
 strcpy(res, t);
 t = strchr(res, ';');
 if(NULL != t)
 {
  *t = '\0';
  l = strlen(res);
 }
 if(l)
 {
  if('\\' == res[l-1] || '/' == res[l - 1])
   res[l - 1] = '\0';
 }
 sprintf(uname, "\\cW-%016llX-%08X.tmp", (unsigned long long)GetCurrentProcessId(), nn++);
 strcat(res, uname);
 return res;
}

#else

char *ctempfile(const char *afile)
{
// make a filename with the same path as argument
// e.g. in the same folder as output file --
// -- it can be slowly on mechanical HDD but looks
// most reliable (output file system should to have
// enough free space but %TMP%/%TEMP% not)
#define MAX_UNAME   (80)
 static unsigned nn = 0;
 size_t len = afile? strlen(afile) : 0;
 char uname[MAX_UNAME + 2] = { 0 };
 char *res = (char *)cmalloc(len + sizeof(uname) + 2, "temp file name");

 *res = 0;

 _snprintf(uname, sizeof(uname) - 2, "CWT-%08X%08X.tmp", (unsigned)GetCurrentProcessId(), nn++);
 uname[sizeof(uname) - 1] = 0;

 if(len)
 {
  size_t ib = len;

  strcpy(res, afile);

  while(ib)
  {
   --ib;
   if('\\' == res[ib] || '/' == res[ib] || ':' == res[ib])
    break;

   res[ib] = '\0';
  }
 }

 strcat(res, uname);

 return res;
}

#endif

/* check file extension
*/
int checkFileExt(const char *fname, const char *ext)
{
 int j = (int)strlen(fname);
 const char *p;

 for(p = fname + j; *p != '.' && j != 0; --j, --p)
  /* nothing */;
 return *p != '.'? 0 : (!_stricmp(p + 1, ext));
}

/* read and check WAV PCM header
*/
void readWavHeader(FILE *fp, unsigned *srate, unsigned *nsamples, unsigned *byteps /* 2 or 3 */)
{
 unsigned char buf[44];
 unsigned bps = 2, ns = 0;

// the code extracted from flac examples (www.flac.org)
// [I agree, that the struct-based approach here is not so good ;))]
 if(fread(buf, 1, 44, fp) != 44)
 {
  error("ERROR: WAV header read error");
 }

 if( // 16 bit
        !memcmp(buf, "RIFF", 4) &&
        !memcmp(&buf[8], "WAVEfmt \020\000\000\000\001\000\002\000", 16) &&
        !memcmp(&buf[32], "\004\000\020\000data", 8))
 {
  bps = 2;
 }
 else if( // legacy 24 bit
        !memcmp(buf, "RIFF", 4) &&
        !memcmp(&buf[8], "WAVEfmt \020\000\000\000\001\000\002\000", 16) &&
        !memcmp(&buf[32], "\006\000\030\000data", 8))
 {
  bps = 3;
 }
 else
 {
  error("ERROR: invalid/unsupported WAVE file,\n"
        "       only 16bps or lagacy 24bps stereo WAVE in canonical form allowed");
 }

 ns = (((((((unsigned)buf[43] << 8) | buf[42]) << 8) | buf[41]) << 8) | buf[40]) / (2 /*ch*/ * bps);
 if(!ns)
 {
  error("ERROR: WAV file has not audio samples");
 }

 *srate = ((((((unsigned)buf[27] << 8) | buf[26]) << 8) | buf[25]) << 8) | buf[24];
 *nsamples = ns;
 *byteps = bps;
}

/* read and convert to double one sample
*/
void readWavSample(FILE *fp, unsigned byteps, double *ls, double *rs)
{
// Windows run on LSB machines only. Integer WAV file always LSB.
 signed short sbuf[2];
 int ibuf[2];
 const char emsg[] = "Input Wave PCM data read error";

 switch(byteps)
 {
  case 2:
   if(fread(sbuf, sizeof(short), 2, fp) != 2)
    error(emsg);

   *ls = (double)sbuf[0];
   *rs = (double)sbuf[1];
   break;

  case 3:
   if(fread(&ibuf[0], 1, 3, fp) != 3 || fread(&ibuf[1], 1, 3, fp) != 3)
    error(emsg);

   *ls = ((double)(ibuf[0] << 8)) / 65536.0;
   *rs = ((double)(ibuf[1] << 8)) / 65536.0;
   break;

  default:
   error("INTERNAL ERROR: Bad sample type");
   break;
 }
}

/* read complex wave (CWAWE) header
*/
void readCwaveHeader(FILE *fp, HCWAVE *hcw)
{
 if(fread(hcw, sizeof(HCWAVE), 1, fp) != 1)
  error("Can't read the Complex data header");

 if(memcmp(&(hcw -> magic[0]), HCW_MAGIC, sizeof(hcw -> magic)) != 0 ||
        hcw -> hsize != sizeof(HCWAVE) ||
        hcw -> version > HCW_VERSION_CUR ||
        hcw -> version == HCW_VERSION_BAD ||
        hcw -> n_channels != 2 ||
        hcw -> n_samples < 2)
  error("Bad or unsupported CWAVE header");

 switch(hcw -> format)
 {
  case HCW_FMT_PCM_DBL64:
  case HCW_FMT_PCM_INT16:
  case HCW_FMT_PCM_INT16_FLT32:
  case HCW_FMT_PCM_FLT32:
   break;
  default:
   error("Bad or unsupported CWAVE data format code '%u'", hcw -> format);
   break;
 }
}

/* write complex wave (CWAWE) header
*/
void writeCwaveHeader(FILE *fp, HCWAVE *hcw)
{
 if(fwrite(hcw, sizeof(HCWAVE), 1, fp) != 1)
  error("Output Complex Wave Header write error");
}

/* read complex data (for CRC check only)
*/
void readComplex(FILE *fp, unsigned t_format, TMP_CRC32 *tcrc)
{
 double buf[4];                 // it has max size for the all formats

 switch(t_format)
 {
  case HCW_FMT_PCM_DBL64:
   if(fread(buf, sizeof(double), 4, fp) != 4)
        error("Input Complex(double) data read error");
   crc32update(buf, sizeof(double) * 4, tcrc);
   break;

  case HCW_FMT_PCM_INT16:
   if(fread(buf, sizeof(short), 4, fp) != 4)
        error("Input Complex(short) data read error");
   crc32update(buf, sizeof(short) * 4, tcrc);
   break;

  case HCW_FMT_PCM_INT16_FLT32:
   if(fread(buf, sizeof(short) + sizeof(float), 2, fp) != 2)
        error("Input Complex(short+float) data read error");
   crc32update(buf, (sizeof(short) + sizeof(float)) * 2, tcrc);
   break;
   
  case HCW_FMT_PCM_FLT32:
   if(fread(buf, sizeof(float), 4, fp) != 4)
    error("Input Complex(float) data read error");
   crc32update(buf, sizeof(float) * 4, tcrc);
   break;

  default:
   error("Internal error 'readComplex()'");
   break;
 }
}

/* write complex data
*/
// Inline helpers
// double -> int16 w/o rounding
static INLINE signed short trunk_dbl(double d, long *n_clips)
{
 int t = (int)d;

 if(t < -32768)
 {
  ++(*n_clips); return -32768;
 }
 else
 {
  if(t > 32767)
  {
   ++(n_clips); return 32767;
  }
  else
  {
   return (signed short)t;
  }
 }
}

// double -> int16 with rounding
static INLINE signed short round_dbl(double d, long *n_clips)
{
 int t;

 if(d < 0.0)
  t = (int)(d - 0.5);
 else
  t = (int)(d + 0.5);

 if(t < -32768)
 {
  ++(*n_clips); return -32768;
 }
 else
 {
  if(t > 32767)
  {
   ++(n_clips); return 32767;
  }
  else
  {
   return (signed short)t;
  }
 }
}

// the main write function
void writeComplex(FILE *fp, double l_re, double l_im,
        double r_re, double r_im, const HCWAVE *hcw,
        long *l_clips, long *r_clips, TMP_CRC32 *tcrc,
        int is16b)
{
 double buf[4];
 short sbuf[4];
 char vbuf[(2 + 4) * 2];        // 2 * (sizeof(short) + sizeof(float))
 float fbuf[4];

 switch(hcw -> format)
 {
  case HCW_FMT_PCM_DBL64:
   buf[0] = l_re;
   buf[1] = l_im;
   buf[2] = r_re;
   buf[3] = r_im;

   crc32update(buf, sizeof(double) * 4, tcrc);

   if(fwrite(buf, sizeof(double), 4, fp) != 4)
        error("Output Complex(double) data write error");
   break;

  case HCW_FMT_PCM_INT16:
   if(hcw -> k_M > 0 && is16b)                  // Hilbert FIR + 16 bit samples
   {
    sbuf[0] = trunk_dbl(l_re, l_clips);         // l_re contain exact value
    sbuf[2] = trunk_dbl(r_re, r_clips);         // r_re contain exact value
   }
   else                                         // direct FFT
   {
    sbuf[0] = round_dbl(l_re, l_clips);         // should be rounded
    sbuf[2] = round_dbl(r_re, r_clips);         // should be rounded
   }
   sbuf[1] = round_dbl(l_im, l_clips);          // should be rounded
   sbuf[3] = round_dbl(r_im, r_clips);          // should be rounded

   crc32update(sbuf, sizeof(short) * 4, tcrc);
   
   if(fwrite(sbuf, sizeof(short), 4, fp) != 4)
    error("Output Complex(short) data write error");
   break;

  case HCW_FMT_PCM_INT16_FLT32:
   if(hcw -> k_M > 0 && is16b)                  // Hilbert FIR + 16 bit samples
   {
    sbuf[0] = trunk_dbl(l_re, l_clips);         // l_re contain exact value
    sbuf[1] = trunk_dbl(r_re, r_clips);         // r_re contain exact value
   }
   else                                         // direct FFT
   {
    sbuf[0] = round_dbl(l_re, l_clips);         // should be rounded
    sbuf[1] = round_dbl(r_re, r_clips);         // should be rounded
   }
   fbuf[0] = (float)l_im;
   fbuf[1] = (float)r_im;

   memcpy(&vbuf[0                                            ], &sbuf[0], sizeof(short));
   memcpy(&vbuf[sizeof(short)                                ], &fbuf[0], sizeof(float));
   memcpy(&vbuf[sizeof(short) + sizeof(float)                ], &sbuf[1], sizeof(short));
   memcpy(&vbuf[sizeof(short) + sizeof(float) + sizeof(short)], &fbuf[1], sizeof(float));

   crc32update(vbuf, (sizeof(short) + sizeof(float)) * 2, tcrc);

   if(fwrite(vbuf, sizeof(short) + sizeof(float), 2, fp) != 2)
    error("Output Complex(short/float) data write error");
   break;

  case HCW_FMT_PCM_FLT32:
   fbuf[0] = (float)l_re;
   fbuf[1] = (float)l_im;
   fbuf[2] = (float)r_re;
   fbuf[3] = (float)r_im;

   crc32update(fbuf, sizeof(float) * 4, tcrc);

   if(fwrite(fbuf, sizeof(float), 4, fp) != 4)
    error("Output Complex(float) data write error");
   break;

  default:
   error("Internal error 'writeComplex()");
   break;
 }
}

/* print CWAVE format information
*/
void PrintCwaveFormat(unsigned cw_format)
{
 printf("-- CWAVE format: ");
 switch(cw_format)
 {
  case HCW_FMT_PCM_DBL64:
   printf("[Re/Im double (2*64 bits) samples, -32768.0..32767.0]");
   break;
 case HCW_FMT_PCM_INT16:
  printf("[Re/Im short (2*16 bits) samples, -32768..32767]");
  break;
 case HCW_FMT_PCM_INT16_FLT32:
  printf("[Re/Im short+float (16+32 bits) samples, -32768..32767]");
  break;
 case HCW_FMT_PCM_FLT32:
  printf("[Re/Im float (2*32 bits) samples, -32768.0..32767.0]");
  break;
 default:
  printf("Unknown complex sample format");
  break;
 }
 printf("\n");
}

/* get number of CPU in the system
*/
int GetCpuNumber(void)
{
 SYSTEM_INFO si;

 GetSystemInfo(&si);
 return (int)(si.dwNumberOfProcessors);
}

/* create auto non-signaling event with check
*/
HANDLE ccreate_event(void)
{
 HANDLE res = CreateEvent(NULL, FALSE, FALSE, NULL);

 if(NULL == res)
  error("Can't create the Windows Event object");
 return res;
}

/* calculate processing time
*/
void CalcTime(time_t tstart)
{
 time_t tot, t;
 unsigned sec;

 time(&tot);
 tot -= tstart;
 sec = (unsigned)tot;
 printf("-- Processing time: ");
 if(0 != (t = tot / (24 * 60 * 60)))    // days
 {
  printf("%u Days ", (unsigned)t);
  tot %= (24 * 60 * 60);
 }
 if(0 != (t = tot / (60 * 60)))         // hours
 {
  printf("%2u Hr ", (unsigned)t);
  tot %= (60 * 60);
 }
 if(0 != (t = tot / 60))                // minutes
 {
  printf("%2u Min ", (unsigned)t);
  tot %= 60;
 }
 printf("%u Sec (%u Sec total)\n", (unsigned)tot, sec); // seconds
}

/* the end...
*/

