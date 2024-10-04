/*
 * The structure of complex wave (CWAVE) file
 * (in general -- for the analitic signal representation)
 * Copyright (C) 2010 Rat and Catcher Tech.
 *
 * This file provide a format information about CWAVE (Complex Audio Wave)
 * files format.
 *
 * Redistribution of this file, with or without modification, are permitted
 * provided that the following conditions are met:
 *
 * - Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer listed
 * in this license in the documentation and/or other materials
 * provided with the distribution.
 *
 * - Neither the name of the copyright holders nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _cwave_h_
#define _cwave_h_

// Version 1
// -- Obsolete -- header of complex wave (cwave) file
typedef struct tagHCWAVE_V1
{
 char magic[8];                 // signature
 unsigned hsize;                // header size in bytes
 unsigned version;              // format version
 unsigned format;               // format of samples
 unsigned n_channels;           // number of channels
 unsigned n_samples;            // number of samples
 unsigned sample_rate;          // sample rate of initial file
 int k_M;                       // Hilbert FIR filter order
 unsigned pad0;                 // NOT USED (par for align k_beta)
 double k_beta;                 // filter parameter
} HCWAVE_V1;

// Version 2
// -- Current -- header of complex wave (cwave) file
typedef struct tagHCWAVE_V2
{
 char magic[8];                 // signature
 unsigned hsize;                // header size in bytes
 unsigned version;              // format version
 unsigned format;               // format of samples
 unsigned n_channels;           // number of channels
 unsigned n_samples;            // number of samples
 unsigned sample_rate;          // sample rate of initial file
 int k_M;                       // Hilbert FIR  filter order (-1 for direct FFT)
 unsigned n_CRC32;              // CRC32 of the data part
 double k_beta;                 // Hilbert FIR filter parameter
} HCWAVE_V2;

typedef HCWAVE_V2 HCWAVE;

// ...and it's field description:
// .magic:
#define HCW_MAGIC               "cPLXwAVE"      /* 8 characters w/o '\0' */
// .version:
#define HCW_VERSION_BAD         (0)             /* mark bad or incomplete file */
#define HCW_VERSION_V1          (1)             /* V1 */
#define HCW_VERSION_V2          (2)             /* V2 */
#define HCW_VERSION_CUR         (HCW_VERSION_V2) /* current version */
// .format:
// not a format
#define HCW_FMT_BAD_FMT         (0xFFFFFFFF)
// double Re(L), Im(L), Re(R), Im(R), [-32768.0..32767.0]
#define HCW_FMT_PCM_DBL64       (0)
// signed short Re(L), Im(L), Re(R), Im(R), [-32768..32767]
#define HCW_FMT_PCM_INT16       (1)
// signed short Re(.), float Im(.), [-32768..32767]/[-32768.0..32767.0]
#define HCW_FMT_PCM_INT16_FLT32 (2)
// float Re(L), Im(L), Re(R), Im(R), [-32768.0..32767.0]
#define HCW_FMT_PCM_FLT32       (3)
// .. no more yet ..

#endif          // def _cwave_h_

/* the end...
*/

