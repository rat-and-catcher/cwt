# cwt -- Complex (audio) Wave Tool

In ./retro folder you can to find old sources of our **old** cwave production tool (cw).
This is Windows only CLI program; and have lots of limitations (input -- 16 bit stereo
WAV w/o any metadata and chunks different from 'fmt ' and 'data'; ASCII charset for
file names and so). Its need FFTW3 (fftw.org) library to compile and work.

./retro/fftw-win contain some fftw-3.3.4 Visual Studio solution as well as our **cw -only**
patch for original fftw-3.3.4 sources (./retro/fftw-win/fftw-3.3.4/kernel/assert.c).
Please note, that the solution and patch are far from ideal, but somewhat works. Note also,
that our assert.c **is cw -only and inacceptable for generic fftw usage**.
To make a build (VS 2013) you should download original fftw-3.3.4.tar.gz from
ftp.fftw.org, unpack it, copy content of our fftw-win/fftw-3.3.4 over original folder
(replace kernel/assert.c) and rebuild fftw-3.3.4-libs/fftw-3.3.4-libs.sln.

./retro/linux folder contain the same code with the only "direct FFT" _(L. Marple Jr.,
Computing the discrete-time analytic signal via FFT, IEEE Transactions on Signal Processing,
47(2001), 2600-2603)_ algorithm for U\*x systems. You, probably, should review and fix our
makefile.

Both implementations are scratch, but (with some fantasy skill) can be classified as
"research quality software".

# ..we do hope to provide better tool in the repo..
