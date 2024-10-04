# cwt -- Complex (audio) Wave Tool

In ./retro folder you can to find old sources of our **old** cwave prodution tool (cw).
This is Windows only CLI program; and have lots of limitations (input -- 16 bit stereo
WAV w/o any metadata and chunks different from 'fmt ' and 'data'; ASCII charset for
file names and so). Its need FFTW3 (fftw.org) library to compile and work.
./retro/linux folder contain the same code with the only "direct FFT" (L. Marple(?))
algorithm for U\*x systems. Both implementations are scratch, but (with some fantasy
skill) can be classified as "research quality software".

# ..we do hope to provide better tool in the repo..
