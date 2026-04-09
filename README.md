# cwt -- Complex (audio) Wave Tool

In ./retro folder you can to find old sources of our **old** cwave production tool (cw).
This is Windows only CLI program; and have lots of limitations (input -- 16 or 24 bit stereo
WAV w/o any metadata and chunks different from 'fmt ' and 'data'; ASCII (or OEM/DOS) charset for
file names and so). Its need FFTW3 (fftw.org) library to compile and work. You shoud to use
cwt???.exe -h for list of program options; its rather big.

We create some (initially for Windows) build / test fftw3 (fftw 3.3.10 now) build / test
infrastructure in the folder ./fftw-3.3.10-build. To build it, you should upack original
fftw-3.3.10 source tarball alongside this folder. You also need cmake (minimum V3.15)
in your path. Windows batch files alongside CMakeLists.txt in the folder could make
your life easy in the way. Note, that them apply some our patch in the build process,
but you can still to build original untouched FFTW. Note, that our projects needs our
patches to compile. If you have Visual Studio 2017, you can simple run win-bcmake.cmd
to build all possible FFTW3 DLL versions; or win-bcmake-static.cmd to buld all possible
FFTW3 static libs, needed for our project(s). If you have another build enviroment,
please review our batch files. Also note about targes names. To test the build,
we adapted for Windows FFTW test infrastructure, see ./fftw-3.3.10-build/tests for detail.

./retro/linux folder contain the same code with the only "direct FFT" _(L. Marple Jr.,
Computing the discrete-time analytic signal via FFT, IEEE Transactions on Signal Processing,
47(2001), 2600-2603)_ algorithm for U\*x systems. You, probably, should review and fix our
makefile.

Both implementations are scratch, but (with some fantasy skill) can be classified as
"research quality software".

Please note, that 32-bit version practically unusable for "direct FFT" algorithm.
Also note, that (1) SSE2 looks optimal for performance / reliability and (2) multithread
run give not so big performance bonus while significantly increase amount of consumed
memory.

# ..we do hope to provide better tool in the repo..
