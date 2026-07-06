# cwt -- Complex (audio) Wave Tool

In ./retro folder you can to find old sources of our **old** cwave production tool (cw).
This is Windows only CLI program; and have lots of limitations (input -- 16 or 24 bit stereo
WAV w/o any metadata and chunks different from `fmt ` and `data`; ASCII (or OEM/DOS) charset for
file names and so). Its need FFTW3 (fftw.org) library to compile and work. You shoud to use
cwt???.exe -h for list of program options; its rather big.

We create some (initially for Windows) build / test fftw3 (fftw 3.3.11 now)
infrastructure in the folder `./fftw-3.3.11-build/`. To build it, you should unpack original
fftw-3.3.11 source tarball alongside this folder. You also need cmake (minimum V3.15)
in your path. Windows batch files alongside `CMakeLists.txt` in the folder could make
your life easy in the way. Note, that they apply some our patch in the build process,
but you can still to build original untouched FFTW. Note, that our projects needs our
patches to compile. If you have Visual Studio 2017, you can simple run win-bcmake.cmd
to build all possible FFTW3 DLL versions; or win-bcmake-static.cmd to build all possible
FFTW3 static libs, needed for our project(s). If you have another build environment,
please review our batch files. Also note about targets names. To test the build,
we adapted for Windows FFTW test infrastructure, see `./fftw-3.3.10-build/tests/` for detail.

`./retro/linux/` folder contain the same code with the only "direct FFT" _(L. Marple Jr.,
Computing the discrete-time analytic signal via FFT, IEEE Transactions on Signal Processing,
47(2001), 2600-2603)_ algorithm for U\*x systems. You, probably, should review and fix our
makefile.

Both implementations are scratch, but (with some fantasy skill) can be classified as
"research quality software".

Please note that 32-bit version practically unusable for "direct FFT" algorithm.
Also note, that (1) SSE2 looks optimal for performance / reliability and (2) multithread
run give not so big performance bonus while significantly increase amount of consumed
memory.

## Notes about FFTW3

For our project we need some DFT. Our choice -- FFTW3. It fast, well tested and uses
widely. But we found the problem. While FFTW run out of memory, it silently call `abort()`.
This can happen in FFT execute functions, and (probably) can happen in planner too.
For example take 60 minutes CD-audio, single channel. It has 158,760,000 samples.
In fftw_complex representation it has about 2.5 GB. It looks rather big, but not so big.
But to make full-complex in-place DFT FFTW3 required 20+GB! It is real example for us,
moreover, memory needs may be even bigger, depending of problem size factorization.
We are sure, that FFTW uses optimal amount of memory for solving the problem. We only want
to have sufficient diagnostics in out of memory situation; unpredictable `abort()` from
library function is not a graceful solution. In the case we want to get possibility
to have a chance to print corresponding message and to exit with predicable code. As
maximum we, of course, want to return to the state which was before failed function call.

As quick-and-dirty solution we replace abort() call in FFTW `fftw_assertion_failed()`
function to the code, which generate memory exception and catch the exception in
our code. This works under Windows and don't works in most Linux distributions.
For the moment this was enough for us. Of course, we see some more sufficient solutions,
but they need some efforts. Probably we will do it, probably not.

Now we have sufficient FFTW3 V3.3.10 and V3.3.11 build environment for Windows-only and workable
workaround of out-of-memory problem; its not a fork.
