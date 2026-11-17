This folder holds small programs for preparing relation sets for NFS_factory, and a driver script that runs a whole
factorization with NFS_factory and msieve.

A typical way to build a relation set: sieve on one side (usually the rational side), check and split the siever output into
numbered files with RelationFilterSplit, remove the duplicate relations of each file with remdups, and convert the files to .bin or
.ab with RelationsToBinaryFormat. batchmsieve.py then factors numbers with that relation set.


RelationFilterSplit
-------------------

RelationFilterSplit <input file> <output name> <polynomial file> <r|a> <smallest prime> <largest prime> <first file number> <last file number>

Checks the relations of a relation file and splits them into numbered files by the largest prime of the side that was sieved, so
that the relation set is in files of a size that NFS_factory can use, in the order of that prime. The purpose is to take sieved
relations for the shared side found with another program, and put them in a format suitable for NFS_factory.

- The input is a text file (the GGNFS/msieve format a,b:rational primes:algebraic primes, as the siever writes it) or a .bin file
  (see RelationsToBinaryFormat), and the output has the same format: <output name><number>.out for text, <output name><number>.bin
  for binary. Converting between the formats is done only by RelationsToBinaryFormat. .ab files hold no factor lists to split by,
  and are refused: split the text or .bin files, and convert the results.
- r or a is the side that was sieved. Its factor list is kept (sorted), and the list of the other side is dropped, which saves
  space: a text line is written as a,b:primes: (r) or a,b::primes (a).
- The primes from the smallest prime up to, but not including, the largest prime are divided evenly over the files numbered from
  the first file to the last file, by the largest prime of each relation. A relation whose largest prime is below the range (or
  that has no primes listed) goes to <output name>_small, and one whose largest prime is at or above it to <output name>_large.
  The primes are given in decimal or as a power of two, such as 2^28.
- A relation is kept if b > 0, a and b are coprime, |a| <= maxa and b <= maxb (if the polynomial file gives them), and every prime
  listed for the sieved side divides its norm on that side, with multiplicity. The polynomial file is read for its polynomials
  (c0, c1, ..., and Y0 and Y1, or m) and maxa/maxb. Lines that are not relations, and other relations that fail these checks, are
  skipped and counted; an incomplete last line (from a siever that was stopped) is skipped.
- The relations are appended to the output files, so that the output of many siever runs can be split into the same files. A run
  that is interrupted and started again appends its relations a second time; remdups removes them. A read or write error stops
  the run with exit code -1.
- The smallest prime is typically the smallest special-q sieved, and the largest prime is determined by lbpr or lpba on the 
  shared side.

The two ranges allow a relation set to be extended. For example, an SNFS-220 set sieved with special-q from 10M to 2^28 with a
lbpa of 2^28 into the files 1 to 360, and later extended with the special-q from 2^28 to 2^29 into the files 361 to 720:

  RelationFilterSplit sieve1.out spairs poly.txt r 10000000 2^28 1 360
  RelationFilterSplit sieve2.out spairs poly.txt r 2^28 2^29 361 720

Should any relations in the raw sieve files have a largest prime outside the range they will be saved to sparis_small and 
spairs_large. Both can be used as relation files of their own, or left out.

It needs MPIR and Polynomial.cpp/h of the NFS_factory folder. On Linux, make utils builds it (see "Compiling the code.txt"). On
Windows, from a Visual Studio command prompt in the main folder (with MPIR in third_party\mpir as for the main project):

  cl /O2 /MT /EHsc /std:c++17 /I NFS_factory /I third_party\mpir\include Utilities\RelationFilterSplit.cpp NFS_factory\Polynomial.cpp /link /LIBPATH:third_party\mpir\lib\x64\Release mpir.lib mpirxx.lib


RelationsToBinaryFormat
-----------------------

RelationsToBinaryFormat <input file> <output file> b [side]
RelationsToBinaryFormat <input file> <output file> a [side]
RelationsToBinaryFormat <input file> <output file> t <side> <polynomial file>

Converts a text relation file to the binary format (.bin, b), a binary file back to text (t), or a text or .bin file to the format
of only a and b (.ab, a), a quarter of the size of .bin, for which NFS_factory rebuilds the factor list of the side that is not
batch checked (see "Relation file formats" in the main Readme).

- side is the side that the relations were sieved on: 0 or r for rational (the default), 1 or a for algebraic. Its factor list is
  the one kept in a .bin file, written by t (as a,b:primes: or a,b::primes), and saved in the header of an .ab file. A .bin file
  does not record its side (NFS_factory takes its primes for the side that it does not batch check), so give the same side when
  converting a .bin file on.
- Before converting, the first 100 relations are checked against the side, so that a swapped side is caught. If the check fails,
  nothing is written and the exit code is -1.
  - Text input (b and a): of the lines with only one factor list (as RelationFilterSplit writes them), most must have it on that
    side; a few relations with an empty list are accepted. Siever output with both lists holds the primes of both sides, and either
    side can be kept.
  - .bin to text (t): the primes must divide the norms on that side for most of the relations, which needs the polynomial file
    (only the polynomials are read from it), since a .bin file does not record its side.
  - .bin to .ab (a): not checked here, but NFS_factory refuses an .ab file whose side is the side that its polynomial file batch
    checks.
- The .bin format holds a in 32 bits, so relations with |a| of 2^31 or more are skipped when converting to it (the .ab format
  holds them); binary files written by older versions of the program can hold such relations with a wrong a.
- Lines that are not relations are skipped and counted; a read or write error stops it with exit code -1.

It needs no libraries (on Windows: cl /O2 /EHsc /std:c++17 Utilities\RelationsToBinaryFormat.cpp).


remdups
-------

remdups (by Greg Childers, distributed with GGNFS and gnfs-lasieve) removes duplicate relations from a relation file. Run it on the
files of RelationFilterSplit, before they are converted, so that the factory runs do not spend time on duplicate relations.
remdups is not included here, since it is under the GPL; it can be found in the contrib/remdups folder of the GGNFS sources (for
example https://github.com/radii/ggnfs/tree/master/contrib/remdups), and is discussed at
https://www.mersenneforum.org/showthread.php?t=13498. remdups.c builds on Windows; remdups4.c is a filter (stdin to stdout) that
needs a POSIX system.


batchmsieve.py
--------------

A stand-alone driver that factors one number by SNFS with a relation set: NFS_factory (the CPU build) or NFS_factory_GPU checks
the relation files for smoothness on the other side, and msieve does the filtering, linear algebra and square root. It is based
on factmsieve.py by Brian Gladman, and needs Python 3.9 or later (gmpy2 is used for the primality test of n if it is installed).
Before running it, there are a few file paths that needs to be set early in the script, to find NFS_factory and Msieve,
for the CPU and GPU versions of the programs. Its usage is as follows:

  python batchmsieve.py NAME -dir FOLDER [options] -gpu 	start NAME (NAME.poly in the current folder) on relations stored 
								in FOLDER on the GPU.
  python batchmsieve.py NAME [options]                		resume NAME, options can be changed from the previous run.
  python batchmsieve.py -h                            		the full help.

- NAME.poly holds n, type: snfs, the polynomial (c0, c1, ... and m, or Y0 and Y1), the sieving parameters (rlim, alim, lpbr,
  lpba, mfbr, mfba, rlambda, alambda), and side: the side that NFS_factory batch checks (1, the algebraic side, for relations
  sieved on the rational side). The limits of the sieved side must match how the relation files were made.
- -dir names a folder of relation files (.bin, .out or .ab; files of the same name in several formats are the same relations, and
  one of them is used), checked in the order of their names (spairs2 before spairs10). It can be given several times.
- -gpu checks with NFS_factory_GPU (one process on a list of all files, which stops when enough relations are found) and runs the
  linear algebra with a CUDA build of msieve; -cpu (the default) runs NFS_factory on one file per process, -threads processes at
  a time. A job can switch between them when it is restarted.
- -minrels sets the relations to collect before msieve's filtering is tried (default 80M for a large prime bound of 29 bits on
  the sieved side, half of it per bit less). If the filtering fails, 10 % more relations are collected and it is tried again. A
  file minrels.txt in the job folder takes precedence, and is read again in every round.
- -threads, -gpunum, -filter (msieve filtering options) and -exclude (relation files to leave out) are the other options. All
  options are saved in NAME.cfg, so that a restart needs none of them.
- The paths of the programs (NFS_FACTORY, NFS_FACTORY_GPU, MSIEVE, MSIEVE_GPU) are set at the top of the script, relative to its
  folder (by default its subfolder Filer) or absolute.

Each step is skipped when its output exists, so running the same command again resumes an interrupted job: the batch checking
(progress in NAME.job.resume and NAME.done), the linear algebra (from msieve's checkpoint), and the square root, which writes a
summary to s<difficulty>-NAME.txt.
