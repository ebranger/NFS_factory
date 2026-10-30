NFS_factory: An implementation of Coppersmith's factorization factory.
Erik Branger


Introduction
------------

NFS_factory is my attempt at implementing Coppersmith's factorization factory method, in a hopefully efficient manner. The program takes pre-computed and stored Number Field Sieve (NFS) relations, that are valid for the shared polynomial to be used by all factorizations, performs a batch smoothness test, identifies the relations and obtains the full factorization ofthe relations, and outputs them. It is designed to take NFS relations in a GGNFS/msieve format, and to output the results in the same format. msieve, or alternatively CADO-NFS, can then perform the post-processing, i.e. the linear algebra and square root, to finish the factorization. The program also requires a polynomial file specifying the polynomial pair of the factorization and reads those in a GGNFS format.

My original inspiration to write this program was the SNFS factorization of 17 Mersenne numbers using the method, demonstrating that performance gains were possible using a batch factorization approach. I have so far used the program to factor about 200 SNFS numbers of difficulty 220 digits sharing a rational polynomial. These factorizations also show that the factoring effort per number using Coppersmiths approach can be reduced by more than 50%, provided that the initial sieving cost can be amortized over sufficiently many numbers.

Developing this code has been more of a programming exercise and an excuse to learn more about the Number Field Sieve and Coppersmith’s factorization factory idea, rather than an attempt to make the fastest or cleanest version possible. However, since I have found no other publicly available implementation of a batch smoothness check for Coppersmith’s method, this code may still be of some interest after all. 

This documentation corresponds to version 0.1 of the program, there is still much work to be done before it will reliably work on any range of numbers of relevance.


Using NFS_factory
-----------------

The first step in using Coppersmiths factorization factory is to sieve for relations and store them. This must be done using another program, and I have mainly used the Franke-Kleinjung lattice siever gnfs-lasieve. Some details regarding this is provided later in this document.

Once enough relations have been stored, it is time for NFS_factory to start working. The program is invoked from the command line taking three arguments:

NFS_factory <polynomial file name> <saved relation file name> <output file name>

The polynomial file describes the polynomial used for the current factorization, in GGNFS format with a few extra options, detailed below. The relation file contains the NFS relations stored by the sieving procedure, and accepts either the GGNFS/msieve relation format, or a binary format I have used to save some space. The program will run until all relations in the file have been processed. All found relations are appended to the output file specified.

The polynomial file contains the polynomial pair and factorization parameters in a GGNFS format, and NFS_factory needs the polynomial pair, rlim/alim, lbpr/lbpa and mfbr/mfba to run. The batch smoothness checking will test the relations for all primes up to rlim/alim, and any remiainder is checked if it is a single large prime if it is less than lbpr/lbpa, or double-large prime if it is less than mfbr/mfba (or 2 lbpr/lbpa, if mfbr/mfba is larger, see "Three large primes" below). One more parameter is required to control which polynomial should be checked for smoothness. Batch smoothness checking is done on the rational side if “side: 0” is specified in the polynomial file, and on the algebraic side if “side: 1” is specified.

A few optional parameters can also be given in the polynomial file:
- format: msieve (default) or cado. With "format: cado" all primes are printed on both sides of each relation, with multiplicity, as CADO-NFS requires. The factor list of the side that is not batch checked is then rebuilt from its norm, and relations where the stored factors do not match the norm are skipped. With "format: msieve" that factor list is printed as it was read, and small primes are left out, since msieve finds them itself.
- print: primes smaller than this are not printed for the batch checked side (default 1000, which msieve handles). Not needed with "format: cado".
- maxa, maxb: relations with |a| larger than maxa, or b larger than maxb, are skipped. There is no limit by default.
- maxrels: stop after this many relations have been read from the relation file. Useful for testing and timing.
- treelevels: the maximum number of levels in the remainder tree (default 20). The number of levels is otherwise chosen from the size of the prime product and the relation values. Each level less halves the memory use of the tree, at some cost in speed.
- 3lp: auto (default) or save. Only used when mfbr/mfba is more than 2 lbpr/lbpa, see "Three large primes" below.
- 3lpfile: the file that "3lp: save" appends the 3LP candidates to. The default is the output file name, followed by a dot, the name of the relation file and ".3lp".
- 3lptreelevels: the maximum number of levels in the remainder tree when running on a 3LP candidate file (default 22). treelevels is not used for such runs.

Three large primes
------------------

With mfbr/mfba larger than 2 lbpr/lbpa (at most 3 lbpr/lbpa), relations with three large primes on the batch checked side are also found. After the batch check against the primes up to rlim/alim, a remainder that is too large for two large primes, but at most mfbr/mfba bits, at least rlim^3/alim^3 and composite, can be the product of three large primes. These remainders are the 3LP candidates. Typically less than 1 % of them are relations, so instead of factoring each of them, they are batch checked in a second stage against the product of all primes below 2^lbpr/2^lbpa. Since that prime product is large (290 Mbit for 2^29 and a degree 4 polynomial, which takes 10 to 20 seconds to set up), the second stage is run when enough candidates have been collected to fill its remainder tree, which with treelevels 20 is about 500000 candidates.

- "3lp: auto" (the default): the second stage runs during the run whenever its tree is full, and on the remaining candidates at the end, so that each run finishes its relation file on its own, and many runs can be started in parallel on different relation files. This needs memory for the second remainder tree (about 550 MB more for lpba 29).
- "3lp: save": the candidates are appended to a candidate file (see 3lpfile) and not checked. Candidate files from many runs can be concatenated, and are checked by running NFS_factory with a candidate file, which must have the extension .3lp, instead of a relation file. Only the second stage is then run, and the found relations are appended to the output file. Use the same polynomial file, so that lpbr/lpba is the same.

The second stage costs less per candidate when its tree is full, and a run on a small relation file may end with a small, and therefore expensive, last batch of candidates. Saving the candidates and checking them together avoids that. A run on a candidate file can also use a larger tree, since such runs are few and do not run in parallel with many others: up to 3lptreelevels (default 22) levels. For lpba 29, 22 levels hold about 2 million candidates, needs about 2 GB of memory, and was about 20 % faster on the second stage than 20 levels.

The candidate file has the binary relation format, with the remainder added to each relation as two 64-bit numbers (low 64 bits first).

The largest candidates are the least likely to be products of three primes below 2^lbpr/2^lbpa, so mfbr/mfba a few bits less than 3 lbpr/lbpa keeps most of the 3LP relations for a much smaller cost. For an SNFS-240 polynomial with lpba 29 and alim 2^29/20, a sample of 8.4 million relations gave:

  mfba   candidates   3LP relations   cost of 3LP (first and second stage)
   87     508350       2543            13.1 Gcycles
   86     452407       2533            12.5 Gcycles
   85     399754       2493            10.3 Gcycles
   84     350156       2377             9.4 Gcycles
   83     303519       2151             8.6 Gcycles
   82     260064       1847
   81     219025       1428

The cost does not include setting up the prime product of the second stage, about 20 Gcycles per run whatever mfba is. For comparison, the whole run without 3LP took about 97 Gcycles and found 19587 relations, about 5 Mcycles per relation. The last few bits cost much more than that per relation found: going from 85 to 87 bits costs more than 50 Mcycles per extra relation, and from 84 to 85 about 8 Mcycles. mfba = 3 lpba - 3 (84 for lpba 29) is a reasonable choice, or 3 lpba - 2 when more relations are worth more CPU time. Without 3LP, use mfbr/mfba = 2 lbpr/lbpa as before.

At the end of a run, the program prints how the time was spent on the different steps, and the peak memory use.

Lines in a GGNFS/msieve relation file that cannot be parsed are skipped and counted, as is an incomplete last line, which a siever that is stopped may leave behind. A binary relation file that ends in the middle of a relation is also detected.

Once NFS_factory has processed and found enough relations, msieve or CADO-NFS is needed for the postprocessing, and this code was made with msieve in mind.

Parameter selection
-------------------

Selecting good parameters is important for the best possible performance of the NFS algorithm. However, I only know of one work where actual parameters have been optimized, for the factorization of 17 Mersenne numbers with a size of 1000 to 1200 bits. These numbers are far too large to be factored without access to a computer cluster and will not be of relevance to an individual having access to a few computers. Thus, if you intend to use this code you may have to start by figuring out good parameters for the factorizations you intend to perform. While I cannot help providing good parameters, I can give a few considerations for some of the parameters.    

The two main parameters to choose is the polynomial degree and the large prime bounds/factor-base limits. By choosing the degree of the algebraic polynomial, it is possible to affect the work required to sieve for relations and for batch-checking them. As an example, consider the factorization of several numbers sharing a rational polynomial, but with different algebraic polynomials. Choosing a low degree for the algebraic polynomials will make it more expensive to sieve for rational relations, and larger factor base will be required on the rational side. However, the norms on the algebraic side will be smaller and more likely to factor over the factor base. As a result, fewer numbers need to be tested for smoothness, making this part faster. It may also be possible to somewhat reduce the factor base bound on the algebraic side, which will further speed up the batch smoothness checking. Thus, this will make the sieving step harder and the batch smoothness step easier. If there are a lot of numbers to factor, this may be worthwhile, since the sieving cost will be amortized across all numbers. This will also ensure that fewer relations need to be stored, which can potentially save a lot of storage memory. 

The size of the factor base needs to be set so that the sieving step can find enough relations to factor all numbers you intend to factor. This will require identifying the worst-case polynomial for your numbers and assessing its smoothness probability when batch-factoring the relations. 


Sieving for relations
---------------------

Although NFS_factory could in principle be used to obtain the relations for the shared polynomial, doing so would be horribly inefficient. A sieve is required to find the relations with reasonable effort. For the factorizations I have done, the relations were obtained using the Franke-Kleinjung lattice siever gnfs-lasieve4I15e. Ideally, this program should be modified to allow sieving over only one polynomial. An alternative that I have used is to run the siever in a non-standard way, which works but with several caveats.

The idea is to use a polynomial pair, where one of the polynomials is the one of interest, and the other one is trivial, such as k*x+l for some small k or l, perhaps around 10. With the two polynomials defined this way, one can calculate the common root, and from that find the number to be factored that this corresponds to, which must be supplied to the siever for it to start. If the algebraic polynomial is the shared one in the factorization factory, the rational polynomial can then be set to something easy when sieving. If the rational polynomial is the shared one, gnfs-lasieve is capable of sieving with a degree-1 “algebraic” polynomial, though I have encountered some strange issues and crashes doing so.

From trial and error, I have found that using a simple polynomial for the side to be ignored of k*x-l for k and l around 10 seems to work best. For too small k and l a lot of special-q are skipped with a “sched_pathology” error and very few relations are found. Also, when sieving with two degree-1 polynomials the siever will sometimes crash on startup, and sometimes it will work, I have found no pattern to this. 

Once the relations are obtained, it is convenient to split the relations into several smaller files, based on the size of the largest prime in the relation. Since NFS_factory will process an entire file before exiting, this will allow the work to be sub-divided into smaller, more manageable pieces. The post-processing can then be done when enough relations have been collected, without having to process the remaining stored relations. The files can also be processed in the order of increasing maximum prime, which effectively uses the lowest large prime bound possible for the shared polynomial. It may also be useful to remove duplicates in each such file, to further reduce the work needed to process one file. 


Implementation details
----------------------

The core of NFS_factory is an implementation of Bernstein’s batch smoothness checking algorithm using remainder trees. It also uses scaled remainder trees, to replace divisions/modulo operations with multiplies, which should further increase the speed. The large integer arithmetic is done using MPIR, which should be the only dependency of the code to external libraries. To calculate the homogeneous polynomial values/norms for the relations, I have recycled some polynomial code I wrote earlier. I’m sure that there are faster alternatives out there for bigint polynomial calculations, but the code works and is not the main speed bottleneck. 

Once the remainder tree has been calculated, any values that are less than mfbr/mfba is checked to see if it is a 1lp or a 2lp relation. The 2lp numbers are split with Pollard rho in native 64-bit arithmetic, and those that rho does not split, or that have more than 63 bits, with a fast tiny quadratic sieve implementation, originally coming from msieve, from a branch with code from Robert Silverman (cofactorize_siqs.c). Every valid relation that is found by NFS_factory is then fully factored with trial division and Pollard rho in native arithmetic up to 128 bits, and Pollard rho with MPIR above that (Smooth_factor.cpp).


Credits
-------

The implementation of the batch smoothness checking has been inspired by the same algorithms used in msieve and CADO-NFS. Some code has also been used coming from msieve, which may have originated from other codes.


A big thanks to all the people of mersenneforum, surprisingly much of what I know about integer factorization comes from your discussions, and the patient explanations you have provided to various questions over the years.


References
----------

The original idea by D. Coppersmith to save relations to be used in later factorizations is his paper: D. Coppersmith, Modifications to the number field sieve, Journal of Cryptology 6(3), 169–180 (1993).

The first successful implementation of Coppersmith’s factorization factory that I know of is found in: T. Kleinjung, J. Bos and A. Lenstra, Mersenne factorization factory, ASIACRYPT 2014

The batch smoothness checking algorithm is described in: D. Bernstein, How to find smooth parts of integers, (2004).

The idea of scaled remainder trees is described in: D. Bernstein, Scaled remainder trees, (2004).

This code has greatly been inspired by and intended to work together with msieve by Jason Papadopoulos. It can be found at https://sourceforge.net/projects/msieve/

The code has also taken some inspiration from CADO-NFS, found at http://cado-nfs.gforge.inria.fr/

The code is also intended to work together with the lattice siever by J. Franke and T. Kleinjung, the source code can be found in the GGNFS code by C. Monico, which can be found at https://sourceforge.net/projects/ggnfs/


