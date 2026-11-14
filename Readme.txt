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

NFS_factory <polynomial file name> <saved relation file name> <output file name> [-rels N] [-threads N] [-gpu N] [-verbose]

The polynomial file describes the polynomial used for the current factorization, in GGNFS format with a few extra options, detailed below. The relation file contains the NFS relations stored by the sieving procedure, in the GGNFS/msieve relation format, or one of two binary formats that save space (see "Relation file formats" below). The program will run until all relations in the file have been processed. All found relations are appended to the output file specified.

Several relation files can be processed in one run, which sets up the prime products only once and keeps filling the same remainder tree across the files:

NFS_factory <polynomial file name> -list <list file name> <output file name> [-rels N] [-threads N] [-gpu N] [-verbose]

The list file names one relation file per line (empty lines and lines starting with # are skipped), and the files are read in that order. With -rels N, the run stops once N relations have been found (counted when their batch is checked, before they are written): it finishes the batches it has read (up to the number of batches of the pipeline option) and the 3LP candidates found so far, so it usually finds more than N. Relations found by an earlier, interrupted run of the same input count towards N.

After each batch, when its relations are checked and written, the position of its last relation is saved in <output file name>.resume: the relation file or list, the file reached, how many lines or binary relations of it are done, the relations found so far, and the size of the output file. If the run is interrupted, running it again with the same relation file or list (and output file) continues from there, and gives the same output as an uninterrupted run; relations written after the saved position are removed from the output file first, since their batch is checked again, and batches that were read but not yet written (with a pipeline) are read again. A run that has done all its files leaves a resume file that says so; remove it to run the same input again. The resume file is replaced in one step, so a run can be stopped at any time. A run does not continue, and stops with an error, if the output file is missing or shorter than when the resume file was written (the relations counted would not all be there), or if the list now names another relation file at the saved position: restore the files, or remove the resume file to start again. With "3lp: save", the size of the candidate file is saved too, and the candidates written after the saved position are removed from it in the same way. 3LP candidates waiting for the second stage when a run is interrupted are not saved, so their relations (about 1 % of the relations of that part of the input) are lost.

At the end of a run, <output file name>.status holds a summary for scripts, one "key value" per line: found (relations found by this run), total (including earlier runs of the same input), target and target_reached, files_done and files, complete (all files done), and error (a relation file could not be read, or the relations found could not be written; the exit code is then -1). The status file of an earlier run is removed when a run starts, so a run that stops before its end (wrong arguments, a polynomial file that is refused, no GPU, a resume file that cannot be used) leaves none, and has exit code -1.

If the output file cannot be written (a full disk, for example), the run stops with an error, and the resume file stays at the last batch whose relations were written, so the run continues from there when it is started again. A read error of a relation file stops the run in the same way. Records with b = 0 in a binary relation file (a part of a file filled with zeros) are skipped and counted.

A typical use is a driver script that runs NFS_factory on a list of relation files with -rels set to the number of relations still needed, and then runs the filtering, linear algebra and square root. With one relation file per run instead, the driver adds up the relations found from the status files.

The GPU build (NFS_factory_gpu, see "Compiling the code.txt") does the product and remainder trees, and the cofactors of the relation values at their leaves, on an NVIDIA GPU, and is meant to be run as one process, on a list of relation files. It checks the batches in a pipeline: while one thread reads the next batch, several batches are in the trees on the GPU at once, the checks of the cofactors of the batches before them run on the CPU, and a writer thread writes their relations, both in the order of the batches, so that the output is the same as without the pipeline. For .ab files (and with "format: cado") the factor lists of the other side of the relations found are rebuilt by a batch check of their norms of its own on the GPU: in trees of 2^17 norms, against the primes up to 2^28, with the cofactors of one or two larger primes split as on the checked side; for the norms with more of them, the cofactor and the rest of the norm are factored apart by Pollard rho. The writer finishes the lines when their tree is done, so relations are written up to about 2^17 relations after they are found. That check needs a prime product of its own, which takes about 8 seconds to set up when the first relations need it. On SNFS-240 (4 .ab files of 139M relations each, with nothing else on the GPU) a run took 46 s, against 57 s when the norms were factored by Pollard rho alone, and 37 s with no factor lists of the other side at all. The number of batches in progress is chosen from the free memory of the GPU (see pipeline below), and the build runs on NVIDIA GPUs of compute capability 7.5 (RTX 20 series) and later. The CPU build runs one relation file per process as before, typically as several processes in parallel.

-threads N sets the number of threads of the checks on the CPU: the cofactors of the relation values (primality tests, Pollard rho), and the factoring of the relations found (default 1 in the CPU build; in the GPU build 8, or the number of hardware threads if that is less). The GPU build does both on the GPU, except tinyqs, relation values of more than 256 bits, and a few rare cases, so it needs few threads: on SNFS-240 (RTX 4070 Ti SUPER) one thread is as fast as eight. The output does not depend on it.

-gpu N selects the GPU of the GPU build, numbered from 0 (default 0); the one used is printed at the start. The GPU build needs a CUDA GPU of compute capability 7.5 or later, and stops with a message saying so if there is none; on such a computer, use the CPU build.

At the start, the program prints the setup: the GPU, the prime products, the size of the batches (relations per batch and tree levels), the batches in progress at once, the threads, and the number and size of the relation files. During the run, one status line shows the time since the reading started, the relation file of the last batch checked, the relations checked so far and per second, and the relations found (with those of the 3LP second stage) and per second. On a console the line is rewritten in place twice a second; when the output goes to a file, it is written as a line of its own once a minute, with a line for each relation file of a list. Messages, such as the 3LP candidates found or a skipped part of a file, are printed above it. -verbose also prints, as earlier versions did, the size of the top of the product tree and the relations found for each batch.

Relation file formats
---------------------

The format of a relation file is chosen from its name:
- .bin: the binary format, (int)a, (uint)b, (uint)number of primes, (uint)primes per relation, all 32-bit. The primes are those of the side that is not batch checked (normally the rational side, shared by all numbers of the relation set), which are printed for the relations found. About 50 bytes per relation for SNFS-240, 84 % of it the primes.
- .ab: only a and b. A header of 16 bytes ("NFAB", then (uint)version 1, (uint)the side whose primes the relations were saved for, and (uint)flags, bit 0 set if the relations were validated), followed by (int64)a, (uint32)b per relation: 12 bytes. The factor list of the side that is not batch checked is rebuilt from its norm for the relations found (on the GPU in the GPU build, along with the factoring of the checked side; the CPU build factors those norms on the CPU, which made a one-thread SNFS-220 run about 30 % slower, and costs more for larger norms). The output is the same as from the .bin file of the same relations. a has 64 bits since the siever's a can reach 2^31 and more, which the .bin format cannot hold.
- .3lp: 3LP candidates (see "Three large primes").
- anything else: the GGNFS/msieve text format, a,b:rational primes:algebraic primes.
The list of a run can mix the formats.

The .ab format is meant for relation sets that are read from a hard disk: for SNFS-240, a .bin file is read at about 225 MB/s from a hard disk, which then limits the GPU build, while the same relations as a .ab file are a quarter of the size, and the GPU is the limit again (about 1.4 times faster with the polynomial of a 2026 SNFS-240 run with 3LP, and about twice as fast with parameters where the GPU checks more relations per second). From an SSD, or the page cache, both formats run at the same speed. The full SNFS-240 relation set takes 0.53 TB as .ab files, instead of 2.2 TB.

Utilities\RelationsToBinaryFormat converts text and .bin files to .ab files. Binary files written by older versions of it can hold relations whose |a| is 2^31 or more with a wrong a (2^32 too large or too small): their factor lists do not match the norms, and NFS_factory leaves such relations out of the output (it checks the saved factor list of every relation it finds against the norm), but cannot find them. Converting such a file to .ab needs the polynomial to find the right a for each relation.

The polynomial file contains the polynomial pair and factorization parameters in a GGNFS format, and NFS_factory needs the polynomial pair, rlim/alim, lbpr/lbpa and mfbr/mfba to run. If the file gives the number n, the polynomials are checked against it: the algebraic polynomial must be 0 at the root m (or -Y0/Y1) of the rational polynomial modulo n, and the run stops if it is not. Without n, the run says that it could not check. The batch smoothness checking will test the relations for all primes up to rlim/alim, and any remiainder is checked if it is a single large prime if it is less than lbpr/lbpa, or double-large prime if it is less than mfbr/mfba (or 2 lbpr/lbpa, if mfbr/mfba is larger, see "Three large primes" below). One more parameter is required to control which polynomial should be checked for smoothness. Batch smoothness checking is done on the rational side if “side: 0” is specified in the polynomial file, and on the algebraic side if “side: 1” is specified.

A few optional parameters can also be given in the polynomial file:
- format: msieve (default) or cado. With "format: cado" all primes are printed on both sides of each relation, with multiplicity, as CADO-NFS requires. The factor list of the side that is not batch checked is then rebuilt from its norm, and relations where the stored factors do not match the norm are skipped. With "format: msieve" that factor list is printed as it was read (rebuilt for .ab files), after checking that its primes divide the norm, and small primes are left out, since msieve finds them itself.
- print: primes smaller than this are not printed for the batch checked side (default 1000, which msieve handles). Not needed with "format: cado".
- maxa, maxb: relations with |a| larger than maxa, or b larger than maxb, are skipped. There is no limit by default.
- maxrels: stop after this many relations have been read from the relation file. Useful for testing and timing.
- treelevels: the maximum number of levels in the remainder tree (default 20). The number of levels is otherwise chosen from the size of the prime product and the relation values. Each level less halves the memory use of the tree, at some cost in speed. The GPU build uses one level more when the memory of the GPU and of the computer still allows six batches in progress at once (for SNFS-240 on a 16 GB GPU: 20 instead of 19 levels, about 4 % faster, about 1.8 GB of GPU memory per batch); treelevels still caps it, and pipeline in the polynomial file turns it off.
- 3lp: auto (default) or save. Only used when mfbr/mfba is more than 2 lbpr/lbpa, see "Three large primes" below.
- 3lpfile: the file that "3lp: save" appends the 3LP candidates to. The default is the output file name, followed by a dot, the name of the relation file and ".3lp".
- 3lptreelevels: the maximum number of levels in the remainder tree when running on a 3LP candidate file (default 22). treelevels is not used for such runs.
- pipeline: the number of batches in progress at once, each with its own tree (on the GPU, about 1 GB each for SNFS-240): one being read, one being checked, and the others in the trees (up to four at once). With the gpu backend the default is as many as fit in 70 % of the free memory of the GPU, up to 6, and in half of the free memory of the computer (about 0.15 GB each for SNFS-240); the choice is printed at the start. 6 (four trees on the GPU at once) is the fastest on an RTX 4070 Ti SUPER, about 7 % faster than 5, and 2 still overlaps the reading and the checks with the GPU. Otherwise the default is 1, which checks each batch in turn.
- backend: how the product and remainder trees are computed. mpir (the default of the CPU build) uses MPIR's products. gpu (the default of the GPU build, and only available there) uses an integer number-theoretic transform on the GPU, from the int_ntt folder. ntt uses the same transform on the CPU, and computes the same numbers as the GPU in every node of the trees; it is much slower than MPIR, and is the reference for testing the GPU version. All of them give exactly the same relations.

Three large primes
------------------

With mfbr/mfba larger than 2 lbpr/lbpa (at most 3 lbpr/lbpa), relations with three large primes on the batch checked side are also found. After the batch check against the primes up to rlim/alim, a remainder that is too large for two large primes, but at most mfbr/mfba bits, at least rlim^3/alim^3 and composite, can be the product of three large primes. These remainders are the 3LP candidates. Typically less than 1 % of them are relations, so instead of factoring each of them, they are batch checked in a second stage against the product of all primes below 2^lbpr/2^lbpa. Since that prime product is large (290 Mbit for 2^29 and a degree 4 polynomial, which takes 10 to 20 seconds to set up), the second stage is run when enough candidates have been collected to fill its remainder tree, which with treelevels 20 is about 500000 candidates.

- "3lp: auto" (the default): the second stage runs during the run whenever its tree is full, and on the remaining candidates at the end, so that each run finishes its relation file on its own, and many runs can be started in parallel on different relation files. This needs memory for the second remainder tree (about 550 MB more for lpba 29).
- "3lp: save": the candidates are appended to a candidate file (see 3lpfile) and not checked. Candidate files from many runs can be concatenated, and are checked by running NFS_factory with a candidate file, which must have the extension .3lp, instead of a relation file. Only the second stage is then run, and the found relations are appended to the output file. Use the same polynomial file, so that lpbr/lpba is the same.

The second stage costs less per candidate when its tree is full, and a run on a small relation file may end with a small, and therefore expensive, last batch of candidates. Saving the candidates and checking them together avoids that. A run on a candidate file can also use a larger tree, since such runs are few and do not run in parallel with many others: up to 3lptreelevels (default 22) levels. For lpba 29, 22 levels hold about 2 million candidates, needs about 2 GB of memory, and was about 20 % faster on the second stage than 20 levels.

The candidate file has the binary relation format, with the remainder added to each relation as two 64-bit numbers (low 64 bits first). Candidates from .ab files are saved without primes, and the factor list of their other side is rebuilt when they are checked. Candidates whose |a| is 2^31 or more do not fit, and are counted as skipped.

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

At the end of a run, the program prints how the time was spent on the different steps, and the peak memory use. Each step is given in Gcycles (of the processor's time stamp counter) and as a share of the time of the whole run. The time of a step is added up over the threads doing it, and with a pipeline several threads work at once (up to four batches in the trees, the reader, and the thread of the checks), so a step can take more than 100 %, and the shares add up to more than 100 %. The waits of the pipeline show what limits the run: the thread of the checks, which takes the batches in order, is the limit when it seldom waits for the trees, and the trees (the GPU) are the limit when the checks wait for them while the reader waits for a free slot.

Lines in a GGNFS/msieve relation file that cannot be parsed are skipped and counted, as is an incomplete last line, which a siever that is stopped may leave behind. A binary or .ab relation file that ends in the middle of a relation is also detected.

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


