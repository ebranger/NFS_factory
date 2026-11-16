# batchmsieve.py (batchmsieve_driver.py) - stand-alone SNFS factorization factory driver for
# NFS_factory (CPU) or NFS_factory_GPU, and msieve
#
# Based on factmsieve.py, Copyright (c) 2010, Brian Gladman, which is a
# conversion of factmsieve.pl, Copyright 2004, Chris Monico.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with
# or without modification, are permitted provided that the
# following conditions are met:
#
#    Redistributions of source code must retain the above
#    copyright notice, this list of conditions and the
#    following disclaimer.
#
#    Redistributions in binary form must reproduce the
#    above copyright notice, this list of conditions and
#    the following disclaimer in the documentation and/or
#    other materials provided with the distribution.
#
#    The names of its contributors may not be used to
#    endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
# THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

"""
Factor one number by SNFS with a factorization factory: relations pre-sieved
on one side (usually the rational side, for the root m of the polynomial) are
batch checked for smoothness on the other side by NFS_factory (CPU) or
NFS_factory_GPU, and msieve does the filtering, linear algebra and square root.

Usage
  python batchmsieve_driver.py NAME -dir FOLDER [options]    start NAME
  python batchmsieve_driver.py NAME [options]                resume NAME
  python batchmsieve_driver.py -h                            this help

  NAME.poly must exist in the current folder; all job files are written there.

Options (all are remembered in NAME.cfg, so a restart needs none of them)
  -dir FOLDER     folder with the pre-sieved relation files: every file in it
                  ending in .bin, .out (text) or .ab (NFS_factory_GPU format)
                  is used, whatever its name. Files with the same name and
                  another ending are the same relations in another format.
                  The files are checked in name order (file2 before file10).
                  Give -dir several times if the files are spread over several
                  folders; on a restart, -dir replaces the folders remembered.
                  The files must belong to NAME.poly (same root m and the
                  pre-sieved side's limits).
  -minrels N      relations to collect before msieve filtering is tried
                  (written to NAME.minrels). Default: 80M for a large prime
                  bound of 29 bits on the pre-sieved side (lpbr, or lpba with
                  side: 0), halved for every bit less: 2^(lpb-29) * 80M.
                  If filtering fails, more relation files are checked and
                  filtering is tried again with 10% of the target more
                  relations. A file minrels.txt in the job folder takes
                  precedence over NAME.minrels (for all jobs there); it is
                  read every round and never deleted.
  -gpu            batch checking with NFS_factory_GPU and linear algebra with
                  msieve_gpu, which stores the matrix only once on the GPU
                  (msieve option single_copy=1, GPU_LA_OPTIONS in the script)
  -cpu            all on the CPU (default for a new job)
  -threads N      NFS_factory processes on the CPU (one relation file each),
                  threads of NFS_factory_GPU and msieve. Default: half the
                  logical processors (usually the physical cores).
  -gpunum N       the GPU to use (default 0)
  -filter "OPTS"  extra msieve filtering options. Default: target_density=120
                  for degree 4 and higher, none for degree 3 (it failed there).
                  -filter "" for none.
  -exclude NAME   a relation file in the -dir folders not to use (its name,
                  with or without the ending; all formats of it are left
                  out). Give it several times for several files; on a
                  restart, -exclude replaces the files remembered.

Switching between -cpu and -gpu when a job is restarted is fine: the job
continues from the relation file where the other one stopped. The files
already checked are listed in NAME.done.

Programs
  NFS_factory (for -cpu), NFS_factory_GPU and msieve_gpu (CUDA build, for
  -gpu), and msieve (always). Their paths are set at the top of this script
  (NFS_FACTORY, NFS_FACTORY_GPU, MSIEVE, MSIEVE_GPU), relative to the folder
  of this script or absolute; by default they are in its subfolder Filer.

NAME.poly
  n: <number>, type: snfs, the coefficients c0..cd, the rational root as m:
  (or Y1: 1 and Y0: -m), skew: (optional), the sieving parameters rlim,
  alim, lpbr, lpba, mfbr, mfba, rlambda, alambda, and side: the side
  NFS_factory batch checks, 1 (algebraic, for relations pre-sieved on the
  rational side) or 0 (rational, for relations pre-sieved on the algebraic
  side). The limits of the pre-sieved side must match how the relation
  files were made.

Steps, each skipped when its output exists, so the same command resumes:
  1. Setup (NAME.dat missing): NAME.ini, NAME.fb and NAME.dat.
  2. Batch checking (NAME.dat.cyc missing) until msieve filtering succeeds.
     The relations are added to NAME.dat; progress is in NAME.job.resume
     and NAME.done.
     CPU: NFS_factory on one relation file per process, a round of -threads
     files at a time; filtering is tried after the round that reaches the
     target, and after a failure once 10% of the target more is found.
     GPU: the relation files left (.ab where there is one) are written to
     NAME.list, and one NFS_factory_GPU process checks them, appending to
     NAME.dat itself (progress in NAME.dat.resume), until the target is
     reached. If filtering fails it continues for 10% of the target more.
  3. Linear algebra (NAME.dat.dep missing): msieve -nc2, or -ncr from the
     checkpoint NAME.dat.chk; with -gpu, msieve_gpu with single_copy=1.
  4. Square root: msieve -nc3, and the summary s<difficulty>-NAME.txt.
"""

import datetime
import glob
import math
import multiprocessing
import os
import platform
import random
import re
import shutil
import subprocess
import sys
import time

# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCRIPT_NAME = os.path.basename(__file__)   # in the help and the log

# Paths of the programs: relative to the folder of this script (here its
# subfolder Filer), or absolute. Set them to where the programs are on this
# computer. '.exe' is added on Windows; '/' works on Windows too.
NFS_FACTORY = 'Filer/NFS_factory'          # CPU build
NFS_FACTORY_GPU = 'Filer/NFS_factory_GPU'  # CUDA build, for -gpu
MSIEVE = 'Filer/msieve'          # CPU build: filtering and square root, and LA with -cpu
MSIEVE_GPU = 'Filer/msieve_gpu'  # CUDA build: linear algebra with -gpu

# msieve_gpu options for the linear algebra: single_copy=1 stores the matrix
# only once on the GPU, so larger matrices fit in its memory.
GPU_LA_OPTIONS = 'single_copy=1'
# When filtering fails, filtering is tried again once this fraction of the
# target more relations has been found (CPU and GPU).
MORE_RELS = 0.10
# Relation target for all jobs in a folder; takes precedence over NAME.minrels.
MINRELS_ALL = 'minrels.txt'
# Default target: DEFAULT_MINRELS relations for a large prime bound of
# DEFAULT_LPB bits on the pre-sieved side, doubling per bit.
DEFAULT_MINRELS = 80000000
DEFAULT_LPB = 29
VERBOSE = True
LARGEP = 3          # large primes per side, for the summary

# NFS_factory sometimes ends with an access violation after its output is
# complete; that exit code is not treated as an error.
IGNORED_EXIT_CODES = (-1073741819, 3221225477)

EXE_SUFFIX = '.exe' if sys.platform.startswith('win') else ''
PRIORITY = 0x00000040 if sys.platform.startswith('win') else 0   # idle priority

# Set from the command line and NAME.cfg in main()
CFG = {}


class Job:
    """File names of one factorization."""

    def __init__(self, name):
        self.name = name
        self.poly = name + '.poly'
        self.cfg = name + '.cfg'
        self.log = name + '.log'
        self.ini = name + '.ini'
        self.fb = name + '.fb'
        self.dat = name + '.dat'
        self.cyc = name + '.dat.cyc'
        self.dep = name + '.dat.dep'
        self.chk = name + '.dat.chk'
        self.resume = name + '.job.resume'
        self.minrels = name + '.minrels'
        self.minrels_all = os.path.join(os.path.dirname(name), MINRELS_ALL)
        self.list = name + '.list'
        self.done = name + '.done'           # relation files already checked
        # NFS_factory_GPU writes to NAME.dat, so its files are named after it
        self.gpu_resume = self.dat + '.resume'
        self.gpu_status = self.dat + '.status'
        self.cpu_out = name + '.T'           # NFS_factory output, NAME.T<j>


# --------------------------------------------------------------------------
# Utilities
# --------------------------------------------------------------------------

def output(job, s, console=True):
    if console:
        print(s, flush=True)
    with open(job.log, 'a') as f:
        f.write(datetime.datetime.today().strftime('%a %b %d %H:%M:%S %Y ') + s + '\n')


def die(job, s):
    output(job, '-> ' + s)
    sys.exit(1)


def delete_file(fn):
    try:
        os.remove(fn)
    except OSError:
        pass


def exe_path(exe):
    """Full path of a program set at the top of the script."""
    p = os.path.normpath(os.path.join(SCRIPT_DIR, exe))
    if EXE_SUFFIX and not p.lower().endswith(EXE_SUFFIX):
        p += EXE_SUFFIX
    return p


def check_programs(job, exes):
    """Exit with a clear message if a program is missing."""
    names = {NFS_FACTORY: 'NFS_FACTORY', NFS_FACTORY_GPU: 'NFS_FACTORY_GPU',
             MSIEVE: 'MSIEVE', MSIEVE_GPU: 'MSIEVE_GPU'}
    missing = [e for e in exes if not os.path.isfile(exe_path(e))]
    if missing:
        for e in missing:
            output(job, '-> Error: {0:s} not found ({1:s} = {2!r})'
                   .format(exe_path(e), names[e], e))
        die(job, 'The paths of the programs must be set {0:s}.'.format(PATHS_SET_IN))


PATHS_SET_IN = ('at the top of {0:s} (NFS_FACTORY, NFS_FACTORY_GPU, MSIEVE, MSIEVE_GPU), '
                'relative to {1:s} or absolute'.format(os.path.abspath(__file__), SCRIPT_DIR))


def start_exe(job, exe, args, cwd=None):
    output(job, '-> {0:s} {1:s}'.format(os.path.basename(exe), ' '.join(args)), console=VERBOSE)
    try:
        return subprocess.Popen([exe_path(exe)] + args, cwd=cwd, creationflags=PRIORITY)
    except OSError as e:
        die(job, 'Error: could not start {0:s} ({1:s}). The paths of the programs are set '
            '{2:s}.'.format(exe_path(exe), str(e), PATHS_SET_IN))


def run_msieve(job, args, gpu=False):
    """Run msieve from its own folder on the files of this job."""
    files = ['-s', job.dat, '-l', job.log, '-i', job.ini, '-nf', job.fb]
    files = [a if a.startswith('-') else os.path.abspath(a) for a in files]
    if gpu:
        files = ['-g', str(CFG['gpunum'])] + files
    exe = MSIEVE_GPU if gpu else MSIEVE
    p = start_exe(job, exe, files + ['-t', str(CFG['threads'])] + args,
                  cwd=os.path.dirname(exe_path(exe)))
    return p.wait()


def count_lines(fn):
    n = 0
    with open(fn, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            n += block.count(b'\n')
    return n


def append_file(src, dst):
    with open(src, 'rb') as i, open(dst, 'ab') as o:
        shutil.copyfileobj(i, o, 1 << 24)


def is_probable_prime(n):
    try:
        import gmpy2
        return gmpy2.is_prime(n, 25)
    except ImportError:
        pass
    if n < 4:
        return n > 1
    d, s = n - 1, 0
    while not d & 1:
        d, s = d >> 1, s + 1
    for _ in range(25):
        x = pow(random.randint(2, n - 2), d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def proc_speed():
    """Processor speed in GHz, 0 if unknown."""
    try:
        if sys.platform.startswith('win'):
            from winreg import OpenKey, QueryValueEx, HKEY_LOCAL_MACHINE
            h = OpenKey(HKEY_LOCAL_MACHINE, 'HARDWARE\\DESCRIPTION\\System\\CentralProcessor\\0')
            return 1e-3 * float(QueryValueEx(h, '~MHz')[0])
        with open('/proc/cpuinfo') as f:
            m = re.search(r'cpu MHz\s+:\s+([0-9.]+)', f.read())
        return 1e-3 * float(m.group(1)) if m else 0.0
    except Exception:
        return 0.0


def read_keys(fn, sep=None):
    """'key value' lines (or 'key: value' with sep=':') of a file, as a dict;
    a key given more than once (dir) gives a list."""
    keys = {}
    if os.path.exists(fn):
        with open(fn) as f:
            for line in f:
                parts = [s.strip() for s in line.strip().split(sep, 1)]
                if parts and parts[0]:
                    val = parts[1] if len(parts) > 1 else ''
                    if parts[0] in keys:
                        old = keys[parts[0]]
                        keys[parts[0]] = (old if isinstance(old, list) else [old]) + [val]
                    else:
                        keys[parts[0]] = val
    return keys


# --------------------------------------------------------------------------
# Settings (command line and NAME.cfg)
# --------------------------------------------------------------------------

OPTIONS = {'-dir': 'dir', '-minrels': 'minrels', '-threads': 'threads', '-gpunum': 'gpunum',
           '-filter': 'filter', '-exclude': 'exclude'}


def parse_args(args):
    """Returns (NAME or None, the options given) or exits with the help."""
    names, opts = [], {}
    i = 0
    while i < len(args):
        a = args[i]
        key = '-' + a.lstrip('-') if a.startswith('-') else a
        if key in ('-h', '-help'):
            print(__doc__.replace('batchmsieve_driver.py', SCRIPT_NAME))
            sys.exit(0)
        if key in ('-gpu', '-cpu'):
            opts['mode'] = key[1:]
        elif key in OPTIONS:
            if i + 1 >= len(args):
                print('Error: {0:s} needs a value'.format(a))
                sys.exit(1)
            i += 1
            if key == '-dir':
                opts.setdefault('dir', []).append(os.path.abspath(args[i]))
            elif key == '-exclude':
                opts.setdefault('exclude', []).append(stem(args[i]))
            elif key in ('-minrels', '-threads', '-gpunum'):
                try:
                    opts[OPTIONS[key]] = int(float(args[i]))   # 80e6 is fine too
                except ValueError:
                    print('Error: {0:s} needs a number'.format(a))
                    sys.exit(1)
            else:
                opts[OPTIONS[key]] = args[i]
        elif a.startswith('-'):
            print('Error: unknown option {0:s}; see python {1:s} -h'.format(a, SCRIPT_NAME))
            sys.exit(1)
        else:
            names.append(a)
        i += 1
    if len(names) != 1:
        print(__doc__.replace('batchmsieve_driver.py', SCRIPT_NAME))
        sys.exit(0 if not args else 1)
    return names[0], opts


def load_settings(job, opts):
    """The settings of NAME.cfg, updated with the options given, saved again."""
    cfg = read_keys(job.cfg, ':')
    for k in ('dir', 'exclude'):
        if isinstance(cfg.get(k), str):
            cfg[k] = [cfg[k]]
    for k in ('threads', 'gpunum'):
        if k in cfg:
            cfg[k] = int(cfg[k])
    changed = {k for k, v in opts.items() if k != 'minrels' and cfg.get(k) != v}
    cfg.update({k: v for k, v in opts.items() if k != 'minrels'})
    cfg.setdefault('mode', 'cpu')
    cfg.setdefault('threads', max(1, multiprocessing.cpu_count() // 2))
    cfg.setdefault('gpunum', 0)
    if not cfg.get('dir'):
        die(job, 'Error: give the folder with the relation files with -dir FOLDER')
    with open(job.cfg, 'w') as f:
        for k in ('mode', 'threads', 'gpunum', 'filter'):
            if k in cfg:
                print('{0:s}: {1}'.format(k, cfg[k]), file=f)
        for d in cfg['dir']:
            print('dir: {0:s}'.format(d), file=f)
        for x in cfg.get('exclude', []):
            print('exclude: {0:s}'.format(x), file=f)
    return cfg, changed


def natural_key(s):
    """Sort key that puts file2 before file10."""
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]


def stem(path):
    """A relation file's name without folder and extension: the same relations
    in .bin, .out and .ab format have the same stem."""
    return os.path.splitext(os.path.basename(path.strip()))[0]


def relation_files(gpu):
    """[(stem, path)] of the relation files (.bin, .out, .ab) in the -dir
    folders, in natural order of the names. For the CPU .out is used rather
    than .bin (it cannot read .ab); NFS_factory_GPU prefers .ab (a converter
    writes its file under that name only once it is complete), then .out,
    then .bin. Other files (such as .ab.tmp) and the -exclude names are ignored."""
    rank = {'.ab': 0, '.out': 1, '.bin': 2} if gpu else {'.out': 1, '.bin': 2}
    exclude = set(CFG.get('exclude', []))
    best = {}
    for d in CFG['dir']:
        if not os.path.isdir(d):
            continue
        for name in os.listdir(d):
            ext = os.path.splitext(name)[1].lower()
            path = os.path.join(d, name)
            if ext not in rank or not os.path.isfile(path):
                continue
            k, r = stem(name), rank[ext]
            if k in exclude:
                continue
            if k not in best or r < best[k][0]:
                best[k] = (r, path)
    return sorted(((k, p) for k, (_, p) in best.items()), key=lambda f: natural_key(f[0]))


def read_done(job):
    """The stems of the relation files already checked (NAME.done)."""
    if not os.path.exists(job.done):
        return set()
    with open(job.done) as f:
        return {line.strip() for line in f if line.strip()}


def add_done(job, stems):
    with open(job.done, 'a') as f:
        for s in stems:
            print(s, file=f)


# --------------------------------------------------------------------------
# Polynomial
# --------------------------------------------------------------------------

REQUIRED_PARAMS = ('rlim', 'alim', 'lpbr', 'lpba', 'mfbr', 'mfba', 'rlambda', 'alambda')


def read_poly(job):
    """Read and check NAME.poly. Returns a dict with n, m, coeffs (list, index =
    power), degree, skew, difficulty and params (the other 'key: value' lines)."""
    params, coeffs = {}, {}
    with open(job.poly) as f:
        for line in f:
            line = re.sub('#.*', '', line).strip()
            if ':' not in line:
                continue
            key, val = (s.strip() for s in line.split(':', 1))
            if re.fullmatch(r'c\d+', key):
                coeffs[int(key[1:])] = int(val)
            elif val:
                params[key] = val

    if 'n' not in params:
        die(job, 'Error: n not given in {0:s}'.format(job.poly))
    n = int(params['n'])
    if params.get('type') != 'snfs':
        die(job, 'Error: the factory only handles SNFS polynomials (type: snfs)')
    if 'm' in params:
        m = int(params['m'])
    elif int(params.get('Y1', 0)) == 1:
        m = -int(params['Y0'])
    else:
        die(job, 'Error: {0:s} needs a rational root m (or Y1: 1 and Y0: -m)'.format(job.poly))
    if not coeffs:
        die(job, 'Error: no polynomial coefficients in {0:s}'.format(job.poly))

    degree = max(coeffs)
    if 'deg' in params and int(params['deg']) != degree:
        die(job, 'Error: deg: {0:s} but the highest coefficient is c{1:d}'
            .format(params['deg'], degree))
    c = [coeffs.get(i, 0) for i in range(degree + 1)]
    if math.gcd(*c) > 1:
        die(job, 'Error: the coefficients have a common factor {0:d}'.format(math.gcd(*c)))
    value = sum(ci * m ** i for i, ci in enumerate(c))
    if value == 0 or value % n:
        die(job, 'Error: f(m) is not a multiple of n')
    if is_probable_prime(n):
        die(job, 'Error: n is probably prime')
    missing = [p for p in REQUIRED_PARAMS if p not in params]
    if missing:
        die(job, 'Error: {0:s} is missing {1:s}'.format(job.poly, ', '.join(missing)))
    if params.get('side') not in ('0', '1'):
        die(job, 'Error: {0:s} needs side: 1 (batch check the algebraic side; relations '
            'pre-sieved on the rational side) or side: 0 (the other way round)'.format(job.poly))

    skew = float(params.get('skew', 0))
    if not skew:
        skew = abs(c[0] / c[degree]) ** (1.0 / degree) if c[0] else 1.0
    return {'n': n, 'm': m, 'coeffs': c, 'degree': degree, 'skew': skew,
            'difficulty': len(str(abs(value))), 'params': params}


def presieved_lpb(poly):
    """('lpbr' or 'lpba', its value) of the pre-sieved side: the side that
    NFS_factory does not batch check (side: 1 means it checks the algebraic
    side, so the relations were pre-sieved on the rational side)."""
    key = 'lpbr' if poly['params']['side'] == '1' else 'lpba'
    return key, int(poly['params'][key])


def default_minrels(poly):
    """DEFAULT_MINRELS at a large prime bound of DEFAULT_LPB bits on the
    pre-sieved side, doubling per bit."""
    return int(DEFAULT_MINRELS * 2.0 ** (presieved_lpb(poly)[1] - DEFAULT_LPB))


def filter_options(poly):
    if 'filter' in CFG:
        return CFG['filter'].split()
    return ['target_density=120'] if poly['degree'] >= 4 else []


# --------------------------------------------------------------------------
# Steps
# --------------------------------------------------------------------------

def setup(job, poly):
    """Start a new factorization: msieve input files and an empty NAME.dat."""
    output(job, '-> Setting up a new factorization')
    for f in ([job.resume, job.done, job.list, job.gpu_resume, job.gpu_status]
              + glob.glob(job.cpu_out + '*')):
        delete_file(f)
    p = poly['params']
    with open(job.ini, 'w') as f:
        print(poly['n'], file=f)
    with open(job.fb, 'w') as f:
        print('N {0:d}'.format(poly['n']), file=f)
        print('SKEW {0:<10.2f}'.format(poly['skew']), file=f)
        for i in range(poly['degree'], -1, -1):
            print('A{0:d} {1:d}'.format(i, poly['coeffs'][i]), file=f)
        print('R1 1', file=f)
        print('R0 {0:d}'.format(-poly['m']), file=f)
        print('FAMAX {0:s}'.format(p['alim']), file=f)
        print('FRMAX {0:s}'.format(p['rlim']), file=f)
        print('SALPMAX {0:d}'.format(2 ** int(p['lpba'])), file=f)
        print('SRLPMAX {0:d}'.format(2 ** int(p['lpbr'])), file=f)
    with open(job.dat, 'w') as f:
        print('N {0:d}'.format(poly['n']), file=f)


def read_resume(job):
    """The progress in NAME.job.resume (the relation files already checked
    are in NAME.done), a dict with
      RELS      relations in NAME.dat (missing: unknown)
      GPU_BASE  relations in NAME.dat when NAME.list was written
      GPU_SIZE  size of NAME.dat then, in bytes
      NEXT      relations needed before filtering is tried again, after it
                failed."""
    res = {}
    if os.path.exists(job.resume):
        with open(job.resume) as f:
            for line in f:
                m = re.match(r'(RELS|GPU_BASE|GPU_SIZE|NEXT):\s*(\d+)', line)
                if m:
                    res[m.group(1)] = int(m.group(2))
    return res


def write_resume(job, res):
    tmp = job.resume + '.tmp'
    with open(tmp, 'w') as f:
        for key, val in res.items():
            print('{0:s}: {1:d}'.format(key, val), file=f)
    os.replace(tmp, job.resume)


def minrels_file(job):
    """minrels.txt of the job folder if it exists, else NAME.minrels."""
    return job.minrels_all if os.path.isfile(job.minrels_all) else job.minrels


def target_relations(job):
    """The number in minrels.txt or NAME.minrels, read every round, so it
    can be edited while the job runs."""
    fn = minrels_file(job)
    with open(fn) as f:
        m = re.search(r'\d+', f.read())
    if not m:
        die(job, 'Error: no number in {0:s}'.format(fn))
    return int(m.group())


def failed_filtering(job, res, target):
    """Filtering failed: try again with MORE_RELS of the target more."""
    res['NEXT'] = res['RELS'] + int(MORE_RELS * target) + 1
    write_resume(job, res)
    output(job, '-> Filtering failed; trying again at {0:d} relations'.format(res['NEXT']))


def filter_relations(job, poly):
    """msieve filtering; True if it succeeded (NAME.dat.cyc written)."""
    # -v: also print msieve's progress (duplicates, singletons, merging) here
    if run_msieve(job, ['-v', '-nc1'] + filter_options(poly)):
        die(job, 'msieve filtering failed')
    return os.path.exists(job.cyc)


def batch_check(job, poly, dirs_changed):
    """Step 2: collect relations until msieve filtering succeeds."""
    res = read_resume(job)
    if 'RELS' not in res:
        output(job, '-> Counting the relations in {0:s}'.format(job.dat))
        res['RELS'] = count_lines(job.dat) - 1
    if os.path.exists(job.list) and 'GPU_BASE' not in res:
        # NAME.job.resume was lost; RELS was counted from NAME.dat above.
        drop_list(job, res, read_keys(job.gpu_resume))
    if os.path.exists(job.list) and (CFG['mode'] == 'cpu' or dirs_changed):
        leave_gpu(job, res)
    if CFG['mode'] == 'gpu':
        gpu_batch_check(job, poly, res)
    else:
        cpu_batch_check(job, poly, res)
    for f in ([job.resume, job.done, job.list, job.gpu_resume, job.gpu_status]
              + glob.glob(job.cpu_out + '*')):
        delete_file(f)


def all_checked(job):
    die(job, 'All relation files have been checked, but filtering has not succeeded yet. '
        'Add relation files with -dir, try filtering by hand, or lower the target in '
        '{0:s}.'.format(minrels_file(job)))


def run_factory(job, chunk):
    """Run NFS_factory on the files in chunk at the same time. Returns the
    output file names, or exits if a process fails."""
    outs = ['{0:s}{1:d}'.format(job.cpu_out, j) for j in range(len(chunk))]
    for out in outs:
        delete_file(out)
    procs = []
    try:
        for (_, path), out in zip(chunk, outs):
            procs.append(start_exe(job, NFS_FACTORY, [job.poly, path, out]))
        rets = [p.wait() for p in procs]
    finally:
        for p in procs:
            if p.poll() is None:
                p.kill()
    for (k, _), out, ret in zip(chunk, outs, rets):
        if ret and ret not in IGNORED_EXIT_CODES:
            die(job, 'NFS_factory failed on {0:s} (exit code {1:d})'.format(k, ret))
        if not os.path.exists(out):
            die(job, 'NFS_factory wrote no output for {0:s}'.format(k))
    return outs


def cpu_batch_check(job, poly, res):
    """Step 2 on the CPU: NFS_factory processes, one relation file each."""
    files = relation_files(gpu=False)
    if not files:
        die(job, 'No relation files (.bin or .out) in {0:s}'.format(', '.join(CFG['dir'])))
    rels = res['RELS']
    done = read_done(job)
    todo = [f for f in files if f[0] not in done]
    output(job, '-> {0:d} relation files, {1:d} left to check; {2:d} relations so far'
           .format(len(files), len(todo), rels))

    while not os.path.exists(job.cyc):
        if not todo:
            all_checked(job)
        chunk, todo = todo[:CFG['threads']], todo[CFG['threads']:]
        output(job, '-> Batch checking input files {0:s}'
               .format(', '.join(k for k, _ in chunk)))
        start_time = time.time()
        outs = run_factory(job, chunk)

        # Add the relations to NAME.dat before recording the progress: if this is
        # interrupted, the round is checked again and msieve removes duplicates.
        for out in outs:
            rels += count_lines(out)
            append_file(out, job.dat)
            delete_file(out)
        add_done(job, [k for k, _ in chunk])
        res['RELS'] = rels
        write_resume(job, res)
        output(job, 'LatSieveTime: {0:d}'.format(int(time.time() - start_time)), console=False)

        target = target_relations(job)
        output(job, 'Found {0:d} relations, {1:2.1f}% of the target ({2:d}).'
               .format(rels, 100.0 * rels / target, target))
        if rels >= max(target, res.get('NEXT', 0)) and not filter_relations(job, poly):
            failed_filtering(job, res, target)


def run_factory_gpu(job, rels):
    """Run NFS_factory_GPU on NAME.list, appending to NAME.dat, until it has
    found rels relations from the list (also counting its earlier runs on
    it). Returns its status."""
    delete_file(job.gpu_status)
    p = start_exe(job, NFS_FACTORY_GPU, [job.poly, '-list', job.list, job.dat,
                                         '-rels', str(rels), '-threads', str(CFG['threads']),
                                         '-gpu', str(CFG['gpunum'])])
    try:
        ret = p.wait()
    finally:
        if p.poll() is None:
            p.kill()
    status = read_keys(job.gpu_status)
    # -1 (4294967295 on Windows): unreadable relation file, bad poly or no GPU
    if (ret and ret not in IGNORED_EXIT_CODES) or not status or status.get('error') == '1':
        die(job, 'NFS_factory_GPU failed (exit code {0:d}{1:s}); running the same command '
            'again resumes it'.format(ret, ', see {0:s}'.format(job.gpu_status) if status else ''))
    return status


def sync_gpu_output(job, res):
    """Bring NAME.dat and RELS in line with the resume file of NFS_factory_GPU:
    the relations up to the output size saved there count. Anything after it
    (the batch of an interrupted run) is cut off, as a resumed run of
    NFS_factory_GPU does itself; before its first saved batch, that is all it
    wrote. Returns its resume keys."""
    fres = read_keys(job.gpu_resume)
    if 'output' in fres:
        end, found = int(fres['output']), int(fres.get('found', 0))
    else:
        end, found = res['GPU_SIZE'], 0
    size = os.path.getsize(job.dat)
    if size < end:
        die(job, '{0:s} is shorter than when the GPU checking started or than {1:s} says. '
            'Delete {2:s}, {1:s} and {3:s} and run again to count the relations in {0:s} '
            'and check all relation files again.'.format(job.dat, job.gpu_resume, job.list,
                                                         job.resume))
    if size > end:
        output(job, '-> Removing {0:d} bytes of an interrupted batch from the end of {1:s}'
               .format(size - end, job.dat))
        with open(job.dat, 'r+b') as f:
            f.truncate(end)
    res['RELS'] = res['GPU_BASE'] + found
    write_resume(job, res)
    return fres


def drop_list(job, res, fres):
    """Stop using NAME.list (and NFS_factory_GPU's files for it): the files of
    the list before the one NFS_factory_GPU was on (fres: its resume keys)
    are added to NAME.done. The CPU, or a new list, continues from that file,
    which is checked again from its start."""
    with open(job.list) as f:
        listed = [line.strip() for line in f if line.strip()]
    if fres.get('complete') == '1':
        n = len(listed)
    else:
        n = int(fres.get('file', 0))
    add_done(job, [stem(p) for p in listed[:n]])
    res.pop('GPU_BASE', None)
    res.pop('GPU_SIZE', None)
    write_resume(job, res)
    for f in (job.list, job.gpu_resume, job.gpu_status):
        delete_file(f)


def leave_gpu(job, res):
    """Keep the relations NFS_factory_GPU has added to NAME.dat and continue
    from the relation file it was on: on the CPU, or with a new list."""
    drop_list(job, res, sync_gpu_output(job, res))
    output(job, '-> {0:d} relation files checked so far; continuing with the others'
           .format(len(read_done(job))))


def gpu_batch_check(job, poly, res):
    """Step 2 on the GPU: one NFS_factory_GPU process on a list of the
    relation files, resumed with a higher target while filtering fails."""
    if not os.path.exists(job.list):
        done = read_done(job)
        files = [f for f in relation_files(gpu=True) if f[0] not in done]
        if not files:
            if done:
                all_checked(job)
            die(job, 'No relation files (.ab, .bin or .out) in {0:s}'
                .format(', '.join(CFG['dir'])))
        # NFS_factory_GPU appends to NAME.dat; a resume file of an earlier list
        # (with the same name) would make it cut NAME.dat to that size.
        delete_file(job.gpu_resume)
        delete_file(job.gpu_status)
        res.update(GPU_BASE=res['RELS'], GPU_SIZE=os.path.getsize(job.dat))
        write_resume(job, res)
        # NFS_factory_GPU knows the list by its name, so it is not changed afterwards.
        with open(job.list, 'w') as f:
            for _, path in files:
                print(os.path.abspath(path), file=f)
        output(job, '-> {0:d} relation files (from {1:s}, {2:d} of them .ab) listed in {3:s}'
               .format(len(files), files[0][0], sum(p.lower().endswith('.ab') for _, p in files),
                       job.list))
    sync_gpu_output(job, res)
    output(job, '-> {0:d} relations so far'.format(res['RELS']))

    while not os.path.exists(job.cyc):
        target = target_relations(job)
        want = max(target, res.get('NEXT', 0))
        status = {}
        if res['RELS'] < want:
            output(job, '-> Batch checking on the GPU until {0:d} relations'.format(want))
            start_time = time.time()
            status = run_factory_gpu(job, want - res['GPU_BASE'])
            sync_gpu_output(job, res)
            output(job, 'LatSieveTime: {0:d}'.format(int(time.time() - start_time)),
                   console=False)
            output(job, 'Found {0:d} relations, {1:2.1f}% of the target ({2:d}); '
                   '{3:s} of {4:s} relation files done.'
                   .format(res['RELS'], 100.0 * res['RELS'] / target, target,
                           status.get('files_done', '?'), status.get('files', '?')))
        if res['RELS'] >= target:
            if filter_relations(job, poly):
                break
            failed_filtering(job, res, target)
        if status.get('complete') == '1':
            all_checked(job)


def linear_algebra(job):
    """Step 3."""
    gpu = CFG['mode'] == 'gpu'
    if os.path.exists(job.chk):
        output(job, '-> Resuming the linear algebra')
        step = '-ncr'
    else:
        output(job, '-> Running the linear algebra')
        step = '-nc2'
    ret = run_msieve(job, [step] + (GPU_LA_OPTIONS.split() if gpu else []), gpu=gpu)
    if ret:
        die(job, 'msieve linear algebra failed (exit code {0:d})'.format(ret))
    if not os.path.exists(job.dep):
        die(job, 'msieve did not write the dependencies ({0:s})'.format(job.dep))


def summary_name(job, poly):
    d, f = os.path.split(job.name)
    return os.path.join(d, 's{0:d}-{1:s}.txt'.format(poly['difficulty'], f))


def read_log(job, n):
    """Factors and statistics from NAME.log."""
    st = {'factors': set(), 'sieve': 0.0, 'relproc': 0.0, 'matrix': 0.0, 'sqrt': 0.0,
          'version': 'Msieve', 'rels': '', 'raw_rels': 0, 'matrix_size': ''}
    times = {'LatSieveTime': 'sieve', 'RelProcTime': 'relproc',
             'BLanczosTime': 'matrix', 'sqrtTime': 'sqrt'}
    with open(job.log, errors='replace') as f:
        for line in f:
            m = re.search(r'(LatSieveTime|RelProcTime|BLanczosTime|sqrtTime):\s*(\d+)', line)
            if m:
                st[times[m.group(1)]] += int(m.group(2)) / 3600.0
            m = re.search(r'(Msieve v\..*?)\s*$', line)
            if m:
                st['version'] = m.group(1)
            m = re.search(r'(\d+) unique relations', line)
            if m:
                st['rels'] = m.group(1) + ' relations'
            m = re.search(r'Found (\d+) relations', line)
            if m:
                st['raw_rels'] = int(m.group(1))
            m = re.search(r'matrix is (\d+) x (\d+)', line)
            if m:
                st['matrix_size'] = m.group(1) + ' x ' + m.group(2)
            # msieve 1.52 writes 'prp72 factor:', 1.54 writes 'p72 factor:'
            m = re.search(r'\b(?:prp|p)\d+ factor: (\d+)', line)
            if m and 1 < int(m.group(1)) < n and n % int(m.group(1)) == 0:
                st['factors'].add(int(m.group(1)))
    return st


def write_summary(job, poly, st):
    n, p = poly['n'], poly['params']
    total = st['sieve'] + st['relproc'] + st['matrix'] + st['sqrt']
    name = summary_name(job, poly)
    with open(name, 'w') as f:
        w = lambda s: print(s, file=f)
        w('Number: {0:s}'.format(job.name))
        w('N = {0:d} ({1:d} digits)'.format(n, len(str(n))))
        w('SNFS difficulty: {0:d} digits.'.format(poly['difficulty']))
        w('Divisors found:')
        for r, d in enumerate(sorted(st['factors']), 1):
            w('r{0:d}={1:d} (pp{2:d})'.format(r, d, len(str(d))))
        w('Version: {0:s}'.format(st['version']))
        w('Total time: {0:1.2f} hours.'.format(total))
        w('Factorization parameters were as follows:')
        with open(job.poly) as pf:
            for line in pf:
                w(line.rstrip('\r\n'))
        w('Batch checking on the {0:s}, {1:d} threads'.format(CFG['mode'].upper(), CFG['threads']))
        w('Factor base limits: {0:s}/{1:s}'.format(p['rlim'], p['alim']))
        w('Large primes per side: {0:d}'.format(LARGEP))
        w('Large prime bits: {0:s}/{1:s}'.format(p['lpbr'], p['lpba']))
        w('Total raw relations: {0:d}'.format(st['raw_rels']))
        w('Relations: {0:s}'.format(st['rels']))
        w('Total batch smoothness checking time: {0:1.2f} hours.'.format(st['sieve']))
        w('Total relation processing time: {0:1.2f} hours.'.format(st['relproc']))
        w('Pruned matrix : {0:s}'.format(st['matrix_size']))
        w('Matrix solve time: {0:1.2f} hours.'.format(st['matrix']))
        w('time per square root: {0:1.2f} hours.'.format(st['sqrt']))
        w('total time: {0:1.2f} hours.'.format(total))
        w(platform.processor())
        w(platform.platform())
        w('processors: {0:d}, speed: {1:.2f}GHz'.format(multiprocessing.cpu_count(), proc_speed()))
    output(job, '-> Factorization summary written to {0:s}'.format(name))


def square_root(job, poly):
    """Step 4."""
    output(job, '-> Running the square root')
    ret = run_msieve(job, ['-nc3'])
    if ret:
        die(job, 'msieve square root failed (exit code {0:d})'.format(ret))
    st = read_log(job, poly['n'])
    if not st['factors']:
        die(job, 'The square root found no factors; see {0:s}'.format(job.log))
    for d in sorted(st['factors']):
        output(job, 'factor: {0:d}'.format(d))
    write_summary(job, poly, st)


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def main():
    global CFG
    name, opts = parse_args(sys.argv[1:])
    job = Job(re.sub(r'\.(poly|fb|n|cfg)$', '', name))
    if not os.path.exists(job.poly):
        print('-> Error: {0:s} does not exist'.format(job.poly))
        sys.exit(1)

    output(job, '-> {0:s}: SNFS factory run of {1:s}'.format(SCRIPT_NAME, job.name))
    CFG, changed = load_settings(job, opts)
    poly = read_poly(job)

    gpu = CFG['mode'] == 'gpu'
    exes = [NFS_FACTORY_GPU, MSIEVE, MSIEVE_GPU] if gpu else [NFS_FACTORY, MSIEVE]
    check_programs(job, exes)
    for d in CFG['dir']:
        if not os.path.isdir(d):
            die(job, 'Error: the relation folder {0:s} does not exist'.format(d))

    if os.path.isfile(job.minrels_all):
        output(job, '-> {0:s} takes precedence over {1:s}: target {2:d}'
               .format(job.minrels_all or MINRELS_ALL, job.minrels, target_relations(job)))
    if 'minrels' in opts or not os.path.exists(job.minrels):
        target = opts.get('minrels') or default_minrels(poly)
        with open(job.minrels, 'w') as f:
            print(target, file=f)
        output(job, '-> Relation target {0:d} written to {1:s}{2:s}'.format(
            target, job.minrels, '' if 'minrels' in opts else
            ' (default for {0:s} {1:d}, the pre-sieved side)'.format(*presieved_lpb(poly))))

    output(job, '-> SNFS difficulty {0:d}, degree {1:d}, m = {2:d}; batch checking and linear '
           'algebra on the {3:s}; relation files in {4:s}'
           .format(poly['difficulty'], poly['degree'], poly['m'], CFG['mode'].upper(),
                   ', '.join(CFG['dir'])))

    if not os.path.exists(job.dat):
        setup(job, poly)
    if not os.path.exists(job.cyc):
        batch_check(job, poly, 'dir' in changed)
    if not os.path.exists(job.dep):
        linear_algebra(job)
    else:
        output(job, '-> {0:s} exists, skipping the linear algebra'.format(job.dep))
    if not os.path.exists(summary_name(job, poly)):
        square_root(job, poly)


if __name__ == '__main__':
    main()
