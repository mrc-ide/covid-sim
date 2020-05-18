#!/usr/bin/env python3
r"""Integration test infrastructure.

Invoke as:

integration-test.py --covidsim <exe> --input <input_dir> \
   --output <output_dir> --popfile <pop file> \
   [--schools <shoolfile>] [--accept]

It will run 3 tests:
a. No intervention, generating network from scratch
b. No intervention, loading network generated in a
c. Intervention, using network generated in a.

It will check outputs from a & b are identical, and that the outputs
match expected checksums.
"""

# WHAT TO DO IF THE TEST FAILS
#
# 1. Clone two copies of the git repo (one called "good" and
#    one called "bad"):
#
#    git clone https://github.com/mrc-ide/covid-sim.git good
#    git clone https://github.com/mrc-ide/covid-sim.git bad
#    cd bad
#    git checkout mybranch
#    cd ..
#    cd good
#    get checkout master
#    cd ..
#
# 2. Run the test in the good and bad clones:
#
#    cd bad/test
#    ./integration-test.py <opts>
#    cd ../..
#    cd good/test
#    ./integration-test.py <opts>
#    cd ../..
#
# 3. Use a file comparison tool, like meld (https://meldmerge.org/) to
#    compare the two directories:
#
#    meld good bad
#
#    You need to look at the changes in the generated .xls files
#    and decide whether the changes are ok.
#
# 4. If the changes are ok, then you can fix the test by updating the
#    checksums file:
#
#    cd bad/test
#    ./integration-test.py <opts> --accept
#    git add *-input/results.checksums.txt
#    git commit -m "Update checksums"
#    git push

# EXPLANATION
#
# The test runs a simulation and compares the generated .xls files to
# make sure that the output has not changed. For the test to pass, the
# output must be byte-for-byte identical to the expected result.
# Since the generated .xls files contain floating point numbers, it is
# possible that a legitimate change to the algorithm could cause the
# test to fail due to very minor numerical differences, like
# 0.0308373797 becoming 0.0308373798. In its current form, the test
# does not have any concept of a "tolerance range", so even such a
# minor numerical difference will cause the test to fail. However,
# this does not mean that test will sometimes fail "randomly".
# Floating point computations are deterministic, so only a change to
# the code can cause the test to fail.
#
# The .xls files are large, so they are not checked into the git
# repo. Instead, the test uses sha512sum to compute a checksum for
# each of the generated files. Those checksums are compared against
# the expected checksums in regressiontest_UK_100th.checksums.txt.
# The disadvantage of using checksums is that they are not helpful
# when you are trying to figure out why the test failed. That's why
# the instructions above recommend cloning two copies of the repo, so
# that you can compare the failing versions of the .xls files against
# good versions.

import argparse
import glob
import gzip
import hashlib
import os
import sys
import shutil
import subprocess


def parse_args():
    """Parse the arguments.

    On exit: Returns the result of calling argparse.parse()
    args.covidsim = name of the CovidSim executable
    args.input = directory with input data
    args.output = directory to put output data in
    args.popfile = tar.gz population file
    args.schools = school file
    args.accept = accept new results
    """
    # Defaults
    script_path = os.path.dirname(os.path.realpath(__file__))
    output_dir = script_path
    src_dir = os.path.join(script_path, os.pardir)

    parser = argparse.ArgumentParser()
    parser.add_argument(
            "--covidsim",
            help="Path to CovidSim binary, if none specified will build")
    parser.add_argument(
            "--input",
            help="Input directory",
            required=True)
    parser.add_argument(
            "--popfile",
            help="Population tar.gz file to use",
            required=True)
    parser.add_argument(
            "--schools",
            help="Schools file to use")
    parser.add_argument(
            "--output",
            help="Output directory, default: current dir",
            default=output_dir)
    parser.add_argument(
            "--srcdir",
            help="Source directory, default = script_dir/../src",
            default=src_dir)
    parser.add_argument(
            "--accept",
            help="Accept the results",
            action="store_true")
    parser.add_argument(
            "--r",
            help="r value to pass to covid-sim, default = 1.5",
            default="1.5")

    return parser.parse_args()


failed = False
args = parse_args()

# Sort out directories
shutil.rmtree(args.output, ignore_errors=True)
os.makedirs(args.output, exist_ok=False)

if args.covidsim is None:
    build_dir = os.path.join(args.output, "build")
    src_dir = os.path.realpath(args.srcdir)

    print("Building CovidSim:\n  Source = {0}\n  Output = {1}".format(
        args.srcdir,
        build_dir))
    os.makedirs(build_dir, exist_ok=False)
    cwd = os.getcwd()
    os.chdir(build_dir)
    subprocess.run(['cmake', args.srcdir])
    subprocess.run(['cmake', '--build', '.'])

    if os.name == 'nt':
        covidsim_exe = os.path.join(build_dir, "Debug", "src", "CovidSim.exe")
    else:
        covidsim_exe = os.path.join(build_dir, "src", "CovidSim")

    os.chdir(cwd)
else:
    covidsim_exe = args.covidsim

print("Using executable: {0}".format(covidsim_exe))

# Check files exist
if not os.path.exists(args.popfile):
    print("Unable to find population file: {0}".format(args.popfile))
    exit(1)
if args.schools is not None and not os.path.exists(args.schools):
    print("Unable to find schools file: {0}".format(args.schools))
    exit(1)
for i in ["pre-params.txt", "input-noint-params.txt",
          "admin-params.txt", "input-params.txt",
          "results.checksums.txt"]:
    f = os.path.join(args.input, i)
    if not os.path.exists(f):
        print("Unable to find input file: {0}".format(f))
        exit(1)

# Population density file in gziped form, text file, and binary file as
# processed by CovidSim
wpop_file = os.path.join(args.output, "pop.txt")

# gunzip wpop fie
with gzip.open(args.popfile, 'rb') as f_in,  open(wpop_file, 'wb') as f_out:
    shutil.copyfileobj(f_in, f_out)

# Run the simulation.
mt_loop = {}
mt_loop[1] = 'st-'
mt_loop[2] = 'mt-'

for threads, prefix in mt_loop.items():

    wpop_bin = os.path.join(args.output, prefix + "pop.bin")
    network_bin = os.path.join(args.output, prefix + "network.bin")
    print('=== RUN 1: No Intervention - Build network {0}:'.format(prefix))
    cmd = [
            covidsim_exe,
            '/c:{0}'.format(threads),
            '/BM:bmp',
            '/PP:' + os.path.join(args.input, 'pre-params.txt'),
            '/P:' + os.path.join(args.input, 'input-noint-params.txt'),
            '/O:' + os.path.join(args.output, prefix + 'results-noint'),
            '/D:' + wpop_file,
            '/M:' + wpop_bin,
            '/A:' + os.path.join(args.input, 'admin-params.txt')
    ]
    if args.schools:
        cmd.extend(["/s:" + args.schools])
    cmd.extend([
            '/S:' + network_bin,
            '/R:' + args.r,
            '98798150',
            '729101',
            '17389101',
            '4797132'
    ])
    print("Command line: " + " ".join(cmd))
    process = subprocess.run(cmd, check=True)

    print('=== RUN 2: No Intervention - Load network:')
    cmd = [
            covidsim_exe,
            '/c:{0}'.format(threads),
            '/BM:bmp',
            '/PP:' + os.path.join(args.input, 'pre-params.txt'),
            '/P:' + os.path.join(args.input, 'input-noint-params.txt'),
            '/O:' + os.path.join(args.output, prefix + 'results-noint-repeat'),
            '/D:' + wpop_bin,
            '/A:' + os.path.join(args.input, 'admin-params.txt')
    ]
    if args.schools:
        cmd.extend(["/s:" + args.schools])
    cmd.extend([
            '/L:' + network_bin,
            '/R:' + args.r,
            '98798150',
            '729101',
            '17389101',
            '4797132'
    ])
    print("Command line: " + " ".join(cmd))
    process = subprocess.run(cmd, check=True)

    print('=== RUN 3: Intervention - Load network:')
    cmd = [
            covidsim_exe,
            '/c:{0}'.format(threads),
            '/BM:bmp',
            '/PP:' + os.path.join(args.input, 'pre-params.txt'),
            '/P:' + os.path.join(args.input, 'input-params.txt'),
            '/O:' + os.path.join(args.output, prefix + 'results-int'),
            '/D:' + wpop_bin,
            '/A:' + os.path.join(args.input, 'admin-params.txt')
    ]
    if args.schools:
        cmd.extend(["/s:" + args.schools])
    cmd.extend([
            '/L:' + network_bin,
            '/R:' + args.r,
            '98798150',
            '729101',
            '17389101',
            '4797132'
    ])
    print("Command line: " + " ".join(cmd))
    process = subprocess.run(cmd, check=True)

    repeat_files_checked = 0
    for fn2 in glob.glob(os.path.join(args.output, prefix + 'results-noint-repeat*.xls')):
        fn1 = fn2.replace('-repeat', '')
        with open(fn1, 'rb') as f:
            dat1 = f.read()
        with open(fn2, 'rb') as f:
            dat2 = f.read()
        # We expect multiple reasonably large files, so only count those that are
        # at least 10 bytes long
        if len(dat1) > 10:
            repeat_files_checked += 1
        if dat1 != dat2:
            print('FAILURE: Contents of ' + fn1 + ' does not match that of ' + fn2)
            failed = True
    if repeat_files_checked < 3:
        print('FAILURE: Not enough repeat files found')
        failed = True

expected_checksums = os.path.join(args.input, 'results.checksums.txt')
actual_checksums = os.path.join(args.output, 'results.checksums.txt')

# Compute SHA-512 checksums for the generated files.
paths = []
sha512sums = []

# Scan output for files to check
for direntry in os.scandir(args.output):
    name = direntry.name

    if name.endswith('.xls'):
        # TSV files
        paths.append((args.output, name))
    elif name.endswith('.ge'):
        # Bitmap files - only check every 100th.
        for direntry2 in os.scandir(direntry.path):
            name2 = direntry2.name
            if name2.endswith('00.bmp'):
                paths.append((direntry.path, name2))

max_filename_len = max(map(len, [filename for dirname, filename in paths]))

for dirname, filename in paths:
    with open(os.path.join(dirname, filename), 'rb') as f:
        dat = f.read()
        sha = hashlib.sha512(dat).hexdigest()
        line = filename + (' ' * (1 + max_filename_len - len(filename))) + sha
        print(line)
        sha512sums.append(line)

sha512sums.sort()
print('New checksums:')
print('\n'.join(sha512sums))
print('end')

# Write the checksums into a file in the temporary directory. (This is
# useful for updating the reference checksums file if the results have
# changed for a legitimate reason.)
with open(actual_checksums, 'wb') as checksums_outfile:
    checksums_outfile.write('\n'.join(sha512sums).encode('utf-8'))

# Read the expected checksums from the reference file.
sha512sums_reference = list(
        filter(None,
               open(expected_checksums, 'rb').read().decode('utf-8')
               .split('\n')))
print('Reference checksums:')
print('\n'.join(sha512sums_reference))
print('end')

# Compare the checksums against the expected values.
if sha512sums == sha512sums_reference:
    print('SUCCESS: checksums match')
else:
    print('FAILURE: checksums do not match')
    len1 = len(sha512sums)
    len2 = len(sha512sums_reference)
    print('length: ' + str(len1))
    print('reference length: ' + str(len2))
    if len1 != len2:
        print('Lengths don\'t match.')
    else:
        for x, y in zip(sha512sums, sha512sums_reference):
            if x != y:
                print('Mismatch:')
                print(str(x))
                print(str(y))
                break

    if args.accept:
        print('Accepting results.')
        with open(actual_checksums, 'rb') as checksums_in, \
                open(expected_checksums, 'wb') as checksums_out:
                shutil.copyfileobj(checksums_in, checksums_out)
    else:
        failed = True

if failed:
    print('TEST FAILED.')
    sys.exit(1)
