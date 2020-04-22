#!/usr/bin/env python3
# Regression test for checking output of CovidSim
# NOTE: This uses test data and is not a run reflective of real-world
# situations.

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
#    ./regressiontest_UK_100th.py
#    cd ../..
#    cd good/test
#    ./regressiontest_UK_100th.py
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
#    python regressiontest_UK_100th.py --accept
#    git add regressiontest_UK_100th.checksums.txt
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

import gzip
import os
import sys
import shutil
import subprocess

testdir='regressiontest_UK_100th'

accept_results = False
if len(sys.argv) > 1 and sys.argv[1] == '--accept':
    accept_results = True

# Portable ../../
updir2 = os.pardir + os.sep + os.pardir + os.sep
datadir = updir2 + "data" + os.sep
updir1 = os.pardir + os.sep

shutil.rmtree(testdir, True)
os.mkdir(testdir)
os.chdir(testdir)
subprocess.check_call(['cmake', '-DCMAKE_CXX_FLAGS=-DNO_WIN32_BM', updir2 + 'src'])
subprocess.check_call(['cmake', '--build', '.'])

if os.name == 'nt':
    covidsim_exe = os.getcwd() + os.sep + 'Debug\CovidSim.exe'
else:
    covidsim_exe = os.getcwd() + os.sep + 'CovidSim'
print(covidsim_exe)

# Population density file in gziped form, text file, and binary file as
# processed by CovidSim
wpop_file_gz = os.path.join(datadir, "populations", "wpop_eur.txt.gz")
wpop_file = "wpop_eur.txt"
wpop_bin = "wpop_test.bin"

# gunzip wpop fie
if not os.path.exists(wpop_file):
    with gzip.open(wpop_file_gz, 'rb') as f_in:
            with open(wpop_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

# Run the simulation.
print('=== Starting building network:')
subprocess.check_call(
    [covidsim_exe, '/c:1',
     '/PP:' +  updir1 + 'preUK_R0=2.0.txt',
     '/P:' + updir1  + 'p_NoInt.txt', '/CLP1:100000',
     '/CLP2:0', '/O:NoInt_R0=2.2', '/D:' + wpop_file, '/M:' + wpop_bin,
     '/A:' + updir1 + 'sample_admin.txt',
     '/S:NetworkUKN_32T_100th.bin', '/R:1.1',
     '98798150', '729101', '17389101', '4797132'
    ])
print('=== Starting running:')
subprocess.check_call(
    [covidsim_exe, '/c:1',
     '/PP:' + updir1 + 'preUK_R0=2.0.txt',
     '/P:' + updir1 + 'p_PC7_CI_HQ_SD.txt', '/CLP1:100',
     '/CLP2:91', '/CLP3:121', '/CLP4:121', '/O:CI_100_91_R0=2.2',
     '/A:' + updir1 + 'sample_admin.txt',
     '/D:' + wpop_bin,
     '/L:NetworkUKN_32T_100th.bin', '/R:1.1',
     '98798150', '729101', '17389101', '4797132'
    ])
print('=== Done')

checksums_filename='regressiontest_UK_100th.checksums.txt'

# Compute SHA-512 checksums for the generated files.
paths = []
sha512sums = []
for direntry in os.scandir():
    name = direntry.name
    if name.endswith('.xls'):
        paths.append((None, name))
paths.append(('CI_100_91_R0=2.2.ge', 'CI_100_91_R0=2.2.00100.bmp'))
max_filename_len = max(map(len, [ filename for dirname, filename in paths ]))

for dirname, filename in paths:
    try:
        # sha512sum (usually installed on Linux, sometimes on Windows)
        output = subprocess.check_output(['sha512sum', '--binary', filename], cwd=dirname)
        sha = output.decode('utf-8').split(' ')[0]
    except FileNotFoundError:
        # certUtil (always? available in modern Windows)
        output = subprocess.check_output(['certUtil', '-hashfile', filename, 'SHA512'], cwd=dirname)
        sha = output.decode('utf-8').splitlines()[1]
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
with open(checksums_filename, 'wb') as checksums_outfile:
    checksums_outfile.write('\n'.join(sha512sums).encode('utf-8'))

# Read the expected checksums from the reference file.
sha512sums_reference = list(filter(None, open(os.pardir + os.sep + checksums_filename, 'rb').read().decode('utf-8').split('\n')))
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

    if accept_results:
        print('Accepting results.')
        with open(checksums_filename, 'rb') as checksums_in:
            with open(os.pardir + os.sep + checksums_filename, 'wb') as checksums_out:
                shutil.copyfileobj(checksums_in, checksums_out)
    else:
        sys.exit(1)
