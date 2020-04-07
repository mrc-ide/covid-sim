#!/usr/bin/env python3
# Regression test for checking output of SpatialSim in US mode
# NOTE: This uses test data and is not a run reflective of real-world
# situations.

# WHAT TO DO IF THE TEST FAILS
#
# 1. Clone two copies of the git repo (one called "good" and
#    one called "bad"):
#
#    git clone https://github.com/mrc-ide/covid-19-spatial-sim.git good
#    git clone https://github.com/mrc-ide/covid-19-spatial-sim.git bad
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
#    ./regressiontest_US_based.py
#    cd ../..
#    cd good/test
#    ./regressiontest_US_based.py
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
#    ./regressiontest_US_based.py --accept
#    git add us-input/results.checksums.txt
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

testdir='regressiontest_US_based'

accept_results = False
if len(sys.argv) > 1 and sys.argv[1] == '--accept':
    accept_results = True

# Sort out directories
script_dir = os.path.dirname(os.path.realpath(__file__))
input_dir = os.path.join(script_dir, 'us-input')
output_dir = os.path.join(script_dir, 'us-output')
src_dir = os.path.join(script_dir, os.pardir, 'src')
data_dir = os.path.join(script_dir, os.pardir, 'data')

shutil.rmtree(output_dir, ignore_errors=True)
os.makedirs(output_dir, exist_ok=False)
os.chdir(output_dir)
subprocess.run(['cmake', src_dir, "-DCOUNTRY=US"])
subprocess.run(['cmake', '--build', '.'])

if os.name == 'nt':
    spacialsim_exe = os.getcwd() + os.sep + 'Debug\SpatialSim.exe'
else:
    spacialsim_exe = os.getcwd() + os.sep + 'SpatialSim'
print(spacialsim_exe)

# Population density file in gziped form, text file, and binary file as
# processed by SpatialSim
wpop_file_gz = os.path.join(data_dir, "populations", "wpop_usacan.txt.gz")
wpop_file = os.path.join(output_dir, "wpop_usacan.txt")
wpop_bin = os.path.join(output_dir, "wpop-test.bin")
wpop_bin_schools = os.path.join(output_dir, "wpop-test-schools.bin")

# gunzip wpop fie
if not os.path.exists(wpop_file):
    with gzip.open(wpop_file_gz, 'rb') as f_in:
            with open(wpop_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

# Run the simulation.
print('=== Starting building network (no schools):')
subprocess.run(
    [spacialsim_exe, '/c:1',
     '/PP:' +  os.path.join(input_dir, 'pre-params.txt'),
     '/P:' + os.path.join(input_dir, 'input-noint-params.txt'),
     '/O:' + os.path.join(output_dir, 'results-noint'),
     '/D:' + wpop_file,
     '/M:' + wpop_bin,
     '/A:' + os.path.join(input_dir, 'admin-params.txt'),
     '/S:' + os.path.join(output_dir, 'network.bin'),
     '/R:1.5',
     '98798150', '729101', '17389101', '4797132'
    ])
print('=== Starting running (no schools):')
subprocess.run(
    [spacialsim_exe, '/c:1',
     '/PP:' +  os.path.join(input_dir, 'pre-params.txt'),
     '/P:' + os.path.join(input_dir, 'input-params.txt'),
     '/O:' + os.path.join(output_dir, 'results-int'),
     '/A:' + os.path.join(input_dir, 'admin-params.txt'),
     '/D:' + wpop_bin,
     '/L:' + os.path.join(output_dir, 'network.bin'),
     '/R:1.5',
     '98798150', '729101', '17389101', '4797132'
    ])
print('=== Starting building network (schools):')
subprocess.run(
    [spacialsim_exe, '/c:1',
     '/PP:' +  os.path.join(input_dir, 'pre-params.txt'),
     '/P:' + os.path.join(input_dir, 'input-noint-params.txt'),
     '/O:' + os.path.join(output_dir, 'results-noint-schools'),
     '/D:' + wpop_file,
     '/M:' + wpop_bin_schools,
     '/A:' + os.path.join(input_dir, 'admin-params.txt'),
     '/S:' + os.path.join(output_dir, 'network-schools.bin'),
     '/s:' + os.path.join(data_dir, 'populations', 'USschools.txt'),
     '/R:1.5',
     '98798150', '729101', '17389101', '4797132'
    ])
print('=== Starting running (no schools):')
subprocess.run(
    [spacialsim_exe, '/c:1',
     '/PP:' +  os.path.join(input_dir, 'pre-params.txt'),
     '/P:' + os.path.join(input_dir, 'input-params.txt'),
     '/O:' + os.path.join(output_dir, 'results-int-schools'),
     '/A:' + os.path.join(input_dir, 'admin-params.txt'),
     '/D:' + wpop_bin_schools,
     '/L:' + os.path.join(output_dir, 'network-schools.bin'),
     '/s:' + os.path.join(data_dir, 'populations', 'USschools.txt'),
     '/R:1.5',
     '98798150', '729101', '17389101', '4797132'
    ])
print('=== Done')

expected_checksums = os.path.join(input_dir, 'results.checksums.txt')
actual_checksums = os.path.join(output_dir, 'results.checksums.txt')

# Compute SHA-512 checksums for the generated .xls files.
sha512sums = []
for direntry in os.scandir():
    name = direntry.name
    if name.endswith('.xls'):
        try:
            # sha512sum (usually installed on Linux, sometimes on Windows)
            sha512sum = subprocess.check_output(['sha512sum', '--binary', name])
        except FileNotFoundError:
            # certUtil (always? available in modern Windows)
            output = subprocess.check_output(['certUtil', '-hashfile', name, 'SHA512'])
            sha512sum = output.splitlines()[1] + b" *" + name.encode() + b"\n"
        print(sha512sum.decode('utf-8'))
        sha512sums.append(sha512sum)

sha512sums.sort()
print('New checksums:')
print(b''.join(sha512sums).decode('utf-8'))
print('end')

# Write the checksums into a file in the temporary directory. (This is
# useful for updating the reference checksums file if the results have
# changed for a legitimate reason.)
with open(actual_checksums, 'wb') as checksums_outfile:
    for sha512sum in sha512sums:
        checksums_outfile.write(sha512sum)

# Read the expected checksums from the reference file.
sha512sums_reference = open(expected_checksums, 'rb').readlines()
print('Reference checksums:')
print(b''.join(sha512sums_reference).decode('utf-8'))
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
        with open(actual_checksums, 'rb') as checksums_in:
            with open(expected_checksums, 'wb') as checksums_out:
                shutil.copyfileobj(checksums_in, checksums_out)
    else:
        sys.exit(1)
