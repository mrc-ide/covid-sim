#!/usr/bin/env python3
"""Regression test for checking output of CovidSim in US mode."""

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

import os
import sys
import subprocess

# Sort out directories
script_dir = os.path.dirname(os.path.realpath(__file__))
input_dir = os.path.join(script_dir, 'us-input')
output_dir = os.path.join(script_dir, 'us-output')
wpop_file_gz = os.path.join(
        script_dir,
        os.pardir,
        "data",
        "populations",
        "wpop_usacan.txt.gz")

print("Invoking integration-test.py")
cmd = [
        sys.executable,
        os.path.join(script_dir, "integration-test.py"),
        "--output", output_dir,
        "--input", input_dir,
        "--popfile", os.path.join(wpop_file_gz)
]

cmd.extend(sys.argv[1:])
print("Command line: " + " ".join(cmd))
process = subprocess.run(cmd, check=True)
