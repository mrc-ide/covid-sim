#!/usr/bin/env python3
r"""Execute a command locally.

This runs a command locally - used by run_many.py CovidSim multiple runner.

    run_local.py --input <input-dir> --output <output-dir> \
            [--lock <lock-dir>] [-j <count]
            <command...>

Options:

    --input <dir> Directory containing input - will be copied elsewhere
    --output <dir> Directory where output ending up in output will be copied
    --lock <dir> Directory to use as the lock dir flag (default
    script_dir/local_lock)
    -j <count> Number of parallel runs to allow at once (default 1)
"""

import argparse
import atexit
import os
import shutil
import subprocess
import sys
import tempfile
import time

def copy_dir(src_dir, dest_dir):
    """Copy src_dir to dest_dir, even if dest_dir exists already."""

    if sys.version_info[0] == 3 and sys.version_info[1] < 8:
        # We have to do something slightly different for versions of Python
        # prior to 3.8.0.
        import distutils.dir_util
        distutils.dir_util.copy_tree(src_dir, dest_dir)
    else:
        shutil.copytree(src_dir, dest_dir, dirs_exist_ok=True)

def parse_args():
    """Parse the arguments.

    On exit: Returns the result of calling argparse.parse()

    args.input is the input directory
    args.output is the directory where output will be stored
    args.j is number of jobs to run at once
    args.lock is lock directory
    args.cmd is the command to run
    """
    parser = argparse.ArgumentParser()
    script_path = os.path.dirname(os.path.realpath(__file__))

    parser.add_argument(
            "--input",
            help="Location of input data",
            required=True)
    parser.add_argument(
            "--output",
            help="Directory to store output data",
            required=True)
    parser.add_argument(
            "-j",
            help="Number of jobs to run at once",
            type=int,
            default=1)
    parser.add_argument(
            "--lock",
            help="Lock directory root",
            default=os.path.join(script_path, "local_lock"))
    parser.add_argument(
            "--keep",
            help="Keep temp directory",
            action="store_true")
    parser.add_argument(
            "cmd",
            help="Command line to run",
            nargs="+")
    args = parser.parse_args()

    return args

args = parse_args()

# Get the lock, need to make lock directory an absolute path as makedirs gets
# confused if it encounters .. or .
args.lock = os.path.abspath(args.lock)
os.makedirs(args.lock, exist_ok=True)
lock_dir = ""
while lock_dir == "":
    for i in range(args.j):
        try:
            lock_dir = os.path.join(args.lock, "job-{0}".format(i))
            os.mkdir(lock_dir)
            break
        except FileExistsError:
            lock_dir = ""

    if lock_dir == "":
        time.sleep(1)

# Delete the lock directory on exit
atexit.register(shutil.rmtree, path=lock_dir, ignore_errors=True)

# Get the temporary directory, and make sure we destroy it after running
run_dir = tempfile.mkdtemp()

# Copy input into run directory
copy_dir(args.input, run_dir)
if not args.keep:
    atexit.register(shutil.rmtree, path=run_dir, ignore_errors=True)
else:
    print("Data files at: {0} for ".format(run_dir) + " ".join(args.cmd))
os.makedirs(os.path.join(run_dir, "output"), exist_ok=True)

# Change to directory and run command
cwd = os.getcwd()
os.chdir(run_dir)
with open(os.path.join(os.curdir, "output", "stdout.txt"), "w") as outf, \
     open(os.path.join(os.curdir, "output", "stderr.txt"), "w") as errf:
    result = subprocess.run(
            args.cmd,
            stdout=outf,
            stderr=errf)
if result.returncode != 0:
    print("Command failed: " + " ".join(args.cmd))

# Copy output tree
os.chdir(cwd)
copy_dir(os.path.join(run_dir, "output"), args.output)

if result.returncode != 0:
    exit(1)
else:
    exit(0)
