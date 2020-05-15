#!/usr/bin/env python3
r"""Do multiple runs locally or remotely.

This Python script provides a framework for doing multiple runs of CovidSim
either locally or remotely.

It should be invoked as follows:

    python3 run_many.py [--covidsim <binary>] [--data <dir>] \
                        [--output <dir>] --runner runner_script... \
                        [--srcdir <srcdir>] config.json

For example:
    python3 run_many.py --runner runner/local.py run_us.json

Options:
    --covidsim Path to CovidSim binary to use, if not specified will build
    one using srcdir as the directory of the source (default parent dir of
    script).

    --output Output directory root.  If not specified will use the
    current directory

    --runner Script to run to actually call CovidSim, use multiple times if
    you need to specify more than one parameter to invoke the runner.

    --data Directory containing the admin units and population data (default
    script dir).

Operand

    config.json: Configuration file.

The runner script does the actual invocation of CovidSim.  It will be invoked
with the following command-line:

    <RUNNER> --input <input-dir> --output <output-dir> <command...>

Its execution is the following steps:

    1. Copy all contents of <input-dir> to the (remote) location where
    execution will take place.

    2. Change to that (remote) location.

    3. Execute <comamnd...>

    4. Copy/move the (remote) contents of ./output/ to <output-dir>.

    5. Exit (0 for success, non-zero for failure).

The JSON configuration takes a set of config options which control what gets
invoked:

{
    "threads" : <NUM> # Number of threads to specify for CovidSim to use
    "num_runs" : <NUM> # Number of runs with different run seeds
    "network_seed0" : <NUM> # Network seed 0
    "network_seed1" : <NUM> # Network seed 1
    "r0" : [<r0>] # R0 values to test
    "geographies" : [<geography>] # Geographies to test
    "param_files" : [<Param file>] # Parameter files to test
}

There will be num_runs runs run for each combination of r0, geography, and
param-file.
"""

import argparse
import gzip
import json
import os
import random
import shutil
import subprocess
import time


def try_remove(f):
    """Try and remove file f.

    No error if removing is impossible.
    """
    try:
        os.remove(f)
    except OSError:
        pass


def parse_args():
    """Parse the arguments.

    On exit: Returns the result of calling argparse.parse()

    args.covidsim is the name of the CovidSim executable
    args.data is the directory with the input data
    args.output is the directory where output will be stored
    args.srcdir is the directory for the source code
    args.runner is the list of arguments to invoke runner
    args.config is the config file
    """
    parser = argparse.ArgumentParser()
    script_path = os.path.dirname(os.path.realpath(__file__))

    # Default values
    data_dir = script_path
    output_dir = os.getcwd()
    src_dir = os.path.join(data_dir, os.pardir)

    parser.add_argument(
            "--covidsim",
            help="Location of CovidSim binary, if none specified will build")
    parser.add_argument(
            "--data",
            help="Directory at root of input data",
            default=data_dir)
    parser.add_argument(
            "--srcdir",
            help="Directory with source in - needed if --covidsim " +
            "isn't specified",
            default=src_dir)
    parser.add_argument(
            "--output",
            help="Directory to store output data",
            default=output_dir)
    parser.add_argument(
            "--runner",
            help="Command line for runner",
            action="append")
    parser.add_argument(
            "config",
            help="Configuration file")
    args = parser.parse_args()

    return args


args = parse_args()

# Lists of places that need to be handled specially
united_states = ["United_States"]
canada = ["Canada"]
usa_territories = ["Alaska", "Hawaii", "Guam", "Virgin_Islands_US",
                   "Puerto_Rico", "American_Samoa"]
nigeria = ["Nigeria"]


def admin_base(geography):
    """Get the admin file basename for a particular geography."""
    return "{0}_admin.txt".format(geography)


def pp_base(geography):
    """Get the basename of the pre-param file for a particular geography."""
    if geography in united_states:
        return "preUS_R0=2.0.txt"
    elif geography in nigeria:
        return "preNGA_R0=2.0.txt"
    else:
        return "preUK_R0=2.0.txt"


# Determine whether we need to build the tool or use a user supplied one:
if args.covidsim is not None:
    exe = args.covidsim
else:
    build_dir = os.path.join(args.output, "build")

    # Ensure we do a clean build
    shutil.rmtree(build_dir, ignore_errors=True)
    os.makedirs(build_dir, exist_ok=False)
    cwd = os.getcwd()
    os.chdir(build_dir)

    # Build
    subprocess.run(['cmake', args.srcdir], check=True)
    subprocess.run(['cmake', '--build', '.'], check=True)

    # Where the exe ends up depends on the OS.
    if os.name == 'nt':
        exe = os.path.join(build_dir, "Debug", "src", "CovidSim.exe")
    else:
        exe = os.path.join(build_dir, "src", "CovidSim")

    os.chdir(cwd)

# Ensure output directory exists
os.makedirs(args.output, exist_ok=True)

# Read the config file
if not os.path.exists(args.config):
    print("Unable to find config file: {0}".format(args.config))
    exit(1)
with open(args.config, 'r') as cf:
    config = json.load(cf)

# If we have no threads variable set use 1:
if config["threads"] is None:
    config["threads"] = 1

if config["network_seed1"] is None:
    print("Need network_seed1 in config")
    exit(1)

if config["network_seed2"] is None:
    print("Need network_seed2 in config")
    exit(1)

if config["num_runs"] is None:
    config["num_runs"] = 1

# Run Seed reference gives the list of run seeds to use for each set of runs.
# We need 2 seeds for each run - so it will be config["num_runs"] * 2 elements
# long.
run_seed_reference = []
random.seed()
for i in range(config["num_runs"] * 2):
    run_seed_reference.append(random.randint(0, 0x7fffffff))

# Generate the input data directories.  We have the following layout:
#  - input-setup/
#  - input/<geography>/
#
# input-setup contains the data needed to run the initial network generation
# run for all geographies.  input-<geography> contains the files needed for
# the runs in a particular geography (including the binary pop data and
# network file).

input_setup = os.path.join(args.output, "input-setup")
os.makedirs(input_setup, exist_ok=True)
os.makedirs(os.path.join(input_setup, "output"), exist_ok=True)

# Copy the executable - use copy2 to preserve exe bit
shutil.copy2(exe, os.path.join(input_setup, "CovidSim.exe"))

# Dictionaries of jobs:
setup_jobs = {}  # Geography -> Setup run, emptied as jobs finish
all_jobs = {}  # ID -> Intervention run

# Run through the setup runs:
for geography in config["geographies"]:
    input_geography = os.path.join(args.output, "input", geography)
    output_geography = os.path.join(args.output, "output", geography, "setup")
    os.makedirs(input_geography, exist_ok=True)
    os.makedirs(output_geography, exist_ok=True)
    os.makedirs(os.path.join(input_geography, "output"), exist_ok=True)

    # Copy the executable
    shutil.copy2(exe, os.path.join(input_geography, "CovidSim.exe"))

    # The admin file to use
    admin_file = os.path.join(args.data, "admin_units", admin_base(geography))

    if not os.path.exists(admin_file):
        print("Unable to find admin file for geography: {0}".format(geography))
        print("Data directory: {0}".format(args.data))
        print("Looked for: {0}".format(admin_file))
        exit(1)

    shutil.copyfile(admin_file,
                    os.path.join(input_setup, admin_base(geography)))
    shutil.copyfile(admin_file,
                    os.path.join(input_geography, admin_base(geography)))

    # Population density file in gziped form, text file, and binary file as
    # processed by CovidSim
    if geography in united_states + canada:
        wpop_file_root = "usacan"
    elif geography in usa_territories:
        wpop_file_root = "us_terr"
    elif geography in nigeria:
        wpop_file_root = "nga_adm1"
    else:
        wpop_file_root = "eur"

    wpop_file_gz = os.path.join(
            args.data,
            "populations",
            "wpop_{0}.txt.gz".format(wpop_file_root))
    wpop_file = os.path.join(
            input_setup,
            "wpop_{0}.txt".format(wpop_file_root))

    if not os.path.exists(wpop_file):
        if not os.path.exists(wpop_file_gz):
            print("Unable to find population file for geography: {0}".
                  format(geography))
            print("Data directory: {0}".format(args.data))
            print("Looked for: {0}".format(wpop_file_gz))
            exit(1)

        # Uncompress
        with gzip.open(wpop_file_gz, 'rb') as fi, open(wpop_file, 'wb') as fo:
            shutil.copyfileobj(fi, fo)

    # We need to run CovidSim to generate the binary - and we don't have that
    # yet.

    # Configure pre-parameter file.  This file doesn't change between runs:

    pp_file = os.path.join(args.data, "param_files", pp_base(geography))
    if not os.path.exists(os.path.join(input_setup, pp_base(geography))):
        if not os.path.exists(pp_file):
            print("Unable to find pre-parameter file")
            print("Data directory: {0}".format(args.data))
            print("Looked for: {0}".format(pp_file))
            exit(1)

        shutil.copyfile(pp_file, os.path.join(input_setup, pp_base(geography)))

    # Copy pre-parameter file into Geography
    shutil.copyfile(pp_file, os.path.join(input_geography, pp_base(geography)))

    # Copy the parameter files
    # We copy all of them into each Geo input directory, but only the first
    # into the input-setup directory.
    for param_file in config["param_files"]:
        if not os.path.exists(param_file):
            print("Unable to find parameter file: {0}".format(param_file))
            exit(1)

        param_base = os.path.basename(param_file)
        if os.path.exists(os.path.join(input_geography, param_base)):
            print("Parameter files with the same basename: {0}".
                  format(param_base))
            exit(1)
        shutil.copyfile(param_file, os.path.join(input_geography, param_base))

    param_base = os.path.basename(config["param_files"][0])
    shutil.copyfile(config["param_files"][0], os.path.join(input_setup,
                    param_base))

    covidsim_cmd = [
        os.path.join(os.path.curdir, "CovidSim.exe"),
        "/c:{0}".format(config["threads"]),
        "/A:" + os.path.join(os.path.curdir, admin_base(geography)),
        "/PP:" + os.path.join(os.path.curdir, pp_base(geography)),
        "/P:" + os.path.join(os.path.curdir, param_base),
        "/O:" + os.path.join(os.path.curdir, "output", "setup"),
        "/D:" + os.path.join(os.path.curdir,
                             "wpop_{0}.txt".format(wpop_file_root)),
        "/M:" + os.path.join(os.path.curdir, "output", "pop.bin"),
        "/S:" + os.path.join(os.path.curdir, "output", "network.bin"),
        "/R:{0}".format(config["r0"][0] / 2),
        "{0}".format(config["network_seed1"]),
        "{0}".format(config["network_seed2"]),
        "17389101",  # Run seed doesn't matter as we're going to drop this
        "4797132"    # run
        ]

    runner_cmd = args.runner + [
                 "--input",
                 input_setup,
                 "--output",
                 output_geography
                 ] + covidsim_cmd

    print("Starting setup job for geography {0}".format(geography))
    print("  Command line: " + " ".join(runner_cmd))
    setup_jobs[geography] = subprocess.Popen(runner_cmd)


def copy_bin_files(geography):
    """Copy the binary files from the setup run to the input places."""
    input_geography = os.path.join(args.output, "input", geography)
    output_geography = os.path.join(args.output, "output", geography, "setup")

    for f in ["pop.bin", "network.bin"]:
        if not os.path.exists(os.path.join(output_geography, f)):
            print("ERROR: Failed to find {1} for geography: {0}".
                  format(geography, f))
            print("Looked in {0}".format(os.path.join(output_geography, f)))
        else:
            shutil.copyfile(os.path.join(output_geography, f),
                            os.path.join(input_geography, f))


def call_runner(geography, r0, param_file, run_seed1, run_seed2, output):
    """Call the runner for a particular set of values."""
    covidsim_cmd = [
        os.path.join(os.path.curdir, "CovidSim.exe"),
        "/c:{0}".format(config["threads"]),
        "/A:" + os.path.join(os.path.curdir, admin_base(geography)),
        "/PP:" + os.path.join(os.path.curdir, pp_base(geography)),
        "/P:" + os.path.join(os.path.curdir, os.path.basename(param_file)),
        "/O:" + os.path.join(os.path.curdir, "output",
                             os.path.basename(output)),
        "/D:" + os.path.join(os.path.curdir, "pop.bin"),
        "/L:" + os.path.join(os.path.curdir, "network.bin"),
        "/R:{0}".format(r0 / 2),
        "{0}".format(config["network_seed1"]),
        "{0}".format(config["network_seed2"]),
        "{0}".format(run_seed1),
        "{0}".format(run_seed2)
        ]

    runner_cmd = args.runner + [
                 "--input",
                 os.path.join(args.output, "input", geography),
                 "--output",
                 output,
                 ] + covidsim_cmd

    print("  Command line: " + " ".join(runner_cmd))
    return subprocess.Popen(runner_cmd)


def start_many_runs(geography):
    """Given a geography do all the runs for it once setup is complete."""
    copy_bin_files(geography)
    for r0 in config["r0"]:
        for param_file in config["param_files"]:
            pf = os.path.basename(param_file)
            run_seeds = run_seed_reference.copy()
            run_id = 0
            while run_seeds:
                run_seed1 = run_seeds.pop()
                run_seed2 = run_seeds.pop()
                output = os.path.join(args.output, "output", geography,
                                      str(r0), pf,
                                      "{0}_{1}_{2}_{3}".format(geography, r0,
                                                               pf, run_id))
                print("Starting run job for geography {0}, r0={1}, "
                      "param_file={2}, run={3}".format(geography, r0,
                                                       param_file, run_id))
                all_jobs[output] = call_runner(geography, r0, param_file,
                                               run_seed1, run_seed2, output)
                run_id += 1


# Wait for the jobs to end, giving status as we go along
print("Waiting for jobs to end:")
total_pending = 1
message_len = 0
while total_pending:
    setup_pending = len(setup_jobs)
    runs_pending = 0
    for job_id, job in all_jobs.items():
        if job.poll() is None:
            runs_pending += 1

    total_pending = setup_pending + runs_pending
    message = "Number of pending jobs: {0} (Setup={1} Runs={2})".format(
            total_pending, setup_pending, runs_pending)
    extra_space = " " * max(0, message_len - len(message))
    print("\r{0}{1}".format(message, extra_space), end="")
    message_len = len(message)
    do_sleep = total_pending > 0

    # When setup jobs complete we need to start their intervention runs
    for geography, process in setup_jobs.items():
        if process.poll() is not None:
            if process.returncode != 0:
                print("\nERROR: Setup job for geography {0} failed".
                      format(geography))
            else:
                print("\nCompleted setup for geography {0}".format(geography))
                start_many_runs(geography)

            del setup_jobs[geography]

            # We assume that if we've started a set of runs we need to update
            # the status immediately
            do_sleep = False
            break

    # Sleep before we poll again
    if do_sleep:
        time.sleep(10)

print("\nAll runs completed")

# Check exit codes:
any_failed = 0
for job_id, job in all_jobs.items():
    if job.returncode != 0:
        print("ERROR: Run failed: {0}".format(job_id))
        any_failed = 1

if any_failed:
    exit(1)
