#!/usr/bin/env python3
r"""Do multiple runs locally or remotely.

This Python script provides a framework for doing multiple runs of CovidSim
either locally or remotely.

It should be invoked as follows:

    python3 run_many.py [--covidsim <binary>] [--data <dir>] \
                        [--output <dir>] --runner runner_script... \
                        [--srcdir <srcdir>] config.json

For example:
    python3 run_many.py --runner runners/run_local.py run_us.json

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

    3. Execute <command...>

    4. Copy/move the (remote) contents of ./output/ to <output-dir>.

    5. Exit (0 for success, non-zero for failure).

The JSON configuration takes a set of config options which control what gets
invoked:

{
    # The following options may all be overriden in the geographies map.
    # All optional unless otherwise specified.
    "pre_param_file": <string> # Pre-param file for the geography
    "pop_density_file": <string> # Population density file
    "threads" : <NUM> # Number of threads to specify for CovidSim to use
    "network_seeds" : [<NUM>, <NUM>] # Network generation seeds
    "r0" : [<r0>] # R0 values to test. Must be specified
    "param_files" : [<Param file>] # Parameter files to test.  Must be
                                   # specified
    "num_runs" : <NUM> # Number of runs with different run seeds
    "run_seeds" : [<NUM>, ...] # 2 * num_run number of seeds for the runs

    "geographies" : {
        <geography> : {
            # The following option defaults to a "sensible" values based
            # on the geography if not specified
            "admin_file": <string> # Admin file for the geography

            # The followng two options default to the value specified at the
            # top level if not specified here, and if not specified at the top
            # level will use a geography dependent default.
            "pre_param_file": <string> # Pre-param file for the geography
            "pop_density_file": <string> # Population density file

            # All the following options override those given at a global
            # level.
            "threads" : <NUM> # Number threads to specify for CovidSim to use
            "network_seeds" : [<NUM>, <NUM>] # Network generation seeds
            "r0" : [<r0>] # R0 values to test.
            "param_files" : [<Param file>] # Parameter files to test.
            "num_runs" : <NUM> # Number of runs with different run seeds
            "run_seeds" : [<NUM>, ...] # 2 * num_run seeds for the runs
        }
        ...
    }
}

There will be num_runs runs for each combination of r0, param-files and
geography.
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


def default_admin_file(data_dir, geography):
    """Get the default admin file for a geography."""
    return os.path.join(data_dir, "admin_units",
                        "{0}_admin.txt".format(geography))


def default_pp_file(data_dir, geography):
    """Get the default pre-param file for a geography."""
    if geography in united_states:
        pp_base = "preUS_R0=2.0.txt"
    elif geography in nigeria:
        pp_base = "preNGA_R0=2.0.txt"
    else:
        pp_base = "preUK_R0=2.0.txt"

    return os.path.join(data_dir, "param_files", pp_base)


def pop_base(geography):
    """Get the basename of the geography GZIP file."""
    if geography in united_states + canada:
        return "wpop_usacan.txt.gz"
    elif geography in usa_territories:
        return "wpop_us_terr.txt.gz"
    elif geography in nigeria:
        return "wpop_nga_adm1.txt.gz"
    else:
        return "wpop_eur.txt.gz"


def default_pop_file(data_dir, geography):
    """Get the default population density file for a geography."""
    return os.path.join(data_dir, "populations", pop_base(geography))


def generate_random_run_seeds(num_runs):
    """Generate a list of 2 * num_runs random seeds.

    We need two run seeds for each run.
    """
    random.seed()
    return [random.randint(0, 0x7fffffff) for i in range(num_runs * 2)]


def validate_run_seeds(config, geography=None):
    """Validate num_runs & run_seeds for a particular geography.

    If geography is None then assume this is the global defaults.
    """
    if "num_runs" not in config and "run_seeds" not in config:
        # Do one run if no values are specified
        config["num_runs"] = 1

    if "run_seeds" not in config:
        config["run_seeds"] = generate_random_run_seeds(config["num_runs"])

    if len(config["run_seeds"]) % 2 != 0:
        print("ERROR: Must specify a multiple of two run seeds in config:")
        print("Specified: {0}".format(config["run_seeds"]))
        if geography is not None:
            print("Geography: {0}".format(geography))
        exit(1)

    if "num_runs" not in config:
        config["num_runs"] = len(config["run_seeds"]) // 2

    if config["num_runs"] * 2 != len(config["run_seeds"]):
        print("WARNING: Number of runs does not match length of run_seeds")
        if geography is not None:
            print("Geography: {0}".format(geography))
        print("{0} * 2 != {1}".format(config["num_runs"],
                                      len(config["run_seeds"])))
        print("Will use length of run_seeds")
        config["num_runs"] = len(config["run_seeds"]) // 2


def read_config(args):
    """Read, validate, and canonicalise the configuration file."""
    if not os.path.exists(args.config):
        print("ERROR: Unable to find config file: {0}".format(args.config))
        exit(1)

    with open(args.config, 'r') as cf:
        config = json.load(cf)

    # If we have no threads variable set use 1:
    if "threads" not in config:
        config["threads"] = 1

    # Need exactly two network seeds
    if "network_seeds" not in config:
        config["network_seeds"] = generate_random_run_seeds(1)

    if len(config["network_seeds"]) != 2:
        print("ERROR: Must specify exactly two network seeds in config:")
        print("Speficied: {0}".format(config["network_seeds"]))
        print("Config file: {0}".format(args.config))
        exit(1)

    validate_run_seeds(config)

    # Ensure we have some geographies
    if "geographies" not in config:
        print("ERROR: No geographies speicified")
        exit(1)

    # Validate individual geographies
    for geography, geography_config in config["geographies"].items():
        if "admin_file" not in geography_config:
            geography_config["admin_file"] = \
                    default_admin_file(args.data, geography)

        if "pre_param_file" not in geography_config:
            if "pre_param_file" not in config:
                geography_config["pre_param_file"] = default_pp_file(args.data,
                                                                     geography)
            else:
                geography_config["pre_param_file"] = config["pre_param_file"]

        if "pop_density_file" not in geography_config:
            if "pop_density_file" not in config:
                geography_config["pop_density_file"] = \
                        default_pop_file(args.data, geography)
            else:
                geography_config["pop_density_file"] = \
                        config["pop_density_file"]

        for i in ["threads", "network_seeds", "r0", "param_files"]:
            if i not in geography_config:
                if i not in config:
                    print("ERROR: Unable to find value for {0}".format(i))
                    print("Geography being considered: {0}".format(geography))
                    exit(1)
                else:
                    geography_config[i] = config[i]

        if "num_runs" not in geography_config and \
           "run_seeds" not in geography_config:
            geography_config["run_seeds"] = config["run_seeds"]
            geography_config["num_runs"] = config["num_runs"]

        validate_run_seeds(geography_config, geography)

    return config


def write_config(fname, config):
    """Write the JSON config to fname."""
    with open(fname, 'w') as cf:
        json.dump(config, cf, indent=2)


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
config = read_config(args)

# Write it out to enable repeatable builds
write_config(os.path.join(args.output, "config.json"), config)

# Generate the input data directories.  We have the following layout of the
# output directory:
#  - <args.output>/
#    - input-setup/ - Input files for setup runs
#    - input/<geography>/ - Input files for "true" geography runs"
#    - output/<geography>/ - Output from each geography run
#      - setup/ - Output from setup run
#      - <r0>/<param>/<run> - Output from a particular run.
#
# We use one input directory for all setup runs because the population density
# files are likely shared, and are large.
#
# We keep a dict of all files we've copied so that we know if we duplicate any
# copying.
input_setup = os.path.join(args.output, "input-setup")
os.makedirs(input_setup, exist_ok=True)

copied_files = {}  # Dest -> source of files we've copied


def do_copy(src, dst):
    """Return True if we should copy src to dst.

    Exits if we have a hard error
    """
    if os.path.exists(dst):
        if dst in copied_files:
            if copied_files[dst] != src:
                print("ERROR: Copying different files to {0}".format(dst))
                print("Previous source: {0}".format(copied_files[dst]))
                print("New source: {0}".format(src))
                exit(1)
            else:
                return False
        else:
            print("WARNING: Overwriting existing file: {0}".format(dst))

    copied_files[dst] = src
    return True


def copy_file(src, dst):
    """Copy src to dst, with additional safety checks.

    If dst exists and we haven't copied it before we overwrite dst, otherwise
    if dst exists and we didn't copy it from src before we error out.
    """
    if do_copy(src, dst):
        print("Copying {0} to {1}".format(src, dst))
        shutil.copy2(src, dst)


def gunzip_file(src, dst):
    """Unzip src into dst, with additional safety checks."""
    if do_copy(src, dst):
        with gzip.open(src, 'rb') as fi, open(dst, 'wb') as fo:
            print("Unzipping {0} to {1}".format(src, dst))
            shutil.copyfileobj(fi, fo)


# Copy the executable
copy_file(exe, os.path.join(input_setup, "CovidSim.exe"))

# Dictionaries of jobs:
setup_jobs = {}  # Geography -> Setup run, emptied as jobs finish
all_jobs = {}  # ID -> Intervention run

# Run through the setup runs:
for geography, geography_config in config["geographies"].items():
    input_geography = os.path.join(args.output, "input", geography)
    output_geography = os.path.join(args.output, "output", geography, "setup")
    os.makedirs(input_geography, exist_ok=True)
    os.makedirs(output_geography, exist_ok=True)

    # Copy the executable
    copy_file(exe, os.path.join(input_geography, "CovidSim.exe"))

    # The admin file to use
    admin_file = geography_config["admin_file"]
    admin_base = os.path.basename(admin_file)

    if not os.path.exists(admin_file):
        print("Unable to find admin file for geography: {0}".format(geography))
        print("Looked for: {0}".format(admin_file))
        exit(1)

    copy_file(admin_file, os.path.join(input_setup, admin_base))
    copy_file(admin_file, os.path.join(input_geography, "admin.txt"))

    # Population density file in gziped form, text file, and binary file as
    # processed by CovidSim
    wpop_file_gz = geography_config["pop_density_file"]
    if not os.path.exists(wpop_file_gz):
        print("Unable to find population file for geography: {0}".
              format(geography))
        print("Looked for: {0}".format(wpop_file_gz))
        exit(1)

    wpop_file_base, ext = os.path.splitext(os.path.basename(wpop_file_gz))
    wpop_file = os.path.join(input_setup, wpop_file_base)
    if ext == '.gz':
        gunzip_file(wpop_file_gz, wpop_file)
    else:
        copy_file(wpop_file_gz, wpop_file)

    # Configure pre-parameter file.  This file doesn't change between runs:

    pp_file = geography_config["pre_param_file"]
    pp_base = os.path.basename(pp_file)
    if not os.path.exists(pp_file):
        print("Unable to find pre-parameter file for geography "
              "{0}".format(geography))
        print("Looked for: {0}".format(pp_file))
        exit(1)

    copy_file(pp_file, os.path.join(input_setup, pp_base))
    copy_file(pp_file, os.path.join(input_geography, "preparam.txt"))

    # Copy the parameter files
    # We copy all of them into each Geo input directory, but only the first
    # into the input-setup directory.
    for param_file in geography_config["param_files"]:
        if not os.path.exists(param_file):
            print("Unable to find parameter file: {0}".format(param_file))
            exit(1)

        param_base = os.path.basename(param_file)
        copy_file(param_file, os.path.join(input_geography, param_base))

    param_base = os.path.basename(geography_config["param_files"][0])
    shutil.copyfile(geography_config["param_files"][0],
                    os.path.join(input_setup, param_base))

    covidsim_cmd = [
        os.path.join(os.path.curdir, "CovidSim.exe"),
        "/c:{0}".format(geography_config["threads"]),
        "/A:" + os.path.join(os.path.curdir, admin_base),
        "/PP:" + os.path.join(os.path.curdir, pp_base),
        "/P:" + os.path.join(os.path.curdir, param_base),
        "/O:" + os.path.join(os.path.curdir, "output", "setup"),
        "/D:" + os.path.join(os.path.curdir, wpop_file_base),
        "/M:" + os.path.join(os.path.curdir, "output", "pop.bin"),
        "/S:" + os.path.join(os.path.curdir, "output", "network.bin"),
        "/R:{0}".format(geography_config["r0"][0] / 2),
        "{0}".format(geography_config["network_seeds"][0]),
        "{0}".format(geography_config["network_seeds"][1]),
        "{0}".format(geography_config["run_seeds"][0]),
        "{0}".format(geography_config["run_seeds"][1])
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


def call_runner(geography, threads, r0, param_file, network_seed1,
                network_seed2, run_seed1, run_seed2, run_id):
    """Call the runner for a particular set of values."""
    param_base = os.path.basename(param_file)
    output = os.path.join(args.output, "output", geography, str(r0),
                          param_base, str(run_id))
    prefix = "{0}_{1}_{2}_{3}".format(geography, r0, param_base, run_id)
    print("Starting run job for geography {0}, r0={1}, param_file={2}, "
          "run={3}".format(geography, r0, param_file, run_id))
    covidsim_cmd = [
        os.path.join(os.path.curdir, "CovidSim.exe"),
        "/c:{0}".format(threads),
        "/A:" + os.path.join(os.path.curdir, "admin.txt"),
        "/PP:" + os.path.join(os.path.curdir, "preparam.txt"),
        "/P:" + os.path.join(os.path.curdir, param_base),
        "/O:" + os.path.join(os.path.curdir, "output", prefix),
        "/D:" + os.path.join(os.path.curdir, "pop.bin"),
        "/L:" + os.path.join(os.path.curdir, "network.bin"),
        "/R:{0}".format(r0 / 2),
        "{0}".format(network_seed1),
        "{0}".format(network_seed2),
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
    all_jobs[prefix] = subprocess.Popen(runner_cmd)


def start_many_runs(geography, geography_config):
    """Given a geography do all the runs for it once setup is complete."""
    copy_bin_files(geography)
    for r0 in geography_config["r0"]:
        for param_file in geography_config["param_files"]:
            run_seeds = geography_config["run_seeds"].copy()
            run_id = 0
            # Hack to iterate over run_seeds two at a time
            for run_seed1, run_seed2 in \
                    zip(*[iter(geography_config["run_seeds"])]*2):
                call_runner(
                        geography,
                        geography_config["threads"],
                        r0,
                        param_file,
                        geography_config["network_seeds"][0],
                        geography_config["network_seeds"][1],
                        run_seed1,
                        run_seed2,
                        run_id)
                run_id += 1


# Wait for the jobs to end, giving status as we go along
print("Waiting for jobs to end:")
total_pending = 1
message_len = 0
any_failed = 0
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
                any_failed = 1
            else:
                print("\nCompleted setup for geography {0}".format(geography))
                start_many_runs(geography, config["geographies"][geography])

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
for job_id, job in all_jobs.items():
    if job.returncode != 0:
        print("ERROR: Run failed: {0}".format(job_id))
        any_failed = 1

if any_failed:
    exit(1)
