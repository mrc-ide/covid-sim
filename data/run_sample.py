#!/usr/bin/env python3

import argparse
import gzip
import multiprocessing
import os
import shutil
import subprocess
import sys

def try_remove(f):
    try:
        os.remove(f)
    except OSError as e:
        pass

def parse_args():
    """Parse the arguments.

    On exit: Returns the result of calling argparse.parse()

    args.spatialsim is the name of the spatial sim executable
    args.datadir is the directory with the input data
    args.paramdir is the directory with the parameters in it
    args.outputdir is the directory where output will be stored
    args.threads is the number of threads to use
    """
    parser = argparse.ArgumentParser()
    try:
        cpu_count = len(os.sched_getaffinity(0))
    except AttributeError:
        # os.sched_getaffinity isn't available
        cpu_count = multiprocessing.cpu_count()
    if cpu_count is None or cpu_count == 0:
        cpu_count = 2

    script_path = os.path.dirname(os.path.realpath(__file__))

    # Default values
    data_dir = script_path
    param_dir = os.path.join(script_path, "param_files")
    output_dir = os.getcwd()
    src_dir = os.path.join(data_dir, os.pardir, "src")

    parser.add_argument(
            "country",
            help="Country to run sample for")
    parser.add_argument(
            "--spatialsim",
            help="Location of SpatialSim binary, if none specified will build")
    parser.add_argument(
            "--datadir",
            help="Directory at root of input data",
            default=script_path)
    parser.add_argument(
            "--paramdir",
            help="Directory with input parameter files",
            default=param_dir)
    parser.add_argument(
            "--srcdir",
            help="Directory with source in - needed if --spatialsim isn't specified",
            default=src_dir)
    parser.add_argument(
            "--outputdir",
            help="Directory to store output data",
            default=output_dir)
    parser.add_argument(
            "--threads",
            help="Number of threads to use",
            default=cpu_count
            )
    args = parser.parse_args()

    return args

args = parse_args()

# Determine whether we need to build the tool or use a user supplied one:
if args.spatialsim is not None:
    exe = args.spatialsim
else:
    build_dir = os.path.join(args.outputdir, "build")

    # Build the right version depending on the country
    if args.country == "United_States" or args.country == "Canada":
        country = "-DCOUNTRY=US"
    else:
        country = "-DCOUNTRY=UK"

    # Ensure we do a clean build
    shutil.rmtree(build_dir, ignore_errors=True)
    os.makedirs(build_dir, exist_ok=False)
    os.chdir(build_dir)

    # Build
    subprocess.run(['cmake', args.srcdir, country], check=True)
    subprocess.run(['cmake', '--build', '.'], check=True)

    # Where the exe ends up depends on the OS.
    if os.name == 'nt':
        exe = os.path.join(build_dir, "Debug", "SpatialSim.exe")
    else:
        exe = os.path.join(build_dir, "SpatialSim")

# Ensure output directory exists
os.makedirs(args.outputdir, exist_ok=True)

# The admin file to use
admin_file = os.path.join(args.datadir, "admin_units",
    "{0}_admin.txt".format(args.country))

if not os.path.exists(admin_file):
    print("Unable to find admin file for country: {0}".format(args.country))
    print("Data directory: {0}".format(args.datadir))
    print("Looked for: {0}".format(admin_file))
    exit(1)

# Population density file in gziped form, text file, and binary file as
# processed by SpatialSim
if args.country == "United_States" or args.country == "Canada":
    wpop_file_gz = os.path.join(args.datadir, "populations", "wpop_usacan.txt.gz")
else:
    wpop_file_gz = os.path.join(args.datadir, "populations", "wpop_eur.txt.gz")

if not os.path.exists(wpop_file_gz):
    print("Unable to find population file for country: {0}".format(args.country))
    print("Data directory: {0}".format(args.datadir))
    print("Looked for: {0}".format(wpop_file_gz))
    exit(1)

wpop_file = os.path.join(args.outputdir, "wpop_eur.txt")
wpop_bin = os.path.join(args.outputdir, "{0}_pop_density.bin".format(args.country))

# gunzip wpop fie
try_remove(wpop_file)
try_remove(wpop_bin)
with gzip.open(wpop_file_gz, 'rb') as f_in:
    with open(wpop_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

# Configure pre-parameter file.  This file doesn't change between runs:
if args.country == "United_States":
    pp_file = os.path.join(args.paramdir, "preUS_R0=2.0.txt")
else:
    pp_file = os.path.join(args.paramdir, "pre_R0=2.0.txt")
if not os.path.exists(pp_file):
    print("Unable to find pre-parameter file")
    print("Param directory: {0}".format(args.paramdir))
    print("Looked for: {0}".format(pp_file))
    exit(1)

# Configure No intervention parameter file.  This is run first
# and provides a baseline
no_int_file = os.path.join(args.paramdir, "p_NoInt.txt")
if not os.path.exists(no_int_file):
    print("Unable to find parameter file")
    print("Param directory: {0}".format(args.paramdir))
    print("Looked for: {0}".format(no_int_file))
    exit(1)

# Configure an intervention (controls) parameter file.
# In reality you will run SpatialSim many times with different parameter
# controls.
control_root = "PC7_CI_HQ_SD"
controls_file = os.path.join(args.paramdir, "p_{0}.txt".format(control_root))
if not os.path.exists(controls_file):
    print("Unable to find parameter file")
    print("Param directory: {0}".format(args.paramdir))
    print("Looked for: {0}".format(controls_file))
    exit(1)

school_file = None
if args.country == "United_States":
    school_file = os.path.join(args.datadir, "populations", "USschools.txt")

    if not os.path.exists(school_file):
        print("Unable to find school file for country: {0}".format(args.country))
        print("Data directory: {0}".format(args.datadir))
        print("Looked for: {0}".format(school_file))
        exit(1)

# Some command_line settings
r = 3.0
rs = r/2

# This is the temporary network that represents initial state of the
# simulation
network_bin = os.path.join(
        args.outputdir,
        "Network_{0}_T{1}_R{2}.bin".format(args.country, args.threads, r))
try_remove(network_bin)

# Run the no intervention sim.  This also does some extra setup which is one
# off for each R.
print("No intervention: {1} {0}".format(r, args.country))
cmd = [
        exe,
        "/c:{0}".format(args.threads),
        "/A:" + admin_file 
]
if school_file:
    cmd.extend(["/s:" + school_file])
cmd.extend([
        "/PP:" + pp_file, # Preparam file
        "/P:" + no_int_file, # Param file
        "/O:" + os.path.join(args.outputdir,
            "{0}_NoInt_R0={1}".format(args.country, r)), # Output
        "/D:" + wpop_file, # Input (this time text) pop density
        "/M:" + wpop_bin, # Where to save binary pop density
        "/S:" + network_bin, # Where to save binary net setup
        "/R:{0}".format(rs),
        "98798150", # These four numbers are RNG seeds
        "729101",
        "17389101",
        "4797132"
        ])
print("Command line: " + " ".join(cmd))
process = subprocess.run(cmd, check=True)

# Parameters that can be chanegd within the run
x = 0
q = 182
qo = q

# In reality loop over this stage for multple parameter files and settings of
# /CLP#:
print("Intervention: {4} {0} {1} {2} {3} {3}".format(r, x, q, qo, control_root))
cmd = [
        exe,
        "/c:{0}".format(args.threads),
        "/A:" + admin_file,
        "/PP:" + pp_file,
        "/P:" + controls_file,
        "/CLP1:{0}".format(x),
        "/CLP2:{0}".format(q),
        "/CLP3:{0}".format(qo),
        "/CLP4:{0}".format(qo),
        "/O:" + os.path.join(args.outputdir,
            "{4}_{0}_{1}_{2}_R0={3}".format(control_root, x, q, r, args.country)),
        "/D:" + wpop_bin, # Binary pop density file (speedup)
        "/L:" + network_bin, # Network to load
        "/R:{0}".format(rs),
        "98798150",
        "729101",
        "17389101",
        "4797132"
        ]
print("Command line: " + " ".join(cmd))
process = subprocess.run(cmd, check=True)
