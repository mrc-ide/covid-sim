# Data Directory and Samples

This directory provides sample input data and an example script to run.

**IMPORTANT**: The parameter files are provided as a sample only and do not
necessarily reflect runs used in published papers.

## Directories

- *populations*: provides population density information

- *admin_units*: specifies the administration units for each country

- *param_files*: Provides sample parameter files

## Running the sample

The file `./run_sample.py` will build and run a sample for a specified geography.

### Usage

```sh
./run_sample.py <GEOGRAPHY>
```

For example:

```sh
./run_sample.py United_Kingdom
```

This will:

1. Build a copy of the CovidSim executable

2. Do a no intervention base run. This run also provides various `*.bin`
   files for the later run.

3. Do a run with controls applied. This is the step to be repeated when
   doing further tests/runs.

The sample script has more command line options available.  See
`./run_sample.py --help` for a list.

## Running many runs

The Python script [./run_many.py](./run_many.py) is a more advanced version of
`./run_sample.py` and provides a framework for doing large numbers of runs.
It generates `CovidSim` command lines and passes them off to a specified
runner program to invoke.

### Usage

```sh
./run_many.py [--covidsim exe] [--srcdir dir] [--data dir] [--output dir] \
              --runner cmd... config
```

The following example will runn the test-config set of jobs on your local
machine, serialising the runs:

```sh
./run_many.py --output output-test --runner ./runners/run_local.py \
              test-config.json
```

### Command line options

 * `--covidsim exe`: Path to executable to use.  If not present will build a
   default configured executable based off the source specified by `--srcdir`.

 * `--srcdir dir`: Path to source directory, only used if building CovidSim.
   Default is the parent directory of the script.

 * `--data dir`: Directory with all the population, admin, and pre-param data
   files.  Default is the directory of the script.

 * `--output dir`: Directory to use to put output in.

 * `--runner cmd`: Command to use to wrap CovidSim execution in.  May be
   specified multiple times to pass commands to the runner.  For instance if
   you do: `--runner ./runners/run_local.py --runner="-j=3"` the local runner
   will run up to three jobs at a time.  See below for more details, and
   details of the runner interface.

 * `config`: JSON configuration of settings to run.  See below for details.

### Config file schema

The config file is a JSON file configuring what variations to run.  The
following is an example giving all the keys possible.

```json
{
  # Number of threads to allow each CovidSim execution to use.  Defaults to 1
  "threads": 4,

  # The two random-number generator seeds to use for the network generation.
  # If specified must be a list of precisely two numbers.
  # Defaults to two random numbers.
  "network_seeds": [ 98798150, 729101 ],

  # Run counts and run seeds.
  # "run_seeds" contains the list of run seeds (two for each run) to use.
  # "num_runs" is the number of runs to do for each combination of settings.
  # If neither is specified then 1 run is done with random run seeds.
  # If num_runs is specified, but not run_seeds then "num_runs" runs are done
  # with random seeds.
  # If run_seeds is specified then that is used to determine the number of
  # runs.
  # "run_seeds" must be a multiple of two elements long.
  "num_runs": 3,
  "run_seeds": [ 1000000, 200000, 30000, 40000, 5000, 6000028 ],

  # Pre-parameter file to use.  Defaults to geography specific file.
  "pre_param_file" : "pre_param.txt"

  # Population density file to use.  If it has extension .gz will be gunzipped
  # before use.  Defaults to geography specific file.
  "pop_density_file" : "pop.txt.gz"

  # Which values of r0 to test.  These are divided by two when passed to
  # CovidSim
  # Must be specified either globally, or for a geography
  "r0": [ 2.8, 3.0, 3.2 ],

  # Which parameter files to execute.  Paths are relative to the current
  # directory.  You cannot have two files with the same basename.
  # Must be specified either globally, or for each geography.
  "param_files": [ "param_files/p_PC7_CI_HQ_SD.txt" ]

  # Which geographies to run through
  "geographies": {
    "GEO_NAME": {
      # A geography can override any of the settings set globally above.
      # If either "num_runs" or "run_seeds" is specified in a geography
      # then both are considered to be overriden.

      # In addition it can specify the admin file for the geography
      # Defaults to a geography specific value.
      "admin_file" : ".../GEO_NAME_admin.txt"
    },
    "GEO2": { ... }
  }
}
```

See [test-config.json](./test-config.json) for an example configuration file.

### Output

The directory specified by --output is populated by as follows:

 * `config.json`: Canonicalised configuration file.  This can be used as input
   to `run_many.py` to reproduce the build exactly.

 * `input-setup/`: Directory used by initial setup runs.  One directory is
   used by all geographies for setup, to help reduce the disk space required.

 * `input/<geography>/`: Directory used by intervention runs for each
   geography.

 * `output/<geography>/`: Top level output directories for each geography.

   * `setup/`: Results of running setup.  The interesting files in here are
     the binary network and population files which are moved to the
     appropriate input directory.

   * `<R>/<param_file>/<run>`:  Results for a particular R0, param file, run
     ID.  This contains the following files:

     * `*.xls`: Output from CovidSim

     * `stderr.txt`: Standard error from CovidSim

     * `stdout.txt`: Standard output from CovidSim

### Available Runners

The following runners are available in the repo currently:

 * runners/run_local.py: Run CovidSim locally

#### run_local.py

This runner runs CovidSim on your local machine.  As well as the options
required by `run_many.py` it also accepts the following options:

 * `-j N`: Run up to `N` invocations of `CovidSim` at once

 * `--lock dir`: Directory to use to store lock files to manage multiple runs.

 * `--keep`: Keep the temporary directory after finishing.  Useful for post
   mortem debugging.

To make use of `-j` specify it as an additional `--runner` argument to
`run_many.py`.  The following example will allow up to three `CovidSim`
executions at once:

```sh
./run_many.py --output output-test --runner ./runners/run_local.py \
              --runner="-j=3" test-config.json
```

### Writing a new Runner

You can write your own runner to make use of any other method of running
CovidSim that you would like, for example via batch job submissions to a
HPC compute cluster.

`./run_many.py` does not do any job control itself and will launch as many
`runner` processes as it can as quickly as it can.  It expects the `runner` to
control when it actually executes the command required, and to exit only when
that command has completed.

#### Command line interface for runners

Runners should support the following command line interface:

```sh
./runner --input dir --output dir cmd...
```

Runners may accept additional options to enable users to control behaviour
(see [`runners/run_local.py`](./runners/run_local.py) for an example).

##### Option descriptions

 * `--input dir`: Directory where input data lives

 * `--output dir`: Directory where output data should go

 * `cmd...`: Command line to execute.

#### Behaviour

The expected behaviour is for the runner to:

 1. Copy input directory to a temporary directory (run directory)

 2. Change to that run directory

 3. Execute the command, redirecting stdout to output/stdout.txt and stderr to
    output/stderr.txt

 4. Copy the contents of "run directory/output" to the output directory.

 5. Exit with the same exit code as the command.

The runner is responsible for tidying up after itself.

In pseudo-shell this looks something like:

```sh
cp -r "$input" "$temp"
cd "$temp"
$cmd
exit=$?
cp -r "$temp/output" "$output"
rm -rf "$temp"
exit $?
```

The runner *must* treat the input directory as read-only, as the same input
directory is used for multiple invocations.  The output directory should be
treated as write only.
