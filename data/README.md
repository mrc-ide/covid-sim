# Data Directory and Samples

This directory provides sample input data and an example script to run.

**IMPORTANT**: The parameter files are provided as a sample only and do not
necessarily reflect runs used in published papers.

## Directories

 * *populations*: provides population density information

 * *admin_units*: specifies the administration units for each country

 * *param_files*: Provides sample parameter files

## Running the sample.

The file `./run_sample.py` will build and run a sample for a specified country.

### Usage:

```sh
./run_sample.py <COUNTRY>
```

For example:

```sh
./run_sample.py United_Kingdom
```

This will:

 1. Build a copy of the SpatialSim executable configured for the appropriate
    country (in this case the UK)

 2. Do a no intervention base run.  This run also provides various `*.bin`
    files for the later run.

 3. Do a run with controls applied.  This is the step to be repeated when
    doing further tests/runs.

The sample script has more command line options available.  See 
`./run_sample.py --help` for a list.
