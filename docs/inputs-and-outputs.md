# The inputs and outputs of the `CovidSim` model

This is WIP. Know something not documented here? Please add and open a PR!

- Table of contents
  - [The geography](#the-geography)
  - [Main command-line arguments](#main-command-line-arguments)
  - [Input files](#input-files)
    - [Parameters](#parameters)
    - [Parameter files](#parameter-files)
    - [Population density file](#population-density-file)
      - [How population density files are produced](#how-population-density-files-are-produced)
    - [School files](#school-files)
  - [Output files](#output-files)
  - [R summary visualisations](#r-summary-visualisations)

## The geography

`CovidSim` simulates disease spread in a geographical region, which in principle
can be at any scale, but in practice is a region or country.

In consequence, the model must be told the geography of a region, such as its
population density, plus other specific information. This information is
specified as a mixture of parameters and input population density files.

## Main command-line arguments

A typical run specifies:

1. Files that contain simulation parameters (the `/A`, `/P` and `/PP` options)
2. A population density file for the country we're simulating (the `/D` option)
3. The name of output files that summarise the results of the simulation (the `/O` option).

```shell
CovidSim
    /O:OutputFilesPrefix
    /P:ParameterFile
    [/NR:NumberOfRealisations]
    [/A:AdminParamFile]
    [/AP:AirTravelFile]
    [/c:NumThreads]
    [/C:PlaceCloseIndepThresh]
    [/CLP[1-6]:ParamOverrideNumber]
    [/d:RegionalDemographyFile]
    [/D:PopulationDensityFile]
    [/I:InterventionFile]
    [/KO:KernelOffsetScale]
    [/KP:KernelPowerScale]
    [/L:NetworkFileToLoad]
    [/LS:SnapshotLoadFile]
    [/M:OutputDensityFile]
    [/PP:PreParameterFile]
    [/R:R0scaling]
    [/s:SchoolFile]
    [/S:NetworkFileToSave]
    [/T:PreControlClusterIdCaseThreshold]
    SetupSeed1 SetupSeed2 RunSeed1 RunSeed2
```
Required arguments:

- `/O` - Output file path prefix for simulation data collection. Output file
  names have the `.xls` extension but use tabular `tsv` data.
  - Example: `/O:./output/NoInt_R0=1`
- `/P` - Intervention parameters for a specific run.
  - Example: `/P:./data/param_files/p_NoInt.txt`
- `SetupSeed1 SetupSeed2` - Random number generator seeds used when initialising
  the model, including creating the network file (large positive integers).
- `RunSeed1 RunSeed2` Random number generator seeds used when running the model.
  These can be varied to do multiple runs with the same network file
  (large positive integers).

Optional Arguments:

- `/NR` - specified the number of simulation realisations (independent runs with the same
parameters) to run at once and average over in the output files.
- `/A` - [Administrative division](./glossary.md#Administrative\ Division) parameter file
  - Example: `/A:./data/admin_units/United_Kingdom_admin.txt`
- `/AP` Air travel data for a specific geography (unused currently)
- `/BM:format`.  Specifies the output bitmap format.  Valid choices are `BMP`, or (when
 available) `PNG` - Default is `PNG` if available, otherwise `BMP`.
- `/c` - Number of parallel threads to use (only used if compiled with OpenMP)
  - Example: `/c:32`
- `/C` - Sets the `P.PlaceCloseIndepThresh` parameter.
- `/CLP[1-6]` - Special parameters that interact with wildcards `#1`, `#2`, etc.
  in the intervention parameter file (and less often the pre-parameter file).
  Wildcard `#n` is replaced by the value of `CLPn`. This is useful to vary parts
  of parameter files at run-time (e.g. to undertake sensitivity analysis)
  without needing to generate entirely new parameter files.
  - Examples: `/CLP1:100000` & `/CLP2:0`
- `/d` - Regional demography file to use.
- `/D` - Population density file for a specific geography (e.g. a country). Can
  be loaded from either the original textual format or a binary format from
  a previous run that used the `/M` option.
  - Examples: `/D:./data/populations/wpop_eur.txt` & `/D:./US_LS2018.bin`
- `/I` - Intervention file. Can be specified more than once.
- `/KO` - Scales the `P.MoveKernelScale` parameter.
- `/KP` - Scales the `P.MoveKernelShape` parameter.
- `/L` - Load a network file saved from a previous run that specified `/S`.
  - Example: `/L:./network_file.bin`
- `/LS` - Load a snapshot file saved by the `/SS` command.
  - Example: `/LS:./snapshot.bin`
- `/M` - Output a population density file to disk
  - Example: `/M:./US_LS2018.bin`
- `/PP` - Transmission and calibration parameter files for a specific run
  - Example: `/PP:./data/param_files/preUS_R0=2.0.txt`
- `/R`. Specifies the basic reproduction number [R0](./glossary.md#R0), as a
  multiplier of 2. This command-line parameter is read into `P.R0scaling` which
  scales the R0 parameter specified in the parameter file. This is useful when
  repeating simulations that *only* vary `R0`). For COVID-19, 1.4 to 1.6 is suitable.
  - Example: `/R:1.6`
- `/s` - School information for a specific geography (currently only used for US).
  - Example: `/s:./data/populations/USschools.txt`
- `/S` - For efficiency, we can run and, as a side-effect, generate a
  [network file](./model-glossary.md#Network-file) that assigns
  [people](./model-glossary.md#People) to [places](./model-glossary.md#Places).
  It may then be re-used for subsequent runs with different input parameters for
  the same geography. ***Note***: this file is non-portable
  - Example: `/S:./network_file.bin`
- `/SS` - Specifies the file and interval at which to save a snapshot when
  running a simulation. The first argument is the number of `P.TimeStep`s that
  should elapse before saving. The second argument is the file to save snapshots
  at.
  - Example: `/SS:100,./snapshot.bin`
- `/T` - Sets the `P.PreControlClusterIdCaseThreshold` parameter.

## Input files

The main inputs files are parameter files and population density files
(for specific geographies).

### Parameters

There are a very large number of parameters to `CovidSim`. This repo is
undergoing active development and rationalisation. The parameters are currently
not self-documenting.

Parameter values are read in from parameter files by function `ReadParams`,
which matches up a parameter description string to the according variable in the
source code. The only method to determine the precise meaning of a specific
parameter is to read the code.

### Parameter files

The parameters are specified in admin, pre-parameter and intervention parameter
files. Both files have the same format.

Admin and pre-parameter files contain parameters whose values are *common* to a
series of runs (i.e. defining geographies and transmission parameters).
Intervention Parameter files group intervention parameters whose values are more
likely to *differ* between a series of runs.

The format is a sequence of:

```
[Description of Parameter]
value
```

If you see multiple numbers below the parameter description, then disregard them.
The simulation uses only the numbers immediately below the parameter description.

An example parameter file is `./data/param_files/p_NoInt.txt`.

### Population density file

A binary geography-specific file used to assign people to cells. Currently these
files are generated and provided by Imperial College.

An example population density file is `./data/populations/wpop_eur.txt`.

The information contained in this file includes:

| longitude | latitude | number of people | country code | admin unit code |
|--|-:|-:|-:|-:|
| -156.68333 | 71.325| 30| 46 | 4602017|
|-156.76666| 71.3 | 1 |46 | 4602017|
| ... | ... | ... | ...| ... |

#### How population density files are produced

Physical geography data: each geography has a shape file (`.shp`) of polygons
and meta-data (`.dbf`) with GPS coordinates. Admin units are a set of polygons.

Human geography data specifies where people live on the same scale as a
`CovidSim`'s [microcell](./model-glossary.md#Microcells) (1/120th of a degree).

Imperial College combines the physical and human data to calculate population
densities per polygon. This process produces the population density file.

A companion to the population density file is a meta-file that maps admin unit
codes to string descriptions (e.g., codes to US state names).

### School files

The first line of a school file has (1 + 2`n`) integer values, where `n` is
the number of school types. The values are:

- Index `0`: The number of types of schools. E.g. a geography might two school
  place types (primary and secondary).
- Index 1 + 2`i`: The total number of schools of type `i`
- Index 2 + 2`i`: The number of age bands in schools of type `i`

E.g., if a geography has 2 school types then the first line of the school file
might be:

`2 100 3 50 4`

representing 2 school types, with 100 of type 0 (which as 3 age classes) and 50
of type 1 (which has 4 age classes).

The remainder of the file has a row per school. E.g.:

| longitude | latitude | place type index | #people in the school | #people in age band 1| # people in age band 2 | ... | # people in age band n |
|-:|-:|-:|-:|-:|-:|-:|-:|
| -156.68333 | 71.325| 0 | 80 | 30 | 46 | ... | 4 |
| -123.32 | 70.35 | 0 | 32 | 23 | 3 | ... | 6 |
| ... | ... | ... | ... | ... | ... | ... | ... |

The place type index for schools is `0`.

## Output files

Simulation output files are produced by each run.
Switches in parameter files can control the precise nature of the outputs
(e.g., at country level, or at admin unit level, or both etc.). E.g.

```
[Do Severity Analysis]
1
```

then `severity.xls` is generated.

A run is extinct if the disease dies out, otherwise a run is non extinct.

Outputs can be averaged over all extinct (`avE` suffix) and non-extinct
(`avNE` suffix) runs. Currently, we are simulating large epidemics that
essentially become deterministic and therefore we focus on `avNE` files.

We pay most attention to `avNE` (average of non-extinct realisations) files.

Below is an incomplete specification of the output file formats.

### `name.avNE.xls`

Contains time-stamped (e.g., daily) statistics for the simulation over the whole country.

| column | meaning |
| ------------- |-------------:|
| t | sample time – specified in the preparam file by Sampling timestep - generally day in 2020 (t=1 -> Jan 1)  |
| S | total number of susceptibles in the population |
| L | total number of latently infected people in the population |
| I | total number of infectious people in the population |
| R | total number of recovered people in the population |
| D | total number of deaths in the population |
| incI | incidence of infections at that timestep |
| incR | incidence of recoveries |
| incFC | incidence of false cases, i.e. false positives |
| incC | incidence of cases |
| incDC | incidence of detected cases |
| incTC | incidence of treated cases |
| incH | incidence of hospitalisations – again, probably can ignore this as was written specifically for the Ebola model and we’re using a  different approach here. |
| cumT | cumulative number of treated cases |
| cumTmax | the maximum number of cumulative treated cases from the runs being averaged over |
| cumTP | cumulative number of privately treated cases |
| cumV | cumulative number of vaccinations |
| cumVmax |  the maximum number of cumulative vaccinations  from the runs being averaged over |
| Extinct | Is the run extinct or not? |
| rmsRad | root mean square radius of infections from seed point |
| maxRad | maximum radius of an infection from the seed point |
| v* | a sequence of columns containing the variance of the above quantities in the same order (excluding the time step) |
| value 1 | Number of non-extinct runs |
| value 2 | Number of extinct runs |
| value 3 | R0 in households |
| value 4 | R0 in places |
| value 5 | R0 of spatial transmission |
| value 6 | Mean peak height |
| value 7 | Variance of peak height |
| value 8 | Mean peak time |
| value 9 | Variance of peak time |

### `name.avNE.adunit.xls`

Contains time-stamped statistics per [admin unit](./model-glossary.md#Admin-unit)
(hopefully with headers matching the codes in a population index file).

| column | meaning |
| ------------- |-------------:|
| t |  time |
| I(admincode) ... | Incidence of infection in each admin unit (the number of columns equals the number of admin units used) |
| C(admincode) ... | Incidence of cases in each admin unit. |
| DC(admincode) ... | Incidence of detected cases in each admin unit |
| T(admincode) ... | Incidence of treated cases in each admin unit |
| value ... | A sequence of column values of the population of each admin unit |

### `name.avNE.age.xls`

| column | meaning |
| ------------- |-------------:|
| t |  time |
| I(age band) ... | incidence of cases in each age band |
| C(age band) ... | incidence of critical cases in each age band |
| D(age band) ... | incidence of deaths in each age band |

### `name.avNE.severity.xls`

Contains statistics on the [prevalence](./model-glossary.md#Prevalence) of the
infection.

| column | meaning |
| ------------- |-------------:|
| t | time |
| PropSchClosed | proportion of schools closed |
| PropSocDist | unknown |
| [mild](./model-glossary.md#Mild) |  total number of mild cases at time t |
| [ILI](./model-glossary.md#ILI) | total number of influenza-like illness cases at time t (assume represents GP demand) |
| [SARI](./model-glossary.md#SARI) | total number of severe acute respiratory illness cases at time t (assume represents hospital demand) |
| [Crit](./model-glossary.md#Crit) | total number of critical cases (assume represents ICU demand) |
| [CritRecov](./model-glossary.md#CritRecovery) | total number of critical cases who are well enough to be out of ICU but still need a hospital bed |
| incMild | incidence of mild cases |
| incILI | incidence of ILI cases |
| incSARI | incidence of SARI cases |
| incCrit | incidence of critical cases |
| incCritRecov | incidence of critical cases still in hospital but no longer requiring ICU |
| incDeath | incidence of death |
| cumMild | cumulative number of mild cases |
| cumILI | cumulative number of ILI cases |
| cumSARI | cumulative number of SARI cases |
| cumCrit | cumulative number of critical cases |
| cumCritRecov | cumulative number of critical cases still in hospital but no longer requiring ICU |
| v* | a sequence of columns containing the variance of the above quantities in the same order (excluding the PropSchClosed, PropSocDist) |

### `name.avNE.severity.adunit.xls`

As per `name.avNE.serverity.xls`, excluding PropSchClosed and PropSocDist, and
with each quantity listed for each admin unit in turn.

<!--
### `name.avNE.adunitVar.xls`

### `name.avNE.controls.xls`

### `name.adunit.xls`

### `name.avNE.country.xls`

### `name.avNE.household.xls`

### `name.avNE.inftype.xls`

### `name.avNE.R0.xls`

### `name.avNE.severity.xls`

### `name.severity.adunit.xls`

### `name.severity.xls`

### `name.xls`
-->

## R summary visualisations

Some [R scripts](../Rscripts) provide basic visualisations of model runs.

If the R software is installed and output files of model runs have been created
in folder `folder`, they can be visualised using the commands

```shell
Rscript Rscripts/PlotsSpatial.R [folder-where-the-data-is]
Rscript Rscripts/CompareScenarios.R [folder-where-the-data-is]
```

This will create `.png`s visualising the data in a new subfolder called `Plots`.
