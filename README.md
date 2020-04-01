# COVID-19 SpatialSim Model

This is the COVID-19 SpatialSim model developed by the MRC Centre for Global
Infectious Disease Analysis hosted at Imperial College, London.

SpatialSim models the instances and severity of COVID-19 infections throughout
a population over time.  It enables modelling of how intervention policies and
healthcare provision affect the spread of COVID-19.

## Status

This model is in active development and subject to significant code changes
to:

 * Enable modelling of more geographies

 * Enable modelling of different intervention scenarios

 * Reduce memory usage

 * Improve performance

## Building

The model is written in C++ and runs on Windows and Linux-based systems.

Running the model for the whole of the UK requires approximately 20GB of RAM.
Other regions will require different amounts of memory (some up to 256GB).

It is strongly recommended to build the model with OpenMP support enabled to
improve performance.

See [build.md](./docs/build.md) for detailed build instructions.

### Testing

The [regressiontest_UK_100th.py](./tests/regressiontest_UK_100th.py) script
provides a basic regression test.

*IMPORTANT*: This script uses test data only and is not a run reflective of
real-world situations.

## Sample Data

The directory [data](./data) contains sample data.

The Python script [run_sample.py](./data/run_sample.py) demonstrates how to invoke
SpatialSim to use this data.  See the [sample README](./data/README.md) for
details on how to run the samples.

## Documentation

Model documentation can be found in the [docs](./docs) directory.  Of
particular interest are:

 * [Model Overview](./docs/model-overview.md)
 * [Model Glossary](./docs/model-glossary.md)
 * [Model Inputs and Outputs](./docs/inputs-and-outputs.md)
 * [R Scripts for Output Visualisation](./docs/inputs-and-outputs.md#r-summary-visualisations)

### Relevant papers

The following papers are relevant to the model.  Please note that some of them
may require a subscription to the appropriate publication to read.

 - https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf
 - https://www.nature.com/articles/nature04795
 - https://www.pnas.org/content/105/12/4639.short

## Copyright and Licensing

The source code for SpatialSim is licensed under the GPLv3, see
[LICENSE.md](LICENSE.md).

It is Copyright Imperial College of Science, Technology and Medicine.

Additional contributions for open-sourcing made by Imperial College of
Science, Technology and Medicine, GitHub Inc, and John Carmack are copyright
the authors.

Licensing details for material from other projects may be found in
[NOTICE.md](NOTICE.md). In summary:

SpatialSim includes code modified from
[RANLIB](https://people.sc.fsu.edu/~jburkardt/c_src/ranlib/ranlib.html) which
is licensed under the LGPLv3.

Sample data in the repository has been derived from the following sources:

WorldPop (www.worldpop.org - School of Geography and Environmental Science,
University of Southampton; Department of Geography and Geosciences, University
of Louisville; Departement de Geographie, Universite de Namur) and Center for
International Earth Science Information Network (CIESIN), Columbia University
(2018). Global High Resolution Population Denominators Project - Funded by The
Bill and Melinda Gates Foundation (OPP1134076).
https://dx.doi.org/10.5258/SOTON/WP00647

WorldPop is licensed under the Creative Commons Attribution 4.0 International
License (CC BY 4.0).  The text of the license can be found at:
https://creativecommons.org/licenses/by/4.0/legalcode

## Contributing, Reporting Issues and Code of Conduct

Due to time pressure on the development team we are not currently accepting
contributions to this repository.  You are free to [fork it]
(https://github.com/mrc-ide/covid-19-spatial-sim/fork).

If you do find issues with the code please raise them in our [Issue
Tracker](https://github.com/mrc-ide/covid-19-spatial-sim/issues).

This repository has a code of conduct which is detailed in
[CODE_OF_CONDUCT.md](./CODE_OF_CONDUCT.md).  When raising an issue in this
repository you agree to abide by the [code of conduct](./CODE_OF_CONDUCT.md).
