# COVID-19 CovidSim Model

This is the COVID-19 CovidSim microsimulation model developed by the MRC Centre
for Global Infectious Disease Analysis hosted at Imperial College, London.

CovidSim models the transmission dynamics and severity of COVID-19 infections
throughout a spatially and socially structured population over time.  It enables
modelling of how intervention policies and healthcare provision affect the
spread of COVID-19. With parameter changes, it can be used to model other
respiratory viruses, such as influenza.

## IMPORTANT NOTES

:warning: This code is released with no support.

:warning: This model is in active development and so parameter name and
behaviours, and output file formats will change without notice.

:warning: The model is stochastic. Multiple runs with different seeds should be
undertaken to see average behaviour. This can now be done easily with the `/NR`
command line parameter. The model code behaves deterministically if run with the
same number of threads enabled and run with the same random number seends.

:warning: As with any mathematical model, it is easy to misconfigure inputs and
therefore get meaningless outputs. The Imperial College COVID-19 team only
endorses outputs it has itself generated.

## Status

This model is in active development and subject to significant code changes
to:

- Enable modelling of more geographies

- Enable modelling of different intervention scenarios

- Improve performance

## Building

The model is written in C++ and runs on Windows and Linux-based systems.

Running the model for the whole of the UK requires approximately 20GB of RAM.
Other regions will require different amounts of memory (some up to 256GB).

It is strongly recommended to build the model with OpenMP support enabled to
improve performance on multi-core processors. 24 to 32 core Xeon systems give
optimal performance for large (e.g. UK, US) populations.

See [build.md](./docs/build.md) for detailed build instructions.

### Testing

From within your build directory do:

```sh
make test
# If you want more progress indication
make test ARGS="-V"
# or
ctest -V
```

*IMPORTANT*: The test scripts use test data only are not runs reflective of
real-world situations.

## Sample Data

The directory [data](./data) contains sample data.

The Python script [run_sample.py](./data/run_sample.py) demonstrates how to
invoke CovidSim to use this data.  See the [sample README](./data/README.md) for
details on how to run the samples.

## Documentation

Model documentation can be found in the [docs](./docs) directory.  Of
particular interest are:

- [Model Overview](./docs/model-overview.md)
- [Model Glossary](./docs/model-glossary.md)
- [Model Inputs and Outputs](./docs/inputs-and-outputs.md)
- [Interventions description](./docs/intervention-description.md)
- [R Scripts for Output Visualisation](./docs/inputs-and-outputs.md#r-summary-visualisations)

Given the entire Imperial College team is working full-time on the COVID-19
response, documentation is currently sparse. More documentation and sample files
will be added as time permits. In the coming few weeks this will include a much
more extensive set of input files to model strategies for exiting lockdown.

### Relevant papers

The following papers are relevant to the model.  Please note that some of them
may require a subscription.

- <https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf>
- <https://www.nature.com/articles/nature04795>
- <https://www.nature.com/articles/nature04017>
- <https://www.pnas.org/content/105/12/4639.short>

## Copyright and Licensing

The source code for CovidSim is licensed under the GPLv3, see
[LICENSE.md](LICENSE.md).

It is Copyright Imperial College of Science, Technology and Medicine. The
lead developers are Neil Ferguson, Gemma Nedjati-Giliani and Daniel Laydon.

Additional contributions for open-sourcing made by Imperial College of
Science, Technology and Medicine, GitHub Inc, and John Carmack are copyright
the authors.

Licensing details for material from other projects may be found in
[NOTICE.md](NOTICE.md). In summary:

CovidSim includes code modified from
[RANLIB](https://people.sc.fsu.edu/~jburkardt/c_src/ranlib/ranlib.html) which
is licensed under the LGPLv3.

Sample data in the repository has been derived from the following sources:

WorldPop (www.worldpop.org - School of Geography and Environmental Science,
University of Southampton; Department of Geography and Geosciences, University
of Louisville; Departement de Geographie, Universite de Namur) and Center for
International Earth Science Information Network (CIESIN), Columbia University
(2018). Global High Resolution Population Denominators Project - Funded by The
Bill and Melinda Gates Foundation (OPP1134076).
<https://dx.doi.org/10.5258/SOTON/WP00647>

WorldPop is licensed under the Creative Commons Attribution 4.0 International
License (CC BY 4.0).  The text of the license can be found at:
<https://creativecommons.org/licenses/by/4.0/legalcode>

## Contributing

Due to time pressure on the development team, we are unable to provide user
support at this time.

If you find issues with the code please raise them in our
[Issue Tracker](https://github.com/mrc-ide/covid-sim/issues).

This repository has a code of conduct which is detailed in
[the code of conduct](./CODE_OF_CONDUCT.md).  When raising an issue in this
repository you agree to abide by the code of conduct.
