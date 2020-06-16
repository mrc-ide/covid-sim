# Overview of the model

This is WIP. Anything wrong or missing? Please improve and open a PR!

## Basic conceptual structure

### Geographical space is divided into cells

This is a spatial model. We divide a country into [cells](./model-glossary.md#Cells)
and [microcells](./model-glossary.md#Microcells) (9x9 microcells to a cell)
which are geolocated.

### People live in households located in cells

[People](./model-glossary.md#people) are allocated according to population
density data (from input files) to cells. People have an age, and other
attributes. People's residence location does not change, but they interact with
people in other cells via places (see below) and via random social interactions
governed by a spatial kernel function.

Each person is modelled using the 'Person' structure. Persons are held in a 1-dimensional array : 'Hosts'.

### People are associated with different civil institutions

People are assigned to [places](./model-glossary.md#Places) (institutions such
as households, offices, schools etc.) that have a geographical location.
[Place groups](./model-glossary.md#Places) which divides places into compartments
(the intent here is that you're less likely to be infected by someone in the same
office but who works on a different floor).

Each place is modelled using the 'Place' structure. Places are held in a 2-dimensional array : 'Places'. 

The first index of Places[][] corresponds to the type of place, the second index picks a specific Place of that type. For example if the index of the 'Hospitals' place type were 5, then Places[5][2] would be the third Hospital (array indices start at 0)

People don't move. Instead the simulation employs spatial mixing probability
distributions (spatial kernels) that control the probability that people in cell
X will infect people in cell Y located in another spatial region.

### Infections spread between people

Infections may be initially seeded in different ways. The simplest way is to
seed according to population density (but seeds can be from specific places,
or randomly etc.)

`InfectSweep` is the main function where infections spread. It loops over people
and transmits infections by calculating a [FOI](./model-glossary.md#FOI)
(force of infection). Infection-spreading is divided into 3 transmission
mechanisms:

- household infections (e.g. between family members)
- place infections (e.g. at work)
- spatial infections (e.g. when travelling around)

Spatial infection models contacts between individuals which have a frequency
which depends upon the distance between home locations (to avoid literally
moving people around cells), modelled using a kernel function that weights
according to both spatial distance and population densities.

## Further information

For more information on the model and associated interventions, please see
the following papers:

 - [Ferguson, Neil, et al. "Report 9: Impact of non-pharmaceutical interventions (NPIs) to reduce COVID19 mortality and healthcare demand." (2020).](https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf)
 - [Ferguson, Neil M., et al. "Strategies for mitigating an influenza pandemic." Nature 442.7101 (2006): 448-452.](https://www.nature.com/articles/nature04795)
 - [Ferguson, Neil M., et al. "Strategies for containing an emerging influenza pandemic in Southeast Asia." Nature 437.7056 (2005): 209-214.](https://www.nature.com/articles/nature04017)
     - [Supplementary information for above](https://static-content.springer.com/esm/art%3A10.1038%2Fnature04017/MediaObjects/41586_2005_BFnature04017_MOESM1_ESM.pdf) - this document provides details of the calibration and spatial modelling
 - [Halloran, M. Elizabeth, et al. "Modeling targeted layered containment of an influenza pandemic in the United States." Proceedings of the National Academy of Sciences 105.12 (2008): 4639-4644.](https://www.pnas.org/content/105/12/4639.short)


Please note that some of the above may require a subscription.
