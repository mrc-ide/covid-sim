# Overview of the model

This is WIP. Anything wrong or missing? Please improve and open a PR!

# Basic conceptual structure

## Geographical space is divided into cells

This is a spatial model. We divide a country into [cells](./model-glossary.md#Cells) and [microcells](./model-glossary.md#Microcells) (9 microcells to a cell). Cells have a location on a geographical map.

## People live in cells

[People](./model-glossary.md#people) are allocated according to population density data (from input files) to cells. People have an age, and other attributes. People do not move between cells.

## People are associated with different civil institutions

People are assigned to [places](./model-glossary.md#Places) (institutions such as households, offices, schools etc.) that have a geographical location. [Place groups](./model-glossary.md#Places) which divides places into compartments (the intent here is that you're less likely to be infected by someone in the same office but who works on a different floor).

People don't move. Instead the simulation employs spatial mixing probability distributions (spatial kernels) that control the probability that people in cell X will infect people in cell Y located in another spatial region.

## Infections spread between people

Infections may be initially seeded in different ways. The simplest way is to seed according to population density (but seeds can be from specific places, or randomly etc.)

`InfectSweep` is the main function where infections spread. It loops over people and transmits infections by calculating a [FOI](./model-glossary.md#FOI) (force of infection). Infection-spreading is divided into 3 transmission mechanisms:
- household infections (e.g. between family members)
- place infections (e.g. at work)
- spatial infections (e.g. when travelling around)

Spatial infection uses spatial kernels to simulate infections across space (to avoid literally moving people around cells). A given cell accumulates a "spatial infectiousness" from the people in it (aggregate summary). This can then spread to other cells (as a kernel function that weights according to both spatial distance and population densities).
