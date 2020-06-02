# CovidSim: how interventions and policies are implemented

See [model overview](./model-overview.md) for further information.

CovidSim models infections through time in a spatially explicit representation
of individuals in a population. Infected people are assumed to infect
susceptible people at household, spatial and place level. Here,
[place](./model-glossary.md#Places) is a generalisation of either primary
schools, secondary schools, universities or workplaces (e.g. offices). The model
includes disease progression for symptomatic individuals, and is parameterised
to reflect current knowledge of Covid-19 disease progression and death.

Interventions are implemented in the code through the parameter and
pre-parameter files (see [input and outputs](./inputs-and-outputs.md)). These
parameters are read in function `ReadParams` in `CovidSim.cpp` and collected in
single structure `PARAMS`, of which there is a single, global instance `P`.

## Intervention overview

The model considers various interventions that alter either the infectiousness
of infected people or the susceptibility of susceptible people (or both).
Fundamentally, interventions are implemented by changing the person-specific
period during which a person is subject to a particular intervention (given in
`#define` functions in `ModelMacros.h`). These values in turn toggle parameters
in functions governing infectiousness and susceptibility in `CalcInfSusc.cpp`,
and in function `InfectSweep` in `Update.cpp`. The interventions broadly fall
into the following categories:

- Case Isolation (CI),
- Household Quarantine (HQ),
- Place Closure (PC),
- Social Distancing (SD).

Example parameter files have been provided for the following scenarios:

1. no interventions (`p_NoInt.txt`)
2. case isolation; household quarantine; social
distancing and place closure (`p_PC7_CI_HQ_SD.txt`).

Units of time `t` are in days, and day 1 is taken to be January 1st 2020 (e.g.
day 76 is March 16th 2020) in output files. The model is calibrated at the start
of runtime to match the user-specified inputs. It will approximately (through an
iterative calibration process) match a user-specified threshold cumulative
epidemic size (either [Number of deaths accumulated before alert] for deaths or
[Number of detected cases needed before outbreak alert triggered] for detected
cases, depending on value of `P.PreControlClusterIdUseDeaths` [Trigger alert on
deaths]) occurring up to the day specified by [Day of year trigger is reached]
(`P.PreControlClusterIdCalTime`). If [Alert trigger starts after interventions]
is set to 1 in the input pre-parameter file, this threshold can occur after the
start of interventions, in which case the start day of interventions,
[Day of year interventions start] (`P.PreIntervIdCalTime`) should be also
specified.

Intervention times in the intervention parameter files are relative to
`P.PreIntervIdCalTime`.
(e.g. if `P.HQuarantineTimeStartBase` / [Household quarantine start time] is set
to 10, then in the start date of household quarantine policy will be ten days
after `P.PreIntervIdCalTime`).

Importantly, the model distinguishes between the period during which a policy is
implemented, and the period during which the actual intervention is implemented.
If a policy is implemented from time `t`, the intervention will not be
implemented until a user-specified threshold has been met by a particular
“trigger” or count (see below). Further, CI and HQ policies have durations at
the level of a person, which denote the length of time that a person would
either self-isolate or be under home quarantine. These person-level durations
are distinct from the duration of policy at the population level. For example,
a policy could be implemented for a year, during which time people would home
quarantine for 14 days if a member of their household becomes sick.

### Intervention triggers

A trigger is defined flexibly, with reference to either the numbers (cumulative
incidence) of confirmed cases, numbers of critical cases that require Intensive
Care Units (ICU beds), or numbers of deaths. Numbers can be either absolute or
per-capita, and the model may be refined so that the number is defined in
relation to ICU bed capacity. The window over which cumulative incidence is
calculated is specified by `P.TriggersSamplingInterval` [Number of sampling
intervals over which cumulative incidence measured for global trigger] days.
Further, triggers may be specified at either admin unit (regional) or at country
(global) level. So too can the start times for when each intervention is
considered. Finally, each intervention has its own threshold value.

Currently, all interventions must use the same epidemiological variable for
their trigger (e.g. all using confirmed cases or all using ICU beds), and
interventions do not differ in whether they consider absolute or per-capita
triggers. Note that triggers can be set to zero, indicating that interventions
are implemented throughout the time window specified by the start time and
duration of that intervention.

The steps below summarise the above process:

- Simulation starts
- Policy implemented
- Trigger (count) surpasses threshold for that intervention
- Intervention implemented
- infectiousness and/or susceptibility parameters turned on for individuals
  affected by intervention

The following interventions are described with reference to particular functions
and parameters in code, and the parameter descriptions in parameter files
(`p_XX.txt`), which are read in through `ReadParams` function in `CovidSim.cpp`.
Functions `TreatSweep` in `Update.cpp` and `RecordSample` in `CovidSim.cpp`
determine whether thresholds have been exceeded.

Notes on variable naming: Parameters with prefix `us` are unsigned short
versions of their namesakes, that have been multiplied by parameter `P.TimeStepsPerDay`
(e.g. `P.usHQuarantineHouseDuration = P.HQuarantineHouseDuration * P.TimeStepsPerDay`).

## CASE ISOLATION (CI)

Policy implemented after time `P.CaseIsolationTimeStartBase`
(“Case isolation start time”) for duration of `P.CaseIsolationPolicyDuration`
days. A value of `P.CaseIsolationTimeStartBase` greater than simulation duration
(`P.SampleTime` “Sampling time” in pre-parameter file) indicates that policy
will not be implemented (e.g. if simulation runs for 720 days and start time set
to 1000 days). If threshold `P.CaseIsolation_CellIncThresh` exceeded, symptomatic
cases stay home for period of `P.CaseIsolationDuration` days, starting after
`P.CaseIsolationDelay` days. Durations and delays can alternatively be specified
at admin unit level.

After infection and upon case detection (`DoDetectCase` function in `Update.cpp`),
a case isolates with probability `P.CaseIsolationProp`.
Their isolation_start_time is set. Isolation stop time formula:
`isolation_start_time + P.usCaseIsolationDelay + P.usCaseIsolationDuration`).
These values are relevant to `HOST_ISOLATED` macro, which returns true or false.
`HOST_ISOLATED` affects every infectiousness function (at person, house, place,
spatial levels) but not susceptibility functions in `CalcInfSusc.cpp`.

If `HOST_ISOLATED`,
place and spatial infectiousness scaled by `P.CaseIsolationEffectiveness`,
whereas household infectiousness scaled by `P.CaseIsolationHouseEffectiveness`.

Additionally, case isolation amends absent start and stop time, in scenarios
where isolated person is a child requiring adult supervision, and where host was
already absent and requires their period of absence amended due to case isolation.
Macro `HOST_ABSENT` relevant in `InfectSweep` in `Sweep.cpp` – all place
infections are contingent on host not being absent from their place.

## HOUSEHOLD/HOME QUARANTINE (HQ)

Policy implemented after time `P.HQuarantineTimeStartBase` (“Household quarantine
start time”) for period of `P.HQuarantinePolicyDuration` days.
Value of `P.HQuarantineTimeStartBase` greater than simulation duration indicates
that policy will not be implemented (e.g. if simulation runs for 720 days and
start time set to 1000 days). If threshold `P.HHQuar_CellIncThresh` exceeded,
people who share household with a symptomatic case stay home for period of
`P.HQuarantineHouseDuration days`, starting after `P.HQuarantineHouseDelay days`.
Durations and delays can alternatively be specified at admin unit level.

Upon case detection (`DoDetectCase` function in `Update.cpp`), `quar_start_time` is set.
A person’s quarantine stop time is `quar_start_time + P.usHQuarantineHouseDuration`.
Compliance is set first at household and then at individual level, both of which
are randomly determined with respective proportions `P.HQuarantinePropHouseCompliant`
and `P.HQuarantinePropIndivCompliant`. Note that 50% compliance means that 50%
of people are 100% compliant, and not that 100% of people are 50% compliant
(although the latter probably more realistic).

Macro `HOST_QUARANTINED` asks whether host is a) compliant with HQ; b) currently
between their quarantine start and stop times. `HOST_QUARANTINED` scales
infectiousness: at household level by `P.HQuarantineHouseEffect`; at spatial
level by `P.HQuarantineSpatialEffect`; and at place level by
`P.HQuarantinePlaceEffect[k]` for place type `k`. `HOST_QUARANTINED` scales
susceptibility: at spatial level by `P.HQuarantineSpatialEffect`; and at place
level by `P. HQuarantinePlaceEffect[k]` for place type `k`.

## SOCIAL DISTANCING (SD)

Policy implemented after time `P.SocDistTimeStartBase` (“Social distancing start
time”) for period of `P.SocDistDuration` days. Value of `P.SocDistTimeStartBase`
greater than simulation duration indicates that policy will not be implemented
(e.g. if simulation runs for 720 days and start time set to 1000 days).
Durations and delays can alternatively be specified at admin unit level.

Social distancing implemented mainly in `TreatSweep`, which determines whether
`P.SocDistCellIncThresh` incidence threshold exceed. If so, then the
susceptibility of all people residing in microcell currently undergoing social
distancing is scaled, either by enhanced or non-enhanced effects (see below).
Infectiousness functions are unaffected by social distancing.

Distinction between social distancing and enhanced social distancing effect:

- Household (`P.SocDistHouseholdEffect` vs. `P.EnhancedSocDistHouseholdEffect`)
- Place (`P.SocDistPlaceEffect[k]` vs. `P.EnhancedSocDistPlaceEffect[k]`, for place type `k`)
- Spatial (`P.SocDistSpatialEffect` vs. `P.EnhancedSocDistSpatialEffect`)

These scale a person’s contact rate at household, place and spatial level.

Both sets of parameters can implemented at the same time, in which case a given
person is compliant with enhanced social distancing (`person.esocdist_comply == 1`)
with age-specific probability `P.EnhancedSocProportionCompliant[age]`
(decided for each person in `SetupModel.cpp`).

A social distancing policy can occur once only (if `P.DoSocDistOnceOnly`), or repeatedly.

## SOCIAL DISTANCING FOR OVER 70s (SDO/SDOL70)

Differences between this and social distancing are in parameter values, rather than source code.

`ReadParams` in `CovidSim.cpp` populates `P.ESocProportionCompliant` array by
scanning first for “Proportion compliant with enhanced social distancing by
age group”, and failing that the non-age-specific “Proportion compliant with
enhanced social distancing”, in which case all values of array same for all age groups.
If this is not found, values default to 0, i.e. zero compliance with enhanced social distancing.

## PLACE CLOSURE (PC)

Policy implemented after time `P.PlaceCloseTimeStartBase` (“Place closure start time”).
Values of `P.PlaceCloseTimeStartBase` greater than simulation duration indicates
that policy will not be implemented (e.g. if simulation runs for 720 days and
start time set to 1000 days). Durations can alternatively be specified at admin
unit level.

`RecordSample` in `CovidSim.cpp` determines whether trigger has exceeded
`P.PlaceCloseCellIncThresh`, and if so sets `P.PlaceCloseTimeStart`.
A subsequent call to `TreatSweep` function in `Sweep.cpp`, determines whether
incidence threshold `P.PlaceCloseCellIncThresh` exceeded. If so, then all places
in microcell marked as closed (and trigger reset), and `DoPlaceClose` in
`Update.cpp` is called. DoPlaceClose changes quantities `place.close_start_time`
and `place.close_stop_time`, which are used in macro function `PLACE_CLOSED`.

Unlike other measures, place closure does not affect infectiousness and
susceptibility functions in `CalcInfSusc.cpp`, but it does affect infectiousness
and susceptibility indirectly in `InfectSweep` function in `Update.cpp` with
calls to `PLACE_CLOSED`. If person has their places closed, their household
infectiousness is scaled by `P.PlaceCloseHouseholdRelContact`. Similarly, infected
people’s spatial infectiousness also scaled by `P.PlaceCloseSpatialRelContact`,
and potential infectee’s susceptibility to spatial infections is also scaled by
`P.PlaceCloseSpatialRelContact`. Further, infected people cannot infect anyone
in their place (and therefore have zero Place infectiousness) when their place
is closed.

Like social distancing, place closure can occur intermittently if parameter
`P.DoPlaceCloseOnceOnly` set to 0.
