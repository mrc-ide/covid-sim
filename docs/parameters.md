## Parameters used for COVID

### Place Closure Related

<table>
  <tr>
    <td colspan="2">[Place closure start time]</td>
  <tr>
</table>

Interventions	Real	1	0	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Start time of place closure intervention - in days after [Day of year interventions start], if set, or after [Number of detected cases needed before outbreak alert triggered] threshold is first met (old report 9 calibration option)

[Place closure second start time]	Interventions	Real	1	100000	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Start time of second period place closure intervention - in days after [Day of year interventions start], if set, or after [Number of detected cases needed before outbreak alert triggered] threshold is first met (old report 9 calibration option). High value means it won't happen. This mechanism for modelling multiple blocks of place closure has been superseded by a more flexible specification of arbitrary multiple periods of interventions, though this format still works

[Place closure in administrative units rather than rings]	Interventions	Integer	1	0	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Not used for COVID - allows school closure to be triggered by total illness related absenteeism in administrative units. Leave set to this value

[Administrative unit divisor for place closure]	Interventions	Integer	1	1	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Not used for COVID - allows school closure to be triggered by total illness related absenteeism in administrative units. Leave set to this value

[Place types to close for admin unit closure (0/1 array)]	Interventions	Integer	[Number of types of places]	1  1  1  0	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Not used for COVID - allows school closure to be triggered by total illness related absenteeism in administrative units. Leave set to this value

[Cumulative proportion of place members needing to become sick for admin unit closure]	Interventions	Real	[Number of types of places]	1  1  1  0	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Not used for COVID - allows school closure to be triggered by total illness related absenteeism in administrative units. Leave set to this value

[Proportion of places in admin unit needing to pass threshold for place closure]	Interventions	Real	1	1	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Not used for COVID - allows school closure to be triggered by total illness related absenteeism in administrative units. Leave set to this value

[Delay to start place closure]	Interventions	Real	1	1	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Delay in days betweeen when place closure policy starts or incidence threshold is met and closure actually occurs

[Duration of place closure]	Interventions	Real	1	91	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Duration of place closure policy in days (unless switched off first via [Trigger incidence per cell for end of place closure])

[Proportion of places remaining open after closure by place type]	Interventions	Real	[Number of types of places]	0  0  0.25  1	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Proportion of each place type which stay open during closure. Places staying open are selected randomly during the initialisation of each run

[Relative household contact rate after closure]	Interventions	Real	1	1.5	Could be varied by +/- 0.5	Multiplicative scaling of household contact rates for individuals affected by place closure

[Relative spatial contact rate after closure]	Interventions	Real	1	1.25	Could be varied by +/- 0.25	Multiplicative scaling of random spatial contact rates for individuals affected by place closure

[Minimum radius for place closure]	Interventions	Real	1	1	No - not used	Not used for COVID - allows school closure to be triggered in rings around local microcells. Leave at 1.

[Place closure incidence threshold]	Interventions	Real	1	0	No - not used	Not used for COVID - allows school closure to be triggered based on incidence in the place itself. Needs to be kept at 0

[Place closure fractional incidence threshold]	Interventions	Real	1	0	No - not used	Not used for COVID - allows school closure to be triggered based on incidence in the place itself. Needs to be kept at 0

[Trigger incidence per cell for place closure]	Interventions	Real	1	200	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	This parameter (which is poorly named) specifies the global or admin unit specific incidence threshold (e.g. of ICU cases) which triggers closure. See [Use global triggers for interventions] and [Use cases per thousand threshold for area controls]

[Trigger incidence per cell for end of place closure]	Interventions	Real	1	50	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	This parameter (which is poorly named) specifies the global or admin unit specific incidence threshold (e.g. of ICU cases) which triggers the end of place closure. See [Use global triggers for interventions] and [Use cases per thousand threshold for area controls]
						
### Household Quarantine

[Household quarantine start time]	Interventions	Real	1	0	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Start time of household quarantine intervention - in days after [Day of year interventions start], if set, or after [Number of detected cases needed before outbreak alert triggered] threshold is first met (old report 9 calibration option)

[Delay to start household quarantine]	Interventions	Real	1	1	Could be varied by +/- 1	Delay in days from when an index case develop symptoms and other household members enter quarantine

[Length of time households are quarantined]	Interventions	Real	1	14	Could be varied by +/- 7	Duration of household quarantine for an individual household

[Duration of household quarantine policy]	Interventions	Real	1	91	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Duration HQ policy is in force

[Relative household contact rate after quarantine]	Interventions	Real	1	2	Could be varied between 1 and 2 (2 is already quite pessimistic)	Multiplicative scaling of household contact rates for a household in HQ

[Residual place contacts after household quarantine by place type]	Interventions	Real	[Number of types of places]	0.25  0.25  0.25  0.25	Could be varied by up to +/-0.25 for the last two entries (universities and workplaces)	Multiplicative scaling of place contact rates for someone in HQ

[Residual spatial contacts after household quarantine]	Interventions	Real	1	0.25	Could be varied by +/- 0.15	Multiplicative scaling of random spatial contact rates for individuals affected by HQ

[Household level compliance with quarantine]	Interventions	Real	1	0.5	Could be increased to 0.75 (0.5 already pessimistic)	Household level compliance with HQ

[Individual level compliance with quarantine]	Interventions	Real	1	1	Could be decreased to 0.75, though contact rate scaling parameters already imply <100% compliance	For a complaint household, what proportion of its members comply with HQ
						

### Case Isolation

[Case isolation start time]	Interventions	Real	1	0	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Start time of case isolation intervention - in days after [Day of year interventions start], if set, or after [Number of detected cases needed before outbreak alert triggered] threshold is first met (old report 9 calibration option)

[Proportion of detected cases isolated]	Interventions	Real	1	0.7	Could be varied by +/-0.2	Proportion of symptomatic cases who self-isolate at home

[Delay to start case isolation]	Interventions	Real	1	1	Could be varied by +/-1	Delay in days from symptom onset to self isolation at home

[Duration of case isolation]	Interventions	Real	1	7	Could be varied by +/-2	Number of days a symptomatic case self isolates for

[Duration of case isolation policy]	Interventions	Real	1	91	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Duration of the case isolation policy

[Residual contacts after case isolation]	Interventions	Real	1	0.25	could be varied by +/-0.2	Proportion of contacts still occuring after self-isolation
						
### Social Distancing

[Social distancing start time]	Interventions	Real	1	0	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Start time of social distancing intervention - in days after [Day of year interventions start], if set, or after [Number of detected cases needed before outbreak alert triggered] threshold is first met (old report 9 calibration option)

[Trigger incidence per cell for social distancing]	Interventions	Real	1	200	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	This parameter (which is poorly named) specifies the global or admin unit specific incidence threshold (e.g. of ICU cases) which triggers social distancing. See [Use global triggers for interventions] and [Use cases per thousand threshold for area controls]

[Trigger incidence per cell for end of social distancing]	Interventions	Real	1	50	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	This parameter (which is poorly named) specifies the global or admin unit specific incidence threshold (e.g. of ICU cases) which triggers the end of social distancing. See [Use global triggers for interventions] and [Use cases per thousand threshold for area controls]

[Duration of social distancing]	Interventions	Real	1	91	Not suited to parametric sensitivity analysis - intended to allow scenario analysis of different intervention strategies	Duration of the social distancing policy (unless switched off via [Trigger incidence per cell for end of social distancing])

[Relative place contact rate given social distancing by place type]	Interventions	Real	[Number of types of places]	1  1  0.75  0.75	Could be varied by up to +/-0.25 for the last two entries (universities and workplaces)	Multiplicative scaling of place contact rates for someone social distancing

[Relative household contact rate given social distancing]	Interventions	Real	1	1.25	Could be varied between 1 and 2 	Multiplicative scaling of household contact rates for someone who is social distancing

[Relative spatial contact rate given social distancing]	Interventions	Real	1	0.25	Could be varied by +/- 0.15	Multiplicative scaling of random spatial contact rates for someone social distancing

[Minimum radius for social distancing]	Interventions	Real	1	1	No - not used	Not used for COVID - allows ssocial distancing to be triggered in rings around local microcells. Leave at 1.

[Proportion compliant with enhanced social distancing]	Interventions	Real	1	0	No - use [Proportion compliant with enhanced social distancing by age group]	Specifies a proportion of people complying with "enhanced" social distancing (typically specified has having higher effectiveness)

[Proportion compliant with enhanced social distancing by age group]	Interventions	Real	17 (5 year age groups up to 80+) if [Age dependent severity delays] set to 1, otherwise 1 mean value for all ages	0  0  0  0  0  0  0  0  0  0  0  0  0  0  0.75  0.75  0.75	Non zero entries could be varied by +/-0.25	If this parameter is present, [Proportion compliant with enhanced social distancing] is ignored. This parameter allows compliance with enhanced social distancing to be specified by 5 year age bands, allowing sheilding to be modelled. Can be included in addition to standard social distancing - those complying with enhanced social distancing have those parameters applied, people not complying with enhanced social distancing have standard social distancing params applied

[Relative place contact rate given enhanced social distancing by place type]	Interventions	Real	[Number of types of places]	0.25  0.25  0.25  0.25	Could be varied by +/- 0.15	Multiplicative scaling of place contact rates for someone complying with enhanced social distancing

[Relative household contact rate given enhanced social distancing]	Interventions	Real	1	1	Could be varied between 1 and 2 	Multiplicative scaling of household contact rates for someone who complying with enhanced social distancing

[Relative spatial contact rate given enhanced social distancing]	Interventions	Real	1	0.25	Could be varied by +/- 0.15	Multiplicative scaling of random spatial contact rates for someone who complying with enhanced social distancing
						
### False Positive

[False positive rate]	Transmission	Real	1	0	No	Not currently used for COVID

[False positive per capita incidence]	Transmission	Real	1	0	No	Not currently used for COVID

[False positive relative incidence by age]	Transmission	Real	2	0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0	No	Not currently used for COVID


## Parameters not currently used in modelling COVID

### Treatment

Parameter name	Class	Type	Length (1=scalar, >1 = vector)	Example value	Sensitivity analysis	Description
[Treatment start time]	Interventions	Real	1	100000	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Relative susceptibility of treated individual]	Interventions	Real	1	0.7	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Relative infectiousness of treated individual]	Interventions	Real	1	0.4	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Proportion of symptomatic cases prevented by treatment]	Interventions	Real	1	0.65	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Proportion of symptomatic cases resulting in death prevented by treatment]	Interventions	Real	1	0	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Delay to treat cell]	Interventions	Real	1	1	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Duration of course of treatment]	Interventions	Real	1	5	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Duration of course of prophylaxis]	Interventions	Real	1	10	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Proportion treated]	Interventions	Real	1	0.9	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Treatment radius]	Interventions	Real	1	0	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Proportion of households of cases treated]	Interventions	Real	1	0	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Proportion of places treated after case detected]	Interventions	Real	[Number of types of places]	0 0 0 0	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Proportion of people treated in targeted places]	Interventions	Real	[Number of types of places]	0 0 0 0	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Only treat mixing groups within places]	Interventions	Integer	1	0	No - not used for COVID	Models antiviral treatment - not currently for COVID
[Maximum number of doses available]	Interventions	Integer	1	60000000	No - not used for COVID	Models antiviral treatment - not currently for COVID
						
### Movement Restrictions

[Movement restrictions start time]	Interventions	Real	1	100000	No - not used for COVID	Models movement restrictions - not currently used for COVID
[Movement restrictions trigger incidence per cell]	Interventions	Integer	1	1000000000	No - not used for COVID	Models movement restrictions - not currently used for COVID
[Delay to start movement restrictions]	Interventions	Real	1	2	No - not used for COVID	Models movement restrictions - not currently used for COVID
[Duration of movement restrictions]	Interventions	Real	1	7	No - not used for COVID	Models movement restrictions - not currently used for COVID
[Residual movements after restrictions]	Interventions	Real	1	0.25	No - not used for COVID	Models movement restrictions - not currently used for COVID
[Minimum radius of movement restrictions]	Interventions	Real	1	5000	No - not used for COVID	Models movement restrictions - not currently used for COVID
[Impose blanket movement restrictions]	Interventions	Integer	1	0	No - not used for COVID	Models movement restrictions - not currently used for COVID
						
### Vaccination

[Vaccination start time]	Interventions	Real	1	100000	No - not used for COVID	Models vaccination - not currently used for COVID
[Duration of household vaccination policy]	Interventions	Real	1	1000	No - not used for COVID	Models vaccination - not currently used for COVID
[Apply mass rather than reactive vaccination]	Interventions	Integer	1	1	No - not used for COVID	Models vaccination - not currently used for COVID
[Priority age range for mass vaccination]	Interventions	Integer	2	50  85  	No - not used for COVID	Models vaccination - not currently used for COVID
[Switch time at which efficacy increases]	Interventions	Real	1	1000	No - not used for COVID	Models vaccination - not currently used for COVID
[Relative susceptibility of vaccinated individual after switch time]	Interventions	Real	1	0.2	No - not used for COVID	Models vaccination - not currently used for COVID
[Relative susceptibility of vaccinated individual]	Interventions	Real	1	0.5	No - not used for COVID	Models vaccination - not currently used for COVID
[Relative infectiousness of vaccinated individual]	Interventions	Real	1	0.5	No - not used for COVID	Models vaccination - not currently used for COVID
[Proportion of symptomatic cases prevented by vaccination]	Interventions	Real	1	0.75	No - not used for COVID	Models vaccination - not currently used for COVID
[Vaccination trigger incidence per cell]	Interventions	Integer	1	0	No - not used for COVID	Models vaccination - not currently used for COVID
[Delay to vaccinate]	Interventions	Real	1	0	No - not used for COVID	Models vaccination - not currently used for COVID
[Delay from vaccination to full protection]	Interventions	Real	1	14	No - not used for COVID	Models vaccination - not currently used for COVID
[Proportion of population vaccinated]	Interventions	Real	1	0.9	No - not used for COVID	Models vaccination - not currently used for COVID
[Vaccination radius]	Interventions	Real	1	1	No - not used for COVID	Models vaccination - not currently used for COVID
[Proportion of households of cases vaccinated]	Interventions	Real	1	0	No - not used for COVID	Models vaccination - not currently used for COVID
[Maximum number of vaccine courses available]	Interventions	Integer	1	1000000000	No - not used for COVID	Models vaccination - not currently used for COVID
[Start time of additional vaccine production]	Interventions	Real	1	10000	No - not used for COVID	Models vaccination - not currently used for COVID
[Rate of additional vaccine production (courses per day)]	Interventions	Integer	1	0	No - not used for COVID	Models vaccination - not currently used for COVID
