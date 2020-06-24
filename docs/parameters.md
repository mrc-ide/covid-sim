# Parameters used for COVID

## Place Closure

<table><tr><td colspan="6"><strong>[Place closure start time]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</strong></td><td>Real, Scalar</td>
  <td><strong>Example:</strong></td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Start time of place closure intervention - in days after [Day of year interventions start], if set, or after [Number of detected cases needed before outbreak alert triggered] threshold is first met (old report 9 calibration option)</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Place closure second start time]<strong></td><tr>
<tr><td><strong>Class</td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>100000</td></tr>
<tr><td><strong>Description:</td><td colspan="5">Start time of second period place closure intervention - in days after [Day of year interventions start], if set, or after [Number of detected cases needed before outbreak alert triggered] threshold is first met (old report 9 calibration option). High value means it won't happen. This mechanism for modelling multiple blocks of place closure has been superseded by a more flexible specification of arbitrary multiple periods of interventions, though this format still works</td></tr>
</table>
<br/>


<table><tr><td colspan="6"><strong>[Administrative unit divisor for place closure]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not used for COVID - allows school closure to be triggered by total illness related absenteeism in administrative units. Leave set to this value</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Place types to close for admin unit closure (0/1 array)]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer Vector, length [Number of types of places]</td>
  <td><strong>Example:</td><td>1  1  1  0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not used for COVID - allows school closure to be triggered by total illness related absenteeism in administrative units. Leave set to this value</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Cumulative proportion of place members needing to become sick for admin unit closure]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real Vector, length [Number of types of places]</td>
  <td><strong>Example:</td><td>1  1  1  0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not used for COVID - allows school closure to be triggered by total illness related absenteeism in administrative units. Leave set to this value</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of places in admin unit needing to pass threshold for place closure]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not used for COVID - allows school closure to be triggered by total illness related absenteeism in administrative units. Leave set to this value</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Delay to start place closure]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Delay in days betweeen when place closure policy starts or incidence threshold is met and closure actually occurs</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Duration of place closure]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>91</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Duration of place closure policy in days (unless switched off first via [Trigger incidence per cell for end of place closure])</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of places remaining open after closure by place type]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real Vector, length [Number of types of places]</td>
  <td><strong>Example:</td><td>0  0  0.25  1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Proportion of each place type which stay open during closure. Places staying open are selected randomly during the initialisation of each run</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative household contact rate after closure]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1.5</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of household contact rates for individuals affected by place closure</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative spatial contact rate after closure]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1.25</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of random spatial contact rates for individuals affected by place closure</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Place closure incidence threshold]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not used for COVID - allows school closure to be triggered based on incidence in the place itself. Needs to be kept at 0</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Place closure fractional incidence threshold]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not used for COVID - allows school closure to be triggered based on incidence in the place itself. Needs to be kept at 0</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Trigger incidence per cell for place closure]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>200</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">This parameter (which is poorly named) specifies the global or admin unit specific incidence threshold (e.g. of ICU cases) which triggers closure. See [Use global triggers for interventions] and [Use cases per thousand threshold for area controls]</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Trigger incidence per cell for end of place closure]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>50</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">This parameter (which is poorly named) specifies the global or admin unit specific incidence threshold (e.g. of ICU cases) which triggers the end of place closure. See [Use global triggers for interventions] and [Use cases per thousand threshold for area controls]</td></tr>
</table>
<br/>


## Household Quarantine

<table><tr><td colspan="6"><strong>[Household quarantine start time]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Start time of household quarantine intervention - in days after [Day of year interventions start], if set, or after [Number of detected cases needed before outbreak alert triggered] threshold is first met (old report 9 calibration option)</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Delay to start household quarantine]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Delay in days from when an index case develop symptoms and other household members enter quarantine</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Length of time households are quarantined]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>14</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Duration of household quarantine for an individual household</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Duration of household quarantine policy]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>91</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Duration HQ policy is in force</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative household contact rate after quarantine]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>2</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of household contact rates for a household in HQ</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Residual place contacts after household quarantine by place type]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Vector of length [Number of types of places]</td>
  <td><strong>Example:</td><td>0.25 0.25 0.25 0.25</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of place contact rates for someone in HQ</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Residual spatial contacts after household quarantine]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.25</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of random spatial contact rates for individuals affected by HQ</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Household level compliance with quarantine]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.5</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Household level compliance with HQ</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Individual level compliance with quarantine]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">For a complaint household, what proportion of its members comply with HQ</td></tr>
</table>
<br/>


## Case Isolation

<table><tr><td colspan="6"><strong>[Case isolation start time]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Start time of case isolation intervention - in days after [Day of year interventions start], if set, or after [Number of detected cases needed before outbreak alert triggered] threshold is first met (old report 9 calibration option)</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of detected cases isolated]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.7</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Proportion of symptomatic cases who self-isolate at home</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Delay to start case isolation]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Delay in days from symptom onset to self isolation at home</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Duration of case isolation]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>7</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Number of days a symptomatic case self isolates for</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Duration of case isolation policy]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>91</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Duration of the case isolation policy</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Residual contacts after case isolation]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.25</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Proportion of contacts still occuring after self-isolation</td></tr>
</table>
<br/>

## Social Distancing

<table><tr><td colspan="6"><strong>[Social distancing start time]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Start time of social distancing intervention - in days after [Day of year interventions start], if set, or after [Number of detected cases needed before outbreak alert triggered] threshold is first met (old report 9 calibration option)</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Trigger incidence per cell for social distancing]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>200</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">This parameter (which is poorly named) specifies the global or admin unit specific incidence threshold (e.g. of ICU cases) which triggers social distancing. See [Use global triggers for interventions] and [Use cases per thousand threshold for area controls]</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Trigger incidence per cell for end of social distancing]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>50</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">This parameter (which is poorly named) specifies the global or admin unit specific incidence threshold (e.g. of ICU cases) which triggers the end of social distancing. See [Use global triggers for interventions] and [Use cases per thousand threshold for area controls]</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Duration of social distancing]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>91</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Duration of the social distancing policy (unless switched off via [Trigger incidence per cell for end of social distancing])</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative place contact rate given social distancing by place type]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Vector of length [Number of types of places]</td>
  <td><strong>Example:</td><td>1  1  0.75  0.75</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of place contact rates for someone social distancing</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative household contact rate given social distancing]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1.25</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of household contact rates for someone who is social distancing</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative spatial contact rate given social distancing]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.25</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of random spatial contact rates for someone social distancing</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Minimum radius for social distancing]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not used for COVID - allows ssocial distancing to be triggered in rings around local microcells. Leave at 1.</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion compliant with enhanced social distancing]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Specifies a proportion of people complying with "enhanced" social distancing (typically specified has having higher effectiveness)</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion compliant with enhanced social distancing by age group]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Vector of length 17 (5 year age groups to 80+)</td>
  <td><strong>Example:</td><td>0  0  0  0  0  0  0  0  0  0  0  0  0  0  0.75  0.75  0.75</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">If [Age dependent severity delays] set to 1, otherwise 1 mean value for all ages. If this parameter is present, [Proportion compliant with enhanced social distancing] is ignored. This parameter allows compliance with enhanced social distancing to be specified by 5 year age bands, allowing sheilding to be modelled. Can be included in addition to standard social distancing - those complying with enhanced social distancing have those parameters applied, people not complying with enhanced social distancing have standard social distancing params applied</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative place contact rate given enhanced social distancing by place type]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Vector of length [Number of types of places]</td>
  <td><strong>Example:</td><td>0.25  0.25  0.25  0.25</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of place contact rates for someone complying with enhanced social distancing</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative household contact rate given enhanced social distancing]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of household contact rates for someone who complying with enhanced social distancing</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative spatial contact rate given enhanced social distancing]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.25</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Multiplicative scaling of random spatial contact rates for someone who complying with enhanced social distancing</td></tr>
</table>
<br/>

## False Positive

<table><tr><td colspan="6"><strong>[False positive rate]</strong></td><tr>
<tr><td><strong>Class</td><td>Transmission</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[False positive per capita incidence]</strong></td><tr>
<tr><td><strong>Class</td><td>Transmission</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[False positive relative incidence by age]</strong></td><tr>
<tr><td><strong>Class</td><td>Transmission</td>
  <td><strong>Type:</td><td>Real, Vector of length 17 (5 year age groups to 80+)</td>
  <td><strong>Example:</td><td>0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not currently used for COVID</td></tr>
</table>
<br/>


# Parameters not currently used in modelling COVID

## Place Closure

<table><tr><td colspan="6"><strong>[Place closure in administrative units rather than rings]<strong></td><tr>
<tr><td><strong>Class</td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer Scalar</td>
  <td><strong>Example:</td><td>100000</td></tr>
<tr><td><strong>Description:</td><td colspan="5">Not used for COVID - allows school closure to be triggered by total illness related absenteeism in administrative units. Leave set to this value</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Minimum radius for place closure]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Not used for COVID - allows school closure to be triggered in rings around local microcells. Leave at 1.</td></tr>
</table>
<br/>


## Treatment


<table><tr><td colspan="6"><strong>[Treatment start time]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>100000</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative susceptibility of treated individual]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.7</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative infectiousness of treated individual]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.4</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of symptomatic cases prevented by treatment]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.65</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of symptomatic cases resulting in death prevented by treatment]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Delay to treat cell]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Duration of course of treatment]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>5</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Duration of course of prophylaxis]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>10</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion treated]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.9</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Treatment radius]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of households of cases treated]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of places treated after case detected]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Vector of length [Number of types of places]</td>
  <td><strong>Example:</td><td>0  0  0  0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of people treated in targeted places]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Vector of length [Number of types of places]</td>
  <td><strong>Example:</td><td>0  0  0  0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Only treat mixing groups within places]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Maximum number of doses available]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer, Scalar</td>
  <td><strong>Example:</td><td>60000000</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models antiviral treatment - not currently for COVID</td></tr>
</table>
<br/>

## Movement Restrictions

<table><tr><td colspan="6"><strong>[Movement restrictions start time]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>100000</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models movement restrictions - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Movement restrictions trigger incidence per cell]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer, Scalar</td>
  <td><strong>Example:</td><td>1000000000</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models movement restrictions - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Delay to start movement restrictions]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>2</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models movement restrictions - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Duration of movement restrictions]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>7</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models movement restrictions - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Residual movements after restrictions]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.25</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models movement restrictions - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Minimum radius of movement restrictions]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>5000</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models movement restrictions - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Impose blanket movement restrictions]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models movement restrictions - not currently used for COVID</td></tr>
</table>
<br/>

## Vaccination

<table><tr><td colspan="6"><strong>[Vaccination start time]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>100000</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Duration of household vaccination policy]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1000</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Apply mass rather than reactive vaccination]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Priority age range for mass vaccination]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer, Vector of length 2 (from, to)</td>
  <td><strong>Example:</td><td>50  85</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Switch time at which efficacy increases]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1000</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative susceptibility of vaccinated individual after switch time]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.2</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative susceptibility of vaccinated individual]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.5</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Relative infectiousness of vaccinated individual]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.5</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of symptomatic cases prevented by vaccination]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.75</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Vaccination trigger incidence per cell]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Delay to vaccinate]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Delay from vaccination to full protection]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>14</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of population vaccinated]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0.9</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Vaccination radius]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>1</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Proportion of households of cases vaccinated]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Maximum number of vaccine courses available]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer, Scalar</td>
  <td><strong>Example:</td><td>1000000000</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Start time of additional vaccine production]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Real, Scalar</td>
  <td><strong>Example:</td><td>10000</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>

<table><tr><td colspan="6"><strong>[Rate of additional vaccine production (courses per day)]</strong></td><tr>
<tr><td><strong>Class</strong></td><td>Interventions</td>
  <td><strong>Type:</td><td>Integer, Scalar</td>
  <td><strong>Example:</td><td>0</td></tr>
<tr><td><strong>Description:</strong></td><td colspan="5">Models vaccination - not currently used for COVID</td></tr>
</table>
<br/>
