# Parameters

WIP - interpretations below are based on reading the SetupModel.cpp file

The parameters class contains:
- configurable options for running the simulation 
- quantitative information used in the simulation
- data produced about the setup process e.g. running totals
- variables used to load/save information
- variable used to manage processor / memory


## Configurable Options
<details>
  <summary>Boolean - on/off options</summary>
  
  A value of '1' for the given variable turns that option 'on'
  
  | Variable | Description |
  | --------------- | --------------- |
  | DoBin |  | 
  | DoHeteroDensity | read population density info from file | 
  | DoAdUnits | use administrative units | 
  | DoAdUnitBoundaries | | 
  | DoPlaces | turn on 'places' functionality | 
  | DoSpecifyPop | Specify number of hosts (given by PopSize variable) | 
  | DoAirports | turn on 'airports' functionality | 
  | DoInitUpdateProbs |  | 
  | DoDigitalContactTracing | use digital contact tracing | 
  | ClusterDigitalContactUsers | cluster users of digital contact tracing by household | 
  | DoSI |  | 
  | DoMassVacc | implement mass vaccination programme prioritised by age | 
  | EnhancedSocDistClusterByHousehold |  | 
  | LocalBeta | | 
  | OutputBitmap | generate bitmap file of results | 
  |  |  | 
  
	
  
</details>

<details>
  <summary>Variables related to Places</summary>
  
  | Variable | Description |
  | --------------- | --------------- |
  | NPlace | number of places | 
  | PlaceTypeNum | number of place types | 
  | HotelPlaceType | index of element in places array corresponding to 'hotel' | 
  | PlaceTypeGroupSizeParam1[] | group size in given place type? | 
  | PlaceTypeTrans[] | transmission rate in given place type | 
  | SymptPlaceTypeContactRate[] | rate of those in given place type displaying symptoms who will use contact tracing? | 
  | SymptPlaceTypeWithdrawalProp[] |  | 
  | PlaceTypePropBetweenGroupLinks[] |  | 
</details>

<details>
  <summary>Variables related to administrative units</summary>
  
  | Variable | Description |
  | --------------- | --------------- |
  | AdunitLevel1Mask |  | 
  | AdunitLevel1Divisor |  | 
  | AdunitLevel1Lookup |  | 
  |  |  | 
</details>



<details>
  <summary>Random seeds</summary>
  
  | Variable | Description |
  | --------------- | --------------- |
  | setupSeed1 |  | 
  | setupSeed2 |  | 
  | nextSetupSeed1 |  | 
  | setupSeed2 |  | 
  
</details>

## Quantitative setup information

<details>
  <summary>Arrays of values for a given age group</summary>
  
  | Variable | Description |
  | --------------- | --------------- |
  | AgeInfectiousness[] | infectiousness values for given age group | 
  | ProportionSymptomatic[] | proportion of those that display symptoms for given age group | 
  | RelativeSpatialContact[] |  | 
  | infectiousness[] |  | 
  | EnhancedSocDistProportionCompliant[] | proportion of given age group likely to comply with social distancing | 
  | ProportionSmartphoneUsersByAge[] | proportion of population of that age that use smart phones - used to determine likelihood of a person of a given age using a smart phone and in turn digital contact tracing  |
</details>


<details>
  <summary>Constant values</summary>
  
  | Variable | Description |
  | --------------- | --------------- |
  | VaccPriorityGroupAge[] | array containing two elements: the min and max age for vaccination priority group | 
  | VaccProp | proportion of population to be vaccinated | 
  | SymptInfectiousnes | infectiousness rate when symptomatic | 
  | InfectiousPeriod | amount of time an infected host is infectious for | 
  | infectious_icdf |  | 
  | SymptSpatialContactRate |  | 
  | LatentToSymptDelay |  |
  | PropPopUsingDigitalContactTracing | proportion of population using digital contact tracing |
  | PopSize | number of hosts |
  |  |  |
  |  |  |

</details>

## Data produced about setup process

<details>
  <summary>Running totals</summary>
  
  | Variable | Description |
  | --------------- | --------------- |
  | NDigitalContactUsers | Count of number of digitabl contact tracing app users | 
  | NDigitalHouseholdUsers | Count of househods using digital contact tracing app | 
  |  |  | 
  |  |  | 
  |  |  | 
</details>

## Variables used to load / save information

<details>
  <summary>Bitmap loading / saving</summary>
  
  | Variable | Description |
  | --------------- | --------------- |
  | SpatialBoundingBox[] | array of four elements - starts large then adjusted to immediately surround x,y | 
  | width | width of calculated spatial bounding box - SpatialBoundingBox[2] - SpatialBoundingBox[0]. for non-heterodensity, width = squre root of population size? | 
  | height | height of calculated spatial bounding box - SpatialBoundingBox[3] - SpatialBoundingBox[1]. for non-heterodensity, width = squre root of population size? | 
  | CountryDivisor |  | 
  | LongitudeCutLine |  | 
  | BinFileLen | number of lines in input file - used for file reading |
  | BinFileBuf | contents of input file - used for file reading |
  | BitmapScale | x-scale |
  | BitmapAspectScale | aspect ratio of bitmap |
  | scalex |  |
  | scaley |  |
  | bwidth | bitmap width |
  | bheight | bitmap height |
  | bheight2 | bheight + 20 to shift above legend |
  |  |  |
  |  |  |
  |  |  |
</details>

<details>
  <summary>Memory / processor management</summary>
  
  | Variable | Description |
  | --------------- | --------------- |
  | KernelShape |  | 
  | KernelScale |  | 
  | KernelP3 |  | 
  | KernelP4 |  | 
  | KernelType |  | 
  | MoveKernelShape |  | 
  | MoveKernelScale |  | 
  | MoveKernelP3 |  | 
  | MoveKernelP4 |  | 
  | MoveKernelType |  | 
  |  |  | 
  |  |  | 
  |  |  | 
</details>


<details>
  <summary></summary>
  
  | Variable | Description |
  | --------------- | --------------- |
  |  |  | 
</details>
