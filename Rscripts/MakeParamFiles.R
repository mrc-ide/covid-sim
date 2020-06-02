#### Author: dlaydon

## Script automates writing parameter files and pre-parameter files, and documents parameters for CovidSim. See covid-sim/docs/inputs-and-outputs.md and covid-sim/docs/intervention-description.md
## Script is a work in progress: if you see a parameter for which you can add to the documentation, please open a PR!
## Broadly, pre-parameter files contain parameters whose values are common to a series of runs. Their values reflect current knowledge of SARS-CoV-2/Covid-19. 
## Parameter files contain parameters for non-pharmaceutical interventions (NPIs), and thus will differ between a series of runs if modelling multiple scenarios.
## Both parameter and pre-parameters have the same format, which is a sequence of:

##		[Description of Parameter]
##		value

## We have endeavored to name the parameters as they appear in the Cpp code.

## == ## == ## == ## == ## == ## == ## == ## == ## == ## == ## == ## == 
## Issues / still to do: 
#		i) 		more comments for parameters; 
#		ii) 	hard-coded parameters in MakePreParamList (any that need changing between countries/populations please add to arguments accordingly)
# 		iii)	more example parameter files 	



NUM_AGE_GROUPS 	= 17
CDF_RES			= 20
options("scipen" = 13) #### set high penalty for scientific display. (so e.g. 10000 is not outputted as 1e+05) 

OutputDir  = "." ## change as appropriate

MakePreParamList = function(DoPlaceCloseOnceOnly = 0, DoSocDistOnceOnly = 0, NumRealisations = 1, NumNonExtinctRealisations = NumRealisations, 
		
		Mean_MildToRecovery 				= 7       ,
		Mean_ILIToRecovery 					= 7       ,
		Mean_ILIToSARI 						= 5       ,
		Mean_ILIToDeath 					= 7       ,
		Mean_SARIToRecovery 				= 1       ,
		Mean_SARIToDeath 					= 1       ,
		Mean_SARIToCritical 				= 1       ,
		Mean_CriticalToCritRecov 			= 1       ,
		Mean_CriticalToDeath 				= 1       ,
		Mean_CritRecovToRecov 				= 1       ,
		MeanTimeToTest 						= 4       ,
		MeanTimeToTestOffset 				= 1       ,
		MeanTimeToTestCriticalOffset 		= 3.3     ,
		MeanTimeToTestCritRecovOffset 		= 9.32    ,
		
		MildToRecovery_icdf 		= c(0	, 0.341579599	, 0.436192391	, 0.509774887	, 0.574196702	, 0.633830053	, 0.690927761	, 0.74691114	, 0.802830695	, 0.859578883	, 0.918015187	, 0.97906363	, 1.043815683	, 1.113669859	, 1.190557274	, 1.277356871	, 1.378761429	, 1.50338422	, 1.670195767	, 1.938414132	, 2.511279379	) , 
		ILIToRecovery_icdf  		= c(0	, 0.341579599	, 0.436192391	, 0.509774887	, 0.574196702	, 0.633830053	, 0.690927761	, 0.74691114	, 0.802830695	, 0.859578883	, 0.918015187	, 0.97906363	, 1.043815683	, 1.113669859	, 1.190557274	, 1.277356871	, 1.378761429	, 1.50338422	, 1.670195767	, 1.938414132	, 2.511279379	) , 
		ILIToSARI_icdf 				= c(0	, 0.341579599	, 0.436192391	, 0.509774887	, 0.574196702	, 0.633830053	, 0.690927761	, 0.74691114	, 0.802830695	, 0.859578883	, 0.918015187	, 0.97906363	, 1.043815683	, 1.113669859	, 1.190557274	, 1.277356871	, 1.378761429	, 1.50338422	, 1.670195767	, 1.938414132	, 2.511279379	) , 
		ILIToDeath_icdf 			= c(0	, 2.257735908	, 3.171065856	, 3.924183798	, 4.608738224	, 5.260437017	, 5.898728066	, 6.53669783	, 7.184755068	, 7.852438367	, 8.549591424	, 9.287408763	, 10.07967529	, 10.94457146	, 11.90769274	, 13.00769447	, 14.3081531	, 15.92655201	, 18.12320384	, 21.71626849	, 29.58154704	) ,
		SARIToRecovery_icdf 		= c(0	, 0.634736097	, 1.217461548	, 1.805695261	, 2.41206761	, 3.044551205	, 3.71010552	, 4.415905623	, 5.170067405	, 5.982314035	, 6.864787504	, 7.833196704	, 8.908589322	, 10.12027655	, 11.51100029	, 13.14682956	, 15.13821107	, 17.69183155	, 21.27093904	, 27.35083955	, 41.35442157	) , 
		SARIToDeath_icdf 			= c(0	, 1.703470233	, 2.39742257	, 2.970367222	, 3.491567676	, 3.988046604	, 4.474541783	, 4.960985883	, 5.455292802	, 5.964726999	, 6.496796075	, 7.06004732	, 7.665014091	, 8.325595834	, 9.061367792	, 9.901900127	, 10.8958347	, 12.133068		, 13.81280888	, 16.56124574	, 22.5803431	) , 
		SARIToCritical_icdf  		= c(0	, 0.108407687	, 0.220267228	, 0.337653773	, 0.46159365	, 0.593106462	, 0.733343356	, 0.88367093	, 1.045760001	, 1.221701998	, 1.414175806	, 1.62669998	, 1.864032461	, 2.132837436	, 2.442868902	, 2.809242289	, 3.257272257	, 3.834402667	, 4.647120033	, 6.035113821	, 9.253953212	) , 
		CriticalToCritRecov_icdf 	= c(0	, 1.308310071	, 1.87022015	, 2.338694632	, 2.76749788	, 3.177830401	, 3.581381361	, 3.986127838	, 4.398512135	, 4.824525291	, 5.270427517	, 5.743406075	, 6.252370864	, 6.809125902	, 7.430338867	, 8.141231404	, 8.983341913	, 10.03350866	, 11.46214198	, 13.80540164	, 18.95469153	) , 
		CriticalToDeath_icdf 		= c(0	, 1.60649128	, 2.291051747	, 2.860938008	, 3.382077741	, 3.880425012	, 4.37026577	, 4.861330415	, 5.361460943	, 5.877935626	, 6.4183471		, 6.991401405	, 7.607881726	, 8.282065409	, 9.034104744	, 9.894486491	, 10.91341144	, 12.18372915	, 13.9113346	, 16.74394356	, 22.96541429	) , 
		CritRecovToRecov_icdf 		= c(0	, 0.133993315	, 0.265922775	, 0.402188416	, 0.544657341	, 0.694774487	, 0.853984373	, 1.023901078	, 1.206436504	, 1.403942719	, 1.619402771	, 1.856711876	, 2.121118605	, 2.419957988	, 2.763950408	, 3.169692564	, 3.664959893	, 4.301777536	, 5.196849239	, 6.7222126		, 10.24997697	) , 
		
		Prop_Mild_ByAge 			= c(0.666244874	,	0.666307235	,	0.666002907	,	0.665309462	,	0.663636419	,	0.660834577	,	0.657465236	,	0.65343285	,	0.650261465	,	0.64478501	,	0.633943755	,	0.625619329	,	0.609080537	,	0.600364976	,	0.5838608	,	0.566553872	,	0.564646465  ) , 
		Prop_ILI_ByAge 				= c(0.333122437	,	0.333153617	,	0.333001453	,	0.332654731	,	0.33181821	,	0.330417289	,	0.328732618	,	0.326716425	,	0.325130732	,	0.322392505	,	0.316971878	,	0.312809664	,	0.304540269	,	0.300182488	,	0.2919304	,	0.283276936	,	0.282323232  ) , 
		Prop_SARI_ByAge 			= c(0.000557744	,	0.000475283	,	0.000877703	,	0.001794658	,	0.004006955	,	0.007711884	,	0.012167229	,	0.017359248	,	0.021140307	,	0.027047193	,	0.03708932	,	0.039871236	,	0.040788928	,	0.027444452	,	0.101605674	,	0.142001415	,	0.150469697  ) , 
		Prop_Critical_ByAge 		= c(7.49444E-05	,	6.38641E-05	,	0.000117937	,	0.000241149	,	0.000538417	,	0.00103625	,	0.001634918	,	0.002491477	,	0.003467496	,	0.005775292	,	0.011995047	,	0.021699771	,	0.045590266	,	0.072008084	,	0.022603126	,	0.008167778	,	0.002560606  ) , 
		CFR_ILI_ByAge 				= rep(0, NUM_AGE_GROUPS),
		CFR_SARI_ByAge 				= c(0.125893251	,	0.12261338	,	0.135672867	,	0.152667869	,	0.174303077	,	0.194187895	,	0.209361731	,	0.224432564	,	0.237013516	,	0.257828065	,	0.290874602	,	0.320763971	,	0.362563751	,	0.390965457	,	0.421151485	,	0.447545892	,	0.482        ) , 
		CFR_Critical_ByAge 			= c(0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896	,	0.5234896    )  

)
{
	### === ### Returns list of parameters for pre-parameter file. 
	### === Parameters are hard-coded except those that may change between countries or model runs. These are given as arguments.
	### === For changes to hard-coded parameters, please add them to function arguments.
	
	PreParamList = list()
	
	PreParamList[["Output every realisation"]] = 1
	PreParamList[["Output bitmap"]] = 0
	PreParamList[["OutputAge"]] = 1
	PreParamList[["OutputSeverityAdminUnit"]] = 1
	PreParamList[["OutputR0"]] = 0
	PreParamList[["OutputControls"]] = 0
	PreParamList[["OutputCountry"]] = 0
	PreParamList[["OutputAdUnitVar"]] = 0
	PreParamList[["OutputHousehold"]] = 0
	PreParamList[["OutputInfType"]] = 0
	PreParamList[["OutputNonSeverity"]] = 0
	PreParamList[["OutputNonSummaryResults"]] = 0
	PreParamList[["Output incidence by administrative unit"]] = 0
	PreParamList[["Output infection tree"]] = 0
	PreParamList[["Only output non-extinct realisations"]] = 0
	PreParamList[["Include administrative units within countries"]] = 1
	PreParamList[["Update timestep"]] = 0.25
	PreParamList[["Equilibriation time"]] = 0
	PreParamList[["Sampling timestep"]] = 1
	PreParamList[["Sampling time"]] = 720
	PreParamList[["Grid size"]] = 0.075
	PreParamList[["Spatial domain for simulation"]] = matrix(c(73, 6.3, 136, 54), nrow = 2, byrow = TRUE)
	PreParamList[["Number of micro-cells per spatial cell width"]] = 9
	PreParamList[["Initial immunity profile by age"]] = rep(0, NUM_AGE_GROUPS)
	PreParamList[["Initial immunity applied to all household members"]] = 1
	PreParamList[["Relative spatial contact rates by age"]] = c(0.6, 0.7, 0.75, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.75, 0.5) ### (POLYMOD, averaging 20-70)
	PreParamList[["Household attack rate"]] = 0.1 ## (Adjusted to be the same as Cauchemez 2004 for R0=1.3.)
	PreParamList[["Household transmission denominator power"]] = 0.8 ## (Cauchemez 2004)
	PreParamList[["Relative transmission rates for place types"]] = c(0.14, 0.14, 0.1, 0.07) ## (School=2 x workplace. This gives Longini AJE 1988 age-specific infection attack rates for R0=1.3. Also comparable with 1957 pandemic attack rates from Chin.)
	PreParamList[["Proportion of between group place links"]] = c(0.25, 0.25, 0.25, 0.25) ## (25% of within-group contacts)
	PreParamList[["Include symptoms"]] = 1
	PreParamList[["Delay from end of latent period to start of symptoms"]] = 0.5 ## (assume average time to symptom onset is half a day)
	PreParamList[["Proportion symptomatic by age group"]] = rep(0.66, NUM_AGE_GROUPS)
	PreParamList[["Symptomatic infectiousness relative to asymptomatic"]] = 1.5
	PreParamList[["Relative rate of random contacts if symptomatic"]] = 0.5
	PreParamList[["Relative level of place attendance if symptomatic"]] = c(0.25, 0.25, 0.5, 0.5)
	PreParamList[["Model symptomatic withdrawal to home as true absenteeism"]] = 1
	PreParamList[["Maximum age of child at home for whom one adult also stays at home"]] = 16
	PreParamList[["Proportion of children at home for whom one adult also stays at home"]] = 1
	PreParamList[["Duration of place absenteeism for cases who withdraw"]] = 7
	PreParamList[["Place close round household"]] = 1
	PreParamList[["Absenteeism place closure"]] = 0
	PreParamList[["Initial number of infecteds"]] = 1000
	PreParamList[["Time when infection rate changes"]] = 30
	PreParamList[["Initial rate of importation of infections"]] = 0
	PreParamList[["Changed rate of importation of infections"]] = 0
	PreParamList[["Length of importation time profile provided"]] = 0
	PreParamList[["Daily importation time profile"]] = c(
			0.00039637	, 0.000422404	, 0.000603232	, 0.000806818	, 0.000977838	, 0.001215654	, 0.001475435	, 0.001759331	, 
			0.002095855	, 0.002476933	, 0.002909161	, 0.003410198	, 0.003982451	, 0.004638805	, 0.005395127	, 0.00626298	,	
			0.007261369	, 0.008410463	, 0.009732184	, 0.011253768	, 0.01300553	, 0.015022328	, 0.017344927	, 0.020019826	, 
			0.023100555	, 0.026649391	, 0.030737285	, 0.035446471	, 0.04087166	, 0.04712199	, 0.054323081	, 0.062619805	,
			0.072178958	, 0.083192875	, 0.095883242	, 0.110505088	, 0.127353055	, 0.146765647	, 0.169134051	, 0.194907414	, 
			0.224605405	, 0.258825037	, 0.298255169	, 0.343688393	, 0.396040694	, 0.456363956	, 0.525872339	, 0.605965353	, 
			0.698253328	, 0.804594637	, 0.927128729	, 1.068321589	, 1.231009953	, 1.418469229	, 1.634476217	, 1.883366252	, 
			2.170149943	, 2.500595049	, 2.881346836	, 3.320054954	, 3.825541988	, 4.407960394	, 5.07901906	, 5.852193458	, 
			6.742999237	, 7.769281222	, 8.951656384	, 10.31379847	, 11.88295306	, 13.690484		, 15.7725322	, 18.17054397	, 
			20.93225706	, 24.11247717	, 27.77414928	, 31.98955641	, 36.84153527	, 42.42484495	, 48.84810841	, 56.23529325	, 
			64.72777172	, 74.48650501	, 85.69401624	, 98.55689323	, 113.307934	, 130.2079818	, 149.5482524	, 171.650606	, 
			196.8670786	, 225.5805289	, 258.1996665	, 295.1479906	, 336.8602628	, 383.7567184	, 436.226154	, 494.5852818	, 
			559.0401758	, 629.6271715	, 706.1466306	, 788.0877032	, 874.5385994	, 964.1277368	, 1054.942274	, 1144.497997	, 
			1229.819016	, 1307.478109	, 1373.849859	, 1425.375658	, 1458.909049	, 1472.112189	, 1463.74025	, 1433.873478	, 
			1383.937547	, 1316.559006	, 1235.24025	, 1143.978092	, 1046.827852	, 947.5742692	, 849.4524822	, 755.0099762	, 
			666.1067926	, 583.9342875	, 509.1260534	, 441.8703774	, 382.0284161	, 329.2333209	, 282.9819783	, 242.6943429	, 
			207.7671547	, 177.6025491	, 151.6341438	, 129.3333359	, 110.223406	, 93.87499488	, 79.90856345	, 67.99027638	, 
			57.82900897	, 49.17215405	, 41.80131863	, 35.5286855	, 30.19272742	, 25.65482116	, 21.79680064	, 18.51741927	, 
			15.73036742	, 13.3620726	, 11.34984238	, 9.640307581	, 8.188006869	, 6.954337378	, 5.906435165	, 5.016354116	, 
			4.260368274	, 3.61826508	, 3.072912654	, 2.609742089	, 2.216371142	, 1.882286966	, 1.598554546	, 1.357591285	, 
			1.152944119	, 0.979149081	, 0.831549648	, 0.70620165	, 0.59974816	, 0.509342763	, 0.43256591	, 0.367362791	, 
			0.311990061	, 0.264964778	, 0.22502858	, 0.19111364	, 0.162311178	, 0.137852214	, 0.117080312	, 0.099440002	, 
			0.084459683	, 0.07173828	, 0.06093513	, 0.051761119	)
	PreParamList[["Reproduction number"]] = 2
	PreParamList[["Power of scaling of spatial R0 with density"]] = 0
	PreParamList[["Include latent period"]] = 1
	PreParamList[["Latent period"]] = 4.59 # - minus half a day to account for infectiousness pre symptom onset
	PreParamList[["Latent period inverse CDF"]] = c(
			0			, 0.098616903	, 0.171170649	, 0.239705594	, 0.307516598,	0.376194441,	0.446827262	, 
			0.520343677	, 0.597665592	, 0.679808341	, 0.767974922	, 0.863671993,	0.968878064,	1.086313899	,	
			1.219915022	, 1.37573215	, 1.563841395	, 1.803041398	, 2.135346254,	2.694118208,	3.964172493	) 
	PreParamList[["Model time varying infectiousness"]] = 1
	PreParamList[["Infectiousness profile"]] = c(
			0.487464241	, 1				, 1.229764827	, 1.312453175	, 1.307955665	, 1.251658756	, 1.166040358	, 1.065716869	, 
			0.960199498	, 0.855580145	, 0.755628835	, 0.662534099	, 0.577412896	, 0.500665739	, 0.432225141	, 0.371729322	, 
			0.318643018 , 0.272340645	, 0.232162632	, 0.19745264	, 0.167581252	, 0.141960133	, 0.120049578	, 0.101361532	,
			0.085459603	, 0.071957123	, 0.060514046	, 0.050833195	, 0.04265624	, 0.035759641	, 0.029950735	, 0.025064045	,
			0.02095788	, 0.017511251	, 0.014621091	, 0.012199802	, 0.010173075	, 0.008477992	, 0.007061366	, 0.005878301	,
			0.00489096	, 0.004067488	, 0.003381102	, 0.00280931	, 0.002333237	, 0.001937064	, 0.001607543	, 0.001333589	,
			0.001105933	, 0.00091683	, 0.000759816	, 0.000629496	, 0.000521372	, 0.000431695	, 0.000357344	, 0.000295719	, 
			0.000244659	)
	PreParamList[["Infectious period"]] = 14
	PreParamList[["Infectious period inverse CDF"]] = c(
			0			, 0.171566836	, 0.424943468, 0.464725594	, 0.50866631	, 0.55773764	, 0.613298069, 
			0.67732916	, 0.752886568	, 0.843151261, 0.895791527	, 0.955973422	, 1.026225109	, 1.110607115, 
			1.216272375	, 1.336349102	, 1.487791911, 1.701882384	, 1.865779085	, 2.126940581	, 2.524164972)
	PreParamList[["k of individual variation in infectiousness"]] = 1
	PreParamList[["Use global triggers for interventions"]] = 1
	PreParamList[["Use admin unit triggers for interventions"]] = 0
	PreParamList[["Number of sampling intervals over which cumulative incidence measured for global trigger"]] = 1000
	PreParamList[["Use cases per thousand threshold for area controls"]] = 0
	PreParamList[["Divisor for per-capita global threshold (default 1000)"]] = 100000
	PreParamList[["Divisor for per-capita area threshold (default 1000)"]] = 1000
	PreParamList[["Trigger alert on deaths"]] = 1
	PreParamList[["Number of deaths accummulated before alert"]] = 10000
	PreParamList[["Number of days to accummulate cases/deaths before alert"]] = 100
	PreParamList[["Day of year trigger is reached"]] = 100
	PreParamList[["Alert trigger starts after interventions"]] = 1
	PreParamList[["Day of year interventions start"]] = 76
	PreParamList[["Proportion of cases detected for treatment"]] = 1
#	PreParamList[["Treatment trigger incidence per cell"]] = 100000
	PreParamList[["Places close only once"]] = DoPlaceCloseOnceOnly
	PreParamList[["Social distancing only once"]] = DoSocDistOnceOnly
	PreParamList[["Use ICU case triggers for interventions"]] = 1
	PreParamList[["Proportion of cases detected before outbreak alert"]] = 1
	PreParamList[["Number of detected cases needed before outbreak alert triggered"]] = 0
	PreParamList[["Number of realisations"]] = NumRealisations
	PreParamList[["Number of non-extinct realisations"]] = NumNonExtinctRealisations
	PreParamList[["Maximum number of cases defining small outbreak"]] = 10000
	PreParamList[["Do one generation"]] = 0
	PreParamList[["Bitmap scale"]] = 60
	PreParamList[["Bitmap y:x aspect scaling"]] = 1.5
	PreParamList[["Calculate spatial correlations"]] = 0
	PreParamList[["Bitmap movie frame interval"]] = 300
	PreParamList[["Record infection events"]] = 0
	PreParamList[["Record infection events per run"]] = 0
	PreParamList[["Max number of infection events to record"]] = 10000000000
	PreParamList[["Limit number of infections"]] = 0
	PreParamList[["Max number of infections"]] = 10000000000
	
	PreParamList[["Do Severity Analysis"]] = 1
	
	PreParamList[["Mean_MildToRecovery"						]] = Mean_MildToRecovery 				
	PreParamList[["Mean_ILIToRecovery"						]] = Mean_ILIToRecovery 					
	PreParamList[["Mean_ILIToSARI"							]] = Mean_ILIToSARI 		
	PreParamList[["Mean_ILIToDeath"							]] = Mean_ILIToDeath 		
	PreParamList[["Mean_SARIToRecovery"						]] = Mean_SARIToRecovery 				
	PreParamList[["Mean_SARIToDeath"						]] = Mean_SARIToDeath 					
	PreParamList[["Mean_SARIToCritical"						]] = Mean_SARIToCritical 				
	PreParamList[["Mean_CriticalToCritRecov"				]] = Mean_CriticalToCritRecov 			
	PreParamList[["Mean_CriticalToDeath"					]] = Mean_CriticalToDeath 				
	PreParamList[["Mean_CritRecovToRecov"					]] = Mean_CritRecovToRecov 				
	
	PreParamList[["MeanTimeToTest"							]] = MeanTimeToTest 						
	PreParamList[["MeanTimeToTestOffset"					]] = MeanTimeToTestOffset 				
	PreParamList[["MeanTimeToTestCriticalOffset"			]] = MeanTimeToTestCriticalOffset 		
	PreParamList[["MeanTimeToTestCritRecovOffset"			]] = MeanTimeToTestCritRecovOffset
	
	PreParamList[["MildToRecovery_icdf"						]] = MildToRecovery_icdf 			
	PreParamList[["ILIToRecovery_icdf"						]] = ILIToRecovery_icdf  		
	PreParamList[["ILIToSARI_icdf"							]] = ILIToSARI_icdf 				
	PreParamList[["ILIToDeath_icdf"							]] = ILIToDeath_icdf 				
	PreParamList[["SARIToRecovery_icdf"						]] = SARIToRecovery_icdf 		
	PreParamList[["SARIToDeath_icdf"						]] = SARIToDeath_icdf 			
	PreParamList[["SARIToCritical_icdf"						]] = SARIToCritical_icdf  		
	PreParamList[["CriticalToCritRecov_icdf"				]] = CriticalToCritRecov_icdf 	
	PreParamList[["CriticalToDeath_icdf"					]] = CriticalToDeath_icdf 		
	PreParamList[["CritRecovToRecov_icdf"					]] = CritRecovToRecov_icdf 		
	
	PreParamList[["Prop_Mild_ByAge"              			]] = Prop_Mild_ByAge 		
	PreParamList[["Prop_ILI_ByAge"               			]] = Prop_ILI_ByAge 		
	PreParamList[["Prop_SARI_ByAge"							]] = Prop_SARI_ByAge 		
	PreParamList[["Prop_Critical_ByAge"						]] = Prop_Critical_ByAge 	
	
	PreParamList[["CFR_ILI_ByAge"							]] = CFR_ILI_ByAge 		
	PreParamList[["CFR_SARI_ByAge"							]] = CFR_SARI_ByAge 		
	PreParamList[["CFR_Critical_ByAge"						]] = CFR_Critical_ByAge 	
	
	PreParamList[["Mean child age gap"]] = 2 
	PreParamList[["Min adult age"]] = 19
	PreParamList[["Max MF partner age gap"]] = 5 
	PreParamList[["Max FM partner age gap"]] = 5 
	PreParamList[["Min parent age gap"]] = 19 
	PreParamList[["Max parent age gap"]] = 44 
	PreParamList[["Max child age"]] = 20 
	PreParamList[["One Child Two Pers Prob"]] = 0.08 
	PreParamList[["Two Child Three Pers Prob"]] = 0.11 
	PreParamList[["One Pers House Prob Old"]] = 0.5 
	PreParamList[["Two Pers House Prob Old"]] = 0.5 
	PreParamList[["One Pers House Prob Young"]] = 0.23 
	PreParamList[["Two Pers House Prob Young"]] = 0.23 
	PreParamList[["One Child Prob Youngest Child Under Five"]] = 0.5 
	PreParamList[["Two Children Prob Youngest Under Five"]] = 	0.0
	PreParamList[["Prob Youngest Child Under Five"]] = 	0
	PreParamList[["Zero Child Three Pers Prob"]] = 	0.25
	PreParamList[["One Child Four Pers Prob"]] = 0.2
	PreParamList[["Young And Single Slope"]] = 0.7
	PreParamList[["Young And Single"]] = 36
	PreParamList[["No Child Pers Age"]] = 44
	PreParamList[["Old Pers Age"]] = 60
	PreParamList[["Three Child Five Pers Prob"]] = 0.5
	PreParamList[["Older Gen Gap"]] = 19
	
	return(PreParamList)
}

PreParamList = MakePreParamList()

MakeParamList = function(
		
		#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
		#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
		#### ==== 			NOTE: Defaults have all interventions enabled. 
		
		#### Interventions are turned on when BOTH start time has passed (and start time + duration has not), AND when trigger threshold exceeded. See https://github.com/mrc-ide/covid-sim/blob/master/docs/intervention-description.md for further discussion.
		#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
		
		NumAdUnits = 3, DoInterventionDelaysByAdUnit = 0, VaryEfficaciesOverTime = 0, 
		
		# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
		# = # = # = # = # = # = # = # = # = # = 	CASE ISOLATION PARAMETERS
		
		##### NOTE: Default is to have case isolation ENABLED
		
		CaseIsolationTimeStartBase = 6,		# a value larger than PreParamList[["Sampling time"]] indicates that case isolation will not start, regardless of trigger threshold(s) CaseIsolation_CellIncThresh / CaseIsolation_CellIncThresh. Otherwise intervention will start once threshold exceeded. 
		CaseIsolationDelay = 1, 			# delay between case detection and case isolation
		CaseIsolationDuration = 7, 			# CaseIsolationDuration = number of days a case will self-isolate. Different from CaseIsolationPolicyDuration = number of days policy is in effect.
		CaseIsolationPolicyDuration = 364, 	# CaseIsolationDuration = number of days a case will self-isolate. Different from CaseIsolationPolicyDuration = number of days policy is in effect.
		CI_Delay_byAdUnit 			= rep(CaseIsolationDelay			, NumAdUnits),  # delay between case detection and case isolation by admin unit
		CI_PolicyDuration_byAdUnit 	= rep(CaseIsolationPolicyDuration	, NumAdUnits), 	# duration of case isolation policy by admin unit.
		
		CaseIsolation_CellIncThresh = 0, CaseIsolationProp = 0.9, CaseIsolationEffectiveness = 0.25, CaseIsolationHouseEffectiveness = 0.5, 
		
		#### Note: numbers here must match Num_CI_ChangeTimes; any times listed in CI_ChangeTimes that are before CaseIsolationTimeStartBase and after CaseIsolationPolicyDuration / CI_PolicyDuration_byAdUnit are irrelevant.
		Num_CI_ChangeTimes = 1, CI_ChangeTimes = c(CaseIsolationTimeStartBase, rep(1000000, Num_CI_ChangeTimes - 1)), ### by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		CI_SpatialAndPlaceEffects_OverTime 		= rep(CaseIsolationEffectiveness			, Num_CI_ChangeTimes), 
		CI_HouseholdEffects_OverTime 			= rep(CaseIsolationHouseEffectiveness		, Num_CI_ChangeTimes),
		CI_Prop_OverTime						= rep(CaseIsolationProp						, Num_CI_ChangeTimes),
		CI_CellIncThresh_OverTime				= rep(CaseIsolation_CellIncThresh			, Num_CI_ChangeTimes),
		
		
		# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
		# = # = # = # = # = # = # = # = # = # = 	HOUSEHOLD QUARANTINE PARAMETERS
		
		##### NOTE: Default is to have household quarantine ENABLED
		
		HQuarantineTimeStartBase = 6, ### a value larger than PreParamList[["Sampling time"]] indicates that case isolation will not start, regardless of trigger threshold(s) CaseIsolation_CellIncThresh / CaseIsolation_CellIncThresh. Otherwise intervention will start once threshold exceeded.  
		HQuarantineDelay = 0, HQuarantinePolicyDuration = 364,
		HQuarantineDelay_byAdUnit 			= rep(HQuarantineDelay			, NumAdUnits), 
		HQuarantinePolicyDuration_byAdUnit 	= rep(HQuarantinePolicyDuration	, NumAdUnits),
		
		HHQuar_CellIncThresh = 0, HQuarantineHouseDuration = 14, 
		HQuarantineHouseEffect = 1.5, HQuarantinePlaceEffect = c(0.25, 0.25, 0.25, 0.25), HQuarantineSpatialEffect = 0.25,
		HQuarantinePropHouseCompliant = 0.75, HQuarantinePropIndivCompliant = 1.0,
		
		#### Note: numbers here must match Num_HQ_ChangeTimes; any times listed in HQ_ChangeTimes that are before HQuarantineTimeStartBase and after HQuarantinePolicyDuration / HQuarantinePolicyDuration_byAdUnit are irrelevant.
		Num_HQ_ChangeTimes = 1, HQ_ChangeTimes = c(HQuarantineTimeStartBase, rep(1000000, Num_HQ_ChangeTimes - 1)) , # by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		HQ_HouseholdEffects_OverTime 	= rep(HQuarantineHouseEffect		, Num_HQ_ChangeTimes),
		HQ_SpatialEffects_OverTime 		= rep(HQuarantineSpatialEffect		, Num_HQ_ChangeTimes),
		HQ_PlaceEffects_OverTime		= matrix(rep(HQuarantinePlaceEffect	, Num_HQ_ChangeTimes), nrow = Num_HQ_ChangeTimes, byrow = TRUE),
		
		HQ_Household_PropComply_OverTime 	= rep(HQuarantinePropHouseCompliant	, Num_HQ_ChangeTimes), # household compliance
		HQ_Individual_PropComply_OverTime 	= rep(HQuarantinePropIndivCompliant	, Num_HQ_ChangeTimes), # individual compliance
		
		HQ_CellIncThresh_OverTime 			= rep(HHQuar_CellIncThresh			, Num_HQ_ChangeTimes), # thresholds
		
		
		# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
		# = # = # = # = # = # = # = # = # = # = 	PLACE CLOSURE PARAMETERS
		
		##### NOTE: Default is to have place closure ENABLED
		
		PlaceCloseTimeStartBase = 7, PlaceCloseTimeStartBase2 = 100000, PlaceCloseDelayMean = 1, PlaceCloseDurationBase = 364, PlaceCloseRadius = 1,  
		PlaceCloseByAdminUnit = 0, PlaceCloseAdminUnitDivisor = 1, PlaceCloseAdunitPlaceTypes = c(1, 1, 1, 0), PlaceCloseCasePropThresh = 1, PlaceCloseAdunitPropThresh = 1,
		PlaceCloseDelayMean_byAdUnit 	= rep(PlaceCloseDelayMean		, NumAdUnits),
		PlaceCloseDuration_byAdUnit 	= rep(PlaceCloseDurationBase	, NumAdUnits),
		PlaceCloseHouseholdRelContact 		= 1.5, 
		PlaceCloseSpatialRelContact 		= 1.25, 
		PropPlacesRemainingOpen_byPlaceType	= c(0, 0, 0.25, 1), 
		PlaceClosePropAttending = c(0.05, 0.05, 0.05, 0.05), ### Partial closure: for places that are closed, what proportion of their members do attend? Set to zero if no partiall admittance. 
		PlaceCloseIncTrig = 0, PlaceCloseFracIncTrig = 0, PlaceCloseCellIncThresh = 0,
		
		#### Note: numbers here must match Num_PC_ChangeTimes; any times listed in PC_ChangeTimes that are before PlaceCloseTimeStartBase and after PlaceCloseDurationBase / PlaceCloseDuration_byAdUnit are irrelevant.
		Num_PC_ChangeTimes = 1,	PC_ChangeTimes = c(PlaceCloseTimeStartBase, rep(1000000, Num_PC_ChangeTimes - 1))		, ### by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		
		PC_IncThresh_OverTime 			= rep(PlaceCloseIncTrig				, Num_PC_ChangeTimes)	,
		PC_FracIncThresh_OverTime 		= rep(PlaceCloseFracIncTrig			, Num_PC_ChangeTimes) 	, 
		PC_CellIncThresh_OverTime 		= rep(PlaceCloseCellIncThresh		, Num_PC_ChangeTimes) 	, 
		PC_HouseholdEffects_OverTime 	= rep(PlaceCloseHouseholdRelContact	, Num_PC_ChangeTimes)	,
		PC_SpatialEffects_OverTime 		= rep(PlaceCloseSpatialRelContact	, Num_PC_ChangeTimes)	,
		
		PC_PlaceEffects_OverTime 		= matrix(rep(PropPlacesRemainingOpen_byPlaceType, Num_PC_ChangeTimes), nrow = Num_PC_ChangeTimes, byrow = TRUE),
		PC_Durs_OverTime				= rep(PlaceCloseDurationBase		, Num_PC_ChangeTimes)	,
		PC_PropAttending_OverTime		= matrix(rep(PlaceClosePropAttending			, Num_PC_ChangeTimes), nrow = Num_PC_ChangeTimes, byrow = TRUE),
		
		# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
		# = # = # = # = # = # = # = # = # = # = 	SOCIAL DISTANCING PARAMETERS
		
		##### NOTE: Default is to have social distancing ENABLED
		
		SocDistTimeStartBase = 6, SocDistDuration = 364, 
		SocDistDelay_byAdUnit 		= rep(0, NumAdUnits),
		SocDistDuration_byAdUnit 	= rep(SocDistDuration, NumAdUnits),
		
		SocDistCellIncThresh = 0, 
		
		SocDistPlaceEffect 		= c(1, 1, 0.5, 0.5)	, EnhancedSocDistPlaceEffect 		= c(0.25, 0.25, 0.25, 0.25), 
		SocDistHouseholdEffect 	= 1.25				, EnhancedSocDistHouseholdEffect 	= 1, 
		SocDistSpatialEffect 	= 0.1				, EnhancedSocDistSpatialEffect 		= 0.25, 
		
		EnhancedSocDistProportionCompliantNonAge = 0, EnhancedSocDistProportionCompliantByAge = rep(0, NUM_AGE_GROUPS),
		
		#### Note: numbers here must match Num_SD_ChangeTimes; any times listed in SD_ChangeTimes that are before SocDistTimeStartBase and after SocDistDuration / SocDistDuration_byAdUnit are irrelevant.
		Num_SD_ChangeTimes = 1, SD_ChangeTimes = c(SocDistTimeStartBase, rep(1000000, Num_SD_ChangeTimes - 1)), ### by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		
		SD_PlaceEffects_OverTime 				= matrix(rep(SocDistPlaceEffect			, Num_SD_ChangeTimes), nrow = Num_SD_ChangeTimes, byrow = TRUE),
		Enhanced_SD_PlaceEffects_OverTime 		= matrix(rep(EnhancedSocDistPlaceEffect	, Num_SD_ChangeTimes), nrow = Num_SD_ChangeTimes, byrow = TRUE),
		SD_HouseholdEffects_OverTime 			= rep(SocDistHouseholdEffect			, Num_SD_ChangeTimes), 
		Enhanced_SD_HouseholdEffects_OverTime	= rep(EnhancedSocDistHouseholdEffect	, Num_SD_ChangeTimes),
		SD_SpatialEffects_OverTime				= rep(SocDistSpatialEffect				, Num_SD_ChangeTimes), 
		Enhanced_SD_SpatialEffects_OverTime		= rep(EnhancedSocDistSpatialEffect		, Num_SD_ChangeTimes), 
		
		SD_CellIncThresh_OverTime 				= rep(SocDistCellIncThresh				, Num_SD_ChangeTimes) 
)
{
	#### Add checks - e.g. number of change times must match number of levels. Add warnings for change times set before and after start and start + duration of policy
	
	#### if change times not provided, default is initialize time to first start time and then arbitrarily large change times (i.e. they won't change). If change times provided, default is to initialize Num_XX_ChangeTimes to equal length XX_ChangeTimes
	if (is.null(CI_ChangeTimes	)) CI_ChangeTimes 	= c(CaseIsolationTimeStartBase			, rep(1000000, Num_CI_ChangeTimes 	- 1)) else Num_CI_ChangeTimes 	= length(CI_ChangeTimes)
	if (is.null(HQ_ChangeTimes	)) HQ_ChangeTimes 	= c(HQuarantineTimeStartBase			, rep(1000000, Num_HQ_ChangeTimes 	- 1)) else Num_HQ_ChangeTimes 	= length(HQ_ChangeTimes)
	if (is.null(PC_ChangeTimes	)) PC_ChangeTimes 	= c(PlaceCloseTimeStartBase				, rep(1000000, Num_PC_ChangeTimes 	- 1)) else Num_PC_ChangeTimes 	= length(PC_ChangeTimes)
	if (is.null(SD_ChangeTimes	)) SD_ChangeTimes 	= c(SocDistTimeStartBase				, rep(1000000, Num_SD_ChangeTimes 	- 1)) else Num_SD_ChangeTimes 	= length(SD_ChangeTimes)
	
	
	ParamList = list()
	
	ParamList[["Include intervention delays by admin unit"]] = DoInterventionDelaysByAdUnit
	ParamList[["Vary efficacies over time"]] = VaryEfficaciesOverTime
	
	# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
	# = # = # = # = # = # = # = # = # = # = 	CASE ISOLATION PARAMETERS
	
	ParamList[["Case isolation start time"]] 					= CaseIsolationTimeStartBase	# Start of case isolation policy.
	ParamList[["Duration of case isolation"]] 					= CaseIsolationDuration  		# CaseIsolationDuration = number of days a case will self-isolate. Different from CaseIsolationPolicyDuration = number of days policy is in effect.
	ParamList[["Duration of case isolation policy"]] 			= CaseIsolationPolicyDuration	# CaseIsolationDuration = number of days a case will self-isolate. Different from CaseIsolationPolicyDuration = number of days policy is in effect.
	ParamList[["Delay to start case isolation"]] 				= CaseIsolationDelay 		    # delay between case detection and case isolation
	
	ParamList[["Delay to case isolation by admin unit"]] 		= CI_Delay_byAdUnit
	ParamList[["Duration of case isolation by admin unit"]] 	= CI_PolicyDuration_byAdUnit   ### delay between case detection and case isolation by admin unit
	
	ParamList[["Case isolation trigger incidence per cell"]] 	= CaseIsolation_CellIncThresh
	ParamList[["Proportion of detected cases isolated"]] 		= CaseIsolationProp
	
	ParamList[["Residual contacts after case isolation"]] 				= CaseIsolationEffectiveness
	ParamList[["Residual household contacts after case isolation"]] 	= CaseIsolationHouseEffectiveness
	ParamList[["Number of change times for levels of case isolation"]] 	= Num_CI_ChangeTimes
	
	ParamList[["Change times for levels of case isolation"]] 					= CI_ChangeTimes #### Note: numbers here must match "Number of change times for levels of case isolation"; that any times listed here that are before "Case isolation start time" and after "Duration of case isolation policy" are irrelevant.
	ParamList[["Residual contacts after case isolation over time"]] 			= CI_SpatialAndPlaceEffects_OverTime
	ParamList[["Residual household contacts after case isolation over time"]] 	= CI_HouseholdEffects_OverTime
	ParamList[["Proportion of detected cases isolated over time"]] 				= CI_Prop_OverTime
	ParamList[["Case isolation trigger incidence per cell over time"]] 			= CI_CellIncThresh_OverTime
	
	
	# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
	# = # = # = # = # = # = # = # = # = # = 	HOUSEHOLD QUARANTINE PARAMETERS
	
	ParamList[["Household quarantine start time"]] 			= HQuarantineTimeStartBase
	ParamList[["Delay to start household quarantine"]] 		= HQuarantineDelay
	ParamList[["Duration of household quarantine policy"]] 	= HQuarantinePolicyDuration
	
	ParamList[["Delay to household quarantine by admin unit"]] 						= HQuarantineDelay_byAdUnit
	ParamList[["Duration of household quarantine by admin unit"]] 					= HQuarantinePolicyDuration_byAdUnit
	
	ParamList[["Household quarantine trigger incidence per cell"]] 					= HHQuar_CellIncThresh
	ParamList[["Length of time households are quarantined"]] 						= HQuarantineHouseDuration
	ParamList[["Relative household contact rate after quarantine"]] 				= HQuarantineHouseEffect
	ParamList[["Residual place contacts after household quarantine by place type"]] = HQuarantinePlaceEffect
	ParamList[["Residual spatial contacts after household quarantine"]] 			= HQuarantineSpatialEffect
	ParamList[["Household level compliance with quarantine"]] 						= HQuarantinePropHouseCompliant
	ParamList[["Individual level compliance with quarantine"]] 						= HQuarantinePropIndivCompliant
	
	ParamList[["Number of change times for levels of household quarantine"]] 					= Num_HQ_ChangeTimes
	ParamList[["Change times for levels of household quarantine"]] 								= HQ_ChangeTimes
	ParamList[["Relative household contact rates over time after quarantine"]] 					= HQ_HouseholdEffects_OverTime
	ParamList[["Residual place contacts over time after household quarantine by place type"]] 	= HQ_PlaceEffects_OverTime
	ParamList[["Residual spatial contacts over time after household quarantine"]] 				= HQ_SpatialEffects_OverTime
	
	ParamList[["Household level compliance with quarantine over time"]] 		= HQ_Household_PropComply_OverTime
	ParamList[["Individual level compliance with quarantine over time"]] 		= HQ_Individual_PropComply_OverTime
	ParamList[["Household quarantine trigger incidence per cell over time"]] 	= HQ_CellIncThresh_OverTime
	
	
	
	
	# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
	# = # = # = # = # = # = # = # = # = # = 	PLACE CLOSURE PARAMETERS
	
	ParamList[["Place closure start time"]] = PlaceCloseTimeStartBase
	ParamList[["Place closure second start time"]] = PlaceCloseTimeStartBase2
	ParamList[["Delay to start place closure"]] = PlaceCloseDelayMean
	ParamList[["Duration of place closure"]] = PlaceCloseDurationBase
	
	ParamList[["Delay to place closure by admin unit"]] = PlaceCloseDelayMean_byAdUnit
	ParamList[["Duration of place closure by admin unit"]] = PlaceCloseDuration_byAdUnit
	
	ParamList[["Place closure in administrative units rather than rings"]] 	= PlaceCloseByAdminUnit
	ParamList[["Administrative unit divisor for place closure"]] 			= PlaceCloseAdminUnitDivisor
	ParamList[["Place types to close for admin unit closure (0/1 array)"]] 	= PlaceCloseAdunitPlaceTypes
	ParamList[["Minimum radius for place closure"]] 						= PlaceCloseRadius
	
	ParamList[["Cumulative proportion of place members needing to become sick for admin unit closure"]] = PlaceCloseCasePropThresh
	ParamList[["Proportion of places in admin unit needing to pass threshold for place closure"]] 		= PlaceCloseAdunitPropThresh
	ParamList[["Proportion of places remaining open after closure by place type"]] 						= PropPlacesRemainingOpen_byPlaceType
	ParamList[["Proportional attendance after closure by place type"]] 									= PlaceClosePropAttending
	
	ParamList[["Relative household contact rate after closure"]] 										= PlaceCloseHouseholdRelContact
	ParamList[["Relative spatial contact rate after closure"]] 											= PlaceCloseSpatialRelContact
	
	ParamList[["Place closure incidence threshold"]] 			= PlaceCloseIncTrig 			## needs to be 0 for global triggers
	ParamList[["Place closure fractional incidence threshold"]] = PlaceCloseFracIncTrig			## needs to be 0 for global triggers or if abs incidence threshold used
	ParamList[["Trigger incidence per cell for place closure"]] = PlaceCloseCellIncThresh 		## change this for global too ###
	
	ParamList[["Number of change times for levels of place closure"]] 							= Num_PC_ChangeTimes
	ParamList[["Change times for levels of place closure"]] 									= PC_ChangeTimes #### Note: numbers here must match "Number of change times for levels of place closure"; that any times listed here that are before "Place closure start time" and after "Duration of place closure" are irrelevant.
	ParamList[["Proportion of places remaining open after closure by place type over time"]] 	= PC_PlaceEffects_OverTime 
	ParamList[["Relative household contact rates over time after place closure"]] 				= PC_HouseholdEffects_OverTime
	ParamList[["Relative spatial contact rates over time after place closure"]] 				= PC_SpatialEffects_OverTime
	ParamList[["Place closure incidence threshold over time"]] 									= PC_IncThresh_OverTime 
	ParamList[["Place closure fractional incidence threshold over time"]] 						= PC_FracIncThresh_OverTime 		
	ParamList[["Trigger incidence per cell for place closure over time"]] 						= PC_CellIncThresh_OverTime 		
	ParamList[["Duration of place closure over time"]] 											= PC_Durs_OverTime #### Note: closure durations longer than interval between change times will be truncated
	ParamList[["Proportional attendance after closure by place type over time" ]] 				= PC_PropAttending_OverTime
	
	# = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = # = 
	# = # = # = # = # = # = # = # = # = # = 	SOCIAL DISTANCING PARAMETERS
	
	ParamList[["Social distancing start time"]] 				= SocDistTimeStartBase
	ParamList[["Duration of social distancing"]] 				= SocDistDuration
	ParamList[["Delay to social distancing by admin unit"]] 	= SocDistDelay_byAdUnit
	ParamList[["Duration of social distancing by admin unit"]] 	= SocDistDuration_byAdUnit
	
	ParamList[["Trigger incidence per cell for social distancing"]] 					= SocDistCellIncThresh
	
	ParamList[["Relative place contact rate given social distancing by place type"]] 			= SocDistPlaceEffect
	ParamList[["Relative household contact rate given social distancing"]] 						= SocDistHouseholdEffect
	ParamList[["Relative spatial contact rate given social distancing"]] 						= SocDistSpatialEffect
	ParamList[["Relative place contact rate given enhanced social distancing by place type"]] 	= EnhancedSocDistPlaceEffect
	ParamList[["Relative household contact rate given enhanced social distancing"]] 			= EnhancedSocDistHouseholdEffect
	ParamList[["Relative spatial contact rate given enhanced social distancing"]] 				= EnhancedSocDistPlaceEffect
	
	ParamList[["Minimum radius for social distancing"]] = 1
	ParamList[["Proportion compliant with enhanced social distancing"]] 				= EnhancedSocDistProportionCompliantNonAge
	ParamList[["Proportion compliant with enhanced social distancing by age group"]] 	= EnhancedSocDistProportionCompliantByAge
	
	### "after change" parameters - now superseded by "over time" parameters. 
	#ParamList[["Delay for change in effectiveness of social distancing"]] 							= 1000
	#ParamList[["Relative place contact rate given social distancing by place type after change"]] 	= c(1, 1, 0.75, 0.75)
	#ParamList[["Relative household contact rate given social distancing after change"]] 			= 1.25
	#ParamList[["Relative spatial contact rate given social distancing after change"]] 				= 0.25
	
	ParamList[["Number of change times for levels of social distancing"]] 	= Num_SD_ChangeTimes #### Must match "Change times for levels of social distancing"
	ParamList[["Change times for levels of social distancing"]] 			= SD_ChangeTimes #### Note: numbers here must match "Number of change times for levels of social distancing"; that any times listed here that are before "Social distancing start time" and after "Duration of social distancing" are irrelevant.
	
	ParamList[["Relative place contact rates over time given social distancing by place type"]] 			= SD_PlaceEffects_OverTime #### Want this to supercede "Relative place contact rate given social distancing by place type". Should be matrix of dimension "Number of change times for levels of social distancing" by Number of place types.
	ParamList[["Relative place contact rates over time given enhanced social distancing by place type"]] 	= Enhanced_SD_PlaceEffects_OverTime
	ParamList[["Relative household contact rates over time given social distancing"]] 						= SD_HouseholdEffects_OverTime ####  Ideally want this to supercede "Relative household contact rate given social distancing" but need to preserved backwards compatibility for now.
	ParamList[["Relative household contact rates over time given enhanced social distancing"]] 				= Enhanced_SD_HouseholdEffects_OverTime
	ParamList[["Relative spatial contact rates over time given social distancing"]] 						= SD_SpatialEffects_OverTime ####  Ideally want this to supercede "Relative spatial contact rate given social distancing" but need to preserved backwards compatibility for now.
	ParamList[["Relative spatial contact rates over time given enhanced social distancing"]] 				= Enhanced_SD_SpatialEffects_OverTime
	ParamList[["Trigger incidence per cell for social distancing over time"]] 								= SD_CellIncThresh_OverTime
	
	
	return(ParamList)
}



WriteParamList = function(ParamList, OutputDir, OutputFileName, PrintToConsole = FALSE)
{
	ListString = ""
	for (ParamNum in 1:length(ParamList))
	{
		ParamNameString = paste0("[", names(ParamList)[ParamNum], "]")
		if (PrintToConsole) cat(paste0(ParamNameString, "\n"))
		ListString 		= paste0(ListString, ParamNameString, "\n")
		if (class (ParamList[[ParamNum]]) == "matrix")
		{
			ParamValueString = ""
			for (row in 1:dim(ParamList[[ParamNum]])[1])
				ParamValueString = paste0(ParamValueString, paste0(ParamList[[ParamNum]][row,], collapse = " "), " \n")
			
		} else 	ParamValueString = paste0(paste0(ParamList[[ParamNum]], collapse = " "), "\n")
		if (PrintToConsole) cat(paste0(ParamValueString, "\n"))
		
		ListString = paste0(ListString, ParamValueString, "\n")
	}
	write.table(ListString, file = paste0(OutputDir, OutputFileName, ".txt"), row.names = F, col.names = F, quote = F, sep = "\t")
}


MakeAndWriteParamList = function(OutputDir, OutputFileName, PrintToConsole = FALSE, ...)
{
	ParamList = MakeParamList(...)
	WriteParamList(ParamList, OutputDir, OutputFileName, PrintToConsole)
}

#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== #### ==== 
#### ==== Examples

#### Write non-US, non-Canada pre-parameter file
PreParamList = MakePreParamList()
WriteParamList(PreParamList, OutputDir = OutputDir, OutputFileName = "pre_R0_2.0")


#### Parameter files corresponding to following set of interventions:

ArbitrarilyLargeNonStartTime = 100000 ## interventions turned off by specifying a start time greater than the simulation duration, i.e. PreParamList[["Sampling time"]]

# i) 	All interventions
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_PC_CI_HQ_SD"	)

# ii) 	No Interventions (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_NoInt", 
		CaseIsolationTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime) 

# iii) 	Case Isolation Only (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_CI", 
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime)

# iv) 	Case Isolation + Home Quarantine (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_CI_HQ", 
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime)

# v) 	Place Closure only (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_PC", 
		CaseIsolationTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		SocDistTimeStartBase 				= ArbitrarilyLargeNonStartTime) 

# vi) 	Social Distancing Only (sufficient to set start times to be larger than simulation time for interventions not implemented).
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_SD", 
		CaseIsolationTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		HQuarantineTimeStartBase 			= ArbitrarilyLargeNonStartTime, 	
		PlaceCloseTimeStartBase 			= ArbitrarilyLargeNonStartTime) 


#### Write parameter files with varying "efficacies" over time. 
MakeAndWriteParamList(OutputDir = OutputDir, OutputFileName = "p_PC_CI_HQ_SDEase_DCT", 
		
		VaryEfficaciesOverTime = 1,
		
		#### Note: numbers here must match Num_SD_ChangeTimes; any times listed in SD_ChangeTimes that are before SocDistTimeStartBase and after SocDistDuration / SocDistDuration_byAdUnit are irrelevant.
		Num_SD_ChangeTimes = 6, SD_ChangeTimes = c(0, 60, 120, 180, 240, 300), ### by default, initialize first time to non-variable time start, then all other times to be arbitrarily large. 
		
		SD_PlaceEffects_OverTime 		= matrix(c(	1, 1, 0.2, 0.2, 
						1, 1, 0.3, 0.3, 
						1, 1, 0.4, 0.4, 
						1, 1, 0.5, 0.5, 
						1, 1, 0.6, 0.6, 
						1, 1, 0.7, 0.7), nrow = 6, byrow = TRUE),
		SD_HouseholdEffects_OverTime 	= c(1.2, 1.2, 1.15, 1.15, 1.1, 1.0), 
		SD_SpatialEffects_OverTime		= c(0.1, 0.2, 0.3, 0.4, 0.5, 0.5)	)	



