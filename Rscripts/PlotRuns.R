#!/usr/bin/env Rscript
# This script creates plots to compare different scenarios.
# Command line arguments: 
#    folder: the relative folder of where to find output files.
# Will create plots in subfolders located at "folder/Plots"


#### Purpose: make plots that Neil outputs in his spreadsheets without having to use spreadsheets. ####

# Requires a set of input runs in subfolders of a base folder
# Puts out a set of plots in the Plots subfolder of the base folder
# Edit the Parameters section to customize

# code sent by d.laydon@imperial.ac.uk
# first refactoring pass wunderalbert@github.com


#### Parameters ####

# Where will the run data be located? The files are in BaseFolder/subfolder_i
# where subfolder_i is one of the possibly several elements of VariousSubFolders
DefaultBaseFolder <- "."
VariousSubFolders <- "."

# Where will the plots be outputted to?
PlotOutputDirRelativeToBaseFolder <- "Plots"

# Date range
Day_0 <- as.Date("2020-01-01")
NumDaysToPlot <- 186

# resolution of output files
PNG_res <- 300

# Which plots should be shown?
OUTPUT_AdminUnitInc		<- TRUE
OUTPUT_SEIR  			<- TRUE
OUTPUT_SevInc			<- TRUE
OUTPUT_AgeInc			<- TRUE

SeverityVariablesToPlot <- c("incSARI", "incCritical", "SARI", "Critical", "incDeath", "cumSARI", "cumCritical")

#### Set up ####
if (!require(sp)) 
{
  do_install <- readline("Need to install package `sp`. Continue (Y/n)? ")
  if (!(do_install == "" | tolower(do_install) == "y")) stop("Cannot continue script without dependencies.")
  install.packages("sp")
  library(sp)
}

# parse command line arguments
args = commandArgs(trailingOnly=TRUE)
BaseFolder 		<- if (length(args > 0)) args[1] else DefaultBaseFolder
PlotOutputDir 	<- file.path(BaseFolder, PlotOutputDirRelativeToBaseFolder)

#### Hardcoded values about the simulation API ####

# variables put out in the severity file
SeverityVariablesToPlot <- c("incSARI", "incCritical", "SARI", "Critical", "incDeath", "cumSARI", "cumCritical")
# file string to look for
Pattern 			    <- ".avNE.severity.xls"


#### Define processing functions ####

CalcPlotWindow 				<- function(EpiCurves, Threshold = 0, PlotFromTo = NULL)
{
	if (!is.null(PlotFromTo) && length(PlotFromTo) != 2) stop("PlotFromTo must be character of length 2")
	# set default
	MinTime <- 1
	MaxTime <- dim(EpiCurves)[1]
	PlotWindowIndices <- MinTime:MaxTime	
	if (is.null(PlotFromTo))
	{
		if (Threshold > 0)
		{
			MinTime <- 1
			MaxTime <- dim(EpiCurves)[1]
			for (timestep in 1:dim(EpiCurves)[1]) if (any(EpiCurves[timestep, ] > Threshold)) {	MinTime <- timestep; break } 
			for (timestep in dim(EpiCurves)[1]:1) if (any(EpiCurves[timestep, ] > Threshold)) {	MaxTime <- timestep; break } 
			if ((MaxTime - MinTime) <= 0)  ## reset in weird cases where exactly one time step is left. Breaks everything. Think of proper solution later.
			{
				MinTime <- 1
				MaxTime <- dim(EpiCurves)[1];
			}
			if (MinTime == 1 & MaxTime == dim(EpiCurves)[1]) 
				warning("CalcPlotWindow: cannot subset, try new threshold")
			PlotWindowIndices <- MinTime:MaxTime		
		}
	} else PlotWindowIndices <- as.numeric(as.Date(PlotFromTo[1]) - Day_0):as.numeric(as.Date(PlotFromTo[2]) - Day_0)
	
	return(PlotWindowIndices)
}
GetVerboseVariableString 	<- function(InfVariableString)
{
	InfVariableStringLong <- sub("inc", "", InfVariableString)
	InfVariableStringLong <- sub("cum", "", InfVariableStringLong)
	return(InfVariableStringLong)
}
GetIncPrevOrCumIncString 	<- function(InfVariableString)
{
	if (length(grep("inc", InfVariableString)) > 0) 
		IncPrevOrCumIncString <- "Incidence" 				else 
	if (length(grep("cum", InfVariableString)) > 0) 
		IncPrevOrCumIncString <- "Cumulative Incidence" 	else 
		IncPrevOrCumIncString <- "Prevalence"
	return(IncPrevOrCumIncString)
}

makeScenarioPlotDir <- function(PlotOutputSubDir, Scenario) 
{
	# make a folder to put plots in and return that folder
	ScenarioPlotDir <- file.path(PlotOutputSubDir, Scenario)
	dir.create(PlotOutputDir	, showWarnings = FALSE)
	dir.create(PlotOutputSubDir	, showWarnings = FALSE)
	dir.create(ScenarioPlotDir	, showWarnings = FALSE)
	stopifnot(dir.exists(ScenarioPlotDir))
	ScenarioPlotDir
}

#### Main processing loop ####
for (SubFolder in VariousSubFolders) 
{
	cur_path		<- file.path(BaseFolder, SubFolder)
	FilesToCheck	<- list.files(path = cur_path, pattern = Pattern)
	Scenarios     	<- sub(Pattern, "", FilesToCheck) 
	
	if (length(FilesToCheck) > 0)
	{
		PlotOutputSubDir <- file.path(PlotOutputDir, SubFolder)
		dir.create(PlotOutputSubDir, showWarnings = FALSE)
	
	}
	cat(paste0("Processing Subfolder: ", SubFolder, " (", length(Scenarios), " scenarios)"))
	for (Scenario in Scenarios) 
	{
		cat(paste0("\nScenario ", Scenario, ": "))
		
		if (OUTPUT_SEIR || OUTPUT_SevInc)
		{
			SeverityFileName 		<- paste0(cur_path, Scenario, ".avNE.severity.xls") 
			if (file.exists(SeverityFileName))
			{
				SeverityResults <- read.delim(SeverityFileName, header = TRUE)
				## clean
				for (col in dim(SeverityResults)[2]:1) if (all(is.na(SeverityResults[,col]))) SeverityResults[,col] = NULL
				Dates <- Day_0 + 0:(dim(SeverityResults)[1]-1)
				
				if (CORRECT_ZERO_SUSCEPTIBLES) 
				{
					NonCalibIndices	= min(which(SeverityResults$S != 0)):dim(SeverityResults)[1]
					SeverityResults = SeverityResults[NonCalibIndices, ] # omit all zero rows at beginning of file.
					Dates 			<- Dates[NonCalibIndices]
				}
				
				### create subfolder so more organised. Don't want empty folders so wastefully include this within if statements.
				ScenarioPlotDir = makeScenarioPlotDir(PlotOutputSubDir, Scenario)
				
				if (OUTPUT_SEIR)
				{
					cat(paste0("SEIR, "))
					##### MAKE SIR plot
					png(file = paste0(ScenarioPlotDir, "SimpleSIR.png"), res = PNG_res, units = "in", width = 7, height = 7)
					plot(Dates, SeverityResults$S, main = paste0(Scenario, "\nSIR/D"),
							col = "green", type = "l", lwd = LWD, ylab = "prevalence", xlab = "", xaxt = "n", 
							ylim = c(0, max(SeverityResults$S)))
					lines(Dates, SeverityResults$I, col = "red"		, lwd = LWD)
					lines(Dates, SeverityResults$R, col = "yellow"	, lwd = LWD)
					lines(Dates, SeverityResults$cumDeath, col = "black"	, lwd = LWD)
					axis.Date(1, at=seq(min(Dates), max(Dates), by = "months"), format="%b", las = 2)
					legend ("topright", legend = c("Susceptible", "Infectious", "Recovered", "Dead"), col = c("green", "red", "yellow", "black"), lty = 1, pch = NA, lwd = LWD, cex = 1)
					dev.off()
				}
				
				if (OUTPUT_SevInc)
				{
					cat(paste0("SeverityInc, "))
					
					### cut off low incidence weeks.
					incVariables = c("incMild", "incILI", "incSARI", "incCritical", "incDeath")
					PlotWindow			= CalcPlotWindow(SeverityResults[,incVariables], Threshold = 1)
					SeverityResults_tmp	= SeverityResults[PlotWindow, incVariables]
					Dates_tmp			= Dates[PlotWindow]
					
					##### MAKE incidence by severity plot
					png(file = paste0(ScenarioPlotDir, "SeverityIncidence.png"), res = PNG_res, units = "in", width = 7, height = 7)
					plot(Dates_tmp, SeverityResults_tmp$incMild, main = paste0(Scenario, "\nincidence by disease severity"),
							col = "pink", type = "l", lwd = LWD, ylab = "incidence", xlab = "", xaxt = "n", 
							ylim = c(0, max(SeverityResults_tmp)))
					lines(Dates_tmp, SeverityResults_tmp$incILI		, col = "palevioletred1", lwd = LWD)
					lines(Dates_tmp, SeverityResults_tmp$incSARI	, col = "orange"		, lwd = LWD)
					lines(Dates_tmp, SeverityResults_tmp$incCritical, col = "red"			, lwd = LWD)
					lines(Dates_tmp, SeverityResults_tmp$incDeath	, col = "black"			, lwd = LWD)
					#axis.Date(1, at=seq(min(Dates), max(Dates), by = "months"), format="%m-%Y", las = 2)
					axis.Date(1, at=seq(min(Dates_tmp), max(Dates_tmp), by = "months"), format="%b", las = 2)
					legend ("topright", legend = c("Mild", "ILI", "SARI", "Critical", "Dead"), col = c("pink", "palevioletred1", "orange", "red", "black"), lty = 1, pch = NA, lwd = LWD, cex = 1)
					dev.off()
					
					#### Do deaths only
					PlotWindow	 = CalcPlotWindow(SeverityResults[c("incDeath")], Threshold = 1)
					png(file = paste0(ScenarioPlotDir, "DeathIncidence.png"), res = PNG_res, units = "in", width = 7, height = 7)
					plot(Dates[PlotWindow], SeverityResults[PlotWindow, c("incDeath")], main = paste0(Scenario, "\nincidence of death"),
							col = "black", type = "l", lwd = LWD, ylab = "incidence", xlab = "", xaxt = "n")
					axis.Date(1, at=seq(min(Dates[PlotWindow]), max(Dates[PlotWindow]), by = "months"), format="%b", las = 2)
					dev.off()
				}
			}
		}
		
		if (OUTPUT_AgeInc) 
		{
			AgeFileName <- file.path(cur_path, paste0(Scenario, ".avNE.age.xls"))
			if (!file.exists(AgeFileName) ) next else
				AgeResults <- read.delim(AgeFileName, header = TRUE)
			#clean
			if (AgeResults[dim(AgeResults)[1], 1] == "dist") AgeResults = AgeResults[-dim(AgeResults)[1], ] ### remove dist calc at tail of xls. 
			
			### create subfolder so more organised. Don't want empty folders so wastefully include this within if statements.      
			ScenarioPlotDir = makeScenarioPlotDir(PlotOutputSubDir, Scenario)
			
			CasesOrDeaths <- c("C", "D")
			AgeBandChar <- function(MinAge, MaxAge, Char = ".") paste0(MinAge, Char, MaxAge)
			
			### make stacked bar chart with indicence of cases or deaths by age group. 
			for (CaseOrDeath in CasesOrDeaths) 
			{
				if (CaseOrDeath == "C") CaseOrDeath_long = "Cases" else if (CaseOrDeath == "D") CaseOrDeath_long = "Deaths"
				cat(paste0(CaseOrDeath_long, " AgeInc, "))
				
				### create weekly incidence by 10 year age group.
				### convert 5 years to 10 years
				CumCases_10yBands <- matrix(nrow = dim(AgeResults)[1], ncol = 9)
				colnames(CumCases_10yBands) <- AgeBandChar(seq(0, 80, by = 10), c(seq(10, 80, by = 10), "85"))
				CumCases_10yBands[,  "0.10"] <- AgeResults[, paste0(CaseOrDeath, "0.5")] + AgeResults[, paste0(CaseOrDeath, "5.10")] 
				CumCases_10yBands[, "10.20"] <- AgeResults[, paste0(CaseOrDeath, "10.15")] + AgeResults[, paste0(CaseOrDeath, "15.20")] 
				CumCases_10yBands[, "20.30"] <- AgeResults[, paste0(CaseOrDeath, "20.25")] + AgeResults[, paste0(CaseOrDeath, "25.30")] 
				CumCases_10yBands[, "30.40"] <- AgeResults[, paste0(CaseOrDeath, "30.35")] + AgeResults[, paste0(CaseOrDeath, "35.40")] 
				CumCases_10yBands[, "40.50"] <- AgeResults[, paste0(CaseOrDeath, "40.45")] + AgeResults[, paste0(CaseOrDeath, "45.50")] 
				CumCases_10yBands[, "50.60"] <- AgeResults[, paste0(CaseOrDeath, "50.55")] + AgeResults[, paste0(CaseOrDeath, "55.60")] 
				CumCases_10yBands[, "60.70"] <- AgeResults[, paste0(CaseOrDeath, "60.65")] + AgeResults[, paste0(CaseOrDeath, "65.70")] 
				CumCases_10yBands[, "70.80"] <- AgeResults[, paste0(CaseOrDeath, "70.75")] + AgeResults[, paste0(CaseOrDeath, "75.80")] 
				CumCases_10yBands[, "80.85"] <- AgeResults[, paste0(CaseOrDeath, "80.85")] 
				
				NumAgeBands = 9
				NumWeeks 	= ceiling(dim(AgeResults)[1]/7)
				
				FirstAgeBand = TRUE 
				DaysToConsider = seq(1, dim(AgeResults)[1], by = 7)
				for (ageband in 1:NumAgeBands)  
				{
					Vec = na.omit(diff(CumCases_10yBands[, ageband], lag = 7)[DaysToConsider])
					if (FirstAgeBand) 
					{
						WeeklyInc = Vec 
						FirstAgeBand = FALSE
						
					}	else WeeklyInc = cbind(WeeklyInc, Vec)
				}
				colnames(WeeklyInc) = colnames(CumCases_10yBands)
				
				### cut off low incidence weeks.
				PlotWindow <- CalcPlotWindow(WeeklyInc, Threshold = 1)
				WeeklyInc <- WeeklyInc[PlotWindow,]
				Dates <- Day_0 + PlotWindow * 7
				Dates <- sub("2019-", "", Dates)
				Dates <- sub("2020-", "", Dates)
				Dates <- sub("2021-", "", Dates)
				rownames(WeeklyInc) <- as.character(sub("2020-", "", Dates))
				
				LWD <- 4
				png(file = file.path(ScenarioPlotDir, 
								paste0("WeeklyIncidence", CaseOrDeath_long, "byAgeGroup.png")), 
						res = PNG_res, units = "in", width = 10, height = 10)
				barplot(t(WeeklyInc), col = rev(bpy.colors(9)), #xaxt = "n",  
						las = 2, main = paste0(Scenario, "\nWeekly ", CaseOrDeath_long, " by age"))
				legend("topright", legend = colnames(CumCases_10yBands), col = rev(bpy.colors(9)), lty = NA, pch = 15, lwd = LWD, cex = 1)
				dev.off()
			}
		}
		
		if (OUTPUT_AdminUnitInc) 
		{
			SeverityVariable = "incCritical"
			cat(paste0("Admin Level Plots: "))
			
			AdUnitFileName <- paste0(cur_path, Scenario, ".avNE.severity.adunit.xls" ) 
			if (!file.exists(AdUnitFileName) ) next else
				AdUnitResults <- read.delim(AdUnitFileName, header = TRUE)
			## clean
			for (col in dim(AdUnitResults)[2]:1) if (all(is.na(AdUnitResults[,col]))) AdUnitResults[,col] = NULL
			Dates 				<- Day_0 + 0:(dim(AdUnitResults)[1]-1)
			
			### create subfolder so more organised. Don't want empty folders so wastefully include this within if statements.
			ScenarioPlotDir = makeScenarioPlotDir(PlotOutputSubDir, Scenario)
			
			for (SeverityVariable in SeverityVariablesToPlot)
			{
				InfVariableStringLong 	= GetVerboseVariableString	(SeverityVariable)
				IncPrevOrCumIncString 	= GetIncPrevOrCumIncString	(SeverityVariable)
				
				## Get admin unit names
				VariablePattern 	= "cumCritRecov_"
				IndicesOfPattern 	= grep(VariablePattern, colnames(AdUnitResults))
				NumAdminUnits 		= length(IndicesOfPattern)
				AdminUnitNames 		= sub(VariablePattern, "", colnames(AdUnitResults)[IndicesOfPattern])
				AdminUnitNames 		= gsub("_", " ", AdminUnitNames)
				
				### select only cols of relevant variables (infection/ case/ Mild/ Critical etc.)
				MinAdUnitColumn = min(grep(SeverityVariable, colnames(AdUnitResults)))
				AdUnitColumns	= colnames(AdUnitResults)[MinAdUnitColumn - 1 + 1:NumAdminUnits]
				
				cat(paste0(InfVariableStringLong, " ", IncPrevOrCumIncString,": "))
				
				MaxNumAdUnitsPerPlot = 20
				NumAdminPlots		 = ceiling(NumAdminUnits / MaxNumAdUnitsPerPlot)
				
				for (AdPlot in 1:NumAdminPlots)
				{
					MinAdUnit 	= (AdPlot - 1) * MaxNumAdUnitsPerPlot + 1
					MaxAdUnit 	= (AdPlot * MaxNumAdUnitsPerPlot)
					if (MaxAdUnit > NumAdminUnits) MaxAdUnit = NumAdminUnits
					Cols 		= bpy.colors(MaxAdUnit - MinAdUnit + 1)
					cat(paste0(" ", MinAdUnit, "-", MaxAdUnit, "\t"))
					
					AdUnitResults_tmp 	= as.matrix(AdUnitResults[,AdUnitColumns][, MinAdUnit:MaxAdUnit])
					Threshold = 1; PlotFromTo = NULL
					
					PlotWindow			= CalcPlotWindow(AdUnitResults_tmp, Threshold = Threshold, PlotFromTo = PlotFromTo)
					AdUnitResults_tmp 	= as.matrix(AdUnitResults_tmp[PlotWindow,])
					Dates_tmp			= Dates[PlotWindow	]
					maxYaxis 			= max(AdUnitResults_tmp)
					
					pngFileName = paste0(ScenarioPlotDir, "AdUnits_", 
							InfVariableStringLong, "_", IncPrevOrCumIncString, "_", MinAdUnit, "_", MaxAdUnit,  ".png")
					
					if (any(AdUnitResults_tmp > 0)) ### don't make empty plots
					{
						png(file = pngFileName, res = PNG_res, units = "in", width = 10, height = 10)
						YLAB = paste(InfVariableStringLong, IncPrevOrCumIncString)
						
						plot(Dates_tmp, AdUnitResults_tmp[,1], 
								main = "", 
								#main = paste0(Scenario, "\n", 
								#		InfVariableStringLong, " ", IncPrevOrCumIncString, " by admin unit (level ", AdLevelToPlot, ")"),
								col = Cols[1], type = "l", lwd = LWD, xaxt = "n", xlab = "",
								ylab = "", ylim = c(0, maxYaxis), cex.axis = 1.7, cex.lab = 1.7)
						mtext(YLAB, side = 2, cex = 2, line = 2.6)
						if (dim(AdUnitResults_tmp)[2] > 1)
							for (Index in 2:(MaxAdUnit - MinAdUnit))
								lines(Dates_tmp, AdUnitResults_tmp[,Index], col = Cols[Index], lwd = LWD)
						axis.Date(1, at = seq(min(Dates_tmp), max(Dates_tmp), by = "months"), 
								format = "%b", las = 2, cex.axis = 1.8)
						
						legend ("topleft", legend = AdminUnitNames, col = Cols, lty = 1, pch = NA, lwd = LWD, cex = 1.3)
						dev.off()
						
					} else cat(paste0("Didn't output ", pngFileName, " as only zeros"))
				}
				cat("\n")
			}
		}
	}
}  

warnings()