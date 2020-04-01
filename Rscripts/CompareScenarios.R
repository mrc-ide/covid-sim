#!/usr/bin/env Rscript
# This script creates plots to compare different scenarios.
# Command line arguments: 
#    folder: the relative folder of where to find output files.
# Will create plots in the "Plots/Comparsons" subfolder of folder

#### Purpose: make plots that Neil outputs in his spreadsheets without having to use spreadsheets. ####

# Requires a set of input runs in a base folder
# Puts out a set of plots in the /Plots/Comparisons subfolder of the base folder
# Edit the Parameters section to customize

# code sent by d.laydon@imperial.ac.uk
# first refactoring pass wunderalbert@github.com




#### Parameters ####

# May later refactor to take inputs from seperate file or command line

# Where will the run data be located?
DefaultInFolder <- "."

# Where will the plots be outputted to?
PlotOutputDirRelativeToBaseFolder <- c("Plots", "Comparisons")

# Date range
Day_0 <- as.Date("2020-01-06")
NumDaysToPlot <- 186

# resolution of output files
PNG_res <- 300

# list the scenarios to compare (or set NULL to take all available)
Scenarios <- NULL

# Which comparison should we run?
# The name will determine title and filename.
ComparisonName 	 <- "Number ICU beds"
# The severity variable is the variable to compare
SeverityVariable <- "incCritical"

# consider the following R value:
R <- 2.2

# default console output width
defaultOutputWidth <- 172L


#### Set up ####

options(width = defaultOutputWidth)

if (!require(sp)){
  do_install <- readline("Need to install package `sp`. Continue (Y/n)? ")
  if (!(do_install == "" | tolower(do_install) == "y")) stop("Cannot continue script without dependencies.")
  install.packages("sp")
  library(sp)
}

# parse command line arguments
args = commandArgs(trailingOnly=TRUE)
InFolder <- if (length(args > 0)) args[1] else DefaultInFolder
CompPlotDir <- do.call(file.path, c(list(InFolder), as.list(PlotOutputDirRelativeToBaseFolder)))

# Create comparison plot directory
dir.create(CompPlotDir, showWarnings = FALSE, recursive = TRUE)
stopifnot(dir.exists(CompPlotDir))


#### Hardcoded values about the simulation API ####

# variables put out in the severity file
SeverityVariables <- c("Mild", "ILI", "SARI", "Critical", "incMild", "incILI", "incSARI", "incCritical", "incDeath")
# severity output from c labelled as incidence of death, but it's actually cumulative deaths. correct here for now but will later fix in c code. 
CORRECT_INC_DEATH <- TRUE #
# filenames to look for
Pattern 			    <- ".avNE.adunitVar.xls"

# where do we find the info we want?
Suffix <-
  if (SeverityVariable %in% SeverityVariables){
    ".avNE.severity.xls" 
  } else {
    ".avNE.xls"
  }

#### Define processing functions ####

CalcPlotWindow <- function(EpiCurves, Threshold = 0, PlotFromTo = NULL) {
  if (!is.null(PlotFromTo) && length(PlotFromTo) != 2) stop("PlotFromTo must be character of length 2")
  # set default
  MinTime <- 1
  MaxTime <- dim(EpiCurves)[1]
  PlotWindow <- MinTime:MaxTime	
  if (is.null(PlotFromTo)) {
    if (Threshold > 0) {
      MinTime <- 1
      MaxTime <- dim(EpiCurves)[1]
      for (timestep in 1:dim(EpiCurves)[1]) if (any(EpiCurves[timestep, ] > Threshold)) {	MinTime = timestep; break } 
      for (timestep in dim(EpiCurves)[1]:1) if (any(EpiCurves[timestep, ] > Threshold)) {	MaxTime = timestep; break } 
      if ((MaxTime - MinTime) <= 0)  ## reset in weird cases where exactly one time step is left. Breaks everything. Think of proper solution later.
      {
        MinTime <- 1
        MaxTime <- dim(EpiCurves)[1];
      }
      PlotWindow <- MinTime:MaxTime		
    }
  } else {
    PlotWindow <- as.numeric(as.Date(PlotFromTo[1]) - Day_0):as.numeric(as.Date(PlotFromTo[2]) - Day_0)
  }
  
  return(PlotWindow)
}

#### Define Plotting Functions ####

CompareEpiCurves <- function(Filenames, SeverityVariable = "incI", ComparisonName = "AllRuns", Names, 
                             Threshold = 0, PlotFromTo = NULL, YLAB = SeverityVariable, LWD = 4, LegendPosition = "topright"){
  if (!is.null(PlotFromTo) && length(PlotFromTo) != 2) stop("PlotFromTo must be of length 2")
  #if (all(SeverityVariable != SeverityVariables)) Suffix = ".avNE.xls" 	else Suffix = ".avNE.severity.xls" 
  
  Index 	<- 1
  Flag    <- TRUE
  for (Index in 1:length(Filenames)) {
    FileName	<- Filenames[Index]
    cat(paste0(Index, "/", length(Filenames), ": Scenario ", sub(InFolder, "", FileName), "\n"))
    Results		<- read.delim(FileName, header = TRUE)
    if (Flag) {
      EpiCurves	<- matrix(nrow = dim(Results)[1], ncol = length(Filenames))
      Flag <- FALSE
    }	
    EpiCurves[, Index]  <- Results[, SeverityVariable]
  }
  
  Cols 		    	<- bpy.colors(length(Filenames))
  Dates 			  <- Day_0 + 0:(dim(Results)[1]-1)
  PlotWindow		<- CalcPlotWindow(EpiCurves, Threshold = Threshold, PlotFromTo = PlotFromTo)
  EpiCurves 		<- EpiCurves[PlotWindow,]
  Dates 		  	<- Dates[PlotWindow	]
  
  maxY <- max(EpiCurves)
  plot(Dates, EpiCurves[,1],
       xaxt = "n", xlab = "",
       main = ComparisonName,
       col = Cols[1], type = "l", lwd = LWD, 
       ylab = YLAB, 
       ylim = c(0, maxY))
  for (Index in 2:length(Filenames))
    lines(Dates, EpiCurves[,Index], col = Cols[Index]	, lwd = LWD)
  axis.Date(1, at = seq(min(Dates), max(Dates), by = "months"), format = "%b", las = 2)
  
  legend (LegendPosition, legend = Names, col = Cols, lty = 1, pch = NA, lwd = LWD, cex = 1)
}


#### Prettify Names for scenarios ####

GetVerboseScenarioName_Single 	= function(AbbreiviatedScenarioName) {
  VerboseScenario <- ""
  
  if (length(grep("NoInt"	, AbbreiviatedScenarioName)) > 0) VerboseScenario = paste0(VerboseScenario, "NoIntervention_"	)
  if (length(grep("MG"	, AbbreiviatedScenarioName)) > 0) VerboseScenario = paste0(VerboseScenario, "MassGatherings_"	)
  if (length(grep("PC"	, AbbreiviatedScenarioName)) > 0) VerboseScenario = paste0(VerboseScenario, "PlaceClosure_"		)
  if (length(grep("CI"	, AbbreiviatedScenarioName)) > 0) VerboseScenario = paste0(VerboseScenario, "CaseIsolation_"	)
  if (length(grep("HQ"	, AbbreiviatedScenarioName)) > 0) VerboseScenario = paste0(VerboseScenario, "HomeQuarantine_"	)
  if (length(grep("SDOL"	, AbbreiviatedScenarioName)) > 0) VerboseScenario = paste0(VerboseScenario, "SocialDistOlder_"	) else
    if (length(grep("SDO"	, AbbreiviatedScenarioName)) > 0) VerboseScenario = paste0(VerboseScenario, "SocialDistOlder_"	) else
      if (length(grep("SD"	, AbbreiviatedScenarioName)) > 0) VerboseScenario = paste0(VerboseScenario, "SocialDist_"		)
      if (length(grep("SC"	, AbbreiviatedScenarioName)) > 0) VerboseScenario = paste0(VerboseScenario, "SchoolClosure_"	)
      
      ### Chop off final "_"
      
      while (substr(VerboseScenario, nchar(VerboseScenario), nchar(VerboseScenario)) == "_") 
        VerboseScenario = substr(VerboseScenario, 1, nchar(VerboseScenario) - 1)
      return(VerboseScenario)
}

GetVerboseScenarioName_Many 	= function(AbbreiviatedScenarioNames) {
  unlist(lapply(AbbreiviatedScenarioNames, GetVerboseScenarioName_Single))
}


#### ==== Get list of scenarios / filenames to compare ####

## this should give every model run in folder (but will be too much info on one plot). 
PatternWithR 		  <- paste0("R0=", R, Pattern)	
FilesToCheck 	    <- list.files(path = InFolder, pattern = PatternWithR)
PossibleScenarios <- sub(PatternWithR, "", FilesToCheck) 

if (0 == length(Scenarios)) {
  # take all possible files
  Scenarios <- PossibleScenarios
  message("Using all scenarios that were found: ", length(PossibleScenarios))
} else {
  message("Found scenarios: ", length(intersect(PossibleScenarios, Scenarios)), " of ", length(Scenarios))
}

Filenames <- file.path(InFolder, paste0(Scenarios, "R0=", R, Suffix))
if (!all(file.exists(Filenames))){
  stop("Aborted plot generation, because the following files were not found: ", Filenames[!file.exists(Filenames)])
}

OutFileName <- file.path(CompPlotDir, paste0(ComparisonName, ".png"))
png(file = OutFileName, res = PNG_res, units = "in", width = 7, height = 7)
CompareEpiCurves(Filenames, SeverityVariable, ComparisonName, 
                 Names = GetVerboseScenarioName_Many(Scenarios), Threshold = 1, YLAB = "ICU beds", LWD = 6)
invisible(dev.off())