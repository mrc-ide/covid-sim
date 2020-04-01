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
Day_0 <- as.Date("2020-01-06")
NumDaysToPlot <- 186

# resolution of output files
PNG_res <- 300

# Which plots should be shown?
OUTPUT_AdminUnitInc  <- TRUE
OUTPUT_SEIR         <- TRUE
OUTPUT_SevInc       <- TRUE
OUTPUT_AgeInc       <- TRUE

# where can spatial data be found?
AdUnitLabelFilename <- "Y:/Europe_Spatial_Data/Fr_It/pop_fr_it_index.txt"

# default console output width
defaultOutputWidth <- 172L






#### Set up ####

options(width = defaultOutputWidth)

if (!require(sp)) {
  do_install <- readline("Need to install package `sp`. Continue (Y/n)? ")
  if (!(do_install == "" | tolower(do_install) == "y")) stop("Cannot continue script without dependencies.")
  install.packages("sp")
  library(sp)
}

# parse command line arguments
args = commandArgs(trailingOnly=TRUE)
BaseFolder <- if (length(args > 0)) args[1] else DefaultBaseFolder
PlotOutputDir <- file.path(BaseFolder, PlotOutputDirRelativeToBaseFolder)

#### Hardcoded values about the simulation API ####

# variables put out in the severity file
SeverityVariables <- c("Mild", "ILI", "SARI", "Critical", "incMild", "incILI", "incSARI", "incCritical", "incDeath")
# severity output from c labelled as incidence of death, but it's actually cumulative deaths. correct here for now but will later fix in c code. 
CORRECT_INC_DEATH <- TRUE #
# filenames to look for
Pattern 			    <- ".avNE.adunitVar.xls"



#### read admin unit data ####

# wunderalbert: I don't have any example of admin unit data, so I'm not touching this except for formatting. 
# However, it looks like there are some hardcoded parameters in here that should be moved.
# Anyways, since we don't have that file, I'll wrap this in a try
if (OUTPUT_AdminUnitInc && is(try({
  str(AdUnitLabelData)
  colnames(AdUnitLabelData) <- c("Number", "Country", "Ad1", "Ad2")
  AdUnitLabelData <- AdUnitLabelData[AdUnitLabelData$Country == "Italy", ]
  levels (AdUnitLabelData$Country) <- droplevels(AdUnitLabelData$Country)
  NumAdmin1s <- length(unique(AdUnitLabelData$Ad1))
  NumAdmin2s <- length(unique(AdUnitLabelData$Ad2))
  MaxNumAdmin_1_Digits <- 2
  MaxNumAdmin_2_Digits <- 2
  CountryCode <- 13
  
  
  AdLevelToPlot = 1
}), "try-error")) {
  warning("Could not find admin unit label data. Can't do OUTPUT_AdminUnitInc")
  OUTPUT_AdminUnitInc <- FALSE
}





#### Define processing functions ####

CalcPlotWindow <- function(EpiCurves, Threshold = 0, PlotFromTo = NULL) {
  if (!is.null(PlotFromTo) && length(PlotFromTo) != 2) stop("PlotFromTo must be character of length 2")
  # set default
  MinTime <- 1
  MaxTime <- dim(EpiCurves)[1]
  PlotWindowIndices <- MinTime:MaxTime  
  if (is.null(PlotFromTo)) {
    if (Threshold > 0) {
      MinTime = 1
      MaxTime = dim(EpiCurves)[1]
      for (timestep in 1:dim(EpiCurves)[1]) if (any(EpiCurves[timestep, ] > Threshold)) {  MinTime = timestep; break } 
      for (timestep in dim(EpiCurves)[1]:1) if (any(EpiCurves[timestep, ] > Threshold)) {  MaxTime = timestep; break } 
      if ((MaxTime - MinTime) <= 0)  {
        ## reset in weird cases where exactly one time step is left. Breaks everything. Think of proper solution later
        MinTime <- 1
        MaxTime <- dim(EpiCurves)[1];
      }
      if (MinTime == 1 & MaxTime == dim(EpiCurves)[1]) 
        warning("CalcPlotWindow: cannot subset, try new threshold")
      PlotWindowIndices <- MinTime:MaxTime    
    }
  } else {
    PlotWindowIndices <- as.numeric(as.Date(PlotFromTo[1]) - Day_0):as.numeric(as.Date(PlotFromTo[2]) - Day_0)
  }
  
  return(PlotWindowIndices)
}

makeScenarioPlotDir <- function(PlotOutputSubDir, Scenario) {
  # make a folder to put plots in and return that folder
  ScenarioPlotDir = file.path(PlotOutputSubDir, Scenario)
  dir.create(PlotOutputDir, showWarnings = FALSE)
  dir.create(PlotOutputSubDir, showWarnings = FALSE)
  dir.create(ScenarioPlotDir, showWarnings = FALSE)
  stopifnot(dir.exists(ScenarioPlotDir))
  ScenarioPlotDir
}



#### Main processing loop ####

for (SubFolder in VariousSubFolders) {
  cur_path       <- file.path(BaseFolder, SubFolder)
  FilesToCheck   <- list.files(path = cur_path, pattern = Pattern)
  Scenarios     <- sub(Pattern, "", FilesToCheck) 
  
  PlotOutputSubDir <- file.path(PlotOutputDir, SubFolder)
  dir.create(PlotOutputSubDir, showWarnings = FALSE)
  
  cat(paste0("Processing Subfolder: ", SubFolder, " (", length(Scenarios), " scenarios)"))
  for (Scenario in Scenarios) {
    cat(paste0("\nScenario ", Scenario, ": "))
    
    if (OUTPUT_SEIR) {
      cat(paste0("SEIR, "))
      
      NonSeverityFileName    <- file.path(cur_path, paste0(Scenario, ".avNE.xls"))
      if (!file.exists(NonSeverityFileName)) next else
        NonSeverityResults   <- read.delim(NonSeverityFileName, header = TRUE)
      ### create subfolder so more organised. Don't want empty folders so wastefully include this within if statements.
      ScenarioPlotDir = makeScenarioPlotDir(PlotOutputSubDir, Scenario)
      Dates       <- Day_0 + 0:(dim(NonSeverityResults)[1]-1)
      
      ##### MAKE SEIR plot
      png(file = file.path(ScenarioPlotDir, "SimpleSEIR.png"), res = PNG_res, units = "in", width = 7, height = 7)
      LWD <- 4
      plot(Dates, NonSeverityResults$S, main = paste0(Scenario, "\nSEIR/D"),
           col = "green", type = "l", lwd = LWD, ylab = "prevalence", xlab = "", xaxt = "n", 
           ylim = c(0, max(NonSeverityResults$S)))
      lines(Dates, NonSeverityResults$L, col = "orange"  , lwd = LWD)
      lines(Dates, NonSeverityResults$I, col = "red"    , lwd = LWD)
      lines(Dates, NonSeverityResults$R, col = "yellow"  , lwd = LWD)
      lines(Dates, NonSeverityResults$D, col = "black"  , lwd = LWD)
      axis.Date(1, at=seq(min(Dates), max(Dates), by = "months"), format="%b", las = 2)
      legend ("topright", legend = c("Susceptible", "Latent", "Infectious", "Recovered", "Dead"), col = c("green", "orange", "red", "yellow", "black"), lty = 1, pch = NA, lwd = LWD, cex = 1)
      dev.off()
    }
    
    if (OUTPUT_SevInc) {
      
      cat(paste0("SeverityInc, "))
      SeverityFileName <- file.path(cur_path, paste0(Scenario, ".avNE.severity.xls"))
      if (!file.exists(SeverityFileName  )) next else
        SeverityResults <- read.delim(SeverityFileName, header = TRUE)
      ## clean
      for (col in dim(SeverityResults)[2]:1) if (all(is.na(SeverityResults[,col]))) SeverityResults[,col] = NULL
      if (CORRECT_INC_DEATH) SeverityResults$incDeath = c(0, diff(SeverityResults$incDeath))
      
      ### create subfolder so more organised. Don't want empty folders so wastefully include this within if statements.
      ScenarioPlotDir <- makeScenarioPlotDir(PlotOutputSubDir, Scenario)
      
      Dates         <- Day_0 + 0:(dim(SeverityResults)[1]-1)
      
      ### cut off low incidence weeks.
      PlotWindow      <- CalcPlotWindow(SeverityResults, Threshold = 1)
      SeverityResults <- SeverityResults[PlotWindow,]
      Dates            <- Dates[PlotWindow]
      
      ##### MAKE incidence by severity plot
      LWD   <- 4
      png(file = file.path(ScenarioPlotDir, "SeverityIncidence.png"), res = PNG_res, units = "in", width = 7, height = 7)
      plot(Dates, SeverityResults$incMild, main = paste0(Scenario, "\nincidence by disease severity"),
           col = "pink", type = "l", lwd = LWD, ylab = "incidence", xlab = "", xaxt = "n", 
           ylim = c(0, max(SeverityResults)))
      lines(Dates, SeverityResults$incILI    , col = "palevioletred1", lwd = LWD)
      lines(Dates, SeverityResults$incSARI  , col = "orangered"    , lwd = LWD)
      lines(Dates, SeverityResults$incCritical, col = "red"      , lwd = LWD)
      lines(Dates, SeverityResults$incDeath  , col = "black"      , lwd = LWD)
      #axis.Date(1, at=seq(min(Dates), max(Dates), by = "months"), format="%m-%Y", las = 2)
      axis.Date(1, at=seq(min(Dates), max(Dates), by = "months"), format="%b", las = 2)
      legend ("topright", legend = c("Mild", "ILI", "SARI", "Critical", "Dead"), col = c("pink", "palevioletred1", "orangered", "red", "black"), lty = 1, pch = NA, lwd = LWD, cex = 1)
      dev.off()
    }
    
    if (OUTPUT_AgeInc) {
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
      for (CaseOrDeath in CasesOrDeaths) {
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
        
        NumWeeks  <- ceiling(dim(AgeResults)[1]/7)
        WeeklyInc <- matrix(nrow = NumWeeks, ncol = 9)
        colnames(WeeklyInc) <- colnames(CumCases_10yBands)
        
        for (col in 1:9)  WeeklyInc[, col] <- diff(CumCases_10yBands[, col], lag = 7)[seq(7, dim(AgeResults)[1], by = 7)]
        if (all(is.na(WeeklyInc[NumWeeks, ] ))) WeeklyInc[NumWeeks, ] = 0
        
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
    
    if (OUTPUT_AdminUnitInc) {
      InfOrCaseString <- "I"
      InfOrCaseString <- "incCritical"
      for (InfOrCaseString in c("I", "C", "incSARI", "incCritical")) {
        Suffix <-
          if (! InfOrCaseString %in% SeverityVariables) {
            ".avNE.adunit.xls"
          } else {
            ".avNE.severity.adunit.xls" 
          }
        
        InfVariableStringLong <- 
          switch(InfOrCaseString, 
                 I = "Infection",
                 C = "Case",
                 InfOrCaseString)
        
        AdUnitFileName <- file.name(cur_path, paste0(Scenario, Suffix))
        if (!file.exists(AdUnitFileName) ) next else
          AdUnitResults <- read.delim(AdUnitFileName, header = TRUE)
        
        ## clean
        for (col in dim(AdUnitResults)[2]:1) if (all(is.na(AdUnitResults[,col]))) AdUnitResults[,col] <- NULL
        
        NumAdminUnits <-
          if (InfOrCaseString %in% c("I", "C")) {
            min(grep("C", colnames(AdUnitResults))) - min(grep("I", colnames(AdUnitResults))) 
          } else {
            length(grep("varincSARI_", colnames(AdUnitResults)))
          }
        
        ### select only cols of relevant variables (infection/ case/ Mild/ Critical etc.)
        WhichAdUnits   <- 1:NumAdminUnits ### change if interested in particular subset of admin units
        NamesAdminUnits <- paste0("ad", WhichAdUnits)
        MinAdUnitColumn <- min(grep(InfOrCaseString, colnames(AdUnitResults)))
        AdUnitColumns  <- colnames(AdUnitResults)[MinAdUnitColumn - 1 + WhichAdUnits]
        AdUnitResults <- AdUnitResults[,AdUnitColumns]
        #colnames(AdUnitResults)
        
        ## get num admin units (of whatever level from data) 
        if (file.exists(AdUnitLabelFilename)) {
          NumAdminUnits <- switch(AdLevelToPlot,
                                  `1` = NumAdmin1s,
                                  `1` =  NumAdmin2s)
        }
        cat(paste0(InfVariableStringLong, " AdminUnitInc, "))
        
        ### create subfolder so more organised. Don't want empty folders so wastefully include this within if statements.
        ScenarioPlotDir = makeScenarioPlotDir(PlotOutputSubDir, Scenario)
        
        Dates             <- Day_0 + 0:(dim(AdUnitResults)[1]-1)
        PlotWindowIndices <- 1:length(Dates)  
        
        if (AdLevelToPlot == 1) {
          ### need to aggregate all columns correponding to region 1, region 2 etc. 
          AdUnitResults_Agg <- matrix(nrow = dim(AdUnitResults)[1], ncol = NumAdmin1s)
          ChoppedColnames <- substr(colnames(AdUnitResults), 1, nchar(colnames(AdUnitResults)) - MaxNumAdmin_2_Digits)
          
          region <- 1
          for (region in 1:NumAdmin1s) {
            ### get columns for this region
            ThisRegionColumns <- grep(unique(ChoppedColnames)[region], colnames(AdUnitResults))
            if (length(ThisRegionColumns) < 2) {
              AdUnitResults_Agg[, region] <- AdUnitResults[, ThisRegionColumns] 
            } else {
              AdUnitResults_Agg[, region] <- rowSums(AdUnitResults[, ThisRegionColumns])
            }
            
          }
        } else AdUnitResults_Agg = AdUnitResults
        
        ## say (these will definitely be country-specific so fix later).
        if (AdLevelToPlot == 1) MaxNumAdUnitsPerPlot = 20 else
          if (AdLevelToPlot == 2) MaxNumAdUnitsPerPlot = 10   
        NumAdminPlots     = ceiling(NumAdminUnits / MaxNumAdUnitsPerPlot)
        
        for (AdPlot in 1:NumAdminPlots) {
          MinAdUnit <- (AdPlot - 1) * MaxNumAdUnitsPerPlot + 1
          MaxAdUnit <- (AdPlot * MaxNumAdUnitsPerPlot)
          if (MaxAdUnit > NumAdminUnits) MaxAdUnit = NumAdminUnits
          Cols <- bpy.colors(MaxAdUnit - MinAdUnit + 1)
          
          AdUnitResults_tmp <- as.matrix(AdUnitResults_Agg[, MinAdUnit:MaxAdUnit])
          Threshold <- 1
          #            if (MaxAdUnit - MinAdUnit == 0) 
          #              PlotWindow = which(AdUnitResults_tmp >= Threshold) else
          PlotWindow <- CalcPlotWindow(AdUnitResults_tmp, Threshold = Threshold)
          AdUnitResults_tmp <- as.matrix(AdUnitResults_tmp[PlotWindow,])
          Dates_tmp <- Dates[PlotWindow  ]
          maxYaxis <- max(AdUnitResults_tmp)
          
          LWD <- 4
          pngFileName <- file.path(ScenarioPlotDir, 
                                   paste0("Inc_AdUnits_Level_", AdLevelToPlot, "_", 
                                          InfVariableStringLong, "_", MinAdUnit, "_", MaxAdUnit,  ".png"))
          
          ### don't make empty plots
          if (any(AdUnitResults_tmp > 0)) {
            png(file = pngFileName, res = PNG_res, units = "in", width = 10, height = 10)
            plot(Dates_tmp, AdUnitResults_tmp[,1], 
                 main = paste0(Scenario, "\n", 
                               InfVariableStringLong, " incidence by admin unit (level ", AdLevelToPlot, ")"),
                 col = Cols[1], type = "l", lwd = LWD, xaxt = "n", xlab = "",
                 ylab = paste0(InfVariableStringLong, " incidence"), ylim = c(0, maxYaxis))
            if (dim(AdUnitResults_tmp)[2] > 1)
              for (Index in 2:(MaxAdUnit - MinAdUnit))
                lines(Dates_tmp, AdUnitResults_tmp[,Index], col = Cols[Index], lwd = LWD)
            axis.Date(1, at = seq(min(Dates_tmp), max(Dates_tmp), by = "months"), 
                      format = "%b", las = 2)
            
            LegLabels <- switch(AdLevelToPlot,
                                `1` = unique(AdUnitLabelData$Ad1)[MinAdUnit:MaxAdUnit],
                                `2` = AdUnitLabelData$Ad2[MinAdUnit:MaxAdUnit],
                                # default
                                LegLabels <- paste0("AdUnit_", MinAdUnit:MaxAdUnit))
            legend ("topright", legend = LegLabels, col = Cols, lty = 1, pch = NA, lwd = LWD, cex = 1)
            dev.off()
            
          } else cat(paste0("Didn't output ", pngFileName, " as only zeros"))
        }
      }
    }
  }
}  

warnings()