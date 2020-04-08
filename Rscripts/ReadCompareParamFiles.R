### Purpose: Extract and compare parameter names and values from (pre-)parameter .txt files
# 
# Author: dlaydon
###############################################################################

rm(list=ls(all=TRUE)) 

#### change as appropriate
CurrentPath 		= "C:/Users/dlaydon/Dropbox (SPH Imperial College)/nCoV_MERS_Covid19/CoVid_19/DJLnCoV_unshared/ParamFiles/"

### issues
# 1) 0 0 0 0 rather than 0\t0\t0\t0: DONE
# 2) Varying number of values after [ParamName] - Mostly done, Very hacky with places and ages. Doesn't work in situations where there are both place and age in parameter description.
# 3) variying number of columns / tabs.
# 4) getting strings into numbers. Should be easy. 
# 5) useless old arrays immediately below ones you want. 
# 6) hashes
#### need to clean NA's from list values. 

ReadParams = function(CurrentPath = CurrentPath, 
		ParamFileName, CleanDuplicateAgeGroups = TRUE, 
		NumAgeGroups = 17, PrintToConsole = FALSE, CleanNAs = TRUE, CleanDuplicatePlaceValues = TRUE, NumPlaceTypes = 4)
{
	## Function cycles through Param Files, and depending on characteristics of line, 
	## either assigns it as a name, single value, or addition to vector.
	Params 	<- read.delim(paste0(CurrentPath, ParamFileName), header = FALSE, stringsAsFactors= FALSE)
	
	ParamList = list(); 
	ListCounter = 0
	LastLineWasAName = FALSE
	
	for (ParamFileLine in 1:dim(Params)[1])
	{
		LineValue = Params[ParamFileLine,]
		class(LineValue)
		
		if (substr(LineValue, 1, 1) == "^" | substr(LineValue, 1, 1) == "=") next
		if (length(grep("\\[", LineValue)) > 0 & length(grep("\\]", LineValue)) > 0) ### then Line denotes a variable name
		{
			ListCounter = ListCounter + 1
			LineValue	= sub("\\[", "", LineValue)
			LineValue	= sub("\\]", "", LineValue)
			ParamName 	= LineValue
			ParamList[[ListCounter]] = NA ### placeholder that will be overwritten in next ParamFileLine
			names(ParamList)[ListCounter] = ParamName
			LastLineWasAName = TRUE
			
		} else 
		{
			if (substr(LineValue, 1, 1) == "#") ### make hashes/#'s into characters, not numerics.  
				ParamList[[ListCounter]] = as.character(LineValue) else
			{
				NumSpacesInString = length(grep(" ", LineValue))
				if (NumSpacesInString == 0) ### if no spaces
				{
					if (LastLineWasAName)  ### i.e. if this is the first value....
						ParamList[[ListCounter]] = as.numeric(as.character(LineValue))	else ## otherwise combine. 
						ParamList[[ListCounter]] = c(ParamList[[ListCounter]], as.numeric(as.character(LineValue)))
					
				} else if (LastLineWasAName) ### this last if statement means you automatically won't consider duplicate 0 0 0 0's. 
					ParamList[[ListCounter]] = as.numeric(unlist(strsplit(LineValue, " ")))
			}
			LastLineWasAName = FALSE
		}
		if (PrintToConsole)
			cat(paste0("Line ", ParamFileLine, " ", Params[ParamFileLine,], "\n"))
	}
	if (CleanNAs) 
		for (param in 1:length(ParamList))
			if (any(is.na(ParamList[[param]])))
				ParamList[[param]] = as.numeric(na.omit(ParamList[[param]]))
	
	if (CleanDuplicateAgeGroups)
		for (param in 1:length(ParamList))
			if ((length(grep("age", names(ParamList)[param])) > 0) 	||	
				(length(grep("Age", names(ParamList)[param])) > 0)	)
				if (length(ParamList[[param]]) > NumAgeGroups) 
					ParamList[[param]] = ParamList[[param]][1:NumAgeGroups]

	if (CleanDuplicatePlaceValues)
		for (param in 1:length(ParamList))
			if ((length(grep("place", names(ParamList)[param])) > 0) 	||	
				(length(grep("Place", names(ParamList)[param])) > 0)	)
				if (length(ParamList[[param]]) > NumPlaceTypes) 
					ParamList[[param]] = ParamList[[param]][1:NumPlaceTypes]
	
	return(ParamList)
}

ComparareLists = function(ParamList1, ParamList2, Name1, Name2, FindDiff = TRUE, CleanNAs = TRUE)
{
	Diffs = list(); 
	Diffs_counter = 0
	P2_checked = rep(FALSE, length(ParamList2))
	p1_index = 1
	p1_index = 12
	for (p1_index in 1:length(ParamList1))
	{
		ParamName = names(ParamList1)[p1_index]
		
		if (ParamName == names(ParamList2)[p1_index]) 
			p2_index = p1_index else
			p2_index = which(names(ParamList2) == ParamName)
		
		if (length(p2_index) > 0) ### if name match found
		{
			P2_checked[p2_index] = TRUE
			### Are values identical? 
			Condition = identical(ParamList1[[p1_index]], ParamList2[[p2_index]])
			
			if (FindDiff) 
				Condition = !Condition ### if looking for differences, negate above condition that looked for similarities
			
			if (Condition)
			{
				Diffs_counter = Diffs_counter + 1
				Diffs[[Diffs_counter]] = ParamList1[[p1_index]]
				
				if (FindDiff) ### only add values twice if looking for differences
				{
					names(Diffs)[Diffs_counter]  = paste(Name1, names(ParamList1)[p1_index])
					
					Diffs_counter = Diffs_counter + 1
					Diffs[[Diffs_counter]] = ParamList2[[p2_index]]
					names(Diffs)[Diffs_counter]  = paste(Name2, names(ParamList2)[p2_index])
					
				} else names(Diffs)[Diffs_counter] = names(ParamList1)[p1_index]
			}
			
		}	else if (FindDiff)	
		{
			Diffs_counter = Diffs_counter + 1
			Diffs[[Diffs_counter]] = ParamList1[[p1_index]]
			names(Diffs)[Diffs_counter]  = paste(Name1, names(ParamList1)[p1_index], "not in", Name2)
		}
		
	}
	##### go through all that haven't been matched. 
	if (FindDiff)
		for (p2_index in 1:length(ParamList2))
			if (!P2_checked[p2_index])
			{
				Diffs_counter = Diffs_counter + 1
				Diffs[[Diffs_counter]] = ParamList2[[p2_index]]
				names(Diffs)[Diffs_counter]  = paste(Name2, names(ParamList2)[p2_index], "not in", Name1)
			}
	if (CleanNAs) 
		for (DiffParam in 1:length(Diffs))
			if (any(is.na(Diffs[[DiffParam]])))
				Diffs[[DiffParam]] = as.numeric(na.omit(Diffs[[DiffParam]]))
	
	return(Diffs)
}

# UK vs USA pre-params
#Diffs = ComparareLists(	
#		ParamList1 = ReadParams(CurrentPath, "preUK_R0=2.0.txt"), Name1 = "UK"	, 
#		ParamList2 = ReadParams(CurrentPath, "preUS_R0=2.0.txt"), Name2 = "USA")
Diffs = ComparareLists(	
		ParamList1 = ReadParams(CurrentPath, "p_NoInt.txt"	), Name1 = "NoInt"	, 
		ParamList2 = ReadParams(CurrentPath, "p_CI.txt"	), Name2 = "CI" 	)
#Diffs = ComparareLists(	
#		ParamList1 = ReadParams(CurrentPath, "p_CI.txt"	), Name1 = "No HH Q", 
#		ParamList2 = ReadParams(CurrentPath, "p_CI_HQ.txt"	), Name2 = "CI_HQ" 	)
#Diffs = ComparareLists(	
#		ParamList1 = ReadParams(CurrentPath, "p_NoInt.txt"	), Name1 = "NoInt"	, 
#		ParamList2 = ReadParams(CurrentPath, "p_MG.txt"	), Name2 = "MG" 	)
#Diffs = ComparareLists(	
#		ParamList1 = ReadParams(CurrentPath, "p_NoInt.txt"	), Name1 = "NoInt", 
#		ParamList2 = ReadParams(CurrentPath, "p_CI_HQ.txt"	), Name2 = "CI_HQ" 	)
#Diffs = ComparareLists(	
#		ParamList1 = ReadParams(CurrentPath, "p_NoInt.txt"	), Name1 = "NoInt"	, 
#		ParamList2 = ReadParams(CurrentPath, "p_CI.txt"	), Name2 = "CI" 	)

Diffs = ComparareLists(	
		ParamList1 = ReadParams(CurrentPath, "p_SD.txt"	), Name1 = "SD", 
		ParamList2 = ReadParams(CurrentPath, "p_SDOL70.txt"), Name2 = "SD > 70" 	)
Similarities = ComparareLists(	FindDiff = FALSE,
		ParamList1 = ReadParams(CurrentPath, "p_SD.txt"	), Name1 = "SD", 
		ParamList2 = ReadParams(CurrentPath, "p_SDOL70.txt"), Name2 = "SD > 70" 	)

# UK vs Italy pre-params
#Diffs = ComparareLists(	 FindDiff = TRUE,
#		ParamList1 = ReadParams(CurrentPath, "preUK_R0=2.0.txt"), Name1 = "UK"	, 
#		ParamList2 = ReadParams(CurrentPath, "preIT_R0=2.0.txt"), Name2 = "Italy")

### example UK vs Italy parameters
#Diffs = ComparareLists(	
#		ParamList1 = ReadParams(CurrentPath, "p_PC7_CI_HQ_SDbad.txt"	), Name1 = "UK SDbad - Neil", 
#		ParamList2 = ReadParams(CurrentPath, "p_Int_Sep07_SD_relaxed50.txt"), Name2 = "IT SD_relaxed50 - Gemma" 	)


#names(ParamList1)
#names(ParamList2)

Diffs
Similarities


