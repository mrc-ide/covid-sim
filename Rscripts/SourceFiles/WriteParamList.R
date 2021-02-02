WriteParamList = function(ParamList, OutputDir, OutputFileName, PrintToConsole = FALSE)
{
	ListString = ""
	for (ParamNum in 1:length(ParamList))
	{
		ParamNameString = paste0("[", names(ParamList)[ParamNum], "]")
		if (PrintToConsole) cat(paste0(ParamNameString, "\n"))
		ListString 		= paste0(ListString, ParamNameString, "\n")
		if (class (ParamList[[ParamNum]]) %in% c("matrix", "data.frame"))
		{
			ParamValueString = ""
			for (row in 1:dim(ParamList[[ParamNum]])[1])
				ParamValueString = paste0(ParamValueString, paste0(ParamList[[ParamNum]][row,], collapse = " "), " \n")
			
		} else 	ParamValueString = paste0(paste0(ParamList[[ParamNum]], collapse = " "), "\n")
		if (PrintToConsole) cat(paste0(ParamValueString, "\n"))
		
		ListString = paste0(ListString, ParamValueString, "\n")
	}
	
	write.table(ListString, file = file.path(OutputDir, OutputFileName), row.names = F, col.names = F, quote = F, sep = "\t")
}
