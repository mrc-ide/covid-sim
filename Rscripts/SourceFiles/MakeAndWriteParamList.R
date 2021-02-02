MakeAndWriteParamList = function(OutputDir, OutputFileName, PrintToConsole = FALSE, ...)
{
	ParamList = MakeParamList(...)
	WriteParamList(ParamList, OutputDir, OutputFileName, PrintToConsole)
}

