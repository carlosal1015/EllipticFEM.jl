
function equationDataAssemble(path::ASCIIString)
	
	include(path)
	
	equationDataGammasInitializeParametrizations(equationData)
	
	return equationData
end

