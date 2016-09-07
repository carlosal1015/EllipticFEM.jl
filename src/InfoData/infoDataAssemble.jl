
function infoDataAssemble(processesStopwatches::Array{ProcessStopwatch,1},meshData::MeshData,lseData::LSEData)
	
	for id = 1:1:length(processesStopwatches)
		processesStopwatches[id].id = id
	end
	
	infoData = InfoData(
		maximum(meshData.elementsMeshsizes),
		size(meshData.elements,2),
		size(meshData.nodes,2),
		size(lseData.M,1),
		(sizeof(Float64)+2*sizeof(Int64))*length(lseData.M.nzval)/1000000,
		sizeof(Float64)*size(lseData.M,1)^2/1000000,
		lseData.solverType,
		processesStopwatches)
	
	return infoData
end
