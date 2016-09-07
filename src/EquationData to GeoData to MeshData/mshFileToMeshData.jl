
function mshFileToMeshData(path::ASCIIString)
	
	filehandle = open(path,"r")
	content = readdlm(filehandle)
	close(filehandle)
	
	i = findfirst(content,"\$Nodes")+2
	nodesNr = content[i-1,1]
	nodes = convert(Array{Float64,2},content[i:1:(i+nodesNr-1),2:1:3]')
	
	i = findfirst(content,"\$Elements")+2
	edgesAndElementsNr = content[i-1,1]
	rowIndicesEdges = findin(content[i:1:(i+edgesAndElementsNr-1),2],1)
	edges = convert(Array{Int64,2},content[i+rowIndicesEdges-1,[6,7]]')
	edgesParentGeoEdgesIds = convert(Array{Int64},content[i+rowIndicesEdges-1,5])
	edgesParentDomainsTypes = []
	edgesParentDomainsIds = []
	edgesParentGammasIds = []
	edgesParentGammasParameterRanges = []'
	edgesVolumes = []
	
	rowIndicesElements = findin(content[i:1:(i+edgesAndElementsNr-1),2],2)
	elements = convert(Array{Int64,2},content[i+rowIndicesElements-1,[6,7,8]]')
	elementsParentSubdomainsIntersectionsIds = convert(Array{Int64},content[i+rowIndicesElements-1,5])
	elementsParentDomainsTypes = []
	elementsParentDomainsIds = []
	elementsDetNablaFTs = []
	elementsNablaFTInv11s = []
	elementsNablaFTInv21s = []
	elementsNablaFTInv12s = []
	elementsNablaFTInv22s = []
	elementsVolumes = []
	elementsCentroids = []'
	elementsMeshsizes = []
	
	iStart = findfirst(content,"\$Periodic")
	if iStart != 0
		
		iStart += 1
		iEnd = findfirst(content,"\$EndPeriodic")-1
		
		rowIndices = (iStart-1) + findin(map(i->sum(map(x->typeof(x)==Int64,content[i,:])),iStart:1:iEnd),2)
		periodicNodes = convert(Array{Int64,2},content[rowIndices,[2,1]]')
		
	else
		periodicNodes = []'
	end
	
	meshData = MeshData(nodes,
		edges,edgesParentGeoEdgesIds,edgesParentDomainsTypes,edgesParentDomainsIds,edgesParentGammasIds,edgesParentGammasParameterRanges,edgesVolumes,
		elements,elementsParentSubdomainsIntersectionsIds,elementsParentDomainsTypes,elementsParentDomainsIds,elementsDetNablaFTs,elementsNablaFTInv11s,elementsNablaFTInv21s,elementsNablaFTInv12s,elementsNablaFTInv22s,elementsVolumes,elementsCentroids,elementsMeshsizes,
		periodicNodes)
	
	return meshData
end

