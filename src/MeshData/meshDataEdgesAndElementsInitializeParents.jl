
function meshDataEdgesInitializeParents(meshData::MeshData,geoData::GeoData)
	
	edgesParentDomainsTypes = Array{ASCIIString}(size(meshData.edges,2))
	edgesParentDomainsIds = Array{Int64}(size(meshData.edges,2))
	edgesParentGammasIds = Array{Int64}(size(meshData.edges,2))
	edgesParentGammasParameterRanges = Array{Float64,2}(3,size(meshData.edges,2))
	for edgeId = 1:1:size(meshData.edges,2)
		
		geoEdge = geoData.geoEdges[meshData.edgesParentGeoEdgesIds[edgeId]]
		
		edgesParentDomainsTypes[edgeId] = geoEdge.parentDomainType
		edgesParentDomainsIds[edgeId] = geoEdge.parentDomainId
		
		edgesParentGammasIds[edgeId] = geoEdge.parentGammaId
		
		edgeNode1Id = meshData.edges[1,edgeId]
		edgeNode2Id = meshData.edges[2,edgeId]
		geoEdgeNode1Id = geoEdge.geoNodesIds[1]
		geoEdgeNode2Id = geoEdge.geoNodesIds[2]
		
		edgeNode1Coords = meshData.nodes[:,edgeNode1Id]
		edgeCentroidCoords = 1/2*(meshData.nodes[:,edgeNode1Id] + meshData.nodes[:,edgeNode2Id])
		edgeNode2Coords = meshData.nodes[:,edgeNode2Id]
		geoEdgeNode1Coords = geoData.geoNodes[geoEdgeNode1Id].coords
		geoEdgeNode2Coords = geoData.geoNodes[geoEdgeNode2Id].coords
		
		quotientEdgeNode1 = norm(edgeNode1Coords-geoEdgeNode1Coords)/norm(geoEdgeNode2Coords-geoEdgeNode1Coords)
		quotientEdgeCentroid = norm(edgeCentroidCoords-geoEdgeNode1Coords)/norm(geoEdgeNode2Coords-geoEdgeNode1Coords)
		quotientEdgeNode2 = norm(edgeNode2Coords-geoEdgeNode1Coords)/norm(geoEdgeNode2Coords-geoEdgeNode1Coords)
		
		edgesParentGammasParameterRanges[:,edgeId] =
			[1-quotientEdgeNode1,1-quotientEdgeCentroid,1-quotientEdgeNode2]*geoEdge.parentGammaParameterRange[1]+
			[quotientEdgeNode1,quotientEdgeCentroid,quotientEdgeNode2]*geoEdge.parentGammaParameterRange[2]
	end
	
	meshData.edgesParentDomainsTypes = edgesParentDomainsTypes
	meshData.edgesParentDomainsIds = edgesParentDomainsIds
	meshData.edgesParentGammasIds = edgesParentGammasIds
	meshData.edgesParentGammasParameterRanges = edgesParentGammasParameterRanges
	
	return nothing
end

function meshDataElementsInitializeParents(meshData::MeshData,geoData::GeoData)
	
	parentGeoEdgesLoopsIds = map(subdomainsIntersection->subdomainsIntersection.parentGeoEdgesLoopId,geoData.subdomainsIntersections)
	
	geoEdgesLoopsParentDomainTypes = map(geoEdgesLoop->geoEdgesLoop.parentDomainType,geoData.geoEdgesLoops)
	geoEdgesLoopsParentDomainIds = map(geoEdgesLoop->geoEdgesLoop.parentDomainId,geoData.geoEdgesLoops)
	
	meshData.elementsParentDomainsTypes = geoEdgesLoopsParentDomainTypes[parentGeoEdgesLoopsIds[meshData.elementsParentSubdomainsIntersectionsIds]]
	meshData.elementsParentDomainsIds = geoEdgesLoopsParentDomainIds[parentGeoEdgesLoopsIds[meshData.elementsParentSubdomainsIntersectionsIds]]
	
	return nothing
end


