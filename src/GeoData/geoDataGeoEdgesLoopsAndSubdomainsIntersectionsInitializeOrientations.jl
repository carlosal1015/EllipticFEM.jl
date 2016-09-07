
function geoDataGeoEdgesLoopsInitializeOrientations(geoData::GeoData)
	
	for geoEdgesLoop in geoData.geoEdgesLoops
		
		geoNodesIds = map(geoEdgeId->geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)],geoEdgesLoop.geoEdgesIds)
		geoNodesCoords = map(geoNode->geoNode.coords,geoData.geoNodes[geoNodesIds])
		
		geoEdgesLoop.orientation = computePolygonOrientation(geoNodesCoords)
	end
	
	return nothing
end

function geoDataSubdomainsIntersectionsInitializeOrientations(geoData::GeoData)
	
	for subdomainsIntersection in geoData.subdomainsIntersections
		
		geoNodesIds = map(geoEdgeId->geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)],subdomainsIntersection.geoEdgesIds)
		geoNodesCoords = map(geoNode->geoNode.coords,geoData.geoNodes[geoNodesIds])
		
		subdomainsIntersection.orientation = computePolygonOrientation(geoNodesCoords)
	end
	
	return nothing
end

