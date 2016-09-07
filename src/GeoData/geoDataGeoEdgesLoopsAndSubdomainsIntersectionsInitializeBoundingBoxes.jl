
function geoDataGeoEdgesLoopsInitializeBoundingBoxes(geoData::GeoData)
	
	for geoEdgesLoop in geoData.geoEdgesLoops
		
		geoNodesIds = map(geoEdge->geoEdge.geoNodesIds[1],geoData.geoEdges[geoEdgesLoop.geoEdgesIds])
		polygonVerticesCoords = map(geoNode->geoNode.coords,geoData.geoNodes[geoNodesIds])
		
		geoEdgesLoop.boundingBox = computePolygonBoundingBox(polygonVerticesCoords)
	end
	
	return nothing
end

function geoDataSubdomainsIntersectionsInitializeBoundingBoxes(geoData::GeoData)
	
	for subdomainsIntersection in geoData.subdomainsIntersections
		
		geoNodesIds = map(geoEdgeId->geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)],subdomainsIntersection.geoEdgesIds)
		polygonVerticesCoords = map(geoNode->geoNode.coords,geoData.geoNodes[geoNodesIds])
		
		subdomainsIntersection.boundingBox = computePolygonBoundingBox(polygonVerticesCoords)
	end
	
	return nothing
end


