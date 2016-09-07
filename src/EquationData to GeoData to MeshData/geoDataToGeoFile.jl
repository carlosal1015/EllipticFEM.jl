

function geoDataToGeoFile(geoData::GeoData,path::ASCIIString)
	
	characteristicLengths = 1
	zCoord = 0
	
	# Assemble contentGeoNodes:
	contentGeoNodes = ""
	for geoNode in geoData.geoNodes
		if geoNode.isActive
			contentGeoNodes *= "Point("*string(geoNode.id)*") = {"*
				string(geoNode.coords[1])*", "*
				string(geoNode.coords[2])*", "*
				string(zCoord)*", "*
				string(characteristicLengths)*"};"*"\n"
		end
	end
	
	# Assemble contentGeoEdges:
	geoNodesIdsMapping = collect(1:1:length(geoData.geoNodes))
	for geoNodesIntersection in geoData.geoNodesIntersections
		masterGeoNodeId = minimum(geoNodesIntersection.geoNodesIds)
		slaveGeoNodeId = maximum(geoNodesIntersection.geoNodesIds)
		geoNodesIdsMapping[slaveGeoNodeId] = masterGeoNodeId
	end
	
	contentGeoEdges = ""
	for geoEdge in geoData.geoEdges
		if geoEdge.isActive
			contentGeoEdges *= "Line("*string(geoEdge.id)*") = {"*
				string(geoNodesIdsMapping[geoEdge.geoNodesIds[1]])*", "*
				string(geoNodesIdsMapping[geoEdge.geoNodesIds[2]])*"};"*"\n"
		end
	end
	
	# Assemble contentSubdomainsIntersections:
	contentSubdomainsIntersectionsLineLoops = ""
	contentSubdomainsIntersectionsPlaneSurfaces = ""
	for subdomainsIntersection in geoData.subdomainsIntersections
		if subdomainsIntersection.isActive
			
			geoEdgesIdsToWrite = copy(subdomainsIntersection.geoEdgesIds)
			for enclosedSubdomainsIntersection in geoData.subdomainsIntersections[subdomainsIntersection.enclosedSubdomainsIntersectionsIds]
				append!(geoEdgesIdsToWrite,enclosedSubdomainsIntersection.geoEdgesIds)
			end
			
			contentSubdomainsIntersectionsLineLoops *= "Line Loop("*string(subdomainsIntersection.id)*") = {"*
				replace(string(geoEdgesIdsToWrite),",",", ")[2:1:(end-1)]*"};"*"\n"
			contentSubdomainsIntersectionsPlaneSurfaces *= "Plane Surface("*string(subdomainsIntersection.id)*") = {"*
				string(subdomainsIntersection.id)*"};"*"\n"
		end
	end
	contentSubdomainsIntersections = contentSubdomainsIntersectionsLineLoops*contentSubdomainsIntersectionsPlaneSurfaces
	
	# Assemble contentPeriodicGeoEdgesPairs:
	contentPeriodicGeoEdgesPairs = ""
	for periodicGeoEdgesPair in geoData.periodicGeoEdgesPairs
		if periodicGeoEdgesPair.isActive
			contentPeriodicGeoEdgesPairs *= "Periodic Line{"*string(periodicGeoEdgesPair.geoEdgesIds[2])*"} = {"*
				string(periodicGeoEdgesPair.geoEdgesIds[1])*"};"*"\n"
		end
	end
	
	# Write contents to file:
	filehandle = open(path,"w")
	print(filehandle,contentGeoNodes*contentGeoEdges*contentSubdomainsIntersections*contentPeriodicGeoEdgesPairs)
	close(filehandle)
	
	return nothing
end



