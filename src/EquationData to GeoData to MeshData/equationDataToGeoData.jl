
function equationDataToGeoData(equationData::EquationData)
	
	allDomains = vcat(equationData.domain,equationData.holes,equationData.subdomains)
	
	# Compute geoNodes and geoEdges data for each domain (without intersections):
	geoNodesPairsIds = map(i->Array{Int64,1}(0),1:1:length(allDomains))
	geoNodesCoords = map(i->Array{Array{Float64,1},1}(0),1:1:length(allDomains))
	geoEdgesParentGammasIds = map(i->Array{Int64,1}(0),1:1:length(allDomains))
	geoEdgesParentGammasParameterRanges = map(i->Array{Array{Float64,1},1}(0),1:1:length(allDomains))
	geoEdgesLoopsBoundingBoxes = Array{Array{Float64,1},1}(length(allDomains))
	for domainIndex = 1:1:length(allDomains)
		for gamma in allDomains[domainIndex].gammas
			for indexGeoNodeInPartition = 1:1:gamma.parametrizationNumberOfPartitions
				s1 = (indexGeoNodeInPartition-1)/(gamma.parametrizationNumberOfPartitions)
				s2 = indexGeoNodeInPartition/(gamma.parametrizationNumberOfPartitions)
				push!(geoNodesPairsIds[domainIndex],0)
				push!(geoNodesCoords[domainIndex],gamma.parametrization(s1))
				push!(geoEdgesParentGammasIds[domainIndex],gamma.id)
				push!(geoEdgesParentGammasParameterRanges[domainIndex],[s1,s2])
			end
		end
		geoEdgesLoopsBoundingBoxes[domainIndex] = computePolygonBoundingBox(geoNodesCoords[domainIndex])
	end
	
	# Compute and sort intersection points data for each domain:
	intersectionPointsDoOccur = falses(length(allDomains))
	intersectionPointsPairsIds = map(i->Array{Int64,1}(0),1:1:length(allDomains))
	intersectionPointsPredecessorsIndices = map(i->Array{Int64,1}(0),1:1:length(allDomains))
	intersectionPointsConvexCombinationParameterValues = map(i->Array{Float64,1}(0),1:1:length(allDomains))
	intersectionPointsCoords = map(i->Array{Array{Float64,1},1}(0),1:1:length(allDomains))
	
	pairsIdsStarter = 1
	for domain1Index = 1:1:length(allDomains)
		for domain2Index = (domain1Index+1):1:length(allDomains)
			polygonsDoIntersect,
			polygon1IntersectionPointsPredecessorsIds,polygon1IntersectionPointsParameterValues,polygon1IntersectionPointsCoords,
			polygon2IntersectionPointsPredecessorsIds,polygon2IntersectionPointsParameterValues,polygon2IntersectionPointsCoords = 
				computeIntersectionPointsBetweenTwoPolygons(
					geoNodesCoords[domain1Index],geoEdgesLoopsBoundingBoxes[domain1Index],
					geoNodesCoords[domain2Index],geoEdgesLoopsBoundingBoxes[domain2Index])
			
			if polygonsDoIntersect
				pairsIds = pairsIdsStarter - 1 + collect(1:1:length(polygon1IntersectionPointsPredecessorsIds))
				pairsIdsStarter += length(polygon1IntersectionPointsPredecessorsIds)
				
				intersectionPointsDoOccur[domain1Index] = true
				append!(intersectionPointsPairsIds[domain1Index],pairsIds)
				append!(intersectionPointsPredecessorsIndices[domain1Index],polygon1IntersectionPointsPredecessorsIds)
				append!(intersectionPointsConvexCombinationParameterValues[domain1Index],polygon1IntersectionPointsParameterValues)
				append!(intersectionPointsCoords[domain1Index],polygon1IntersectionPointsCoords)
				
				intersectionPointsDoOccur[domain2Index] = true
				append!(intersectionPointsPairsIds[domain2Index],pairsIds)
				append!(intersectionPointsPredecessorsIndices[domain2Index],polygon2IntersectionPointsPredecessorsIds)
				append!(intersectionPointsConvexCombinationParameterValues[domain2Index],polygon2IntersectionPointsParameterValues)
				append!(intersectionPointsCoords[domain2Index],polygon2IntersectionPointsCoords)
			end
		end
	end
	
	for domainIndex = 1:1:length(allDomains)
		if intersectionPointsDoOccur[domainIndex]
			permutation1 = sortperm(intersectionPointsConvexCombinationParameterValues[domainIndex])
			intersectionPointsPairsIds[domainIndex] = intersectionPointsPairsIds[domainIndex][permutation1]
			intersectionPointsPredecessorsIndices[domainIndex] = intersectionPointsPredecessorsIndices[domainIndex][permutation1]
			intersectionPointsConvexCombinationParameterValues[domainIndex] = intersectionPointsConvexCombinationParameterValues[domainIndex][permutation1]
			intersectionPointsCoords[domainIndex] = intersectionPointsCoords[domainIndex][permutation1]
			
			permutation2 = sortperm(intersectionPointsPredecessorsIndices[domainIndex])
			intersectionPointsPairsIds[domainIndex] = intersectionPointsPairsIds[domainIndex][permutation2]
			intersectionPointsPredecessorsIndices[domainIndex] = intersectionPointsPredecessorsIndices[domainIndex][permutation2]
			intersectionPointsConvexCombinationParameterValues[domainIndex] = intersectionPointsConvexCombinationParameterValues[domainIndex][permutation2]
			intersectionPointsCoords[domainIndex] = intersectionPointsCoords[domainIndex][permutation2]
		end
	end
	
	# Insert intersection points data into geoNodes and geoEdges data:
	
	# TODO: Set threshold properly.
	threshold = 0.000001
	
	for domainIndex = 1:1:length(allDomains)
		if intersectionPointsDoOccur[domainIndex]
			for i = length(intersectionPointsPredecessorsIndices[domainIndex]):-1:1
				
				predecessorIndex = intersectionPointsPredecessorsIndices[domainIndex][i]
				if abs(intersectionPointsConvexCombinationParameterValues[domainIndex][i]-0.0) <= threshold
					geoNodesPairsIds[domainIndex][predecessorIndex] = intersectionPointsPairsIds[domainIndex][i]
				else
					insert!(geoNodesPairsIds[domainIndex],predecessorIndex+1,intersectionPointsPairsIds[domainIndex][i])
					insert!(geoNodesCoords[domainIndex],predecessorIndex+1,intersectionPointsCoords[domainIndex][i])
					insert!(geoEdgesParentGammasIds[domainIndex],predecessorIndex+1,geoEdgesParentGammasIds[domainIndex][predecessorIndex])
					insert!(geoEdgesParentGammasParameterRanges[domainIndex],predecessorIndex+1,
						[evaluateConvexCombination(
						geoEdgesParentGammasParameterRanges[domainIndex][predecessorIndex][1],
						geoEdgesParentGammasParameterRanges[domainIndex][predecessorIndex][2],
						intersectionPointsConvexCombinationParameterValues[domainIndex][i]),
						-1.0])
				end
			end
		end
	end
	
	for domainIndex = 1:1:length(allDomains)
		if intersectionPointsDoOccur[domainIndex]
			for i = 1:1:length(geoEdgesParentGammasParameterRanges[domainIndex])
				if i <= length(geoEdgesParentGammasParameterRanges[domainIndex]) - 1
					if geoEdgesParentGammasIds[domainIndex][i+1] == geoEdgesParentGammasIds[domainIndex][i]
						geoEdgesParentGammasParameterRanges[domainIndex][i][2] = geoEdgesParentGammasParameterRanges[domainIndex][i+1][1]
					else
						geoEdgesParentGammasParameterRanges[domainIndex][i][2] = 1.0
					end
				else
					geoEdgesParentGammasParameterRanges[domainIndex][i][2] = 1.0
				end
			end
		end
	end
	
	# Build geoNodes, geoEdges, geoEdgesLoops and geoNodesIntersections from geoNodes and geoEdges data:
	geoNodesTotalNumber = sum(map(domainGeoNodesCoords->length(domainGeoNodesCoords),geoNodesCoords))
	geoEdgesTotalNumber = geoNodesTotalNumber
	geoEdgesLoopsTotalNumber = length(geoNodesCoords)
	geoNodesIntersectionsTotalNumber = maximum(vcat([0],map(intersectionPointsPairsIdsPerDomain->maximum(intersectionPointsPairsIdsPerDomain),intersectionPointsPairsIds[intersectionPointsDoOccur])))
	
	geoNodes = Array{GeoNode}(geoNodesTotalNumber)
	geoEdges = Array{GeoEdge}(geoEdgesTotalNumber)
	geoEdgesLoops = Array{GeoEdgesLoop}(geoEdgesLoopsTotalNumber)
	geoNodesIntersections = map(i->GeoNodesIntersection(i,Array{Int64,1}(0)),1:1:geoNodesIntersectionsTotalNumber)
	
	geoNodeId = 1
	geoEdgeId = 1
	geoEdgesLoopId = 1
	for domainIndex = 1:1:length(allDomains)
		
		firstInLoopGeoNodeId = geoNodeId
		firstInLoopGeoEdgeId = geoEdgeId
		
		for i = 1:1:length(geoNodesCoords[domainIndex])
			geoNodes[geoNodeId] = GeoNode(geoNodeId,geoNodesCoords[domainIndex][i],true)
			geoEdges[geoEdgeId] = GeoEdge(geoEdgeId,[geoNodeId,geoNodeId+1],
				allDomains[domainIndex].domainType,allDomains[domainIndex].id,
				geoEdgesParentGammasIds[domainIndex][i],geoEdgesParentGammasParameterRanges[domainIndex][i],
				true)
			
			geoNodesPairId = geoNodesPairsIds[domainIndex][i]
			if geoNodesPairId != 0
				push!(geoNodesIntersections[geoNodesPairId].geoNodesIds,geoNodeId)
			end
			
			geoNodeId += 1
			geoEdgeId += 1
		end
		
		# Overwrite last edge's geoNodesIds[2]:
		geoEdges[geoEdgeId-1].geoNodesIds[2] = firstInLoopGeoNodeId
		
		# Add current geoEdgeLoop:
		geoEdgesLoops[geoEdgesLoopId] = GeoEdgesLoop(geoEdgesLoopId,collect(firstInLoopGeoEdgeId:1:(geoEdgeId-1)),allDomains[domainIndex].domainType,allDomains[domainIndex].id,"",[],[])
		geoEdgesLoopId += 1
	end
	
	# Assemble geoData:
	geoData = GeoData(geoNodes,geoEdges,geoEdgesLoops,geoNodesIntersections,[],[])
	
	# Initialize additional fields:
	geoDataGeoEdgesLoopsInitializeOrientations(geoData)
	geoDataGeoEdgesLoopsInitializeBoundingBoxes(geoData)
	geoDataGeoEdgesLoopsInitializeSubsetRelationsToGeoEdgesLoops(geoData)
	
	geoDataGeoNodesAndGeoEdgesInitializeNonIsActives(geoData)
	
	geoDataInitializeSubdomainsIntersections(geoData)
	geoDataSubdomainsIntersectionsInitializeOrientations(geoData)
	geoDataSubdomainsIntersectionsInitializeBoundingBoxes(geoData)
	geoDataSubdomainsIntersectionsInitializeSubsetRelationsToSubdomainsIntersections(geoData)
	geoDataSubdomainsIntersectionsInitializeSubsetRelationsToGeoEdgesLoops(geoData)
	geoDataSubdomainsIntersectionsInitializeParentGeoEdgesLoopIds(geoData)
	geoDataSubdomainsIntersectionsInitializeEnclosedSubdomainsIntersectionsIds(geoData)
	
	geoDataInitializePeriodicGeoEdgesPairs(equationData,geoData)
	
	return geoData
end


