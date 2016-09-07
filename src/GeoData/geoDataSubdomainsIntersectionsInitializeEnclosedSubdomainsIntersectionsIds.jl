
using LightGraphs

function geoDataSubdomainsIntersectionsInitializeEnclosedSubdomainsIntersectionsIds(geoData::GeoData)
	
	graph = DiGraph(length(geoData.subdomainsIntersections))
	
	edgeWeights = zeros(length(geoData.subdomainsIntersections),length(geoData.subdomainsIntersections))
	for subdomainsIntersection in geoData.subdomainsIntersections
		for subdomainsIntersectionTmpId = 1:1:length(subdomainsIntersection.subsetRelationsToSubdomainsIntersections)
			if in(subdomainsIntersection.subsetRelationsToSubdomainsIntersections[subdomainsIntersectionTmpId],["touchingSuperset","strictSuperset"])
				add_edge!(graph,subdomainsIntersection.id,subdomainsIntersectionTmpId)
				edgeWeights[subdomainsIntersection.id,subdomainsIntersectionTmpId] = -1
			end
		end
	end
	
	if is_cyclic(graph)
		error("This should not happen. Graph contains cycle.")
	else
		dijkstraState = dijkstra_shortest_paths(graph,1,edgeWeights)
		
		for subdomainsIntersection in geoData.subdomainsIntersections[2:1:end]
			parentId = dijkstraState.parents[subdomainsIntersection.id]
			if parentId != 0
				push!(geoData.subdomainsIntersections[parentId].enclosedSubdomainsIntersectionsIds,subdomainsIntersection.id)
			end
		end
	end
	
	return nothing
end


