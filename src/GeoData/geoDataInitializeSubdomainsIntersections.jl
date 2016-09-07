
using LightGraphs

function geoDataInitializeSubdomainsIntersections(geoData::GeoData)
	
	connectedGeoEdgesLoopsBundles = geoDataSubdomainsIntersectionsComputeConnectedGeoEdgesLoopsBundles(geoData)
	
	subdomainsIntersections = Array{SubdomainsIntersection}(0)
	
	subdomainsIntersectionId = 1
	for connectedGeoEdgesLoopsBundle in connectedGeoEdgesLoopsBundles
		
		if length(connectedGeoEdgesLoopsBundle)==1
			
			for geoEdgesLoop in geoData.geoEdgesLoops[connectedGeoEdgesLoopsBundle]
				
				subdomainsIntersection = SubdomainsIntersection(
					subdomainsIntersectionId,
					geoEdgesLoop.id,
					1,
					(geoEdgesLoop.orientation=="positive" ? copy(geoEdgesLoop.geoEdgesIds) : -(geoEdgesLoop.geoEdgesIds[end:(-1):1])),
					"",
					[],
					[],
					[],
					[],
					-1,
					in(geoEdgesLoop.parentDomainType,["domain","subdomain"]))
				
				push!(subdomainsIntersections,subdomainsIntersection)
				subdomainsIntersectionId += 1
			end
			
		elseif length(connectedGeoEdgesLoopsBundle)>=2
			
			remainingGeoEdgesBools = Array{Array{Bool,1},1}(length(geoData.geoEdgesLoops))
			for geoEdgesLoopTmp in geoData.geoEdgesLoops[connectedGeoEdgesLoopsBundle]
				remainingGeoEdgesBools[geoEdgesLoopTmp.id] = map(geoEdge->geoEdge.isActive,geoData.geoEdges[geoEdgesLoopTmp.geoEdgesIds])
			end
			
			outsideStarterGeoEdgeHasBeenSet = false
			outsideStarterGeoEdgesLoopId = -1
			outsideStarterGeoEdgesLoopGeoEdgesIdsIndex = -1
			
			for geoEdgesLoop in geoData.geoEdgesLoops[connectedGeoEdgesLoopsBundle]
				
				counterMax = sum(remainingGeoEdgesBools[geoEdgesLoop.id])
				for counter = 1:1:(counterMax+1)
					
					starterGeoEdgesLoopGeoEdgesIdsIndex = findfirst(remainingGeoEdgesBools[geoEdgesLoop.id])
					if starterGeoEdgesLoopGeoEdgesIdsIndex != 0
						
						subdomainsIntersection = SubdomainsIntersection(
							subdomainsIntersectionId,
							geoEdgesLoop.id,
							starterGeoEdgesLoopGeoEdgesIdsIndex,
							[],
							"",
							[],
							[],
							[],
							[],
							-1,
							in(geoEdgesLoop.parentDomainType,["domain","subdomain"]))
							
						geoDataSubdomainsIntersectionsComputeGeoEdgesIds(geoData,connectedGeoEdgesLoopsBundle,subdomainsIntersection,"leftmost")
						
						push!(subdomainsIntersections,subdomainsIntersection)
						subdomainsIntersectionId += 1
						
						for geoEdgesLoopTmp in geoData.geoEdgesLoops[connectedGeoEdgesLoopsBundle]
							indices = findin((geoEdgesLoopTmp.orientation=="positive" ? 1 : -1)*geoEdgesLoopTmp.geoEdgesIds,subdomainsIntersection.geoEdgesIds)
							remainingGeoEdgesBools[geoEdgesLoopTmp.id][indices] = false
						end
						
						if !outsideStarterGeoEdgeHasBeenSet
							subsetRelationsStarterGeoEdge = repmat([""],length(geoData.geoEdgesLoops))
							for geoEdgesLoopTmp in geoData.geoEdgesLoops[connectedGeoEdgesLoopsBundle]
								starterGeoEdgeId = geoEdgesLoop.geoEdgesIds[starterGeoEdgesLoopGeoEdgesIdsIndex]
								subsetRelationStarterGeoEdge,closeGeoEdgeIdStarterGeoEdge =
									geoDataComputeSubsetRelationGeoEdgeGeoEdgesIds(geoData,geoData.geoEdges[starterGeoEdgeId],geoEdgesLoopTmp.boundingBox,geoEdgesLoopTmp.geoEdgesIds)
								
								subsetRelationsStarterGeoEdge[geoEdgesLoopTmp.id] = subsetRelationStarterGeoEdge
							end
							
							if (in("touchingOutside",subsetRelationsStarterGeoEdge)|in("strictOutside",subsetRelationsStarterGeoEdge)) & 
								!(in("touchingInside",subsetRelationsStarterGeoEdge)|in("strictInside",subsetRelationsStarterGeoEdge))
								
								outsideStarterGeoEdgesLoopId = geoEdgesLoop.id
								outsideStarterGeoEdgesLoopGeoEdgesIdsIndex = starterGeoEdgesLoopGeoEdgesIdsIndex
								
								outsideStarterGeoEdgeHasBeenSet = true
							end
						end
						
					else
						break
					end
					
					if counter == (counterMax+1)
						warn("There are more subdomainsIntersections than active geoEdges in geoEdgesLoop "*string(geoEdgesLoop.id)*".")
					end
					
				end
			end
			
			if outsideStarterGeoEdgeHasBeenSet
				if !in(1,connectedGeoEdgesLoopsBundle)
					
					subdomainsIntersection = SubdomainsIntersection(
						subdomainsIntersectionId,
						outsideStarterGeoEdgesLoopId,
						outsideStarterGeoEdgesLoopGeoEdgesIdsIndex,
						[],
						"",
						[],
						[],
						[],
						[],
						-1,
						false)
					
					geoDataSubdomainsIntersectionsComputeGeoEdgesIds(geoData,connectedGeoEdgesLoopsBundle,subdomainsIntersection,"rightmost")
					
					push!(subdomainsIntersections,subdomainsIntersection)
					subdomainsIntersectionId += 1
				end
			else
				warn("connectedGeoEdgesLoopsBundle "*string(connectedGeoEdgesLoopsBundle)*" has no starter geoEdge on outside boundary.")
			end
			
		end
	end
	
	geoData.subdomainsIntersections = subdomainsIntersections
	
	return nothing
end

function geoDataSubdomainsIntersectionsComputeConnectedGeoEdgesLoopsBundles(geoData::GeoData)
	
	graph = Graph(length(geoData.geoEdgesLoops))
	
	for geoEdgesLoop in geoData.geoEdgesLoops
		for geoEdgesLoopTmpId = 1:1:length(geoEdgesLoop.subsetRelationsToGeoEdgesLoops)
			if in(geoEdgesLoop.subsetRelationsToGeoEdgesLoops[geoEdgesLoopTmpId],["touchingSubset","touchingSuperset","intersecting","touchingDisjoint"])
				add_edge!(graph,geoEdgesLoop.id,geoEdgesLoopTmpId)
			end
		end
	end
	
	connectedGeoEdgesLoopsBundles = connected_components(graph)
	
	return connectedGeoEdgesLoopsBundles
end

function geoDataSubdomainsIntersectionsComputeGeoEdgesIds(geoData::GeoData,connectedGeoEdgesLoopsBundle::Array{Int64},subdomainsIntersection::SubdomainsIntersection,leftmostRightmostChooser::ASCIIString)
	
	# Compute lookup tables from geoData.geoNodesIntersections:
	intersectionsGeoNodesIds = Array{Int64,1}(0)
	nextGeoEdgesLoopsIdsLookupTable = Array{Int64,1}(0)
	nextGeoEdgesLoopsGeoEdgesIdsIndicesLookupTable = Array{Array{Int64,1},1}(0)
	nextDirectionsLookupTable = Array{Array{ASCIIString,1},1}(0)
	for geoNodesIntersection in geoData.geoNodesIntersections
		for k = 1:1:2
			currentGeoNodeId = geoNodesIntersection.geoNodesIds[k]
			nextGeoNodeId = geoNodesIntersection.geoNodesIds[(k==1 ? 2 : 1)]
			nextGeoEdgesLoopId = findfirst(map(geoEdgesLoop->in(nextGeoNodeId,geoEdgesLoop.geoEdgesIds),geoData.geoEdgesLoops))
			
			push!(intersectionsGeoNodesIds,currentGeoNodeId)
			push!(nextGeoEdgesLoopsIdsLookupTable,nextGeoEdgesLoopId)
			index = findfirst(geoData.geoEdgesLoops[nextGeoEdgesLoopId].geoEdgesIds,nextGeoNodeId)
			push!(nextGeoEdgesLoopsGeoEdgesIdsIndicesLookupTable,[moduloIndex(index-1,length(geoData.geoEdgesLoops[nextGeoEdgesLoopId].geoEdgesIds)),index])
			push!(nextDirectionsLookupTable,["descendingIds","ascendingIds"])
		end
	end
	permutation = sortperm(intersectionsGeoNodesIds)
	intersectionsGeoNodesIds = intersectionsGeoNodesIds[permutation]
	nextGeoEdgesLoopsIdsLookupTable = nextGeoEdgesLoopsIdsLookupTable[permutation]
	nextGeoEdgesLoopsGeoEdgesIdsIndicesLookupTable = nextGeoEdgesLoopsGeoEdgesIdsIndicesLookupTable[permutation]
	nextDirectionsLookupTable = nextDirectionsLookupTable[permutation]
	
	# Assemble subdomainsIntersection.geoEdgesIds:
	geoEdgesLoopId = subdomainsIntersection.starterGeoEdgesLoopId
	geoEdgesLoopGeoEdgesIdsIndex = subdomainsIntersection.starterGeoEdgesLoopGeoEdgesIdsIndex
	direction = (geoData.geoEdgesLoops[geoEdgesLoopId].orientation=="positive" ? "ascendingIds" : "descendingIds")
	
	counterMax = length(geoData.geoEdges)
	for counter = 1:1:(counterMax+2)
		
		if counter==(counterMax+2)
			warn("There are more subdomainsIntersection.geoEdgesIds than the total number of geoData.geoEdges.")
		end
		
		geoEdgeId = geoData.geoEdgesLoops[geoEdgesLoopId].geoEdgesIds[geoEdgesLoopGeoEdgesIdsIndex]
		
		# Break the loop or set next geoEdgesLoopGeoEdgesIdsIndex:
		if (length(subdomainsIntersection.geoEdgesIds)>=1)&
			(geoEdgeId==geoData.geoEdgesLoops[subdomainsIntersection.starterGeoEdgesLoopId].geoEdgesIds[subdomainsIntersection.starterGeoEdgesLoopGeoEdgesIdsIndex])
			break
		else
			# Add the current geoEdgeId to the final collection:
			push!(subdomainsIntersection.geoEdgesIds,(direction=="ascendingIds" ? geoEdgeId : -geoEdgeId))
			
			# Set the candidates for the next geoEdge:
			nextGeoEdgesLoopsIds = [geoEdgesLoopId]
			nextGeoEdgesLoopsGeoEdgesIdsIndices = [moduloIndex(geoEdgesLoopGeoEdgesIdsIndex + (direction=="ascendingIds" ? 1 : -1),length(geoData.geoEdgesLoops[geoEdgesLoopId].geoEdgesIds))]
			nextDirections = [direction]
			
			intersectionGeoNodeId = geoData.geoEdges[geoEdgeId].geoNodesIds[(direction=="ascendingIds" ? 2 : 1)]
			lookupTableIndex = findfirst(intersectionsGeoNodesIds,intersectionGeoNodeId)
			if lookupTableIndex != 0
				append!(nextGeoEdgesLoopsIds,repmat([nextGeoEdgesLoopsIdsLookupTable[lookupTableIndex]],2))
				append!(nextGeoEdgesLoopsGeoEdgesIdsIndices,nextGeoEdgesLoopsGeoEdgesIdsIndicesLookupTable[lookupTableIndex])
				append!(nextDirections,nextDirectionsLookupTable[lookupTableIndex])
			end
			
			# Choose the next geoEdge from candidates:
			if length(nextGeoEdgesLoopsIds) == 1
				nextArraysIndex = 1
			elseif length(nextGeoEdgesLoopsIds) >= 2
				
				orientedAngles = Array{Float64}(nextGeoEdgesLoopsIds)
				for i = 1:1:length(nextGeoEdgesLoopsIds)
					geoEdgeCoords = -(geoData.geoNodes[geoData.geoEdges[geoEdgeId].geoNodesIds[(direction=="ascendingIds" ? 2 : 1)]].coords - 
						geoData.geoNodes[geoData.geoEdges[geoEdgeId].geoNodesIds[(direction=="ascendingIds" ? 1 : 2)]].coords)
					
					nextGeoEdgeId = geoData.geoEdgesLoops[nextGeoEdgesLoopsIds[i]].geoEdgesIds[nextGeoEdgesLoopsGeoEdgesIdsIndices[i]]
					nextGeoEdgeCoords = geoData.geoNodes[geoData.geoEdges[nextGeoEdgeId].geoNodesIds[(nextDirections[i]=="ascendingIds" ? 2 : 1)]].coords - 
						geoData.geoNodes[geoData.geoEdges[nextGeoEdgeId].geoNodesIds[(nextDirections[i]=="ascendingIds" ? 1 : 2)]].coords
					
					orientedAngles[i] = orientedAngle(geoEdgeCoords,nextGeoEdgeCoords)
				end
				
				if leftmostRightmostChooser == "leftmost"
					nextArraysIndex = indmax(orientedAngles)	# Chooses leftmost geoEdge.
				elseif leftmostRightmostChooser == "rightmost"
					nextArraysIndex = indmin(orientedAngles)	# Chooses rightmost geoEdge.
				end
				
			end
			
			geoEdgesLoopId = nextGeoEdgesLoopsIds[nextArraysIndex]
			geoEdgesLoopGeoEdgesIdsIndex = nextGeoEdgesLoopsGeoEdgesIdsIndices[nextArraysIndex]
			direction = nextDirections[nextArraysIndex]
			
		end
		
	end
	
	return nothing
end


