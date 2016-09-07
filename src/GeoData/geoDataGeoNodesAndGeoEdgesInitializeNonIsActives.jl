
function geoDataGeoNodesAndGeoEdgesInitializeNonIsActives(geoData::GeoData)
	
	for geoEdgesLoop in geoData.geoEdgesLoops
		if geoEdgesLoop.parentDomainType == "subdomain"
			
			geoEdgesLoopIsStrictSupersetOfDomain = (geoEdgesLoop.subsetRelationsToGeoEdgesLoops[1]=="strictSuperset")
			geoEdgesLoopAndDomainAreStrictDisjoint = (geoEdgesLoop.subsetRelationsToGeoEdgesLoops[1]=="strictDisjoint")
			geoEdgesLoopIsStrictSubsetOfAHole = any(map(geoEdgesLoopTmp->(geoEdgesLoopTmp.parentDomainType=="hole")&(geoEdgesLoop.subsetRelationsToGeoEdgesLoops[geoEdgesLoopTmp.id]=="strictSubset"),geoData.geoEdgesLoops))
			geoEdgesLoopInterferesWithDomainOrAHole = any(map(geoEdgesLoopTmp->in(geoEdgesLoopTmp.parentDomainType,["domain","hole"])&in(geoEdgesLoop.subsetRelationsToGeoEdgesLoops[geoEdgesLoopTmp.id],["touchingSubset","touchingSuperset","intersecting"]),geoData.geoEdgesLoops))
			
			if geoEdgesLoopIsStrictSupersetOfDomain|geoEdgesLoopAndDomainAreStrictDisjoint|geoEdgesLoopIsStrictSubsetOfAHole
				
				for geoEdge in geoData.geoEdges[geoEdgesLoop.geoEdgesIds]
					geoEdge.isActive = false
					geoData.geoNodes[geoEdge.geoNodesIds[1]].isActive = false
					geoData.geoNodes[geoEdge.geoNodesIds[2]].isActive = false
				end
				
			elseif geoEdgesLoopInterferesWithDomainOrAHole
				
				for geoEdge in geoData.geoEdges[geoEdgesLoop.geoEdgesIds]
					
					geoEdgeIsActive = true
					
					subsetRelationGeoEdge,closeGeoEdgeId = geoDataComputeSubsetRelationGeoEdgeGeoEdgesIds(geoData,geoEdge,geoData.geoEdgesLoops[1].boundingBox,geoData.geoEdgesLoops[1].geoEdgesIds)
					if in(subsetRelationGeoEdge,["touchingInside","strictInside"])
						
						for geoEdgesLoopTmp in geoData.geoEdgesLoops
							if geoEdgesLoopTmp.parentDomainType == "hole"
								subsetRelationGeoEdge2,closeGeoEdgeId2 = geoDataComputeSubsetRelationGeoEdgeGeoEdgesIds(geoData,geoEdge,geoEdgesLoopTmp.boundingBox,geoEdgesLoopTmp.geoEdgesIds)
								
								if in(subsetRelationGeoEdge2,["touchingInside","strictInside"])
									geoEdgeIsActive = false
									break
								end
							end
						end
						
					elseif in(subsetRelationGeoEdge,["oneOfGeoEdges","closeToOneOfGeoEdges","touchingOutside","strictOutside"])
						geoEdgeIsActive = false
					end
					
					if !geoEdgeIsActive
						geoEdge.isActive = false
						geoData.geoNodes[geoEdge.geoNodesIds[1]].isActive = false
						geoData.geoNodes[geoEdge.geoNodesIds[2]].isActive = false
					end
					
				end
				
			end
			
		end
	end
	
	# Assure that isActives of geoNodes of subdomains at intersections are set inactive:
	for geoNodesIntersection in geoData.geoNodesIntersections
		inactiveGeoNodeId = maximum(geoNodesIntersection.geoNodesIds)
		geoData.geoNodes[inactiveGeoNodeId].isActive = false
	end
	
	return nothing
end


