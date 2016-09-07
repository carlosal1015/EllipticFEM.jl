

function geoDataComputeSubsetRelationGeoNodeGeoEdgesIds(geoData::GeoData,geoNode::GeoNode,geoEdgesIdsBoundingBox::Array{Float64},geoEdgesIds::Array{Int64})
	
	closeGeoNodeId = -1
	if !pointIsInBoundingBox(geoNode.coords,geoEdgesIdsBoundingBox)
		subsetRelationGeoNode = "strictOutside"
	else
		
		# TODO: Set threshold properly:
		threshold = 0.000001
		
		isOneOfGeoNodes = false
		isCloseToOneOfGeoNodes = false
		for geoEdgeId = geoEdgesIds
			
			geoEdgeGeoNode1 = geoData.geoNodes[geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)]]
			
			if geoNode.id == geoEdgeGeoNode1.id
				isOneOfGeoNodes = true
				closeGeoNodeId = geoEdgeGeoNode1.id
				break
			elseif norm(geoNode.coords-geoEdgeGeoNode1.coords) <= threshold
				isCloseToOneOfGeoNodes = true
				closeGeoNodeId = geoEdgeGeoNode1.id
				break
			end
		end
		
		if isOneOfGeoNodes
			subsetRelationGeoNode = "oneOfGeoNodes"
		elseif isCloseToOneOfGeoNodes
			subsetRelationGeoNode = "closeToOneOfGeoNodes"
		else
			
			# Perform Jordan point in polygon test:
			geoNodesIds = map(geoEdgeId->geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)],geoEdgesIds)
			subsetRelationGeoNode = performJordanPointInPolygonTest(geoNode.coords,map(geoNode->geoNode.coords,geoData.geoNodes[geoNodesIds]))
		end
		
	end
	
	return subsetRelationGeoNode,closeGeoNodeId
	# "oneOfGeoNodes", "closeToOneOfGeoNodes", "strictInside", "strictOutside"
end

function geoDataComputeSubsetRelationGeoEdgeGeoEdgesIds(geoData::GeoData,geoEdge::GeoEdge,geoEdgesIdsBoundingBox::Array{Float64},geoEdgesIds::Array{Int64})
	
	closeGeoEdgeId = -1
	if in(geoEdge.id,abs(geoEdgesIds))
		
		subsetRelationGeoEdge = "oneOfGeoEdges"
		closeGeoEdgeId = geoEdge.id
		
	else
		
		# Compute geoDataComputeSubsetRelationGeoNodeGeoEdgesIds for both geoNodes of geoEdge:
		subsetRelationGeoEdgeGeoNode1,closeGeoNodeIdGeoEdgeGeoNode1 = 
			geoDataComputeSubsetRelationGeoNodeGeoEdgesIds(geoData,geoData.geoNodes[geoEdge.geoNodesIds[1]],geoEdgesIdsBoundingBox,geoEdgesIds)
		subsetRelationGeoEdgeGeoNode2,closeGeoNodeIdGeoEdgeGeoNode2 = 
			geoDataComputeSubsetRelationGeoNodeGeoEdgesIds(geoData,geoData.geoNodes[geoEdge.geoNodesIds[2]],geoEdgesIdsBoundingBox,geoEdgesIds)
		
		# Deduce subsetRelationGeoEdge from subsetRelationGeoEdgeGeoNode1 and subsetRelationGeoEdgeGeoNode2:
		if in(subsetRelationGeoEdgeGeoNode1,["oneOfGeoNodes","closeToOneOfGeoNodes"])
			if in(subsetRelationGeoEdgeGeoNode2,["oneOfGeoNodes","closeToOneOfGeoNodes"])
				
				isCloseToOneOfGeoEdges = false
				for geoEdgeTmp in geoData.geoEdges[abs(geoEdgesIds)]
					if sort([closeGeoNodeIdGeoEdgeGeoNode1,closeGeoNodeIdGeoEdgeGeoNode2])==sort(geoEdgeTmp.geoNodesIds)
						
						isCloseToOneOfGeoEdges = true
						closeGeoEdgeId = geoEdgeTmp.id
						break
					end
				end
				
				if isCloseToOneOfGeoEdges
					subsetRelationGeoEdge = "closeToOneOfGeoEdges"
				else
					geoEdgeGeoNodeMidpoint = GeoNode(-1,(1/2)*(geoData.geoNodes[geoEdge.geoNodesIds[1]].coords + geoData.geoNodes[geoEdge.geoNodesIds[2]].coords),false)
					subsetRelationGeoEdgeGeoNodeMidpoint,closeGeoNodeIdGeoEdgeGeoNodeMidpoint = 
						geoDataComputeSubsetRelationGeoNodeGeoEdgesIds(geoData,geoEdgeGeoNodeMidpoint,geoEdgesIdsBoundingBox,geoEdgesIds)
					if in(subsetRelationGeoEdgeGeoNodeMidpoint,["oneOfGeoNodes","closeToOneOfGeoNodes"])
						subsetRelationGeoEdge = ""
						warn("This should not happen. geoEdge.id = "*string(geoEdge.id)*", geoEdgesIds = "*string(geoEdgesIds)*".")
					elseif subsetRelationGeoEdgeGeoNodeMidpoint == "strictInside"
						subsetRelationGeoEdge = "touchingInside"
					elseif subsetRelationGeoEdgeGeoNodeMidpoint == "strictOutside"
						subsetRelationGeoEdge = "touchingOutside"
					end
				end
				
			elseif subsetRelationGeoEdgeGeoNode2 == "strictInside"
				subsetRelationGeoEdge = "touchingInside"
			elseif subsetRelationGeoEdgeGeoNode2 == "strictOutside"
				subsetRelationGeoEdge = "touchingOutside"
			end
		elseif subsetRelationGeoEdgeGeoNode1 == "strictInside"
			if in(subsetRelationGeoEdgeGeoNode2,["oneOfGeoNodes","closeToOneOfGeoNodes"])
				subsetRelationGeoEdge = "touchingInside"
			elseif subsetRelationGeoEdgeGeoNode2 == "strictInside"
				subsetRelationGeoEdge = "strictInside"
			elseif subsetRelationGeoEdgeGeoNode2 == "strictOutside"
				subsetRelationGeoEdge = ""
				warn("This should not happen. geoEdge.id = "*string(geoEdge.id)*", geoEdgesIds = "*string(geoEdgesIds)*".")
			end
		elseif subsetRelationGeoEdgeGeoNode1 == "strictOutside"
			if in(subsetRelationGeoEdgeGeoNode2,["oneOfGeoNodes","closeToOneOfGeoNodes"])
				subsetRelationGeoEdge = "touchingOutside"
			elseif subsetRelationGeoEdgeGeoNode2 == "strictInside"
				subsetRelationGeoEdge = ""
				warn("This should not happen. geoEdge.id = "*string(geoEdge.id)*", geoEdgesIds = "*string(geoEdgesIds)*".")
			elseif subsetRelationGeoEdgeGeoNode2 == "strictOutside"
				subsetRelationGeoEdge = "strictOutside"
			end
		end
		
	end
	
	return subsetRelationGeoEdge,closeGeoEdgeId
	# "oneOfGeoEdges", "closeToOneOfGeoEdges", "touchingInside", "strictInside", "touchingOutside", "strictOutside"
end

function geoDataComputSubsetRelationGeoEdgesIdsGeoEdgesIds(geoData::GeoData,
	geoEdgesIds1Type::ASCIIString,geoEdgesIds1Id::Int64,geoEdgesIds1BoundingBox::Array{Float64},geoEdgesIds1::Array{Int64},
	geoEdgesIds2Type::ASCIIString,geoEdgesIds2Id::Int64,geoEdgesIds2BoundingBox::Array{Float64},geoEdgesIds2::Array{Int64})
	
	if ((geoEdgesIds1Type==geoEdgesIds2Type)&(geoEdgesIds1Id==geoEdgesIds2Id)) | 
		(geoEdgesIds1==geoEdgesIds2)
		subsetRelationGeoEdgesIds1 = "identical"
		subsetRelationGeoEdgesIds2 = "identical"
	else
		
		if !boundingBoxesDoIntersect(geoEdgesIds1BoundingBox,geoEdgesIds2BoundingBox)
			subsetRelationGeoEdgesIds1 = "strictDisjoint"
			subsetRelationGeoEdgesIds2 = "strictDisjoint"
		else
			
			# Compute subsetRelationGeoEdge for all geoEdges of geoEdgesIds1:
			subsetRelationsGeoEdges = Array{ASCIIString}(length(geoEdgesIds1))
			for index = 1:1:length(geoEdgesIds1)
				subsetRelationGeoEdge,closeGeoEdgeId = geoDataComputeSubsetRelationGeoEdgeGeoEdgesIds(
					geoData,geoData.geoEdges[abs(geoEdgesIds1[index])],geoEdgesIds2BoundingBox,geoEdgesIds2)
				
				subsetRelationsGeoEdges[index] = subsetRelationGeoEdge
			end
			
			# Deduce subsetRelationGeoEdgesIds1 and subsetRelationGeoEdgesIds2 from subsetRelationsGeoEdges:
			if all(subsetRelationsGeoEdges.=="strictInside")
				
				subsetRelationGeoEdgesIds1 = "strictSubset"
				subsetRelationGeoEdgesIds2 = "strictSuperset"
				
			elseif all(subsetRelationsGeoEdges.=="strictOutside")
				
				geoNodesIds2 = map(geoEdgeId->geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)],geoEdgesIds2)
				somePointInsideGeoEdgesIds2Coords = computeSomePointInPolygon(geoEdgesIds2BoundingBox,map(geoNode->geoNode.coords,geoData.geoNodes[geoNodesIds2]))
				
				subsetRelationSomePointInsideGeoEdgesIds2,closeGeoNodeIdSomePointInsideGeoEdgesIds2 = 
					geoDataComputeSubsetRelationGeoNodeGeoEdgesIds(geoData,GeoNode(-1,somePointInsideGeoEdgesIds2Coords,false),geoEdgesIds1BoundingBox,geoEdgesIds1)
				if in(subsetRelationSomePointInsideGeoEdgesIds2,["oneOfGeoNodes","closeToOneOfGeoNodes"])
					subsetRelationGeoEdgesIds1 = ""
					subsetRelationGeoEdgesIds2 = ""
					warn("This should not happen.")
				elseif subsetRelationSomePointInsideGeoEdgesIds2 == "strictInside"
					subsetRelationGeoEdgesIds1 = "strictSuperset"
					subsetRelationGeoEdgesIds2 = "strictSubset"
				elseif subsetRelationSomePointInsideGeoEdgesIds2 == "strictOutside"
					subsetRelationGeoEdgesIds1 = "strictDisjoint"
					subsetRelationGeoEdgesIds2 = "strictDisjoint"
				end
				
			elseif all(map(subsetRelation->in(subsetRelation,["oneOfGeoEdges","closeToOneOfGeoEdges"]),subsetRelationsGeoEdges))
				
				geoNodesIds1 = map(geoEdgeId->geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)],geoEdgesIds1)
				somePointInsideGeoEdgesIds1Coords = computeSomePointInPolygon(geoEdgesIds1BoundingBox,map(geoNode->geoNode.coords,geoData.geoNodes[geoNodesIds1]))
				
				subsetRelationSomePointInsideGeoEdgesIds1,closeGeoNodeIdSomePointInsideGeoEdgesIds1 = 
					geoDataComputeSubsetRelationGeoNodeGeoEdgesIds(geoData,GeoNode(-1,somePointInsideGeoEdgesIds1Coords,false),geoEdgesIds2BoundingBox,geoEdgesIds2)
				if in(subsetRelationSomePointInsideGeoEdgesIds1,["oneOfGeoNodes","closeToOneOfGeoNodes"])
					subsetRelationGeoEdgesIds1 = ""
					subsetRelationGeoEdgesIds2 = ""
					warn("This should not happen.")
				elseif subsetRelationSomePointInsideGeoEdgesIds1 == "strictInside"
					subsetRelationGeoEdgesIds1 = "touchingSubset"
					subsetRelationGeoEdgesIds2 = "touchingSuperset"
				elseif subsetRelationSomePointInsideGeoEdgesIds1 == "strictOutside"
					subsetRelationGeoEdgesIds1 = "touchingDisjoint"
					subsetRelationGeoEdgesIds2 = "touchingDisjoint"
				end
				
			else
				
				if in("touchingInside",subsetRelationsGeoEdges)
					if in("touchingOutside",subsetRelationsGeoEdges)
						subsetRelationGeoEdgesIds1 = "intersecting"
						subsetRelationGeoEdgesIds2 = "intersecting"
					else
						subsetRelationGeoEdgesIds1 = "touchingSubset"
						subsetRelationGeoEdgesIds2 = "touchingSuperset"
					end
				else
					geoNodesIds2 = map(geoEdgeId->geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)],geoEdgesIds2)
					somePointInsideGeoEdgesIds2Coords = computeSomePointInPolygon(geoEdgesIds2BoundingBox,map(geoNode->geoNode.coords,geoData.geoNodes[geoNodesIds2]))
					
					subsetRelationSomePointInsideGeoEdgesIds2,closeGeoNodeIdSomePointInsideGeoEdgesIds2 = 
						geoDataComputeSubsetRelationGeoNodeGeoEdgesIds(geoData,GeoNode(-1,somePointInsideGeoEdgesIds2Coords,false),geoEdgesIds1BoundingBox,geoEdgesIds1)
					if in(subsetRelationSomePointInsideGeoEdgesIds2,["oneOfGeoNodes","closeToOneOfGeoNodes"])
						subsetRelationGeoEdgesIds1 = ""
						subsetRelationGeoEdgesIds2 = ""
						warn("This should not happen.")
					elseif subsetRelationSomePointInsideGeoEdgesIds2 == "strictInside"
						subsetRelationGeoEdgesIds1 = "touchingSuperset"
						subsetRelationGeoEdgesIds2 = "touchingSubset"
					elseif subsetRelationSomePointInsideGeoEdgesIds2 == "strictOutside"
						subsetRelationGeoEdgesIds1 = "touchingDisjoint"
						subsetRelationGeoEdgesIds2 = "touchingDisjoint"
					end
				end
				
			end
			
		end
		
	end
	
	return subsetRelationGeoEdgesIds1, subsetRelationGeoEdgesIds2
	# "identical", "touchingSubset", "strictSubset", "touchingSuperset", "strictSuperset", "intersecting", "touchingDisjoint", "strictDisjoint"
end


