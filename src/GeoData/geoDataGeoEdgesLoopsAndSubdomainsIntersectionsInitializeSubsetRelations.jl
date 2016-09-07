
function geoDataGeoEdgesLoopsInitializeSubsetRelationsToGeoEdgesLoops(geoData::GeoData)
	
	for geoEdgesLoop in geoData.geoEdgesLoops
		geoEdgesLoop.subsetRelationsToGeoEdgesLoops = Array{ASCIIString}(length(geoData.geoEdgesLoops))
	end
	
	for geoEdgesLoop1 in geoData.geoEdgesLoops
		for geoEdgesLoop2 in geoData.geoEdgesLoops[geoEdgesLoop1.id:1:end]
			
			geoEdgesLoop1SubsetRelation,geoEdgesLoop2SubsetRelation = 
				geoDataComputSubsetRelationGeoEdgesIdsGeoEdgesIds(geoData,
					"geoEdgesLoop",geoEdgesLoop1.id,geoEdgesLoop1.boundingBox,geoEdgesLoop1.geoEdgesIds,
					"geoEdgesLoop",geoEdgesLoop2.id,geoEdgesLoop2.boundingBox,geoEdgesLoop2.geoEdgesIds)
			
			geoEdgesLoop1.subsetRelationsToGeoEdgesLoops[geoEdgesLoop2.id] = geoEdgesLoop1SubsetRelation
			geoEdgesLoop2.subsetRelationsToGeoEdgesLoops[geoEdgesLoop1.id] = geoEdgesLoop2SubsetRelation
		end
	end
	
	return nothing
end

function geoDataSubdomainsIntersectionsInitializeSubsetRelationsToSubdomainsIntersections(geoData::GeoData)
	
	for subdomainsIntersection in geoData.subdomainsIntersections
		subdomainsIntersection.subsetRelationsToSubdomainsIntersections = Array{ASCIIString}(length(geoData.subdomainsIntersections))
	end
	
	for subdomainsIntersection1 in geoData.subdomainsIntersections
		for subdomainsIntersection2 in geoData.subdomainsIntersections[subdomainsIntersection1.id:1:end]
			
			subdomainsIntersection1SubsetRelation,subdomainsIntersection2SubsetRelation = 
				geoDataComputSubsetRelationGeoEdgesIdsGeoEdgesIds(geoData,
					"subdomainsIntersection",subdomainsIntersection1.id,subdomainsIntersection1.boundingBox,subdomainsIntersection1.geoEdgesIds,
					"subdomainsIntersection",subdomainsIntersection2.id,subdomainsIntersection2.boundingBox,subdomainsIntersection2.geoEdgesIds)
			
			subdomainsIntersection1.subsetRelationsToSubdomainsIntersections[subdomainsIntersection2.id] = subdomainsIntersection1SubsetRelation
			subdomainsIntersection2.subsetRelationsToSubdomainsIntersections[subdomainsIntersection1.id] = subdomainsIntersection2SubsetRelation
		end
	end
	
	return nothing
end

function geoDataSubdomainsIntersectionsInitializeSubsetRelationsToGeoEdgesLoops(geoData::GeoData)
	
	for subdomainsIntersection in geoData.subdomainsIntersections
		subdomainsIntersection.subsetRelationsToGeoEdgesLoops = Array{ASCIIString}(length(geoData.geoEdgesLoops))
	end
	
	for subdomainsIntersection in geoData.subdomainsIntersections
		for geoEdgesLoop in geoData.geoEdgesLoops
			
			subdomainsIntersectionSubsetRelation,geoEdgesLoopSubsetRelation = 
				geoDataComputSubsetRelationGeoEdgesIdsGeoEdgesIds(geoData,
					"subdomainsIntersection",subdomainsIntersection.id,subdomainsIntersection.boundingBox,subdomainsIntersection.geoEdgesIds,
					"geoEdgesLoop",geoEdgesLoop.id,geoEdgesLoop.boundingBox,geoEdgesLoop.geoEdgesIds)
			
			subdomainsIntersection.subsetRelationsToGeoEdgesLoops[geoEdgesLoop.id] = subdomainsIntersectionSubsetRelation
		end
	end
	
	return nothing
end

