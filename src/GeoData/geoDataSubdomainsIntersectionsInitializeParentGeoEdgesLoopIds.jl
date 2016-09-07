
function geoDataSubdomainsIntersectionsInitializeParentGeoEdgesLoopIds(geoData::GeoData)
	
	for subdomainsIntersection in geoData.subdomainsIntersections
		if subdomainsIntersection.isActive
			
			index1 = findfirst(map(geoEdgesLoop->
				(geoEdgesLoop.parentDomainType=="subdomain")&
				in(subdomainsIntersection.subsetRelationsToGeoEdgesLoops[geoEdgesLoop.id],["identical","touchingSubset","strictSubset"]),
				geoData.geoEdgesLoops))
			
			if index1 != 0
				subdomainsIntersection.parentGeoEdgesLoopId = index1
			else
				
				index2 = findfirst(map(geoEdgesLoop->
					(geoEdgesLoop.parentDomainType=="hole")&
					in(subdomainsIntersection.subsetRelationsToGeoEdgesLoops[geoEdgesLoop.id],["identical","touchingSubset","strictSubset"]),
					geoData.geoEdgesLoops))
				
				if index2 != 0
					subdomainsIntersection.parentGeoEdgesLoopId = index2
				else
					subdomainsIntersection.parentGeoEdgesLoopId = 1
				end
			end
		end
	end
	
	return nothing
end


