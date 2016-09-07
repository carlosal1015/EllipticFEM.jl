
function geoDataInitializePeriodicGeoEdgesPairs(equationData::EquationData,geoData::GeoData)
	
	mastersGeoEdgesIds_all = Array{Int64,1}(0)
	slavesGeoEdgesIds_all = Array{Int64,1}(0)
	for domain in vcat(equationData.domain,equationData.holes,equationData.subdomains)
		for periodicGammasPair in domain.periodicGammasPairs
			
			masterGammaId = periodicGammasPair.gammasIds[1]
			slaveGammaId = periodicGammasPair.gammasIds[2]
			
			mastersGeoEdgesIds = find(map(
				geoEdge->(geoEdge.parentDomainType==domain.domainType)&(geoEdge.parentDomainId==domain.id)&(geoEdge.parentGammaId==masterGammaId),
				geoData.geoEdges))
			
			if periodicGammasPair.useSameOrientation
				slavesGeoEdgesIds = find(map(
					geoEdge->(geoEdge.parentDomainType==domain.domainType)&(geoEdge.parentDomainId==domain.id)&(geoEdge.parentGammaId==slaveGammaId),
					geoData.geoEdges))
			else
				slavesGeoEdgesIds = -flipdim(find(map(
					geoEdge->(geoEdge.parentDomainType==domain.domainType)&(geoEdge.parentDomainId==domain.id)&(geoEdge.parentGammaId==slaveGammaId),
					geoData.geoEdges)),1)
			end
			
			if length(mastersGeoEdgesIds) == length(slavesGeoEdgesIds)
				append!(mastersGeoEdgesIds_all,mastersGeoEdgesIds)
				append!(slavesGeoEdgesIds_all,slavesGeoEdgesIds)
			else
				warn("length(mastersGeoEdgesIds) != length(slavesGeoEdgesIds)")
			end
			
		end
	end
	
	periodicGeoEdgesPairs = Array{PeriodicGeoEdgesPair}(length(mastersGeoEdgesIds_all))
	for i = 1:1:length(periodicGeoEdgesPairs)
		periodicGeoEdgesPairs[i] = PeriodicGeoEdgesPair(i,[mastersGeoEdgesIds_all[i],slavesGeoEdgesIds_all[i]],true)
	end
	
	geoData.periodicGeoEdgesPairs = periodicGeoEdgesPairs
	
	return nothing
end
