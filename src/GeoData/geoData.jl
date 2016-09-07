
# The geometry data defined in the Gammas of equationData are discretized and stored in the
# type GeoData. Then geoData is used to build a mesh. The type GeoData mimicks the respective
# types used by Gmsh, but adds some additional information.

type GeoNode
	id::Int64
	coords::Array{Float64}
	isActive::Bool
end

type GeoEdge
	id::Int64
	geoNodesIds::Array{Int64}
	
	parentDomainType::ASCIIString
	parentDomainId::Int64
	parentGammaId::Int64
	parentGammaParameterRange::Array{Float64}
	
	isActive::Bool
end

type GeoEdgesLoop
	id::Int64
	geoEdgesIds::Array{Int64}
	
	parentDomainType::ASCIIString
	parentDomainId::Int64
	
	orientation::ASCIIString	# "positive" (counter-clockwise), "negative" (clockwise)
	boundingBox::Array{Float64}		# [xMin,xMax,yMin,yMax]
	
	subsetRelationsToGeoEdgesLoops::Array{ASCIIString}		# "identical", "touchingSubset", "strictSubset", "touchingSuperset", "strictSuperset", "intersecting", "touchingDisjoint", "strictDisjoint"
end

type GeoNodesIntersection
	id::Int64
	geoNodesIds::Array{Int64}
end

type SubdomainsIntersection
	id::Int64
	
	starterGeoEdgesLoopId::Int64
	starterGeoEdgesLoopGeoEdgesIdsIndex::Int64
	geoEdgesIds::Array{Int64}
	
	orientation::ASCIIString	# "positive" (counter-clockwise), "negative" (clockwise)
	boundingBox::Array{Float64}		# [xMin,xMax,yMin,yMax]
	
	enclosedSubdomainsIntersectionsIds::Array{Int64}
	
	subsetRelationsToSubdomainsIntersections::Array{ASCIIString}	# "identical", "touchingSubset", "strictSubset", "touchingSuperset", "strictSuperset", "intersecting", "touchingDisjoint", "strictDisjoint"
	subsetRelationsToGeoEdgesLoops::Array{ASCIIString}		# "identical", "touchingSubset", "strictSubset", "touchingSuperset", "strictSuperset", "intersecting", "touchingDisjoint", "strictDisjoint"
	parentGeoEdgesLoopId::Int64
	
	isActive::Bool
end

type PeriodicGeoEdgesPair
	id::Int64
	geoEdgesIds::Array{Int64}		# [masterId,slaveId]
	isActive::Bool
end

type GeoData
	geoNodes::Array{GeoNode}
	geoEdges::Array{GeoEdge}
	geoEdgesLoops::Array{GeoEdgesLoop}
	
	geoNodesIntersections::Array{GeoNodesIntersection}
	subdomainsIntersections::Array{SubdomainsIntersection}
	
	periodicGeoEdgesPairs::Array{PeriodicGeoEdgesPair}
end

include("./geoDataComputeSubsetRelations.jl")
include("./geoDataDisplayAndPlot.jl")
include("./geoDataGeoEdgesLoopsAndSubdomainsIntersectionsInitializeBoundingBoxes.jl")
include("./geoDataGeoEdgesLoopsAndSubdomainsIntersectionsInitializeOrientations.jl")
include("./geoDataGeoEdgesLoopsAndSubdomainsIntersectionsInitializeSubsetRelations.jl")
include("./geoDataGeoNodesAndGeoEdgesInitializeNonIsActives.jl")
include("./geoDataInitializePeriodicGeoEdgesPairs.jl")
include("./geoDataInitializeSubdomainsIntersections.jl")
include("./geoDataSubdomainsIntersectionsInitializeEnclosedSubdomainsIntersectionsIds.jl")
include("./geoDataSubdomainsIntersectionsInitializeParentGeoEdgesLoopIds.jl")
include("./geometricAlgorithms.jl")


