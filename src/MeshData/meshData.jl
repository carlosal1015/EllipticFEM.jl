
type MeshData
	
	# Stores all data of the triangulation of the domain \Omega.
	
	nodes::Array{Float64,2}
	
	edges::Array{Int64,2}
	edgesParentGeoEdgesIds::Array{Int64}
	edgesParentDomainsTypes::Array{ASCIIString}
	edgesParentDomainsIds::Array{Int64}
	edgesParentGammasIds::Array{Int64}
	edgesParentGammasParameterRanges::Array{Float64,2}
	edgesVolumes::Array{Float64}
	
	elements::Array{Int64,2}
	elementsParentSubdomainsIntersectionsIds::Array{Int64}
	elementsParentDomainsTypes::Array{ASCIIString}
	elementsParentDomainsIds::Array{Int64}
	elementsDetNablaFTs::Array{Float64}
	elementsNablaFTInv11s::Array{Float64}
	elementsNablaFTInv21s::Array{Float64}
	elementsNablaFTInv12s::Array{Float64}
	elementsNablaFTInv22s::Array{Float64}
	elementsVolumes::Array{Float64}
	elementsCentroids::Array{Float64,2}
	elementsMeshsizes::Array{Float64}
	
	periodicNodes::Array{Int64,2}
end

include("./meshDataDisplayAndPlot.jl")
include("./meshDataEdgesAndElementsInitializeParents.jl")
include("./meshDataElementsInitializeNablaFTsVolumesCentroidsAndMeshsizes.jl")


