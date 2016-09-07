
include("./geoDataToGeoFile.jl")
include("./geoFileToMshFile.jl")
include("./mshFileToMeshData.jl")

function geoDataToMeshData(geoData::GeoData,intendedMeshsize::Float64)
	
	filepath = "./Gmsh saved files/gmshData"
	
	# Write .geo file:
	geoDataToGeoFile(geoData,filepath*".geo")
	
	# Run gmsh:
	geoFileToMshFile(filepath*".geo",intendedMeshsize)
	
	# Read .msh file:
	meshData = mshFileToMeshData(filepath*".msh")
	
	# Initialize additional meshData fields:
	meshDataElementsInitializeNablaFTsVolumesCentroidsAndMeshsizes(meshData)
	meshDataEdgesInitializeParents(meshData,geoData)
	meshDataElementsInitializeParents(meshData,geoData)
	
	return meshData
end


