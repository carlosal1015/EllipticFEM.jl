
function geoFileToMshFile(path::ASCIIString,intendedMeshsize::Float64)
	readall(`gmsh $path -2 -clmax $intendedMeshsize`)
	return nothing
end
