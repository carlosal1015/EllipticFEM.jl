
function equationDataGammasRotateGammaNodes(centerCoords::Array{Float64,1},angle::Float64,gammaNodes::Array{GammaNode,1})
	
	for gammaNode in gammaNodes
		gammaNode.coords = rotate2DPoint(centerCoords,angle,gammaNode.coords)
	end
	
	return nothing
end
