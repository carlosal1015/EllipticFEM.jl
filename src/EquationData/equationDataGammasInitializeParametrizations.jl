
function equationDataGammasInitializeParametrizations(equationData::EquationData)
	for subdomain in vcat(equationData.domain,equationData.holes,equationData.subdomains)
		for gamma in subdomain.gammas
			
			if gamma.parametrizationType == "gammaLine"
				gamma.parametrization = equationDataGammasDefineGammaLineParametrization(
					subdomain.gammaNodes[gamma.parametrizationGammaNodesIds[1]],
					subdomain.gammaNodes[gamma.parametrizationGammaNodesIds[2]])
			elseif gamma.parametrizationType == "gammaArc"
				gamma.parametrization = equationDataGammasDefineGammaArcParametrization(
					subdomain.gammaNodes[gamma.parametrizationGammaNodesIds[1]],
					subdomain.gammaNodes[gamma.parametrizationGammaNodesIds[2]],
					gamma.parametrizationCoefficients[1])
			elseif gamma.parametrizationType == "gammaParabolaArc"
				gamma.parametrization = equationDataGammasDefineGammaParabolaArcParametrization(
					subdomain.gammaNodes[gamma.parametrizationGammaNodesIds[1]],
					subdomain.gammaNodes[gamma.parametrizationGammaNodesIds[2]],
					gamma.parametrizationCoefficients[1],
					gamma.parametrizationCoefficients[2])
			end
			
		end
	end
	
	return nothing
end

function equationDataGammasDefineGammaLineParametrization(gammaNode1::GammaNode,gammaNode2::GammaNode)
	
	gamma = (s) -> evaluateConvexCombination(gammaNode1.coords,gammaNode2.coords,s)
	
	return gamma
end

function equationDataGammasDefineGammaArcParametrization(gammaNodeCenter::GammaNode,gammaNode1::GammaNode,angle::Float64)
	
	gamma = (s) -> rotate2DPoint(gammaNodeCenter.coords,angle*s,gammaNode1.coords)
	
	return gamma
end

function equationDataGammasDefineGammaParabolaArcParametrization(gammaNodeFocus::GammaNode,gammaNodeVertex::GammaNode,angleMin::Float64,angleMax::Float64)
	
	gamma = (s) -> gammaNodeFocus.coords + 2*(gammaNodeVertex.coords-gammaNodeFocus.coords) + 
		tan(evaluateConvexCombination(angleMin,angleMax,s)/2)*hcat([0,-1],[1,0])*2*(gammaNodeFocus.coords-gammaNodeVertex.coords) + 
		(1/(2*norm(gammaNodeVertex.coords-gammaNodeFocus.coords)))*(tan(evaluateConvexCombination(angleMin,angleMax,s)/2)^2 + norm(gammaNodeVertex.coords-gammaNodeFocus.coords))*2*(gammaNodeFocus.coords-gammaNodeVertex.coords)
	
	return gamma
end


