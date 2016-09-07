
function lseDataSolve(equationData::EquationData,meshData::MeshData,lseData::LSEData)
	
	if in(lseData.solverType,["standardSolverElliptic","saddlepointSolverElliptic"])
		
		lseData.solutionsTimestamps = [0]
		if lseData.solverType == "standardSolverElliptic"
			lseData.solutions = Vector{Float64}[lseDataSolveStandardSystemElliptic(lseData.M,lseData.r)]
		elseif lseData.solverType == "saddlepointSolverElliptic"
			lseData.solutions = Vector{Float64}[lseDataSolveSaddlepointSystemElliptic(lseData.M,lseData.integralValues,lseData.r,0.0)]
		end
		
	elseif lseData.solverType == "solverParabolic"
		
		numberOfTimesteps = round(Int64,(equationData.tMax-0)/lseData.deltaT)
		
		solutionsTimestamps = Array{Float64,1}(numberOfTimesteps+1)
		solutionsTimestamps[1] = 0
		
		solutions = Array{Array{Float64,1},1}(numberOfTimesteps+1)
		uInitialValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"uInitial",0.0)
		uInitialValuesAtNodes = valuesAtElementsCentroidsToValuesAtNodes(meshData,uInitialValuesAtCentroids)
		solutions[1] = uInitialValuesAtNodes
		
		for i = 2:1:(numberOfTimesteps+1)
			
			t = (lseData.deltaT)*(i-1)
			
			println("t = "*string(round(t,2)))
			
			lseDataUpdate(lseData,equationData,meshData,t)
			
			solutionsTimestamps[i] = t
			solutions[i] = lseDataSolveSystemParabolic(lseData.M,lseData.MMass,lseData.r,solutions[i-1],lseData.deltaT)
		end
		
		lseData.solutionsTimestamps = solutionsTimestamps
		lseData.solutions = solutions
		
	elseif lseData.solverType == "solverHyperbolic"
		
		numberOfTimesteps = round(Int64,(equationData.tMax-0)/lseData.deltaT)
		
		solutionsTimestamps = Array{Float64,1}(numberOfTimesteps+1)
		solutionsTimestamps[1] = 0
		solutionsTimestamps[2] = lseData.deltaT
		
		solutions = Array{Array{Float64,1},1}(numberOfTimesteps+1)
		uInitialValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"uInitial",0.0)
		uInitialValuesAtNodes = valuesAtElementsCentroidsToValuesAtNodes(meshData,uInitialValuesAtCentroids)
		solutions[1] = uInitialValuesAtNodes
		uTInitialValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"uTInitial",0.0)
		uTInitialValuesAtNodes = valuesAtElementsCentroidsToValuesAtNodes(meshData,uTInitialValuesAtCentroids)
		solutions[2] = uInitialValuesAtNodes + lseData.deltaT*uTInitialValuesAtNodes
		
		for i = 3:1:(numberOfTimesteps+1)
			
			t = (lseData.deltaT)*(i-1)
			
			println("t = "*string(round(t,2)))
			
			lseDataUpdate(lseData,equationData,meshData,t)
			
			solutionsTimestamps[i] = t
			solutions[i] = lseDataSolveSystemHyperbolic(lseData.M,lseData.MMass,lseData.r,solutions[i-1],solutions[i-2],lseData.deltaT)
		end
		
		lseData.solutionsTimestamps = solutionsTimestamps
		lseData.solutions = solutions
		
	end
	
	return nothing
end

function lseDataSolveStandardSystemElliptic(M::SparseMatrixCSC{Float64,Int64},r::Array{Float64,1})
	
	solution = M\r
	
	return solution
end

function lseDataSolveSaddlepointSystemElliptic(M::SparseMatrixCSC{Float64,Int64},N::Array{Float64,1},r::Array{Float64,1},s::Float64)
	
	# Solves the saddlepoint system
	#
	# hcat([M,N^T],[N,0])*[solution,lambda] = [r,s],
	#
	# where N \in \R^{n x 1}, s \in \R^1 and lambda \in \R^1.
	# (In this case one can solve for lambda directly!)
	
	MInvRAndMInvN = M\hcat(r,N)
	MInvR = collect(MInvRAndMInvN[:,1])
	MInvN = collect(MInvRAndMInvN[:,2])
	
	lambda = (dot(MInvR,N) - s)/dot(MInvN,N)
	
	solution = MInvR - lambda*MInvN
	
	return solution
end

function lseDataSolveSystemParabolic(M::SparseMatrixCSC{Float64,Int64},MMass::SparseMatrixCSC{Float64,Int64},r::Array{Float64,1},solutionCurrent::Array{Float64,1},deltaT::Float64)
	
	# Implicit Euler algorithm:
	# (\partial_t u)(t+\deltaT) = (1/deltaT)*[u(t+deltaT) - u(t)] + O(deltaT)
	
	solutionNext = (MMass + deltaT*M)\(deltaT*r + MMass*solutionCurrent)
	
	return solutionNext
end

function lseDataSolveSystemHyperbolic(M::SparseMatrixCSC{Float64,Int64},MMass::SparseMatrixCSC{Float64,Int64},r::Array{Float64,1},solutionCurrent::Array{Float64,1},solutionOld::Array{Float64,1},deltaT::Float64)
	
	# Implicit Euler algorithm:
	# (\partial_{tt} u)(t+\deltaT) = (1/\deltaT^2)*[u(t+\deltaT) - 2*u(t) + u(t-\deltaT)] + O(deltaT)
	
	solutionNext = (MMass + deltaT^2*M)\(deltaT^2*r + MMass*(2*solutionCurrent-solutionOld))
	
	return solutionNext
end

