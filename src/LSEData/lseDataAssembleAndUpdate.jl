
function lseDataAssemble(equationData::EquationData,meshData::MeshData,deltaT::Float64)
	
	# This function assembles the system matrix M and the right-hand side r
	# of the LSE occurring in the P1-FEM (piecewise linear test functions on the given meshData) 
	# for the equation described by equationData.
	
	# Determine types of given boundary data:
	boundaryDataAvailableNeumann = any(map(tmpDomain->any(map(gamma->(gamma.boundaryDataType=="Neumann"),tmpDomain.gammas)),vcat(equationData.domain,equationData.holes)))
	boundaryDataAvailableDirichlet = any(map(tmpDomain->any(map(gamma->(gamma.boundaryDataType=="Dirichlet"),tmpDomain.gammas)),vcat(equationData.domain,equationData.holes)))
	boundaryDataAvailablePeriodic = any(map(tmpDomain->any(map(gamma->(gamma.boundaryDataType=="periodic"),tmpDomain.gammas)),vcat(equationData.domain,equationData.holes)))
	
	# Flags for quickly computable parts of LSE assembly:
	doQuickAssemblyMA = all(map(tmpDomain->tmpDomain.isIdentityA,vcat(equationData.domain,equationData.subdomains)))
	doQuickAssemblyMb = all(map(tmpDomain->tmpDomain.isZeroB,vcat(equationData.domain,equationData.subdomains)))
	doQuickAssemblyMc = all(map(tmpDomain->tmpDomain.isZeroC,vcat(equationData.domain,equationData.subdomains)))
	doQuickAssemblyRF = all(map(tmpDomain->tmpDomain.isZeroF,vcat(equationData.domain,equationData.subdomains)))
	doQuickAssemblyRFDiv = all(map(tmpDomain->tmpDomain.isZeroFDiv,vcat(equationData.domain,equationData.subdomains)))
	doQuickAssemblyRNeumann = all(map(tmpDomain->all(map(gamma->(gamma.boundaryDataType!="Neumann")|gamma.isZeroBoundaryDataFunction,tmpDomain.gammas)),vcat(equationData.domain,equationData.holes)))
	# doQuickAssemblyRDirichlet not available.
	
	# Flags for omittable parts of LSE assembly:
	doAssembleDirichletNodesIds = boundaryDataAvailableDirichlet
	doAssembleGDirichletValuesAtDirichletNodes = boundaryDataAvailableDirichlet
	doAssembleRDirichlet = boundaryDataAvailableDirichlet
	doEnforceDirichletBoundaryConditionsInM = boundaryDataAvailableDirichlet
	doEnforceDirichletBoundaryConditionsInR = boundaryDataAvailableDirichlet
	doEnforcePeriodicBoundaryConditionsInM = boundaryDataAvailablePeriodic
	doEnforcePeriodicBoundaryConditionsInR = boundaryDataAvailablePeriodic
	doAssembleMMass = in(equationData.equationType,["parabolic","hyperbolic"])
	doAssembleIntegralValues = !boundaryDataAvailableDirichlet
	
	# Assemble phiHats:
	phiHats = lseDataAssemblePhiHats()
	
	# Assemble MA:
	MA = lseDataAssembleMA(doQuickAssemblyMA,equationData,meshData,phiHats,0.0)
	
	# Assemble Mb:
	Mb = lseDataAssembleMb(doQuickAssemblyMb,equationData,meshData,phiHats,0.0)
	
	# Assemble Mc:
	Mc = lseDataAssembleMc(doQuickAssemblyMc,equationData,meshData,phiHats,0.0)
	
	# Assemble MAbc:
	MAbc = lseDataAssembleMAbc(MA,Mb,Mc)
	
	# Assemble M:
	M = lseDataAssembleM(MAbc)
	
	# Assemble rF:
	rF = lseDataAssembleRF(doQuickAssemblyRF,equationData,meshData,phiHats,0.0)
	
	# Assemble rFDiv:
	rFDiv = lseDataAssembleRFDiv(doQuickAssemblyRFDiv,equationData,meshData,phiHats,0.0)
	
	# Assemble rNeumann:
	rNeumann = lseDataAssembleRNeumann(doQuickAssemblyRNeumann,equationData,meshData,0.0)
	
	# Assemble dirichletNodesIds,gDirichletValuesAtDirichletNodes and rDirichlet:
	if doAssembleDirichletNodesIds|doAssembleGDirichletValuesAtDirichletNodes
		dirichletNodesIds,gDirichletValuesAtDirichletNodes = 
			lseDataAssembleDirichletNodesIdsAndGDirichletValuesAtDirichletNodes(equationData,meshData,0.0)
	else
		dirichletNodesIds = []
		gDirichletValuesAtDirichletNodes = []
	end
	if doAssembleRDirichlet
		rDirichlet = lseDataAssembleRDirichlet(meshData,MAbc,dirichletNodesIds,gDirichletValuesAtDirichletNodes)
	else
		rDirichlet = Array{Float64,1}(0)
	end
	
	# Assemble r:
	r = lseDataAssembleR(rF,rFDiv,rNeumann,rDirichlet)
	
	# Enforce Dirichlet boundary conditions in system matrix and right-hand side:
	if doEnforceDirichletBoundaryConditionsInM
		enforceDirichletBoundaryConditionsInM(M,dirichletNodesIds)
	end
	if doEnforceDirichletBoundaryConditionsInR
		enforceDirichletBoundaryConditionsInR(r,dirichletNodesIds,gDirichletValuesAtDirichletNodes)
	end
	
	# Enforce periodic boundary conditions in system matrix and right-hand side:
	if doEnforcePeriodicBoundaryConditionsInM
		enforcePeriodicBoundaryConditionsInM(M,meshData)
	end
	if doEnforcePeriodicBoundaryConditionsInR
		enforcePeriodicBoundaryConditionsInR(r,meshData)
	end
	
	# Assemble MMass:
	if doAssembleMMass
		MMass = lseDataAssembleMMass(meshData,phiHats)
	else
		MMass = spzeros(0,0)
	end
	
	# Assemble integralValues:
	if doAssembleIntegralValues
		integralValues = lseDataAssembleIntegralValues(meshData,phiHats)
	else
		integralValues = []
	end
	
	# Infer necessary solver type:
	if equationData.equationType == "elliptic"
		if boundaryDataAvailableDirichlet
			solverType = "standardSolverElliptic"
		else
			solverType = "saddlepointSolverElliptic"
		end
	elseif equationData.equationType == "parabolic"
		solverType = "solverParabolic"
	elseif equationData.equationType == "hyperbolic"
		solverType = "solverHyperbolic"
	end
	
	# Assemble lseData:
	lseData = LSEData(phiHats,
		MA,Mb,Mc,MAbc,M,
		rF,rFDiv,rNeumann,dirichletNodesIds,gDirichletValuesAtDirichletNodes,rDirichlet,r,
		MMass,integralValues,
		solverType,
		deltaT,[],[])
	
	return lseData
end

function lseDataUpdate(lseData::LSEData,equationData::EquationData,meshData::MeshData,t::Float64)
	
	# Flags for omittable parts of LSE update:
	doUpdateMA = any(map(tmpDomain->!tmpDomain.isConstantInTimeA,vcat(equationData.domain,equationData.subdomains)))
	doUpdateMb = any(map(tmpDomain->!tmpDomain.isConstantInTimeB,vcat(equationData.domain,equationData.subdomains)))
	doUpdateMc = any(map(tmpDomain->!tmpDomain.isConstantInTimeC,vcat(equationData.domain,equationData.subdomains)))
	doUpdateMAbc = doUpdateMA|doUpdateMb|doUpdateMc
	doUpdateM = doUpdateMAbc
	doUpdateRF = any(map(tmpDomain->!tmpDomain.isConstantInTimeF,vcat(equationData.domain,equationData.subdomains)))
	doUpdateRFDiv = any(map(tmpDomain->!tmpDomain.isConstantInTimeFDiv,vcat(equationData.domain,equationData.subdomains)))
	doUpdateRNeumann = any(map(tmpDomain->any(map(gamma->(gamma.boundaryDataType=="Neumann")&(!gamma.isConstantInTimeBoundaryDataFunction),tmpDomain.gammas)),vcat(equationData.domain,equationData.holes)))
	doUpdateDirichletNodesIds = false
	doUpdateGDirichletValuesAtDirichletNodes = any(map(tmpDomain->any(map(gamma->(gamma.boundaryDataType=="Dirichlet")&(!gamma.isConstantInTimeBoundaryDataFunction),tmpDomain.gammas)),vcat(equationData.domain,equationData.holes)))
	doUpdateRDirichlet = doUpdateMAbc|doUpdateDirichletNodesIds|doUpdateGDirichletValuesAtDirichletNodes
	doUpdateR = doUpdateRF|doUpdateRFDiv|doUpdateRNeumann|doUpdateRDirichlet
	
	doEnforceDirichletBoundaryConditionsInM = doUpdateM|doUpdateDirichletNodesIds
	doEnforceDirichletBoundaryConditionsInR = doUpdateR|doUpdateDirichletNodesIds|doUpdateGDirichletValuesAtDirichletNodes
	doEnforcePeriodicBoundaryConditionsInM = doUpdateM
	doEnforcePeriodicBoundaryConditionsInR = doUpdateR
	
	# Update MA:
	if doUpdateMA
		lseData.MA = lseDataAssembleMA(false,equationData,meshData,lseData.phiHats,t)
	end
	
	# Update Mb:
	if doUpdateMb
		lseData.Mb = lseDataAssembleMb(false,equationData,meshData,lseData.phiHats,t)
	end
	
	# Update Mc:
	if doUpdateMc
		lseData.Mc = lseDataAssembleMc(false,equationData,meshData,lseData.phiHats,t)
	end
	
	# Update MAbc:
	if doUpdateMAbc
		lseData.MAbc = lseDataAssembleMAbc(lseData.MA,lseData.Mb,lseData.Mc)
	end
	
	# Update M:
	if doUpdateM
		lseData.M = lseDataAssembleM(lseData.MAbc)
	end
	
	# Update rF:
	if doUpdateRF
		lseData.rF = lseDataAssembleRF(false,equationData,meshData,lseData.phiHats,t)
	end
	
	# Update rFDiv:
	if doUpdateRFDiv
		lseData.rFDiv = lseDataAssembleRFDiv(false,equationData,meshData,lseData.phiHats,t)
	end
	
	# Update rNeumann:
	if doUpdateRNeumann
		lseData.rNeumann = lseDataAssembleRNeumann(false,equationData,meshData,t)
	end
	
	# Update dirichletNodesIds and gDirichletValuesAtDirichletNodes:
	if doUpdateDirichletNodesIds|doUpdateGDirichletValuesAtDirichletNodes
		lseData.dirichletNodesIds,lseData.gDirichletValuesAtDirichletNodes = 
			lseDataAssembleDirichletNodesIdsAndGDirichletValuesAtDirichletNodes(equationData,meshData,t)
	end
	
	# Update rDirichlet:
	if doUpdateRDirichlet
		lseData.rDirichlet = lseDataAssembleRDirichlet(meshData,lseData.MAbc,lseData.dirichletNodesIds,lseData.gDirichletValuesAtDirichletNodes)
	end
	
	# Update r:
	if doUpdateR
		lseData.r = lseDataAssembleR(lseData.rF,lseData.rFDiv,lseData.rNeumann,lseData.rDirichlet)
	end
	
	# Enforce Dirichlet boundary conditions in system matrix and right-hand side:
	if doEnforceDirichletBoundaryConditionsInM
		enforceDirichletBoundaryConditionsInM(lseData.M,lseData.dirichletNodesIds)
	end
	if doEnforceDirichletBoundaryConditionsInR
		enforceDirichletBoundaryConditionsInR(lseData.r,lseData.dirichletNodesIds,lseData.gDirichletValuesAtDirichletNodes)
	end
	
	# Enforce periodic boundary conditions in system matrix and right-hand side:
	if doEnforcePeriodicBoundaryConditionsInM
		enforcePeriodicBoundaryConditionsInM(lseData.M,meshData)
	end
	if doEnforcePeriodicBoundaryConditionsInR
		enforcePeriodicBoundaryConditionsInR(lseData.r,meshData)
	end
	
	return nothing
end

################################################

function lseDataAssemblePhiHats()
	
	PhiHat0 = (1/6)*[1,1,1]
	PhiHatX = (1/2)*[-1,1,0]
	PhiHatY = (1/2)*[-1,0,1]
	PhiHat00 = (1/24)*hcat([2,1,1],[1,2,1],[1,1,2])
	PhiHatX0 = (1/6)*hcat([-1,1,0],[-1,1,0],[-1,1,0])
	PhiHatY0 = (1/6)*hcat([-1,0,1],[-1,0,1],[-1,0,1])
	PhiHat0X = (1/6)*hcat([-1,-1,-1],[1,1,1],[0,0,0])
	PhiHat0Y = (1/6)*hcat([-1,-1,-1],[0,0,0],[1,1,1])
	PhiHatXX = (1/2)*hcat([1,-1,0],[-1,1,0],[0,0,0])
	PhiHatXY = (1/2)*hcat([1,-1,0],[0,0,0],[-1,1,0])
	PhiHatYX = (1/2)*hcat([1,0,-1],[-1,0,1],[0,0,0])
	PhiHatYY = (1/2)*hcat([1,0,-1],[0,0,0],[-1,0,1])
	
	phiHats = PhiHats(
		PhiHat0,
		PhiHatX,
		PhiHatY,
		PhiHat00,
		PhiHatX0,
		PhiHatY0,
		PhiHat0X,
		PhiHat0Y,
		PhiHatXX,
		PhiHatXY,
		PhiHatYX,
		PhiHatYY)
	
	return phiHats
end

################################################

function lseDataAssembleMA(doQuickAssemblyMA::Bool,equationData::EquationData,meshData::MeshData,phiHats::PhiHats,t::Float64)
	
	# Assembles the Matrix
	#
	# MA := (\int_{\Omega} \nabla v_M^T*A*\nabla v_N dx)_{M,N \in \Nodes},
	#
	# where v_M \in \S^{1,1}(\Elements) is the piecewise linear hat function belonging to the node M.
	
	# We implement the following formula (1-point approximation of A):
	#
	# MA = sum_{i,j=1}^{3} \sum_{T \in \Elements} abs(\det\nabla F_T)*
	#		[A11Tilde(T)*(\hat{\Phi}_{xx})_{ij} + A12Tilde(T)*(\hat{\Phi}_{xy})_{ij} +
	#		A21Tilde(T)*(\hat{\Phi}_{yx})_{ij} + A22Tilde(T)*(\hat{\Phi}_{yy})_{ij}]*
	#		e_{N_i(T)}*e_{N_j(T)}^T,
	#
	# where
	#
	# \hat{\Phi}_{...} \in \R^{3x3} ... integrals over the reference basis functions on the reference element,
	# e_{N_i(T)},e_{N_j(T)} \in \R^{#nodes} ... euclidean unit vectors,
	#
	# C_T \in \R^2 ... element centroid,
	# Aij(C_T) \in \R ... the value of the function A_ij at the point C_T,
	# Fij(T) := ((\nabla F_T)^{-1})_{ij} \in \R ... the entries of the inverse of the Jacobian matrix,
	# 		of the reference-element-mapping F_T : \hat{T} --> T,
	#
	# A11Tilde(T) := A11(C_T)*F11(T)^2 + A12(C_T)*F11(T)*F12(T) + A21(C_T)*F11(T)*F12(T) + A22(C_T)*F12(T)^2,
	# A12Tilde(T) := A11(C_T)*F11(T)*F21(T) + A12(C_T)*F11(T)*F22(T) + A21(C_T)*F12(T)*F21(T) + A22(C_T)*F12(T)*F22(T),
	# A21Tilde(T) := A11(C_T)*F11(T)*F21(T) + A12(C_T)*F12(T)*F21(T) + A21(C_T)*F11(T)*F22(T) + A22(C_T)*F12(T)*F22(T),
	# A22Tilde(T) := A11(C_T)*F21(T)^2 + A12(C_T)*F21(T)*F22(T) + A21(C_T)*F21(T)*F22(T) + A22(C_T)*F22(T)^2.
	
	# Compute all necessary values:
	if doQuickAssemblyMA
		
		A11TildeValuesAtCentroids = (meshData.elementsNablaFTInv11s).^2 + (meshData.elementsNablaFTInv12s).^2
		A12TildeValuesAtCentroids = (meshData.elementsNablaFTInv11s).*(meshData.elementsNablaFTInv21s) + 
			(meshData.elementsNablaFTInv12s).*(meshData.elementsNablaFTInv22s)
		A21TildeValuesAtCentroids = A12TildeValuesAtCentroids
		A22TildeValuesAtCentroids = (meshData.elementsNablaFTInv21s).^2 + (meshData.elementsNablaFTInv22s).^2
		
	else
		
		A11ValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"A11",t)
		A21ValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"A21",t)
		A12ValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"A12",t)
		A22ValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"A22",t)
		
		A11TildeValuesAtCentroids = 
			A11ValuesAtCentroids.*((meshData.elementsNablaFTInv11s).^2) + 
			A12ValuesAtCentroids.*(meshData.elementsNablaFTInv11s).*(meshData.elementsNablaFTInv12s) + 
			A21ValuesAtCentroids.*(meshData.elementsNablaFTInv11s).*(meshData.elementsNablaFTInv12s) + 
			A22ValuesAtCentroids.*((meshData.elementsNablaFTInv12s).^2)
		
		A12TildeValuesAtCentroids = 
			A11ValuesAtCentroids.*(meshData.elementsNablaFTInv11s).*(meshData.elementsNablaFTInv21s) + 
			A12ValuesAtCentroids.*(meshData.elementsNablaFTInv11s).*(meshData.elementsNablaFTInv22s) + 
			A21ValuesAtCentroids.*(meshData.elementsNablaFTInv12s).*(meshData.elementsNablaFTInv21s) + 
			A22ValuesAtCentroids.*(meshData.elementsNablaFTInv12s).*(meshData.elementsNablaFTInv22s)
		
		A21TildeValuesAtCentroids = 
			A11ValuesAtCentroids.*(meshData.elementsNablaFTInv11s).*(meshData.elementsNablaFTInv21s) + 
			A12ValuesAtCentroids.*(meshData.elementsNablaFTInv12s).*(meshData.elementsNablaFTInv21s) + 
			A21ValuesAtCentroids.*(meshData.elementsNablaFTInv11s).*(meshData.elementsNablaFTInv22s) + 
			A22ValuesAtCentroids.*(meshData.elementsNablaFTInv12s).*(meshData.elementsNablaFTInv22s)
		
		A22TildeValuesAtCentroids = 
			A11ValuesAtCentroids.*((meshData.elementsNablaFTInv21s).^2) + 
			A12ValuesAtCentroids.*(meshData.elementsNablaFTInv21s).*(meshData.elementsNablaFTInv22s) + 
			A21ValuesAtCentroids.*(meshData.elementsNablaFTInv21s).*(meshData.elementsNablaFTInv22s) + 
			A22ValuesAtCentroids.*((meshData.elementsNablaFTInv22s).^2)
		
	end
	
	MAIs = Array{Int64,2}(9,size(meshData.elements,2))
	MAJs = Array{Int64,2}(9,size(meshData.elements,2))
	MAVs = Array{Float64,2}(9,size(meshData.elements,2))
	for i = 1:1:3
		for j = 1:1:3
			MAIs[3*(i-1)+j,:] = meshData.elements[i,:]
			MAJs[3*(i-1)+j,:] = meshData.elements[j,:]
			MAVs[3*(i-1)+j,:] = abs(meshData.elementsDetNablaFTs).*(
				A11TildeValuesAtCentroids*(phiHats.PhiHatXX[i,j]) + 
				A12TildeValuesAtCentroids*(phiHats.PhiHatXY[i,j]) + 
				A21TildeValuesAtCentroids*(phiHats.PhiHatYX[i,j]) + 
				A22TildeValuesAtCentroids*(phiHats.PhiHatYY[i,j]))
		end
	end
	
	# Assemble sparse matrix:
	# If an index tuple [i,j] occurs multiple times, then sparse() sums up the corresponding values.
	MA = sparse(collect(MAIs),collect(MAJs),collect(MAVs),size(meshData.nodes,2),size(meshData.nodes,2))
	
	return MA
end

function lseDataAssembleMb(doQuickAssemblyMb::Bool,equationData::EquationData,meshData::MeshData,phiHats::PhiHats,t::Float64)
	
	# Assembles the Matrix
	#
	# Mb := (\int_{\Omega} (b^T*\nabla v_N)*v_M dx)_{M,N \in \Nodes},
	#
	# where v_M \in \S^{1,1}(\Elements) is the piecewise linear hat function belonging to the node M.
	
	# We implement the following formula (1-point approximation of b):
	#
	# Mb = sum_{i,j=1}^{3} \sum_{T \in \Elements} abs(\det\nabla F_T)*
	#		[b1Tilde(T)*(\hat{\Phi}_{x0})_{ji} + b2Tilde(T)*(\hat{\Phi}_{y0})_{ji}]*
	#		e_{N_i(T)}*e_{N_j(T)}^T,
	#
	# where
	#
	# \hat{\Phi}_{...} \in \R^{3x3} ... integrals over the reference basis functions on the reference element,
	# e_{N_i(T)},e_{N_j(T)} \in \R^{#nodes} ... euclidean unit vectors,
	#
	# C_T \in \R^2 ... element centroid,
	# bi(C_T) \in \R ... the value of the function b_i at the point C_T,
	# Fij(T) := ((\nabla F_T)^{-1})_{ij} \in \R ... the entries of the inverse of the Jacobian matrix,
	# 		of the reference-element-mapping F_T : \hat{T} --> T,
	#
	# b1Tilde(T) := b1(C_T)*F11(T)^2 + b2(C_T)*F12(T),
	# b2Tilde(T) := b1(C_T)*F21(T)^2 + b2(C_T)*F22(T).
	
	if doQuickAssemblyMb
		
		Mb = spzeros(size(meshData.nodes,2),size(meshData.nodes,2))
		
	else
		
		# Compute all necessary values:
		b1ValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"b1",t)
		b2ValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"b2",t)
		
		b1TildeValuesAtCentroids = 
			b1ValuesAtCentroids.*(meshData.elementsNablaFTInv11s) + 
			b2ValuesAtCentroids.*(meshData.elementsNablaFTInv12s)
		
		b2TildeValuesAtCentroids = 
			b1ValuesAtCentroids.*(meshData.elementsNablaFTInv21s) + 
			b2ValuesAtCentroids.*(meshData.elementsNablaFTInv22s)
		
		MbIs = Array{Int64,2}(9,size(meshData.elements,2))
		MbJs = Array{Int64,2}(9,size(meshData.elements,2))
		MbVs = Array{Float64,2}(9,size(meshData.elements,2))
		for i = 1:1:3
			for j = 1:1:3
				MbIs[3*(i-1)+j,:] = meshData.elements[i,:]
				MbJs[3*(i-1)+j,:] = meshData.elements[j,:]
				MbVs[3*(i-1)+j,:] = abs(meshData.elementsDetNablaFTs).*(
					b1TildeValuesAtCentroids*(phiHats.PhiHatX0[j,i]) + 
					b2TildeValuesAtCentroids*(phiHats.PhiHatY0[j,i]))
			end
		end
		
		# Assemble sparse matrix:
		# If an index tuple [i,j] occurs multiple times, then sparse() sums up the corresponding values.
		Mb = sparse(collect(MbIs),collect(MbJs),collect(MbVs),size(meshData.nodes,2),size(meshData.nodes,2))
		
	end
	
	return Mb
end

function lseDataAssembleMc(doQuickAssemblyMc::Bool,equationData::EquationData,meshData::MeshData,phiHats::PhiHats,t::Float64)
	
	# Assembles the Matrix
	#
	# Mc := (\int_{\Omega} c*v_N*v_M dx)_{M,N \in \Nodes},
	#
	# where v_M \in \S^{1,1}(\Elements) is the piecewise linear hat function belonging to the node M.
	
	# We implement the following formula (1-point approximation of c):
	#
	# Mc = sum_{i,j=1}^{3} \sum_{T \in \Elements} abs(\det\nabla F_T)*
	#		[c(T)*(\hat{\Phi}_{00})_{ji}]*e_{N_i(T)}*e_{N_j(T)}^T,
	#
	# where
	#
	# \hat{\Phi}_{...} \in \R^{3x3} ... integrals over the reference basis functions on the reference element,
	# e_{N_i(T)},e_{N_j(T)} \in \R^{#nodes} ... euclidean unit vectors,
	#
	# C_T \in \R^2 ... element centroid,
	# c(C_T) \in \R ... the value of the function c at the point C_T,
	# Fij(T) := ((\nabla F_T)^{-1})_{ij} \in \R ... the entries of the inverse of the Jacobian matrix,
	# 		of the reference-element-mapping F_T : \hat{T} --> T.
	
	if doQuickAssemblyMc
		
		Mc = spzeros(size(meshData.nodes,2),size(meshData.nodes,2))
		
	else
		
		# Compute all necessary values:
		cValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"c",t)
		
		McIs = Array{Int64,2}(9,size(meshData.elements,2))
		McJs = Array{Int64,2}(9,size(meshData.elements,2))
		McVs = Array{Float64,2}(9,size(meshData.elements,2))
		for i = 1:1:3
			for j = 1:1:3
				McIs[3*(i-1)+j,:] = meshData.elements[i,:]
				McJs[3*(i-1)+j,:] = meshData.elements[j,:]
				McVs[3*(i-1)+j,:] = abs(meshData.elementsDetNablaFTs).*(
					cValuesAtCentroids*(phiHats.PhiHat00[j,i]))
			end
		end
		
		# Assemble sparse matrix:
		# If an index tuple [i,j] occurs multiple times, then sparse() sums up the corresponding values.
		Mc = sparse(collect(McIs),collect(McJs),collect(McVs),size(meshData.nodes,2),size(meshData.nodes,2))
		
	end
	
	return Mc
end

function lseDataAssembleMAbc(MA::SparseMatrixCSC{Float64,Int64},Mb::SparseMatrixCSC{Float64,Int64},Mc::SparseMatrixCSC{Float64,Int64})
	
	MAbc = MA + Mb + Mc
	
	return MAbc
end

function lseDataAssembleM(MAbc::SparseMatrixCSC{Float64,Int64})
	
	M = copy(MAbc)
	
	return M
end

################################################

function lseDataAssembleRF(doQuickAssemblyRF::Bool,equationData::EquationData,meshData::MeshData,phiHats::PhiHats,t::Float64)
	
	# Assembles the vector
	#
	# rF := (\int_{\Omega} f*v_M dx)_{M \in \Nodes},
	#
	# where v_M \in \S^{1,1}(\Elements) is the piecewise linear hat function belonging to the node M.
	
	# We implement the following formula (1-point approximation of f):
	#
	# rF = sum_{i=1}^{3} \sum_{T \in \Elements} abs(\det\nabla F_T)*
	#		[f(C_T)*(\hat{\Phi}_0)_i]*e_{N_i(T)},
	#
	# where
	#
	# \hat{\Phi}_{...} \in \R^3 ... integrals over the reference basis functions on the reference element,
	# e_{N_i(T)} \in \R^{#nodes} ... euclidean unit vectors,
	#
	# C_T \in \R^2 ... element centroid,
	# f(C_T) \in \R ... the value of the function f at the point C_T.
	
	if doQuickAssemblyRF
		
		rF = zeros(Float64,size(meshData.nodes,2))
		
	else
		
		# Compute all necessary values:
		fValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"f",t)
		
		rFIs = Array{Int64,2}(3,size(meshData.elements,2))
		rFVs = Array{Float64,2}(3,size(meshData.elements,2))
		for i = 1:1:3
			rFIs[i,:] = meshData.elements[i,:]
			rFVs[i,:] = abs(meshData.elementsDetNablaFTs).*fValuesAtCentroids*(phiHats.PhiHat0[i])
		end
		
		# Assemble sparse vector:
		# If an index i occurs multiple times, then sparsevec() sums up the corresponding values.
		rF = collect(full(sparsevec(collect(rFIs),collect(rFVs),size(meshData.nodes,2))))
		
	end
	
	return rF
end

function lseDataAssembleRFDiv(doQuickAssemblyRFDiv::Bool,equationData::EquationData,meshData::MeshData,phiHats::PhiHats,t::Float64)
	
	# Assembles the vector
	#
	# rFDiv := (\int_{\Omega} fDiv^T*v_M dx)_{M \in \Nodes},
	#
	# where v_M \in \S^{1,1}(\Elements) is the piecewise linear hat function belonging to the node M.
	
	# We implement the following formula (1-point approximation of fDiv):
	#
	# rFDiv = sum_{i=1}^{3} \sum_{T \in \Elements} abs(\det\nabla F_T)*
	#		[fDiv1Tilde(T)*(\hat{\Phi}_x)_i + fDiv2Tilde(T)*(\hat{\Phi}_y)_i]*
	#		e_{N_i(T)},
	#
	# where
	#
	# \hat{\Phi}_{...} \in \R^{3x3} ... integrals over the reference basis functions on the reference element,
	# e_{N_i(T)} \in \R^{#nodes} ... euclidean unit vectors,
	#
	# C_T \in \R^2 ... element centroid,
	# fDivi(C_T) \in \R ... the value of the function fDiv_i at the point C_T,
	# Fij(T) := ((\nabla F_T)^{-1})_{ij} \in \R ... the entries of the inverse of the Jacobian matrix,
	# 		of the reference-element-mapping F_T : \hat{T} --> T,
	#
	# fDiv1Tilde(T) := fDiv1(C_T)*F11(T)^2 + fDiv2(C_T)*F12(T),
	# fDiv2Tilde(T) := fDiv1(C_T)*F21(T)^2 + fDiv2(C_T)*F22(T).
	
	if doQuickAssemblyRFDiv
		
		rFDiv = zeros(Float64,size(meshData.nodes,2))
		
	else
		
		# Compute all necessary values:
		fDiv1ValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"fDiv1",t)
		fDiv2ValuesAtCentroids = computeFunctionValuesAtElementsCentroids(equationData,meshData,"fDiv2",t)
		
		fDiv1TildeValuesAtCentroids = 
			fDiv1ValuesAtCentroids.*(meshData.elementsNablaFTInv11s) + 
			fDiv2ValuesAtCentroids.*(meshData.elementsNablaFTInv12s)
		
		fDiv2TildeValuesAtCentroids = 
			fDiv1ValuesAtCentroids.*(meshData.elementsNablaFTInv21s) + 
			fDiv2ValuesAtCentroids.*(meshData.elementsNablaFTInv22s)
		
		rFDivIs = Array{Int64,2}(3,size(meshData.elements,2))
		rFDivVs = Array{Float64,2}(3,size(meshData.elements,2))
		for i = 1:1:3
			rFDivIs[i,:] = meshData.elements[i,:]
			rFDivVs[i,:] = abs(meshData.elementsDetNablaFTs).*(
				fDiv1TildeValuesAtCentroids*(phiHats.PhiHatX[i]) + 
				fDiv2TildeValuesAtCentroids*(phiHats.PhiHatY[i]))
		end
		
		# Assemble sparse vector:
		# If an index i occurs multiple times, then sparsevec() sums up the corresponding values.
		rFDiv = collect(full(sparsevec(collect(rFDivIs),collect(rFDivVs),size(meshData.nodes,2))))
		
	end
	
	return rFDiv
end

function lseDataAssembleRNeumann(doQuickAssemblyRNeumann::Bool,equationData::EquationData,meshData::MeshData,t::Float64)
	
	# Assembles the vector
	#
	# rNeumann := (\int_{\Gamma_Neumann} gNeumann*v_M dx)_{M \in \Nodes},
	#
	# where v_M \in \S^{1,1}(\Elements) is the piecewise linear hat function belonging to the node M.
	
	# We implement the following formula (1-point approximation of Neumann data):
	#
	# rNeumann = sum_{E \in neumannEdges} 1/2*volume(E)*gNeumann(C(E))*sum_{i=1}^{2} e_{N_i(E)},
	#
	# where
	#
	# gNeumann: \Gamma_Neumann --> \R ... Neumann data (encoded in equationData.domain.gammas[...].boundaryDataFunction),
	# C(E) \in \R^2 ... edge centroid,
	# volume(E) \in \R ... surface volume of the edge,
	# e_{N_i(E)} \in \R^{#nodes} ... euclidean unit vector.
	
	if doQuickAssemblyRNeumann
		
		rNeumann = zeros(Float64,size(meshData.nodes,2))
		
	else
		
		# Compute all necessary values:
		gNeumannValuesAtCentroids = computeFunctionValuesAtNeumannEdgesCentroids(equationData,meshData,t)
		
		rNeumannIs = Array{Int64,2}(2,size(meshData.edges,2))
		rNeumannVs = Array{Float64,2}(2,size(meshData.edges,2))
		for i = 1:1:2
			rNeumannIs[i,:] = meshData.edges[i,:]
			rNeumannVs[i,:] = 1/2*(meshData.edgesVolumes).*gNeumannValuesAtCentroids
		end
		
		# Assemble sparse vector:
		# If an index i occurs multiple times, then sparsevec() sums up the corresponding values.
		rNeumann = collect(full(sparsevec(collect(rNeumannIs),collect(rNeumannVs),size(meshData.nodes,2))))
		
	end
	
	return rNeumann
end

function lseDataAssembleDirichletNodesIdsAndGDirichletValuesAtDirichletNodes(equationData::EquationData,meshData::MeshData,t::Float64)
	
	dirichletNodesIds,gDirichletValuesAtDirichletNodes = computeFunctionValuesAtDirichletNodes(equationData,meshData,t)
	
	return dirichletNodesIds,gDirichletValuesAtDirichletNodes
end

function lseDataAssembleRDirichlet(meshData::MeshData,MAbc::SparseMatrixCSC{Float64,Int64},dirichletNodesIds::Array{Int64,1},gDirichletValuesAtDirichletNodes::Array{Float64,1})
	
	# Assembles the vector
	#
	# rDirichlet := (\int_{\Omega} \nabla v_M^T*A*\nabla gDirichlet + (b^T*\nabla gDirichlet)*v_M + c*gDirichlet*v_M dx)_{M \in \Nodes},
	#
	# where v_M \in \S^{1,1}(\Elements) is the piecewise linear hat function belonging to the node M.
	
	# We implement the following formula (exact Dirichlet data, not an approximation):
	#
	# rDirichlet = MAbc*sum_{N \in dirichletNodes} gDirichlet(N)*e_N,
	#
	# where
	#
	# MAbc \in \R^{#nodes x #nodes} ... the already computed system matrix.
	
	# If an index i occurs multiple times, then sparsevec() sums up the corresponding values.
	rDirichlet = MAbc*collect(full(sparsevec(dirichletNodesIds,gDirichletValuesAtDirichletNodes,size(meshData.nodes,2))))
	
	return rDirichlet
end

function lseDataAssembleR(rF::Array{Float64,1},rFDiv::Array{Float64,1},rNeumann::Array{Float64,1},rDirichlet::Array{Float64,1})
	
	r = rF + rFDiv + rNeumann
	
	if !isempty(rDirichlet)
		r += -rDirichlet
	end
	
	return r
end

################################################

function enforceDirichletBoundaryConditionsInM(M::SparseMatrixCSC{Float64,Int64},dirichletNodesIds::Array{Int64,1})
	
	# Enforce the dirichlet boundary conditions in M:
	
	# Delete rows and columns and set diagonal entries to 1.
	M[dirichletNodesIds,:] = 0
	M[:,dirichletNodesIds] = 0
	M[sub2ind(size(M),dirichletNodesIds,dirichletNodesIds)] = 1
	
	return nothing
end

function enforceDirichletBoundaryConditionsInR(r::Array{Float64,1},dirichletNodesIds::Array{Int64,1},gDirichletValuesAtDirichletNodes::Array{Float64,1})
	
	# Enforce the dirichlet boundary conditions in r:
	
	r[dirichletNodesIds] = gDirichletValuesAtDirichletNodes
	
	return nothing
end

function enforcePeriodicBoundaryConditionsInM(M::SparseMatrixCSC{Float64,Int64},meshData::MeshData)
	
	# Enforce the periodic boundary conditions in M:
	
	if !isempty(meshData.periodicNodes)
		
		mastersIds = collect(meshData.periodicNodes[1,:])
		slavesIds = collect(meshData.periodicNodes[2,:])
		slaveIdToMasterIdMapping = collect(full(sparsevec(slavesIds,mastersIds,size(meshData.nodes,2))))
		
		MIs,MJs = ind2sub(size(M),find(M))
		MIsMapped = slaveIdToMasterIdMapping[MIs]
		MIs[(MIsMapped.!=0)] = MIsMapped[(MIsMapped.!=0)]
		MVs = nonzeros(M)
		M[:,:] = sparse(MIs,MJs,MVs,size(M,1),size(M,2))
		
		M[sub2ind(size(M),slavesIds,mastersIds)] = 1
		M[sub2ind(size(M),slavesIds,slavesIds)] = -1
	end
	
	return nothing
end

function enforcePeriodicBoundaryConditionsInR(r::Array{Float64,1},meshData::MeshData)
	
	# Enforce the periodic boundary conditions in r:
	
	if !isempty(meshData.periodicNodes)
		
		mastersIds = collect(meshData.periodicNodes[1,:])
		slavesIds = collect(meshData.periodicNodes[2,:])
		slaveIdToMasterIdMapping = collect(full(sparsevec(slavesIds,mastersIds,size(meshData.nodes,2))))
		
		rIs = collect(1:1:size(meshData.nodes,2))
		rIsMapped = slaveIdToMasterIdMapping[rIs]
		rIs[(rIsMapped.!=0)] = rIsMapped[(rIsMapped.!=0)]
		rVs = r
		r[:] = collect(full(sparsevec(rIs,rVs,size(meshData.nodes,2))))
	end
	
	return nothing
end

################################################

function lseDataAssembleMMass(meshData::MeshData,phiHats::PhiHats)
	
	# Assembles the Matrix
	#
	# MMass := (\int_{\Omega} v_N*v_M dx)_{M,N \in \Nodes},
	#
	# where v_M \in \S^{1,1}(\Elements) is the piecewise linear hat function belonging to the node M.
	
	# We implement the following formula (exact values, not an approximation):
	#
	# MMass = sum_{i,j=1}^{3} \sum_{T \in \Elements} abs(\det\nabla F_T)*
	#		[(\hat{\Phi}_{00})_{ji}]*e_{N_i(T)}*e_{N_j(T)}^T,
	#
	# where
	#
	# \hat{\Phi}_{...} \in \R^{3x3} ... integrals over the reference basis functions on the reference element,
	# e_{N_i(T)},e_{N_j(T)} \in \R^{#nodes} ... euclidean unit vectors.
	
	# Compute all necessary values:
	MMassIs = Array{Int64,2}(9,size(meshData.elements,2))
	MMassJs = Array{Int64,2}(9,size(meshData.elements,2))
	MMassVs = Array{Float64,2}(9,size(meshData.elements,2))
	for i = 1:1:3
		for j = 1:1:3
			MMassIs[3*(i-1)+j,:] = meshData.elements[i,:]
			MMassJs[3*(i-1)+j,:] = meshData.elements[j,:]
			MMassVs[3*(i-1)+j,:] = abs(meshData.elementsDetNablaFTs)*(phiHats.PhiHat00[j,i])
		end
	end
	
	# Assemble sparse matrix:
	# If an index tuple [i,j] occurs multiple times, then sparse() sums up the corresponding values.
	MMass = sparse(collect(MMassIs),collect(MMassJs),collect(MMassVs),size(meshData.nodes,2),size(meshData.nodes,2))
	
	return MMass
end

function lseDataAssembleIntegralValues(meshData::MeshData,phiHats::PhiHats)
	
	# Assembles the vector
	#
	# integralValues := (\int_{\Omega} v_M dx)_{M \in \Nodes},
	#
	# where v_M \in \S^{1,1}(\Elements) is the piecewise linear hat function belonging to the node M.
	
	# We implement the following formula (exact values, not an approximation):
	#
	# integralValues = sum_{i=1}^{3} \sum_{T \in \Elements} abs(\det\nabla F_T)*[(\hat{\Phi}_0)_i]*e_{N_i(T)},
	#
	# where
	#
	# \hat{\Phi}_{...} \in \R^3 ... integrals over the reference basis functions on the reference element,
	# e_{N_i(T)} \in \R^{#nodes} ... euclidean unit vectors.
	
	# Compute all necessary values:
	integralValuesIs = Array{Int64,2}(3,size(meshData.elements,2))
	integralValuesVs = Array{Float64,2}(3,size(meshData.elements,2))
	for i = 1:1:3
		integralValuesIs[i,:] = meshData.elements[i,:]
		integralValuesVs[i,:] = abs(meshData.elementsDetNablaFTs)*(phiHats.PhiHat0[i])
	end
	
	# Assemble sparse vector:
	# If an index i occurs multiple times, then sparsevec() sums up the corresponding values.
	integralValues = collect(full(sparsevec(collect(integralValuesIs),collect(integralValuesVs),size(meshData.nodes,2))))
	
	return integralValues
end

################################################

function computeFunctionValuesAtElementsCentroids(equationData::EquationData,meshData::MeshData,coefficientChooser::ASCIIString,t::Float64)
	
	functionValuesAtCentroids = zeros(Float64,size(meshData.elements,2))
	
	for domain in vcat(equationData.domain,equationData.subdomains)
		
		if coefficientChooser == "uInitial"
			doRun = !domain.isZeroUInitial
			coefficientFunction = domain.uInitial
		elseif coefficientChooser == "uTInitial"
			doRun = !domain.isZeroUTInitial
			coefficientFunction = domain.uTInitial
		elseif coefficientChooser == "A11"
			doRun = true
			coefficientFunction = domain.A[1,1]
		elseif coefficientChooser == "A21"
			doRun = true
			coefficientFunction = domain.A[2,1]
		elseif coefficientChooser == "A12"
			doRun = true
			coefficientFunction = domain.A[1,2]
		elseif coefficientChooser == "A22"
			doRun = true
			coefficientFunction = domain.A[2,2]
		elseif coefficientChooser == "b1"
			doRun = !domain.isZeroB
			coefficientFunction = domain.b[1]
		elseif coefficientChooser == "b2"
			doRun = !domain.isZeroB
			coefficientFunction = domain.b[2]
		elseif coefficientChooser == "c"
			doRun = !domain.isZeroC
			coefficientFunction = domain.c
		elseif coefficientChooser == "f"
			doRun = !domain.isZeroF
			coefficientFunction = domain.f
		elseif coefficientChooser == "fDiv1"
			doRun = !domain.isZeroFDiv
			coefficientFunction = domain.fDiv[1]
		elseif coefficientChooser == "fDiv2"
			doRun = !domain.isZeroFDiv
			coefficientFunction = domain.fDiv[2]
		end
		
		if doRun
			
			insertionIndices = (meshData.elementsParentDomainsTypes.==domain.domainType)&(meshData.elementsParentDomainsIds.==domain.id)
			
			if (equationData.equationType=="elliptic")|in(coefficientChooser,["uInitial","uTInitial"])
				functionValuesAtCentroids[insertionIndices] += map(coefficientFunction,
					collect(meshData.elementsCentroids[1,insertionIndices]),
					collect(meshData.elementsCentroids[2,insertionIndices]))
			elseif in(equationData.equationType,["parabolic","hyperbolic"])
				functionValuesAtCentroids[insertionIndices] += map(coefficientFunction,
					collect(meshData.elementsCentroids[1,insertionIndices]),
					collect(meshData.elementsCentroids[2,insertionIndices]),
					t*ones(Float64,sum(insertionIndices)))
			end
		end
		
	end
	
	return functionValuesAtCentroids
end

function computeFunctionValuesAtNeumannEdgesCentroids(equationData::EquationData,meshData::MeshData,t::Float64)
	
	gNeumannValuesAtCentroids = zeros(Float64,size(meshData.edges,2))
	
	for domain in vcat(equationData.domain,equationData.holes)
		for gamma in domain.gammas
			if (gamma.boundaryDataType == "Neumann")&(!gamma.isZeroBoundaryDataFunction)
				
				insertionIndices =
					(meshData.edgesParentDomainsTypes.==domain.domainType)&
					(meshData.edgesParentDomainsIds.==domain.id)&
					(meshData.edgesParentGammasIds.==gamma.id)
				
				if equationData.equationType == "elliptic"
					gNeumannValuesAtCentroids[insertionIndices] += map(
						gamma.boundaryDataFunction,
						collect(meshData.edgesParentGammasParameterRanges[2,insertionIndices]))
				elseif in(equationData.equationType,["parabolic","hyperbolic"])
					gNeumannValuesAtCentroids[insertionIndices] += map(
						gamma.boundaryDataFunction,
						collect(meshData.edgesParentGammasParameterRanges[2,insertionIndices]),
						t*ones(Float64,sum(insertionIndices)))
				end
				
			end
		end
	end
	
	return gNeumannValuesAtCentroids
end

function computeFunctionValuesAtDirichletNodes(equationData::EquationData,meshData::MeshData,t::Float64)
	
	dirichletEdgesIds = falses(size(meshData.edges,2))
	gDirichletValuesAtDirichletNode1s = zeros(Float64,size(meshData.edges,2))
	gDirichletValuesAtDirichletNode2s = zeros(Float64,size(meshData.edges,2))
	
	for domain in vcat(equationData.domain,equationData.holes)
		for gamma in domain.gammas
			if gamma.boundaryDataType == "Dirichlet"
				
				insertionIndices =
					(meshData.edgesParentDomainsTypes.==domain.domainType)&
					(meshData.edgesParentDomainsIds.==domain.id)&
					(meshData.edgesParentGammasIds.==gamma.id)
				
				dirichletEdgesIds[insertionIndices] = true
				
				if !gamma.isZeroBoundaryDataFunction
					if equationData.equationType == "elliptic"
						
						gDirichletValuesAtDirichletNode1s[insertionIndices] += map(
							gamma.boundaryDataFunction,
							collect(meshData.edgesParentGammasParameterRanges[1,insertionIndices]))
						gDirichletValuesAtDirichletNode2s[insertionIndices] += map(
							gamma.boundaryDataFunction,
							collect(meshData.edgesParentGammasParameterRanges[3,insertionIndices]))
						
					elseif in(equationData.equationType,["parabolic","hyperbolic"])
						
						gDirichletValuesAtDirichletNode1s[insertionIndices] += map(
							gamma.boundaryDataFunction,
							collect(meshData.edgesParentGammasParameterRanges[1,insertionIndices]),
							t*ones(Float64,sum(insertionIndices)))
						gDirichletValuesAtDirichletNode2s[insertionIndices] += map(
							gamma.boundaryDataFunction,
							collect(meshData.edgesParentGammasParameterRanges[3,insertionIndices]),
							t*ones(Float64,sum(insertionIndices)))
						
					end
				end
				
			end
		end
	end
	
	dirichletNodesIds_nonUnique = vcat(
		collect(meshData.edges[1,dirichletEdgesIds]),
		collect(meshData.edges[2,dirichletEdgesIds]))
	gDirichletValuesAtDirichletNodes_nonUnique = vcat(
		gDirichletValuesAtDirichletNode1s[dirichletEdgesIds],
		gDirichletValuesAtDirichletNode2s[dirichletEdgesIds])
	
	uniquePerm = myUniqueperm(dirichletNodesIds_nonUnique)
	
	dirichletNodesIds = dirichletNodesIds_nonUnique[uniquePerm]
	gDirichletValuesAtDirichletNodes = gDirichletValuesAtDirichletNodes_nonUnique[uniquePerm]
	
	return dirichletNodesIds,gDirichletValuesAtDirichletNodes
end

################################################

function valuesAtNodesToValuesAtElementsCentroids(meshData::MeshData,valuesAtNodes::Array{Float64,1},phiHats::PhiHats)
	
	# Implements the action of the L^2(\Omega)-orthogonal-projection \Pi : L^2(\Omega) --> \S^{0,0}(\Elements)
	# on functions \S^{1,1}(\Elements).
	
	valuesAtCentroids = zeros(Float64,size(meshData.elements,2))
	for i = 1:1:3
		valuesAtCentroids += 2*(phiHats.PhiHat0[i])*valuesAtNodes[collect(meshData.elements[i,:])]
	end
	
	return valuesAtCentroids
end

function valuesAtElementsCentroidsToValuesAtNodes(meshData::MeshData,valuesAtElementsCentroids::Array{Float64,1})
	
	# Implements the action of the Clement-quasi-interpolation-operator \Pi : L^2(\Omega) --> \S^{1,1}(\Elements)
	# on functions \S^{0,0}(\Elements).
	
	nodesPatchesVolumesIs = Array{Int64,2}(3,size(meshData.elements,2))
	nodesPatchesVolumesVs = Array{Float64,2}(3,size(meshData.elements,2))
	for i = 1:1:3
		nodesPatchesVolumesIs[i,:] = meshData.elements[i,:]
		nodesPatchesVolumesVs[i,:] = meshData.elementsVolumes
	end
	nodesPatchesVolumes = collect(full(sparsevec(collect(nodesPatchesVolumesIs),collect(nodesPatchesVolumesVs),size(meshData.nodes,2))))
	
	valuesAtNodesIs = Array{Int64,2}(3,size(meshData.elements,2))
	valuesAtNodesVs = Array{Float64,2}(3,size(meshData.elements,2))
	for i = 1:1:3
		valuesAtNodesIs[i,:] = meshData.elements[i,:]
		valuesAtNodesVs[i,:] = (meshData.elementsVolumes)./(nodesPatchesVolumes[collect(meshData.elements[i,:])]).*(valuesAtElementsCentroids)
	end
	valuesAtNodes = collect(full(sparsevec(collect(valuesAtNodesIs),collect(valuesAtNodesVs),size(meshData.nodes,2))))
	
	return valuesAtNodes
end

################################################

function calculateL2Norm(meshData::MeshData,u::Array{Float64},phiHats::PhiHats)
	
	# We implement the following formula (exact values, not an approximation):
	#
	# ||u||_{L^2(\Omega)}^2 = sum_{i,j=1}^{3} \sum_{T \in \Elements} abs(\det\nabla F_T)*
	#		u(N_i(T))*u(N_j(T))*[(\hat{\Phi}_{00})_{ij}],
	#
	# where
	#
	# \hat{\Phi}_{...} \in \R^{3x3} ... integrals over the reference basis functions on the reference element.
	
	uL2NormSquared = 0
	for i = 1:1:3
		for j = 1:1:3
			uL2NormSquared += sum(abs(meshData.elementsDetNablaFTs).*u[collect(meshData.elements[i,:])].*u[collect(meshData.elements[j,:])]*(phiHats.PhiHat00[j,i]))
		end
	end
	uL2Norm = sqrt(uL2NormSquared)
	
	return uL2Norm
end

function calculateH1Seminorm(meshData::MeshData,u::Array{Float64},phiHats::PhiHats)
	
	# We implement the following formula (exact values, not an approximation):
	#
	# |u|_{H^1(\Omega)}^2 = sum_{i,j=1}^{3} \sum_{T \in \Elements} abs(\det\nabla F_T)*
	#		u(N_i(T))*u(N_j(T))*
	#		[A11Tilde(T)*(\hat{\Phi}_{xx})_{ij} + A12Tilde(T)*(\hat{\Phi}_{xy})_{ij} +
	#		A21Tilde(T)*(\hat{\Phi}_{yx})_{ij} + A22Tilde(T)*(\hat{\Phi}_{yy})_{ij}],
	#
	# where
	#
	# \hat{\Phi}_{...} \in \R^{3x3} ... integrals over the reference basis functions on the reference element,
	#
	# Fij(T) := ((\nabla F_T)^{-1})_{ij} \in \R ... the entries of the inverse of the Jacobian matrix,
	# 		of the reference-element-mapping F_T : \hat{T} --> T,
	#
	# A11Tilde(T) := F11(T)^2 + F12(T)^2,
	# A12Tilde(T) := F11(T)*F21(T) + F12(T)*F22(T),
	# A21Tilde(T) := F11(T)*F21(T) + F12(T)*F22(T),
	# A22Tilde(T) := F21(T)^2 + F22(T)^2.
	
	# Compute all necessary values:
	A11TildeValuesAtCentroids = (meshData.elementsNablaFTInv11s).^2 + (meshData.elementsNablaFTInv12s).^2
	A12TildeValuesAtCentroids = (meshData.elementsNablaFTInv11s).*(meshData.elementsNablaFTInv21s) + 
		(meshData.elementsNablaFTInv12s).*(meshData.elementsNablaFTInv22s)
	A21TildeValuesAtCentroids = A12TildeValuesAtCentroids
	A22TildeValuesAtCentroids = (meshData.elementsNablaFTInv21s).^2 + (meshData.elementsNablaFTInv22s).^2
	
	# Add all values:
	uH1SeminormSquared = 0
	for i = 1:1:3
		for j = 1:1:3
			uH1SeminormSquared += sum(abs(meshData.elementsDetNablaFTs).*
				u[collect(meshData.elements[i,:])].*u[collect(meshData.elements[j,:])].*
				(A11TildeValuesAtCentroids*(phiHats.PhiHatXX[i,j]) + 
				A12TildeValuesAtCentroids*(phiHats.PhiHatXY[i,j]) + 
				A21TildeValuesAtCentroids*(phiHats.PhiHatYX[i,j]) + 
				A22TildeValuesAtCentroids*(phiHats.PhiHatYY[i,j])))
		end
	end
	uH1Seminorm = sqrt(uH1SeminormSquared)
	
	return uH1Seminorm
end

function calculateIntegralValue(meshData::MeshData,u::Array{Float64},phiHats::PhiHats)
	
	integralValues = lseDataAssembleIntegralValues(meshData,phiHats)
	uIntegralValue = dot(u,integralValues)
	
	return uIntegralValue
end

################################################

