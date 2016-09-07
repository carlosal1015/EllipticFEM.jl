
type PhiHats
	# Defines the integrals over the reference basis functions \phi_i on the
	# 2D reference element \hat{T}.
	#
	# \hat{T} := {(x,y) \in \R^2 : x \in (0,1), y \in (0,1-x)}
	# \hat{\phi}_1(x,y) := 1-x-y
	# \hat{\phi}_2(x,y) := x
	# \hat{\phi}_3(x,y) := y
	#
	# \hat{\Phi}_0 := (\int_{\hat{T}} \phi_i dx)_{i=1}^{3} \in \R^3
	# \hat{\Phi}_x := (\int_{\hat{T}} \partial_x\phi_i dx)_{i=1}^{3} \in \R^3
	# \hat{\Phi}_y := (\int_{\hat{T}} \partial_y\phi_i dx)_{i=1}^{3} \in \R^3
	# \hat{\Phi}_{00} := (\int_{\hat{T}} \phi_i*\phi_j dx)_{i,j=1}^{3} \in \R^{3x3}
	# \hat{\Phi}_{x0} := (\int_{\hat{T}} (\partial_x \phi_i)*\phi_j dx)_{i,j=1}^{3} \in \R^{3x3}
	# \hat{\Phi}_{y0} := (\int_{\hat{T}} (\partial_y \phi_i)*\phi_j dx)_{i,j=1}^{3} \in \R^{3x3}
	# \hat{\Phi}_{0x} := (\int_{\hat{T}} \phi_i*(\partial_x \phi_j) dx)_{i,j=1}^{3} \in \R^{3x3}
	# \hat{\Phi}_{0y} := (\int_{\hat{T}} \phi_i*(\partial_y \phi_j) dx)_{i,j=1}^{3} \in \R^{3x3}
	# \hat{\Phi}_{xx} := (\int_{\hat{T}} (\partial_x \phi_i)*(\partial_x \phi_j) dx)_{i,j=1}^{3} \in \R^{3x3}
	# \hat{\Phi}_{xy} := (\int_{\hat{T}} (\partial_x \phi_i)*(\partial_y \phi_j) dx)_{i,j=1}^{3} \in \R^{3x3}
	# \hat{\Phi}_{yx} := (\int_{\hat{T}} (\partial_y \phi_i)*(\partial_x \phi_j) dx)_{i,j=1}^{3} \in \R^{3x3}
	# \hat{\Phi}_{yy} := (\int_{\hat{T}} (\partial_y \phi_i)*(\partial_y \phi_j) dx)_{i,j=1}^{3} \in \R^{3x3}
	
	PhiHat0::Array{Float64,1}
	PhiHatX::Array{Float64,1}
	PhiHatY::Array{Float64,1}
	PhiHat00::Array{Float64,2}
	PhiHatX0::Array{Float64,2}
	PhiHatY0::Array{Float64,2}
	PhiHat0X::Array{Float64,2}
	PhiHat0Y::Array{Float64,2}
	PhiHatXX::Array{Float64,2}
	PhiHatXY::Array{Float64,2}
	PhiHatYX::Array{Float64,2}
	PhiHatYY::Array{Float64,2}
end

type LSEData
	
	# Stores system matrices and right hand sides of the LSE.
	
	phiHats::PhiHats
	
	MA::SparseMatrixCSC{Float64,Int64}
	Mb::SparseMatrixCSC{Float64,Int64}
	Mc::SparseMatrixCSC{Float64,Int64}
	MAbc::SparseMatrixCSC{Float64,Int64}
	M::SparseMatrixCSC{Float64,Int64}
	
	rF::Array{Float64,1}
	rFDiv::Array{Float64,1}
	rNeumann::Array{Float64,1}
	dirichletNodesIds::Array{Int64,1}
	gDirichletValuesAtDirichletNodes::Array{Float64,1}
	rDirichlet::Array{Float64,1}
	r::Array{Float64,1}
	
	MMass::SparseMatrixCSC{Float64,Int64}
	integralValues::Array{Float64,1}
	
	solverType::ASCIIString
	
	deltaT::Float64
	solutionsTimestamps::Array{Float64,1}
	solutions::Array{Array{Float64,1},1}
end

include("./lseDataAssembleAndUpdate.jl")
include("./lseDataPlot.jl")
include("./lseDataSolve.jl")
