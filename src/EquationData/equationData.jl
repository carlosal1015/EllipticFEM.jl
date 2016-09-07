
# We are solving the elliptic equation
#
# -div(A(x,y)*nabla u) + b(x,y)^T*nabla u + c(x,y)*u = f(x,y) - div(fDiv(x,y)) in \Omega,
# u = gDirichlet on \Gamma_Dirichlet,
# (A(x,y)*\nabla u)^T*n = gNeumann on \Gamma_Neumann,
# u periodic on \Gamma_periodic,
#
# respectively the parabolic equation
#
# \partial_t u - div(A(x,y,t)*nabla u) + b(x,y,t)^T*nabla u + c(x,y,t)*u = f(x,y,t) - div(fDiv(x,y,t)) in \Omega,
# u = gDirichlet on \Gamma_Dirichlet,
# (A(x,y,t)*\nabla u)^T*n = gNeumann on \Gamma_Neumann,
# u periodic on \Gamma_periodic,
# u(t=0) = uInitial
#
# respectively the hyperbolic equation
#
# \partial_{tt} u - div(A(x,y,t)*nabla u) + b(x,y,t)^T*nabla u + c(x,y,t)*u = f(x,y,t) - div(fDiv(x,y,t)) in \Omega,
# u = gDirichlet on \Gamma_Dirichlet,
# (A(x,y,t)*\nabla u)^T*n = gNeumann on \Gamma_Neumann,
# u periodic on \Gamma_periodic,
# u(t=0) = uInitial,
# (\partial_t u)(t=0) = uTInitial,
#
# where the coefficients A,b,c,f,fDiv,uInitial,uTInitial may depend on space (x,y) and time t and can be defined for each subdomain of \Omega individually.
# The boundary data gDirichlet and gNeumann may also depend on space and time and can be defined for each part of the boundary \partial\Omega individually.
# The dimensions are given as follows:
#
# A(x,y,t) = hcat([A11,A21],[A12,A22]) \in \R^{2x2},
# b(x,y,t) = [b1,b2] \in \R^2,
# c(x,y,t) \in \R,
# f(x,y,t) \in \R,
# fDiv(x,y,t) = [fDiv1,fDiv2] \in \R^2 (assuming fDiv = 0 on \Gamma_Neumann and fDiv periodic on \Gamma_periodic),
# gDirichlet(x,y,t) \in \R,
# gNeumann(x,y,t) \in \R,
# uInitial(x,y) \in \R,
# uTInitial(x,y) \in \R.

type GammaNode
	id::Int64
	
	coords::Array{Float64}
end

type Gamma
	id::Int64
	
	parametrizationType::ASCIIString	# "gammaLine", "gammaArc", "gammaParabolaArc"
	parametrizationGammaNodesIds::Array{Int64}
	parametrizationCoefficients::Array{Float64}		# e.g.: angle for "gammaArc"
	parametrizationNumberOfPartitions::Int64
	parametrization::Function
	
	boundaryDataType::ASCIIString		# "Dirichlet", "Neumann", "periodic"
	boundaryDataFunction::Function
	
	isZeroBoundaryDataFunction::Bool
	isConstantInTimeBoundaryDataFunction::Bool
end

type PeriodicGammasPair
	id::Int64
	gammasIds::Array{Int64}
	useSameOrientation::Bool
end

type Domain
	id::Int64
	domainType::ASCIIString		# "domain", "hole", "subdomain"
	gammaNodes::Array{GammaNode}
	gammas::Array{Gamma}
	periodicGammasPairs::Array{PeriodicGammasPair}
	
	uInitial::Function
	uTInitial::Function
	
	A::Array{Function,2}
	b::Array{Function}
	c::Function
	f::Function
	fDiv::Array{Function}
	
	isZeroUInitial::Bool
	isZeroUTInitial::Bool
	
	isIdentityA::Bool
	isZeroB::Bool
	isZeroC::Bool
	isZeroF::Bool
	isZeroFDiv::Bool
	
	isConstantInTimeA::Bool
	isConstantInTimeB::Bool
	isConstantInTimeC::Bool
	isConstantInTimeF::Bool
	isConstantInTimeFDiv::Bool
end

type EquationData
	equationType::ASCIIString		# "elliptic", "parabolic", "hyperbolic"
	tMax::Float64
	
	domain::Domain
	holes::Array{Domain}
	subdomains::Array{Domain}
end

include("./equationDataAssemble.jl")
include("./equationDataGammasInitializeParametrizations.jl")
include("./equationDataGammasRotateGammaNodes.jl")


