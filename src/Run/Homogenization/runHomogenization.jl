
include("../../EquationData/equationData.jl")
include("../../GeoData/geoData.jl")
include("../../MeshData/meshData.jl")
include("../../EquationData to GeoData to MeshData/equationDataToGeoData.jl")
include("../../EquationData to GeoData to MeshData/geoDataToMeshData.jl")
include("../../LSEData/lseData.jl")
include("../../Miscellaneous/miscellaneous.jl")
include("../../InfoData/infoData.jl")

##############################################

infoDataProcessesStopwatches = Array{ProcessStopwatch,1}(0)
tic()

##############################################

println("Assemble equationData...")
tic()
equationData = equationDataAssemble("./homogenizationEquationData.jl")
push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Assemble equationData",toq(),true))

##############################################

println("Assemble geoData...")
tic()
geoData = equationDataToGeoData(equationData)
push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Assemble geoData",toq(),true))

##############################################

println("Plot geoData...")
tic()
# geoDataDisplay(geoData)
geoDataPlot(geoData)
push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Plot geoData",toq(),true))

##############################################

intendedMeshsize = 0.025

println("Assemble meshData...")
tic()
meshData = geoDataToMeshData(geoData,intendedMeshsize)
push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Assemble meshData",toq(),true))

##############################################

println("Plot meshData...")
tic()
# displayMeshData(meshData)
if intendedMeshsize >= 0.15
	meshDataPlotSlow(meshData)
else
	meshDataPlot(meshData)
end
push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Plot meshData",toq(),true))

##############################################
##############################################

domainFDiv1 = (x,y) -> -equationData.domain.A[1,1](x,y)
domainFDiv2 = (x,y) -> -equationData.domain.A[2,1](x,y)
equationData.domain.fDiv = [domainFDiv1,domainFDiv2]
for subdomain in equationData.subdomains
	subdomainsFDiv1 = (x,y) -> -subdomain.A[1,1](x,y)
	subdomainsFDiv2 = (x,y) -> -subdomain.A[2,1](x,y)
	subdomain.fDiv = [subdomainsFDiv1,subdomainsFDiv2]
end

##############################################

deltaT = 0.0

println("Assemble lseData...")
tic()
lseData = lseDataAssemble(equationData,meshData,deltaT)
push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Assemble lseData",toq(),true))

##############################################

println("Solve lseData...")
tic()
lseDataSolve(equationData,meshData,lseData)
push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Solve lseData",toq(),true))

##############################################

chi1 = copy(lseData.solutions[1])
r1 = copy(lseData.r)

##############################################

# println("Plot lseData...")
# tic()
# lseDataPlot(meshData,lseData)
# push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Plot lseData",toq(),true))

##############################################
##############################################

domainFDiv1 = (x,y) -> -equationData.domain.A[1,2](x,y)
domainFDiv2 = (x,y) -> -equationData.domain.A[2,2](x,y)
equationData.domain.fDiv = [domainFDiv1,domainFDiv2]
for subdomain in equationData.subdomains
	subdomainsFDiv1 = (x,y) -> -subdomain.A[1,2](x,y)
	subdomainsFDiv2 = (x,y) -> -subdomain.A[2,2](x,y)
	subdomain.fDiv = [subdomainsFDiv1,subdomainsFDiv2]
end

##############################################

deltaT = 0.0

println("Assemble lseData...")
tic()
lseData = lseDataAssemble(equationData,meshData,deltaT)
push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Assemble lseData",toq(),true))

##############################################

println("Solve lseData...")
tic()
lseDataSolve(equationData,meshData,lseData)
push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Solve lseData",toq(),true))

##############################################

chi2 = copy(lseData.solutions[1])
r2 = copy(lseData.r)

##############################################

println("Plot lseData...")
tic()
lseDataPlot(meshData,lseData)
push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Plot lseData",toq(),true))

##############################################
##############################################

# Compute effective A (only valid for symmetric A):
A11ValuesAtNodes = valuesAtElementsCentroidsToValuesAtNodes(meshData,
	computeFunctionValuesAtElementsCentroids(equationData,meshData,"A11",0.0))
A21ValuesAtNodes = valuesAtElementsCentroidsToValuesAtNodes(meshData,
	computeFunctionValuesAtElementsCentroids(equationData,meshData,"A21",0.0))
A12ValuesAtNodes = valuesAtElementsCentroidsToValuesAtNodes(meshData,
	computeFunctionValuesAtElementsCentroids(equationData,meshData,"A12",0.0))
A22ValuesAtNodes = valuesAtElementsCentroidsToValuesAtNodes(meshData,
	computeFunctionValuesAtElementsCentroids(equationData,meshData,"A22",0.0))

A11Effective = calculateIntegralValue(meshData,A11ValuesAtNodes,lseData.phiHats) + dot(chi1,-r1)
A21Effective = calculateIntegralValue(meshData,A21ValuesAtNodes,lseData.phiHats) + dot(chi1,-r2)
A12Effective = calculateIntegralValue(meshData,A12ValuesAtNodes,lseData.phiHats) + dot(chi2,-r1)
A22Effective = calculateIntegralValue(meshData,A22ValuesAtNodes,lseData.phiHats) + dot(chi2,-r2)
AEffective = hcat([A11Effective,A21Effective],[A12Effective,A22Effective])

println("")
println("AEffective = ")
display(AEffective)

##############################################

push!(infoDataProcessesStopwatches,ProcessStopwatch(-1,"Total",toq(),true))

##############################################

infoData = infoDataAssemble(infoDataProcessesStopwatches,meshData,lseData)
infoDataDisplay(infoData)



