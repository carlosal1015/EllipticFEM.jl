
# Stores some basic information about runtimes and sizes of the LSE.

type ProcessStopwatch
	id::Int64
	processName::ASCIIString
	processRuntime::Float64
	isDisplayed::Bool
end

type InfoData
	
	meshsize::Float64
	numberOfElements::Int64
	numberOfNodes::Int64
	systemSize::Int64
	sparseSystemMatrixStorage::Float64
	estimatedFullSystemMatrixStorage::Float64
	solverType::ASCIIString
	
	processesStopwatches::Array{ProcessStopwatch}
end

include("./infoDataAssemble.jl")
include("./infoDataDisplay.jl")
