
function infoDataDisplay(infoData::InfoData)
	
	println("")
	println("Meshsize: "*string(round(infoData.meshsize,5)))
	println("Number of elements: "*string(infoData.numberOfElements))
	println("Number of nodes: "*string(infoData.numberOfNodes))
	println("System size: "*string(infoData.systemSize)*" x "*string(infoData.systemSize))
	println("Sparse system matrix storage: "*string(round(infoData.sparseSystemMatrixStorage,2))*" MBytes")
	println("(Estim. full stiffness matrix storage: "*string(round(infoData.estimatedFullSystemMatrixStorage,2))*" MBytes)")
	println("Used solver: "*infoData.solverType)
	println("")
	
	for processStopwatch in infoData.processesStopwatches
		if processStopwatch.isDisplayed
			println("ProcessStopwatch "*string(processStopwatch.id)*": "*
				processStopwatch.processName*": "*
				string(round(processStopwatch.processRuntime,2))*" s")
		end
	end
	
	return nothing
end
