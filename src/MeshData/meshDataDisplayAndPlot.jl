
using PyPlot

function meshDataDisplay(meshData::MeshData)
	
	println("")
	println("---------------------------------------")
	println("meshData.nodes = ")
	display(meshData.nodes)
	println("")
	
	println("---------------------------------------")
	println("meshData.edges = ")
	display(meshData.edges)
	println("")
	
	println("meshData.edgesParentGeoEdgesIds' = ")
	display(meshData.edgesParentGeoEdgesIds')
	println("")
	
	println("meshData.edgesParentDomainsTypes' = ")
	display(meshData.edgesParentDomainsTypes')
	println("")
	
	println("meshData.edgesParentDomainsIds' = ")
	display(meshData.edgesParentDomainsIds')
	println("")
	
	println("meshData.edgesParentGammasIds' = ")
	display(meshData.edgesParentGammasIds')
	println("")
	
	println("meshData.edgesParentGammasParameterRanges = ")
	display(meshData.edgesParentGammasParameterRanges)
	println("")
	
	println("meshData.edgesVolumes' = ")
	display(meshData.edgesVolumes')
	println("")
	
	println("---------------------------------------")
	println("meshData.elements = ")
	display(meshData.elements)
	println("")
	
	println("meshData.elementsParentSubdomainsIntersectionsIds' = ")
	display(meshData.elementsParentSubdomainsIntersectionsIds')
	println("")
	
	println("meshData.elementsParentDomainsTypes' = ")
	display(meshData.elementsParentDomainsTypes')
	println("")
	
	println("meshData.elementsParentDomainsIds' = ")
	display(meshData.elementsParentDomainsIds')
	println("")
	
	println("meshData.elementsDetNablaFTs' = ")
	display(meshData.elementsDetNablaFTs')
	println("")
	
	println("meshData.elementsNablaFTInv11s' = ")
	display(meshData.elementsNablaFTInv11s')
	println("")
	
	println("meshData.elementsNablaFTInv21s' = ")
	display(meshData.elementsNablaFTInv21s')
	println("")
	
	println("meshData.elementsNablaFTInv12s' = ")
	display(meshData.elementsNablaFTInv12s')
	println("")
	
	println("meshData.elementsNablaFTInv22s' = ")
	display(meshData.elementsNablaFTInv22s')
	println("")
	
	println("meshData.elementsVolumes' = ")
	display(meshData.elementsVolumes')
	println("")
	
	println("meshData.elementsCentroids = ")
	display(meshData.elementsCentroids)
	println("")
	
	println("meshData.elementsMeshsizes' = ")
	display(meshData.elementsMeshsizes')
	println("")
	
	println("---------------------------------------")
	println("meshData.periodicNodes = ")
	display(meshData.periodicNodes)
	println("")
	
	println("---------------------------------------")
	println("")
	
	return nothing
end

function meshDataPlot(meshData::MeshData)
	
	doSaveFigure = true
	
	# Define nodes plot style:
	nodesVisible = false
	
	nodesColor = "red"
	nodesMarkerType = "o"
	nodesMarkerSize = 3
	
	# Define edges plot style:
	edgesVisible = true
	
	edgesColor = "black"
	edgesLineStyle = "-"
	edgesLineWidth = 2

	edgesColorSubdomains = "blue"
	edgesLineWidthSubdomains = 2

	# Define elements plot style:
	elementsVisible = true
	
	elementsColor = "blue"
	elementsLineStyle = "-"
	elementsLineWidth = 1

	# Plot:
	figure("MeshData")
	clf()
	
	if elementsVisible
		triplot(collect(meshData.nodes[1,:]),collect(meshData.nodes[2,:]),meshData.elements'-1,
			color=elementsColor,
			linewidth=elementsLineWidth,
			marker="",
			figure=figure("MeshData"))
	end
	
	if edgesVisible
		edgesIds = findin(meshData.edgesParentDomainsTypes,["domain","hole"])
		if !isempty(edgesIds)
			triplot(collect(meshData.nodes[1,:]),collect(meshData.nodes[2,:]),vcat(meshData.edges[1,edgesIds],meshData.edges[2,edgesIds],meshData.edges[2,edgesIds])'-1,
				color=edgesColor,
				linewidth=edgesLineWidth,
				marker="",
				figure=figure("MeshData"))
		end
		
		edgesIds = find(meshData.edgesParentDomainsTypes.=="subdomain")
		if !isempty(edgesIds)
			triplot(collect(meshData.nodes[1,:]),collect(meshData.nodes[2,:]),vcat(meshData.edges[1,edgesIds],meshData.edges[2,edgesIds],meshData.edges[2,edgesIds])'-1,
				color=edgesColorSubdomains,
				linewidth=edgesLineWidthSubdomains,
				marker="",
				figure=figure("MeshData"))
		end
	end
	
	if nodesVisible
		plot(collect(meshData.nodes[1,:]),collect(meshData.nodes[2,:]),
			color=nodesColor,
			markeredgecolor=nodesColor,
			linestyle="",
			marker=nodesMarkerType,
			markersize=nodesMarkerSize,
			figure=figure("MeshData"))
	end
	
	axis("equal")
	title("MeshData")
	xlabel("x")
	ylabel("y")
	
	xMin,xMax = xlim()
	xMean = (xMin+xMax)/2
	xlim(xMean + 1.5*(xMin-xMean),xMean + 1.5*(xMax-xMean))
	
	if doSaveFigure
		savefig("./Saved images/meshData.png")
	end
	
	return nothing
end

function meshDataPlotSlow(meshData::MeshData)
	
	doSaveFigure = true
	
	# Define nodes plot style:
	nodesVisible = true
	
	nodesColor = "red"
	nodesMarkerType = "o"
	nodesMarkerSize = 5
	
	nodesLabelsVisible = false
	nodesLabelsOffset = 0.005
	
	nodesLabelsContentId = true
	nodesLabelsContentCoords = false
	
	# Define edges plot style:
	edgesVisible = true
	
	edgesColor = "black"
	edgesLineStyle = "-"
	edgesLineWidth = 2
	
	edgesLabelsVisible = false
	edgesLabelsOffset = 0.005
	
	edgesLabelsContentId = true
	edgesLabelsContentNodesIds = false
	edgesLabelsContentParentGeoEdgesIds = false
	edgesLabelsContentParentDomainsTypes = false
	edgesLabelsContentParentDomainsIds = false
	edgesLabelsContentParentGammasIds = false
	edgesLabelsContentParentGammasParameterRanges = false
	edgesLabelsContentVolumes = false
	
	edgesColorSubdomains = "blue"
	edgesLineWidthSubdomains = 2
	edgesLabelsVisibleSubdomains = edgesLabelsVisible
	
	# Define elements plot style:
	elementsVisible = true
	
	elementsColor = "blue"
	elementsLineStyle = "-"
	elementsLineWidth = 1
	
	elementsLabelsVisible = false
	elementsLabelsOffset = -0.015
	
	elementsLabelsContentId = true
	elementsLabelsContentNodesIds = false
	elementsLabelsContentParentSubdomainsIntersectionsIds = false
	elementsLabelsContentParentDomainsTypes = false
	elementsLabelsContentParentDomainsIds = false
	elementsLabelsContentDetNablaFTs = false
	elementsLabelsContentNablaFTInv11s = false
	elementsLabelsContentNablaFTInv21s = false
	elementsLabelsContentNablaFTInv12s = false
	elementsLabelsContentNablaFTInv22s = false
	elementsLabelsContentVolumes = false
	elementsLabelsContentCentroids = false
	elementsLabelsContentMeshsizes = false
	
	# Plot:
	figure("MeshData")
	clf()
	
	hasBeenPlotted = falses(size(meshData.nodes,2),size(meshData.nodes,2))
	if edgesVisible
		for edgeId = 1:1:size(meshData.edges,2)
			if in(meshData.edgesParentDomainsTypes[edgeId],["domain","hole"])
				currentEdgeColor = edgesColor
				currentEdgeLineWidth = edgesLineWidth
				currentEdgeLabelVisible = edgesLabelsVisible
			else
				currentEdgeColor = edgesColorSubdomains
				currentEdgeLineWidth = edgesLineWidthSubdomains
				currentEdgeLabelVisible = edgesLabelsVisibleSubdomains
			end
			
			nodeId1 = meshData.edges[1,edgeId]
			nodeId2 = meshData.edges[2,edgeId]
			
			# Plot data:
			edgePositionXs = [meshData.nodes[1,nodeId1],meshData.nodes[1,nodeId2]]
			edgePositionYs = [meshData.nodes[2,nodeId1],meshData.nodes[2,nodeId2]]
			plot(edgePositionXs,edgePositionYs,
				color=currentEdgeColor,
				linestyle=edgesLineStyle,
				linewidth=currentEdgeLineWidth,
				figure=figure("MeshData"))
			
			hasBeenPlotted[nodeId1,nodeId2] = true
			hasBeenPlotted[nodeId2,nodeId1] = true
			
			if currentEdgeLabelVisible
				
				# Build label content:
				labelContent = ""
				if edgesLabelsContentId
					labelContent *= string(edgeId)*"|"
				end
				if edgesLabelsContentNodesIds
					labelContent *= string(collect(meshData.edges[:,edgeId]))*"|"
				end
				if edgesLabelsContentParentGeoEdgesIds
					labelContent *= string(meshData.edgesParentGeoEdgesIds[edgeId])*"|"
				end
				if edgesLabelsContentParentDomainsTypes
					labelContent *= string(meshData.edgesParentDomainsTypes[edgeId])*"|"
				end
				if edgesLabelsContentParentDomainsIds
					labelContent *= string(meshData.edgesParentDomainsIds[edgeId])*"|"
				end
				if edgesLabelsContentParentGammasIds
					labelContent *= string(meshData.edgesParentGammasIds[edgeId])*"|"
				end
				if edgesLabelsContentParentGammasParameterRanges
					labelContent *= string(round(collect(meshData.edgesParentGammasParameterRanges[:,edgeId]),2))*"|"
				end
				if edgesLabelsContentVolumes
					labelContent *= string(round(meshData.edgesVolumes[edgeId],2))*"|"
				end
				labelContent = labelContent[1:1:(end-1)]
				
				# Build label position:
				labelPositionX = sum(edgePositionXs)/2 + edgesLabelsOffset
				labelPositionY = sum(edgePositionYs)/2 + edgesLabelsOffset
				
				# Plot label:
				text(labelPositionX,labelPositionY,
					labelContent,
					color=currentEdgeColor,
					figure=figure("MeshData"))
			end
		end
	end
	
	if elementsVisible
		for elementId = 1:1:size(meshData.elements,2)
			for pairingId = 1:1:3
				if pairingId==1
					nodeId1 = meshData.elements[1,elementId]
					nodeId2 = meshData.elements[2,elementId]
				elseif pairingId==2
					nodeId1 = meshData.elements[2,elementId]
					nodeId2 = meshData.elements[3,elementId]
				elseif pairingId==3
					nodeId1 = meshData.elements[3,elementId]
					nodeId2 = meshData.elements[1,elementId]
				end
				
				if !hasBeenPlotted[nodeId1,nodeId2]
					# Plot data:
					edgePositionXs = [meshData.nodes[1,nodeId1],meshData.nodes[1,nodeId2]]
					edgePositionYs = [meshData.nodes[2,nodeId1],meshData.nodes[2,nodeId2]]
					plot(edgePositionXs,edgePositionYs,
						color=elementsColor,
						linestyle=elementsLineStyle,
						linewidth=elementsLineWidth,
						figure=figure("MeshData"))
					
					hasBeenPlotted[nodeId1,nodeId2] = true
					hasBeenPlotted[nodeId2,nodeId1] = true
				end
			end
			
			if elementsLabelsVisible
				
				# Build label content:
				labelContent = ""
				if elementsLabelsContentId
					labelContent *= string(elementId)*"|"
				end
				if elementsLabelsContentNodesIds
					labelContent *= string(collect(meshData.elements[:,elementId]))*"|"
				end
				if elementsLabelsContentParentSubdomainsIntersectionsIds
					labelContent *= string(meshData.elementsParentSubdomainsIntersectionsIds[elementId])*"|"
				end
				if elementsLabelsContentParentDomainsTypes
					labelContent *= string(meshData.elementsParentDomainsTypes[elementId])*"|"
				end
				if elementsLabelsContentParentDomainsIds
					labelContent *= string(meshData.elementsParentDomainsIds[elementId])*"|"
				end
				if elementsLabelsContentDetNablaFTs
					labelContent *= string(round(meshData.elementsDetNablaFTs[elementId],2))*"|"
				end
				if elementsLabelsContentNablaFTInv11s
					labelContent *= string(round(meshData.elementsNablaFTInv11s[elementId],2))*"|"
				end
				if elementsLabelsContentNablaFTInv21s
					labelContent *= string(round(meshData.elementsNablaFTInv21s[elementId],2))*"|"
				end
				if elementsLabelsContentNablaFTInv12s
					labelContent *= string(round(meshData.elementsNablaFTInv12s[elementId],2))*"|"
				end
				if elementsLabelsContentNablaFTInv22s
					labelContent *= string(round(meshData.elementsNablaFTInv22s[elementId],2))*"|"
				end
				if elementsLabelsContentVolumes
					labelContent *= string(round(meshData.elementsVolumes[elementId],2))*"|"
				end
				if elementsLabelsContentCentroids
					labelContent *= string(round(collect(meshData.elementsCentroids[:,elementId]),2))*"|"
				end
				if elementsLabelsContentMeshsizes
					labelContent *= string(round(meshData.elementsMeshsizes[elementId],2))*"|"
				end
				labelContent = labelContent[1:1:(end-1)]
				
				# Build label position:
				nodeId1 = meshData.elements[1,elementId]
				nodeId2 = meshData.elements[2,elementId]
				nodeId3 = meshData.elements[3,elementId]
				labelPositionX = sum([meshData.nodes[1,nodeId1],meshData.nodes[1,nodeId2],meshData.nodes[1,nodeId3]])/3 + elementsLabelsOffset
				labelPositionY = sum([meshData.nodes[2,nodeId1],meshData.nodes[2,nodeId2],meshData.nodes[2,nodeId3]])/3 + elementsLabelsOffset
				
				# Plot label:
				text(labelPositionX,labelPositionY,
					labelContent,
					color=elementsColor,
					figure=figure("MeshData"))
			end
		end
	end
	
	if nodesVisible
		for nodeId = 1:1:size(meshData.nodes,2)
			# Plot data:
			nodePositionX = meshData.nodes[1,nodeId]
			nodePositionY = meshData.nodes[2,nodeId]
			plot(nodePositionX,nodePositionY,
				color=nodesColor,
				markeredgecolor=nodesColor,
				marker=nodesMarkerType,
				markersize=nodesMarkerSize,
				figure=figure("MeshData"))
				
			if nodesLabelsVisible
				
				# Build label content:
				labelContent = ""
				if nodesLabelsContentId
					labelContent *= string(nodeId)*"|"
				end
				if nodesLabelsContentCoords
					labelContent *= string(round(collect(meshData.nodes[:,nodeId]),2))*"|"
				end
				labelContent = labelContent[1:1:(end-1)]
				
				# Build label position:
				labelPositionX = nodePositionX + nodesLabelsOffset
				labelPositionY = nodePositionY + nodesLabelsOffset
				
				# Plot label:
				text(labelPositionX,labelPositionY,
					labelContent,
					color=nodesColor,
					figure=figure("MeshData"))
			end
		end
	end
	
	axis("equal")
	title("MeshData")
	xlabel("x")
	ylabel("y")
	
	xMin,xMax = xlim()
	xMean = (xMin+xMax)/2
	xlim(xMean + 1.5*(xMin-xMean),xMean + 1.5*(xMax-xMean))
	
	if doSaveFigure
		savefig("./Saved images/meshData.png")
	end
	
	return nothing
end


