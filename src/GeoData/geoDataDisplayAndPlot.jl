
using PyPlot

function geoDataDisplay(geoData::GeoData)
	
	println("")
	println("---------------------------------------")
	println("geoData.geoNodes = ")
	display(geoData.geoNodes)
	println("")
	
	println("---------------------------------------")
	println("geoData.geoEdges = ")
	display(geoData.geoEdges)
	println("")
	
	println("---------------------------------------")
	println("geoData.geoEdgesLoops = ")
	display(geoData.geoEdgesLoops)
	println("")
	
	println("---------------------------------------")
	println("geoData.geoNodesIntersections = ")
	display(geoData.geoNodesIntersections)
	println("")
	
	println("---------------------------------------")
	println("geoData.subdomainsIntersections = ")
	display(geoData.subdomainsIntersections)
	println("")
	
	println("---------------------------------------")
	println("geoData.periodicGeoEdgesPairs = ")
	display(geoData.periodicGeoEdgesPairs)
	println("")
	
	println("---------------------------------------")
	println("")
	
	return nothing
end

function geoDataPlot(geoData::GeoData)
	
	doSaveFigure = true
	
	# Define geoNodes plot style:
	geoNodesVisible = true
	
	geoNodesColor = "red"
	geoNodesMarkerType = "o"
	geoNodesMarkerSize = 5
	
	geoNodesNonIsActivesVisible = true
	
	geoNodesLabelsVisible = false
	geoNodesLabelsOffset = 0.005
	geoNodesLabelsContentId = true
	geoNodesLabelsContentCoords = false
	geoNodesLabelsContentIsActive = false
	
	# Define geoEdges plot style:
	geoEdgesVisible = true
	
	geoEdgesColor = "black"
	geoEdgesLineStyle = "-"
	geoEdgesLineWidth = 2
	
	geoEdgesNonIsActivesVisible = true
	
	geoEdgesLabelsVisible = false
	geoEdgesLabelsOffset = 0.005
	geoEdgesLabelsContentId = true
	geoEdgesLabelsContentGeoNodesIds = false
	geoEdgesLabelsContentParentDomainType = false
	geoEdgesLabelsContentParentDomainId = false
	geoEdgesLabelsContentParentGammaId = false
	geoEdgesLabelsContentParentGammaParameterRange = false
	geoEdgesLabelsContentIsActive = false
	
	geoEdgesColorSubdomains = "blue"
	geoEdgesLineWidthSubdomains = 2
	
	# Define geoEdgesLoops plot style:
	geoEdgesLoopsVisible = true
	
	geoEdgesLoopsColor = geoEdgesColor
	geoEdgesLoopsColorSubdomains = geoEdgesColorSubdomains
	geoEdgesLoopsAlpha = 0.1
	
	geoEdgesLoopsBoundingBoxVisible = false
	geoEdgesLoopsBoundingBoxLinestyle = "--"
	geoEdgesLoopsBoundingBoxLinewidth = 1
	
	geoEdgesLoopsLabelsVisible = false
	geoEdgesLoopsLabelsOffset = 0.005
	geoEdgesLoopsLabelsContentId = true
	geoEdgesLoopsLabelsContentGeoEdgesIds = false
	geoEdgesLoopsLabelsContentParentDomainType = false
	geoEdgesLoopsLabelsContentParentDomainId = false
	geoEdgesLoopsLabelsContentOrientation = false
	geoEdgesLoopsLabelsContentBoundingBox = false
	geoEdgesLoopsLabelsContentSubsetRelationsToGeoEdgesLoops = false
	
	# Define subdomainsIntersections plot style:
	subdomainsIntersectionsVisible = false
	
	subdomainsIntersectionsColor = [0.0,1.0,0.0]
	subdomainsIntersectionsLabelsColor = 0.5*subdomainsIntersectionsColor
	subdomainsIntersectionsAlpha = 0.3
	
	subdomainsIntersectionsNonIsActivesVisible = true
	
	subdomainsIntersectionsBoundingBoxVisible = false
	subdomainsIntersectionsBoundingBoxLinestyle = "--"
	subdomainsIntersectionsBoundingBoxLinewidth = 1
	
	subdomainsIntersectionsLabelsVisible = false
	subdomainsIntersectionsLabelsOffset = 0.005
	subdomainsIntersectionsLabelsContentId = true
	subdomainsIntersectionsLabelsContentStarterGeoEdgesLoopId = false
	subdomainsIntersectionsLabelsContentStarterGeoEdgesLoopGeoEdgesIdsIndex = false
	subdomainsIntersectionsLabelsContentGeoEdgesIds = false
	subdomainsIntersectionsLabelsContentOrientation = false
	subdomainsIntersectionsLabelsContentBoundingBox = false
	subdomainsIntersectionsLabelsContentEnclosedSubdomainsIntersectionsIds = false
	subdomainsIntersectionsLabelsContentSubsetRelationsToSubdomainsIntersections = false
	subdomainsIntersectionsLabelsContentSubsetRelationsToGeoEdgesLoops = false
	subdomainsIntersectionsLabelsContentParentGeoEdgesLoopId = false
	subdomainsIntersectionsLabelsContentIsActive = false
	
	# Plot:
	figure("GeoData")
	clf()
	
	if geoEdgesLoopsVisible
		for geoEdgesLoop in geoData.geoEdgesLoops
			
			# Choose current color:
			if in(geoEdgesLoop.parentDomainType,["domain","hole"])
				currentSubdomainsIntersectionColor = geoEdgesLoopsColor
			elseif geoEdgesLoop.parentDomainType == "subdomain"
				currentSubdomainsIntersectionColor = geoEdgesLoopsColorSubdomains
			end
			
			# Plot boundingBox:
			if geoEdgesLoopsBoundingBoxVisible
				boundingBoxXs = [geoEdgesLoop.boundingBox[1],geoEdgesLoop.boundingBox[2],geoEdgesLoop.boundingBox[2],geoEdgesLoop.boundingBox[1],geoEdgesLoop.boundingBox[1]]
				boundingBoxYs = [geoEdgesLoop.boundingBox[3],geoEdgesLoop.boundingBox[3],geoEdgesLoop.boundingBox[4],geoEdgesLoop.boundingBox[4],geoEdgesLoop.boundingBox[3]]
				plot(boundingBoxXs,boundingBoxYs,
					color=currentSubdomainsIntersectionColor,
					linestyle=geoEdgesLoopsBoundingBoxLinestyle,
					linewidth=geoEdgesLoopsBoundingBoxLinewidth,
					figure=figure("GeoData"))
			end
			
			# Plot data:
			geoEdgesLoopPositionXs = map(geoNode->geoNode.coords[1],geoData.geoNodes[geoEdgesLoop.geoEdgesIds])
			geoEdgesLoopPositionYs = map(geoNode->geoNode.coords[2],geoData.geoNodes[geoEdgesLoop.geoEdgesIds])
			fill(geoEdgesLoopPositionXs,geoEdgesLoopPositionYs,
				color=currentSubdomainsIntersectionColor,
				alpha=geoEdgesLoopsAlpha,
				figure=figure("GeoData"))
			
			if geoEdgesLoopsLabelsVisible
				
				# Build label content:
				labelContent = ""
				if geoEdgesLoopsLabelsContentId
					labelContent *= string(geoEdgesLoop.id)*"|"
				end
				if geoEdgesLoopsLabelsContentGeoEdgesIds
					labelContent *= string(geoEdgesLoop.geoEdgesIds)*"|"
				end
				if geoEdgesLoopsLabelsContentParentDomainType
					labelContent *= string(geoEdgesLoop.parentDomainType)*"|"
				end
				if geoEdgesLoopsLabelsContentParentDomainId
					labelContent *= string(geoEdgesLoop.parentDomainId)*"|"
				end
				if geoEdgesLoopsLabelsContentOrientation
					labelContent *= string(geoEdgesLoop.orientation)*"|"
				end
				if geoEdgesLoopsLabelsContentBoundingBox
					labelContent *= string(round(geoEdgesLoop.boundingBox,2))*"|"
				end
				if geoEdgesLoopsLabelsContentSubsetRelationsToGeoEdgesLoops
					labelContent *= string(geoEdgesLoop.subsetRelationsToGeoEdgesLoops)*"|"
				end
				labelContent = labelContent[1:1:(end-1)]
				
				# Build label position:
				labelPositionX = geoEdgesLoopPositionXs[1] + geoEdgesLoopsLabelsOffset
				labelPositionY = geoEdgesLoopPositionYs[1] + geoEdgesLoopsLabelsOffset
				
				# Plot label:
				text(labelPositionX,labelPositionY,
					labelContent,
					color=currentSubdomainsIntersectionColor,
					figure=figure("GeoData"))
			end
			
		end
	end
	
	if subdomainsIntersectionsVisible
		for subdomainsIntersection in geoData.subdomainsIntersections
			if (subdomainsIntersectionsNonIsActivesVisible ? true : subdomainsIntersection.isActive)
				
				# Plot boundingBox:
				if subdomainsIntersectionsBoundingBoxVisible
					boundingBoxXs = [subdomainsIntersection.boundingBox[1],subdomainsIntersection.boundingBox[2],subdomainsIntersection.boundingBox[2],subdomainsIntersection.boundingBox[1],subdomainsIntersection.boundingBox[1]]
					boundingBoxYs = [subdomainsIntersection.boundingBox[3],subdomainsIntersection.boundingBox[3],subdomainsIntersection.boundingBox[4],subdomainsIntersection.boundingBox[4],subdomainsIntersection.boundingBox[3]]
					plot(boundingBoxXs,boundingBoxYs,
						color=subdomainsIntersectionsColor,
						linestyle=subdomainsIntersectionsBoundingBoxLinestyle,
						linewidth=subdomainsIntersectionsBoundingBoxLinewidth,
						figure=figure("GeoData"))
				end
				
				# Plot data:
				subdomainsIntersectionPositionXs = map(geoEdgeId->geoData.geoNodes[geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)]].coords[1],subdomainsIntersection.geoEdgesIds)
				subdomainsIntersectionPositionYs = map(geoEdgeId->geoData.geoNodes[geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)]].coords[2],subdomainsIntersection.geoEdgesIds)
				for enclosedSubdomainsIntersection in geoData.subdomainsIntersections[subdomainsIntersection.enclosedSubdomainsIntersectionsIds]
					enclosedSubdomainsIntersectionPositionXs = map(geoEdgeId->geoData.geoNodes[geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)]].coords[1],enclosedSubdomainsIntersection.geoEdgesIds)
					enclosedSubdomainsIntersectionPositionYs = map(geoEdgeId->geoData.geoNodes[geoData.geoEdges[abs(geoEdgeId)].geoNodesIds[(geoEdgeId>=0 ? 1 : 2)]].coords[2],enclosedSubdomainsIntersection.geoEdgesIds)
					if enclosedSubdomainsIntersection.orientation == subdomainsIntersection.orientation
						append!(subdomainsIntersectionPositionXs,
							vcat(subdomainsIntersectionPositionXs[1],
							enclosedSubdomainsIntersectionPositionXs[1],
							enclosedSubdomainsIntersectionPositionXs[end:(-1):1]))
						append!(subdomainsIntersectionPositionYs,
							vcat(subdomainsIntersectionPositionYs[1],
							enclosedSubdomainsIntersectionPositionYs[1],
							enclosedSubdomainsIntersectionPositionYs[end:(-1):1]))
					else
						append!(subdomainsIntersectionPositionXs,
							vcat(subdomainsIntersectionPositionXs[1],
							enclosedSubdomainsIntersectionPositionXs,
							enclosedSubdomainsIntersectionPositionXs[1]))
						append!(subdomainsIntersectionPositionYs,
							vcat(subdomainsIntersectionPositionYs[1],
							enclosedSubdomainsIntersectionPositionYs,
							enclosedSubdomainsIntersectionPositionYs[1]))
					end
				end
				fill(subdomainsIntersectionPositionXs,subdomainsIntersectionPositionYs,
					color=subdomainsIntersectionsColor,
					alpha=subdomainsIntersectionsAlpha,
					edgecolor="none",
					figure=figure("GeoData"))
				
				if subdomainsIntersectionsLabelsVisible
					
					# Build label content:
					labelContent = ""
					if subdomainsIntersectionsLabelsContentId
						labelContent *= string(subdomainsIntersection.id)*"|"
					end
					if subdomainsIntersectionsLabelsContentStarterGeoEdgesLoopId
						labelContent *= string(subdomainsIntersection.starterGeoEdgesLoopId)*"|"
					end
					if subdomainsIntersectionsLabelsContentStarterGeoEdgesLoopGeoEdgesIdsIndex
						labelContent *= string(subdomainsIntersection.starterGeoEdgesLoopGeoEdgesIdsIndex)*"|"
					end
					if subdomainsIntersectionsLabelsContentGeoEdgesIds
						labelContent *= string(subdomainsIntersection.geoEdgesIds)*"|"
					end
					if subdomainsIntersectionsLabelsContentOrientation
						labelContent *= string(subdomainsIntersection.orientation)*"|"
					end
					if subdomainsIntersectionsLabelsContentBoundingBox
						labelContent *= string(round(subdomainsIntersection.boundingBox,2))*"|"
					end
					if subdomainsIntersectionsLabelsContentEnclosedSubdomainsIntersectionsIds
						labelContent *= string(subdomainsIntersection.enclosedSubdomainsIntersectionsIds)*"|"
					end
					if subdomainsIntersectionsLabelsContentSubsetRelationsToSubdomainsIntersections
						labelContent *= string(subdomainsIntersection.subsetRelationsToSubdomainsIntersections)*"|"
					end
					if subdomainsIntersectionsLabelsContentSubsetRelationsToGeoEdgesLoops
						labelContent *= string(subdomainsIntersection.subsetRelationsToGeoEdgesLoops)*"|"
					end
					if subdomainsIntersectionsLabelsContentParentGeoEdgesLoopId
						labelContent *= string(subdomainsIntersection.parentGeoEdgesLoopId)*"|"
					end
					if subdomainsIntersectionsLabelsContentIsActive
						labelContent *= string(subdomainsIntersection.isActive)*"|"
					end
					labelContent = labelContent[1:1:(end-1)]
					
					# Build label position:
					labelPositionX = subdomainsIntersectionPositionXs[1] + subdomainsIntersectionsLabelsOffset
					labelPositionY = subdomainsIntersectionPositionYs[1] + subdomainsIntersectionsLabelsOffset
					
					# Plot label:
					text(labelPositionX,labelPositionY,
						labelContent,
						color=subdomainsIntersectionsLabelsColor,
						figure=figure("GeoData"))
				end
			end
		end
	end
	
	if geoEdgesVisible
		for geoEdge in geoData.geoEdges
			if (geoEdgesNonIsActivesVisible ? true : geoEdge.isActive)
				
				# Choose current color, linewidth, label visibility:
				if in(geoEdge.parentDomainType,["domain","hole"])
					currentGeoEdgeColor = geoEdgesColor
					currentGeoEdgeLineWidth = geoEdgesLineWidth
				elseif geoEdge.parentDomainType == "subdomain"
					currentGeoEdgeColor = geoEdgesColorSubdomains
					currentGeoEdgeLineWidth = geoEdgesLineWidthSubdomains
				end
				
				# Plot data:
				geoEdgePositionXs = [geoData.geoNodes[geoEdge.geoNodesIds[1]].coords[1],geoData.geoNodes[geoEdge.geoNodesIds[2]].coords[1]]
				geoEdgePositionYs = [geoData.geoNodes[geoEdge.geoNodesIds[1]].coords[2],geoData.geoNodes[geoEdge.geoNodesIds[2]].coords[2]]
				plot(geoEdgePositionXs,geoEdgePositionYs,
					color=currentGeoEdgeColor,
					linestyle=geoEdgesLineStyle,
					linewidth=currentGeoEdgeLineWidth,
					figure=figure("GeoData"))
				
				if geoEdgesLabelsVisible
					
					# Build label content:
					labelContent = ""
					if geoEdgesLabelsContentId
						labelContent *= string(geoEdge.id)*"|"
					end
					if geoEdgesLabelsContentGeoNodesIds
						labelContent *= string(geoEdge.geoNodesIds)*"|"
					end
					if geoEdgesLabelsContentParentDomainType
						labelContent *= string(geoEdge.parentDomainType)*"|"
					end
					if geoEdgesLabelsContentParentDomainId
						labelContent *= string(geoEdge.parentDomainId)*"|"
					end
					if geoEdgesLabelsContentParentGammaId
						labelContent *= string(geoEdge.parentGammaId)*"|"
					end
					if geoEdgesLabelsContentParentGammaParameterRange
						labelContent *= string(round(geoEdge.parentGammaParameterRange,2))*"|"
					end
					if geoEdgesLabelsContentIsActive
						labelContent *= string(geoEdge.isActive)*"|"
					end
					labelContent = labelContent[1:1:(end-1)]
					
					# Build label position:
					labelPositionX = sum(geoEdgePositionXs)/2 + geoEdgesLabelsOffset
					labelPositionY = sum(geoEdgePositionYs)/2 + geoEdgesLabelsOffset
					
					# Plot label:
					text(labelPositionX,labelPositionY,
						labelContent,
						color=currentGeoEdgeColor,
						figure=figure("GeoData"))
				end
			end
		end
	end
	
	if geoNodesVisible
		for geoNode in geoData.geoNodes
			if (geoNodesNonIsActivesVisible ? true : geoNode.isActive)
				# Plot data:
				geoNodePositionX = geoNode.coords[1]
				geoNodePositionY = geoNode.coords[2]
				plot(geoNodePositionX,geoNodePositionY,
					color=geoNodesColor,
					markeredgecolor=geoNodesColor,
					marker=geoNodesMarkerType,
					markersize=geoNodesMarkerSize,
					figure=figure("GeoData"))
				
				if geoNodesLabelsVisible
					# Build label content:
					labelContent = ""
					if geoNodesLabelsContentId
						labelContent *= string(geoNode.id)*"|"
					end
					if geoNodesLabelsContentCoords
						labelContent *= string(round(geoNode.coords,2))*"|"
					end
					if geoNodesLabelsContentIsActive
						labelContent *= string(geoNode.isActive)*"|"
					end
					labelContent = labelContent[1:1:(end-1)]
					
					# Build label position:
					labelPositionX = geoNodePositionX + geoNodesLabelsOffset
					labelPositionY = geoNodePositionY + geoNodesLabelsOffset
					
					# Plot label:
					text(labelPositionX,labelPositionY,
						labelContent,
						color=geoNodesColor,
						figure=figure("GeoData"))
				end
			end
		end
	end
	
	axis("equal")
	title("GeoData")
	xlabel("x")
	ylabel("y")
	
	xMin,xMax = xlim()
	xMean = (xMin+xMax)/2
	xlim(xMean + 1.5*(xMin-xMean),xMean + 1.5*(xMax-xMean))
	
	if doSaveFigure
		savefig("./Saved images/geoData.png")
	end
	
	return nothing
end


