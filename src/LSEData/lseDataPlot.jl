
function lseDataPlot(meshData::MeshData,lseData::LSEData)
	
	dimensionChooser = "2D"
	
	doToggleFullscreen = false
	doRedrawCanvas = true
	doSaveFigures = false
	doResizeFrames = true
	doConvertToVideo = true
	doPlayVideo = false
	doRemoveFrames = true
	
	doPlotColorbar = true
	
	# Compute color limits:
	colorLimits = [minimum(map(solution->minimum(solution),lseData.solutions)),
		maximum(map(solution->maximum(solution),lseData.solutions))]
	
	figure("Solution")
	clf()
	
	# Toggle window fullscreen:
	if doToggleFullscreen
		get_current_fig_manager()[:full_screen_toggle]()
	end
	
	if dimensionChooser == "2D"
		
		plotHandle = tripcolor(collect(meshData.nodes[1,:]),collect(meshData.nodes[2,:]),meshData.elements'-1,
			lseData.solutions[1],
			cmap=ColorMap("jet"),
			vmin=colorLimits[1],
			vmax=colorLimits[2],
			figure=figure("Solution"),
			shading="gouraud")
		
	elseif dimensionChooser == "3D"
		
		plotHandle = plot_trisurf(collect(meshData.nodes[1,:]),collect(meshData.nodes[2,:]),meshData.elements'-1,
			lseData.solutions[1],
			cmap=ColorMap("jet"),
			vmin=colorLimits[1],
			vmax=colorLimits[2],
			edgecolor="None",
			figure=figure("Solution"),
			shade=true)
		
		zlim(colorLimits)
		
	end
	
	if doPlotColorbar
		colorbar(plotHandle)
	end
	
	axis("equal")
	gca()[:set_title]("Solution",loc="left")
	xlabel("x")
	ylabel("y")
	
	for i = 1:1:length(lseData.solutions)
		
		println("t = "*string(round(lseData.solutionsTimestamps[i],2)))
		
		# Set facecolors according to solutions values:
		if dimensionChooser == "2D"
			faceColors = plotHandle[:to_rgba](lseData.solutions[i])
			plotHandle[:set_facecolors](faceColors)
		end
		
		# Set z-data according to solutions values:
		if dimensionChooser == "3D"
			
			# TODO: Gouraud shading does not work in 3D ?
# 			faceColors = plotHandle[:to_rgba](lseData.solutions[i])
# 			plotHandle[:set_facecolors](faceColors)
# 			plotHandle[:set_array](lseData.solutions[i])
			
			clf()
			plotHandle = plot_trisurf(collect(meshData.nodes[1,:]),collect(meshData.nodes[2,:]),meshData.elements'-1,
				lseData.solutions[i],
				cmap=ColorMap("jet"),
				vmin=colorLimits[1],
				vmax=colorLimits[2],
				edgecolor="None",
				figure=figure("Solution"))
			zlim(colorLimits)
			
		end
		
		# Set title:
		gca()[:set_title]("t = "*string(round(lseData.solutionsTimestamps[i],2))*" s",loc="left")
		
		# Make changes visible in current figure:
		if doRedrawCanvas
			gcf()[:canvas][:draw]()
		end
		
		# Save figure:
		if doSaveFigures
			savefig("./Saved images/frame"*dec(i,3)*".png")
		end
	end
	
	# Toggle window fullscreen:
	if doToggleFullscreen
		get_current_fig_manager()[:full_screen_toggle]()
	end
	
	# Resize frames:
	if doResizeFrames
		
		modeChooser = "evenWidthAndHeight"
		
		width = 812
		height = 612
		
		cd("./Saved images")
		if modeChooser == "fixedWidth"
			readall(`avconv -loglevel quiet -i frame%03d.png -vf scale=$width:-1 frame%03d.png`)
		elseif modeChooser == "fixedHeight"
			readall(`avconv -loglevel quiet -i frame%03d.png -vf scale=-1:$height frame%03d.png`)
		elseif modeChooser == "fixedWidthAndHeight"
			readall(`avconv -loglevel quiet -i frame%03d.png -vf scale=$width:$height frame%03d.png`)
		elseif modeChooser == "evenWidthAndHeight"
			readall(`avconv -loglevel quiet -i frame%03d.png -vf scale=trunc(iw/2)*2:trunc(ih/2)*2 frame%03d.png`)
		end
		cd("..")
	end
	
	# Convert frames to a video:
	if doConvertToVideo
		cd("./Saved images")
		readall(`avconv -loglevel quiet -i frame%03d.png video.mp4 -y`)
		cd("..")
	end
	
	# Play the video:
	if doPlayVideo
		# TODO: Remove warnings and errors.
		cd("./Saved images")
		readall(`vlc video.mp4 --quiet`)
		cd("..")
	end
	
	# Remove frames:
	if doRemoveFrames
		# TODO: Remove frames without loop.
		cd("./Saved images")
		allFilesNames = readdir()
		framesFilesNames = allFilesNames[map(fileName->fileName[1:1:5]=="frame",allFilesNames)]
		for fileName in framesFilesNames
			rm(fileName)
		end
		cd("..")
	end
	
	return nothing
end

