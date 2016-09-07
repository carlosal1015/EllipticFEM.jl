
function evaluateConvexCombination(a::Float64,b::Float64,s::Float64)
	return (1-s)*a + s*b
end

function evaluateConvexCombination(a::Array{Float64,1},b::Array{Float64,1},s::Float64)
	return (1-s)*a + s*b
end

function rotate2DPoint(centerCoords::Array{Float64,1},angle::Float64,pointCoords::Array{Float64,1})
	return centerCoords + hcat([cos(angle),sin(angle)],[-sin(angle),cos(angle)])*(pointCoords-centerCoords)
end

function valueIsInInterval(value::Float64,intervalCoords::Array{Float64})
	return (intervalCoords[1]<=value)&(value<=intervalCoords[2])
end

function intervalsDoIntersect(interval1Coords::Array{Float64},interval2Coords::Array{Float64})
	return valueIsInInterval(interval1Coords[1],interval2Coords)|
		valueIsInInterval(interval1Coords[2],interval2Coords)|
		valueIsInInterval(interval2Coords[1],interval1Coords)|
		valueIsInInterval(interval2Coords[1],interval1Coords)
end

function cross2D(a::Array{Float64},b::Array{Float64})
	return a[1]*b[2] - a[2]*b[1]
end

function orientedAngle(a::Array{Float64},b::Array{Float64})
	# Computes the oriented angle from vector a to vector b.
	# E.g.: orientedAngle([1.0,0.0],[0.0,1.0]) == +pi/2.
	
	angle = acos(projectValueIntoInterval(dot(b,a)/(norm(b)*norm(a)),[-1.0,1.0]))
	
	if cross2D(a,b) < 0
		angle = 2*pi - angle
	end
	
	return angle
end

#######################################################################

function intersectLineWithLine(A::Array{Float64},a::Array{Float64},B::Array{Float64},b::Array{Float64})
	
	# TODO: Set threshold properly.
	thresholdDeterminant = 0.000001
	thresholdFourierCoefficient = 0.000001
	
	matrix = hcat(a,-b)
	
	if abs(det(matrix)-0) <= thresholdDeterminant
		if abs(dot(B-A,[-a[2],a[1]])/(norm(a)^2)) <= thresholdFourierCoefficient
			
			alpha = dot(B-A,a)/(norm(a)^2)		# (A + alpha*a == B + 0*b)
			beta = dot(A-B,b)/(norm(b)^2)		# (B + beta*b == A + 0*a)
			intersectionPoint1Coords = A
			intersectionPoint2Coords = B
			
			intersectionRelation = "collinear"
			parameters = Array{Array{Float64,1},1}(2)
			parameters[1] = [0,beta]
			parameters[2] = [alpha,0]
			intersectionPointsCoords = Array{Array{Float64,1},1}(2)
			intersectionPointsCoords[1] = intersectionPoint1Coords
			intersectionPointsCoords[2] = intersectionPoint2Coords
			
		else
			intersectionRelation = "parallel"
			parameters = Array{Array{Float64,1},1}(0)
			intersectionPointsCoords = Array{Array{Float64,1},1}(0)
		end
	else
		alphaBeta = matrix\(B-A)		# (A + alphaBeta[1]*a == B + alphaBeta[2]*b)
		intersectionPointCoords = A + alphaBeta[1]*a
		
		intersectionRelation = "intersecting"
		parameters = Array{Array{Float64,1},1}(1)
		parameters[1] = alphaBeta
		intersectionPointsCoords = Array{Array{Float64,1},1}(1)
		intersectionPointsCoords[1] = intersectionPointCoords
	end
	
	return intersectionRelation,parameters,intersectionPointsCoords
end

#######################################################################

function computePolygonBoundingBox(polygonVerticesCoords::Array{Array{Float64,1},1})
	
	polygonVerticesXs = map(polygonVertexCoords->polygonVertexCoords[1],polygonVerticesCoords)
	polygonVerticesYs = map(polygonVertexCoords->polygonVertexCoords[2],polygonVerticesCoords)
	
	polygonBoundingBox = [minimum(polygonVerticesXs),maximum(polygonVerticesXs),minimum(polygonVerticesYs),maximum(polygonVerticesYs)]
	
	return polygonBoundingBox
end

function pointIsInBoundingBox(A::Array{Float64},boundingBox::Array{Float64})
	
	# TODO: Set threshold properly.
	threshold = 0.000001
	
	if valueIsInInterval(A[1],boundingBox[1:1:2]+threshold*[-1,1])&
		valueIsInInterval(A[2],boundingBox[3:1:4]+threshold*[-1,1])
		
		isInside = true
	else
		isInside = false
	end
	
	return isInside
end

function boundingBoxesDoIntersect(boundingBox1::Array{Float64},boundingBox2::Array{Float64})
	
	# TODO: Set threshold properly.
	threshold = 0.000001
	
	if intervalsDoIntersect(boundingBox1[1:1:2]+threshold*[-1,1],boundingBox2[1:1:2]+threshold*[-1,1])&
		intervalsDoIntersect(boundingBox1[3:1:4]+threshold*[-1,1],boundingBox2[3:1:4]+threshold*[-1,1])
		
		doIntersect = true
	else
		doIntersect = false
	end
	
	return doIntersect
end

function computePolygonSignedArea(polygonVerticesCoords::Array{Array{Float64,1},1})
	
	sum = 0
	for i = 1:1:length(polygonVerticesCoords)
		sum += cross2D(polygonVerticesCoords[i],polygonVerticesCoords[moduloIndex(i+1,length(polygonVerticesCoords))])
	end
	
	signedArea = 1/2*sum
	
	return signedArea
end

function computePolygonOrientation(polygonVerticesCoords::Array{Array{Float64,1},1})
	
	signedArea = computePolygonSignedArea(polygonVerticesCoords)
	orientation = (signedArea>=0 ? "positive" : "negative")
	
	return orientation
end

function computePolygonCenterOfMass(polygonVerticesCoords::Array{Array{Float64,1},1})
	
	sumX = 0
	sumY = 0
	for i = 1:1:length(polygonVerticesCoords)
		sumX += (polygonVerticesCoords[i][1] + polygonVerticesCoords[moduloIndex(i+1,length(polygonVerticesCoords))][1])*
			cross2D(polygonVerticesCoords[i],polygonVerticesCoords[moduloIndex(i+1,length(polygonVerticesCoords))])
		sumY += (polygonVerticesCoords[i][2] + polygonVerticesCoords[moduloIndex(i+1,length(polygonVerticesCoords))][2])*
			cross2D(polygonVerticesCoords[i],polygonVerticesCoords[moduloIndex(i+1,length(polygonVerticesCoords))])
	end
	
	centerOfMassCoords = (1/(6*computePolygonSignedArea(polygonVerticesCoords)))*[sumX,sumY]
	
	return centerOfMassCoords
end

function computeIsolatedValue(values::Array{Float64},intervalCoords::Array{Float64},periodicBoundaries::Bool)
	
	valuesSorted = sort(values)
	if !periodicBoundaries
		valuesAndBoundariesSorted = vcat(intervalCoords[1],valuesSorted,intervalCoords[2])
	else
		valuesAndBoundariesSorted = vcat(valuesSorted,intervalCoords[2]+(valuesSorted[1]-intervalCoords[1]))
	end
	
	(maxDifference,maxIndex) = findmax(valuesAndBoundariesSorted[2:1:end] - valuesAndBoundariesSorted[1:1:(end-1)])
	isolatedValue = valuesAndBoundariesSorted[maxIndex] + (1/2)*maxDifference
	
	if periodicBoundaries
		isolatedValue = mod(isolatedValue-intervalCoords[1],intervalCoords[2]-intervalCoords[1]) + intervalCoords[1]
	end
	
	return isolatedValue
end

function computeIntersectionPointsBetweenRayAndPolygon(rayPointCoords::Array{Float64},rayDirectionVectorCoords::Array{Float64},polygonVerticesCoords::Array{Array{Float64,1},1})
	
	# TODO: Do we need a threshold here?
	
	intersectionPointsCoords = Array{Array{Float64,1},1}(0)
	
	for i = 1:1:length(polygonVerticesCoords)
		intersectionRelationCurrentEdge,parametersCurrentEdge,intersectionPointsCoordsCurrentEdge = intersectLineWithLine(
			rayPointCoords,
			rayDirectionVectorCoords,
			polygonVerticesCoords[i],
			polygonVerticesCoords[moduloIndex(i+1,length(polygonVerticesCoords))] - polygonVerticesCoords[i])
		
		if intersectionRelationCurrentEdge == "intersecting"
			if (0<=parametersCurrentEdge[1][1]) & valueIsInInterval(parametersCurrentEdge[1][2],[0.0,1.0])
				push!(intersectionPointsCoords,intersectionPointsCoordsCurrentEdge[1])
			end
		elseif intersectionRelationCurrentEdge == "parallel"
			# Do nothing.
		elseif intersectionRelationCurrentEdge == "collinear"
			if (0<=parametersCurrentEdge[1][1])
				# Do nothing.
			end
			warn("This should not happen. ray and geoEdge nearly parallel. "*
				"rayPointCoords = "*string(round(rayPointCoords,3))*
				", rayDirectionVectorCoords = "*string(round(rayDirectionVectorCoords,3))*
				", polygonVerticesCoords[i] = "*string(round(polygonVerticesCoords[i],3))*".")
		end
	end
	
	return intersectionPointsCoords
end

function computeRayDirectionVectorCoordsForJordanTest(pointCoords::Array{Float64},polygonVerticesCoords::Array{Array{Float64,1},1})
	
	angles1 = Array{Float64}(length(polygonVerticesCoords))
	angles2 = Array{Float64}(length(polygonVerticesCoords))
	angles3 = Array{Float64}(length(polygonVerticesCoords))
	for i = 1:1:length(polygonVerticesCoords)
		angle = orientedAngle([1.0,0.0],
			polygonVerticesCoords[moduloIndex(i+1,length(polygonVerticesCoords))]-polygonVerticesCoords[i])
		angles1[i] = angle
		angles2[i] = mod(angle-pi,2*pi)
		
		angles3[i] = orientedAngle([1.0,0.0],polygonVerticesCoords[i]-pointCoords)
	end
	
	rayDirectionVectorAngle = computeIsolatedValue(vcat(angles1,angles2,angles3),[0.0,2*pi],true)
	rayDirectionVectorCoords = [cos(rayDirectionVectorAngle),sin(rayDirectionVectorAngle)]
	
	return rayDirectionVectorCoords
end

function performJordanPointInPolygonTest(pointCoords::Array{Float64},polygonVerticesCoords::Array{Array{Float64,1},1})
	
	rayDirectionVectorCoords = computeRayDirectionVectorCoordsForJordanTest(pointCoords,polygonVerticesCoords)
	intersectionPointsCoords = computeIntersectionPointsBetweenRayAndPolygon(pointCoords,rayDirectionVectorCoords,polygonVerticesCoords)
	subsetRelation = (isodd(length(intersectionPointsCoords)) ? "strictInside" : "strictOutside")
	
	return subsetRelation
end

function computeSomePointInPolygon(polygonBoundingBox::Array{Float64},polygonVerticesCoords::Array{Array{Float64,1},1})
	
	# TODO: Add bounding box checks.
	
	centerOfMassCoords = computePolygonCenterOfMass(polygonVerticesCoords)
	subsetRelation = performJordanPointInPolygonTest(centerOfMassCoords,polygonVerticesCoords)
	
	if subsetRelation == "strictInside"
		
		somePointCoords = centerOfMassCoords
		
	elseif subsetRelation == "strictOutside"
		
		rayPointCoords = [polygonBoundingBox[1]-0.1*(polygonBoundingBox[2]-polygonBoundingBox[1]),
			computeIsolatedValue(map(polygonVertexCoords->polygonVertexCoords[2],polygonVerticesCoords),polygonBoundingBox[3:1:4],false)]
		rayDirectionVectorCoords = [1.0,0.0]
		intersectionPointsCoords = computeIntersectionPointsBetweenRayAndPolygon(rayPointCoords,rayDirectionVectorCoords,polygonVerticesCoords)
		
		distances = Array{Float64}(length(intersectionPointsCoords))
		for i = 1:1:length(intersectionPointsCoords)
			distances[i] = norm(intersectionPointsCoords[i]-rayPointCoords)
		end
		indices = sortperm(distances)
		
		somePointCoords = (1/2)*(intersectionPointsCoords[indices[1]] + intersectionPointsCoords[indices[2]])
		
	end
	
	return somePointCoords
end

function computeIntersectionPointsBetweenTwoPolygons(polygon1VerticesCoords::Array{Array{Float64,1},1},polygon1BoundingBox::Array{Float64,1},polygon2VerticesCoords::Array{Array{Float64,1},1},polygon2BoundingBox::Array{Float64,1})
	
	# TODO: Set threshold properly.
	threshold = 0.000001
	
	polygonsDoIntersect = false
	
	polygon1IntersectionPointsPredecessorsIds = Array{Int64,1}(0)
	polygon1IntersectionPointsParameterValues = Array{Float64,1}(0)
	polygon1IntersectionPointsCoords = Array{Array{Float64,1}}(0)
	polygon2IntersectionPointsPredecessorsIds = Array{Int64,1}(0)
	polygon2IntersectionPointsParameterValues = Array{Float64,1}(0)
	polygon2IntersectionPointsCoords = Array{Array{Float64,1}}(0)
	
	if boundingBoxesDoIntersect(polygon1BoundingBox,polygon2BoundingBox)
		
		for i = 1:1:length(polygon1VerticesCoords)
			for j = 1:1:length(polygon2VerticesCoords)
				
				A = polygon1VerticesCoords[i]
				a = polygon1VerticesCoords[moduloIndex(i+1,length(polygon1VerticesCoords))] - polygon1VerticesCoords[i]
				B = polygon2VerticesCoords[j]
				b = polygon2VerticesCoords[moduloIndex(j+1,length(polygon2VerticesCoords))] - polygon2VerticesCoords[j]
				
				intersectionRelation,parameters,intersectionPointsCoords = intersectLineWithLine(A,a,B,b)
				
				if intersectionRelation == "intersecting"
					if valueIsInInterval(parameters[1][1],[0.0,1.0]+threshold*[-1,-1])&
						valueIsInInterval(parameters[1][2],[0.0,1.0]+threshold*[-1,-1])
						
						push!(polygon1IntersectionPointsPredecessorsIds,i)
						push!(polygon1IntersectionPointsParameterValues,parameters[1][1])
						push!(polygon1IntersectionPointsCoords,intersectionPointsCoords[1])
						
						push!(polygon2IntersectionPointsPredecessorsIds,j)
						push!(polygon2IntersectionPointsParameterValues,parameters[1][2])
						push!(polygon2IntersectionPointsCoords,intersectionPointsCoords[1])
						
						polygonsDoIntersect = true
						
					end
				end
				
			end
		end
		
	end
	
	return polygonsDoIntersect,
		polygon1IntersectionPointsPredecessorsIds,polygon1IntersectionPointsParameterValues,polygon1IntersectionPointsCoords,
		polygon2IntersectionPointsPredecessorsIds,polygon2IntersectionPointsParameterValues,polygon2IntersectionPointsCoords
end

###############################################

# For debugging purposes:

function randomPolygon(numberOfVertices::Int64,maySelfIntersect::Bool)
	
	polygonVerticesCoords = Array{Array{Float64,1},1}(numberOfVertices)
	
	if maySelfIntersect
		
		for i = 1:1:length(polygonVerticesCoords)
			polygonVerticesCoords[i] = rand(2)
		end
		
	else
		
		counterMax = 100
		for counter = 1:1:counterMax
			for i = 1:1:length(polygonVerticesCoords)
				polygonVerticesCoords[i] = rand(2)
			end
			
			if !polygonIsSelfIntersecting(polygonVerticesCoords)
				println("counter (end) = "*string(counter))
				break
			end
			
			if counter==counterMax
				warn("Polygon is self-intersecting.")
			end
		end
		
	end
	
	return polygonVerticesCoords
end

function polygonIsSelfIntersecting(polygonVerticesCoords::Array{Array{Float64,1},1})
	
	isSelfIntersecting = false
	breakOuterLoop = false
	for i = 1:1:length(polygonVerticesCoords)
		for j = (i+2):1:(length(polygonVerticesCoords) + (i==1 ? -1 : 0))
			
			A = polygonVerticesCoords[i]
			a = polygonVerticesCoords[moduloIndex(i+1,length(polygonVerticesCoords))] - polygonVerticesCoords[i]
			B = polygonVerticesCoords[j]
			b = polygonVerticesCoords[moduloIndex(j+1,length(polygonVerticesCoords))] - polygonVerticesCoords[j]
			
			intersectionRelation,parameters,intersectionPointCoords = intersectLineWithLine(A,a,B,b)
			
			# TODO: Do we need a threshold here?
			if intersectionRelation == "intersecting"
				if valueIsInInterval(parameters[1][1],[0.0,1.0])&
					valueIsInInterval(parameters[1][2],[0.0,1.0])
					
					isSelfIntersecting = true
					breakOuterLoop = true
					break
				end
			end
		end
		if breakOuterLoop
			break
		end
	end
	
	return isSelfIntersecting
end

function plotPolygonOutline(polygonVerticesCoords::Array{Array{Float64,1},1})
	
	figure("Polygon")
	clf()
	axis("equal")
	
	xs = map(polygonVertexCoords->polygonVertexCoords[1],polygonVerticesCoords)
	ys = map(polygonVertexCoords->polygonVertexCoords[2],polygonVerticesCoords)
	plot(vcat(xs,xs[1]),vcat(ys,ys[1]),color="black",linewidth=2,figure=figure("Polygon"))
	
	for i = 1:1:length(polygonVerticesCoords)
		text(polygonVerticesCoords[i][1],polygonVerticesCoords[i][2],
			string(i),
			color="black",
			figure=figure("Polygon"))
	end
	
	return nothing
end

function plotPolygonBoundingBox(polygonBoundingBox::Array{Float64})
	
	xs = [polygonBoundingBox[1],polygonBoundingBox[2],polygonBoundingBox[2],polygonBoundingBox[1]]
	ys = [polygonBoundingBox[3],polygonBoundingBox[3],polygonBoundingBox[4],polygonBoundingBox[4]]
	plot(vcat(xs,xs[1]),vcat(ys,ys[1]),color="black",linewidth=1,figure=figure("Polygon"))
	
	return nothing
end

function plotPolygonSomePoint(somePointCoords::Array{Float64})
	
	plot(somePointCoords[1],somePointCoords[2],color="yellow",marker="o",markersize=5,figure=figure("Polygon"))
	
	return nothing
end

function plotPolygonRayAndIntersectionPoints(rayPointCoords::Array{Float64},rayDirectionVectorCoords::Array{Float64},intersectionPointsCoords::Array{Array{Float64,1},1})
	
	xs = rayPointCoords[1]+[0,rayDirectionVectorCoords[1]]
	ys = rayPointCoords[2]+[0,rayDirectionVectorCoords[2]]
	plot(xs,ys,color="red",linewidth=1,figure=figure("Polygon"))
	
	plot(rayPointCoords[1],rayPointCoords[2],color="red",marker="o",markersize=10,figure=figure("Polygon"))
	
	for i = 1:1:length(intersectionPointsCoords)
		plot(intersectionPointsCoords[i][1],intersectionPointsCoords[i][2],color="green",marker="o",markersize=10,figure=figure("Polygon"))
		text(intersectionPointsCoords[i][1],intersectionPointsCoords[i][2],
			string(i),
			color="green",
			figure=figure("Polygon"))
	end
	
	return nothing
end


