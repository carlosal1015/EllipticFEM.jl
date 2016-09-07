##############################################
# coefficients in a (sub-)domain:
##############################################

uInitial = (x,y,t) -> 0.0
uTInitial = (x,y,t) -> 0.0

A11 = (x,y,t) -> 1.0
A21 = (x,y,t) -> 0.0
A12 = (x,y,t) -> 0.0
A22 = (x,y,t) -> 1.0
A = hcat([A11,A21],[A12,A22])
b1 = (x,y,t) -> 0.0
b2 = (x,y,t) -> 0.0
b = [b1,b2]
c = (x,y,t) -> 0.0
f = (x,y,t) -> 0.0
fDiv1 = (x,y,t) -> 0.0
fDiv2 = (x,y,t) -> 0.0
fDiv = [fDiv1,fDiv2]

isZeroUInitial = true
isZeroUTInitial = true

isIdentityA = true
isZeroB = true
isZeroC = true
isZeroF = true
isZeroFDiv = true

isConstantInTimeA = true
isConstantInTimeB = true
isConstantInTimeC = true
isConstantInTimeF = true
isConstantInTimeFDiv = true

domain = Domain(1,"domain",gammaNodes,gammas,[],
	uInitial,uTInitial,A,b,c,f,fDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# rectangle:
##############################################

x0 = 0.0
y0 = 0.0
dx = 1.0
dy = 1.0

gammaNode1 = GammaNode(1,[x0,y0]+[0.0,0.0])
gammaNode2 = GammaNode(2,[x0,y0]+[dx,0.0])
gammaNode3 = GammaNode(3,[x0,y0]+[dx,dy])
gammaNode4 = GammaNode(4,[x0,y0]+[0.0,dy])
gammaNodes = [gammaNode1,gammaNode2,gammaNode3,gammaNode4]
gamma1 = Gamma(1,"gammaLine",[gammaNode1.id,gammaNode2.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaLine",[gammaNode2.id,gammaNode3.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma3 = Gamma(3,"gammaLine",[gammaNode3.id,gammaNode4.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma4 = Gamma(4,"gammaLine",[gammaNode4.id,gammaNode1.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2,gamma3,gamma4]

##############################################
# L-shape:
##############################################

x0 = 0.0
y0 = 0.0
dx1 = 0.5
dx2 = 0.5
dy1 = 0.5
dy2 = 0.5

gammaNode1 = GammaNode(1,[x0,y0] + [0.0,0.0])
gammaNode2 = GammaNode(2,[x0,y0] + [dx1+dx2,0.0])
gammaNode3 = GammaNode(3,[x0,y0] + [dx1+dx2,dy1])
gammaNode4 = GammaNode(4,[x0,y0] + [dx1,dy1])
gammaNode5 = GammaNode(5,[x0,y0] + [dx1,dy1+dy2])
gammaNode6 = GammaNode(6,[x0,y0] + [0.0,dy1+dy2])
gammaNodes = [gammaNode1,gammaNode2,gammaNode3,gammaNode4,gammaNode5,gammaNode6]
gamma1 = Gamma(1,"gammaLine",[gammaNode1.id,gammaNode2.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaLine",[gammaNode2.id,gammaNode3.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma3 = Gamma(3,"gammaLine",[gammaNode3.id,gammaNode4.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma4 = Gamma(4,"gammaLine",[gammaNode4.id,gammaNode5.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma5 = Gamma(5,"gammaLine",[gammaNode5.id,gammaNode6.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma6 = Gamma(6,"gammaLine",[gammaNode6.id,gammaNode1.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2,gamma3,gamma4,gamma5,gamma6]

##############################################
# L-pipe:
##############################################

x0 = 0.0
y0 = 0.0
dx = 0.5
dy = 0.5
radiusMax = 0.5
width = 0.25

gammaNode1 = GammaNode(1,[x0,y0] + [0.0,0.0])
gammaNode2 = GammaNode(2,[x0,y0] + [dx,0.0])
gammaNode3 = GammaNode(3,[x0,y0] + [dx+radiusMax,radiusMax])
gammaNode4 = GammaNode(4,[x0,y0] + [dx+radiusMax,dy+radiusMax])
gammaNode5 = GammaNode(5,[x0,y0] + [dx+radiusMax-width,dy+radiusMax])
gammaNode6 = GammaNode(6,[x0,y0] + [dx+radiusMax-width,radiusMax])
gammaNode7 = GammaNode(7,[x0,y0] + [dx,width])
gammaNode8 = GammaNode(8,[x0,y0] + [0.0,width])
gammaNodeCenter = GammaNode(9,[x0,y0] + [dx,radiusMax])
gammaNodes = [gammaNode1,gammaNode2,gammaNode3,gammaNode4,gammaNode5,gammaNode6,gammaNode7,gammaNode8,gammaNodeCenter]
gamma1 = Gamma(1,"gammaLine",[gammaNode1.id,gammaNode2.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaArc",[gammaNodeCenter.id,gammaNode2.id],[pi/2],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma3 = Gamma(3,"gammaLine",[gammaNode3.id,gammaNode4.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma4 = Gamma(4,"gammaLine",[gammaNode4.id,gammaNode5.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma5 = Gamma(5,"gammaLine",[gammaNode5.id,gammaNode6.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma6 = Gamma(6,"gammaArc",[gammaNodeCenter.id,gammaNode6.id],[-pi/2],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma7 = Gamma(7,"gammaLine",[gammaNode7.id,gammaNode8.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma8 = Gamma(8,"gammaLine",[gammaNode8.id,gammaNode1.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8]

##############################################
# U-shape:
##############################################

x0 = 0.0
y0 = 0.0
dx1 = 0.25
dx2 = 0.5
dx3 = 0.25
dy1 = 0.25
dy2 = 0.75

gammaNode1 = GammaNode(1,[x0,y0] + [0.0,0.0])
gammaNode2 = GammaNode(2,[x0,y0] + [dx1+dx2+dx3,0.0])
gammaNode3 = GammaNode(3,[x0,y0] + [dx1+dx2+dx3,dy1+dy2])
gammaNode4 = GammaNode(4,[x0,y0] + [dx1+dx2,dy1+dy2])
gammaNode5 = GammaNode(5,[x0,y0] + [dx1+dx2,dy1])
gammaNode6 = GammaNode(6,[x0,y0] + [dx1,dy1])
gammaNode7 = GammaNode(7,[x0,y0] + [dx1,dy1+dy2])
gammaNode8 = GammaNode(8,[x0,y0] + [0.0,dy1+dy2])
gammaNodes = [gammaNode1,gammaNode2,gammaNode3,gammaNode4,gammaNode5,gammaNode6,gammaNode7,gammaNode8]
gamma1 = Gamma(1,"gammaLine",[gammaNode1.id,gammaNode2.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaLine",[gammaNode2.id,gammaNode3.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma3 = Gamma(3,"gammaLine",[gammaNode3.id,gammaNode4.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma4 = Gamma(4,"gammaLine",[gammaNode4.id,gammaNode5.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma5 = Gamma(5,"gammaLine",[gammaNode5.id,gammaNode6.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma6 = Gamma(6,"gammaLine",[gammaNode6.id,gammaNode7.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma7 = Gamma(7,"gammaLine",[gammaNode7.id,gammaNode8.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma8 = Gamma(8,"gammaLine",[gammaNode8.id,gammaNode1.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8]

##############################################
# full circle:
##############################################

x0 = 0.0
y0 = 0.0
radius = 0.5

gammaNodeCenter = GammaNode(1,[x0,y0] + [0.0,0.0])
gammaNode1 = GammaNode(2,[x0,y0] + [radius,0.0])
gammaNodes = [gammaNodeCenter,gammaNode1]
gamma1 = Gamma(1,"gammaArc",[gammaNodeCenter.id,gammaNode1.id],[2*pi],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1]

##############################################
# circular sector:
##############################################

x0 = 0.0
y0 = 0.0
radius = 1.0
angle = pi/2

gammaNodeCenter = GammaNode(1,[x0,y0] + [0.0,0.0])
gammaNode1 = GammaNode(2,[x0,y0] + [radius,0.0])
gammaNode2 = GammaNode(3,rotate2DPoint([x0,y0],angle,[x0,y0] + [radius,0.0]))
gammaNodes = [gammaNodeCenter,gammaNode1,gammaNode2]
gamma1 = Gamma(1,"gammaLine",[gammaNodeCenter.id,gammaNode1.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaArc",[gammaNodeCenter.id,gammaNode1.id],[angle],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma3 = Gamma(3,"gammaLine",[gammaNode2.id,gammaNodeCenter.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2,gamma3]

##############################################
# turbine:
##############################################

x0 = 0.0
y0 = 0.0
dx1 = 0.5
dx2 = 0.5
width = 0.25
radius = 0.25

intersectionsDx = sqrt(radius^2 - (width/2)^2)
intersectionsAngle = pi - 2*asin((width/2)/radius)
gammaNode1 = GammaNode(1,[x0,y0] + [0.0,0.0])
gammaNode2 = GammaNode(2,[x0,y0] + [dx1-intersectionsDx,0.0])
gammaNode3 = GammaNode(3,[x0,y0] + [dx1+intersectionsDx,0.0])
gammaNode4 = GammaNode(4,[x0,y0] + [dx1+dx2,0.0])
gammaNode5 = GammaNode(5,[x0,y0] + [dx1+dx2,width])
gammaNode6 = GammaNode(6,[x0,y0] + [dx1+intersectionsDx,width])
gammaNode7 = GammaNode(7,[x0,y0] + [dx1-intersectionsDx,width])
gammaNode8 = GammaNode(8,[x0,y0] + [0.0,width])
gammaNodeCenter = GammaNode(9,[x0,y0] + [dx1,width/2])
gammaNodes = [gammaNode1,gammaNode2,gammaNode3,gammaNode4,gammaNode5,gammaNode6,gammaNode7,gammaNode8,gammaNodeCenter]
gamma1 = Gamma(1,"gammaLine",[gammaNode1.id,gammaNode2.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaArc",[gammaNodeCenter.id,gammaNode2.id],[intersectionsAngle],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma3 = Gamma(3,"gammaLine",[gammaNode3.id,gammaNode4.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma4 = Gamma(4,"gammaLine",[gammaNode4.id,gammaNode5.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma5 = Gamma(5,"gammaLine",[gammaNode5.id,gammaNode6.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma6 = Gamma(6,"gammaArc",[gammaNodeCenter.id,gammaNode6.id],[intersectionsAngle],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma7 = Gamma(7,"gammaLine",[gammaNode7.id,gammaNode8.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma8 = Gamma(8,"gammaLine",[gammaNode8.id,gammaNode1.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8]

##############################################
# H-shape:
##############################################

x0 = 0.0
y0 = 0.0
dx1 = 0.45
dx2 = 0.1
dx3 = 0.45
dy1 = 0.4
dy2 = 0.2
dy3 = 0.4

gammaNode1 = GammaNode(1,[x0,y0] + [0.0,0.0])
gammaNode2 = GammaNode(2,[x0,y0] + [dx1,0.0])
gammaNode3 = GammaNode(3,[x0,y0] + [dx1,dy1])
gammaNode4 = GammaNode(4,[x0,y0] + [dx1+dx2,dy1])
gammaNode5 = GammaNode(5,[x0,y0] + [dx1+dx2,0.0])
gammaNode6 = GammaNode(6,[x0,y0] + [dx1+dx2+dx3,0.0])
gammaNode7 = GammaNode(7,[x0,y0] + [dx1+dx2+dx3,dy1+dy2+dy3])
gammaNode8 = GammaNode(8,[x0,y0] + [dx1+dx2,dy1+dy2+dy3])
gammaNode9 = GammaNode(9,[x0,y0] + [dx1+dx2,dy1+dy2])
gammaNode10 = GammaNode(10,[x0,y0] + [dx1,dy1+dy2])
gammaNode11 = GammaNode(11,[x0,y0] + [dx1,dy1+dy2+dy3])
gammaNode12 = GammaNode(12,[x0,y0] + [0.0,dy1+dy2+dy3])
gammaNodes = [gammaNode1,gammaNode2,gammaNode3,gammaNode4,gammaNode5,gammaNode6,gammaNode7,gammaNode8,gammaNode9,gammaNode10,gammaNode11,gammaNode12]
gamma1 = Gamma(1,"gammaLine",[gammaNode1.id,gammaNode2.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaLine",[gammaNode2.id,gammaNode3.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma3 = Gamma(3,"gammaLine",[gammaNode3.id,gammaNode4.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma4 = Gamma(4,"gammaLine",[gammaNode4.id,gammaNode5.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma5 = Gamma(5,"gammaLine",[gammaNode5.id,gammaNode6.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma6 = Gamma(6,"gammaLine",[gammaNode6.id,gammaNode7.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma7 = Gamma(7,"gammaLine",[gammaNode7.id,gammaNode8.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma8 = Gamma(8,"gammaLine",[gammaNode8.id,gammaNode9.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma9 = Gamma(9,"gammaLine",[gammaNode9.id,gammaNode10.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma10 = Gamma(10,"gammaLine",[gammaNode10.id,gammaNode11.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma11 = Gamma(11,"gammaLine",[gammaNode11.id,gammaNode12.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma12 = Gamma(12,"gammaLine",[gammaNode12.id,gammaNode1.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8,gamma9,gamma10,gamma11,gamma12]

##############################################
# loudspeaker:
##############################################

x0 = 0.0
y0 = 0.0
focus = 0.5
dx = 1.0

angle = 2*atan(sqrt(dx))
dy = 4*focus*sqrt(dx)
gammaNodeVertex = GammaNode(1,[x0,y0])
gammaNodeFocus = GammaNode(2,[x0,y0] + [focus,0.0])
gammaNode1 = GammaNode(3,[x0,y0] + [dx,-dy/2])
gammaNode2 = GammaNode(4,[x0,y0] + [dx,dy/2])
gammaNodes = [gammaNodeVertex,gammaNodeFocus,gammaNode1,gammaNode2]
gamma1 = Gamma(1,"gammaParabolaArc",[gammaNodeFocus.id,gammaNodeVertex.id],[-angle,angle],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaLine",[gammaNode1.id,gammaNode2.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2]

##############################################
# lens:
##############################################

x0 = 0.0
y0 = 0.0
focus = 1.0
dx = 0.5

angle = 2*atan(sqrt(dx))
dy = 4*focus*sqrt(dx)
gammaNodeVertexLeft = GammaNode(1,[x0,y0])
gammaNodeFocusLeft = GammaNode(2,[x0,y0] + [focus,0.0])
gammaNodeVertexRight = GammaNode(3,[x0,y0] + [2*dx,0.0])
gammaNodeFocusRight = GammaNode(4,[x0,y0] + [2*dx-focus,0.0])
gammaNodes = [gammaNodeVertexLeft,gammaNodeFocusLeft,gammaNodeVertexRight,gammaNodeFocusRight]
gamma1 = Gamma(1,"gammaParabolaArc",[gammaNodeFocusLeft.id,gammaNodeVertexLeft.id],[-angle,angle],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaParabolaArc",[gammaNodeFocusRight.id,gammaNodeVertexRight.id],[-angle,angle],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2]

##############################################
# tracks for pcb (printed circuit board):
##############################################

function computeGammaNodesCoordsOfAxesParallelTrack(coordsStart::Array{Float64,1},thickness::Float64,directions::Array{ASCIIString,1},distances::Array{Float64,1})
	
	gammaNodesCoordsCenterOfTrack = Array{Array{Float64,1},1}(length(directions)+1)
	gammaNodesCoordsCenterOfTrack[1] = copy(coordsStart)
	for i = 1:1:length(directions)
		
		if directions[i] == "north"
			directionVectorCenter = [0,1]
		elseif directions[i] == "east"
			directionVectorCenter = [1,0]
		elseif directions[i] == "south"
			directionVectorCenter = [0,-1]
		elseif directions[i] == "west"
			directionVectorCenter = [-1,0]
		end
		
		gammaNodesCoordsCenterOfTrack[i+1] = gammaNodesCoordsCenterOfTrack[i] + distances[i]*directionVectorCenter
	end
	
	gammaNodesCoordsLefthandBoundaryOfTrack = Array{Array{Float64,1},1}(length(gammaNodesCoordsCenterOfTrack))
	gammaNodesCoordsRighthandBoundaryOfTrack = Array{Array{Float64,1},1}(length(gammaNodesCoordsCenterOfTrack))
	for i = 1:1:length(gammaNodesCoordsCenterOfTrack)
		
		if i == 1
			
			if directions[1] == "north"
				directionVectorLeft = [-1,0]
				directionVectorRight = [1,0]
			elseif directions[1] == "east"
				directionVectorLeft = [0,1]
				directionVectorRight = [0,-1]
			elseif directions[1] == "south"
				directionVectorLeft = [1,0]
				directionVectorRight = [-1,0]
			elseif directions[1] == "west"
				directionVectorLeft = [0,-1]
				directionVectorRight = [0,1]
			end
			
		elseif i == length(gammaNodesCoordsCenterOfTrack)
			
			if directions[end] == "north"
				directionVectorLeft = [-1,0]
				directionVectorRight = [1,0]
			elseif directions[end] == "east"
				directionVectorLeft = [0,1]
				directionVectorRight = [0,-1]
			elseif directions[end] == "south"
				directionVectorLeft = [1,0]
				directionVectorRight = [-1,0]
			elseif directions[end] == "west"
				directionVectorLeft = [0,-1]
				directionVectorRight = [0,1]
			end
			
		else
			
			if directions[i-1] == "north"
				if directions[i] == "north"
					directionVectorLeft = [-1,0]
					directionVectorRight = [1,0]
				elseif directions[i] == "east"
					directionVectorLeft = [-1,1]
					directionVectorRight = [1,-1]
				elseif directions[i] == "south"
					error("Inconsistent directions given.")
				elseif directions[i] == "west"
					directionVectorLeft = [-1,-1]
					directionVectorRight = [1,1]
				end
			elseif directions[i-1] == "east"
				if directions[i] == "north"
					directionVectorLeft = [-1,1]
					directionVectorRight = [1,-1]
				elseif directions[i] == "east"
					directionVectorLeft = [0,1]
					directionVectorRight = [0,-1]
				elseif directions[i] == "south"
					directionVectorLeft = [1,1]
					directionVectorRight = [-1,-1]
				elseif directions[i] == "west"
					error("Inconsistent directions given.")
				end
			elseif directions[i-1] == "south"
				if directions[i] == "north"
					error("Inconsistent directions given.")
				elseif directions[i] == "east"
					directionVectorLeft = [1,1]
					directionVectorRight = [-1,-1]
				elseif directions[i] == "south"
					directionVectorLeft = [1,0]
					directionVectorRight = [-1,0]
				elseif directions[i] == "west"
					directionVectorLeft = [1,-1]
					directionVectorRight = [-1,1]
				end
			elseif directions[i-1] == "west"
				if directions[i] == "north"
					directionVectorLeft = [-1,-1]
					directionVectorRight = [1,1]
				elseif directions[i] == "east"
					error("Inconsistent directions given.")
				elseif directions[i] == "south"
					directionVectorLeft = [1,-1]
					directionVectorRight = [-1,1]
				elseif directions[i] == "west"
					directionVectorLeft = [0,-1]
					directionVectorRight = [0,1]
				end
			end
			
		end
		
		gammaNodesCoordsLefthandBoundaryOfTrack[i] = gammaNodesCoordsCenterOfTrack[i] + (1/2)*thickness*directionVectorLeft
		gammaNodesCoordsRighthandBoundaryOfTrack[i] = gammaNodesCoordsCenterOfTrack[i] + (1/2)*thickness*directionVectorRight
	end
	
	return gammaNodesCoordsCenterOfTrack,gammaNodesCoordsLefthandBoundaryOfTrack,gammaNodesCoordsRighthandBoundaryOfTrack
end

coordsStart = [0.0,0.0]
thickness = 0.1
directions = ["east","north","west","south"]
distances = [1.0,1.0,1.0,0.5]

gammaNodesCoordsCenterOfTrack,gammaNodesCoordsLefthandBoundaryOfTrack,gammaNodesCoordsRighthandBoundaryOfTrack = 
	computeGammaNodesCoordsOfAxesParallelTrack(coordsStart,thickness,directions,distances)

gammaNodes = Array{GammaNode,1}(length(gammaNodesCoordsRighthandBoundaryOfTrack)+length(gammaNodesCoordsLefthandBoundaryOfTrack))
for i = 1:1:length(gammaNodesCoordsRighthandBoundaryOfTrack)
	gammaNodes[i] = GammaNode(i,gammaNodesCoordsRighthandBoundaryOfTrack[i])
end
for i = 1:1:length(gammaNodesCoordsLefthandBoundaryOfTrack)
	gammaNodes[i+length(gammaNodesCoordsRighthandBoundaryOfTrack)] = GammaNode(i+length(gammaNodesCoordsRighthandBoundaryOfTrack),gammaNodesCoordsLefthandBoundaryOfTrack[length(gammaNodesCoordsLefthandBoundaryOfTrack)+1-i])
end
gammas = Array{Gamma,1}(length(gammaNodes))
for i = 1:1:length(gammaNodes)
	gammas[i] = Gamma(i,"gammaLine",[i,moduloIndex(i+1,length(gammaNodes))],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
end

##############################################

