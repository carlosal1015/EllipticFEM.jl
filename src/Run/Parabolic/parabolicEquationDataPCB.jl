
##############################################

equationType = "parabolic"
tMax = 10.0

##############################################

conductivityDomain = 0.00025
conductivityTracks = 1.0
fanTOn = 5.0
fanStrength = 0.1

##############################################
# domain:
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

uInitial = (x,y,t) -> 0.0
uTInitial = (x,y,t) -> 0.0

A11 = (x,y,t) -> conductivityDomain
A21 = (x,y,t) -> 0.0
A12 = (x,y,t) -> 0.0
A22 = (x,y,t) -> conductivityDomain
A = hcat([A11,A21],[A12,A22])
b1 = (x,y,t) -> (t <= fanTOn ? 0.0 : fanStrength)
b2 = (x,y,t) -> 0.0
b = [b1,b2]
c = (x,y,t) -> 0.0
f = (x,y,t) -> 0.0
fDiv1 = (x,y,t) -> 0.0
fDiv2 = (x,y,t) -> 0.0
fDiv = [fDiv1,fDiv2]

isZeroUInitial = true
isZeroUTInitial = true

isIdentityA = false
isZeroB = false
isZeroC = true
isZeroF = true
isZeroFDiv = true

isConstantInTimeA = true
isConstantInTimeB = false
isConstantInTimeC = true
isConstantInTimeF = true
isConstantInTimeFDiv = true

domain = Domain(1,"domain",gammaNodes,gammas,[],
	uInitial,uTInitial,A,b,c,f,fDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# holes:
##############################################

holes = []

##############################################
# subdomain 1:
##############################################

x0 = 0.4
y0 = 0.4
dx = 0.2
dy = 0.2

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

uInitial = (x,y,t) -> 0.0
uTInitial = (x,y,t) -> 0.0

A11 = (x,y,t) -> 1.0
A21 = (x,y,t) -> 0.0
A12 = (x,y,t) -> 0.0
A22 = (x,y,t) -> 1.0
A = hcat([A11,A21],[A12,A22])
b1 = (x,y,t) -> (t <= fanTOn ? 0.0 : fanStrength)
b2 = (x,y,t) -> 0.0
b = [b1,b2]
c = (x,y,t) -> 0.0
f = (x,y,t) -> 1.0
fDiv1 = (x,y,t) -> 0.0
fDiv2 = (x,y,t) -> 0.0
fDiv = [fDiv1,fDiv2]

isZeroUInitial = true
isZeroUTInitial = true

isIdentityA = true
isZeroB = false
isZeroC = true
isZeroF = false
isZeroFDiv = true

isConstantInTimeA = true
isConstantInTimeB = false
isConstantInTimeC = true
isConstantInTimeF = true
isConstantInTimeFDiv = true

subdomain1 = Domain(1,"subdomain",gammaNodes,gammas,[],
	uInitial,uTInitial,A,b,c,f,fDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# subdomains 2-5:
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

##############################################
# subdomain 2:
##############################################

epsilon = 0.001
coordsStart = [0.6+epsilon,0.5]
thickness = 0.05
directions = ["east","north","west","south"]
distances = [0.2,0.3,0.2,0.1]

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

uInitial = (x,y,t) -> 0.0
uTInitial = (x,y,t) -> 0.0

A11 = (x,y,t) -> conductivityTracks
A21 = (x,y,t) -> 0.0
A12 = (x,y,t) -> 0.0
A22 = (x,y,t) -> conductivityTracks
A = hcat([A11,A21],[A12,A22])
b1 = (x,y,t) -> (t <= fanTOn ? 0.0 : fanStrength)
b2 = (x,y,t) -> 0.0
b = [b1,b2]
c = (x,y,t) -> 0.0
f = (x,y,t) -> 0.0
fDiv1 = (x,y,t) -> 0.0
fDiv2 = (x,y,t) -> 0.0
fDiv = [fDiv1,fDiv2]

isZeroUInitial = true
isZeroUTInitial = true

isIdentityA = false
isZeroB = false
isZeroC = true
isZeroF = true
isZeroFDiv = true

isConstantInTimeA = true
isConstantInTimeB = false
isConstantInTimeC = true
isConstantInTimeF = true
isConstantInTimeFDiv = true

subdomain2 = Domain(2,"subdomain",gammaNodes,gammas,[],
	uInitial,uTInitial,A,b,c,f,fDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# subdomain 3:
##############################################

coordsStart = [0.5,0.6+epsilon]
thickness = 0.05
directions = ["north","east","south","west"]
distances = [0.3,0.4,0.6,0.2]

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

uInitial = (x,y,t) -> 0.0
uTInitial = (x,y,t) -> 0.0

A11 = (x,y,t) -> conductivityTracks
A21 = (x,y,t) -> 0.0
A12 = (x,y,t) -> 0.0
A22 = (x,y,t) -> conductivityTracks
A = hcat([A11,A21],[A12,A22])
b1 = (x,y,t) -> (t <= fanTOn ? 0.0 : fanStrength)
b2 = (x,y,t) -> 0.0
b = [b1,b2]
c = (x,y,t) -> 0.0
f = (x,y,t) -> 0.0
fDiv1 = (x,y,t) -> 0.0
fDiv2 = (x,y,t) -> 0.0
fDiv = [fDiv1,fDiv2]

isZeroUInitial = true
isZeroUTInitial = true

isIdentityA = false
isZeroB = false
isZeroC = true
isZeroF = true
isZeroFDiv = true

isConstantInTimeA = true
isConstantInTimeB = false
isConstantInTimeC = true
isConstantInTimeF = true
isConstantInTimeFDiv = true

subdomain3 = Domain(3,"subdomain",gammaNodes,gammas,[],
	uInitial,uTInitial,A,b,c,f,fDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# subdomain 4:
##############################################

coordsStart = [0.4-epsilon,0.5]
thickness = 0.05
directions = ["west","north","west","south","east"]
distances = [0.1,0.4,0.2,0.8,0.85]

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

uInitial = (x,y,t) -> 0.0
uTInitial = (x,y,t) -> 0.0

A11 = (x,y,t) -> conductivityTracks
A21 = (x,y,t) -> 0.0
A12 = (x,y,t) -> 0.0
A22 = (x,y,t) -> conductivityTracks
A = hcat([A11,A21],[A12,A22])
b1 = (x,y,t) -> (t <= fanTOn ? 0.0 : fanStrength)
b2 = (x,y,t) -> 0.0
b = [b1,b2]
c = (x,y,t) -> 0.0
f = (x,y,t) -> 0.0
fDiv1 = (x,y,t) -> 0.0
fDiv2 = (x,y,t) -> 0.0
fDiv = [fDiv1,fDiv2]

isZeroUInitial = true
isZeroUTInitial = true

isIdentityA = false
isZeroB = false
isZeroC = true
isZeroF = true
isZeroFDiv = true

isConstantInTimeA = true
isConstantInTimeB = false
isConstantInTimeC = true
isConstantInTimeF = true
isConstantInTimeFDiv = true

subdomain4 = Domain(4,"subdomain",gammaNodes,gammas,[],
	uInitial,uTInitial,A,b,c,f,fDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# subdomain 5:
##############################################

coordsStart = [0.5,0.4-epsilon]
thickness = 0.05
directions = ["south","west","south","east"]
distances = [0.1,0.3,0.1,0.725]

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

uInitial = (x,y,t) -> 0.0
uTInitial = (x,y,t) -> 0.0

A11 = (x,y,t) -> conductivityTracks
A21 = (x,y,t) -> 0.0
A12 = (x,y,t) -> 0.0
A22 = (x,y,t) -> conductivityTracks
A = hcat([A11,A21],[A12,A22])
b1 = (x,y,t) -> (t <= fanTOn ? 0.0 : fanStrength)
b2 = (x,y,t) -> 0.0
b = [b1,b2]
c = (x,y,t) -> 0.0
f = (x,y,t) -> 0.0
fDiv1 = (x,y,t) -> 0.0
fDiv2 = (x,y,t) -> 0.0
fDiv = [fDiv1,fDiv2]

isZeroUInitial = true
isZeroUTInitial = true

isIdentityA = false
isZeroB = false
isZeroC = true
isZeroF = true
isZeroFDiv = true

isConstantInTimeA = true
isConstantInTimeB = false
isConstantInTimeC = true
isConstantInTimeF = true
isConstantInTimeFDiv = true

subdomain5 = Domain(5,"subdomain",gammaNodes,gammas,[],
	uInitial,uTInitial,A,b,c,f,fDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# subdomains:
##############################################

subdomains = [subdomain1,subdomain2,subdomain3,subdomain4,subdomain5]

##############################################

equationData = EquationData(equationType,tMax,
	domain,holes,subdomains)

##############################################



