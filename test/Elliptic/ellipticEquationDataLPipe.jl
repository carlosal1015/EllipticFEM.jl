
##############################################

equationType = "elliptic"
tMax = 0.0

##############################################
# domain:
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
gamma4 = Gamma(4,"gammaLine",[gammaNode4.id,gammaNode5.id],[],1,()->[],"Neumann",(s,t)->0.0,true,true)
gamma5 = Gamma(5,"gammaLine",[gammaNode5.id,gammaNode6.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma6 = Gamma(6,"gammaArc",[gammaNodeCenter.id,gammaNode6.id],[-pi/2],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma7 = Gamma(7,"gammaLine",[gammaNode7.id,gammaNode8.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma8 = Gamma(8,"gammaLine",[gammaNode8.id,gammaNode1.id],[],1,()->[],"Neumann",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2,gamma3,gamma4,gamma5,gamma6,gamma7,gamma8]

uInitial = (x,y) -> 0.0
uTInitial = (x,y) -> 0.0

A11 = (x,y) -> 1.0
A21 = (x,y) -> 0.0
A12 = (x,y) -> 0.0
A22 = (x,y) -> 1.0
A = hcat([A11,A21],[A12,A22])
bStrength = 100.0
b1 = (x,y) -> (valueIsInInterval(x,[x0,x0+dx]) ? bStrength : 0.0) + 
	(valueIsInInterval(x,[x0+dx,x0+dx+radiusMax])&valueIsInInterval(y,[y0,y0+radiusMax]) ? -bStrength*(y-(y0+radiusMax))/norm([x,y]-([x0,y0]+[dx,radiusMax])) : 0.0) + 
	(valueIsInInterval(y,[y0+radiusMax,y0+radiusMax+dy]) ? 0.0 : 0.0)
b2 = (x,y) -> (valueIsInInterval(x,[x0,x0+dx]) ? 0.0 : 0.0) + 
	(valueIsInInterval(x,[x0+dx,x0+dx+radiusMax])&valueIsInInterval(y,[y0,y0+radiusMax]) ? bStrength*(x-(x0+dx))/norm([x,y]-([x0,y0]+[dx,radiusMax])) : 0.0) + 
	(valueIsInInterval(y,[y0+radiusMax,y0+radiusMax+dy]) ? bStrength : 0.0)
b = [b1,b2]
c = (x,y) -> 0.0
f = (x,y) -> 0.0
fDiv1 = (x,y) -> 0.0
fDiv2 = (x,y) -> 0.0
fDiv = [fDiv1,fDiv2]

isZeroUInitial = true
isZeroUTInitial = true

isIdentityA = true
isZeroB = false
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
# hole 1:
##############################################

x0Hole = 0.85
y0Hole = 0.5
dxHole = 0.05
dyHole = 0.3

gammaNode1 = GammaNode(1,[x0Hole,y0Hole]+[0.0,0.0])
gammaNode2 = GammaNode(2,[x0Hole,y0Hole]+[dxHole,0.0])
gammaNode3 = GammaNode(3,[x0Hole,y0Hole]+[dxHole,dyHole])
gammaNode4 = GammaNode(4,[x0Hole,y0Hole]+[0.0,dyHole])
gammaNodes = [gammaNode1,gammaNode2,gammaNode3,gammaNode4]
gamma1 = Gamma(1,"gammaLine",[gammaNode1.id,gammaNode2.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaLine",[gammaNode2.id,gammaNode3.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma3 = Gamma(3,"gammaLine",[gammaNode3.id,gammaNode4.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma4 = Gamma(4,"gammaLine",[gammaNode4.id,gammaNode1.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2,gamma3,gamma4]

uInitial = (x,y) -> 0.0
uTInitial = (x,y) -> 0.0

A11 = (x,y) -> 1.0
A21 = (x,y) -> 0.0
A12 = (x,y) -> 0.0
A22 = (x,y) -> 1.0
A = hcat([A11,A21],[A12,A22])
b1 = (x,y) -> 0.0
b2 = (x,y) -> 0.0
b = [b1,b2]
c = (x,y) -> 0.0
f = (x,y) -> 0.0
fDiv1 = (x,y) -> 0.0
fDiv2 = (x,y) -> 0.0
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

hole1 = Domain(1,"hole",gammaNodes,gammas,[],
	uInitial,uTInitial,A,b,c,f,fDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# holes:
##############################################

holes = [hole1]

##############################################
# subdomain 1:
##############################################

x0Subdomain = 0.1
y0Subdomain = 0.125
radius = 0.05

gammaNodeCenter = GammaNode(1,[x0Subdomain,y0Subdomain] + [0.0,0.0])
gammaNode1 = GammaNode(2,[x0Subdomain,y0Subdomain] + [radius,0.0])
gammaNodes = [gammaNodeCenter,gammaNode1]
gamma1 = Gamma(1,"gammaArc",[gammaNodeCenter.id,gammaNode1.id],[2*pi],10,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1]

uInitial = (x,y) -> 0.0
uTInitial = (x,y) -> 0.0

A11 = (x,y) -> 1.0
A21 = (x,y) -> 0.0
A12 = (x,y) -> 0.0
A22 = (x,y) -> 1.0
A = hcat([A11,A21],[A12,A22])
b1 = (x,y) -> 0.0
b2 = (x,y) -> 0.0
b = [b1,b2]
c = (x,y) -> 0.0
f = (x,y) -> 1.0
fDiv1 = (x,y) -> 0.0
fDiv2 = (x,y) -> 0.0
fDiv = [fDiv1,fDiv2]

isZeroUInitial = true
isZeroUTInitial = true

isIdentityA = true
isZeroB = true
isZeroC = true
isZeroF = false
isZeroFDiv = true

isConstantInTimeA = true
isConstantInTimeB = true
isConstantInTimeC = true
isConstantInTimeF = true
isConstantInTimeFDiv = true

subdomain1 = Domain(1,"subdomain",gammaNodes,gammas,[],
	uInitial,uTInitial,A,b,c,f,fDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# subdomains:
##############################################

subdomains = [subdomain1]

##############################################

equationData = EquationData(equationType,tMax,
	domain,holes,subdomains)

##############################################



