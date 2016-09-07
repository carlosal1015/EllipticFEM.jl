
##############################################

equationType = "hyperbolic"
tMax = 10.0

##############################################

waveLength = 1.0
waveSpeed = 0.5
frequency = waveSpeed/waveLength

##############################################
# domain:
##############################################

x0 = 0.0
y0 = 0.0
focus = 0.5
dx = 5.0

angle = 2*atan(sqrt(dx))
dy = 4*focus*sqrt(dx)
gammaNodeVertex = GammaNode(1,[x0,y0])
gammaNodeFocus = GammaNode(2,[x0,y0] + [focus,0.0])
gammaNode1 = GammaNode(3,[x0,y0] + [dx,-dy/2])
gammaNode2 = GammaNode(4,[x0,y0] + [dx,dy/2])
gammaNodes = [gammaNodeVertex,gammaNodeFocus,gammaNode1,gammaNode2]
gamma1 = Gamma(1,"gammaParabolaArc",[gammaNodeFocus.id,gammaNodeVertex.id],[-angle,angle],50,()->[],"Dirichlet",(s,t)->0.0,true,true)
gamma2 = Gamma(2,"gammaLine",[gammaNode1.id,gammaNode2.id],[],1,()->[],"Dirichlet",(s,t)->0.0,true,true)
gammas = [gamma1,gamma2]

uInitial = (x,y,t) -> 0.0
uTInitial = (x,y,t) -> 0.0

A11 = (x,y,t) -> waveSpeed^2
A21 = (x,y,t) -> 0.0
A12 = (x,y,t) -> 0.0
A22 = (x,y,t) -> waveSpeed^2
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

isIdentityA = false
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
# hole 1:
##############################################

x0 = 0.5
y0 = 0.0
radius = 0.1

gammaNodeCenter = GammaNode(1,[x0,y0] + [0.0,0.0])
gammaNode1 = GammaNode(2,[x0,y0] + [radius,0.0])
gammaNodes = [gammaNodeCenter,gammaNode1]
gamma1 = Gamma(1,"gammaArc",[gammaNodeCenter.id,gammaNode1.id],[2*pi],10,()->[],"Dirichlet",(s,t)->(1/2)*(1 + sin(2*pi*frequency*t)),false,false)
gammas = [gamma1]

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

hole1 = Domain(1,"hole",gammaNodes,gammas,[],
	uInitial,uTInitial,A,b,c,f,fDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# holes:
##############################################

holes = [hole1]

##############################################
# subdomains:
##############################################

subdomains = []

##############################################

equationData = EquationData(equationType,tMax,
	domain,holes,subdomains)

##############################################



