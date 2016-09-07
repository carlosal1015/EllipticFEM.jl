
##############################################

function assembleCentersAndRadiiDeterministic()
	
	xs = [0.0,0.25,0.5,0.75,1.0]
	ys = [0.5,0.5,0.5,0.5,0.5]
	radii = [0.1,0.1,0.1,0.1,0.1]
	
	centersAndRadii = hcat(xs,ys,radii)'
	
	return centersAndRadii
end

function assembleCentersAndRadiiRandom(intendedNumberOfInclusions::Int64)
	
	# TODO: Set threshold properly:
	threshold1 = 0.000001
	threshold2 = 0.000001
	
	xMin = 0
	xMax = 1
	yMin = 0
	yMax = 1
	radiusMin = 0.1
	radiusMax = 0.1
	
	xs = Array{Float64}(0)
	ys = Array{Float64}(0)
	radii = Array{Float64}(0)
	
	counterMax = 100
	for counter = 1:1:counterMax
		
		# Choose random center coordinates:
		x = xMin + rand()*(xMax-xMin)
		y = yMin + rand()*(yMax-yMin)
		r = radiusMin + rand()*(radiusMax-radiusMin)
		
		xsCurrent = x + (xMax-xMin)*[-1,0,1,-1,0,1,-1,0,1]
		ysCurrent = y + (yMax-yMin)*[-1,-1,-1,0,0,0,1,1,1]
		radiiCurrent = r*ones(9)
		isActiveCurrent = falses(9)
		for i = 1:1:9
			if valueIsInInterval(xsCurrent[i],[xMin,xMax] + radiiCurrent[i]*[-1,1] + threshold1*[1,-1])&
				valueIsInInterval(ysCurrent[i],[yMin,yMax] + radiiCurrent[i]*[-1,1] + threshold1*[1,-1])
				
				isActiveCurrent[i] = true
			end
		end
		
		# Check if inclusions overlap:
		doPush = true
		breakOuterLoop = false
		for i = 1:1:length(xs)
			for j = 1:1:length(xsCurrent)
				if isActiveCurrent[j]
					if norm([xs[i],ys[i]]-[xsCurrent[j],ysCurrent[j]]) <= radii[i]+radiiCurrent[j]+threshold2
						doPush = false
						breakOuterLoop = true
						break
					end
				end
			end
			if breakOuterLoop
				break
			end
		end
		
		# Add current inclusions:
		if doPush
			append!(xs,xsCurrent[isActiveCurrent])
			append!(ys,ysCurrent[isActiveCurrent])
			append!(radii,radiiCurrent[isActiveCurrent])
		end
		
		if length(xs) >= intendedNumberOfInclusions
			println("stopped at counter = "*string(counter))
			break
		end
		
		if counter == counterMax
			println("stopped at counterMax = "*string(counter))
		end
	end
	
	centersAndRadii = hcat(xs,ys,radii)'
	
	return centersAndRadii
end

##############################################

equationType = "elliptic"
tMax = 0.0

##############################################

# intendedNumberOfInclusions = 10
# centersAndRadii = assembleCentersAndRadiiRandom(intendedNumberOfInclusions)
centersAndRadii = assembleCentersAndRadiiDeterministic()

AValueDomain = 1.0
AValueInclusions = 10.0

##############################################
# domain:
##############################################

domainWidth = 1.0
domainHeight = 1.0

x0 = 0.0
y0 = 0.0
dx = 1.0
dy = 1.0

gammaNode1 = GammaNode(1,[0,0])
gammaNode2 = GammaNode(2,[domainWidth,0])
gammaNode3 = GammaNode(3,[domainWidth,domainHeight])
gammaNode4 = GammaNode(4,[0,domainHeight])
gammaNodes = [gammaNode1,gammaNode2,gammaNode3,gammaNode4]
gamma1 = Gamma(1,"gammaLine",[gammaNode1.id,gammaNode2.id],[],1,()->[],"periodic",(s)->0.0,true,true)
gamma2 = Gamma(2,"gammaLine",[gammaNode2.id,gammaNode3.id],[],1,()->[],"periodic",(s)->0.0,true,true)
gamma3 = Gamma(3,"gammaLine",[gammaNode3.id,gammaNode4.id],[],1,()->[],"periodic",(s)->0.0,true,true)
gamma4 = Gamma(4,"gammaLine",[gammaNode4.id,gammaNode1.id],[],1,()->[],"periodic",(s)->0.0,true,true)
gammas = [gamma1,gamma2,gamma3,gamma4]

periodicGammasPairs =
	[PeriodicGammasPair(1,[1,3],false),
	PeriodicGammasPair(2,[2,4],false)]

uInitial = (x,y) -> 0.0
uTInitial = (x,y) -> 0.0

domainA11 = (x,y) -> AValueDomain
domainA21 = (x,y) -> 0.0
domainA12 = (x,y) -> 0.0
domainA22 = (x,y) -> AValueDomain
domainA = hcat([domainA11,domainA21],[domainA12,domainA22])
b1 = (x,y) -> 0.0
b2 = (x,y) -> 0.0
b = [b1,b2]
c = (x,y) -> 0.0
f = (x,y) -> 0.0
domainFDiv1 = (x,y) -> 0.0
domainFDiv2 = (x,y) -> 0.0
domainFDiv = [domainFDiv1,domainFDiv2]

isZeroUInitial = true
isZeroUTInitial = true

isIdentityA = false
isZeroB = true
isZeroC = true
isZeroF = true
isZeroFDiv = false

isConstantInTimeA = true
isConstantInTimeB = true
isConstantInTimeC = true
isConstantInTimeF = true
isConstantInTimeFDiv = true

domain = Domain(1,"domain",gammaNodes,gammas,periodicGammasPairs,
	uInitial,uTInitial,domainA,b,c,f,domainFDiv,
	isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
	isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)

##############################################
# holes:
##############################################

holes = []

##############################################
# subdomains:
##############################################

subdomains = Array{Domain}(size(centersAndRadii,2))
for i = 1:1:size(centersAndRadii,2)
	
	x = centersAndRadii[1,i]
	y = centersAndRadii[2,i]
	r = centersAndRadii[3,i]
	
	numberOfPoints = 16

	gammaNodeCenter = GammaNode(1,[x,y])
	gammaNode1 = GammaNode(2,[x+r,y])
	gammaNodes = [gammaNodeCenter,gammaNode1]
	gamma1 = Gamma(1,"gammaArc",[gammaNodeCenter.id,gammaNode1.id],[2*pi],numberOfPoints,()->[],"Dirichlet",(s)->0.0,true,true)
	gammas = [gamma1]
	
	uInitial = (x,y) -> 0.0
	uTInitial = (x,y) -> 0.0

	subdomainsA11 = (x,y) -> AValueInclusions
	subdomainsA21 = (x,y) -> 0.0
	subdomainsA12 = (x,y) -> 0.0
	subdomainsA22 = (x,y) -> AValueInclusions
	subdomainsA = hcat([subdomainsA11,subdomainsA21],[subdomainsA12,subdomainsA22])
	b1 = (x,y) -> 0.0
	b2 = (x,y) -> 0.0
	b = [b1,b2]
	c = (x,y) -> 0.0
	f = (x,y) -> 0.0
	subdomainsFDiv1 = (x,y) -> 0.0
	subdomainsFDiv2 = (x,y) -> 0.0
	subdomainsFDiv = [subdomainsFDiv1,subdomainsFDiv2]

	isZeroUInitial = true
	isZeroUTInitial = true

	isIdentityA = false
	isZeroB = true
	isZeroC = true
	isZeroF = true
	isZeroFDiv = false

	isConstantInTimeA = true
	isConstantInTimeB = true
	isConstantInTimeC = true
	isConstantInTimeF = true
	isConstantInTimeFDiv = true

	subdomains[i] = Domain(i,"subdomain",gammaNodes,gammas,[],
		uInitial,uTInitial,subdomainsA,b,c,f,subdomainsFDiv,
		isZeroUInitial,isZeroUTInitial,isIdentityA,isZeroB,isZeroC,isZeroF,isZeroFDiv,
		isConstantInTimeA,isConstantInTimeB,isConstantInTimeC,isConstantInTimeF,isConstantInTimeFDiv)
	
end

##############################################

equationData = EquationData(equationType,tMax,
	domain,holes,subdomains)

##############################################


