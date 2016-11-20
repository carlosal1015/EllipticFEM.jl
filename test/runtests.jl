
using Base.Test

include("./Elliptic/runElliptic.jl")

# testing equationData:
@test equationData.equationType == "elliptic"
@test !isempty(equationData.subdomains)
@test !isempty(equationData.holes)

# testing geoData:
@test !isempty(geoData.geoNodes)
@test !isempty(geoData.geoEdges)
@test !isempty(geoData.geoEdgesLoops)

# testing meshData:
@test !isempty(meshData.nodes)
@test !isempty(meshData.edges)
@test !isempty(meshData.elements)

# testing lseData:
@test !isempty(lseData.M)
@test !isempty(lseData.r)
@test !isempty(lseData.solutions[1])

