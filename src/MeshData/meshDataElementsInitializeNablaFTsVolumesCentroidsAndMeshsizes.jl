
function meshDataElementsInitializeNablaFTsVolumesCentroidsAndMeshsizes(meshData::MeshData)
	
	# Edges volumes:
	node1sIds = collect(meshData.edges[1,:])
	node2sIds = collect(meshData.edges[2,:])
	node1sXs = collect(meshData.nodes[1,node1sIds])
	node1sYs = collect(meshData.nodes[2,node1sIds])
	node2sXs = collect(meshData.nodes[1,node2sIds])
	node2sYs = collect(meshData.nodes[2,node2sIds])

	meshData.edgesVolumes = sqrt((node2sXs-node1sXs).^2 + (node2sYs-node1sYs).^2)

	# Elements volumes, \nabla FT, centroids and meshsizes:
	node1sIds = collect(meshData.elements[1,:])
	node2sIds = collect(meshData.elements[2,:])
	node3sIds = collect(meshData.elements[3,:])
	node1sXs = collect(meshData.nodes[1,node1sIds])
	node1sYs = collect(meshData.nodes[2,node1sIds])
	node2sXs = collect(meshData.nodes[1,node2sIds])
	node2sYs = collect(meshData.nodes[2,node2sIds])
	node3sXs = collect(meshData.nodes[1,node3sIds])
	node3sYs = collect(meshData.nodes[2,node3sIds])
	
	meshData.elementsDetNablaFTs = (node2sYs-node1sYs).*(node1sXs-node3sXs) - (node2sXs-node1sXs).*(node1sYs-node3sYs)
	meshData.elementsNablaFTInv11s = (node3sYs-node1sYs)./(meshData.elementsDetNablaFTs)
	meshData.elementsNablaFTInv21s = -(node2sYs-node1sYs)./(meshData.elementsDetNablaFTs)
	meshData.elementsNablaFTInv12s = -(node3sXs-node1sXs)./(meshData.elementsDetNablaFTs)
	meshData.elementsNablaFTInv22s = (node2sXs-node1sXs)./(meshData.elementsDetNablaFTs)
	meshData.elementsVolumes = (1/2)*abs(meshData.elementsDetNablaFTs)
	meshData.elementsCentroids = vcat(1/3*(meshData.nodes[1,node1sIds] + meshData.nodes[1,node2sIds] + meshData.nodes[1,node3sIds]),
		1/3*(meshData.nodes[2,node1sIds] + meshData.nodes[2,node2sIds] + meshData.nodes[2,node3sIds]))
	meshData.elementsMeshsizes = collect(maximum(hcat(
		sqrt((node2sXs-node1sXs).^2 + (node2sYs-node1sYs).^2),
		sqrt((node3sXs-node2sXs).^2 + (node3sYs-node2sYs).^2),
		sqrt((node1sXs-node3sXs).^2 + (node1sYs-node3sYs).^2)),
		2))

	return nothing
end
