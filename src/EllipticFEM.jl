module EllipticFEM

using PyPlot

export Mesh,read_mesh,assemble,solve,plot_solution,calculate_norm

type Mesh
    #2xnrnodes array with [x;y] coords in each column
    nodes::Array{Float64, 2}

    #4xnredges array with [node1; node2; phys_prop; geom_prop] in each column
    edges::Array{Int64, 2}

    #5xnrelements array with [node1; node2; node3; phys_prop; geom_prop]
    elements::Array{Int64, 2}
end
 

#############################################################################
function read_mesh(path)

    in = open(path, "r")
    lines = readdlm(in)
    i = 1
    @assert lines[i] == "\$MeshFormat"
    i += 1
    s = lines[i]
    @assert s[1] == 2.2
    i += 1
    @assert lines[i] == "\$EndMeshFormat"
    i += 1

    if lines[i] == "\$PhysicalNames"
        i += 1
        nr_names = lines[i]
        i += nr_names + 2 #Skip names
    end

    #Save each node as type Node and collect nodes into nodes array
    @assert lines[i] == "\$Nodes"
    i += 1
    nr_nodes = lines[i]
    nodes = zeros(2, nr_nodes)
    @assert lines[i+1:i+nr_nodes,4]==vec(zeros(nr_nodes,1))
    nodes[:,:]=lines[i+1:i+nr_nodes,2:3]'
    i += nr_nodes+1
    @assert lines[i] == "\$EndNodes"
    i += 1

    #Save each element and edge
    @assert lines[i] == "\$Elements"
    i += 1
    nr_all_elements = lines[i]
    pos_edges = findin(lines[i+1:i+nr_all_elements,2],1)
    nr_edges = length(pos_edges)
    edges = Array{Int64}(4,nr_edges)
    pos_elements = findin(lines[i+1:i+nr_all_elements,2],2)
    nr_elements = length(pos_elements)
    elements = Array{Int64}(5,nr_elements)
    @assert lines[i+pos_edges,3] == vec(ones(nr_edges,1)*2)
    edges[:,:] = lines[i+pos_edges,[6,7,4,5]]'
    @assert lines[i+pos_elements,3] == vec(ones(nr_elements,1)*2)
    elements[:,:] = lines[i+pos_elements,[6,7,8,4,5]]'

    if nr_all_elements != nr_edges+nr_elements
            error("Unknown type of elements used! Please only use edges and triangular elements!")
    end

    i+=nr_all_elements+1
    @assert lines[i] == "\$EndElements"
    i += 1

    ## Read next section if available.
    if i <= size(lines,1)
    warn("There are lines in the mesh file that are not needed and cannot be taken into account.")
    end
    close(in)
    Mesh(nodes, edges, elements)
end

###############################################################################
function assemble(mesh::Mesh, A::Float64, f::Function, bddata::Array)
#assemble_mass creates the first-order mass matrix for the Poisson equation
# - \nabla \cdot (A \nabla u) = f with boundary conditions consisting of both Dirichlet and Neumann parts on a domain
#[-L/2,L/2]^2
    #data acquiry
    nq = size(mesh.nodes,2)
    nme = size(mesh.elements,2)

    ph_prop_11 = ones(Float64,2,nme)*A
   
    #uu = 2xnr_elems array of node2-node3 with order [x;y]
    #vv = 2xnr_elems array of node3-node1
    #ww = 2xnr_elems array node1-node2
    #ar = determinant of hcat(ones(3), nodes_x, nodes_y)

    #corner point indices
    a1 = vec(mesh.elements[1,:])
    a2 = vec(mesh.elements[2,:])
    a3 = vec(mesh.elements[3,:])

    # generate stiffness matrix
    q1 = mesh.nodes[:,a1]
    q2 = mesh.nodes[:,a2]
    q3 = mesh.nodes[:,a3]
    uu = q2-q3
    vv = q3-q1
    ww = q1-q2
    ar = (uu[1,:].*vv[2,:]-uu[2,:].*vv[1,:])./2


    #stiffness matrix assembly
    ar4 = abs(ar.*4)
    Kg = zeros(9,nme)
    Kg[1,:] = sum(uu.*ph_prop_11.*uu,1)./ar4
    Kg[2,:] = sum(vv.*ph_prop_11.*uu,1)./ar4
    Kg[3,:] = sum(ww.*ph_prop_11.*uu,1)./ar4
    Kg[5,:] = sum(vv.*ph_prop_11.*vv,1)./ar4
    Kg[6,:] = sum(ww.*ph_prop_11.*vv,1)./ar4
    Kg[9,:] = sum(ww.*ph_prop_11.*ww,1)./ar4
    Kg[[4,7,8],:] = Kg[[2,3,6],:]
    Ig = vec(mesh.elements[[1, 2, 3, 1, 2, 3, 1, 2, 3],:])
    Jg = vec(mesh.elements[[1, 1, 1, 2, 2, 2, 3, 3, 3],:])
    Kg = vec(Kg[:])

    #load vector assembly
    #this uses types uu, vv, ww, ar from above, for basis
    #calculate center of triangles, rhs is vector of function
    #values at centers
 	centroid_x = (mesh.nodes[1,a1] + mesh.nodes[1,a2] + mesh.nodes[1,a3])./3
   	centroid_y = (mesh.nodes[2,a1] + mesh.nodes[2,a2] + mesh.nodes[2,a3])./3
    rhs = map(f,centroid_x,centroid_y)

    Irhs = vec([a1 a2 a3])
    Jrhs = vec(ones(Int64,length(Irhs)))
    rhs = [rhs rhs rhs]
    rhsint = rhs.*abs([ar ar ar])./3 
    Krhs = vec(rhsint)


	# implement boundary conditions
	# find out which BC applies where
	neu_bd = div(findin(bddata,'N'),3) + 1
	diri_bd = div(findin(bddata,'D'),3) + 1
	#check if all edges are covered
	if length(neu_bd) + length(diri_bd) != size(bddata,2)
		error("Could not understand boundary file. Not all parts of the boundary have BC assigned.")
	end
		
	# Begin Neumann Part
	if neu_bd != []
   		neu_nodes1 = []
   		neu_nodes2 = []
   		neu_fct = []
   		for j = 1:size(neu_bd,1)
   			neu_edges = findin(mesh.edges[4,:],bddata[1,neu_bd[j]])
   			append!(neu_nodes1,vec(mesh.edges[1,neu_edges]))
   			append!(neu_nodes2,vec(mesh.edges[2,neu_edges]))
   			append!(neu_fct,vec(map(bddata[3,neu_bd[j]],(mesh.nodes[1,vec(mesh.edges[1,neu_edges])]+mesh.nodes[1,vec(mesh.edges[2,neu_edges])])/2,(mesh.nodes[2,vec(mesh.edges[1,neu_edges])]+mesh.nodes[2,vec(mesh.edges[2,neu_edges])])/2)))
   		end

   		# The flow is given by the surface integral over the respective edges.
   		# We use the formula int_dE gn*eta ds ~ |E|/2*g(center_x,center_y)
   		
   		##Calculate the flows at the edges.
   		neumanndata = sqrt((mesh.nodes[1,neu_nodes1] - mesh.nodes[1,neu_nodes2]).^2 + (mesh.nodes[2,neu_nodes1] - mesh.nodes[2,neu_nodes2]).^2)'.*neu_fct/2

    	#The vector [neu_nodes1;neu_nodes2] together with the repeated neumanndata adds all the Neumann data to the correct nodes.
  		append!(Irhs,vec([neu_nodes1;neu_nodes2]))
   		Jrhs = vec(ones(Int64,length(Irhs)))
   		append!(Krhs,repmat(vec(neumanndata),2))

	end

	# Begin Dirichlet Part
	diri_nodes = []
	diri_values = []
	for j = 1:size(diri_bd,1)
		diri_edges = findin(mesh.edges[4,:],bddata[1,diri_bd[j]])
		append!(diri_nodes,unique(vec(mesh.edges[1:2,diri_edges])))
       	append!(diri_values,map(node->bddata[3,diri_bd[j]](mesh.nodes[1,node],mesh.nodes[2,node]),diri_nodes[end-length(diri_edges):end]))
    end

    diri_data = unique([diri_nodes diri_values],1)
    gD = map(Float64,diri_data[:,2])
    n = diri_data[:,1]

	# Prepare for changing stiffness matrix entries
    rows = findin(Ig, n)
    cols = findin(Jg, n)

    # bdM gives the matrix necessary for subtraction from RHS
    bdM = sparse(vcat(Ig[rows],Ig[cols]),vcat(Jg[rows],Jg[cols]),vcat(Kg[rows],Kg[cols]),nq,nq)  
    # Change stiffness matrix entries

    Kg[rows] = 0
    Kg[cols] = 0
    append!(Ig, n)
    append!(Jg, n)
    append!(Kg, ones(Int64,length(n)))

	# create matrix
    M = sparse(Ig[:], Jg[:], Kg[:], nq, nq)
       
    # Subtract the columns from the right side.
    b = full(sparse(Irhs,Jrhs,Krhs,nq,1))  

    b -= bdM[:, n] * gD
    b[n] = gD

    return M, b
end

###############################################################################
function solve(meshpath::AbstractString,A::Float64,f::Function,bddata::Array)
# reads the mesh, assembles stiffness matrix and right hand side, and solves the 
# matrix equation
	mesh = read_mesh(meshpath)
	M,b=assemble(mesh,A,f,bddata)
	u=M\b
	
	return mesh,u
end

###############################################################################
function plot_solution(mesh::Mesh, u::Array; fig::Int = 2)
#plot_solution plots a solution u solved on the mesh "mesh"
    x = vec(mesh.nodes[1,:])
    y = vec(mesh.nodes[2,:])
    l = length(mesh.elements)
    triangles = mesh.elements[1:3,:]'-1
    figure(fig)
    hold(true)
    clf()
    plot_trisurf(x,y,triangles,vec(u),cmap=ColorMap("jet"))
    hold(false)
end


###############################################################################
function calculate_norm(mesh,u)
# Calculate the 2-norm of the solution.

	a1 = vec(mesh.elements[1,:])
    a2 = vec(mesh.elements[2,:])
    a3 = vec(mesh.elements[3,:])

    q1 = mesh.nodes[:,a1]
    q2 = mesh.nodes[:,a2]
    q3 = mesh.nodes[:,a3]
    uu = q2-q3
    vv = q3-q1
    ww = q1-q2
    ar = (uu[1,:].*vv[2,:] - uu[2,:].*vv[1,:])./2

    sol = ((u[a1] + u[a2] + u[a3])./3)
    sol2 = ones(size(sol,1),1)
    norm = sqrt(sum(abs(ar'.*sol.^2)))
    return norm
end

end # module
