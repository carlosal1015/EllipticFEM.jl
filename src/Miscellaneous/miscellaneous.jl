
function myUniqueperm(a::Array{Int64,1})
	# Computes the indices of the unique entries in a.
	# Therefore unique(a) == a[myUniqueperm(a)].
	
	if !isempty(a)
		Pi = sortperm(a)
		aSorted = a[Pi]
		aSortedUniqueIndices = vcat([1],find(aSorted[2:1:end]-aSorted[1:1:(end-1)])+1)
		
		uniqueperm = sort(Pi[aSortedUniqueIndices])
	else
		uniqueperm = Array{Int64,1}(0)
	end
	
	return uniqueperm
end

function moduloIndex(index::Int64,n::Int64)
	return mod(index-1,n)+1
end

function projectValueIntoInterval(value::Float64,intervalCoords::Array{Float64})
	return max(intervalCoords[1],min(value,intervalCoords[2]))
end



