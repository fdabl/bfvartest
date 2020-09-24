using Distributions, StatsBase, Random, DataStructures


struct PartitionDistribution{T} <: Distribution{Multivariate, Discrete} where T <: Integer
	partitions::T
	vertices::T
	function PartitionDistribution(partitions::T, vertices::T) where T <: Integer
		vertices > 1 || throw(DomainError(vertices, "condition vertices > 1 is violated"))
		-1 < partitions < vertices || throw(DomainError(partitions, "condition 0 < partitions < vertices - 1 is violated")) 
		new{T}(partitions, vertices)
	end
end

Base.length(x::PartitionDistribution{T}) where T<:Integer = x.vertices
Base.eltype(x::PartitionDistribution{T}) where T<:Integer = T

# function _rand_indices!(::AbstractRNG, d::PartitionDistribution{T}, u::AbstractMatrix) where T <: Integer

#     # everything disconnected
#     iszero(d.partitions) && return u

#     # everything connected
#     if d.partitions + one(T) == d.vertices
#         # everything is connected so just return everything connected to 1
#         o = one(T)
#         for j in axes(u, 2)
#             count = one(T)
#             for i in axes(u, 1)
#                 u[i, j] = [o, o + count]
#                 count += o
#             end
#         end
#         return u
#     end

#     vertices = T.(collect(1:d.vertices))
#     for j in axes(u, 2)
#         omitted = Set{T}()
#         current_view = @view vertices[(1:end) .∉ Ref(omitted)]
#         for i in axes(u, 1)

#             # sample 2 edges to be equal -- both 1, 2 and 2, 1 are possible
#             idx = sort(sample(current_view, 2; replace = false))
#             u[i, j] = idx
#             # remove the vertex with larger index from the view (i.e., "shrink" a dimension)
#             push!(omitted, idx[2])
#             current_view = @view vertices[(1:end) .∉ Ref(omitted)]

#         end
#     end
#     return u
# end

function Distributions._rand!(::AbstractRNG, d::PartitionDistribution{T}, u::AbstractVector) where T <: Integer
	if iszero(d.partitions)
		u[:] .= 1:d.vertices
		return u
	end

	# everything connected
	if d.partitions + one(T) == d.vertices 
		fill!(u, one(T))
		return u
	end

	vertices = T.(collect(1:d.vertices))
	omitted = Set{T}()
	# buckets = Vector{Set{T}}() # <- should this be DataStructures.IntDisjointSets?
	buckets = DataStructures.IntDisjointSets(d.vertices)
	current_view = @view vertices[(1:end) .∉ Ref(omitted)]

	# new approach with DataStructures -- guarantees that first element is always a 1
	for i in 1:d.partitions

		idx = sort(sample(current_view, 2; replace = false))
		union!(buckets, idx[1], idx[2])
		push!(omitted, idx[2])
		current_view = @view vertices[(1:end) .∉ Ref(omitted)]

	end
	# empty!(omitted)
	parents = buckets.parents
	uniq = sort(unique(parents))
	for i in eachindex(uniq)
		u[parents .== uniq[i]] .= i
	end

	# old approach
	# for i in 1:d.partitions
	# 	buckets = DataStructures.IntDisjointSets(d.vertices)
	# 	current_view = @view vertices[(1:end) .∉ Ref(omitted)]
# 
		# # sample 2 edges to be equal -- both 1, 2 and 2, 1 are possible (and equally likely)
		# idx = sort(sample(current_view, 2; replace = false))
		# if isempty(buckets)
		# 	push!(buckets, Set(idx))
		# else
		# 	found = false
		# 	for b in buckets
		# 		if idx[1] ∈ b
		# 			found = true
		# 			push!(b, idx[2])
		# 		elseif idx[2] ∈ b
		# 			found = true
		# 			push!(b, idx[1])
		# 		end
		# 		# there should be a break if it's found!
		# 	end
		# 	if !found
		# 		push!(buckets, Set(idx))
		# 	end
		# end
		# # remove the vertex with larger index from the view (i.e., "shrink" a dimension)
		# push!(omitted, idx[2])
		# current_view = @view vertices[(1:end) .∉ Ref(omitted)]

	# end
	# u[:] .= zero(T)
	# for (i, b) in enumerate(buckets)
	# 	u[collect(b)] .= i
	# end
	# idx = findall(iszero, u[:])
	# n_buckets = length(buckets)
	# u[idx] .= n_buckets + 1:n_buckets + length(idx)
	return u
end

function Distributions._rand!(::AbstractRNG, d::PartitionDistribution{T}, u::AbstractMatrix) where T <: Integer

	# println(size(u))
	println(typeof(u))

	# everything disconnected
	if iszero(d.partitions)
		# println("fully disconnected")
		for j in axes(u, 2)
			u[:, j] .= 1:d.vertices
		end
		return u
	end

	# everything connected
	if d.partitions + one(T) == d.vertices 
		# println("fully connected")
		# everything is connected so just return everything connected to 1
		fill!(u, one(T))
		return u
	end

	vertices = T.(collect(1:d.vertices))
	omitted = Set{T}()
	# buckets = Vector{Set{T}}()
	for j in axes(u, 2)
		buckets = DataStructures.IntDisjointSets(d.vertices)
		current_view = @view vertices[(1:end) .∉ Ref(omitted)]

		# new approach with DataStructures -- guarantees that first element is always a 1
		for i in 1:d.partitions

			idx = sort(sample(current_view, 2; replace = false))
			union!(buckets, idx[1], idx[2])
			push!(omitted, idx[2])
			current_view = @view vertices[(1:end) .∉ Ref(omitted)]

		end
		empty!(omitted)
		parents = buckets.parents
		uniq = sort(unique(parents))
		for i in eachindex(uniq)
			u[parents .== uniq[i], j] .= i
		end

		# old approach
		# for i in 1:d.partitions

		# 	# sample 2 edges to be equal -- both 1, 2 and 2, 1 are possible (and equally likely)
		# 	idx = sort(sample(current_view, 2; replace = false))
		# 	if isempty(buckets)
		# 		push!(buckets, Set(idx))
		# 	else
		# 		found = false
		# 		for b in buckets
		# 			if idx[1] ∈ b
		# 				found = true
		# 				push!(b, idx[2])
		# 			elseif idx[2] ∈ b
		# 				found = true
		# 				push!(b, idx[1])
		# 			end
		# 		end
		# 		if !found
		# 			push!(buckets, Set(idx))
		# 		end
		# 	end
		# 	# remove the vertex with larger index from the view (i.e., "shrink" a dimension)
		# 	push!(omitted, idx[2])
		# 	current_view = @view vertices[(1:end) .∉ Ref(omitted)]

		# end
		# u[:, j] .= zero(T)
		# for (i, b) in enumerate(buckets)
		# 	u[collect(b), j] .= i
		# end
		# idx = findall(iszero, u[:, j])
		# n_buckets = length(buckets)
		# u[idx, j] .= n_buckets + 1:n_buckets + length(idx)
		# empty!(omitted)
		# empty!(buckets)
	end
	return u

end

function Distributions.logpdf(d::PartitionDistribution{T}, x::Vector{U}) where {T <: Integer, U <: Integer}
	# doesn't matter for now? we could count this though
	invalidx = !iszero(length(unique(u)) - d.vertices + d.partitions)
	invalidx && return -Inf
	return 0.0
	# return 0 
end

# tests
# d = PartitionDistribution(3, 6)
# rand(d, 5)
# rand(PartitionDistribution(2, 4), 1)
# rand(PartitionDistribution(3, 4), 1)
# length(d)
#
function count_probs(d::PartitionDistribution, nmax::Int = 10_000)

	probs = zeros(Float64, d.vertices, d.vertices)
	for i in 1:nmax
		u = rand(d, 1)
		for j in eachindex(u)
			idx = j .+ findall(==(u[j]), u[j+1:end])
			probs[idx, j] .+= 1.0
		end
	end
	return probs ./ nmax

end

# d = PartitionDistribution(1, 4)
# u = Vector{Int}(undef, length(d))
# Distributions._rand!(Distributions.GLOBAL_RNG, d, u)
# rand(d, 1)

# d = PartitionDistribution(3, 5)
# u = rand(d, 1)
# length(unique(u)) - d.vertices + d.partitions

#= Tests
count_probs(PartitionDistribution(0, 4))
count_probs(PartitionDistribution(1, 4))
count_probs(PartitionDistribution(2, 4))
count_probs(PartitionDistribution(3, 4))


d = PartitionDistribution(0, 4)
rand(d, 10)

d = PartitionDistribution(1, 4)
rand(d, 10)

d = PartitionDistribution(2, 4)
rand(d, 10)

d = PartitionDistribution(3, 4)
rand(d, 10)

=#

# d = PartitionDistribution(1, 5)
# rand(d, 10)