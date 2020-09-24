using Turing
include("julia/PartitionDistribution.jl")

function compute_probs(samps)
	n_samps, n_groups = size(samps)
	probs = zeros(Float64, n_groups, n_groups)
	for row in eachrow(samps)
		for j in eachindex(row)
			idx = j .+ findall(==(row[j]), row[j+1:end])
			probs[idx, j] .+= 1.0
		end
	end
	return probs ./ n_samps
end

function print_probs(dimension, equal_indices)
	
	for d in sort(unique(dimension))
		idx = findall(==(d), dimension)
		println("times dimension == $(d): $(length(idx)) times")
		probs = compute_probs(view(equal_indices, idx, :))
		println("probs:")
		display(probs)
	end
end

@model function prior(k::Int, α::Float64 = 1.0, β::Float64 = 1.0)

	θ ~ Beta(α, β)
	d ~ Binomial(k - 1, θ)
	equal_indices ~ PartitionDistribution(d, k)

end

spl = Gibbs(HMC(0.2, 3, :θ), PG(20, :d), PG(20, :equal_indices))
samples = sample(prior(3), spl, 10_000)

dimension = vec(Int.(samples[:d]))
equal_indices = Int.(samples[samples.name_map.parameters[2:4]].value.data)
print_probs(dimension, equal_indices)

#= 
	TODO: 
		figure out why arent all edges connected when d = 0 and no edges connected when d = 2

=#