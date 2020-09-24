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

	d ~ BetaBinomial(k - 1, α, β)
	equal_indices ~ PartitionDistribution(d, k)
	if !isfinite(logpdf(PartitionDistribution(d, k), equal_indices))
		# Exit the model evaluation early
		Turing.@addlogprob! -Inf
		return
	end

end

n_groups = 4
samples = sample(prior(n_groups), PG(20, :d, :equal_indices), 10_000)

dimension = vec(Int.(samples[:d]))
equal_indices = Int.(samples[samples.name_map.parameters[2:1+n_groups]].value.data)
print_probs(dimension, equal_indices)
