using Turing
include("julia/PartitionDistribution.jl")

@model function bfvar_gamma_eq(n::VecFloat, b::VecFloat, α::VecFloat, ::Type{T} = Float64) where {T}

	k = length(n)
	d ~ BetaBinomial(k - 1, 1, 1)
	equal_indices ~ PartitionDistribution(d, k)

	τ ~ Gamma(0.001, 0.001)
	gammas = Vector{T}(undef, k)
	for i in 1:k
		gammas[i] ~ Gamma(α[i], 1)
	end
	ρ = Vector{T}(undef, k)
	for i in 1:k
		ρ[i] = mean(gammas[equal_indices .== equal_indices[i]])
	end
	ρ ./= sum(ρ)
	Turing.@addlogprob! myloglikelihood(n, b, ρ, τ)

end
myloglikelihood
get_eq_ind_nms(samples) = filter(x->startswith(string(x), "equal_indices"), samples.name_map.parameters)

function compute_post_prob_eq(samples)
	eq_ind_nms = get_eq_ind_nms(samples)
	s = size(samples[eq_ind_nms])
	samps = reshape(samples[eq_ind_nms].value.data, :, s[2])
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

sds = Float64[1, 1, 1]
k   = length(sds)
ns  = 1000 .* ones(Int, k)

ss = (sds .* ((ns .- 1) ./ ns)).^2
n = (ns .- 1) ./ 2
b = ns .* ss
α = ones(Float64, length(ss))

spl = Gibbs(HMC(0.2, 3, :τ, :gammas), PG(20, :d, :equal_indices))
model_gamma_eq = bfvar_gamma_eq(n, b, α)
samples_gamma_eq = sample(model_gamma_eq, spl, 10_000)

StatsBase.countmap(vec(samples_gamma_eq[:d].data))
mean(samples_gamma_eq[:d])
compute_post_prob_eq(samples_gamma_eq)

dimension = vec(Int.(samples_gamma_eq[:d]))
eq_ind_nms = filter(x->startswith(string(x), "equal_indices"), samples.name_map.parameters)
equal_indices = Int.(samples[eq_ind_nms].value.data)
print_probs(dimension, equal_indices)


logpdf(BetaBinomial(2, 1, 1), 0.01)
function logpdf_bb(k, n, α, β)
	SpecialFunctions.logabsbinomial(n, k)[1] + SpecialFunctions.lbeta(k + α, n - k + β) - SpecialFunctions.lbeta(α, β)
end

logpdf.(BetaBinomial(2, 1, 1), 0:2)
logpdf_bb.(0:2, 2, 1, 1)

SpecialFunctions.logabsbinomial(2, 0)[1]
SpecialFunctions.lbeta(0 + 1, 2 - 0 + 1)
SpecialFunctions.lbeta(1, 1)