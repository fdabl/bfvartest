using Turing # 0.14.1 last time I checked

#=

	TODO:

		write model using gamma variables

=#

# consider SVector!
const VecFloat = Vector{Float64}
# const VecInt = Vector{Int}

function myloglikelihood(n, b, ρ, τ)

	prec = ρ * τ * length(n)
	out = 
		-logpdf(Gamma(0.001, 0.001), τ) +
		-log(τ) +
		n' * log.(prec) +
		-0.5 * (prec' * b)
	return out
end

@model function bfvar_rho(n::VecFloat, b::VecFloat, α::VecFloat)

	τ ~ Gamma(0.001, 0.001)
	ρ ~ Dirichlet(α)
	Turing.@addlogprob! myloglikelihood(n, b, ρ, τ)

end

@model function bfvar_gamma(n::VecFloat, b::VecFloat, α::VecFloat, ::Type{T} = Float64) where {T}

	τ ~ Gamma(0.001, 0.001)
	gammas = Vector{T}(undef, length(n))
	for i in eachindex(gammas)
		gammas[i] ~ Gamma(α[i], 1)
	end
	ρ = gammas ./ sum(gammas)
	Turing.@addlogprob! myloglikelihood(n, b, ρ, τ)

end

function samples_2_array(samples, has_gamma::Bool = false)
	# transform turing array to useful samples
	samps = samples[samples.name_map.parameters].value.data[:, :, 1]
	k = size(samps)[2] - 1
	if has_gamma
		gammas = @view samps[:, 1:k]
		rhos = similar(gammas)
		for i in axes(rhos, 1)
			rhos[i, :] = gammas[i, :] ./ sum(gammas[i, :])
		end
	else
		rhos = @view samps[:, 1:k]
	end
	taus = @view samps[:, k + 1]
	@assert all(x->(abs(x - 1)) < √eps(), sum(rhos, dims = 2))
	precs = similar(rhos)
	for i in axes(rhos, 1)
		precs[i, :] = (taus[i] * k) .* rhos[i, :]
	end 
	stds = 1 ./ sqrt.(precs)
	out = hcat(samps, precs, stds)
	means = Dict{Symbol, Vector{Float64}}(
		:ρ => vec(mean(rhos, dims = 1)),
		:τ => [mean(taus)],
		:σ => vec(mean(stds, dims = 1)),
	)
	return means, out
end

sds = Float64[1, 2, 3]
k   = length(sds)
ns  = 100 .* ones(Int, k)

ss = (sds .* ((ns .- 1) ./ ns)).^2
n = (ns .- 1) ./ 2
b = ns .* ss
α = ones(Float64, length(ss))

model_rho = bfvar_rho(n, b, α)
samples_rho = sample(model_rho, HMC(0.05, 10), 10_000)

means_rho, samples2_rho = samples_2_array(samples_rho)
means_rho

model_gamma = bfvar_gamma(n, b, α)
samples_gamma = sample(model_gamma, HMC(0.05, 10), 10_000)

means_gamma, samples2_gamma = samples_2_array(samples_gamma, true)
means_gamma