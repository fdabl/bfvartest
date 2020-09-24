using Turing

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

@model function bfvar(n::VecFloat, b::VecFloat, α::VecFloat)

	τ ~ Gamma(0.001, 0.001)
	ρ ~ Dirichlet(α)
	Turing.@addlogprob! myloglikelihood(n, b, ρ, τ)

end

function samples_2_array(samples)
	# transform turing array to useful samples
	samps = samples[samples.name_map.parameters].value.data[:, :, 1]
	k = size(samps)[2] - 1
	rhos = @view samps[:, 1:k]
	taus = @view samps[:, k + 1]
	@assert all(isone, sum(rhos, dims = 2))
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

model = bfvar(n, b, α)
samples = sample(model, HMC(0.05, 10), 10_000)

means, samples2 = samples_2_array(samples)
means