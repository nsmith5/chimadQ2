module Dynamics

include("Constants.jl")

# Export time stepper
export step, volume, free_energy

# Name scalar field and complex scalar fields

typealias field	Array{Float64, 3}
typealias ffield Array{Complex128, 3}

# Auxiliary Arrays for time stepping

const Φ   = Array(Complex128, Nx>>1 + 1, Ny, Nz)	# 𝐅[ϕ]
const nl  = Array(Float64, Nx, Ny, Nz)				# non-linear term
const fnl = Array(Complex128, Nx>>1 + 1, Ny, Nz)	# 𝐅[non-linear term]

# Construct fft plan

FFTW.set_num_threads(8)		
const plan = plan_rfft(nl, [1,2,3], flags=FFTW.MEASURE) # Real->half complex on axis 1,2,3

# Methods

function nonlinear_term(ϕ::field, nl::field)
	# Compute non-linear term for equation of motion:
	#  	nl(x) = -6 h₀ϕ(x)² + 4h₀ϕ(x)³ + ϵ*σ*h(ϕ)
	for i in 1:Nx*Ny*Nz
		@inbounds begin
			nl[i] = h₀ * ϕ[i]^2 * (4.0 * ϕ[i] - 6.0)
		end
	end
	return
end


function step(ϕ::field)
	
	# Propagate ϕ(t) -> ϕ(t + Δt)
	
	A_mul_B!(Φ, plan, ϕ)
	nonlinear_term(ϕ, nl)
	A_mul_B!(fnl, plan, nl)
	
	for i in 1:Nz*Ny*(Nx>>1 + 1)
		@inbounds begin
			pref = 1.0/(1.0 - dt*M_ϕ*∇²[i] * (2*h₀ - κ * ∇²[i]))
			Φ[i] = pref * (Φ[i] + dt * M_ϕ * ∇²[i] * fnl[i])
		end
	end
	A_ldiv_B!(ϕ, plan, Φ)
	
	return
end	

function free_energy(ϕ::field)
	f = ϕ -> ϕ^2*(ϕ-1)^2
	bulk = sum([h₀*f(ϕᵢ)*dx*dy*dz for ϕᵢ in ϕ])
	∇²Φ = similar(Φ)
	∇²ϕ = similar(ϕ)
	A_mul_B!(Φ, plan, ϕ)
	for i in 1:(Nx>>1 + 1)*Ny*Nz
		@inbounds ∇²Φ[i] = ∇²[i]*Φ[i]
	end
	A_ldiv_B!(∇²ϕ, plan, ∇²Φ)
	surface = sum([-0.5*κ*ϕ[i]*∇²ϕ[i]*dx*dy*dz for i in length(ϕ)])
	return bulk + surface
end

function volume(ϕ)
	return sum(ϕ)*dx*dy*dz
end

end 



