### Constants for part (a)

const h₀ = 2.4e8	# [J/m³] Barrier term
const κ  = 3e-9		# [J/m]	Surface energy 

const Lx = 0.3e-6	# [m] Physical length simulation
const Ly = 0.3e-6	
const Lz = 0.3e-6	

const Nx = 128		# Number of grid points
const Ny = 128
const Nz = 128

const dx = Lx/Nx	# [m] Spacing between grid points
const dy = Ly/Ny
const dz = Lz/Nz

const kx0 = 2π/Lx	# [1/m] Spatial freq. of fundamental 
const ky0 = 2π/Ly	
const kz0 = 2π/Lz

const dt = dx		# [s] 'Time' step size

const M_ϕ = 0.001 	# [1/s] 'Mobility'

function k²(i,j,k)
	# Computes |k⃗|² 
	kx² = i < Nx>>1 ? (kx0*(i-1))^2 : (kx0*(i-1-Nx))^2
	ky² = j < Ny>>1 ? (ky0*(j-1))^2 : (ky0*(j-1-Ny))^2
	kz² = k < Nz>>1 ? (kz0*(k-1))^2 : (kz0*(k-1-Nz))^2
	return kx² + ky² + kz² 
end

const ∇² = [-k²(i,j,k) for i in 1:Nx>>1 + 1, j in 1:Ny, k in 1:Nz] 

