using Dynamics
using JLD

D = Dynamics

const phi = zeros(D.Nx, D.Ny, D.Nz)

r = 0.1e-6 			# radius of ball is	0.1μm
L = (4.2e-21)^(1/3)

for z in 1:D.Nz
	for y in 1:D.Ny
		for x in 1:D.Nx
			if (D.dx*(x-1-D.Nx>>1))^2 < (L/2)^2 &&
			   (D.dy*(y-1-D.Ny>>1))^2 < (L/2)^2 &&
			   (D.dz*(z-1-D.Nz>>1))^2 < (L/2)^2
				phi[x, y, z] = 1.0
			end
		end
	end
end					


