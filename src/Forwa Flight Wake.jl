# key words in forward flight
V #forward velocity
  #thus the rotor blade loading and motion are periodic with a fundamental
  #frequency equal to the rotor speed Ω
A = π*R^2 #the area of a circle with 1/2 diameter
          #also: Cdi = CL^2/(π*A*R)
vi #the induced velocity of the helicopter rotor in high-speed forward flight
vi = T/(2*ρ*A*V)
AR = R/c = (N/π)*σ #aspect ratio should be high
U = sqrt(uT^2+uP^2) #
Φ = atan(uP/UT)
α = θ - Φ #the aerodynamic angle-of-attack of the blade element
L = 1/2*ρ*U^2*c*cl
D = 1/2*ρ*U^2*c*cd
Fz = L*cos(Φ)-D*sin(Φ)
Fx = L*sin(Φ)+D*cos(Φ)
dT = N*Fz*dr
dQ = N*Fx*r*dr
dP = Ω*dQ = N*Fx*Ω*r*dr
