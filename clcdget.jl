# This file is to get the Cl and Cd of Blade Element


#the excel interpolation method
function cla(α,ma)
  # NACA0012
  x = α/π*180
  if ma<0.1
    cla = 2e-17*x^6 + 1e-10*x^5 - 8e-13*x^4 - 4e-06*x*3 + 1e-08*x^2 + 0.0351*x + 1e-06
  elseif ma<0.2
    cla = 9e-18*x^6 + 1e-10*x^5 - 5e-13*x^4 - 4e-06*x^3 + 6e-09*x^2 + 0.0351*x - 3e-06
  elseif ma<0.3
    cla = 2e-17*x^6 + 1e-10*x^5 - 8e-13*x^4 - 4e-06*x^3 + 1e-08*x^2 + 0.0351*x + 1e-06
  elseif ma<0.4
    cla = -1e-15*x^6 + 1e-10*x^5 + 6e-11*x^4 - 4e-06*x^3 - 8e-07*x^2 + 0.0349*x + 0.0025
  elseif ma<0.5
    cla = -3e-17*x^6 + 1e-10*x^5 + 1e-12*x^4 - 4e-06*x^3 - 3e-09*x^2 + 0.0349*x - 0.0001
  elseif ma<0.6
    cla = -6e-17*x^6 + 1e-10*x^5 + 2e-12*x^4 - 4e-06*x^3 - 9e-09*x^2 + 0.0347*x - 0.0002
  elseif ma<0.7
    cla = -6e-17*x^6 + 1e-10*x^5 + 2e-12*x^4 - 4e-06*x^3 - 9e-09*x^2 + 0.0345*x - 0.0002
  elseif ma<0.8
    cla = -6e-17*x^6 + 1e-10*x^5 + 2e-12*x^4 - 4e-06*x^3 - 9e-09*x^2 + 0.0345*x - 0.0002
  elseif ma<0.9
    cla = 2e-16*x^6 + 9e-11*x^5 - 7e-12*x^4 - 4e-06*x^3 + 4e-08*x^2 + 0.0313*x + 0.0004
  elseif ma<1.0
    cla = 1e-16*x^6 + 9e-11*x^5 - 5e-12*x^4 - 4e-06*x^3 - 4e-09*x^2 + 0.03*x + 0.0007
  else
    cla = 1e-16*x^6 + 9e-11*x^5 - 5e-12*x^4 - 4e-06*x^3 - 7e-09*x^2 + 0.03*x + 0.0007
  end
  return cla
end

#the simple interpolation method
function liftcoffi(α,Ma)
  if abs(α)>(20.0*π/180.0)
    Cl = 0.0
  else
    Cl = 5.73*α+0.085
  end
  return Cl
end

function dragcoffi(α,Ma)
  if abs(α)<(3.0*π/180.0)
    Cd = 0.0076
  elseif abs(α)<(5.0*π/180.0)
    Cd = 0.008
  elseif abs(α)<(6.0*π/180.0)
    Cd = 0.0086
  elseif abs(α)<(7.0*π/180.0)
    Cd = 0.0093
  elseif abs(α)<(8.0*π/180.0)
    Cd = 0.0103
  elseif abs(α)<(9.0*π/180.0)
    Cd = 0.012
  elseif abs(α)<(10.0*π/180.0)
    Cd = 0.0146
  elseif abs(α)<(11.0*π/180.0)
    Cd = 0.0177
  elseif abs(α)<(12.0*π/180.0)
    Cd = 0.0215
  elseif abs(α)<(13.0*π/180.0)
    Cd = 0.0266
  elseif abs(α)<(14.0*π/180.0)
    Cd = 0.055
  elseif abs(α)<(15.0*π/180.0)
    Cd = 0.125
  elseif abs(α)<(16.0*π/180.0)
    Cd = 0.21
  else
    Cd = 0.5
  end

  Cd = Cd*1.05  #The real blade correction
  return Cd
end

function clcdget(α_aero,vall_r)
  Ma = Array{Float64}(NR,Nbe)
  for k in 1:NR #take Mach Number into consideration
    for i in 1:Nbe
      Ma[k,i] = norm(vall_r[k,i])/v_sound
    end
  end
  clift = Array{Float64}(NR,Nbe)
  cdrag = Array{Float64}(NR,Nbe)
  for k in 1:NR #solve the local clcd
    for i in 1:Nbe
      clift[k,i] = cla(α_aero[k,i],Ma[k,i])
      cdrag[k,i] = dragcoffi(α_aero[k,i],Ma[k,i])
    end
  end
  return clift,cdrag,Ma
end
