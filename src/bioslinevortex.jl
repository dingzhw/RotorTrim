# This file solve the line vortex Bio-Savar induced velocity

function bslinev(gama,pointa,pointb,pointiv,t)
  xm = pointiv[1]
  ym = pointiv[2]
  zm = pointiv[3]
  x1 = pointa[1]
  y1 = pointa[2]
  z1 = pointa[3]
  x2 = pointb[1]
  y2 = pointb[2]
  z2 = pointb[3]
  xl = x2-x1
  yl = y2-y1
  zl = z2-z1
  xr1 = xm-x1
  yr1 = ym-y1
  zr1 = zm-z1
  xr2 = xm-x2
  yr2 = ym-y2
  zr2 = zm-z2
  lvec = [xl,yl,zl]
  r1 = [xr1,yr1,zr1]
  r2 = [xr2,yr2,zr2]
  cosα = cosvecang(r1,lvec)
  cosβ = cosvecang(r2,lvec)
  if cosα<=1e-3||cosβ<=1e-3
    return [0.0,0.0,0.0]
  else
    sinα = sqrt(abs(1-cosα^2))
    h = abs(norm(r1)*sinα)
    vvec = [yl*zr1-zl*yr1,
            xl*zr1-zl*xr1,
            xl*yr1-yl*xr1] # the vector direction of induced velocity
    if norm(vvec)<=1e-3
      return [0.0,0.0,0.0]
    else
      vvec_ = vvec/norm(vvec)
    end
    v_θ = vcore(gama,h,t,3)[1]
    vind = v_θ*(cosα-cosβ)*vvec_
    return vind
  end
end

function vcore(Γ_avg,h,t,n=3)  # the vortex core model
  # v_θ(r) = Γ*r/(2*π*(rc^(2*n)+R^(2*n))^(1/n))
  # n is a Int VARIABLE
  # r is the distance between center of vortex to the calculated point
  # rc is the radius of vortex core
  # rc = 1.12*√(4*δ*ν*t)
  # ν is the coefficient of kinematic viscosity
  # t is the wake vortex age time
  # t = ζ/Ω (ζ is the wake vortex age angle)
  # δ = 1+a1*Γ_avg/ν
  # Γ_avg is the mean value of Γ
  # Γ_avg = 1/Nmax*sum(Γ)
  # a1 is an empirical parameter
  # a1 = 0.1

  a1 = 0.1
  δ = 1+a1*abs(Γ_avg)/ν # here I could refer to Tan J's codes for correction
  rc = 0.00855*√(δ*t)
  v_θ = Γ_avg*h/(2*π*(rc^(2*n)+h^(2*n))^(1/n))
  return v_θ,rc
end

function bslinevnocore(gama,pointa,pointb,pointiv)  #Bio-Savar calculate for Bound Vortex
    xm = pointiv[1]
    ym = pointiv[2]
    zm = pointiv[3]
    x1 = pointa[1]
    y1 = pointa[2]
    z1 = pointa[3]
    x2 = pointb[1]
    y2 = pointb[2]
    z2 = pointb[3]
    xl = x2-x1
    yl = y2-y1
    zl = z2-z1
    xr1 = xm-x1
    yr1 = ym-y1
    zr1 = zm-z1
    xr2 = xm-x2
    yr2 = ym-y2
    zr2 = zm-z2
    lvec = [xl,yl,zl]
    r1 = [xr1,yr1,zr1]
    r2 = [xr2,yr2,zr2]
    cosα = cosvecang(r1,lvec)
    cosβ = cosvecang(r2,lvec)
    if abs(cosα)<=1e-3||abs(cosβ)<=1e-3
      return 0.0
    else
      sinα = sqrt(abs(1-cosα^2))
      h = abs(norm(r1)*sinα)
      vvec = [yl*zr1-zl*yr1,
              xl*zr1-zl*xr1,
              xl*yr1-yl*xr1] # the vector direction of induced velocity
      if norm(vvec)<=1e-3
        return 0.0
      else
        vvec_ = vvec/norm(vvec)
      end
      v_θ = gama/(4*π*h)
      vind = v_θ*(cosα-cosβ)*vvec_
      return vind
    end
  end
