# 这是一个风洞配平程序

function uisol(uitmp,θ7,θlat,θlon)  # 均匀入流求解力和力矩的流程
  # 用于雅可比矩阵中求平均周期力
  ψxx = 0.0
  iter = 1
  sumf = 0.0
  blat = 0.0
  blon = 0.0
  θ0 = Array{Float64}(NR,Nbe)
  θcp = zeros(Float64,NR)
  for k in 1:NR
    θcp[k] = θ7
    for i in 1:Nbe
      θ0[k,i] = θcp[k]+θtw*(rb[k,i]/R-0.7)
    end
  end

  while iter<=npsi
    vind_r = uitmp[1]
    vall_r = uitmp[2]
    beta = uitmp[4]
    dbeta = uitmp[5]
    blat = blat+uitmp[4]*sin(ψxx)
    blon = blon+uitmp[4]*cos(ψxx)

    # This file is used for get angle-of-attack
    aoatmp = aoaget(vall_r,ψxx,θ0,θlat,θlon,beta)
    α_aero = aoatmp[1]
    θ = aoatmp[2]

    # This file is to get the Cl and Cd of Blade Element
    clcdtmp = clcdget(α_aero,vall_r)
    Cl = clcdtmp[1]
    Cd = clcdtmp[2]

    # This file calculate the whole aerodynamic force of the rotor in the hub
    rftmp = rotorforce(ψxx,vall_r,α_aero,θ,Cl,Cd)
    sumf = sumf+rftmp[3]
    Mbeta_aero = rftmp[5][1]

    # The file calculate the uniform induced velocity
    uitmp = uniforminflow(ψxx,θcp,θlat,θlon,beta,dbeta,Mbeta_aero)

    ψxx = ψxx+dψ
    iter = iter+1
  end
  sumf = sumf/npsi
  blat = blat/npsi*2
  blon = blon/npsi*2

  return sumf,blat,blon
end

function yagb(uitmp,θcp,θlat,θlon,epsi=1/180*pi)
  # 前向差分雅可比矩阵

  # 后向差分的雅可比矩阵

  # 中心差分雅可比矩阵
  thecp = θcp+epsi/2
  uisoltmp = uisol(uitmp,thecp,θlat,θlon)
  rotfor = uisoltmp[1]
  blatfor = uisoltmp[2]
  blonfor = uisoltmp[3]

  thecp = θcp-epsi/2
  uisoltmp = uisol(uitmp,thecp,θlat,θlon)
  rotbac = uisoltmp[1]
  blatbac = uisoltmp[2]
  blonbac = uisoltmp[3]

  dthecp = (rotfor-rotbac)/epsi
  dtheblat = (blatfor-blatbac)/epsi
  dtheblon = (blonfor-blonbac)/epsi

  lat = θlat+epsi/2
  uisoltmp = uisol(uitmp,θcp,lat,θlon)
  rotfor = uisoltmp[1]
  blatfor = uisoltmp[2]
  blonfor = uisoltmp[3]

  lat = θlat-epsi/2
  uisoltmp = uisol(uitmp,θcp,lat,θlon)
  rotbac = uisoltmp[1]
  blatbac = uisoltmp[2]
  blonbac = uisoltmp[3]

  dlat = (rotfor-rotbac)/epsi
  dlatblat = (blatfor-blatbac)/epsi
  dlatblon = (blonfor-blonbac)/epsi

  lon = θlon+epsi/2
  uisoltmp = uisol(uitmp,θcp,θlat,lon)
  rotfor = uisoltmp[1]
  blatfor = uisoltmp[2]
  blonfor = uisoltmp[3]

  lon = θlon-epsi/2
  uisoltmp = uisol(uitmp,θcp,θlat,lon)
  rotbac = uisoltmp[1]
  blatbac = uisoltmp[2]
  blonbac = uisoltmp[3]

  dlon = (rotfor-rotbac)/epsi
  dlonblat = (blatfor-blatbac)/epsi
  dlonblon = (blonfor-blonbac)/epsi

  Myg = [dthecp  dlat  dlon;
          dtheblat  dlatblat  dlonblat;
          dtheblon  dlatblon  dlonblon]

  return Myg
end


function trimwt(uitmp,rot,beta_lat,beta_lon,θcp,θlat,θlon,epsi=1/180*π) # 纵横向挥舞配平法
  fzt = abs(rot*cos(αs)-T) # 力的配平条件
  epsi = 1/180*π # 雅可比求解的变化小量
  # 求解雅可比矩阵
  # 输入量：θ0,θ_lat,θ_lon
  # 输出量：rot,beta_lat,beta_lon

  if fzt<10&&abs(beta_lat)<1e-3&&abs(beta_lon)<1e-3 # 配平力和挥舞
  # if abs(beta_lat)<1e-3&abs(beta_lon)<1e-3 # 仅配平挥舞
    return true,θ0,θlat,θlon,θcp
  else
    Myg = yagb(uitmp,θcp[1],θlat,θlat,epsi) # 时变雅可比矩阵
    dc = inv(Myg)*[T-rot,-beta_lat,-beta_lon]
    for k = 1:NR
      θcp[k] = θcp[k]+dc[1]
      θlat[k] = θlat[k]+dc[2]
      θlon[k] = θlon[k]+dc[3]
      for i in 1:Nbe
        θ0[k,i] = θcp[k]+θtw*(rb[k,i]/R-0.7)
      end
    end
  end

  return false,θ0,θlat,θlon,θcp
end
