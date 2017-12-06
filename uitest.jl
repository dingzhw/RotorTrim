# The file test the uniform inflow condition

include("const.jl")
include("mathfunctions.jl")
include("flapreponse.jl")
include("uniforminflow.jl")
include("aoaget.jl")
include("clcdget.jl")
include("rotorforce.jl")
include("trim_windtu.jl")

tic() # 程序起始标志
print("===计算初始化===")

ψ = 0.0
beta_lat = 0.0
beta_lon = 0.0
θcp = zeros(Float64,NR)
θlat = zeros(Float64,NR)
θlon = zeros(Float64,NR)
θ0 = Array{Float64}(NR,Nbe)
for k in 1:NR
  θcp[k] = θ7
  θlat[k] = thelat
  θlon[k] = thelon
  for i in 1:Nbe
    θ0[k,i] = θcp[k]+θtw*(rb[k,i]/R-0.7)
  end
end
index = 1
betanow = βang0
betaxnow = dβ0
rot = 0.0
blat = 0.0
blon = 0.0
uitmp = uniforminflow(ψ,θcp,θlat,θlon,betanow,betaxnow)

while index<=60*npsi # 如果计算了100圈还没收敛就结束计算
  vind_r = uitmp[1]
  vall_r = uitmp[2]
  beta = uitmp[4]
  if index%npsi==0
    print("===当前挥舞角$(beta/π*180)°===\n")
  end
  dbeta = uitmp[5]
  blat = blat+uitmp[4]*sin(ψ)
  blon = blon+uitmp[4]*cos(ψ)

  # This file is used for get angle-of-attack
  aoatmp = aoaget(vall_r,ψ,θ0,θlat,θlon,beta)
  α_aero = aoatmp[1]
  θ = aoatmp[2]

  # This file is to get the Cl and Cd of Blade Element
  clcdtmp = clcdget(α_aero,vall_r)
  Cl = clcdtmp[1]
  Cd = clcdtmp[2]

  # This file calculate the whole aerodynamic force of the rotor in the hub
  rftmp = rotorforce(ψ,vall_r,α_aero,θ,Cl,Cd)
  rot = rot+rftmp[3]
  Mbeta_aero = rftmp[5][1]
  # if index%npsi==0
  #   print("===当前挥舞力矩：$(Mbeta_aero)N·m===\n")
  # end
  # print("$(θ/π*180)\n")
  # print("$(α_aero/π*180)\n")
  # print("\n")
  # The file calculate the uniform induced velocity
  uitmp = uniforminflow(ψ,θcp,θlat,θlon,beta,dbeta,0.0)

  ψ = ψ+dψ

  # 此处要输入转过一周进行配平的条件，今天来不及了明天完成，作此标志 12/2/2017
  if index%npsi==0
    rot = rot/npsi
    print("===当前拉力：$(rot)N===\n")
    beta_lat = blat/npsi
    beta_lon = blon/npsi

    # rot = 0.0
    # blat = 0.0
    # blon = 0.0
    ψ = 0.0

    # 此处开始进行配平
    trimtmp = trimwt(uitmp,rot,beta_lat,beta_lon,θcp,θlat,θlon)
    if trimtmp[1]
      # print(trimtmp)
      print("配平总距：$(trimtmp[5]*180/π)\n")
      print("配平横向变距：$(trimtmp[3]*180/π)\n")
      print("配平纵向变距：$(trimtmp[4]*180/π)\n")
      break
    else
      θcp = trimtmp[5]
      θ0 = trimtmp[2]
      print("配平总距：$(trimtmp[5]*180/π)\n")
      θlat = trimtmp[3]
      θlon = trimtmp[4]
      rot = 0.0
      blat = 0.0
      blon = 0.0
    end
  end


  index = index+1
end

print("===计算收敛===\n")
print("===总计算步数$(index-1)===\n")
toc() # 程序结束
