# The main file of the project

include("const.jl")
include("wakemath.jl")
include("bioslinevortex.jl")
include("flapreponse.jl")
include("uniforminflow.jl")
include("freewake.jl")
include("aoaget.jl")
include("clcdget.jl")
include("rotorforce.jl")

# 桨叶安装角和总距角初始化
θcp = Array{Float64}(NR)  #总距角
θroot = Array{Float64}(NR)  #桨根安装角
tw_θ = Array{Float64}(NR) #桨叶扭转
θ_lon = Array{Float64}(NR)  #纵向周期变距
θ_lat = Array{Float64}(NR)  #横向周期变距

for k in 1:NR
  θcp[k] = 10.0*π/180
  θroot[k] = 50.0*π/180
  tw_θ[k] = 1.0*π/180
  θ_lon[k] = 0.0
  θ_lat[k] = 0.0
end

ψ = 0.0
index = 1
sumforce = 0.0

# The file calculate the uniform induced velocity
v_res = uniforminflow(ψ)
v_ind = v_res[1]
v_all = v_res[2]

#begin :: Calculate the rigid wake
riwake = rigidwake(v_res[3],"RIGIDWAKE.PLT")
wp = riwake[1]
ζ0 = riwake[2]
Nwp = trunc(Int,length(wp)/(NR*Nbe))
#end :: Rigid wake calculated

α_aero = aoaget(v_all,ψ)[1]
Cl = clcdget(α_aero,v_all)[1]

while ψ<=(20*π)

  wp_gen = wpgen(ψ,ζ0,wp,Nwp)
  wp = wp_gen[1]
  Nwp = trunc(Int,length(wp)/(NR*Nbe))
  ζ0 = wp_gen[2]
  fwake = freewake(wp,Nwp,ψ,ζ0,Cl,v_all)
  gama_wp = fwake[2]
  wp = fwake[1]
  vbelm = vbe(ψ,ζ0,gama_wp,wp)
  v_all = vbelm[2]

  # This file is used for get angle-of-attack

  α_aero = aoaget(v_all,ψ)[1]
  θ = aoaget(v_all,ψ)[2]

  # This file is to get the Cl and Cd of Blade Element

  Cl = clcdget(α_aero,v_all)[1]
  Cd = clcdget(α_aero,v_all)[2]

  # This file calculate the whole aerodynamic force of the rotor in the hub
  sumforce += rotorforce(v_all,α_aero,θ,Cl,Cd,
                         pwd()*"\\force_record\\bladeforce_$(index).txt")[2]

    ψ += dψ
    index += 1
end

wpplt = open("wakepoints.PLT","w")
write(wpplt,"  TITLE = 'RIGID WAKE OF ROTOR'\n")
for k in 1:NR
  for j in 1:Nbe
    write(wpplt,"  VARIABLES = 'X','Y','Z'\n")
    write(wpplt,"  ZONE I = $(Nwp)\tF = POINT\n")
    for i in 1:Nwp
      write(wpplt,
            "\t$(wp[i,j,k][1])\t$(wp[i,j,k][2])\t$(wp[i,j,k][3])\n")
    end
  end
end
close(wpplt)

print("The thrust of the rotor is : $(sumforce/index)")
