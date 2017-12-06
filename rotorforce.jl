# This file calculate the whole aerodynamic force of the rotor in the hub

function beforce(vall_r,chord,Cl,Cd,α_aero,θ,dr,beta=0.0)
  # vall_r = Array{Vector}(NR,Nbe)
  # for k in 1:NR # 将总速度由旋转坐标系转化到叶素坐标系
  #   for i in 1:Nbe
  #     vall_r[k,i] = rotobe(vall_r[k,i],beta,θ[k,i])
  #   end
  # end
  L = abs(1/2*ρ*norm(vall_r)^2*chord*Cl)
  D = abs(1/2*ρ*norm(vall_r)^2*chord*Cd)
  if α_aero>=0
    L_be = [0,cos(α_aero+π/2),-sin(α_aero+π/2)] # 叶素当地升力单位矢量
  else
    L_be = [0,cos(α_aero-π/2),-sin(α_aero-π/2)]
  end
  D_be = [0,cos(α_aero),-sin(α_aero)] # 叶素当地阻力单位矢量
  L_ro = betoro(L_be,beta,θ) # 旋转坐标系叶素升力单位矢量
  D_ro = betoro(D_be,beta,θ) # 旋转坐标系叶素阻力单位矢量
  Fbez = (L*L_ro[3]+D*D_ro[3])*dr # 旋转坐标系叶素z向力
  Fbey = (L*L_ro[2]+D*D_ro[2])*dr # 旋转坐标系叶素y向力
  return Fbez,Fbey
end

function rotorforce(ψ,vall_r,α_aero,θ,Cl,Cd,filename=pwd()*"\\force_record\\fhub.txt")
  # filebf = open(filename,"w")
  # write(filebf,"List  ZForce  YForce\n")

  Fz_r = 0.0 # 旋转坐标系轴向力
  Fy_r = 0.0 # 旋转坐标系y向力
  MQ = 0.0 # 扭矩
  Mβ_aero = zeros(Float64,NR) # 初始化气动力挥舞力矩
  for k in 1:NR
    # write(filebf,"$(k)\t\n")
    for i in 1:Nbe
      fbe = beforce(vall_r[k,i],ch[k,i],Cl[k,i],Cd[k,i],
                      α_aero[k,i],θ[k,i],dr[k,i])
      fbez = fbe[1]
      fbey = fbe[2]
      MQ += fbey*rb[k,i]
      Mβ_aero[k] += fbez*(rb[k,i]-ecut*R)
      Fz_r += fbez
      Fy_r += fbey
      # write(filebf,"$(i)  $(Fzr)  $(Fyr)\n")
    end
  end

  # write(filebf,"The Whole Force of Z is : $(Fzr)\n")
  # write(filebf,"The Whole Force of Y is : $(Fyr)\n")
  #
  # close(filebf)
  Fy_s = Fy_r*cos(ψ)
  Fx_s = -Fy_r*sin(ψ)
  Fz_s = Fz_r

  return Fx_s,Fy_s,Fz_s,MQ,Mβ_aero
end
