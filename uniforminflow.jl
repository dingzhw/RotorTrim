# The file calculate the uniform induced velocity
# 给出均匀入流模型

function uniforminflow(ψ,θcp,θlat,θlon,betanow,betaxnow,Mbeta_aero=0.0)
  vind_s = [0.0,0.0,-T/(2*ρ*A*norm(v_air))] # 固定坐标系诱导速度初值
  vall_s = v_air + vind_s # 固定坐标系总的滑流速度初值
  # lmdaui = -T/(2*ρ*A*norm(vall_s))/(Ω*R) # 前飞均匀诱导速度系数
  vindui_s = [0.0,0.0,-T/(2*ρ*A*norm(vall_s))]  # 第一步迭代的固定坐标系均匀入流值

  while true  # 迭代求解均匀入流
    if norm(vindui_s-vind_s)<=1e-3
      break
    else
      vind_s = vindui_s
      vall_s = v_air + vind_s
      vindui_s = [0.0,0.0,-T/(2*ρ*A*norm(vall_s))]
    end
  end
  lmdaui = vindui_s[3]/(Ω*R) # 收敛前飞均匀诱导速度系数

  vinduitmp_r = Array{Vector}(NR,Nbe)
  valluitmp_r = Array{Vector}(NR,Nbe)

  # # 准静态挥舞解方法
  # β = sflapre(ψ,vall_s[3]/(Ω*R),θcp,θlat,θlon)[1]
  # dβ = sflapre(ψ,vall_s[3]/(Ω*R),θcp,θlat,θlon)[2]
  # vβ_β = [0.0,0.0,-dβ] # 挥舞坐标系下挥舞引起的气流速度分量
  # vβ_r = betatoro(vβ_β,β)


  if Mbeta_aero==0.0 # 入流场无挥舞初始化
    vβ_r = [0.0,0.0,0.0]
    β = 0.0
    dβ = 0.0
  else
    # 时间步进挥舞解方法
    betatmp = dflapre(betanow,betaxnow,Mbeta_aero)
    β = betatmp[1]
    dβ = betatmp[2]
    vβ_β = [0.0,0.0,-dβ] # 挥舞坐标系下挥舞引起的气流速度分量
    vβ_r = betatoro(vβ_β,β)
  end

  for k in 1:NR
    ψk = ψ+(k-1)*2*π/NR
    for i in 1:Nbe
      # 将速度从固定坐标系转为旋转坐标系
      vinduitmp_r[k,i] = systoro(vindui_s,ψk)
      vOmega_r = [0,-Ω*rb[k,i],0] # 旋转坐标系下旋转导致的气流分量
      valluitmp_r[k,i] = systoro(v_air,ψk)+vinduitmp_r[k,i]+vOmega_r+vβ_r*(rb[k,i]-ecut*R)
    end
  end

  return vinduitmp_r,valluitmp_r,lmdaui,β,dβ
end
