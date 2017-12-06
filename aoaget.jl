# This file is used for get angle-of-attack
# be carful of the coordination

function aoaget(vall_r,ψ,θ0,θ_lat,θ_lon,beta=0.0)
#The solution is directly solving the angle between
#pitch angle vector and wind angle vector.

  α_aero = Array{Float64}(NR,Nbe)
  θ = Array{Float64}(NR,Nbe)

  for k in 1:NR
    ψk = ψ+(k-1)*2*π/NR
    for i in 1:Nbe
      v_per = vall_r[k,i][3]
      v_tan = vall_r[k,i][2]
      v_rad = vall_r[k,i][1]
      θ[k,i] = θ0[k,i]+θ_lat[k]*cos(ψk)+θ_lon[k]*sin(ψk) #pitch angle

      wind_v = [0,v_tan,v_per] #wind velocity
      wind_ = wind_v/norm(wind_v) #unit vector of wind velocity
      wind_tmp = rotobe(wind_,beta,θ[k,i])
      wind_be = [wind_tmp[2],wind_tmp[3]] # 叶素当地二维气流速度
      α_aero[k,i] = aoaang(wind_be) # 最终的气动迎角解
    end
  end
  return α_aero,θ #,dF_be
end
