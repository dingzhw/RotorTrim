# This file solve the flap reponses of blades
function sflapre(ψ,λ_α,θ_cp,θ_lat,θ_lon) #staticflapre
  # 其中θ应当输入平均值，例如mean(θ_cp)
  γ_ = (ρ*5.73*(0.06)*R^4)/Iβ
  μ = μ_air
  F0 = 1+μ^2/2
  F_B1 = μ
  F_λ = 1
  A0 = 2*μ
  A_B1 = 1+3*μ^2/4
  A_λ = μ
  A_a1 = 1-μ^2/4
  B_β0 = μ
  B_A1 = 1+μ^2/4
  B_b1 = 1+μ^2/4
  γ_B = γ_/2*B_β0/B_b1
  θ0_ = θcp[1]
  θ_lat_ = θ_lat[1]
  θ_lon_ = θ_lon[1]
  Matrix_β = [γ_/2*F0 0 γ_/2*F_B1 γ_/2*F_λ;
              A0/A_a1 0 A_B1/A_a1 A_λ/A_a1;
              γ_B*F0 -B_A1/B_b1  γ_B*F_B1 γ_B*F_λ]*[θ0_,θ_lat_,θ_lon_,λ_α]
  β0 = Matrix_β[1]
  β_lon = Matrix_β[2]
  β_lat = Matrix_β[3]
  β = β0-β_lon*cos(ψ)-β_lat*sin(ψ)
  dβ = Ω*(β_lon*sin(ψ)-β_lat*cos(ψ))
  ddβ = Ω*(β_lon*cos(ψ)+β_lat*cos(ψ))
  return β,dβ,ddβ,β_lon,β_lat,β0
end

function dflapre(betanow,betaxnow,Mbeta_aero)	#动态挥舞响应
  # 输入：当前挥舞角、挥舞角速度
  # 输出：步进挥舞角、挥舞角速度
	# 常规参数
	Mcen = ecut*R*Ω^2*m_*1/2*(R-ecut*R)^2
	MG = m_*g*1/2*(R-ecut*R)^2
  nite = 10 # 迭代步数
  h = dψ/nite # 龙格库塔求解间隔
  ite = 1
  betanext = 0.0
  betaxnext = 0.0

	# 挥舞方程
	funcbetax(xxbeta) = 1/(Iβ*Ω^2)*(Mbeta_aero-MG*cos(xxbeta)-
											Mcen*sin(xxbeta))-
											cos(xxbeta)*sin(xxbeta)

	# 龙格库塔求解
  while ite<=nite
  	k1 = betaxnow
  	kx1 = funcbetax(betanow)
  	k2 = betaxnow+h/2*k1
  	kx2 = funcbetax(betanow+h/2*kx1)
  	k3 = betaxnow+h/2*k2
  	kx3 = funcbetax(betanow+h/2*kx2)
  	k4 = betaxnow+h*k3
  	kx4 = funcbetax(betanow+h*kx3)
  	betanext = betanow + h/6*(k1+2*k2+2*k3+k4)
  	betaxnext = betaxnow + h/6*(kx1+2*kx2+2*kx3+kx4)
    betanow = betanext
    betaxnow = betaxnext
    ite = ite+1
  end

	return betanext,betaxnext
end
