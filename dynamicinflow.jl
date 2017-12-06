# This is a version of Dynamic Inflow from the C++ version
# The auther is Dingzw and all rignts reserved
# The starting time is 9/23/2017
# The restarting time is 11/11/2017

# Dynamic Inflow Codes for Rotor

function staticdynamicinflow(v,Mx,My) #Steady Pitt-peter dynamic Inflow
  # v::前飞来流速度
  # Mx,My::桨毂力矩

  L = Array{Float64}(3,3) #初始化动态入流求解矩阵

  # 前飞来流速度的三个分量
  μ1 = v[1]/(Ω*R)
  μ2 = v[2]/(Ω*R)
  μ3 = v[3]/(Ω*R)

  CT = T/fnonc
  C1 = -Mx/mnonc
  C2 = My/mnonc

  λ_m = uniforminflow(ψ)[3]
  λ = λ_m - μ3
  μ = √(μ1^2+μ2^2)

  if √(μ1^2+μ2^2)<0.05  #Hovering
    λ_0 = λ_m
    λ_c = 0.0
    λ_s = 0.0
  else  #Forward Flight
    salfa = sin(atan(abs(λ)/μ))
    B = 15.0*π*√((1.0-salfa)/(1.0+salfa))/64.0
    D = 4.0*salfa/(1.0+salfa)
    E = 4.0/(1.0+salfa)
    sdelta = μ2/μ
    cdelta = μ1/μ
    VT = √(λ^2+μ^2)
    V = (μ^2+(2.0*λ_m-μ3)*(λ_m-μ3))/VT

    L[1,1] = 0.5/VT
    L[1,2] = -B*sdelta/V
    L[1,3] = -B*cdelta/V
    L[2,1] = B*sdelta/VT
    L[2,2] = (E*cdelta^2+D*sdelta^2)/V
    L[2,3] = (D-E)*sdelta*cdelta/V
    L[3,1] = B*cdelta/VT
    L[3,2] = (D-E)*cdelta*sdelta/V
    L[3,3] = (E*cdelta*cdelta+D*sdelta^2)/V

    λ_0 = L[1,1]*CT+L[1,2]*C1-L[1,3]*C2
    λ_s = L[2,1]*CT+L[2,2]*C1-L[2,3]*C2
    λ_c = L[3,1]*CT+l[3,2]*C1-L[3,3]*C2
  end

  # 有量纲的轴向入流
  v_axi = [λ_0*Ω*R,λ_c*Ω*R,λ_s*Ω*R]

  return v_axi

end
