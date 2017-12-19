  # Clear
  const p = 0
  const q = 0
  const r = 0
  const θ7 = 10
  const A1s = 0
  const B1s = 0
  const ρ = 1.225
  const G = 280
  const m = 280
  const Ix = 26
  const Iy = 113
  const Iz = 89
  const Ixy = -1.8
  const Iyz = 0.1
  const Izx = 4.3
  const g = 9.8
  const x = 1.64
  const y = 0
  const z = 0.42

  const R = 2.54
  const Ω = 76.8 # rad/s
  # Ω = 734rpm
  const Ib = 7.806
  const Mb = 5.588
  const a00 = 7.19
  const alpha00 = -0.75
  const b = 0.178
  const n = 2
  const eroot = 0
  const deltw = 0
  const theta1 = -10
  const sigma = 0.045
  const alpha00 = -0.75
  const X_mr = x-1.65
  const Y_mr = 0-y
  const Z_mr = z-1.24
  # 尾桨原始参数
  const theta1_T = 0
  const n = 2
  const ΩT = 470
  const sigmaT = 0.11
  # 4488r/min
  const RT = 0.38
  const a00T = 6.05
  const alpha00T = 0
  const bT = 0.068
  const K_QT = 1
  const X_T = x-4.65
  const Y_T = -0.12-y
  const Z_T = z-0.473
  #平尾原始参数
  const K_MH = 1
  const K_QH = 1
  # v_01 = v_0
  const alpha_sf = 0
  const phi_H = -4
  const A_H = 0.3
  const X_H = x-3.81
  const Y_H = 0
  const Z_H = z-0.4073
  #垂尾原始参数
  #机身原始参数
  const C_DF = 0.2
  const C_LF = 0.02
  const C_SF = 0.06
  const C_MxF = -0.001
  const C_MyF = 0.001
  const C_MzF = -0.01
  const A_F = 0.72
  const l_F = 4.908
  const X_F = 0
  const Y_F = 0
  const Z_F = 0

  Trzb = [0 Z_mr -Y_mr; -Z_mr 0 X_mr; Y_mr -X_mr 0]
  u = [20.0,0.0,0.0]+Trzb*[p,q,r]
  Vx = u[1]
  Vy = u[2]
  Vz = u[3]
  q1 = 1/2*ρ*π*R^2*(Ω*R)^2
  μ = sqrt(Vx*Vx+Vy*Vy)/Ω/R
  lamda00 = -Vz/Ω/R
  theta0 = θ7+7
  if abs(Vx)<=1e-3&&abs(Vy)<=1e-3
      alphas = π/2
  else
      alphas = atan(Vz/sqrt(Vx*Vx+Vy*Vy))
  end

  lamda_b = (ρ*a00*b*R^4)/Ib
  vO0 = 0.000001
  vOK = 1
  vOT = 0.000001
  eM = 0.92*a00*sigma*(theta0/57.3*(1/3+μ^2/2)+
                          theta1/57.3*(1+μ^2)/4-
                            (vOT-lamda00)/2-μ*B1s/57.3/2)       # CT旋翼拉力系数
  fM = 4*sqrt(μ^2+(vOT-lamda00)^2)
  FvOT = vOT-eM/fM
  i = 0
  while abs(FvOT)>0.000001
      if i != 0
          if FvOT>0
              vOK = vOT
          else
              vO0 = vOT
          end
      end
      i = 1
      aM = 0.92*a00*sigma*(theta0/57.3*
                            (1/3+μ^2/2)+theta1/57.3*
                              (1+μ^2)/4-(vO0-lamda00)/2-μ*B1s/57.3/2)
      bM = 4*sqrt(μ^2+(vO0-lamda00)^2)
      cM = 0.92*a00*sigma*(theta0/57.3*
                            (1/3+μ^2/2)+theta1/57.3*
                              (1+μ^2)/4-(vOK-lamda00)/2-μ*B1s/57.3/2)
      dM = 4*sqrt(μ^2+(vOK-lamda00)^2)
      FvO0 = vO0-aM/bM
      FvOK = vOK-cM/dM
      vOT = vOK-FvOK*(vOK-vO0)/(FvOK- FvO0)
      eM = 0.92*a00*sigma*(theta0/57.3*
                            (1/3+μ^2/2)+theta1/57.3*
                              (1+μ^2)/4-(vOT-lamda00)/2-μ*B1s/57.3/2)
      fM = 4*sqrt(μ^2+(vOT-lamda00)^2)
      FvOT = vOT-eM/fM
  end

  v_0 = vOT
  v_01 = v_0*Ω*R
  C_T = eM
  beta0 = 3/8*lamda_b*C_T/a00/sigma-1.5*9.8*R/((Ω*R)^2)
  beta00 = beta0*57.3
  a1 = (2*μ/(1-0.5*μ*μ)*
        (4/3*theta0/57.3+theta1/57.3+μ*alphas-v_0/Ω/R)-
          (1+1.5*μ*μ)*B1s/57.3/(1-0.5*μ*μ))
  a11 = a1*57.3
  b1 = 1/(1+0.5*μ*μ)*(4/3*μ*beta0+v_0/Ω/R)+A1s/57.3
  b11 = b1*57.3
  xt = R
  xs = 0.1*R
  rstep = (xt-xs)/20
  r1 = Array{Float64}(21)
  for j in 1:21
    r1[j] = xs+(j-1)*rstep
  end
  # r1 = xs:rstep:xt
  ψ = linspace(0,2*π,12)
  Ts1 = 0
  Hs1 = 0
  Ss1 = 0
  M1 = 0
  P1 = 0
  Cd1 = Array{Float64}(11,20)
  DQ1 = Array{Float64}(11,20)
  DSs1 = Array{Float64}(11,20)
  Ss11 = Array{Float64}(11)
  Hs11 = Array{Float64}(11)
  for i = 1:11
      Ts = 0
      Hs = 0
      Ss = 0
      M = 0
      P = 0
      sp = sin(ψ[i])
      cp = cos(ψ[i])
      beta = beta0-a1*cp-b1*sp
      beta_d = Ω*(a1*sp-b1*cp)
      beta_dd = Ω*Ω*(a1*cp+b1*sp)
      for j = 1:20
          rad = r1[j]
          theta = θ7+(rad/R-0.7)*theta1-A1s*cp-B1s*sp
          theta = theta/57.3
          ut = Ω*R*(rad/R+μ*sp)
          up = -(Ω*R*(-lamda00-μ*beta*cp)-v_0*Ω*R*(1+rad/R*cp)-rad*beta_d)
          alpha = atan(up/ut)
          alpha11 = theta-alpha-alpha00/57.3
  #         mach = sqrt(up^2+ut^2)/340
  #         CLCD0 = CLCD(alpha11,mach)
  #         CL0 = CLCD0 [1]
  #         CD0 = CLCD0[2]
          CL = a00*alpha11
          Cd = 0.008-0.003*CL+0.01*CL*CL
          Cd1[i,j] = Cd
          q2 = 0.5*ρ*(up^2+ut^2)
          Dy = q2*CL*b
          Dx = q2*Cd*b
          DT = Dy*cos(alpha)-Dx*sin(alpha)
          DQ = Dx*cos(alpha)+Dy*sin(alpha)
          DQ1[i,j] = DQ
          DTs = DT*cos(beta)
          DHs = DQ*sin(ψ[i])-DT*sin(beta)*cos(ψ[i])
          DSs = -DQ*cos(ψ[i])-DT*sin(beta)*sin(ψ[i])
          DSs1[i,j] = DSs
          DM = DQ*rad
          DP = DQ*rad*Ω
          Ts = Ts+DTs*rstep
          Hs = Hs+DHs*rstep
          Ss = Ss+DSs*rstep
          M = M+DM*rstep
          P = P+DP*rstep
      end
      Ts1 = Ts1+Ts
      Hs1 = Hs1+Hs
      Ss1 = Ss1+Ss
      Ss11[i] = Ss
      Hs11[i] = Hs
      M1 = M1+M
      P1 = P1+P
  end
  Ts = Ts1*n/11
  Hs = Hs1*n/11
  Ss = Ss1*n/11
  M = M1*n/11
  P = P1*n/11
  Mx = 0.5*n*eroot*Mb*Ω^2*b1
  My = 0.5*n*eroot*Mb*Ω^2*a1
  Trzb1 = [cos(deltw) -sin(deltw) 0; 0 0 1; sin(deltw) cos(deltw) 0]
  #旋翼侧滑角
  betaw = asin(Vy/sqrt(Vx*Vx+Vy*Vy))
  Ts = Ts
  Hs = Hs*cos(betaw)+Ss*sin(betaw)
  Ss = -Hs*sin(betaw)+Ss*cos(betaw)
  # v_0 = v_0*Ω*R
  beta0 = beta0*57.3
  a1 = a1*57.3
  b1 = b1*57.3
  # C_T = C_T0
  F = Trzb1*[-Hs,-Ts,-Ss]
  M = Trzb1*[Mx,M,My]+-1*Trzb*F
  FM = [F,M]
  FX = FM
  Cd1 = Cd1
  DQ1 = DQ1
  DSs1 = DSs1
  Ss11 = Ss11
  Hs11 = Hs11
