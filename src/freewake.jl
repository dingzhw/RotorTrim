#This file is the Free-Wake code

function rigidwake(λ,filename)
  pointfile = open(filename,"w")
  write(pointfile,"  TITLE = 'RIGID WAKE OF ROTOR'\n")
  ktip = 0.88 #Default Tip Vortex ScrollUp Position
  # dt = dψ/Ω
  Nrw = Int64(2.0*π/dψ)
  rwps = Array{Vector}(Nrw,Nb,NR)
  rbp = Array{Float64}(NR,Nb)
  ζ0 = Array{Float64}(Nrw,Nb,NR)
  # Rtip = Array{Float64}(NR)
  for k in 1:NR
    # Rtip[k] = ktip*R
    for j in 1:Nb
      write(pointfile,"  VARIABLES = 'X','Y','Z'\n")
      write(pointfile,"  ZONE I = $(Nrw)\tF = POINT\n")
      rbp[k,j] = ecut+R/Nb*(j-1)
      for i in 1:Nrw
        ψ = (i-1)*dψ
        ψk = ψ+(k-1)*2*π/NR
        ζ0[i,j,k] = -ψ
        rwps[i,j,k] = [rbp[k,j]*cos(ψk)+μ*R*ψk;
                      -rbp[k,j]*sin(ψk);
                      λ*R*ψ]
        write(pointfile,
              "\t$(rwps[i,j,k][1])\t$(rwps[i,j,k][2])\t$(rwps[i,j,k][3])\n")
      end
    end
  end
  close(pointfile)
  return rwps,ζ0
end

function Gama(Cl,V∞,chord)
  gama = abs(0.5*Cl*V∞*chord)
  return gama
end

function freewake(wp0,Nwp,ψ,ζ0,Cl,v_all)
  #wake parameters
  ktip = 0.88
  dt = dψ/Ω
  rcut = 2*R
  # Nwp = Int64(2*π/dψ)
  bp = Array{Vector}(NR,Nb+1)
  for k in 1:NR
    ψk = ψ+2*π/NR*(k-1)
    for i in 1:(Nb+1)
      bp[k,i] = [(ecut+R/NR*(i-1))*cos(ψk),(ecut+R/NR*(i-1))*sin(ψk),0]
    end
  end

  Rtip = Array{Float64}(NR)
  for k in 1:NR
    Rtip[k] = ktip*R
  end

  #begin :: solve the Γ of wakepoints
  gama_bw = Array{Float64}(NR,Nb)
  for k in 1:NR
    for i in 1:Nb
      gama_bw[k,i] = Gama(Cl[k,i],v_all[k,i][2],c[k,i])
    end
  end
  #end :: the Γ of wakeponits solved

  #begin :: solve the Γ of wakepoints
  gama_wp = Array{Float64}(NR,Nb+1)
  for k in 1:NR
    gama_wp[k,1] = 0-gama_bw[k,1]
    gama_wp[k,Nb+1] = gama_bw[k,Nb]-0.0
    for i in 2:Nb
      gama_wp[k,i] = gama_bw[k,i-1]-gama_bw[k,i]
    end
  end
  #end :: the Γ of wakepoints solved


  # gama_ = sum(gama_bw)/NR # average of whole gama
  # vgamaind = Array{Vector}(Nb,Nwp)

  #begin :: initialize rigid wake
  # wp0 = rigidwake(uniforminflow(ψ)[3],"RIGIDWAKE.PLT")[1]
  #end :: rigid wake initialization

  #begin :: calculate bound vortex induced influnce on wakepoints
  vgamaind = Array{Vector}(Nwp,Nb,NR)
  for k in 1:NR
    for l in 1:Nb
      for i in 1:Nwp
        vgind = [0.0,0.0,0.0]
        for m in 1:NR
          for n in 1:Nb
            vgind += bslinevnocore(gama_bw[m,n],bp[m,n],bp[m,n+1],wp[i,l,k])
          end
        end
        vgamaind[i,l,k] = vgind
      end
    end
  end
  #end :: bound vortex induced influnce solving

  #begin :: calculate wake vortex influnce on wakepoints
  vwakeind = Array{Vector}(Nwp,Nb,NR)
  for k in 1:NR
    for l in 1:Nb
      for i in 1:Nwp
        vwvind = [0.0,0.0,0.0]
        for m in 1:NR
          for n in 1:Nb
            for r in 1:(Nwp-1)
              t_wp = abs(ψ-ζ0[r,n,m])/Ω
              vwvind += bslinev(gama_wp[m,n],
                                wp[r,n,m],wp[r+1,n,m],wp[i,l,k],t_wp)
              # vwvind += bslinevnocore(gama_wp[m,n],
              #                   wp[r,n,m],wp[r+1,n,m],wp[i,l,k])
            end
          end
        end
        vwakeind[i,l,k] = vwvind
      end
    end
  end
  #end :: wake vortex influnce on wakepoints calculated

  #begin :: wakepoints move with its local velocity
  vwpsind = Array{Vector}(Nwp,Nb,NR)
  vwps = Array{Vector}(Nwp,Nb,NR)
  dwps = Array{Vector}(Nwp,Nb,NR)
  for k in 1:NR
    for i in 1:Nb
      for j in 1:Nwp
        vwpsind[j,i,k] = vwakeind[j,i,k]+vgamaind[j,i,k]
        vwps[j,i,k] = vwpsind[j,i,k]+v_for
        dwps[j,i,k] = vwps[j,i,k]*dt
        wp[j,i,k] = wp[j,i,k]+dwps[j,i,k]
        ζ0[j,i,k] = ζ0[j,i,k]+dψ
      end
    end
  end
  #end :: wakepoints moved

  return wp,gama_wp
end

function wpgen(ψ,ζ0,wp,Nwp)  #new wakepoints generate
  l = 1
  if norm(ψ)<=1e-3
    return wp,ζ0
  else
    wpsadd = Array{Vector}(l,Nb,NR)
    ζadd = Array{Float64}(l,Nb,NR)
    rbp = Array{Float64}(NR,Nb)
    for k in 1:NR
      ψk = ψ+2*pi/NR*(k-1)
      for i in 1:Nb
        for l in 1
          rbp[k,i] = ecut+R/Nb*(i-1)
          wpsadd[l,i,k] = [rbp[k,i]*cos(ψk),
                            rbp[k,i]*sin(ψk),
                            0.0]
          ζadd[l,i,k] = ψ
        end
      end
    end
    newwp = maprepend(wp,wpsadd,Nwp,Nb,NR)
    newζ0 = maprepend(ζ0,ζadd,Nwp,Nb,NR)
    return newwp,newζ0
  end
end

function vbe(ψ,ζ0,gama_wp,wp) #calculate velocity of blade elemnts(both sys and ro coordinations)
  bp = Array{Vector}(NR,Nb+1)
  for k in 1:NR
    ψk = ψ+2*π/NR*(k-1)
    for i in 1:(Nb+1)
      bp[k,i] = [(ecut+R/NR*(i-1))*cos(ψk),(ecut+R/NR*(i-1))*sin(ψk),0]
    end
  end
  vb_ind = Array{Vector}(NR,Nb)
  for k in 1:NR
    for l in 1:Nb
      vbind = [0.0,0.0,0.0]
      for m in 1:NR
        for n in 1:Nb
          for r in 1:(Nwp-1)
            t = abs(ψ-ζ0[r,n,m])/Ω
            vbind += bslinev(gama_wp[m,n],wp[r,n,m],wp[r+1,n,m],bp[k,l],t)
          end
        end
      end
      vb_ind[k,l] = vbind
    end
  end

  rbp = Array{Float64}(NR,Nb) #radius of balde element points
  v_ro = Array{Vector}(NR,Nb)
  vb_sys = Array{Vector}(NR,Nb)
  vb_ro = Array{Vector}(NR,Nb)
  for k in 1:NR
    ψk = ψ+2*π/NR*(k-1)
    for i in 1:Nb
      rbp[k,i] = ecut+R/Nb*(i-1)
      v_ro[k,i] = [Ω*rbp[k,i]*sin(ψk),
                    -Ω*rbp[k,i]*cos(ψk),
                    0.0]
      vb_sys[k,i] = v_ro[k,i]+v_for+vb_ind[k,i]
      vb_ro[k,i] = systoro(vb_sys[k,i],ψ)
    end
  end
  return vb_sys,vb_ro
end
