# This script is written for calculation the Hover Result of Helicopter Rotor

@everywhere function ui2hover(ψ,rb,T_hover=T,A_hover=A)
    # this function is for solving uniform induced velocity and total velocity of
    # blade elements in rotation coordination in hover condition
    vind_hover = [0.,0.,-sqrt(T_hover/2/ρ/A_hover)] # 动量理论求解悬停诱导速度

    # ---Initialize Induced Velocity Array & Total Velocity Array---
    vindhover_r = Array{Vector}(NR,Nbe)
    vallhover_r = Array{Vector}(NR,Nbe)

    # ---Calculate Induced Velocity and Total Velocity in rotation coordination---
    for k in 1:NR
        ψk = ψ+(k-1)*2*π/NR
        for i in 1:Nbe
            vindhover_r[k,i] = systoro(vind_hover,ψk)
            vOmega_r = [0.,-Ω*rb[k,i],0.]
            vallhover_r[k,i] = vindhover_r[k,i]+vOmega_r
        end
    end

    return vindhover_r,vallhover_r,vind_hover[3]
end

@everywhere function betaang(θcp,θlon,θlat,λ,γ)
    # calculate current flap angle
    beta0 = γ/2*θcp+γ/2*λ
    betalon = θlon
    betalat = -θlat
    return beta0,betalon,betalat
end

@everywhere function binaryse(rot, θcp, dθ, T_hover=T, simθ=10.)
    if abs(rot-T_hover)<simθ
        return true, θcp, dθ
    elseif rot<T_hover
        θcp = θcp-dθ
        return false, θcp, dθ/2.
    else
        θcp = θcp+dθ
        return false, θcp, dθ/2.
    end
end

@everywhere function uitest2hover(θ7,thelat,thelon,rb)
    # this function is for running the whole hovering calculation of rotor

    ψ = 0.0
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
    betalon = thelon
    betalat = -thelat
    index = 1
    trimind = 1
    rot = 0.
    MQ  = 0.
    dθ  = 0.2*π/180
    trimrot = zeros(Float64,2)

    # ---Initialize Induced Velocity Field---
    uitmp = ui2hover(ψ,rb)

    while true  # trim iteration
        vindhover_r = uitmp[1]
        vallhover_r = uitmp[2]
        λ           = uitmp[3]

        # ---calculate flap coefficient
        betatmp = betaang(θcp[1],thelon,thelat,λ,(ρ*8.2*(0.062)*R^4)/Iβ)
        β0      = betatmp[1]
        βlon    = betatmp[2]
        βlat    = betatmp[3]

        # ---calculate current flap angle
        betaval = β0-βlon*cos(ψ)-βlat*sin(ψ)

        # ---calculate the a.o.a of every blade element---
        aoatmp = aoaget(vallhover_r,ψ,θ0,θlat,θlon,betaval)
        α_aero = aoatmp[1]
        θ      = aoatmp[2]

        # ---calculate the lift&drag coefficient of each blade element
        clcdtmp = clcdget(α_aero,vallhover_r)
        cl = clcdtmp[1]
        cd = clcdtmp[2]

        # calculate the whole aerodynamic force and moment of the rotor in hover
        rftmp = rotorforce(ψ,vallhover_r,α_aero,θ,cl,cd)
        rot = rot+rftmp[3]
        MQ  = MQ+rftmp[4]
        ψ += dψ

        # judge if trim condition is satisfied
        if index%npsi==0
            rot = rot/npsi
            if (trimind-1)%3!=0
                trimrot[(trimind-1)%3] = rot
            else
                trimrot0 = rot
            end
            MQ  = MQ/npsi
            print("---The Trust now is : $(rot) N\n")
            if abs(rot-T)<10||index>=1000*npsi
                print("---The collective Pitch is : $(θcp/π*180)\n")
                break
            else
                print("\n$(trimind)\n")
                if trimind%3 == 1
                    for k in 1:NR
                        θcp[k] = θcp[k]+dθ
                        for i in 1:Nbe
                            θ0[k,i] = θcp[k]+θtw*(rb[k,i]/R-0.7)
                        end
                    end
                elseif trimind%3 == 2
                    for k in 1:NR
                        θcp[k] = θcp[k]-dθ
                        for i in 1:Nbe
                            θ0[k,i] = θcp[k]+θtw*(rb[k,i]/R-0.7)
                        end
                    end
                else
                    for k in 1:NR
                        θcp[k] = θcp[k]+(T-trimrot0)*2*dθ/(trimrot[1]-trimrot[2])
                        for i in 1:Nbe
                            θ0[k,i] = θcp[k]+θtw*(rb[k,i]/R-0.7)
                        end
                    end
                end
            end

            ψ = 0.  # Initialize the azimuth every circle
            trimind += 1
            rot = 0.
            MQ  = 0.
            # # ---run the binary search---
            # bsetmp = binaryse(rot, θcp[1], dθ)
            # if bsetmp[1]||index==100*npsi
            #     print("---The collective Pitch is : $(θcp/π*180)\n")
            #     break
            # else
            #     dθ = bsetmp[3]
            #     print("\n$(bsetmp)\n\n")
            #     for k in 1:NR
            #         θcp[k] = bsetmp[2]
            #         for i in 1:Nbe
            #             θ0[k,i] = θcp[k]+θtw*(rb[k,i]/R-0.7)
            #         end
            #     end
            #     rot = 0.
            #     MQ  = 0.
            # end
        end

        uitmp = ui2hover(ψ,rb)
        index += 1
    end


    if index>=1000*npsi
        return false,abs(MQ*Ω*10)
    else
        return true,abs(MQ*Ω)
    end
end
