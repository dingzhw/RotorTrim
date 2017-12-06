# 本代码主要是针对此项目制定一个通用命名标准

# 桨叶的物理描述
# 本项目中旋翼转向为逆时针旋转，也就是俗称右旋旋翼
# 旋翼一般都是沿径向扭转的，且在实际分析中多为负扭转，Δθ=θtw*r，常规直升机一般为负扭转
R = 10  #the rotor radius; the length of the blade, measured from the rotation to tip
Omega = 66  #the rotor ratational speed or angular velocity
rho = 1.25  #air density
psi = 0 #azimuth angle of the blade; for constant rotational speed, ψ = Ω*t
r = 0.1 #radial location on the blade, measured from the center of ratation(r = 0)
ch = 0.1 #blade chord
NR = 2 #number of blades
Nbe = 10 #number of blade elements
m = 0.3 #blade mass per unit length
Ib = 10 #∫₀ᴿ(r^2*m); moment of inertia of the blade about the center of rotation
A = π*R^2 #rotor disk area
sigma = Nch/(π*R) #rotor solidity
gama = rho*cla*ch*R^4/Ib #blade lock number (洛克数代表了气动力和惯性力的比值)

# 桨叶气动参数
cla = 0.6 #blade section two-dimensional lift curve slope
alpha = 0.5 #blade section angle-of-attack
Ma = 0.3 #blade section Mach number

# 桨叶运动参数
# 在本项目中，桨叶运动被认为是绕桨毂的刚体运动
beta = 0.5 #blade flap angle
zeta = 0.5 #blade lag angle
theta = 0.5 #blade pitch angle
thetatw = -5 #blade twist
beta = beta0+betalon*cos(psi)+betalat*sin(psi)
theta = theta0+thetalat*cos(psi)+thetalon*cos(psi)+thetatw*r/R

# 旋翼迎角和速度
# 下标s表示旋翼固定坐标系
alphai = 0.5 #rotor disk plane incidence angle(低头为正)
vair = 20 #rotor or helicopter velocity with respect to the air
vind_s = 10 #rotor induced velocity, normal to the disk plane
mu = vair*cos(alphai)/(Omega*R) #rotor advance ratio
muz = vair*sin(alphai)/(Omega*R) #normal velocity ratio
lambda = (vair*sin(alphai)+vind_s)/(Omega*R) #rotor inflow ratio
lambdai = vind_s/(Omega*R) #induced inflow ratio

# 旋翼力和力矩
T = 100 #rotor thrust
H = 10 #rotor drag force in the disk plane
Y = 1 #rotor side force in the disk plane
Q = 5 #rotor shaft torque
P = 10 #rotor shaft power(P=Ω*Q)
CT = T/(rho*A*(Omega*R)^2)
CH = H/(rho*A*(Omega*R)^2)
CY = Y/(rho*A*(Omega*R)^2)
CQ = Q/(rho*A*(Omega*R)^2*R)
CP = P/(rho*A*(Omega*R)^3) #CP=CQ

#  旋翼桨盘平面
# tpp tip-path plane
# nfp no-feathering plane
# hp hub plane
# cp control plane
#
#  各坐标系下标
# s system(固定坐标系)
# r rotation(旋转坐标系)
# β beta(挥舞坐标系)
# be blade element(叶素当地坐标系)

# 其他下标说明
# all means total
# x means x coordination
# y means y coordination
# z means z coordination
# aero表示该项是由气动导致的
# ui表示该项是均匀入流
# ind代表诱导速度
