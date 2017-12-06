function maprepend(a,b,i,j,k,l::Int64=1)
  # 多维矢量数组的新插值方法
  # a replace of prepend! for multi-dim vectors
  # i,j,k is the dimension of Array a
  # l,j,k is the dimension of Array b

  if typeof(a)==Array{Array{T,1} where T,3}
    c = Array{Vector}(i+l,j,k)
  elseif typeof(a)==Array{Float64,3}
    c = Array{Float64}(i+l,j,k)
  else
    print("Maprepend return an error message!")
    print("The value input could not be handled!")
    return false
  end

  for m in 1:l
    for n in 1:j
      for r in 1:k
        c[m,n,r] = b[m,n,r]
      end
    end
  end

  for m in (l+1):(i+l)
    for n in 1:j
      for r in 1:k
        c[m,n,r] = a[m-l,n,r]
      end
    end
  end

  return c
end

# 坐标系转换矩阵
function systoro(m1,ψ) #固定坐标系转为旋转坐标系
  m2 = zeros(Float64,3)
  stor = [cos(ψ) sin(ψ) 0;
          -sin(ψ) cos(ψ) 0;
          0 0 1]
  # m2[1] = cos(ψ)*m1[1]+sin(ψ)*m1[2]
  # m2[2] = -sin(ψ)*m1[1]+cos(ψ)*m1[2]
  # m2[3] = m1[3]
  m2 = stor*m1
  return m2
end

function rotosys(m2,ψ) #旋转坐标系转为固定坐标系
  m1 = zeros(Float64,3)
  rtos = inv([cos(ψ) sin(ψ) 0;
              -sin(ψ) cos(ψ) 0;
              0 0 1])
  m1 = rtos*m2
  return m1
end

function rotobeta(mro,beta) # 旋转坐标系到挥舞坐标系
  mbeta = zeros(Float64,3)
  rtb = [cos(beta)  sin(beta) 0;
          0 1 1;
          -sin(beta)  0 cos(beta)]
  mbeta = rtb*mro
  return mbeta
end

function betatoro(mbeta,beta) # 挥舞坐标系到旋转坐标系
  mro = zeros(Float64,3)
  btr = inv([cos(beta)  sin(beta) 0;
              0 1 1;
              -sin(beta)  0 cos(beta)])
  mro = btr*mbeta
  return mro
end

function betatobe(mbeta,θ) # 从挥舞坐标系转换到叶素坐标系
  mbe = zeros(Float64,3)
  betobe = [1 0 0;
            0 -cos(θ) -sin(θ);
            0 sin(θ)  -cos(θ)]
  mbe = betobe*mbeta
  return mbe
end

function rotobe(mr,beta,θ) # 从旋转坐标系转化到叶素坐标系
  mbe = zeros(Float64,3)
  rtob = [cos(beta)   -cos(θ)*sin(beta) -sin(beta)sin(θ);
          0   -cos(θ) -sin(θ);
          -sin(beta)   cos(beta)*sin(θ)  -cos(beta)*cos(θ)]
  mbe = rtob*mr
  return mbe
end

function betoro(mbe,beta,θ) # 从叶素坐标系转换到旋转坐标系
  mr = zeros(Float64,3)
  btor = inv([cos(beta)   -cos(θ)*sin(beta) -sin(beta)sin(θ);
              0   -cos(θ) -sin(θ);
              -sin(beta)   cos(beta)*sin(θ) -cos(beta)*cos(θ)])
  mr = btor*mbe
  return mr
end

# 矢量夹角函数
function aoaang(vec_,x_=[1,0]) # 二维坐标系下，任意矢量与x轴夹角(-180deg~180deg)
#  x_ = [1,0]
  if vec_[2]>=x_[2] # 判断矢量与x轴夹角方向
    return -acos(dot(vec_,x_)/(norm(vec_)*norm(x_)))
  else
    return acos(dot(vec_,x_)/(norm(vec_)*norm(x_)))
  end
end

function cosvecang(a,b) #计算夹角余弦值
  if norm(a)<=1e-3||norm(b)<=1e-3
    return 0.0
  elseif norm(a+b)<1e-3||norm(a-b)<1e-3
    return cos(π)
  elseif norm(a-b)<1e-3
    return cos(0.0)
  else
    return dot(a,b)/(norm(a)*norm(b))
  end
end
