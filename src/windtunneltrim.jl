# 这是一个风洞配平程序
# 配平结果要求旋翼没有侧向力
# ——————要求旋翼没有俯仰和滚装力矩
# ——————要求旋翼的桨尖平面与旋翼轴垂直
# 配平结果不要求配平旋翼拉力
# 最后根据一般性的风洞培平要求，决定只通过挥舞量配平操纵量

# 程序首次编辑时间为11/11/2017

function wttrim()
  jacobi = [A0/A_a1 0 A_B1/A_a1;
              γ_B*F0 -B_A1/B_b1  γ_B*F_B1]
  soljaco = inv(jacobi)
  # delthe_lon = 0.1  #纵向周期变距步进
  # delthe_lat = 0.1  #横向周期变距步进

  while abs(Fhub[2])>1e-3 | abs(Mx)>1e-3 | abs(My)>1e-3
    the_lon = the_lon+delthe_lon
    the_lat = the_lat+delthe_lat
    aoaget(v_all,ψ,the_lon,the_lat)
  end

  return the_lon,the_lat

end 
