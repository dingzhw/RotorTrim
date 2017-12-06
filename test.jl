include("const.jl")
include("flapreponse.jl")

a = 0.0
b = 0.0
c = 1
i = 1
while true
  dftmp = dflapre(a,b,c*100)
  if i%2==0
    c = -1
  end
  print("$(a*180/π)\t")
  print("$(b*10)\n")
  if abs(a)>π/2
    print("Error")
    break
  end
  a = dftmp[1]
  b = dftmp[2]
  i = i+1
end
