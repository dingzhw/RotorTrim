# The script s written for test moudle using

module soltest

importall Base

export functest,roro 

# ---
function functest(x::Float64,y::Float64)
 return x^y
end

type roro
 x::Int
 sinxx
end

function sinxx(x)
 return sin(x)
end

end
