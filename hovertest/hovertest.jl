# this script is used for testing hover functions

include(pwd()*"/src/const.jl")
include(pwd()*"/src/mathfunctions.jl")
include(pwd()*"/src/aoaget.jl")
include(pwd()*"/src/clcdget.jl")
include(pwd()*"/src/rotorforce.jl")
include(pwd()*"/hovertest/hoverfunc.jl")

hovertmp = uitest2hover(2.*π/180,0.,2.*π/180,rb)
print("the result is : $(hovertmp)")
