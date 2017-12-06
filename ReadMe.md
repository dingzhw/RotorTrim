> The most important thing from C++ to Julia is to change
> the index of the "ARRAYS"
> Must be done before the first run


# Change Log
 1. I have transfer the C++ file into a Julia version
 2. It seems nice and I would like to make it OPEN SOURCE
 3. 暂时看来，整体模块化的工作初见成效，各个模块**可维护性**
      和**可修改性**还是不错的，后续将在保证整体框架合理准确的情况下继续添加模块和修正已有模块[Time Now is 11:14 PM 11/11/2017]

# Details
  1. ~~Now I am handling with the v_all~~
  2. ~~Vortex Core Model is to be taken into consideration~~
  3. ~~Assume that if a vortex is in the line of a vortex line, then the
     induced velocity of the vortex by the vortex line is zero.~~
  4. ~~Wake Points Vector Array Fomat is WP::Array{Vector}(Nwp,Nb,NR)~~
  5. Today is 10/26/2017
    1. The blade element partation is still to be adjusted.
    2. The output of force is too many.
    3. The freewake must be cut.
  6. Today is 11/11/2017
    1. 后续任务首先是需要在rotorforce的文档中将桨毂俯仰力矩、滚装力矩以及
        旋翼的扭矩计算方式列出来，这项任务应当最晚在**下周一**完成
    2. 今天尝试了一些求解非线性系统的方法，最后发现还是只有雅克比矩阵是可以拿来用的
    3. 所有最后还是决定采用风洞配平，并且，只考虑配平挥舞，也就是通过改变操纵量，
        最后实现一阶纵向挥舞和横向挥舞都为**零**
    3. 下一步就是将动态入流的模型和配平的方法结合起来
