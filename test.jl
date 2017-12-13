addprocs(4)

function test()
  a = randn(1000)
  @parallel (+) for i=1:100000
    fetch(a[randi(end)])
  end
end
