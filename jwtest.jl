@time let
mat = rand(3000,3000)
mat += mat'
data = eigen(mat)
end

return nothing
