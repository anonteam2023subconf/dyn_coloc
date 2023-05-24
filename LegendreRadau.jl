include("LegendreRadau_Coefs.jl")
function compute_diffMatrix(tau, skip_first=false)
    n = size(tau)[1]
    prod = zeros(n,n)
    eval0 = zeros(n) .+ 1
    i_min= 1
    if skip_first
        i_min=2
    end
    for i=i_min:n
        prod_ = 1. 
        for k=i_min:n
            if k != i
                prod_ *= (tau[i]-tau[k])
            end
        end

        for k=i_min:n
            if k != i
                eval0[i] /= (tau[i]-tau[k])
                eval0[i] *= (0. -tau[k])
            end
        end
        
        for j=1:n
            for l=i_min:n
                if l!=i
                    partial_prod = 1.
                    for k=i_min:n
                        if k != i
                            if l != k
                                partial_prod *= (tau[j]-tau[k])
                            end
                        end
                    end
                    prod[i,j] += partial_prod / prod_
                end
            end
        end
    end
    return prod,eval0
end


tau_degree =Dict()
w_degree =Dict()

for a=1:size(tau_w)[1]
    tau_degree[string(size(tau_w[a][1])[1]-1)] = reverse(1 .- tau_w[a][1][2:end])/2.
    w_degree[string(size(tau_w[a][1])[1]-1)] = reverse((tau_w[a])[2][2:end] ./sum(tau_w[a][2]))
end

diff_matrix_degree = Dict()
diff_matrix_degree_skip_ = Dict()
diff_matrix_inv_ = Dict()
for (key, value) in tau_degree
    diff_matrix, eval0 = compute_diffMatrix(value)
    diff_matrix_degree[key]= diff_matrix
    diff_matrix_degree_skip_[key] = compute_diffMatrix(value, true)[1]
    A_ = transpose(diff_matrix) .+ 0.
    A_[1, :] .= 0.
    A_[1, 1] = 1
    diff_matrix_inv_[key] = inv(A_)
end

function get_diff_matrix()
    return diff_matrix_degree, diff_matrix_inv_, tau_degree, w_degree
end

function compute_lagrange_interpo(t_a::Float64, t_b::Float64, t_is::Vector{Float64}, degree::Int64, transpose::Bool=false)

    ans = zeros(degree, size(t_is)[1])
    rescaled = (t_is .- t_a) ./ (t_b - t_a)
    if transpose
        ans = zeros(size(t_is)[1], degree)
    end
    tau= tau_degree[string(degree)]
    for i=1:degree
      prod_ = 1. 
      prods_i = ones(size(t_is)[1])
      for k=1:degree
          if k != i
              prod_ *= (tau[i]-tau[k])
              prods_i = prods_i .* (rescaled .- tau[k])
          end
      end
      if transpose
        ans[:,i] += prods_i ./ prod_
      else
        ans[i,:] += prods_i ./ prod_
      end
    end
    return ans
  end