function ar1_sample_mean(N, beta, alpha, s)
    sm = 0.0
    x = beta / (1 - alpha)
    for i in 1:N 
        sm += x
        x = beta + alpha * x + s * randn()
    end
    return sm / N
end

N = 10000000
beta = 1.0
alpha = 0.9
s = 1.0
tic()
result = ar1_sample_mean(N, beta, alpha, s)
println("mean = $result")
toc()

