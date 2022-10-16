using Plots
using Graphs
using Random

include("../main.jl")

Random.seed!(0)

k = 100
n = 1000

G0 = erdos_renyi(k, 0.05)
P = [[] for _ in 1:k]

for i in 1:n
    push!(P[i%k+1], i)
end

G = SimpleGraph(n)

for Pi in 1:k
    N = neighbors(G0, Pi)
    for Pj in N
        for u in P[Pi]
            for v in P[Pj]
                add_edge!(G, u, v)
            end
        end
    end
end

stable = []
q_stable = []

println(G)

push!(stable, length(refine_fixpoint(G)))
push!(q_stable, length(refine_fixpoint(G, eps=2.0)))
for i in [1 2 4 8 16 32 64 128 256]
    for _ in 1:i
        add_edge!(G, rand(1:n), rand(1:n))
    end
    push!(stable, length(refine_fixpoint(G, early_stop=1000)))
    push!(q_stable, length(refine_fixpoint(G, eps=4.0, early_stop=500)))
end

plot([stable, q_stable], labels=["stable" "q-stable, q=4"], legend=:topleft,
    markershape=:auto, ylim=(0, 1100), size=(425, 200))
xticks!([1:64;], ["Original", "1", "2", "4", "8", "16", "32", "64", "128", "256"])
xlabel!("Random edges added")
ylabel!("Colors needed")
hline!([1000], label="# of graph vertices", linestyle=:dash)

savefig("robustness.pdf")
