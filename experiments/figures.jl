using Graphs
using GraphPlot
using Colors
using Compose
import Cairo, Fontconfig
using Random: MersenneTwister, shuffle!

rng = MersenneTwister(1904)

function lookup(P)
    lookup = Dict{Int,Int}()
    sizehint!(lookup, sum(length(x) for x in P))
    for (color, nodes) in enumerate(P)
        for x in nodes
            lookup[x] = color
        end
    end
    return lookup
end

G = smallgraph(:karate)
P_q = refine_fixpoint(G, early_stop=6)
P_s = refine_fixpoint(G, early_stop=100)

C_q = distinguishable_colors(length(P_q), [RGB(1, 1, 1), RGB(0, 0, 0)], dropseed=true)
C_s = distinguishable_colors(length(P_s), lchoices=range(20, stop=75, length=20))
shuffle!(rng, C_q)
shuffle!(rng, C_s)

X_q = lookup(P_q)
X_s = lookup(P_s)

NC_q = [C_q[X_q[i]] for i in 1:nv(G)]
NC_s = [C_s[X_s[i]] for i in 1:nv(G)]

x_loc, y_loc = spring_layout(G, 2 .* rand(rng, nv(G)) .- 1.0,
    2 .* rand(rng, nv(G)) .- 1.0; MAXITER=35, INITTEMP=2.0, C=5.0)

nodelabel = collect(1:nv(G))
plot_q = gplot(G, x_loc, y_loc, nodefillc=NC_q, nodelabel=nodelabel, nodelabelc="white")
plot_s = gplot(G, x_loc, y_loc, nodefillc=NC_s, nodelabel=nodelabel, nodelabelc="white")

draw(PDF("karate_q_stable.pdf", 14cm, 7cm), plot_q)
draw(PDF("karate_stable.pdf", 14cm, 7cm), plot_s)
