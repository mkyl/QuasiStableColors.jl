module datasets

using DelimitedFiles
using CSV
using Graphs
using SparseArrays

DBLP = "Datasets/DBLP.tsv"
AIRPORTS = "Datasets/routes.csv"
ENRON = "Datasets/enron.tsv"
ASTROPHYSICS = "Datasets/astrophysics.tsv"
FACEBOOK = "Datasets/facebook.tsv"
TWITTER = "Datasets/twitter.tsv"
EPINIONS = "Datasets/epinions.tsv"
DEEZER = "Datasets/deezer.tsv"
TSUKUBA0 = "Datasets/BVZ-tsukuba0.max"
TSUKUBA2 = "Datasets/BVZ-tsukuba2.max"
SIMCELLS = "Datasets/simcells.sd2.max"
CELLS = "Datasets/cells.sd2.max"
VENUS0 = "Datasets/BVZ-venus0.max"
VENUS1 = "Datasets/BVZ-venus1.max"
SAWTOOTH0 = "Datasets/BVZ-sawtooth0.max"
SAWTOOTH1 = "Datasets/BVZ-sawtooth1.max"

function tsv_graph(name)
    M = readdlm(name, '\t', Int, '\n', comments=true, comment_char='#')
    n = max(maximum(M[:, 1]), maximum(M[:, 2]))
    G = SimpleGraph{Int}(n)
    for r in eachrow(M)
        u, v = r[1], r[2]
        add_edge!(G, u, v)
    end
    return G
end

enron = () -> tsv_graph(ENRON)
astrophysics = () -> tsv_graph(ASTROPHYSICS)
facebook = () -> tsv_graph(FACEBOOK)
twitter = () -> tsv_graph(TWITTER)
epinions = () -> tsv_graph(EPINIONS)
deezer = () -> tsv_graph(DEEZER)

function dblp()
    M = readdlm(DBLP, '\t', Int, '\n')
    n = max(maximum(M[:, 1]), maximum(M[:, 2]))
    G = SimpleGraph{Int}(n)
    for r in eachrow(M)
        u, v = r[1], r[2]
        add_edge!(G, u, v)
    end
    return G
end

function airports_labelled()
    f = CSV.File(AIRPORTS, header=false)
    V = Set(f.Column3) ∪ Set(f.Column5)
    lookup = Dict((x, i) for (i, x) in enumerate(V))

    G = SimpleDiGraph(length(V))
    for r in f
        u, v = r.Column3, r.Column5
        add_edge!(G, lookup[u], lookup[v])
    end
    return G, lookup
end

function openflight()
    f = CSV.File(AIRPORTS, header=false)
    V = Set(f.Column3) ∪ Set(f.Column5)
    lookup = Dict((x, i) for (i, x) in enumerate(V))

    G = SimpleGraph(length(V))
    for r in f
        u, v = r.Column3, r.Column5
        add_edge!(G, lookup[u], lookup[v])
    end
    return G
end

function DIMACS(file)
    # see http://lpsolve.sourceforge.net/5.5/DIMACS_maxf.htm
    local G::SimpleDiGraph{Int}
    local s::Int
    local t::Int
    I, J, V = [], [], []
    for line in eachline(file)
        if length(line) == 0
            continue
        end
        if line[1] == 'p'
            _, _, n, _ = split(line)
            G = SimpleDiGraph(parse(Int, n))
        elseif line[1] == 'n'
            _, node, type = split(line)
            if type == "s"
                s = parse(Int, node)
            elseif type == "t"
                t = parse(Int, node)
            end
        elseif line[1] == 'a'
            _, u, v, c = split(line)
            if parse(Int, c) == 0
                continue
            end

            push!(I, parse(Int, u))
            push!(J, parse(Int, v))
            push!(V, parse(Int, c))
            add_edge!(G, parse(Int, u), parse(Int, v))
        end
    end
    C::SparseMatrixCSC{Int,Int} = sparse(I, J, V, nv(G), nv(G))
    return G, C, s, t
end

tsukuba0 = () -> DIMACS(TSUKUBA0)

tsukuba2 = () -> DIMACS(TSUKUBA2)

simcells = () -> DIMACS(SIMCELLS)

cells = () -> DIMACS(CELLS)

venus0 = () -> DIMACS(VENUS0)

venus1 = () -> DIMACS(VENUS1)

sawtooth0 = () -> DIMACS(SAWTOOTH0)

sawtooth1 = () -> DIMACS(SAWTOOTH1)

end


