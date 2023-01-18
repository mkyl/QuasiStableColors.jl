module datasets

using DelimitedFiles
using CSV
using Graphs
using DataDeps

register(DataDep("OpenFlights",
    "OpenFlights.org Route database, see https://openflights.org/data.html",
    "https://github.com/jpatokal/openflights/raw/1b649733a4d35f7668d901b936361ad4d3b7c6ac/data/routes.dat",
    "bd373706238134f619c624c606dccc74c05c2582a977c489c81de501735f2390",
))

register(DataDep("DBLP",
    "DBLP computer science co-authorship graph, see https://snap.stanford.edu/data/com-DBLP.html",
    "https://snap.stanford.edu/data/bigdata/communities/com-dblp.ungraph.txt.gz",
    "9eb0bd30312ddd04e2624f7c36c0983a2e99b116f0385be5a7fce6d6170f4cb3",
    post_fetch_method=unpack,
))

register(DataDep("Karate",
    "Zachary's Karate Club social network, see https://en.wikipedia.org/wiki/Zachary%27s_karate_club",
    "http://konect.cc/files/download.tsv.ucidata-zachary.tar.bz2",
    "37ad42a41e95a09c1d1e5de5ce9eae3a53d8170c6c4af6c4a0121885cafc2cef",
    post_fetch_method=unpack,
))

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

function tsv_graph(path)
    M = readdlm(path, '\t', Int, '\n', comments=true, comment_char='#')
    n = max(maximum(M[:, 1]), maximum(M[:, 2]))
    G = SimpleGraph{Int}(n)
    for r in eachrow(M)
        u, v = r[1], r[2]
        add_edge!(G, u, v)
    end
    return G
end

karate = () -> tsv_graph(datadep"Karate/ucidata-zachary/out.ucidata-zachary")
enron = () -> tsv_graph(ENRON)
astrophysics = () -> tsv_graph(ASTROPHYSICS)
facebook = () -> tsv_graph(FACEBOOK)
twitter = () -> tsv_graph(TWITTER)
epinions = () -> tsv_graph(EPINIONS)
deezer = () -> tsv_graph(DEEZER)

function dblp()
    path = datadep"DBLP/com-dblp.ungraph.txt"
    M = readdlm(path, '\t', Int, '\n', comments=true, comment_char='#')
    n = max(maximum(M[:, 1]), maximum(M[:, 2]))
    G = SimpleGraph{Int}(n)
    for r in eachrow(M)
        u, v = r[1], r[2]
        add_edge!(G, u, v)
    end
    return G
end

function airports_labelled()
    path = datadep"OpenFlights/routes.dat"
    f = CSV.File(path, header=false)
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
    path = datadep"OpenFlights/routes.dat"
    f = CSV.File(path, header=false)
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


