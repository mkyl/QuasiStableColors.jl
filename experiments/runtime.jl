using Plots

# collected from logs
l = ["qap15, |P|=10", "qap15, |P|=100", "sawtooth0, |P|=15", "sawtooth0, |P|=35", "cells, |P|=15",
    "cells, |P|=35", "astrophysics, |P|=5", "astrophysics, |P|=200"]
x = [1.14035, 9.822, 0.458346, 1.4518, 38.34, 115.57, 0.065539, 26.79]
y = [0.000693, 0.002232, 0.000030, 0.0002, 0.000038, 0.000134, 0.85, 12.721]

bar([((x ./ (x + y)) + y ./ (x + y)) (x ./ (x + y))], labels=["Solving" "Coloring"],
    xticks=(1:8, l), xrotation=15, size=(400, 250), legend=:bottomleft)
savefig("breakdown.pdf")
