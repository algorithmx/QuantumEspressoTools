using Plots

spl(x) = parse.(Float64,split(x," ",keepempty=false))

function plot_ph(d; mc=:red, ms=1)
    dd = hcat(d...)
    X = dd[1,:]
    plt = scatter()
    for i=2:size(dd,1)
        plt = scatter!(X,dd[i,:], markersize=ms,markercolor=mc, legend=nothing)
    end
    return plt
end


function plot_ph!(plt, d; mc=:red, ms=1)
    dd = hcat(d...)
    X = dd[1,:]
    for i=2:size(dd,1)
        plt = scatter!(X,dd[i,:], markersize=ms,markercolor=mc, legend=nothing)
    end
    return plt
end


function plot_ph_diff(d1,d2; mc=:red, ms=10, cutoff=0.5)
    dd = hcat(d1...)
    ddd = hcat(d2...)
    X = dd[1,:]
    plt = scatter()
    for i=2:size(dd,1)
        s1 = abs.(dd[i,:].-ddd[i,:])
        ss = [min(cutoff,x) for x in s1]
        plt = scatter!(X,dd[i,:], markersize=ms.*ss,markercolor=mc, legend=nothing)
    end
    return plt
end


##

raw = [readlines("/home/dabajabaza/Downloads/ReO3_$(i)bar/ReO3_MT_PBE.frq.gp") for i=1:5]



data = [spl.(r) for r in raw]


pp = plot_ph(data[1]; mc=:black, ms=1)
pp = plot_ph!(pp, data[2]; mc=:red, ms=1)
pp = plot_ph!(pp, data[3]; mc=:blue, ms=1)
pp = plot_ph!(pp, data[4]; mc=:green, ms=1)
pp = plot_ph!(pp, data[5]; mc=:orange, ms=1)

##

PS = "MT_PW_LDA"
d1 = spl.(readlines("/home/dabajabaza/Downloads/ReO3_1bar/ReO3_$(PS).frq.gp")) ;
d2 = spl.(readlines("/home/dabajabaza/Downloads/ReO3_450bar/ReO3_$(PS).frq.gp")) ;
d3 = spl.(readlines("/home/dabajabaza/Downloads/ReO3_250bar/ReO3_$(PS).frq.gp")) ;

pp = plot_ph_diff(d1, d3; mc=:black, ms=8, cutoff=2.0)
pp = plot_ph_diff(d1, d2; mc=:red,   ms=8, cutoff=2.0)

