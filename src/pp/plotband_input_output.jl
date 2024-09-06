#=
cat > plotband.in <<EOF
al.dyn.freq
0 400
freq.plot
freq.ps
0.0
50.0 0.0
EOF

cat > plotband.in <<EOF
alas.freq
0 600
freq.plot
freq.ps
0.0
50.0 0.0
EOF


cat > plotband.in <<EOF
bn.freq
0 1650
freq.disp.plot
freq.disp.ps
0.0
50.0 0.0
EOF


cat > plotband.in << EOF
fe.band
0 50
ciao
EOF

=#

function plotband_input(config::Dict)
    inp = [
        "$(config[:freq])",
        "$(config[:freq_min_max_cm_1][1])   $(config[:freq_min_max_cm_1][2])",
        "$(config[:plot])",
        "$(config[:ps])",
        "$(config[:E_Fermi])",
        "$(config[:dE])  0.0",
        "ciao",
        ""
    ]
    return  inp
end
