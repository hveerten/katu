filename = sprintf("../data/lepton_data_%04d.tsv", @ARG1)

norm = "$2/($3+$4+$5+$6+$7)"

set xrange [1:2e11]
set yrange [1e-30:1e2]

plot filename u 1:2, \
     "" u 1:(@norm*$3) t "Inj", \
     "" u 1:(@norm*$5) t "PP",  \
     "" u 1:(@norm*$6) t "BH",  \
     "" u 1:(@norm*$7) t "MD"
