filename = sprintf("hadron_data_%04d.tsv", @ARG1)

norm = "$2/($4+$5+$6)"

set xrange [1:1e8]
set yrange [1e-30:1e2]

plot filename u 1:2, \
     "" u 1:(@norm*$4) t "Inj", \
     "" u 1:(@norm*$5) t "MR", \
     "" u 1:(@norm*$6) t "D"
