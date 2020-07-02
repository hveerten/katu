filename = sprintf("../data/photon_data_%04d.tsv", @ARG1)

norm = "$1*$1*$3/($4+$5+$6+$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17)"

set xrange [1e-12:2e11]
set yrange [1e-15:1e5]
set logscale xy

plot filename u 1:($1*$1*$3) w l, \
           "" u 1:(@norm* $4) t "Inj", \
           "" u 1:(@norm* $5) t "e^- sync",  \
           "" u 1:(@norm* $6) t "e^+ sync",  \
           "" u 1:(@norm* $7) t "p^+ sync",  \
           "" u 1:(@norm*$14) t "IC up", \
           "" u 1:(@norm*$15) t "IC Down", \
           "" u 1:(@norm*$16) t "Pi0", \
           "" u 1:(@norm*$17) t "PA"
