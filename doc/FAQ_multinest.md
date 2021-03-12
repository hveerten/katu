# Frequently Asked Questions - Wrappers

Here I try to gather some of the questions that you might have while using
`multinest`'s wrapper. If you have an issue that is not reflected here, please
submit it to the relevant place (links available at the FAQ for Katu)


1. The wrapper crashes in a random way, what can I do?

Remember that this wrapper needs to have `multinest` installed and available.
In particular, remember that the location of the library *must* have been added
to `$LD_LIBRARY_PATH` for the cases where `multinest` has been compilated
from source.

Having discarded that, that there are dependencies that you have forgotten 
and Python errors, please check the relevant question for wrappers in general.

If the problem still persist, please open a new issue detailing what you have
done and attaching all the relevant code.


2. The wrapper takes a really long time to even start, what can I do?

The first thing that `multinest` does is to find a number of points to start
with (by default 400) and after that it will output something every ~10
replacements. In order to make sure that there is no problem with `Katu`,
follow the instructions on detecting if `Katu` is stuck. If not, then you
just need to give everything more time, it should be working fine.


3. The only output I get is warnings about the log likelihood not being finite,
    what can I do?

This warning is harmless (even if it floods the output) and it is telling you
that `multinest` has tried a point that is forbidden by the log likely function.

That said, it could point out to some disparity between your `prior` function
(the function that translates from hypercube coordinates to physical
coordinates) and the `lnprior` function (the function that evaluates the prior
probability of the point). In those cases, the user should fix and make those
two functions match.

In other cases, it could be the result of using hard cutt-offs for upper and
lower limits. In those cases there is nothing that we can do but to ignore
the output.


4. The wrapper is taking a really long time to finish, what can I do?

`multinest` is set up so that it will finish eventually, and the condition
is that the variation of the log likelihood caused by the best point is less
than 0.5.

If you want to check how `multinest` is doing, besides checking its output one
can also analyze the files it writes. The length of <prefix>.txt is a good
indicator of it, the longer the better in general. When the run is finishing,
the log likelihood of the last live point should tend towards that of the
best point so you can run:

`awk -v CONVFMT='%.4g' '{ for (i=1; i<=NF; ++i) $i+=0; print }' data_phys_live.points | sort -g -k 9 | column -t | (head -n 5; tail -n 5)`

(modify it to your file name and the number of dimensions).

This will output the worst and best five points. The value of the worst
log likelihood (top - second from the right) should be not much below the best
(bottom - second from the right)

The other option is to get the whole set of points (from <prefix>\ev.dat and
<prefix>\_phys\_live.points) and make the whole analysis oneself.
