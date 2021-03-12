# Frequently Asked Questions - Wrappers

Here I try to gather some of the questions that you might have while using
`emcee`'s wrapper. If you have an issue that is not reflected here, please
submit it to the relevant place (links available at the FAQ for Katu)


1. The wrapper crashes in a random way, what can I do?

First, make sure that you are using `emcee` version '2.2.1'. There is *no*
support for the latest version of `emcee` as that changes the way of accessing
its internal members and the wrapper relies on doing it.

Having discarded that, that there are dependencies that you have forgotten 
and Python errors, please check the relevant question for wrappers in general.

If the problem still persist, please open a new issue detailing what you have
done and attaching all the relevant code.


2. The wrapper takes a really long time to even start, what can I do?

Unfortunately, when `emcee` checks for updating the position of a walker,
it does a comparison between the values of the log likelihood at both points.
In the case where one of them is `-infinity` (AKA invalid point), the
comparison evaluates to false and the new point is **never** accepted. In order
to avoid this problem, before starting `emcee` the wrapper checks the
position of every walker and makes sure that it starts from a valid set of
parameters.

In the cases where this search does not finish, the user should modify the
starting point and the spread of the initial walkers. See question 4.


3. The wrapper is taking a really long time to finish, what can I do?

The first point to solving this problem is for the user to check and understand
how `emcee` works and to understand the stopping condition that we have
added to the wrapper.

By default, `emcee`'s wrapper is set to stop when the relative change in the
value of the medians and the 16 and 84 percentiles betwen the last 50 and the
next 50 blocks of data is less than 1e-3.

In theory, this is a very stringeng requirement and in practice it is almost
never achieved. The variation depends, generally, on the shape of the log
likelihood and the number of walkers, being harder to achieve the lower the
number of walkers.

This means that, in general, the walker never finishes working and that it
should be stopped by the user after an analysis of the output. Please, refer
to question 5 for advices on how to stop the wrapper.


4. Where should I start the walkers from, how much spread?

This is a _very_ complicated question that has no easy answer as it depends
heavily on the shape of the log likelihood function: number of maxima,
behaviour near the maxima and the edges of the parameter space...

In general, my advice is the following: A 'pre-analysis' with `Katu` should
be done in order to find a good starting position (read a high log likelihood)
and make `emcee` start from there. For the cases where this is not possible
or the user can not find a good starting position, then the best is usually
to let the wrapper do all the work by starting at some point that we
know is valid and with a wide spread and wait for a few steps. Usually, the
log likelihood will improve quickly and you can switch so that the new run
starts at the new best point. After some cycles you should be good to go for
the definite run.


5. When should I stop the wrapper?

As explained in question 3, the wrapper has a very stringent stopping criteria
that in practice it is almost never achieved. This leaves us with manually
stopping the wrapper as the way of, well... stopping it.

There is no exact number of steps that `emcee` needs to run before we are
sure that we have good statistical data, so the stopping has to be eye balled.
During the runs, the wrapper will output to the terminal the best set of
parameters along with the median and the first sigma deviation of *the whole*
run. In the limit of infinite steps that should converge towards fixed numbers.
When the variation is very small, one can guess that the run could be
finishing.

Either way, the best way of checking the status of the run is to analyze the
whole data. Each step, the wrapper creates a numpy file with the whole chain
of walkers and their positions that can be analyzed.

We know, from the general behaviour of MCMC that there is going to be a first
set of steps that correspond to a burn-in and then the 'real' data starts.
The best way of detecting this is to use the `plot_bin_evolution` function
from `utils/analyzer_utils.py` with something like 16 or so bins (more if you
have more data). What you should get is that for every dimension there is a
first drift and then some stabilization. Note the number of the latest bin
where the drift ends for all the dimensions and calculate the number of 'good'
points. If it is bigger than a few thousands then you should be good to go
and to get good statistical information.

Of course, this is just eye-balling things and your problem might need other
kind of treatment for the output. In that case, you can modify the given
code to suit your needs or ask for help by opening an issue.
