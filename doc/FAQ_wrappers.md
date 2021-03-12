# Frequently Asked Questions - Wrappers

Here I try to gather some of the questions that you might have while using
any of the wrappers that I have added to support to. If you have an issue
that is not reflected here and is also not at the corresponding FAQ for
the wrapper, please submit it to the relevant place (links available at the
FAQ for Katu)


1. The wrapper crashes in a random way, what can I do?

First, make sure that you have followed the instructions that I have supplied
in the `Readme` about how to use the wrappers. Then, make sure that they are
not a problem from your part (as in you have placed some files in the wrong
place or have modified the wrappers and had introduced some bug). After
discarding these kind of problems and/or the cause of the problem is `Katu`,
please test `Katu` with the configuration file that has made the wrapper crash
and refer to the `Katu` FAQ. If you are still having problems, please open
a new issue with the relevant files (wrapper, data and the configuration)
and I will try to reproduce the issue.


2. The wrapper seems stuck and does nothing, what can I do?

Depending on the problem, it might be the case that the problem is that it
is taking time for the wrapper to do whatever it needs to do. This could
happen for `emcee` if you are using a very big number of walkers and for
`multinest` when generating the first set of points.

Either way, in these cases it is good to make sure that it is not `Katu` the
part that is stuck. For this open some monitoring program such as `top` or
`ps` and see how much time has spent `Katu` trying one configuration. If it
is much higher than what it normally takes it might be stuck somewhere. Please
refer to `Katu`'s FAQ to diagnose the problem.

In these cases, it could be a signal that the tested configuration is 'strange'
in some whay. It might be interesting for the user to remove this part of
the configuration space by modifying the wrappers.


3. The Wrapper is taking forever to finish, what can I do?

Please, refer to the documentation of each wrapper, as they all work
differently and have their own set of issues.


4. Can you add support for XYZ in ABC's wrapper?
   Can you add support for XYZ statistical package?

Miaube. Open an issue stating what you need and why and I will consider it.
