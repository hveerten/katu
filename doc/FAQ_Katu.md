# Frequently Asked Questions - Katu

Here I try to gather some of the questions that you can have while using Katu.
Hopefully it will be as detailed as possible.

Remember that the issues can/should be left at both:
 - github https://github.com/hveerten/katu/issues
 - gitlab https://gitlab.com/Noxbru/Katu/-/issues


1. Where did Katu came from?

Katu, in its basic form, is the result of Bruno Jiménez's PhD. All the code has
been written by him and all the praises/blames are his.


2. Where can I learn a bit more about Katu's internals?

Besides reading the code, the best place to learn about Katu is the publication
where we present it to the public <fill with link to it at some point>. In it,
I have tried to be as thorough as possible in explaining the general 
characteristics of Katu and the physical processes it simulates.

Note that it reflects the state of Katu at the time of writting the article,
and it might be the case that Katu has evolved beyond it.


3. Why have you done/not done XYZ?

The general answer is a mixture of time-effort-reward. It takes time to
research whatever I want to add to Katu, effort to implement and debug it to
make sure that it works as it should, and this has to be balanced by the
potential reward.

As an example, there are certain processes that we do not fully simulate or that
we only do in a first order approximation because it is not worth going
further.


4. Are you going to add XYZ to Katu?

Maybe, maybe not. It will depend on whether I need the new feature and/or if
many people ask for it.

If you need something added to Katu, please open a new issue stating what
you need, what your use case is and any reference that supports your claim
or I can use to add the requested feature.


5. There are many numerical parameters in the configuration file, any
    guidelines for setting them?

`Katu` has a lot of parameters and most of them can be set either by physical
reasonings or are set by the wrappers that the user must have written. Still,
there are numerical parameters that must be set.

- `dt`: corresponds to the initial value for the time step and should be set to
    a small value compared to `dt_max`. (E.g. `dt` ~ `dt_max` / 100)
- `dt_max`: corresponds to the maximum value for the time step. By default,
    `Katu` takes the minium value between the user supplied one and the
    light crossing time. This is usually enough to avoid oscillations while
    being fast. In some edge cases, this might not be enough and a lower value
    can be set.

- `tol`: Should be set to a low value, but the user **MUST** make sure that it is
    low enough by analysing the evolution of its value with time.
    In our experiments, a value of 1e-6 is low enough.


- `size`: Corresponds to the size of either of the given distributions. The
    size must be set as a compromise between the desired accuracy and the
    allowed running time and memory requirements.
    As a rough approximation, the heaviest algorithms/tables scale as the
    cube of size (square photon/electron times electron/photon, respectively)
    and must be taken into account.
    By some experimentation, we have found that: n\_photon ~ 128 and
    n\_electron ~ n\_proton ~ 160 gives good results while being fast enough.


6. Is it possible to start `Katu` from given populations?

Actually, it is! `Katu` does support to save and load populations with the
`state_{save,load}_state_from_file` functions. Although these functions might
break at random times as I do not use them much and sometimes I forget to
update them when I modify any of the code related to populations. Note
that there are no examples for this and that you will have to write your
own code.



7. Katu has crashed/shown some error, what should I do?

First, make sure that all the examples (except for that of a sphere inside
an sphere) crash in the same (or similar way). Then, open an issue copying
the output of Katu and attach the configuration file that created the problem.
I will try to reproduce the issue locally and help you with it.


8. Katu seems to be stuck somewhere/does not finish, what should I do?

Even though we have tried to make Katu resistant to some problems that make
the kinetic equation to not converge numerically, there are still weird cases
where this is possible. Please, test your configuration with the
`example_simple` example and note the evolution of the different statistics
it outputs.

In the case where these statistics reach a low value (1e-8 or so) and do not go
lower from it, it might be the case that numerically it is the lowest
achievable value for the tolerance. Unfortunatelly that means that with the
current algorithms for calculating the tolerance we can not go further. If you
*trully* need to achieve a lower tolerance, open a new issue with the
configuration and stating *very* clearly why you need it.

In the case where the statistics oscillate. These are weird cases that even
if we have tried hard to avoid can still happen. The best way to solve them
are to reduce the value of the maximum allowed value of dt in the configuration
file (`dt_max`). Note that by default, `Katu` clamps it to the light crossing
timescale, so it must be set to a value lower than it.
