# Hybrid systems

## Introduction


Our running example is the *bouncing ball* model; although it is a very hybrid automaton, it can be used to introduce the main notions involved in hybrid systems reachability.


## Clocked linear dynamics

So far we have focused on transitions that involve "spatial" variables.
If the system under consideration has transitions governed by time variables,
i.e. by variables whose dynamics are of the form ``t' = 1``, then decoupling the
spatial variables with the clock variables gives a computational advantage.
We refer to [[HG19]].
