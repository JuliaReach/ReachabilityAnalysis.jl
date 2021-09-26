## Domain splitting

A common technique to reduce wrapping effects is to split the set of initial
states. If an initial-value problem has been setup with an *array of sets*, then
the flowpipe starting from each initial set scomputed in parallel, using Julia's
built-in multithreaded support.


!!! note
    To turn off multithreading, pass the `multithreaded=false` option flag to
    `solve` method. It is `true` by default.


!!! note
    To change the number of threads being used, change the `THREADS` flag in . . .
