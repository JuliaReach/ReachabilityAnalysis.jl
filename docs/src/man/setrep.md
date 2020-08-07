# Set representations

## Polyhedra





## Reach-sets

## Flowpipes

## Taylor models


f(x) = -6x^3 + (13/3)x^2 + (31/3)x
dom = -3.5 .. 3.5

plot(f, -3.5, 3.5, lab="f", xlab="x")

x = Taylor1(5)
set_taylor1_varname("x")
f(x)


rem = 0 .. 0
x0 = 0.0
dom = -3.5 .. 3.5
tm = TaylorModel1(f(x), rem, x0, dom)
