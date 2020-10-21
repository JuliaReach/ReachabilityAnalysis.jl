## Overview

We organize the models by the type of nonlinearities (if there are some), and
whether they are purely continuous or present discrete transitions, i.e. hybrid systems.
We have added a column with the associated scientific domain, and another column with
the number of state variables. Roughly speaking, a higher number of state variables
usually correspnds to problems which are harder to solve, altough strictly speaking,
this usually depends on the property to be verified.

Column `P.V.` refers to the cases
where the example is presented with at lest one instance with parameter variation.

## Linear continuous

|Name|Area|State dim.|
|----|------|---------|
|Damped oscillator|Physics|2|
|Building|Mechanical Engineering|48|
|International Space Station|Aerospace Engineering|270|
|Modified Nodal Analysis 1|Electronics|1002|
|Modified Nodal Analysis 2|Electronics|10913|
|Heat PDE|Physics|125, 1000, 8000, 125000|

## Linear hybrid

|Name|Area|State dim.|
|----|------|---------|
|Amplifier circuit|Electronic Engineering|2|
|Electromechanic break|Electronic Engineering|
|Gearbox|Mechanical Engineering||
|Platoon|Autonomous Driving||
|Powertrain|Mechanical Engineering||
|Transmission line|Power Systems Stability||

## Nonlinear continuous

|Name|Area|State dim.|
|----|------|---------|
|Brusselator||||
|Laub-Loomis|Molecular Biology|7|
|Lorenz system|||
|Lotka-Volterra|||
|Production-Destruction|Electrical engineering||
|Quadrotor|||
|SEIR Model|||
|Van der Pol|||

## Nonlinear hybrid

|Name|Area|State dim.|
|----|------|---------|
|Lotka-Volterra w/crossing|Biology, Nonlinear physics||
|Spacecraft|||
