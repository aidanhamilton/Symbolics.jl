using Symbolics
using Test

@variables x,y

A = [0 1 
     -u[2]^2/u[1]^2 + u[1] 2*u[2]/u[1] ]

eigevals(A) 