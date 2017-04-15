# DSR
A Dynamic Structural Response solver based on Julia

基于Julia语言的结构动力响应求解器

## Simple Usage

```julia
using DSR
T, zta = 1.0, 0.05
omg = 2.0*pi/T
t1, dt = 40, 0.005
s = SDF(omg, zta)
p = (t) -> sin(2.0*π/T*t)
r = response(s,p,t1,dt)
```