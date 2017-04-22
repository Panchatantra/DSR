# DSR
A Dynamic Structural Response solver based on Julia

基于Julia语言的结构动力响应求解器

## Simple Usage

For a classic Single-Degree-of-Freedom system under resonant harmonic load, the equation of motion is

![eq1](http://www.sciweavers.org/tex2img.php?eq=%5Cddot%7Bu%7D%20%2B%202%20%5Czeta%20%5Comega%20%5Cdot%7Bu%7D%20%2B%20%5Comega%5E2%20u%20%3D%20%5Csin%20%28%5Comega%20t%29%20&bc=White&fc=Black&im=png&fs=12&ff=modern&edit=0)

By using DSR, the response can be solved as follows. 

```julia
using DSR

T, zta = 1.0, 0.05
omg = 2.0*pi/T
t1, dt = 40, 0.005
s = SDF(omg, zta)
p = (t) -> sin(2.0*π/T*t)
r = response(s,p,t1,dt)
```

With the help of PyPlot, the resonant response can be drawn:

```julia
using PyPlot

figure("Resonant Response",(12,4))
plot(r[:,1], r[:,2])
grid(true)
xlabel("Time [s]")
ylabel("Displacement [m]")
show()
```

![resonant_response.png](resonant_response.png "Resonant Response")