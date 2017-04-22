using PyPlot
# plt[:style][:use]("ggplot")
plt[:rcParams]["font.family"][1] = "serif"
# include("../src/DSR.jl")
using DSR

function test()

    T, zta = 1.0, 0.05
    omg = 2.0*pi/T
    t1, dt = 20.0, 0.005
    s = SDF(omg, zta)
    p = (t) -> sin(2.0*π/T*t)
    r = response(s,p,t1,dt)

    figure("Resonant Response",(9,3))
    plot(r[:,1], r[:,2])
    xlabel("Time [s]")
    ylabel("Displacement [m]")
    grid(true)
    # savefig("resonant_response.png"; dpi=72, bbox_inches="tight")
    show()
end

test();