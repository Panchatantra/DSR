using PyPlot
# plt[:style][:use]("classic")
# plt[:rcParams]["font.family"][1] = "serif"
# include("../src/DSR.jl")
using DSR

function test1()
    # A classic Single-Degree-of-Freedom system under resonant harmonic load
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
    # savefig("../Docs/resonant_response.png"; dpi=96, bbox_inches="tight")
    show()
end


function test2()
    # A classic Single-Degree-of-Freedom system under resonant harmonic load
    T, zta = 1.0, 0.05
    omg = 2.0*pi/T
    xi, alpha = 0.015, 0.2
    s = SDF_VD(omg, zta, xi, alpha)

    a_g = readtxt("EQ-S-1.txt")
    dt = 0.005
    nstep = 9000
    r  = response(s,  a_g, dt, dt, nstep)

    figure("Seismic Response of SDF with VD",(9,9))
    subplot(3,1,1)
    plot(r[:,1], r[:,2])
    xlabel("\$ t [\\rm{s}] \$")
    ylabel("\$ u [\\rm{m}] \$")
    grid(true)

    subplot2grid((3,1),(1,0),rowspan=2,colspan=1)
    plot(r[:,2],r[:,5],label="Inherent Damping")
    plot(r[:,2],r[:,6],label="Nonlinear Damper")
    xlabel("\$ u [\\rm{m}] \$")
    ylabel("\$ f_d/m [ \\rm{m/s^{2}} ] \$")
    grid(true)
    legend()
    # savefig("../Docs/seismic_response_VD.png"; dpi=96, bbox_inches="tight")
    show()
end


test1();
# test2();