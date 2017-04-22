using PyPlot
# plt[:style][:use]("ggplot")
plt[:rcParams]["font.family"][1] = "serif"
# include("../src/DSR.jl")
using DSR

function test()

    T, zta = 1.0, 0.02
    xi, alpha = 0.008, 0.2
    kpa = 1.0
    phi = 0.005
    omg = 2.0*pi/T
    dt = 0.005
    n = 4001

    # t = linspace(0,(n-1)*dt,n)
    # p = (t) -> sin(2.0*π/T*t) # 正弦激励

    s0 = DSR.SDF(omg, zta)
    # s  = DSR.SDF_MYD_BW(omg, zta, kpa, phi, 0.02)
    # s.n = 10
    # s  = DSR.SDF_MYD_EP(omg, zta, kpa, phi, 0.05)
    # s  = DSR.SDF_FD(omg, zta, phi)
    # s  = DSR.SDF_VD(omg, zta, xi, alpha)

    mu, xi, kpa = 0.050, 0.02, 0.085
    s  = DSR.SDF_PID_VD(omg, zta, mu, xi, kpa)
    # s  = DSR.SDF_MAXWELL_VD(omg, zta, xi, kpa, alpha)
    # s  = DSR.SDF_MAXWELL_VD(omg, zta, xi*1000.0^(1.0-alpha), kpa, alpha)

    # r = DSR.response(s, p, t)

    a_g = DSR.readtxt("EQ-S-1.txt")
    dt_g = 0.005
    ndiv = 1
    dt = dt_g/ndiv
    nstep = convert(Int64,round(length(a_g)*ndiv*1.1))

    r0 = DSR.response(s0, a_g, dt_g, dt, nstep)
    r  = DSR.response(s,  a_g, dt_g, dt, nstep)

    figure(1,(9,7.5))
    subplot(3,1,1)
    plot(r0[:,1],r0[:,2])
    plot(r[:,1], r[:,2])
    grid(true)

    # subplot(3,1,2)
    # plot(t, r[:,2])
    # plot(t,r0[:,2])
    # grid(true)
    # subplot(3,1,3)
    # plot(t, r[:,3])
    # plot(t,r0[:,3])
    # grid(true)

    subplot2grid((3,1),(1,0),rowspan=2,colspan=1)
    plot(r[:,2],r[:,5])
    # plot(r[:,2],r[:,6])
    plot(r[:,7],r[:,8])
    grid(true)

    # show()
end

test();