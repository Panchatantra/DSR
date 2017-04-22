type SDF_FD <: SDFDynamicSystem
    omg::Float64
    zta::Float64
    phi::Float64
    u0::Float64
    v0::Float64

    SDF_FD(omg::Float64, zta::Float64, phi::Float64) = new(omg, zta, phi, 0.0, 0.0)
end


function systemEquation(s::SDF_FD, p::Function)
    function dydt(y::AbstractArray{Float64,1}, t::Float64)
        [y[2]; p(t)-2.0*s.zta*s.omg*y[2]-s.phi*g*sign(y[2])-s.omg*s.omg*y[1]]
    end
    dydt
end


function response(s::SDF_FD, p::Function, t::AbstractArray{Float64,1})
    y0 = [s.u0; s.v0]

    f = systemEquation(s, p)

    y = ode4(f, y0, t)

    u = y[:,1]
    v = y[:,2]
    fdz = 2.0*s.zta*s.omg*v
    fdp = s.phi*g*sign(v)
    a = map(p,t)-fdz-fdp-s.omg*s.omg*u

    [u v a fdz fdp]
end