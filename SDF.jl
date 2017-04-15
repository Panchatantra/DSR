type SDF <: SDFDynamicSystem
    omg::Float64
    zta::Float64
    u0::Float64
    v0::Float64

    SDF(omg::Float64, zta::Float64) = new(omg, zta, 0.0, 0.0)
end


function SDF(m::Float64, c::Float64, k::Float64)
    omg = sqrt(k/m)
    zta = c/(2.0*m*omg)
    SDF(omg, zta)
end


function systemEquation(s::SDF, p::Function)
    function dydt(y::AbstractArray{Float64,1}, t::Float64)
        [y[2]; p(t)-2.0*s.zta*s.omg*y[2]-s.omg*s.omg*y[1]]
    end
    dydt
end


function response(s::SDF, p::Function, t::AbstractArray{Float64,1})
    y0 = [s.u0; s.v0]

    f = systemEquation(s, p)

    y = ode4(f, y0, t)

    u = y[:,1]
    v = y[:,2]
    fdz = 2.0*s.zta*s.omg*v
    a = map(p,t)-fdz-s.omg*s.omg*u

    [u v a fdz zeros(length(t))]
end