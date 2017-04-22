type SDF_MYD_BW <: SDFDynamicSystem
    omg::Float64
    zta::Float64
    kpa::Float64
    phi::Float64
    alpha::Float64
    u0::Float64
    v0::Float64
    z0::Float64
    n::Int64
    bta::Float64

    SDF_MYD_BW(omg::Float64, zta::Float64, kpa::Float64, phi::Float64, alpha::Float64, n::Int64, bta::Float64) = new(omg, zta, kpa, phi, alpha, 0.0, 0.0, 0.0, n, bta)
    SDF_MYD_BW(omg::Float64, zta::Float64, kpa::Float64, phi::Float64, alpha::Float64) = new(omg, zta, kpa, phi, alpha, 0.0, 0.0, 0.0, 5, 0.5)
    SDF_MYD_BW(omg::Float64, zta::Float64, kpa::Float64, phi::Float64) = new(omg, zta, kpa, phi, 0.0, 0.0, 0.0, 0.0, 5, 0.5)
end


function systemEquation(s::SDF_MYD_BW, p::Function)
    function dydt(y::AbstractArray{Float64,1}, t::Float64)
        u = y[1]
        v = y[2]
        z = y[3]
        omg2 = s.omg*s.omg
        dvdt = p(t)-2.0*s.zta*s.omg*v-(s.alpha*s.kpa*omg2*u+(1.0-s.alpha)*s.phi*g*z)-omg2*u
        dzdt = v - (s.bta*abs(v)*abs(z)^(s.n-1)*z + (1.0-s.bta)*v*abs(z)^(s.n))
        dzdt *= s.kpa*omg2/(s.phi*g)
        [v; dvdt; dzdt]
    end
    dydt
end


function response(s::SDF_MYD_BW, p::Function, t::AbstractArray{Float64,1})
    y0 = [s.u0; s.v0; s.z0]

    f = systemEquation(s, p)

    y = ode4(f, y0, t)

    u = y[:,1]
    v = y[:,2]
    z = y[:,3]
    omg2 = s.omg*s.omg
    fdz = 2.0*s.zta*s.omg*v
    fdp = s.alpha*s.kpa*omg2*u+(1.0-s.alpha)*s.phi*g*z
    a = map(p,t)-fdz-fdp-s.omg*s.omg*u

    [u v a fdz fdp]
end