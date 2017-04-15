type SDF_VD <: SDFDynamicSystem
    omg::Float64
    zta::Float64
    xi::Float64
    alpha::Float64
    u0::Float64
    v0::Float64

    SDF_VD(omg::Float64, zta::Float64, xi::Float64, alpha::Float64) = new(omg, zta, xi, alpha, 0.0, 0.0)
    SDF_VD(omg::Float64, zta::Float64, xi::Float64) = new(omg, zta, xi, 1.0, 0.0, 0.0)
end


function systemEquation(s::SDF_VD, p::Function)
    function dydt(y::AbstractArray{Float64,1}, t::Float64)
        if s.alpha == 1.0
            dvdt = p(t)-2.0*(s.zta + s.xi)*s.omg*y[2]-s.omg*s.omg*y[1]
        else
            dvdt = p(t)-2.0*s.zta*s.omg*y[2]-2.0*s.xi*s.omg*sign(y[2])*abs(y[2])^s.alpha-s.omg*s.omg*y[1]
        end
        [y[2]; dvdt]
    end
    dydt
end


function response(s::SDF_VD, p::Function, t::AbstractArray{Float64,1})
    y0 = [s.u0; s.v0]

    f = systemEquation(s, p)

    y = ode4(f, y0, t)

    u = y[:,1]
    v = y[:,2]
    fdz = 2.0*s.zta*s.omg*v
    fdx = 2.0*s.xi*s.omg*sign(v).*abs(v).^s.alpha
    a = map(p,t)-fdz-fdx-s.omg*s.omg*u

    [u v a fdz fdx]
end