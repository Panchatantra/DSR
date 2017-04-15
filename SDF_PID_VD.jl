type SDF_PID_VD <: SDFDynamicSystem
    omg::Float64
    zta::Float64
    mu::Float64
    xi::Float64
    kpa::Float64
    alpha::Float64
    u0::Float64
    v0::Float64
    ud0::Float64
    vd0::Float64

    SDF_PID_VD(omg::Float64, zta::Float64, mu::Float64, xi::Float64, kpa::Float64, alpha::Float64) = new(omg, zta, mu, xi, kpa, alpha, 0.0, 0.0, 0.0, 0.0)
    SDF_PID_VD(omg::Float64, zta::Float64, mu::Float64, xi::Float64, kpa::Float64) = new(omg, zta, mu, xi, kpa, 1.0, 0.0, 0.0, 0.0, 0.0)
end


function systemEquation(s::SDF_PID_VD, p::Function)
    function dydt(y::AbstractArray{Float64,1}, t::Float64)
        u, v, ud, vd = y
        omg2 = s.omg*s.omg
        dvdt = p(t)-2.0*s.zta*s.omg*v-omg2*u-s.kpa*omg2*(u-ud)
        
        if s.alpha == 1.0
            dvddt = (s.kpa*omg2*(u-ud) - 2.0*s.xi*s.omg*vd)/s.mu
        else
            dvddt = (s.kpa*omg2*(u-ud) - 2.0*s.xi*s.omg*sign(vd)*abs(vd)^s.alpha)/s.mu
        end
        [v; dvdt; vd; dvddt]
    end
    dydt
end

function response(s::SDF_PID_VD, p::Function, t::AbstractArray{Float64,1})
    y0 = [s.u0; s.v0; s.ud0; s.vd0]

    f = systemEquation(s, p)

    y = ode4(f, y0, t)

    u = y[:,1]
    v = y[:,2]
    ud = y[:,3]
    vd = y[:,4]
    fdz = 2.0*s.zta*s.omg*v
    fdx = s.kpa*s.omg*s.omg*(u-ud)
    fdd = 2.0*s.xi*s.omg*sign(vd).*abs(vd).^s.alpha
    a = map(p,t)-fdz-fdx-s.omg*s.omg*u

    [u v a fdz fdx ud fdd]
end