type SDF_MYD_EP <: SDFDynamicSystem
    omg::Float64
    zta::Float64
    kpa::Float64
    phi::Float64
    alpha::Float64
    u0::Float64
    v0::Float64
    hv::AbstractArray{Float64,1}

    SDF_MYD_EP(omg::Float64, zta::Float64, kpa::Float64, phi::Float64, alpha::Float64) = new(omg, zta, kpa, phi, alpha, 0.0, 0.0, zeros(5))
    SDF_MYD_EP(omg::Float64, zta::Float64, kpa::Float64, phi::Float64) = new(omg, zta, kpa, phi, 0.0, 0.0, 0.0, zeros(5))
end


function fd(s::SDF_MYD_EP, u::Float64; update::Bool=false)
    up, dup, fp, kp, status = s.hv[1:5]
    k0 = s.kpa*s.omg*s.omg
    k1 = s.alpha*k0
    fy = s.phi*g
    uy = fy/k0
    du = u - up
    f_try = fp+du*kp
    bup = fy+k1*(u-uy)
    bdn = -fy+k1*(u+uy)

    if du == 0.0
        kc = kp
        f = fp
    elseif status == 0
        if du > 0.0
            if f_try > bup
                f = bup
                kc = k1
                status = 1
            else
                kc = k0
                f = f_try
            end
        else
            if f_try < bdn
                f = bdn
                kc = k1
                status = 1
            else
                kc = k0
                f = f_try
            end
        end
    elseif status == 1
        if dup*du > 0
            f = f_try
            kc = k1
        else
            kc = k0
            f = fp + kc*du
            status = 0
            if f > bup
                f = bup
                kc = k1
                status = 1
            elseif f < bdn
                f = bdn
                kc = k1
                status = 1
            end
        end
    end
    if update
        s.hv[1:5] = [u, du, f, kc, status]
    end
    f
end


function systemEquation(s::SDF_MYD_EP, p::Function)
    function dydt(y::AbstractArray{Float64,1}, t::Float64)
        u = y[1]
        v = y[2]
        omg2 = s.omg*s.omg
        dvdt = p(t)-2.0*s.zta*s.omg*v-fd(s,u)-omg2*u
        [v; dvdt]
    end
    dydt
end


function response(s::SDF_MYD_EP, p::Function, t::AbstractArray{Float64,1})
    y0 = [s.u0; s.v0]

    f = systemEquation(s, p)
    
    updateHistVar = (y) -> fd(s, y[1]; update=true)
    y = ode4hv(f, y0, t, updateHistVar)

    u = y[:,1]
    v = y[:,2]
    fdz = 2.0*s.zta*s.omg*v
    fdp = y[:,3]
    a = map(p,t)-fdz-fdp-s.omg*s.omg*u

    [u v a fdz fdp]
end