module DSR

export SDF, SDF_VD, SDF_FD, SDF_MYD_BW, SDF_MAXWELL_VD, response

const g = 9.8


function ode4(f::Function, y0::AbstractArray{Float64,1}, t::AbstractArray{Float64,1})
    # solve ODE dy/dt = f(y,t) using 4th-order Runge Kutta Method
    nstep = length(t)
    neq   = length(y0)

    y  = zeros(nstep, neq)
    y[1,:] = y0

    yn = copy(y0)

    for n in 2:nstep
        tn = t[n-1]
        h = t[n] - tn
        h2 = h*0.5

        k1 = f(yn      , tn   )
        k2 = f(yn+h2*k1, tn+h2)
        k3 = f(yn+h2*k2, tn+h2)
        k4 = f(yn+h*k3 , tn+h )

        yn += h/6.0*(k1+2.0*k2+2.0*k3+k4)
        y[n,:] = yn
    end

    y
end


function ode4hv(f::Function, y0::AbstractArray{Float64,1}, t::AbstractArray{Float64,1}, updateHistVar::Function)
    # solve ODE dy/dt = f(y,t,hv) using 4th-order Runge Kutta Method
    nstep = length(t)
    neq   = length(y0)

    y  = zeros(nstep, neq)
    y[1,:] = y0

    yn = copy(y0)
    fd = zeros(nstep)

    for n in 2:nstep
        tn = t[n-1]
        h = t[n] - tn
        h2 = h*0.5

        k1 = f(yn      , tn   )
        k2 = f(yn+h2*k1, tn+h2)
        k3 = f(yn+h2*k2, tn+h2)
        k4 = f(yn+h*k3 , tn+h )

        yn += h/6.0*(k1+2.0*k2+2.0*k3+k4)
        y[n,:] = yn
        fd[n] = updateHistVar(yn)
    end

    [y fd]
end


function readtxt(fn::String; skiprows::Int64=0, delim::Char=' ')
    data = readdlm(fn, delim; skipstart=skiprows)
    if size(data,2) == 1
        return data[:,1]
    elseif size(data,2) == 2
        return data[:,1], data[:,2]
    else
        return data
    end
end


function ds(t::Float64,dt::Float64,a_g::AbstractArray{Float64,1})
    # Convert a decrete signal into a continuous one via linear interpolation
    ind = convert(Int64, floor(t/dt)) + 1
    if ind == length(a_g)
        res = a_g[ind]
    elseif ind > length(a_g)
        res = 0.0
    elseif floor(t/dt)==ceil(t/dt)
        res = a_g[ind]
    else
        al = a_g[ind]
        ar = a_g[ind + 1]
        k = (ar - al) / dt
        return al+k*(t-(ind-1)*dt)
    end

    res
end

abstract DynamicSystem
abstract SDFDynamicSystem <: DynamicSystem


include("SDF.jl")
include("SDF_VD.jl")
include("SDF_FD.jl")
include("SDF_MYD_BW.jl")
include("SDF_MYD_EP.jl")
include("SDF_MAXWELL_VD.jl")
include("SDF_PID_VD.jl")


function response(s::SDFDynamicSystem, p::Function, t0::Float64, t1::Float64, dt::Float64)
    t = t0:dt:t1
    [t response(s, p, t)]
end


function response(s::SDFDynamicSystem, p::Function, t1::Float64, dt::Float64)
    t = 0.0:dt:t1
    [t response(s, p, t)]
end


function response(s::SDFDynamicSystem, p::Function, dt::Float64, nstep::Int64)
    t = 0.0:dt:dt*(nstep-1)
    [t response(s, p, t)]
end


function response(s::SDFDynamicSystem, a_g::AbstractArray{Float64,1}, t::AbstractArray{Float64,1})  
    dt_g = t[2] - t[1]
    p = (t) -> ds(t, dt_g, -a_g)
    [t response(s, p, t)]
end


function response(s::SDFDynamicSystem, a_g::AbstractArray{Float64,1}, dt_g::Float64, t::AbstractArray{Float64,1})  
    p = (t) -> ds(t, dt_g, -a_g)
    [t response(s, p, t)]
end


function response(s::SDFDynamicSystem, a_g::AbstractArray{Float64,1}, dt_g::Float64, t0::Float64, t1::Float64, dt::Float64)  
    p = (t) -> ds(t, dt_g, -a_g)
    response(s, p, t0, t1, dt)
end


function response(s::SDFDynamicSystem, a_g::AbstractArray{Float64,1}, dt_g::Float64, t1::Float64, dt::Float64)  
    p = (t) -> ds(t, dt_g, -a_g)
    response(s, p, t1, dt)
end


function response(s::SDFDynamicSystem, a_g::AbstractArray{Float64,1}, dt_g::Float64, dt::Float64, nstep::Int64)  
    p = (t) -> ds(t, dt_g, -a_g)
    response(s, p, dt, nstep)
end

function response(s::SDFDynamicSystem, a_g::AbstractArray{Float64,1}, t1::Float64, dt::Float64)  
    p = (t) -> ds(t, dt, -a_g)
    response(s, p, t1, dt)
end


function response(s::SDFDynamicSystem, a_g::AbstractArray{Float64,1}, dt::Float64, nstep::Int64)  
    p = (t) -> ds(t, dt, -a_g)
    response(s, p, dt, nstep)
end


end