# supertype for a system of ODEs. Different ODE systems can be defined
# by creating subtypes, which should implement in-place versions of
# rhs and (optionally) jac, representing the differential equations and
# jacobian respectively. They should also implement "state_dims", giving the
# phase space dimension for a particular instance
abstract type ODESys end;

# dummy functions to define the ODESys interface
function rhs(dx, x, p::ODESys, t=0.0) ::Nothing end
function jac(J, x, p::ODESys, t=0.0) ::Nothing end
function state_size(sys::ODESys) ::Tuple end

# don't attempt broadcast over an ODESys
# facilitates broadcasting with "." over methods that accept ODESys alongisde
# other arguments you DO want to boradcast over
broadcastable(s::ODESys) = Ref(s)

# generic out-of-place versions of rhs & jac, for convience
function rhs(x, p::ODESys, t=0.0)
    dxdt = similar(x)
    rhs(dxdt, x, p, t)
    return dxdt
end

function jac(x, p::ODESys, t=0.0)
    n = length(x)
    J = similar(x, n, n)
    jac(J, x, p, t)
    return J
end
