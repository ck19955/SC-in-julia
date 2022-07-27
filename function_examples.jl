function ode_first_order(t, u)
    return u
end


function ode_second_order(t, u)
    u_array = [u[2], -u[1]]
    return u_array
end

function exact_second_order(t, u)
    x = u[0]
    y = u[1]
    u_array = [x*cos(t) + y*sin(t), -x*sin(t) + y*cos(t)]
    return u_array
end

function exponential(t, u)
    return exp(t)
end

function pred_prey(t, u_values, a, b, d)
    x = u_values[1]
    y = u_values[2]
    u_array = [x*(1-x) - (a*x*y)/(d+x), b*y*(1-(y/x))]
    return u_array
end

function hopf_bif(t, u_values, β)
    x = u_values[1]
    y = u_values[2]
    u_array = [β*x - y - x*(x^2 + y^2), x + β*y - y*(x^2 + y^2)]
    return u_array
end

function exact_hopf_bif(t, u_values, β, θ)
    u_array = [(β^0.5)*cos(t + θ), (β^0.5)*sin(t + θ)]
    return u_array
end

function alg_cubic(t, x, c)
    return x^3 - x + c
end

function mod_hopf_bif(t, u_values, β)
    x = u_values[1]
    y = u_values[2]
    u_array = [β*x - y + x*(x^2 + y^2) - x*(x^2 + y^2)^2,
                        x + β*y + y*(x^2 + y^2) - y*(x^2 + y^2)^2]
    return u_array
end

function u_I(x, l)
    # initial temperature distribution
    y = sin(π*x/l)
    return y
end

function alternate_u_I(x, l)
    # initial temperature distribution
    y = (sin(π*x/l))^2
    return y
end

function u_exact(x, t, k, l)
    # the exact solution to the temperature equation
    y = exp(-k*(π^2/l^2)*t)*sin(π*x/l)
    return y
end

function p(t)
    return 5
end

function q(t)
    return 3
end