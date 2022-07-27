using Plots
include("function_examples.jl")


function euler_step(t, x, step_size, ode)
    x0 = x + step_size * ode(t, x)
    return x0
end


function rk4(t, x, step_size, ode)
    """
    Executes a single step of the 4th Order Runge Kutta method for given value, t_n
    Parameters:
    ----------
        t_n : float
            The value of the independent variable
        x_n : numpy array
            The value of the dependant variable(s)
        step_size : float
            The step-size of the rk4 step to be executed
        ode : function
            The ODE for which the rk4 step predicts
        args : numpy array
            The parameters of the ODE
    Returns:
    ----------
        x : numpy array
            The new value of the dependant variable(s) after an rk4 step
"""
    k1 = ode(t, x)
    k2 = ode(t + step_size/2, x + k1*(step_size/2))
    k3 = ode(t + step_size/2, x + k2*(step_size/2))
    k4 = ode(t + step_size, x + k3*step_size)
    x = x + ((k1 + 2*k2 + 2*k3 + k4)/6)*step_size
    return x
end


function solve_to(t_0, t_end, x_0, deltaT_max, method, ode)
    """
    Solves an ODE for an array of values of the independent variable, t.
    Parameters:
    ----------
        t_0 : float
            The initial value of the independent variable
        t_end : float
            The final value of the independent variable
        x_0 : numpy array
            The value of the dependant variable(s)
        deltaT_max : float
            The maximum value the step-size can take when performing a step of a numerical method
        ode : function
            The ODE to be integrated
        args : numpy array
            The parameters of the ODE
    Returns:
    ----------
        x : numpy array
            The final value of the dependant variable(s)
    """

    x = x_0
    t = t_0
    while t < t_end
        if t + deltaT_max <= t_end
            x = method(t, x, deltaT_max, ode)
            t += deltaT_max
        else
            deltaT_max = t_end - t  # Reduces deltaT_max to fit next iteration perfectly
            # This coding decision was chosen to ensure the final time value, t_end is accounted for
        end
    end
    return x
end



function solve_ode(t_values, x_0, deltaT_max, method, ode)
    """
    Solves an ODE from an initial value, t_0 to a final value t_end. However, the difference
    between the two values must be smaller than deltaT_max.
    Parameters:
    ----------
        t_values : array
            values of the independent variable
        x_0 : numpy array
            The initial value of the dependant variable(s)
        deltaT_max : float
            The maximum value the step-size can take when performing a step of a numerical method
        method : function
            The one-step integration method to be used on the ODE
        ode : function
            The ODE to be integrated
        args : numpy array
            The parameters of the ODE
    Returns:
    ----------
        x : numpy array
            An array of dependant variable values for the integrated ODE
    """

    x_values = Array{Float64}(undef, (length(t_values), length(x_0)))   # Create a list of solution x values

    for j in 1:(length(x_0))
        x_values[1, j] = x_0[j]
    end
    for i in 1:(length(t_values)-1)
        x_values[i+1, :] .= solve_to(t_values[i], t_values[i + 1], x_values[i, :], deltaT_max, method, ode)
    end
    return x_values
end

t_values = LinRange(1, 10, 100)
x_values = solve_ode(t_values, [3, 4], 0.01, rk4, ode_second_order)

plot(t_values, x_values[:, 1])