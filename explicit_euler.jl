module explicit_euler
#explicit euler method with adaptive time step control 
function euler_adaptive(f, t0, t_end, y0, epsilon; mu=nothing)
    # Set initial values
    t = t0
    y = y0
    h = (t_end - t0) / 1000  # initial time step size
    tvals = [t0]
    yvals = [y0']

    while t < t_end
        # Calculate solution at next time step using explicit Euler method
        y_n_plus_1 = y + h * f(t, y, mu)

        # Estimate local error using difference between explicit Euler and Runge-Kutta methods
        y_rk4 = y + h * f(t, y, mu)
        for i in 1:3
            k = h * f(t + (i-1)/2 * h, y_rk4, mu)
            y_rk4 += k/3
        end
        k4 = h * f(t + h, y_rk4, mu)
        err = maximum(abs.((k4/2 - h * f(t + h/2, y + h/2 * f(t, y, mu), mu))))[:]

        # Adjust time step size
        if err < epsilon
            # Accept time step
            t += h
            y = y_n_plus_1
            tvals = [tvals; t]
            yvals = [yvals; y']
        else
            # Reduce time step size and try again
            h /= 2
        end

        # Increase time step size if error is small enough
        if err < epsilon/2
            h *= 2
        end
    end

    return tvals', yvals'
end
end
