using LinearAlgebra, Roots

function implicit_euler_adaptive_step_doubling(f, y0, t0, tf, rtol, dt_max, dt_min)
    #dt_max is max time step, dt_min is the minimum, rest is implied, rtol is relative error tolerance
    # Set initial values
    t = t0
    y = y0
    dt = dt_max

    while t < tf
        # Compute derivative at current time and state
        f_current = f(y, t)

        # Use Newton's method to solve for the next state
        function implicit_func(y_new)
            return y_new - y - dt * f(y_new, t+dt)
        end
        y_new = similar(y)
        for i in eachindex(y)
            y_new[i] = find_zero(implicit_func, y[i], Order1(), verbose=false)
        end

        # Compute error and update time step
        e = norm(y_new - y) / (rtol * (norm(y) + norm(y_new)))
        if e > 1
            dt = 0.9 * dt * (1/e)^0.5
        elseif e < 0.5 || dt < dt_min
            dt = min(dt_max, 2*dt)
            t += dt
            y = y_new
        else
            dt = min(dt_max, 2*dt)
            t += dt
            y = y_new
        end
    end

    return y
end
