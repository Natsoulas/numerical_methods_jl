module implicit_euler

using LinearAlgebra, Roots

using LinearAlgebra, Roots

function implicit_euler_adaptive_step_doubling(f, y0, t0, tf, rtol, dt_max, dt_min, mu=nothing)
    # Set initial values
    t = t0
    y = copy(y0)
    dt = dt_max

    # Initialize output arrays
    ts = [t0]
    ys = [y0']

    while t < tf
        # Compute derivative at current time and state
        f_current = f(y, t, mu)

        # Use Newton's method to solve for the next state
        function implicit_func(y_new)
            return y_new - y - dt * f(y_new, t+dt, mu)
        end
        y_new = find_zero(implicit_func, y, Order1(), verbose=false)

        # Compute error and update time step
        e = norm(y_new - y) / (rtol * (norm(y) + norm(y_new)))
        if e > 1
            dt = 0.9 * dt * (1/e)^0.5
        elseif e < 0.5 || dt < dt_min
            dt = min(dt_max, 2*dt)
            t += dt
            y = copy(y_new)
        else
            dt = min(dt_max, 2*dt)
            t += dt
            y = copy(y_new)
        end

        # Update output arrays
        push!(ts, t)
        push!(ys, y')
    end

    return ts, ys
end
end