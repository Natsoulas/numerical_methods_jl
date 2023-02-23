module explicit_euler
#explicit euler method with adaptive time step control 
function euler_adaptive(f, t0, t_end, y0, epsilon, mu)
    # Set initial values
    t = t0
    y = y0
    h = (t_end - t0) / 1000  # initial time step size
    tvals = [t0]
    yvals = [y[1];y[2]]'

    while t < t_end
        # Calculate solution at next time step using explicit Euler method
        y_n_plus_1 = y + h .* f(t, y, mu)

        # Estimate local error using difference between explicit Euler and midpoint methods
        y_mid = y + (h/2) .* f(t, y, mu)
        y_mid_n_plus_1 = y_mid + (h/2) .* f(t + h/2, y_mid, mu)
        err = (y_mid_n_plus_1 - y_n_plus_1)'

        # Adjust time step size
        if maximum(abs, err) < epsilon
            # Accept time step
            t += h
            y = y_n_plus_1
            tvals = [tvals; t]
            yvals = [yvals; y']
        end

        # Increase/decrease time step size if error is small/large enough
        if maximum(abs, err) < epsilon/2
            h *= 2
        elseif maximum(abs, err) > epsilon && t + h < t_end
            h /= 2
        end
    end

    return tvals, yvals
end


end
