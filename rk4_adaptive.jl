module rk4
# RK4 method with adaptive time step control
function rk4_adaptive(f, t0, tf, y0, tol, mu=nothing)
    t = [t0]
    y = [y0]
    h = (tf - t0) / 100 # initial time step size
    while t[end] < tf
        k1 = h * f(t[end], y[end], mu)
        k2 = h * f(t[end] + h/2, y[end] + k1/2, mu)
        k3 = h * f(t[end] + h/2, y[end] + k2/2, mu)
        k4 = h * f(t[end] + h, y[end] + k3, mu)
        y1 = y[end] + (k1 + 2*k2 + 2*k3 + k4)/6
        k1 = h * f(t[end] + h, y1, mu)
        k2 = h * f(t[end] + h/2, y1 + k1/2, mu)
        k3 = h * f(t[end] + h/2, y1 + k2/2, mu)
        k4 = h * f(t[end] + h, y1 + k3, mu)
        y2 = y[end] + (k1 + 2*k2 + 2*k3 + k4)/6
        err = norm(y2 - y1) # error estimate
        if err < tol
            t_new = t[end] + h
            y_new = y2
            push!(t, t_new)
            push!(y, y_new)
        end
        h = 0.8 * h * (tol/err)^(1/5) # adjust time step size
    end
    return t, y
end



end