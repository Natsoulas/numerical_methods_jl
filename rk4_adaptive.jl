# RK4 method with adaptive time step control
function rk4_adaptive(f, t0, tf, x0, er_tol)
    #f is function to be integrated/derivative: f(t,x)
    #er_tol is tolerance for error for adaptive time step control
    t = t0
    x = x0
    h = (tf - t0) / 100 # initial time step size
    while t < tf
        k1 = h * f(t, x)
        k2 = h * f(t + h/2, x + k1/2)
        k3 = h * f(t + h/2, x + k2/2)
        k4 = h * f(t + h, x + k3)
        x1 = x + (k1 + 2*k2 + 2*k3 + k4)/6
        k1 = h * f(t+h, x1)
        k2 = h * f(t + h/2, x1 + k1/2)
        k3 = h * f(t + h/2, x1 + k2/2)
        k4 = h * f(t + h, x1 + k3)
        x2 = x + (k1 + 2*k2 + 2*k3 + k4)/6
        err = norm(x2 - x1) # error estimate
        if err < er_tol
            t += h
            x = x2
        end
        h = 0.8 * h * (er_tol/err)^(1/5) # adjust time step size
    end
    return x
end
