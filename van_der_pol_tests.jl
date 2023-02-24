import Pkg
Pkg.add("Plots")
include("explicit_euler.jl")
#import implicit_euler
#import rk4
#using implicit_euler
#using rk4

#test numerical methods on Van Der Pol Problem for different values of mu (1, 3, 10, 100)
using Plots
function VanDerPol(t,Y,μ)
    # VANDERPOL Implementation of the Van der Pol model
    # Syntax: Ydot = VanDerPol(t,Y,mu)
    #Ydot = zeros(2,1)
    Ydot_1 = Y[2]
    Ydot_2 = μ*(1-Y[1]*Y[1])*Y[2]-Y[1]
    Ydot = [Ydot_1; Ydot_2]
    return Ydot
end    
μ = 10
#stiff is when it has a high mu value
y0 = [2.0; 0.0]
tol = 1.0E-6

#try explicit euler test
T,Y = explicit_euler.euler_adaptive(VanDerPol,0.0,5.0*μ,y0,tol,μ)


display(plot(T,Y[:,1]))


display(plot(T,Y[:,2]))


display(plot(Y[:,1],Y[:,2]))

