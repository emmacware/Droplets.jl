##################################################
# Condensation  example
##################################################
using Droplets
using DifferentialEquations
using Plots

# the differential equations for condensation from "condensation.jl" take radius as the argument and 
#have the following options for parameters:
#drdtcondensation1 p = (a,b,S,M,denom)
#drdtcondensation2 p = (M,m,T,qv,P)
#drdtcondensation3 p = (M,m,T,Senv)

#where S is environmental saturation, M is the mass of the solute,
#m is the molecular weight of the solute, T is the temperature, 
#qv is the mixing ration, and P is the pressure



#here is a test example with drdtcondensation1

####SETTINGS####
constants_ = Constants{Float32}()

T = 273 # Temperature in K
S = 1.01 #Saturation
smallS = 0.95 #Saturation < Sactivation for these settings
M = 1e-16 #Mass of solute
a = 3.3*10^(-7)/T #(m)
m = 58.44 #molecular weight of NaCL
b = 4.3 *2/m/1e6 #m^3 for NaCL
denom = FK(T,constants_)+FD(T,constants_) #denominator for the Köhler equation

#parameters for ODE problem
pcondense = (a,b,S,M,denom)
peq = (a,b,smallS,M,denom)
Rc = 7.4e-8
tspan = (0.00001, 100.0)

rODE = ODEProblem{false}(drdtcondensation1,Rc,tspan,pcondense)
rnew = solve(rODE,Rosenbrock23())
R = solve(rODE,ImplicitEuler())

req = ODEProblem{false}(drdtcondensation1,Rc,tspan,peq)
reqRosenbrock = solve(req,Rosenbrock23())
reqEuler = solve(req,ImplicitEuler())


Köhler(r) = 1+a/r-(b*M)/r^3
rarray = 10 .^ range(-7.3, stop=-5, length=100)


############GRAPH#########



p1 = plot(rarray*1e6,Köhler.(rarray),xscale=:log10,label="Köhler")
p1 = scatter!(R.u*1e6,Köhler.(R.u),xscale=:log10,label="ImplicitEuler")
p1 = scatter!(rnew.u*1e6,Köhler.(rnew.u),xscale=:log10,label="Rosenbrock23")
p1 = title!("Köhler Curve")
p1 = xlabel!("Radius (μm)")
p1 = ylabel!("Saturation")
p1 = Köhler(6e-8)
rmax = sqrt(3*b*1e-16/a)
Sact = Köhler(rmax)

#annotate the plot
p1 = plot!([5e-2,10],[S,S],color="black",label="Senv")
p1 = plot!([Rc*1e6,Rc*1e6],[0.97,.99],color="red",label="Ri")
p1 = plot!([rmax*1e6,rmax*1e6],[0.995,1.01],color="blue",label="Ract")
p1 = plot!([R.u[end]*1e6,R.u[end]*1e6],[0.995,1.01],color="orange",label="Euler"*string(tspan[2])*"seconds")
p1 = plot!([rnew.u[end]*1e6,rnew.u[end]*1e6],[0.995,1.01],color="green",label="Rosenbrock"*string(tspan[2])*"seconds")
p1 = plot!([0.05,rmax*1e6],[Sact,Sact],color="pink",label="Sact")

p2 = scatter(R.t,R.u*1e6,color="orange",xscale=:log10,label="ImplicitEuler,  "*string(length(R.t))*"steps")
p2 = scatter!(rnew.t,rnew.u*1e6,color="green",xscale=:log10,label="Rosenbrock23,  "*string(length(R.t))*"steps")
p2 = plot!(legend=:topleft)
p2 = plot!([50,200],[R.u[end]*1e6,R.u[end]*1e6],color="orange",label="Euler"*string(tspan[2])*"seconds")
p2 = plot!([50,200],[rnew.u[end]*1e6,rnew.u[end]*1e6],color="green",label="Rosenbrock"*string(tspan[2])*"seconds")
p2 = plot!([0.00001,200],[rmax*1e6,rmax*1e6],color="blue",label="Ract")
p2 = annotate!(1e-3, 1.5, text("S="*string(S)*",Ri=6e-8", 10, :left))
p2 = annotate!(8,R.u[end]*1e6, text(string(round((R.u[end]*1e6),digits=2)), 10, :left,:orange))
p2 = annotate!(8,rnew.u[end]*1e6, text(string(round((rnew.u[end]*1e6),digits=2)), 10, :left,:green))
p2 = title!("Condensation Growth")
p2 = xlabel!("Time (s)")
p2 = ylabel!("Radius (μm)")

p3 = plot(rarray*1e6,Köhler.(rarray),xscale=:log10,label="Köhler")
p3 = scatter!(reqEuler.u*1e6,Köhler.(reqEuler.u),xscale=:log10,label="ImplicitEuler")
p3 = scatter!(reqRosenbrock.u*1e6,Köhler.(reqRosenbrock.u),xscale=:log10,label="Rosenbrock23")
p3 = title!("Köhler Curve")
p3 = xlabel!("Radius (μm)")
p3 = ylabel!("Saturation")
rmax = sqrt(3*b*1e-16/a)
Sact = Köhler(rmax)
p3 = plot!([5e-2,10],[smallS,smallS],color="black",label="Senv")
p3 = plot!([Rc*1e6,Rc*1e6],[0.97,.99],color="red",label="Ri")
p3 = plot!([rmax*1e6,rmax*1e6],[0.995,1.01],color="blue",label="Ract")
p3 = plot!([reqRosenbrock.u[end]*1e6,reqRosenbrock.u[end]*1e6],[0.94,0.96],color="green",label="Rosenbrock"*string(tspan[2])*"seconds")
p3 = plot!([reqEuler.u[end]*1e6,reqEuler.u[end]*1e6],[0.94,0.96],color="orange",label="Euler"*string(tspan[2])*"seconds")
p3 = plot!([0.05,rmax*1e6],[Sact,Sact],color="pink",label="Sact")

p4 = scatter(reqEuler.t,reqEuler.u*1e6,color="orange",xscale=:log10,label="ImplicitEuler,  "*string(length(reqEuler.t))*"steps")
p4 = scatter!(reqRosenbrock.t,reqRosenbrock.u*1e6,color="green",xscale=:log10,label="Rosenbrock23,  "*string(length(reqRosenbrock.t))*"steps")
p4 = plot!(legend=:topleft)
p4 = plot!([50,200],[reqRosenbrock.u[end]*1e6,reqRosenbrock.u[end]*1e6],color="green",label="Rosenbrock"*string(tspan[2])*"seconds")
p4 = plot!([50,200],[reqEuler.u[end]*1e6,reqEuler.u[end]*1e6],color="orange",label="Euler"*string(tspan[2])*"seconds")
p4 = plot!([0.00001,200],[rmax*1e6,rmax*1e6],color="blue",label="Ract")
p4 = annotate!(1e-3, 1.5, text("S="*string(smallS)*",Ri=6e-8", 10, :left))
p4 = title!("Condensation Growth")
p4 = xlabel!("Time (s)")
p4 = ylabel!("Radius (μm)")



title1 = plot(title = "T ="*string(T)*", Ms(NaCL) ="*string(M), grid = false, showaxis = false)
title2 = plot(title = "Sact<Senv", grid = false, showaxis = false)
title3 = plot(title = "Sact>Senv", grid = false, showaxis = false)
plot(title1,title2,p1,p2,title3,p3,p4,layout=@layout([A{0.01h}; B{0.01h};[C D];E{0.01h};[F G]]),size=(1000,1000),left_margin=5Plots.mm,right_margin=5Plots.mm,top_margin=1Plots.mm,bottom_margin=5Plots.mm)#,dpi=300)

# savefig("Kohlertime.pdf")#,dpi=300)
