import numpy as np
import pylab as plt
from math import pi, log, exp, sqrt, log
from scipy.integrate import odeint
from scipy.constants import inch
from scipy.integrate import quad

import matplotlib
matplotlib.rcParams.update({'font.size': 20})


# controls:
# rtime, Ea*, block_size, P1

gamma = 1.4
R_ = 8.31  # universal gas constant J/mol/degK
M = 0.029  # molecular weight of air kg/mol
R = R_/M   # specific gas constant J/kg/degK


rho0 = 15.5 # kg/m3
p0 = 1e5
gamma = 1.4
p1 = 1e6

# %precision %.3g
e0 = p0/rho0/(gamma-1)
print "e0", e0
t0 = p0*M/rho0/R_
print "t0", t0
t1 = p1*M/rho0/R_
print "t1", t1

Q = p1/rho0/(gamma-1) - e0
print "Q", Q

rad = 1.5*inch/2.0
h = w = l = 0.25
rock_rho = 2200.0
explosive_volume = h*pi*rad**2
explosive_mass = explosive_volume * rho0
rock_mass = (h*l*w - explosive_volume)* rock_rho


explosive_volume_atmo = (p1*explosive_volume**gamma/p0)**(1/gamma)
print "expansion factor", explosive_volume_atmo/explosive_volume

p = lambda v : max(0, p1 * explosive_volume**gamma/v**gamma-p0)
work = quad(p, explosive_volume, explosive_volume_atmo)[0]
block_vel = sqrt(2*work/rock_mass)
print block_vel, "m/s"



Ea = R_ * t0 * 5.0

end_time = 1.2e-1
t = [0.0] + np.logspace(-5, np.log10(end_time), 1000).tolist()
#res = odeint(f,0,t)
rt = [1e-4, 1e-3, 2e-3]
if True:
    plt.ylabel("Pressure [MPa]")
    plt.xlabel("Time [s]")
    my_A = -log(0.5)/rt[0]
    f = lambda x,t : my_A*(1-x)*exp(-Ea/((gamma-1)*M*(e0+x*Q)))
    res = odeint(f,0,t)
    press = (e0+Q*res)*rho0*(gamma-1)
    plt.plot(t,press/1e6,linewidth=2)
    my_A = -log(0.5)/rt[1]
    f = lambda x,t : my_A*(1-x)*exp(-Ea/((gamma-1)*M*(e0+x*Q)))
    res = odeint(f,0,t)
    press = (e0+Q*res)*rho0*(gamma-1)
    plt.plot(t,press/1e6,linewidth=2)
    my_A = -log(0.5)/rt[2]
    f = lambda x,t : my_A*(1-x)*exp(-Ea/((gamma-1)*M*(e0+x*Q)))
    res = odeint(f,0,t)
    press = (e0+Q*res)*rho0*(gamma-1)
    plt.plot(t,press/1e6,linewidth=2)
    plt.legend(("Rise time {:.1e} s".format(rt[0]),
                "Rise time {:.1e} s".format(rt[1]),
                "Rise time {:.1e} s".format(rt[2])),loc=4)
    plt.show()

def deriv(x,t):
    # x[0] is x position
    # x[1] is x' vel
    # x[2] is ev # specific internal energy
    # x[3] lambda
    pos       = x[0]
    vel       = x[1]
    e         = x[2]
    lam       = x[3]

    # define density and ev
    volume = pi*(rad+pos)**2*h
    density = explosive_mass/(volume)
    pressure = e * density * (gamma-1)
    temp = pressure/density/R
    dlamdt = A*(1-lam)*exp(-Ea/R_/temp)
    dvdt = pi*h*(2*rad + 2*pos)*vel
    dedt = dlamdt*Q - pressure/explosive_mass*dvdt
    lp = max((pressure-p0),0)
    #lp = pressure
    #if lp<p0: lp=0.0
    #area = 4*sqrt(2)*rad
    area = pi*2*(rad+pos)*h
    acceleration = area *lp /rock_mass
    return [vel,
            acceleration,
            dedt,
            dlamdt]
A=1
def model(rtime, end_time=0.5e-1):
    global A
    A = -log(0.5)/rtime

    time = [0.0] + np.logspace(-5, np.log10(end_time), 1000).tolist()
    x0 = [0.0, 0.0, e0, 0.0]
    res  = odeint(deriv, x0, time)
    #print res
    x,v,e,lam = res.T
    volume = pi*(rad+x)**2*h
    density = explosive_mass/(volume)
    pressure = (e) * density * (gamma-1)
    # print pressure[-1]*M/density[-1]/R_   # end temperature
    pressure[pressure<p0] = p0
    return time, x,v,pressure,lam, volume

# _tl  = np.load("time.npy")
# _uyl = np.load("vel.npy")
# _al  = np.load("ap.npy")
# _pl  = np.load("pl.npy")


time0, x0, v0, pressure0, lam0, vol0 = model(rt[0])
time1, x1, v1, pressure1, lam1, vol1 = model(rt[1])
time2, x2, v2, pressure2, lam2, vol2 = model(rt[2])


plt.loglog(vol0, pressure0, linewidth=2)
plt.loglog(vol1, pressure1, linewidth=2)
plt.loglog(vol2, pressure2, linewidth=2)
plt.ylabel("Log$_{10}$ Pressure [Pa]")
plt.xlabel("Log$_{10}$ Volume [m$^3$]")
plt.legend(("Rise time {:.1e} s".format(rt[0]),
            "Rise time {:.1e} s".format(rt[1]),
            "Rise time {:.1e} s".format(rt[2])),loc=1)

plt.show()


my_A = -log(0.5)/rt[0]
f = lambda x,t : my_A*(1-x)*exp(-Ea/((gamma-1)*M*(e0+x*Q)))
res0 = odeint(f,0,time0)
press0 = (e0+Q*res0)*rho0*(gamma-1)

my_A = -log(0.5)/rt[1]
f = lambda x,t : my_A*(1-x)*exp(-Ea/((gamma-1)*M*(e0+x*Q)))
res1 = odeint(f,0,time0)
press1 = (e0+Q*res1)*rho0*(gamma-1)

my_A = -log(0.5)/rt[2]
f = lambda x,t : my_A*(1-x)*exp(-Ea/((gamma-1)*M*(e0+x*Q)))
res2 = odeint(f,0,time0)
press2 = (e0+Q*res2)*rho0*(gamma-1)

plt.subplot(311)
plt.plot(time0,v0)
plt.plot(time1,v1)
plt.plot(time2,v2)
#plt.plot(_tl+2e-4, _uyl)
plt.axhline(block_vel, linestyle="--", color="black")
plt.ylabel("Block vel [m/s]")
plt.subplot(312)
#plt.ylim(None,P1)
plt.plot(time0,pressure0/1e6)
plt.plot(time1,pressure1/1e6)
plt.plot(time2,pressure2/1e6)
plt.plot(time0,press0/1e6, "--", color="blue")
plt.plot(time1,press1/1e6, "--", color="green")
plt.plot(time2,press2/1e6, "--", color="red")
#plt.plot(_tl+2e-4,_pl/1e6)
plt.ylabel("Pressure [MPa]")
plt.subplot(313)
plt.ylabel("$\lambda$ []")
plt.plot(time0,lam0)
plt.plot(time1,lam1)
plt.plot(time2,lam2)
plt.plot(time0,res0, "--", color="blue")
plt.plot(time1,res1, "--", color="green")
plt.plot(time2,res2, "--", color="red")

plt.xlabel("Time [s]")
plt.show()
