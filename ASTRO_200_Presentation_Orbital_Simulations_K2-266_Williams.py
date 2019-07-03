
# coding: utf-8

# In[1]:


import numpy as np
import astropy
import scipy
import math 
import matplotlib.pyplot as plt
from scipy import constants
from astropy.io import ascii
from astropy import cosmology 
from astropy import units as u 
from astropy.units import astrophys as astro
from vpython import * 
from vpython import graph 


# In[2]:


#window for the progrom if it's run on a computer 
scene = display(title = 'Planets  Orbiting K2 - 266', width = 800, 
                height = 800, range = (25000,25000,25000), center = (1000,0,0)) 


# In[3]:


#getting file from directory 
filename = 'C://Users/Janel/Documents/ASTRO_200_Python/K2-266.txt'

#reading data from the file into an ascii table 
data = ascii.read(filename)



#indexing raw values for mathematical calculation 
planet_letter= data['pl_letter'] #planet letter 
#print(planet_letter)
P = data['pl_orbper'] #planets orbital period in days 
#print(P)
i_deg = data['pl_orbincl'] #planets inclincation 
#print(i_deg)
pl_massj= data['pl_bmassj'] #planets mass in jupiter units
#print(pl_massj)
pl_rad= data['pl_radj'] #planet radius in jupiter radii
#print(pl_rad)
St_dist= data['st_dist'] #distance of star from us in parsecs
#print(St_dist)
St_M= data['st_mass'] #stellar mass in solar mass 
#print(St_M)
St_R= data['st_rad'] #stellar radius in solar radii
#print(St_R)
e = data['pl_orbeccen'] #orbital eccentricity 
#print(e)
semi_ax = data['pl_orbsmax'] #semi - major axis of each planet in AU 
#print(semi_ax)
orb_peri = data['pl_orbper']
V_star = data['st_rad']
#print(orb_peri)


# In[4]:


#defining constants 
G = scipy.constants.gravitational_constant *1e-9 #gravitational constant from scipy database #changed from m^3 to kilometers^3
#print(G)
pi = scipy.constants.pi # pi constant from data bate 
#print(pi)
pi_2 = pi**2 #pi squared 
#print(pi_2)
au_km = (1.496e8) # 1au in kilometers 
#print(au_km)
Mjup_kg= (1.8981872e27) #jupiter masses in kilograms, conversion factor 
#print(Mjup_kg)
Msol_kg = (1.9884754e30) #solar masses in kilograms, conversion factor 
#print(Msol_kg)
Rjup_km = (71492.0) #jupiter radii in kilometers 
#print(Rjup_km)
Rsol_km = (695700.0) #star radii in kilometers 
#print(Rsol_km)
#change from integer to real number 


# In[5]:


#indexing raw extracted raw data 
Mstar = St_M[0] * Msol_kg  #multiplying indexed raw stellar mass data by solar constant 
#print(Mstar)
Rstar = St_R[0]*Rsol_km #indexing Rstar value which is in km 
#print(Rstar)
v_star = V_star[0]
Mp = pl_massj*Mjup_kg #multiplying indexed raw planet masses by jupiter masses to get kg 
#print(Mp)
Rp = pl_rad*Rjup_km #multiplying indexed raw planet radii by jupiter radians to get km 
#print(Rp)
ap = semi_ax*au_km #multiplying indexed raw planet obital radius by 1au in km to get the ap from au to km
#print(ap)
i_rad = (i_deg)*(pi/180)
#print(i_rad)
Wop = orb_peri

#print(Wop)

#print(ap)
#print(planet_letter)

#------------Planet Masses in kg -----------# 
#bmass = Mp[0:1]
bmass = Mp[0]
#print(bmass)
cmass = Mp[1]
#print(cmass)
dmass = Mp[2]
#print(dmass)
emass = Mp[3]
#print(emass)

print(bmass)
print(cmass)
print(dmass)
#----------Planet Period in days  ----------#
P_b = P[0]
#print(P_b)
P_c = P[1]
#print(P_c)
P_d = P[2]
#print(P_d)
P_e = P[3]
#print(P_e)


#--------Planet Orbital Semimajor axis in km ------#
ap_b = ap[0]
#print(ap_b)
ap_c = ap[1]
#print(ap_c)
ap_d = ap[2]
#print(ap_d)
ap_e = ap[3]
#print(ap_e)


#---------------Planet Radii in km ---------------#

Rb = Rp[0]
#print(Rb)
Rc = Rp[1]
#print(Rc)
Rd = Rp[2]
#print(Rd)
Re = Rp[3]
#print(Re)


#--------------Planet Inclincation in radians -------------#
i_b = i_rad[0]
#print(i_b)
i_c = i_rad[1]
#print(i_c)
i_d = i_rad[2]
#print(i_d)
i_e = i_rad[3]
#print(i_e)

#-------------Planet Orbital Eccentricity ---------------#
e_b = .044
#print(e_b)
e_c = e[1]
#print(e_c)
e_d = e[2]
#print(e_d)
e_e = e[3]
#print(e_e)

#----------Planet Long. Periastrian angle in deg----------------# 
Wob = 88.00 #abitrarily chosen due to missing values
#print(Wob)
Woc = Wop[1]
Wod = Wop[2]
#print(Wod)
Woe = Wop[3]
#print(Woe)

#----------------FORCE on the planets------------------#
F_b = (G*(bmass*Mstar))/(ap_b**2)
#print(F_b)
F_c = (G*(cmass*Mstar))/(ap_c**2)
#print(F_c)
F_d = (G*(dmass*Mstar))/(ap_d**2)
#print(F_d)
F_e = (G*(emass*Mstar))/(ap_e**2)
#print(F_e)


#------------Angular Velocity of the planets, the initial velocity-----------------#
wb = math.sqrt(F_b/(bmass*ap_b))
wc = math.sqrt(F_c/(cmass*ap_c))
wd = math.sqrt(F_d/(dmass*ap_d))
we = math.sqrt(F_e/(emass*ap_e))

#------------Initial Velocities of the planets-----------------#
vb = wb*ap_b
vc = wc*ap_c
vd = wd*ap_d
ve = we*ap_e



#print(vb)
#print(vc)
#print(vd)
#print(ve)


# In[6]:


#this defines the inital position, size, and velocity of the objects 

rsc = 10

star = sphere(pos = vector(0,0,0), radius = Rstar, color = color.red, make_trail = True ) # star's starting position at the origin
star.velocity= vector(0,0,0)

b = sphere(pos = vector(ap_b,0,0), radius = Rb*rsc, color = color.orange, make_trail = True ) # star's starting position at the origin of the simulation
b.velocity = vector(0,vb,0)

c = sphere(pos = vector(ap_c,0,0), radius = Rc*rsc, color = color.green, make_trail = True ) # star's starting position at the origin ''
c.velocity = vector(0,vc,0)

d = sphere(pos = vector(ap_d,0,0), radius = Rd*rsc, color = color.blue, make_trail = True ) # star's starting position at the origin ''
d.velocity = vector(0,vd,0)

e = sphere(pos = vector(ap_e,0,0), radius = Re*rsc, color = color.purple, make_trail = True ) # star's starting position at the origin ''
e.velocity = vector(0,ve,0)


# In[7]:


#---------Initializing Graphics Curve -----------------------# 
f1 = gcurve(color=color.orange)	# a graphics curve
f2 = gcurve(color = color.green)
f3 = gcurve(color = color.blue)
f4 = gcurve(color = color.purple)
f5 = gcurve(color = color.black)

t = 0 


# In[8]:


t = 0 #inital time 
dt = 10 #timestep
while t <= 1e50:
    rate(1000)
    dist_b = mag(b.pos)                           
    unvec_b = (b.pos - star.pos)/dist_b
    fgrav_b_mag = -(G*(bmass*Mstar))/(dist_b**2)
    fgrav_b = fgrav_b_mag * unvec_b
    b.velocity = b.velocity + (np.divide(fgrav_b,bmass))*dt 
    b.pos = b.pos + b.velocity*dt 
    
#    f1.plot(t,b.pos.y)
    
    dist_c = mag(c.pos)                           
    unvec_c = (c.pos - star.pos)/dist_c
    fgrav_c_mag = -(G*(cmass*Mstar))/(dist_c**2)
    fgrav_c = fgrav_c_mag * unvec_c
    c.velocity = c.velocity + (np.divide(fgrav_c,cmass))*dt 
    c.pos = c.pos + c.velocity*dt 
    
#    f2.plot(t,c.pos.y)
    
    dist_d = mag(d.pos)
    unvec_d = (d.pos - star.pos)/dist_d
    fgrav_d_mag = -(G*(dmass*Mstar))/(dist_d**2)
    fgrav_d = fgrav_d_mag * unvec_d 
    d.velocity = d.velocity + (np.divide(fgrav_d,dmass))*dt 
    d.pos = d.pos + d.velocity*dt 
    
#    f3.plot(t,d.pos.y)
    
    dist_e = mag(e.pos)
    unvec_e = (e.pos - star.pos)/dist_e
    fgrav_e_mag = -(G*(emass*Mstar))/(dist_e**2)
    fgrav_e = fgrav_e_mag*unvec_e
    e.velocity = e.velocity + (np.divide(fgrav_e,emass))*dt 
    e.pos = e.pos + e.velocity*dt

#    f4.plot(t,e.pos.y)
    
    mratio_b = bmass/Mstar
    mratio_c = cmass/Mstar
    mratio_d = dmass/Mstar
    mratio_e = emass/Mstar
    
    star.pos.y = star.pos.y + (b.velocity.y*mratio_b + c.velocity.y*mratio_c + d.velocity.y*mratio_d + e.velocity.y*mratio_e)*dt
    star.pos.x = star.pos.x + (b.velocity.x*mratio_b + c.velocity.x*mratio_c + d.velocity.x*mratio_d + e.velocity.x*mratio_e)*dt
    star.pos.z = star.pos.z + (b.velocity.z*mratio_b + c.velocity.z*mratio_c + d.velocity.z*mratio_d + e.velocity.z*mratio_e)*dt
    
    #print(star.pos.y)
    #plotting the y compnent of the velocity over time 
    f1.plot(t,b.velocity.y)
    f2.plot(t,c.velocity.y)
    f3.plot(t,d.velocity.y)
    f4.plot(t,e.velocity.y)   
#    f5.plot(t,star.pos)
    f5.plot(t,star.pos.y)
    
    #plotting the magnitiudes of the gravitional force on the planets over time 
    t+=dt

