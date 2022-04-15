import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

## Variables
L = 150                      #Length (m)
H = 150                       #Height (m)
W = 150                        #Width (m)
                #Temperature at the bottem on the Area (C)
T_Surface = 0                  #Temperature of the surface (C)
T_plume = 1500                 #Temperature of the Plume (C)
k = 0.5                       #Thermal Diffusivity (m2/s)


## Numerical Parameters
nx = 101                       #Gridpoints (x)
ny = 101                       #Gridpoints (y)
nt = 50                    #Total Seconds (s)
dt = 1                         #Timestep (s)
dx = L/(nx-1)                  #Distance between grids (x)
dy = H/(ny-1)                  #Distance between grids (y)
x = np.linspace(-L/2,L/2,nx)   #(x) Coordinates
y = np.linspace(-H,0,ny)       #(y) Coordinates
xm,ym = np.meshgrid(x,y)       #(x,y) Meshgrid


# Checking if the stability condition is met
Limit = min([pow(dx,2),pow(dy,2)])/(4*k)
if Limit < dt:
    sys.exit('Error due to Time step being too large')
    
# Creating T
T = np.empty((int(nt),nx,ny))
T.fill(T_Surface)

ind = np.argwhere(abs(xm[0,:])<=(W/2))
T[0,0,ind] = T_plume

## Compute Temperatures
def Temperature(T):
    for t in  range(0,nt-1):
        for i in range(1,nx-1):
            for j in range(1,ny-1):
                S_x = (k*dt)/pow(dx,2)
                S_y = (k*dt)/pow(dy,2)
                T[t+1,i,j] = T[t,i,j]\
                    + S_x*(T[t,i,j+1] - 2*T[t,i,j] + T[t,i,j-1])\
                    + S_y*(T[t,i+1,j] - 2*T[t,i,j] + T[t,i-1,j])
                    
    return T

T = Temperature(T)                    

def Plot(T,t):
    plt.clf()
    plt.title(f"Time = {t:.1f} (s)")
    plt.contourf(xm,ym,T,40,cmap='hot')
    plt.colorbar()
    plt.xlabel('X axis')
    plt.ylabel('Y axis')
                    
    return plt

def animate(t):
  Plot(T[t], t)
  
anim = animation.FuncAnimation(plt.figure(), animate, interval=dt, frames=nt, repeat=False)
anim.save("heat_equation_solution.gif")
    

        

                    
                    
    
                    
            