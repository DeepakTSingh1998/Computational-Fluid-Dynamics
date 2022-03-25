# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt 

## Initial conditons ##
u = 2 # velcoity
time = 1 # total time
length = 5 # total length of space domain
length_cos = 1 # length of initial cos stop
Nodes = 100
Nodes_cos = int(np.round(Nodes - (length-length_cos)*(Nodes/length))) # Nodes required to cos inpulse
x = np.linspace(0,length,Nodes)
dt = time/(Nodes-1)
dx = length/(Nodes-1)

# Cos impulse
phi_0 = np.zeros(Nodes)
phi_cos = np.cos(np.linspace((np.pi/2),-(np.pi/2),Nodes_cos))
phi = phi_0.copy()
phi[:len(phi_cos)] += phi_cos

# alpha 
alpha = (u*dt)/(2*dx)
CLF = alpha*2

## Alpha matrix ##
Matrix = np.zeros((Nodes,Nodes))
# Boundary Conditons
Matrix[0,0] = 1
Matrix[Nodes-1,Nodes-1] = 1

for i in range(0,Nodes-2):
    Matrix[i+1,i] = -alpha
    Matrix[i+1,i+1] = 1
    Matrix[i+1,i+2] = alpha
    
Matrix_inv = np.linalg.inv(Matrix) # Inverse of Matrix

## Main Process ##
t = 0 #intial velocity
phi_new = phi
while t < time:
    phi_new = np.dot(Matrix_inv,phi_new)
    
    t += dt
    plt.clf()
    plt.plot(x,phi_new)
    plt.grid(True)
    plt.xlabel('Distance')
    plt.ylabel('Phi')
    plt.xlim(0,length)
    plt.ylim(0,2)
    plt.title('1D Scaler Convection t =' + str(t))
    plt.pause(0.1)

