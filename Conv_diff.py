# Import libraries
import numpy as np
import matplotlib.pyplot as plt


# Local to Global stiffness matrx
def Gather(IN1,OUT1,IND,LOCAL,type_array):
    if type_array == 1: #if == to matrix
        for i in range(0,len(IN1)):
            for j in range(0,len(IN1)):
                LOCAL[OUT1[i],OUT1[j]]  += IND[IN1[i],IN1[j]]
    else: #if not matrix, then vector
        for i in range(0,len(IN1)):
            LOCAL[OUT1[i]] += IND[IN1[i]]
    return LOCAL


# Python Matrix function (Minus 1 in dictionary)
def Matrix(DICT):
    for i in range(0,len(DICT)):
        DICT[i] = DICT[i]-1
    return DICT 


# Inputs 
u = 0
k = 1
L = 2
Nodes = 150
h = L/(Nodes-1)
Pe = (u*h)/(2*k)
Q = 0
print('The Pecklet number is '  + str(np.round(Pe,3)))
ndof = 1


# Boundary conditions
phi = []
for i in range(0,Nodes):
    phi.append("Free") 
phi[0] = 0 #phi @ 0
phi[Nodes-1] = 1 #phi @ L


# Stiffness and loadvectors
K = np.zeros((Nodes,Nodes))
F = np.zeros((Nodes))
for i in range(0,Nodes-1):
    ke = (u/2)*np.array(([-1,1],[-1,1]))+(k/h)*np.array(([1,-1],[-1,1]))
    fe = (Q/2)*np.array(([1],[1]))
    IN = Matrix({0:1,1:2}) 
    OUT = {0:1,1:2}
    for j in range(0,len(OUT)):
        OUT[j] = ndof*i+j
    K = Gather(IN,OUT,ke,K,1)
    F = Gather(IN,OUT,fe,F,2)

A = K
G = F
c = 0
phi_loc = []


# Condense matrix
for i in range(0,Nodes):
    if phi[i] == "Free":
        phi_loc.append(i)
    else:
        G += -phi[i]*A[:,i-c]
        A = np.delete(A,i-c,axis=0)
        A = np.delete(A,i-c,axis=1)
        G = np.delete(G,i-c,axis=0)
        c += 1


# Solve missing phi values
phiv = np.dot(np.linalg.inv(A),G)

c = 0
for i in range(0,Nodes):
    if i == phi_loc[c]:
        phi[i] = phiv[c]
        c += 1
        if c == len(phi_loc):
            c = 0
    else:
        phi[i] = phi[i]


# Post Processing
x = np.linspace(0,L,Nodes)
plt.plot(x,phi)
plt.xlabel('x')
plt.ylabel('phi')
plt.title('1D Convection Diffusion, Pe = '+str(Pe))
plt.grid(True)
        

        