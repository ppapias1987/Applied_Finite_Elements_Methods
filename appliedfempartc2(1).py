#!/usr/bin/env python
# coding: utf-8

# In[67]:


from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import random


# In[68]:


#Create a mesh and define function space
mesh = Mesh("/home/chris/Desktop/Computational Science/Courses/Applied Finite Element Methods/Project PartC/sweden.xml.gz")


# In[69]:


# Construct the finite element space
P1 = FiniteElement("Lagrange",mesh.ufl_cell(),1)
TH = P1 * P1 * P1  # Taylor Hood Element
V = FunctionSpace(mesh,TH) # Creates space for u, v, w


# In[70]:


#Define parameters:
# Define parameters:
T = 1200
dt = 0.5
delta1 = Constant(1.0)
delta2 = Constant(1.0)
delta3 = Constant(1.0)
alpha = Constant(0.4)
beta = Constant(1)
gamma = Constant(0.8)
zeta = Constant(2.0)
L_0 = Constant(0.4)
l = Constant(0.6)
m = Constant(0.12)
r=random.uniform(0,1)


# In[71]:


# Class representing the initial conditions
class InitialConditions(UserExpression):
    def eval(self,values,x):
        if x[1]<100:
            values[0]=5/1000*r
            values[1]=1/2*(1 - r)
            values[2]=1/4 + 1/2*r
        else :
            values[0]=1/100
            values[1]=1/100
            values[2]=1/100


            
            
          
    def value_shape(self):
        return (3,)


# In[72]:


# Define Test and Trial Functions

#Trial-Initial Conditions
indata = InitialConditions(degree=2)
trial_0=TrialFunction(V) 
trial_0=interpolate(indata, V)

#Test
k=TestFunction(V)

#Trial
trial=TrialFunction(V)

#For 3-PDEs System

#Trial_Initial
uz=trial_0[0]
vz=trial_0[1]
wz=trial_0[2]



#Test
k1=k[0]
k2=k[1]
k3=k[2]

#Trial
u=trial[0]
v=trial[1]
w=trial[2]





# In[73]:


# #Plot solution at t=0

# Mutualists (t=0)
# plot(trial_0[0])
# plt.title("mutualists,t="+str(0))
# plt.savefig("u"+str(0)+".png")
# plt.show()
# plt.clf()

# # Preys (t=0)
# plot(trial_0[1])
# plt.title("preys,t="+str(0))
# plt.savefig("v"+str(0)+".png")
# plt.show()
# plt.clf()

# # Predators (t=0)
# plot(trial_0[2])
# plt.title("predators,t="+str(0))
# plt.savefig("w"+str(0)+".png")
# plt.show()
# plt.clf()


# In[74]:


#Define Variational Problem - System of PDEs
F=(((u-uz)/dt)*k1)*dx \
  +(((v-vz)/dt)*k2)*dx \
  +(((w-wz)/dt)*k3)*dx \
  +(delta1*inner(grad(u+uz),grad(k1)))*0.5*dx \
  +(delta2*inner(grad(v+vz),grad(k2)))*0.5*dx \
  +(delta3*inner(grad(w+wz),grad(k3)))*0.5*dx \
  -(alpha*((u+uz)*0.5)*k1)*dx \
  -(beta*((v+vz)*0.5)*k2)*dx \
  +(gamma*((w+wz)*0.5)*k3)*dx \
  +(alpha*uz**2)/(L_0+l*vz)*k1*dx\
  +beta*(vz**2)*k2*dx \
  +vz*wz/(alpha + vz + m*uz)*k2*dx \
  -(zeta*(wz*vz)/(alpha+vz+m*uz))*k3*dx



# In[75]:


a= lhs(F) 
L=rhs(F)
trial=Function(V)


# In[76]:


upop=[assemble(uz*dx)] 
vpop=[assemble(vz*dx)]
wpop=[assemble(wz*dx)]

time=[]





# In[77]:


t=0
while t < T:
    time.append(t)
    t=t+1

    solve(a==L,trial) 
    trial_0.assign(trial)

    utemp=assemble(trial_0[0]*dx)
    vtemp=assemble(trial_0[1]*dx)
    wtemp=assemble(trial_0[2]*dx)

    upop.append(utemp)
    vpop.append(vtemp)
    wpop.append(wtemp)

     # Plot
    
    # if int(t)%100==0: 

        #Mutualists
        # plot(trial_0[0])
        # plt.title("mutualists,t="+str(t))
        # plt.savefig("u"+str(t)+".png")
        # plt.show()
        # plt.clf()

        # #Preys
    
        # plot(trial_0[1])
        # plt.title("preys,t="+str(t))
        # plt.savefig("v"+str(t)+".png")
        # plt.show()
        # plt.clf()

    #     #Predators
    
    #     plot(trial_0[2])
    #     plt.title("preys,t="+str(t))
    #     plt.savefig("w"+str(t)+".png")
    #     plt.show()
    #     plt.clf()



      





        
        
     

        
       



    

    





# In[80]:


# vpop=vpop[:-1]
# wpop=wpop[:-1]
# upop=upop[:-1]

plt.plot(time,upop,label="mutualists")
plt.plot(time,vpop, label="preys")
plt.plot(time,wpop, label="predators")
plt.legend()
plt.show
plt.savefig("u,v"+str(t)+".png")
plt.clf


# #Phase portrait (preys-predators)
# plt.plot(vpop,wpop, label="preys-predators")
# plt.legend()
# plt.show
# plt.savefig("v,w"+".png")
# plt.clf

# Phase portrait (mutualists-preys-predators)
plt.plot(upop,vpop,wpop)
plt.title("mutualists-preys-predators")
plt.legend()
plt.grid()
plt.show
plt.savefig("u,v,w"+".png")
plt.clf




# In[79]:


# Mutualists (t=1000)
# plot(trial_0[0])
# plt.title("mutualists,t="+str(1000))
# plt.savefig("u"+str(t)+".png")
# plt.show()
# plt.clf()

# # Preys (t=1000)
# plot(trial_0[1])
# plt.title("preys,t="+str(1000))
# plt.savefig("v"+str(t)+".png")
# plt.show()
# plt.clf()


# In[ ]:




