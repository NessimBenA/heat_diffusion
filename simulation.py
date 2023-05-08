import numpy as np
import matplotlib.pyplot as plt
import math

# Parameters
m = 100 # m : distance where it is diffused
time_step = 1 # ro
space_step = 1 # delta
viscosite=0.02
prandlt=0.71
coeffision=viscosite/prandlt
relaxation_time = 1/(3*viscosite + 0.5) #wm 
relaxation_time2 = 1/(3*coeffision + 0.5) #wd
cs=3 #1/cs^2 de nos calculs
constantev0 = 0.2
c=1

U0 = 1 # U_0
iterations = 20000 # number of time step simulated
rho=time_step/space_step

# Init
p = 5*np.ones((m+1,m+1)) # p(x,t)
u0 = p*(4/9) #f0
u1 = p*(1/9) #f1
u2 = p*(1/9)
u3 = p*(1/9)
u4 = p*(1/9)
u5 = p*(1/36)
u6 = p*(1/36)
u7 = p*(1/36)
u8 = p*(1/36)

x = np.zeros((m+1,m+1)) # X(x,t)
w0=(4/9)*x #g1
w1=(1/9)*x #g2
w2=(1/9)*x
w3=(1/9)*x
w4=(1/9)*x
w5=(1/36)*x
w6=(1/36)*x
w7=(1/36)*x
w8=(1/36)*x


v1 = np.zeros((m+1,m+1)) #Vitesse vx
v1[0,:]=constantev0
v2 = np.zeros((m+1,m+1)) #vitesse vy
tk=([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36])
c0=[0,0]
c1=[1,0]
c2=[0,1]
c3=[-1,0]
c4=[0,-1]
c5=[1,1]
c6=[-1,1]
c7=[-1,-1]
c8=[1,-1]

def functioncalculf(cx): #fonction pour calculer feq
    if cx==0:
        return tk[cx]* p * (1 + cs * (c0[0] * v1 + c0[1] * v2) + (1/2)*((cs**2) * (c0[0] * v1 + c0[1] * v2)**2) - (1/2 * cs * (c0[0] * v1 + c0[1]*v2)))
    if cx==1:
        return tk[cx]* p * (1 + cs * (c1[0] * v1 + c1[1] * v2) + (1/2)*((cs**2) * (c1[0] * v1 + c1[1] * v2)**2) - (1/2 * cs * (c1[0] * v1 + c1[1]*v2)))
    if cx==2:
        return tk[cx]* p * (1 + cs * (c2[0] * v1 + c2[1] * v2) + (1/2)*((cs**2) * (c2[0] * v1 + c2[1] * v2)**2) - (1/2 * cs * (c2[0] * v1 + c2[1]*v2)))
    if cx==3:
        return tk[cx]* p * (1 + cs * (c3[0] * v1 + c3[1] * v2) + (1/2)*((cs**2) * (c3[0] * v1 + c3[1] * v2)**2) - (1/2 * cs * (c3[0] * v1 + c3[1]*v2)))
    if cx==4:
        return tk[cx]* p * (1 + cs * (c4[0] * v1 + c4[1] * v2) + (1/2)*((cs**2) * (c4[0] * v1 + c4[1] * v2)**2) - (1/2 * cs * (c4[0] * v1 + c4[1]*v2)))
    if cx==5:
        return tk[cx]* p * (1 + cs * (c5[0] * v1 + c5[1] * v2) + (1/2)*((cs**2) * (c5[0] * v1 + c5[1] * v2)**2) - (1/2 * cs * (c5[0] * v1 + c5[1]*v2)))
    if cx==6:
        return tk[cx]* p * (1 + cs * (c6[0] * v1 + c6[1] * v2) + (1/2)*((cs**2) * (c6[0] * v1 + c6[1] * v2)**2) - (1/2 * cs * (c6[0] * v1 + c6[1]*v2)))
    if cx==7:
        return tk[cx]* p * (1 + cs * (c7[0] * v1 + c7[1] * v2) + (1/2)*((cs**2) * (c7[0] * v1 + c7[1] * v2)**2) - (1/2 * cs * (c7[0] * v1 + c7[1]*v2)))
    if cx==8:
        return tk[cx]* p * (1 + cs * (c8[0] * v1 + c8[1] * v2) + (1/2)*((cs**2) * (c8[0] * v1 + c8[1] * v2)**2) - (1/2 * cs * (c8[0] * v1 + c8[1]*v2)))




# Process temperature diffusion simulation
for i in range(iterations+1):
    
    u0eq = functioncalculf(0)
    u1eq = functioncalculf(1)
    u2eq = functioncalculf(2)
    u3eq = functioncalculf(3)
    u4eq = functioncalculf(4)
    u5eq = functioncalculf(5)
    u6eq = functioncalculf(6)
    u7eq = functioncalculf(7)
    u8eq = functioncalculf(8)
    
    # Collision
    u0 = u0*(1-relaxation_time) + relaxation_time*u0eq
    u1 = u1*(1-relaxation_time) + relaxation_time*u1eq
    u2 = u2*(1-relaxation_time) + relaxation_time*u2eq
    u3 = u3*(1-relaxation_time) + relaxation_time*u3eq
    u4 = u4*(1-relaxation_time) + relaxation_time*u4eq
    u5 = u5*(1-relaxation_time) + relaxation_time*u5eq
    u6 = u6*(1-relaxation_time) + relaxation_time*u6eq
    u7 = u7*(1-relaxation_time) + relaxation_time*u7eq
    u8 = u8*(1-relaxation_time) + relaxation_time*u8eq
    
    # Streaming
    u1[:,1:] = u1[:,:-1]
    u2[:-1,:] = u2[1:,:]
    u3[:,:-1] = u3[:,1:]
    u4[1:,:] = u4[:-1,:]
    u5[:-1,1:] = u5[1:,:-1]
    u6[:-1,:-1] = u6[1:,1:]
    u7[1:,:-1] = u7[:-1,1:]
    u8[1:,1:] = u8[:-1,:-1]
    
   
    # Border conditions
    u1[:,0] = u3[:,0]
    u5[:,0] = u7[:,0]
    u8[:,0] = u6[:,0]


    u3[:,m] = u1[:,m]
    u7[:,m] = u5[:,m]
    u6[:,m] = u8[:,m]


    u2[m,:] = u4[m,:]
    u5[m,:] = u7[m,:]
    u6[m,:] = u8[m,:]


    p[0,:] = u0[0,:] + u1[0,:] + u3[0,:] + 2*(u2[0,:] + u6[0,:] + u5[0,:])
    u4[0,:] = u2[0,:]
    u7[0,:] = u5[0,:] - (p[0,:]*constantev0)/6
    u8[0,:] = u6[0,:] + (p[0,:]*constantev0)/6

    p = u0 + u1 + u2 + u3 + u4 + u5 + u6 + u7 + u8 

  
    v1=(u1*c1[0] + u2*c2[0] + u3*c3[0] + u4*c4[0] + u5*c5[0] + u6*c6[0] + u7*c7[0] + u8*c8[0])/p
    v2=(u1*c1[1] + u2*c2[1] + u3*c3[1] + u4*c4[1] + u5*c5[1] + u6*c6[1] + u7*c7[1] + u8*c8[1])/p


    
    #on passe sur g 
    w0eq= (4/9)*x
    w1eq= (1/9)*x*(1+3*v1)
    w2eq= (1/9)*x*(1+3*v2)
    w3eq= (1/9)*x*(1-3*v1)
    w4eq= (1/9)*x*(1-3*v2)
    w5eq= (1/36)*x*(1+3*v1+3*v2)
    w6eq= (1/36)*x*(1-3*v1+3*v2)
    w7eq= (1/36)*x*(1-3*v1-3*v2)
    w8eq= (1/36)*x*(1+3*v1-3*v2)
    
    # Collision
    w0 = w0*(1-relaxation_time2) + relaxation_time2*w0eq
    w1 = w1*(1-relaxation_time2) + relaxation_time2*w1eq
    w2 = w2*(1-relaxation_time2) + relaxation_time2*w2eq
    w3 = w3*(1-relaxation_time2) + relaxation_time2*w3eq
    w4 = w4*(1-relaxation_time2) + relaxation_time2*w4eq
    w5 = w5*(1-relaxation_time2) + relaxation_time2*w5eq
    w6 = w6*(1-relaxation_time2) + relaxation_time2*w6eq
    w7 = w7*(1-relaxation_time2) + relaxation_time2*w7eq
    w8 = w8*(1-relaxation_time2) + relaxation_time2*w8eq
    
    # Streaming
    w1[:,1:] = w1[:,:-1]
    w2[:-1,:] = w2[1:,:]
    w3[:,:-1] = w3[:,1:]
    w4[1:,:] = w4[:-1,:]
    w5[:-1,1:] = w5[1:,:-1]
    w6[:-1,:-1] = w6[1:,1:]
    w7[1:,:-1] = w7[:-1,1:]
    w8[1:,1:] = w8[:-1,:-1]

    
    # Border conditions
    # Left border
    w1[:,0] = -w3[:,0]
    w5[:,0] = -w7[:,0]
    w8[:,0] = -w6[:,0]
    w0[:,0]=0

    # right border
    w3[:,m] = -w1[:,m]
    w6[:,m] = -w8[:,m]
    w7[:,m] = -w5[:,m]
    w0[:,m]=0

    # top border
    w4[0,:] = (U0/9) + (U0/9) - w2[0,:]
    w7[0,:] = (U0/36) + (U0/36) - w5[0,:]
    w8[0,:] = (U0/36) + (U0/36) - w6[0,:]
    w0[0,:]=U0


    # Left border
    w1[-1,:] = w1[-2,:]
    w2[-1,:] = w2[-2,:]
    w3[-1,:] = w3[-2,:]
    w4[-1,:] = w4[-2,:]
    w5[-1,:] = w5[-2,:]
    w6[-1,:] = w6[-2,:]
    w7[-1,:] = w7[-2,:]
    w8[-1,:] = w8[-2,:]
    w0[-1,:] = w0[-2,:]

    x=w1+w2+w3+w4+w5+w6+w7+w8+w0
    


# Plot
fig = plt.figure(figsize=(12,4))
plt.imshow(x, cmap="plasma")
plt.colorbar()
plt.legend()
plt.show()