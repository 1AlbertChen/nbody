import sys
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

N =int(sys.argv[1])     # Number of bodies
M= 1 # Total Mass
a= 1 # Plummer Radius
X, Y, Z, VX, VY, VZ, MASS = range(7)



def getV(r): #Velocity Distribution
    while(1==1):
        c=np.random.rand()
        d=0.1*np.random.rand()
        if((c*c)*np.power(1-c*c, 7/2)<d):
            continue
        else:
            return (c*np.sqrt(2)*np.power(1+r*r, -1/4))
        
        
    
def main():
    bodies = np.zeros((N, 7))
    for i in np.arange(N):
        while(1==1):
            r=np.power(np.random.rand(),1/3) #Randomize radius
            plum = (3*np.power(1+(r*r)/(a*a),-5/2))/(4*np.pi*np.power(a,3)) #Plummer density equation
            if (plum<np.power(np.random.rand(),1/3)): #Rejection sampling
                continue
            else:
                v=getV(r)
                #Randomize location vector direction
                Phi=2*np.pi*np.random.rand() 
                Theta=np.arccos((2*np.random.rand())-1)
                #Randomize velocity vector direction
                VPhi=2*np.pi*np.random.rand()
                VTheta=np.arccos((2*np.random.rand())-1)
                
                bodies[i, X] = r * np.sin(Theta) * np.cos (Phi)
                bodies[i, Y] = r * np.sin(Theta) * np.sin (Phi)
                bodies[i, Z] = r * np.cos(Theta) 
            
                bodies[i, VX] =v * np.sin(VTheta) * np.cos (VPhi)
                bodies[i, VY] =v * np.sin(VTheta) * np.sin (VPhi)
                bodies[i, VZ] =v * np.cos(VTheta) 

                #Uniform mass distribution
                bodies[i, MASS]=M/N

    bodies = pd.DataFrame(bodies)
    bodies.columns = 'x', 'y', 'z', 'vx', 'vy', 'vz', 'mass'
    bodies.to_csv('products_globular/init.csv', index=False)


if __name__ == '__main__':
    
    main()
