import sys
sys.path.insert(0, './nuclide-data')
from nuclide_data import *
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import numpy as np
import uncertainties as unc

massP = 938.27208816 # Mev/c^2
massN = 939.56542052 # Mev/c^2
toMev = 931.493614838934 # ua -> Mev

def fit(A,Z):
    # how many points [left]
    nL = 2
    # how many points [right]
    nR = 2

    # retrieve a few isotones given A
    zi=[]
    b=[]

    # get nuclides [left]
    # and retrieve masses
    for i in range(1,nL+1):
       zi.append(Z-i)
       sperimentalMass = weight(Z-i,A)*toMev
       theoreticalMass = (Z-i)*massP + (A-(Z-i))*massN
       #b.append((theoreticalMass-sperimentalMass)/A)
       b.append(sperimentalMass)

    # get nuclides [right]
    # and retrieve masses
    for i in range(0,nR+1):
       zi.append(Z+i)
       sperimentalMass = weight(Z+i,A)*toMev
       theoreticalMass =  (Z+i)*massP + (A-(Z+i))*massN
       #b.append((theoreticalMass-sperimentalMass)/A)
       b.append(sperimentalMass)
    return (zi,b)

def SEMF(Z,A):
    AV=15.76
    AS=17.81
    AC=0.711
    AA=23.702
    #b=AV*A-AS*A**(2/3) - AC*(Z**2)/A**(1/3) - AA*((A-Z -Z)**2)/A
    b=(AV-AA)*A - AS*A**(2./3) + 4.*AA*Z - (AC*A**(-1./3)+4*AA*A**(-1))*Z**2
    return   (Z*massP + (A-Z)*massN) - b

def main():
    #print(nuc(44,101))
    #print(weight(44,101))
    X=fit(101,44)[0]
    Y=fit(101,44)[1]
    fig=plt.figure()
    plt.plot(X,Y,'.')
    plt.plot(np.linspace(40,50),SEMF(np.linspace(40,50),101))
    plt.show()
    #print(fit(101,44))

# Call the main function when running the script
if __name__ == "__main__":
    main()
