#import sys
#sys.path.insert(0, './nuclide-data')
#from nuclide_data import *
from nuclideData.nuclide_data import *
import matplotlib.pyplot as plt
import matplotlib.colors as cl
from scipy.optimize import curve_fit
import numpy as np
import uncertainties as unc
import re

massP = 938.27208816 # Mev/c^2
massN = 939.56542052 # Mev/c^2
toMev = 931.493614838934 # ua -> Mev

def fetch_data(A,Z):
    # how many points [left]
    nL = 2
    # how many points [right]
    nR = 2

    Y=[]
    Y_e=[]
    X=[]

    print(Z)
    # get nuclides [left]
    for i in range(1,nL+1):
       X.append(Z-i)
       # get mass + error (as strings)
       txt = str(nuc(Z-i,A)['weight'])
       # split mass and error
       Msp = re.split("/",txt)[0][:-1]
       Msp_e = re.split("/",txt)[1][1:]
       # convert in MeV
       sperimentalMass = float(Msp)*toMev
       sperimentalMass_e = float(Msp_e)*toMev
       #theoreticalMass = (Z-i)*massP + (A-(Z-i))*massN
       #b.append((theoreticalMass-sperimentalMass)/A) Y.append(sperimentalMass)
       # store data
       Y.append(sperimentalMass)
       Y_e.append(sperimentalMass_e)

    # get nuclides [right]
    for i in range(0,nR+1):
       X.append(Z+i)
       txt = str(nuc(Z+i,A)['weight'])
       Msp = re.split("/",txt)[0][:-1]
       Msp_e = re.split("/",txt)[1][1:]
       sperimentalMass = float(Msp)*toMev
       sperimentalMass_e = float(Msp_e)*toMev
       #sperimentalMass = weight(Z+i,A)*toMev
       #theoreticalMass =  (Z+i)*massP + (A-(Z+i))*massN
       #b.append((theoreticalMass-sperimentalMass)/A)
       Y.append(sperimentalMass)
       Y_e.append(sperimentalMass_e)
    return (X,Y,Y_e)

def SEMF(Z,A,AV,AS,AC, AA):
    #b=AV*A-AS*A**(2/3) - AC*(Z**2)/A**(1/3) - AA*((A-Z -Z)**2)/A
    #b=(AV-AA)*A - AS*A**(2./3) + 4.*AA*Z - (AC*A**(-1./3)+4.*AA*A**(-1))*Z**2
    b=AV*A-AS*A**(2/3) - AC*Z*(Z-1)*A**(-1/3) - AA*(A-2*Z)**2/A
    return (Z*massP + (A-Z)*massN) - b

def main():
    #print(nuc(44,101))
    #print(weight(44,101))
    zed = 44 # Z di partenza
    ei = 101 # A di partenza

    # get data
    data = fetch_data(ei,zed)
    X=data[0]
    Y=data[1]
    Y_e=data[2]

    # fit
    par,cov=curve_fit(lambda Z,AV,AS,AC,AA: SEMF(Z,ei,AV,AS,AC,AA),
                      X,Y,
                      sigma=Y_e, absolute_sigma=True,
                      #p0= (15.6,17.23, 0.7, 23.28),
                      bounds=([14,16,0,17],[17,19,2,25]),)
    # plot
    fig=plt.figure()
    plt.errorbar(X,Y,Y_e,fmt='.')
    plt.plot(np.linspace(zed-zed*10/100,zed+zed*10/100),SEMF(np.linspace(zed-zed*10/100,zed+zed*10/100),ei, *par))
    plt.show()

    # print results
    print(par)
    print(cov)
    print(X)
    print(Y)

# Call the main function when running the script
if __name__ == "__main__":
    main()
