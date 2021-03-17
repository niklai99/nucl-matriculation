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
    nL = 3
    # how many points [right]
    nR = 3

    Y=[]
    Y_e=[]
    X=[]

    #print(Z)

    # get nuclides [left]
    for i in range(1,nL+1):

       X.append(Z-i)
       # get mass + par_error (as strings)
       txt = str(nuc(Z-i,A)['weight'])
       # split mass and par_error
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
       # get mass + par_error (as strings)
       txt = str(nuc(Z+i,A)['weight'])
       # split mass and par_error
       Msp = re.split("/",txt)[0][:-1]
       Msp_e = re.split("/",txt)[1][1:]
       # convert in MeV
       sperimentalMass = float(Msp)*toMev
       sperimentalMass_e = float(Msp_e)*toMev

       #sperimentalMass = weight(Z+i,A)*toMev
       #theoreticalMass =  (Z+i)*massP + (A-(Z+i))*massN
       #b.append((theoreticalMass-sperimentalMass)/A)

       Y.append(sperimentalMass)
       Y_e.append(sperimentalMass_e)

    # return Z, MASS, MASS_ERROR
    return X, Y, Y_e


def SEMF(Z,A,AV,AS,AC,AA):

    #b=AV*A-AS*A**(2/3) - AC*(Z**2)/A**(1/3) - AA*((A-Z -Z)**2)/A
    #b=(AV-AA)*A - AS*A**(2./3) + 4.*AA*Z - (AC*A**(-1./3)+4.*AA*A**(-1))*Z**2

    b = AV*A - AS*A**(2/3) - AC*Z*(Z-1)*A**(-1/3) - AA*(A-2*Z)**2/A

    return (Z*massP + (A-Z)*massN) - b


def make_fit(X, Y, Y_e, A_start):

    # fit
    par, cov=curve_fit(
                        lambda Z,AV,AS,AC,AA: SEMF(Z,A_start,AV,AS,AC,AA),
                        xdata = X, ydata = Y,
                        sigma=Y_e, absolute_sigma=True,
                        #p0= (15.6,17.23, 0.7, 23.28),
                        bounds=([14,16,0,17],[17,19,2,25]),
                  ) 

    residuals = []
    chi2 = 0
    for i in range(len(X)):
        residuals.append(Y[i] - SEMF(X[i], A_start, *par))
        chi2 += (residuals[i] / Y_e[i])**2

    # get parameters error from cov matrix
    par_error = []

    for i in range(len(par)):
        try:
            par_error.append(np.absolute(cov[i][i])**0.5)
        except:
            par_error.append( 0.00 )

    fit_par = par
    fit_err = np.array(par_error)

    return fit_par, fit_err, chi2


def make_plot(X, Y, Y_e, A_start, par, par_e, chi2):

    # figure config
    fig=plt.figure(figsize=(10,6))
    ax = fig.add_subplot(1, 1, 1)

    # set plot range limits 
    XMIN = np.amin(X) - 1
    XMAX = np.amax(X) + 1

    # text template for parameters
    textA  = 'A  = ' + format(A_start, '1.0f')
    textAV = 'AV = ' + format(par[0], '1.1f') + ' +/- ' + format(par_e[0],'1.1f')
    textAS = 'AS = ' + format(par[1], '1.1f') + ' +/- ' + format(par_e[1],'1.1f')
    textAC = 'AC = ' + format(par[2], '1.4f') + ' +/- ' + format(par_e[2],'1.4f')
    textAA = 'AA = ' + format(par[3], '1.3f') + ' +/- ' + format(par_e[3],'1.3f')
    chisq  = '$\chi^{2}$ / ndf = ' + format(chi2, '1.0f') + ' / ' + format(len(X) - len(par), '1.0f')

    # plot data
    ax.errorbar(x = X, y = Y, yerr = Y_e, marker = '.', markersize = 10,
                elinewidth=1, color = '#000000', linewidth=0, capsize=2, label = 'Data')

    # plot fit
    ax.plot(np.arange(XMIN, XMAX, 0.01), SEMF(np.arange(XMIN, XMAX, 0.01), A_start, *par), 
            color = '#FF4B00', linewidth = 2, linestyle = '-', label = 'Fit')

    # display texts
    ax.text(0.40, 0.50, textA + '\n\n' + textAV + '\n' + textAS + '\n' + textAC + '\n' + textAA + '\n\n' + chisq, fontsize = 14, color = '#000000', transform = ax.transAxes)

    # plot title
    ax.set_title('SEMF FIT')

    # axis labels
    ax.set_xlabel('Z')
    ax.set_ylabel('Mass (MeV)')

    # axis range
    ax.set_xlim(XMIN, XMAX)

    # do not use fancy notations for axis ticks
    ax.ticklabel_format(useOffset=False)

    return fig, ax


def main():

    #print(nuc(44,101))
    #print(weight(44,101))

    Z_start = 40 # Z di partenza
    A_start = 91 # A di partenza

    # get data
    X, Y, Y_e = fetch_data(A_start,Z_start)

    par, par_e, chi2 = make_fit(X, Y, Y_e, A_start)

    fig, ax = make_plot(X, Y, Y_e, A_start, par, par_e, chi2)

    # print results
    print(par)
    print(par_e)

    #print(X)
    #print(Y)

    plt.show()
    
    return


# Call the main function when running the script
if __name__ == "__main__":
    main()
