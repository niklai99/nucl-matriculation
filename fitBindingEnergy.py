from nuclideData.nuclide_data import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import re
from operator import itemgetter
from itertools import groupby

massP = 938.27208816 # Mev/c^2
massP_e = 0.00000029 # Mev/c^2
massN = 939.56542052 # Mev/c^2
massN_e = 0.00000054 # Mev/c^2
toMev = 931.493614838934 # ua -> Mev/c^2

ZMAX = 119


def getIsotopes():

    data = []

    # loop over Z
    for Z in range(0,ZMAX):

        iso = isotopes[Z]

        # loop over isotopes for current Z
        for j in range(0,len(iso)):

            # get N from A of current isotope
            A = iso[j]
            N = iso[j] - Z
        
            data.append( (Z, A) )

    # take second element for sort
    def takeSecond(elem):
        return elem[1]
    
    # sort by mass number A
    data.sort(key=takeSecond)
    
    return data


def computeB(data):

    par = []

    for i in range(0, len(data)):

        if 'weight' in nuc(data[i][0], data[i][1]):
            if nuc(data[i][0], data[i][1])['weight'] != '':
    
                # get sperimental mass
                txt = str( nuc(data[i][0], data[i][1])['weight'] )

                # string splitting
                Msp = re.split("/",txt)[0][:-1]
                Msp_e = re.split("/",txt)[1][1:]

                # sperimental mass value and error
                sperimentalMass = float(Msp)*toMev
                sperimentalMass_e = float(Msp_e)*toMev

                # compute theoretical mass 
                theoreticalMass = data[i][0]*massP + (data[i][1]-data[i][0])*massN

                # compute theoretical mass error
                theoreticalMass_e = np.sqrt( (data[i][0]*massP_e)**2 + ((data[i][1]-data[i][0])*massN_e)**2 )

                # compute B/A (MeV per nucleon)
                BA = (theoreticalMass - sperimentalMass) / data[i][1]
                
                # compute B/A error
                BA_e = np.sqrt( theoreticalMass_e**2 + sperimentalMass_e**2 ) / data[i][1]

                # append data
                par.append( (data[i][0], data[i][1], BA, BA_e) )
    
    return par


def computeMean(par):

    # remove duplicates from As
    As = list(set(x[1] for x in par))
    # group Zs energy by A
    Zs = [[x for x,y,z,w in g] for k, g in  groupby(par,key=itemgetter(1))]
    # group binding energy by A
    BAs = [[z for x,y,z,w in g] for k, g in  groupby(par,key=itemgetter(1))]
    # group weights of binding energy by A
    W = [[w**-2 for x,y,z,w in g] for k, g in  groupby(par,key=itemgetter(1))]

    BA_avg = []
    BA_avg_e = []
    Z_avg = []

    # compute weighted average and standard deviation
    for i in range(len(BAs)):
        Z_avg.append( np.average( Zs[i] ) )
        w_avg = np.average( BAs[i], weights=W[i] )
        BA_avg.append(w_avg)
        w_avg_e = np.average( (np.array(BAs[i]) - w_avg)**2, weights=W[i] )
        BA_avg_e.append(w_avg_e)

    return As, Z_avg, BA_avg, BA_avg_e


def BINDING(A, AV, AS, AC, AA):

    B = AV*A[0] - AS*A[0]**(2/3) - AC*A[1]*(A[1]-1)*A[0]**(-1/3) - AA*(A[0]-2*A[1])**2/A[0]

    return B/A[0]


def makeFit(X, Zs, Y, Y_e):
    
    # fit
    par, cov=curve_fit(
                        BINDING,
                        xdata = [X, Zs], ydata = Y,
                        sigma=Y_e, absolute_sigma=True,
                        #p0= (15.6,17.23, 0.7, 23.28),
                        #bounds=([14,16,0,17],[17,19,2,25])
                    )

    # residuals and chi2
    residuals = np.array(Y) - BINDING([np.array(X), np.array(Zs)], *par)
    chi2 = np.sum((residuals / np.array(Y_e))**2)

    # get parameters error from cov matrix
    par_error = np.zeros(len(par))

    for i in range(len(par)):
        try:
            par_error[i] = cov[i][i]**0.5
        except:
            par_error[i] = 0.00

    print(par)
    print(par_error)
    
    return par, par_error, chi2


def makeFig(A, BoverA, par, par_e, chi2):

    # figure config
    fig=plt.figure(figsize=(10,6))
    ax = fig.add_subplot(1, 1, 1)

    # set plot range limits 
    XMIN = np.amin(A) - 1
    XMAX = np.amax(A) + 1
    YMIN = np.amin(BoverA) - 0.5
    YMAX = np.amax(BoverA) + 0.5


    # axis range
    ax.set_xlim(XMIN, XMAX)
    ax.set_ylim(YMIN, YMAX)

    # plot title
    ax.set_title('Binding Energy Trend')

    # axis labels
    ax.set_xlabel('Mass number A')
    ax.set_ylabel('B/A (MeV per nucleon)')

    # do not use fancy notations for axis ticks
    ax.ticklabel_format(useOffset=False)

    # text template for parameters
    textAV = '$A_V$ = ' + format(par[0], '1.5f') + ' $\pm$ ' + format(par_e[0],'1.5f') + ' MeV'
    textAS = '$A_S$ = ' + format(par[1], '1.5f') + ' $\pm$ ' + format(par_e[1],'1.5f') + ' MeV'
    textAC = '$A_C$ = ' + format(par[2], '1.5f') + ' $\pm$ ' + format(par_e[2],'1.5f') + ' MeV'
    textAA = '$A_A$ = ' + format(par[3], '1.5f') + ' $\pm$ ' + format(par_e[3],'1.5f') + ' MeV'
    chisq  = '$\chi^{2}$ / ndf = ' + format(chi2, '1.0f') + ' / ' + format(len(A) - len(par), '1.0f')

    # display texts
    ax.text(0.70, 0.70, textAV + '\n' + textAS + '\n' + textAC + '\n' + textAA + '\n\n' + chisq, fontsize = 18, color = '#000000', transform = ax.transAxes)

    return fig, ax


def plotData(ax, A, BoverA, BoverA_e):
    
    # make the plot
    ax.errorbar(x = A, y = BoverA, yerr = BoverA_e, marker = '.', markersize = 10, elinewidth=1, color = '#000000', linewidth=0, capsize=2, label = 'Data', zorder = 0)

    return


def plotFit(ax, A, Z, par, FITSLICE):

    XMAX = np.amax(A) + 1

    ax.plot(
        np.arange(A[FITSLICE], XMAX, 0.1), 
        BINDING([ np.arange(A[FITSLICE], XMAX, 0.1), np.arange( Z[FITSLICE], ZMAX, 0.1 * (ZMAX-Z[FITSLICE]) / (XMAX-A[FITSLICE]) ) ], *par), 
        color = '#FF4B00', linewidth = 3, linestyle = '-', label = 'Fit', zorder = 1
        ) 

    return 


def main():

    # get data
    data = getIsotopes()

    # compute binding energy 
    par = computeB(data)

    # compute weighted average for isotopes with same As
    As, Z_avg, BA_avg, BA_avg_e = computeMean(par)

    # starting A for plotting
    PLOTSLICE = 20

    # starting A for fitting
    FITSLICE = 50

    # slices
    A_plot = As[PLOTSLICE:]
    BA_plot = BA_avg[PLOTSLICE:]
    BA_e_plot = BA_avg_e[PLOTSLICE:]

    A_fit = As[FITSLICE:]
    Z_fit = Z_avg[FITSLICE:]
    BA_fit = BA_avg[FITSLICE:]
    BA_e_fit = BA_avg_e[FITSLICE:]


    # perform the fit
    fit_par, fit_par_e, chi2 = makeFit(A_fit, Z_fit, BA_fit, BA_e_fit)

    # make the figure 
    fig, ax = makeFig(A_plot, BA_plot, fit_par, fit_par_e, chi2)

    # plot data
    plotData(ax, A_plot, BA_plot, BA_e_plot)

    #plot fit function
    plotFit(ax, As, Z_avg, fit_par, FITSLICE)

    plt.show()

    return



if __name__ == "__main__":
    main()