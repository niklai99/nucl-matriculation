from nuclideData.nuclide_data import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import numpy as np
import re
from operator import itemgetter
from itertools import groupby
from statistics import mean

massP = 938.27208816 # Mev/c^2
massP_e = 0.00000029 # Mev/c^2
massN = 939.56542052 # Mev/c^2
massN_e = 0.00000054 # Mev/c^2
toMev = 931.493614838934 # ua -> Mev/c^2

Zmax = 119


def getIsotopes():

    data = []

    # loop over Z
    for Z in range(0,Zmax):

        iso = isotopes[Z]

        # loop over isotopes for current Z
        for j in range(0,len(iso)):

            # get N from A of current isotope
            A = iso[j]
            N = iso[j] - Z
           
            # # save only stable isotopes
            # if(nuc(Z,A)['stable']==True):
            #     # store data
            #     data.append((Z, A))

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
    # group binding energy by A
    BAs = [[z for x,y,z,w in g] for k, g in  groupby(par,key=itemgetter(1))]
    # BAs_e = [[w for x,y,z,w in g] for k, g in  groupby(par,key=itemgetter(1))]
    # group weights of binding energy by A
    W = [[w**-2 for x,y,z,w in g] for k, g in  groupby(par,key=itemgetter(1))]

    BA_avg = []
    BA_avg_e = []

    # compute weighted average and standard deviation
    for i in range(len(BAs)):
        w_avg = np.average( BAs[i], weights=W[i] )
        BA_avg.append(w_avg)
        w_avg_e = np.average( (np.array(BAs[i]) - w_avg)**2, weights=W[i] )
        BA_avg_e.append(w_avg_e)

    return As, BA_avg, BA_avg_e


def makePlot(A, BoverA, BoverA_e):

    # figure config
    fig=plt.figure(figsize=(10,6))
    ax = fig.add_subplot(1, 1, 1)

    # set plot range limits 
    XMIN = np.amin(A) - 1
    XMAX = np.amax(A) + 1
    YMIN = np.amin(BoverA) - 1
    YMAX = np.amax(BoverA) + 1

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

    # make the plot
    ax.errorbar(x = A, y = BoverA, yerr = BoverA_e, marker = '.', markersize = 10, elinewidth=1, color = '#000000', linewidth=0, capsize=2, label = 'Data')

    return fig, ax


def main():

    data = getIsotopes()
    # print(data)

    par = computeB(data)
    # print(par)
    
    # Z = [x[0] for x in par]
    # A = [x[1] for x in par]
    # BoverA = [x[2] for x in par]
    # BoverA_e = [x[3] for x in par]
    # print(BoverA_e)

    As, BA_avg, BA_avg_e = computeMean(par)

    fig, ax = makePlot(As, BA_avg, BA_avg_e)

    plt.show()

    return



if __name__ == "__main__":
    main()