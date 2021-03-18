import matplotlib.pyplot as plt
import matplotlib.colors as cl
import numpy as np
import uncertainties as unc
#import sys
#sys.path.insert(0, './nuclide-data')
#from nuclide_data import *
from nuclideData.nuclide_data import *


def getData():

    unstableData=[] 
    stableData=[]
    Zmax=119

    # loop over Z
    for Z in range(0,Zmax):

        iso = isotopes[Z]

        # loop over isotopes for current Z
        for j in range(0,len(iso)):
            
            # get N from A of current isotope
            A = iso[j]
            N = iso[j] - Z

            if (nuc(Z,A)['stable']==False):

                # get decay time 
                decTime = nuc(Z, A)['half-life']
                # store data
                unstableData.append((N,Z,decTime))

            # fix "inf" for stable nuclides
            if(nuc(Z,A)['stable']==True):   
                # store data
                stableData.append((N,Z))
    
    return unstableData, stableData


def figureSetup(data):

    # figure and axes
    fig = plt.figure(figsize=(12, 6), dpi = 100)
    ax = fig.add_subplot(1, 1, 1)

    ax.set_title('Nuclides Chart')
    ax.set_xlabel('N')
    ax.set_ylabel('Z')
    ax.set_xlim(0, data[-1][0])
    ax.set_ylim(0, data[-1][1])

    return fig, ax


def main():
    
    unstableData, stableData = getData()

    unstableX=[]
    unstableY=[]
    unstableZ=[]

    stableX=[]
    stableY=[]

    for i in range(len(unstableData)):
        # store N as X
        unstableX.append(unstableData[i][0])
        # store Z as Y
        unstableY.append(unstableData[i][1])
        # store decay time as Z
        unstableZ.append(unstableData[i][2])

    for i in range(len(stableData)):
        # store N as X
        stableX.append(stableData[i][0])
        # store Z as Y
        stableY.append(stableData[i][1])

    # create figure and axes
    fig, ax = figureSetup(unstableData)

    # scatterplot
    unstablePlot = ax.scatter(unstableX, unstableY, c=np.array(unstableZ), cmap='rainbow', norm=cl.LogNorm(), s=5, marker = 's')
    stablePlot = ax.scatter(stableX, stableY, color='black', s=5, marker = 's')

    # add X=Y diagonal dashed line
    ax.plot(np.linspace(0, unstableData[-1][1]), np.linspace(0, unstableData[-1][1]), '--', color='black', alpha = 0.3)

    # add the colorbar on the right
    fig.colorbar(unstablePlot, label = 'Decay time')

    plt.show()


# Call the main function when running the script
if __name__ == "__main__":
    main()