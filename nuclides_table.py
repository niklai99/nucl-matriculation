import sys
sys.path.insert(0, './nuclide-data')
from nuclide_data import *
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import numpy as np
import uncertainties as unc


def getData():

    data=[] 
    Zmax=119

    # loop over Z
    for Z in range(0,Zmax):

        iso = isotopes[Z]

        # loop over isotopes for current Z
        for j in range(0,len(iso)):
            
            # get N from A of current isotope
            A = iso[j]
            N = iso[j] - Z

            # get decay time 
            decTime = nuc(Z, A)['half-life']

            # fix "inf" for stable nuclides
            if(nuc(Z,A)['stable']==True):
                decTime=10e26
                
            # store data
            data.append((N,Z,decTime))
    
    return data


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
    
    data = getData()

    X=[]
    Y=[]
    Z=[]

    for i in range(len(data)):
        # store N as X
        X.append(data[i][0])
        # store Z as Y
        Y.append(data[i][1])
        # store decay time as Z
        Z.append(data[i][2])

    # create figure and axes
    fig, ax = figureSetup(data)

    # scatterplot
    plot = ax.scatter(X, Y, c=np.array(Z), cmap='rainbow', norm=cl.LogNorm(), s=5)

    # add X=Y diagonal dashed line
    ax.plot(np.linspace(0, data[-1][1]), np.linspace(0, data[-1][1]), '--', color='black', alpha = 0.3)

    # add the colorbar on the right
    fig.colorbar(plot, label = 'Decay time')

    plt.show()


# Call the main function when running the script
if __name__ == "__main__":
    main()