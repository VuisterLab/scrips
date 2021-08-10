__info__    = "A python script to compare molecular Smiles."
__author__  = "Luca Mureddu"
__date__    = "2020-01-08 16:23:25 +0000 (Wed, Jan 08, 2020)"
__required__= "python 3, numpy, pandas, matplotlib, scipy"

import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
import scipy.stats as ss


def correlateProfiles(array1, array2, engine='correlation'):

    """
    :param array1 and array2:  curves to compare. 1d or nd array
    :param engine: 'correlation' or 'spearmanr'.
    :return: correlation value (float 0 to 1)
    - if nd concatenate the single curves one after the other for array1 and 2 to create two single curves to be correlated.
    - If not of same lenght they'll be interpolated.
    - Apply np.corrcoef or spearmanr
    """
    corrcoef = None
    c1, c2 = array1, array2
    if array1.ndim >1 and array2.ndim >1:
        c1, c2 = np.concatenate(array1), np.concatenate(array2)
    if len(c1) != len(c2):  # we interpolate first to get same xs/ys length
        points = int(np.mean([len(c1),len(c2)]))
        x, x2 =  np.arange(len(c1)), np.arange(len(c2))
        xnew, xnew2 = np.linspace(min(x), max(x), num=points), np.linspace(min(x2), max(x2), num=points)
        f, f2 = interp1d(x, c1), interp1d(x2, c2)
        c1, c2 = f(xnew), f2(xnew2)
    if engine == 'correlation':
        corrcoef = np.corrcoef(c1, c2)[0][1]  # Gives back the correlation matrix for the two arrays
    if engine == 'spearmanr':
        spearmanr = ss.spearmanr(c1, c2)  # Gives the spearman correlation for the two arrays
        corrcoef = spearmanr.correlation
        pvalue = spearmanr.pvalue
    return corrcoef, c1, c2

def _interpolate(curves):
    """
    make multiple curve of same length by linear interpolating the missing points of the shortest curves.
    """
    from scipy.interpolate import interp1d
    ys = []
    counts = [len(i) for i in curves]
    points = int(np.mean(counts))
    for i in curves:
        y = np.arange(len(i))
        ynew = np.linspace(min(y), max(y), num=points)
        f = interp1d(y, i)
        c = f(ynew)
        ys.append(c)
    return ys

def compareSimils(AZDfile, ABTfile, VRSfile):
    '''
    Create figure 4A for the  Frontiers Journal paper 2021 .
    Plot the SIMILARITIES of three case studies.
    params: the excel files containing the smiles and other chemical properties. Located in this directory.
    '''

    AZD =  pd.read_excel(AZDfile)
    ABT = pd.read_excel(ABTfile)
    VRS =  pd.read_excel(VRSfile)

    yAZD = AZD.similarities
    yABT = ABT.similarities
    yVRS = VRS.similarities


    ny = yAZD, yABT, yVRS = _interpolate([yAZD, yABT, yVRS])
    cAZD_ABT = correlateProfiles(yAZD, yABT, engine='correlation')[0]
    cAZD_VRS = correlateProfiles(yAZD, yVRS, engine='correlation')[0]
    cABT_VRS = correlateProfiles(yABT, yVRS, engine='correlation')[0]

    print('cAZD_ABT', cAZD_ABT, )
    print('cAZD_VRS', cAZD_VRS, )
    print('cABT_VRS', cABT_VRS, )

    ax = plt.plot(ny[0], label='AZD')
    ax = plt.plot(ny[1], label='ABT')
    ax = plt.plot(ny[2], label='VRS')
    plt.xlabel('Development Steps (interpolated)')
    plt.ylabel('Tanimoto Coefficient (Similarities)')
    plt.title('Fragments to Drugs similarities')
    plt.legend()
    plt.show()



if __name__ == '__main__':

    localPath = os.getcwd()
    AZDfile = os.path.join(localPath, 'smilesAndOthers_AZD.xlsx')
    ABTfile = os.path.join(localPath, 'smilesAndOthers_Venetoclax.xlsx')
    VRSfile = os.path.join(localPath, 'smilesAndOthers_s64315.xlsx')
    compareSimils(AZDfile, ABTfile, VRSfile)