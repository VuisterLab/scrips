
__info__    = "A script to compute the molecular fingerprint of smiles and calculate the Tanimoto coefficient."
__author__  = "Luca Mureddu"
__date__    = "2020-01-08 16:23:25 +0000 (Wed, Jan 08, 2020)"
__required__= "python 3, openbabel, pandas"


from openbabel import pybel
import os
import pandas as pd


def getSMILESsimilarities(filePath, smilesColumn='smiles', similarityColumn='similarity'):
    """
    filePath: an excel file containing SMILES in a column named smiles
    return a dataframe with an extra column "similarities"

    Compute the  molecular fingerprint from the smiles and calculate the Tanimoto coefficient.
    The  Tanimoto coefficient is used as a value of similarity between the first smiles and all others in the dataframe.
    A value of 1 ==>  same molecules. Lower the value, lower the similarity.

    References
    ===========
    1)  https://en.wikipedia.org/wiki/Jaccard_index
    2)  Bajusz, D., Rácz, A. & Héberger, K. Why is Tanimoto index an appropriate choice for fingerprint-based
        similarity calculations?. J Cheminform 7, 20 (2015). https://doi.org/10.1186/s13321-015-0069-3
    """

    df = pd.read_excel(filePath)
    smiles = df.get(smilesColumn)
    if smiles is None:
        raise RuntimeWarning('No smiles column found in %s.' %filePath)

    referenceSmiles = smiles.values[0]
    refMolPyBel = pybel.readstring("smi", referenceSmiles)

    similarities = []
    for targetSmiles in smiles.values:
        molTarget = pybel.readstring("smi", targetSmiles)
        simil = refMolPyBel.calcfp() | molTarget.calcfp()
        similarities.append(simil)
    df[similarityColumn] = similarities
    return df



if __name__ == '__main__':

    localPath = os.getcwd()
    AZDfile = os.path.join(localPath, 'smilesAndOthers_AZD.xlsx')
    ABTfile = os.path.join(localPath, 'smilesAndOthers_Venetoclax.xlsx')
    VRSfile = os.path.join(localPath, 'smilesAndOthers_s64315.xlsx')

    for filePath in [AZDfile, ABTfile, VRSfile]:
        df = getSMILESsimilarities(filePath)
        df.to_excel(filePath)
    print('Done. Similarity scores updated inside the original files.')
