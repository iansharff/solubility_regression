import math

import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from tools import helpers as h

def main():
    COLUMN_NAMES = ['iupac', 'log_solubility', 'log_solubility_pred', 'SMILES']
    MOLS_PER_PAGE = 50

    data = pd.read_csv('data/delaney.csv', names=COLUMN_NAMES, header=0)
    molecules = [Chem.MolFromSmiles(smi) for smi in data.SMILES]
    
    for i, mol in enumerate(molecules):

        # Add _Name property to molecule
        mol.SetProp('_Name', data.iupac.iloc[i])
        # Compute 2D coordinates of molecule
        _ = AllChem.Compute2DCoords(mol)


    n_pages = math.ceil(len(molecules) / MOLS_PER_PAGE)
    print(n_pages)

    for i in range(n_pages):
        start = i * MOLS_PER_PAGE
        end = start + MOLS_PER_PAGE
        mols_on_page = molecules[start:end]
        img = Draw.MolsToGridImage(
            mols_on_page,
            molsPerRow=5,
            legends=[mol.GetProp('_Name') for mol in mols_on_page],
        )
        filepath = f"images/molecules/Page_{i:04d}.o.png"
        img.save(filepath)


if __name__ == '__main__': main()