import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

RANDOM_SEED = 42

def main():
    delaney = pd.read_csv('data/raw_data/delaney.csv', usecols=['Compound ID', 'SMILES'])
    delaney_mols = [get_molecule(name, smiles) for name, smiles in delaney.itertuples(index=False)]
    write_molecules('./data/outputs/delaney.sdf', delaney_mols)

    # wang = pd.read_excel(
    #     'data/raw_data/wang.xls',
    #     sheet_name=None,
    #     header=2,
    #     names=['Name', 'SLN']
    # )
    # wang_full = pd.concat(wang.values(), ignore_index=True).dropna()



def get_molecule(name, smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp('_Name', name)
    
    mol_H = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_H, randomSeed=RANDOM_SEED)
    return mol_H

def write_molecules(filepath, molecules):
    writer =  Chem.SDWriter(filepath)
    for m in molecules:
        writer.write(m)
    writer.close()


if __name__ == '__main__':
    main()