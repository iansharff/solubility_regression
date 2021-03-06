import pickle

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

def get_feature_array(mols):
    """
    Return an pd.DataFrame of molecule properties given an array (or array-like) of molecule objects

    Parameters
    ----------
    mols: array-like, array of molecule objects
    
    Returns
    ----------
    mol_features: pd.DataFrame of molecule features
    """
    entries = [get_predictors(mol) for mol in  mols]
    mol_features = pd.DataFrame(data=entries, dtype=float)
    return mol_features

    
def get_predictors(mol):
    """Return a dictionary of properties of an RDKit molecule object to be used as model predictors"""
    
    molecular_weight = Descriptors.MolWt(mol)
    oct_water_partition_coefficient = Descriptors.MolLogP(mol)
    num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
    aromatic_proportion = get_aromatic_proportion(mol)

    entry = {
        'MW': molecular_weight,
        'cLogP': oct_water_partition_coefficient,
        'RB': num_rotatable_bonds,
        'AP': aromatic_proportion
    }

    return entry


def get_aromatic_proportion(mol):
    """Return the calculated aromatic proportion of a molecule"""

    are_aromatic = sum([mol.GetAtomWithIdx(i).GetIsAromatic() for i in range(mol.GetNumAtoms())])
    are_heavy = Descriptors.HeavyAtomCount(mol)
    return are_aromatic / are_heavy


def save(obj, filepath):
    """Write an object to a PKL file at `filepath`"""
    with open(filepath, 'wb') as pkl:
        pickle.dump(obj, filepath)
    
def load(filepath):
    """Return a pickled object at `filepath`"""
    pkl = open(filepath, 'rb')
    obj = pickle.load(pkl)
    pkl.close()
    return obj