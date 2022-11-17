import csv
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem import Descriptors,FilterCatalog,rdqueries,RDConfig
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcFractionCSP3
# from rdkit.Chem.rdchem import GetNumHeavyAtoms
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

import pandas as pd
from tqdm import tqdm

df = pd.DataFrame(columns=["SMILES", "Formula", "Total Molecular Weight","Isomers", "Total Isomers in SMILES", "Number of Atoms","Number of Aromatic Atoms","Fraction csp3","Rotatable Bonds","NumHAcceptors", "NumHDonors","Total Molar Refractivity","tPSA","LogP","Log S (delancy equation)","Hetro Atoms","Number of Rings","Total Carbon","Synthetic Accessibility"])

def logs_delancy_equation(LogP, Molecular_Weight, rotatable_bonds, Aromatic_proportion):
	return 0.16 - (0.63 * LogP) - (0.0062 * Molecular_Weight) + (0.066 * rotatable_bonds) - (0.74 * Aromatic_proportion)

i = 1
with open('temp.smi','r') as csvfile:
	csvreader = csv.reader(csvfile)
	for each_row in csvreader: #tqdm(csvreader):
		chem = each_row[0]
		mol = Chem.MolFromSmiles(chem)
		formula = CalcMolFormula(mol)
		Molecular_Weight = Descriptors.ExactMolWt(mol)
				# mol1 = Chem.AddHs(mol)
				# try:
				# 	embed_check = AllChem.EmbedMolecule(mol1,maxAttempts=5)
				# # print(embed_check)
				# finally:
				# 	if embed_check == -1 :
				# 		embed_check=AllChem.EmbedMolecule(mol1,useRandomCoords=True)
				# print(embed_check)

				# volume = AllChem.ComputeMolVolume(mol1)
		number_of_atoms = mol.GetNumAtoms()
		isomers = tuple(EnumerateStereoisomers(mol))
		total_isomers_in_smiles = ""
		for smi in sorted(Chem.MolToSmiles(x, isomericSmiles=True) for x in isomers):
			total_isomers_in_smiles = smi + " , " + total_isomers_in_smiles

		tPSA = Descriptors.TPSA(mol)
		LogP = Descriptors.MolLogP(mol)
		H_bond_doner = Chem.Lipinski.NumHDonors(mol)
		H_bond_acceptors = Chem.Lipinski.NumHAcceptors(mol)
		aromatic_atoms = len(list(mol.GetAromaticAtoms()))
		heavy_atoms = Chem.rdchem.GetNumHeavyAtoms(mol)
		Aromatic_proportion = aromatic_atoms / heavy_atoms
		molecular_refractivity = Chem.Crippen.MolMR(mol)
		number_of_rings = Descriptors.rdMolDescriptors.CalcNumRings(mol)
		number_of_hetro_atoms = Descriptors.rdMolDescriptors.CalcNumHeteroatoms(mol)
		rotatable_bonds = Chem.Lipinski.NumRotatableBonds(mol)
		carbon = Chem.rdqueries.AtomNumEqualsQueryAtom(6)
		number_of_carbon = len(mol.GetAtomsMatchingQuery(carbon))
		logs_delancy = logs_delancy_equation(LogP, Molecular_Weight, rotatable_bonds, Aromatic_proportion)
		print("logs_delancy",logs_delancy)
		df.loc[i] = [chem, formula,Molecular_Weight, len(isomers),total_isomers_in_smiles, number_of_atoms, aromatic_atoms, CalcFractionCSP3(mol),rotatable_bonds, H_bond_acceptors, H_bond_doner, molecular_refractivity, tPSA, LogP, logs_delancy,number_of_hetro_atoms,number_of_rings,number_of_carbon,sascorer.calculateScore(mol)]
		i = i + 1

df.to_csv('physicochemical_property.csv')
