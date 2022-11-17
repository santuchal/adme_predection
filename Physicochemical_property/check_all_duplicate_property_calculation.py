import csv
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, QED, rdMolDescriptors
from rdkit.Chem.rdMolDescriptors import CalcAUTOCORR2D, CalcAUTOCORR3D, CalcAsphericity, CalcCrippenDescriptors, CalcExactMolWt
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from rdkit.Chem import Descriptors,FilterCatalog,rdqueries,RDConfig
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcFractionCSP3, CalcHallKierAlpha
from rdkit.Chem.rdMolDescriptors import CalcKappa1, CalcKappa2, CalcKappa3, CalcLabuteASA
# from rdkit.Chem.rdchem import GetNumHeavyAtoms
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

import pandas as pd
from tqdm import tqdm

# df = pd.DataFrame(columns=["SMILES", "Formula", "QED_MW", "QED_ALogP", "QED_HBA", "QED_HBD", "QED_PSA", "QED_ROTB", "QED_AROM", "QED_ALERTS","QED_mol_weight_max","QED_mol_weight_mean","QED_mol_weight_none","QED_default_mol"])#"Total Molecular Weight","Isomers", "Total Isomers in SMILES", "Number of Atoms","Number of Aromatic Atoms","Fraction csp3","Rotatable Bonds","NumHAcceptors", "NumHDonors","Total Molar Refractivity","tPSA","LogP","Log S (delancy equation)","Hetro Atoms","Number of Rings","Total Carbon","Synthetic Accessibility"])

# i = 1
# with open('temp.smi','r') as csvfile:
# 	csvreader = csv.reader(csvfile)
# 	for each_row in tqdm(csvreader):
# 		chem = each_row[0]
# 		mol = Chem.MolFromSmiles(chem)
# 		formula = CalcMolFormula(mol)
# 		QED_MW, QED_ALogP, QED_HBA, QED_HBD, QED_PSA, QED_ROTB, QED_AROM, QED_ALERTS = QED.properties( mol )		
# 		QED_mol_weight_max = QED.weights_max(mol)
# 		QED_mol_weight_mean = QED.weights_mean(mol)
# 		QED_mol_weight_none = QED.weights_none(mol)
# 		QED_default_mol = QED.default(mol)
# 		df.loc[i] = [chem, formula, QED_MW, QED_ALogP, QED_HBA, QED_HBD, QED_PSA, QED_ROTB, QED_AROM, QED_ALERTS,QED_mol_weight_max,QED_mol_weight_mean,QED_mol_weight_none,QED_default_mol ]#Molecular_Weight, len(isomers),total_isomers_in_smiles, number_of_atoms, aromatic_atoms, CalcFractionCSP3(mol),rotatable_bonds, H_bond_acceptors, H_bond_doner, molecular_refractivity, tPSA, LogP, logs_delancy,number_of_hetro_atoms,number_of_rings,number_of_carbon,sascorer.calculateScore(mol)]
# 		i = i + 1

# df.to_csv('./output/all_QED_physicochemical_property.csv')


df = pd.DataFrame(columns=["SMILES", "Formula", "Autocorrelation_descriptors_vector_2d","Autocorrelation_descriptors_vector_3d","CrippenDescriptors","ExactMolWt","FractionCSP3","HallKierAlpha", "CalcKappa1", "CalcKappa2", "CalcKappa3", "LabuteASA"])

i = 1
with open('temp.smi','r') as csvfile:
	csvreader = csv.reader(csvfile)
	for each_row in tqdm(csvreader):
		chem = each_row[0]
		mol = Chem.MolFromSmiles(chem)
		formula = CalcMolFormula(mol)
		Autocorrelation_descriptors_vector_2d = CalcAUTOCORR2D(mol)
		Autocorrelation_descriptors_vector_3d = 0
		# Autocorrelation_descriptors_vector_3d = CalcAUTOCORR3D(mol)	
		CalcCrippenDescriptor = CalcCrippenDescriptors(mol) #2-tuple with the Wildman-Crippen logp,mr values
		ExactMolWt = CalcExactMolWt(mol) # molecule’s exact molecular weight
		FractionCSP3 = CalcFractionCSP3(mol) # fraction of C atoms that are SP3 hybridized
		HallKierAlpha = CalcHallKierAlpha(mol)
		CalcKappa1_1 = CalcKappa1(mol)
		CalcKappa2_1 = CalcKappa2(mol)
		CalcKappa3_1 = CalcKappa3(mol)
		LabuteASA = CalcLabuteASA(mol) # Labute ASA value for a molecule
		print(CalcLabuteASA(mol))
		df.loc[i] = [chem,formula,Autocorrelation_descriptors_vector_2d,Autocorrelation_descriptors_vector_3d,CalcCrippenDescriptor,ExactMolWt,FractionCSP3,HallKierAlpha, CalcKappa1_1, CalcKappa2_1, CalcKappa3_1, LabuteASA]
		i = i + 1

df.to_csv('./output/all_molecular_descriptors_physicochemical_property.csv')