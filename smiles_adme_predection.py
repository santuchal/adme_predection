import csv
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from rdkit import Chem
from rdkit.Chem import Descriptors,FilterCatalog,rdqueries,RDConfig
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcFractionCSP3
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

import pandas as pd
from tqdm import tqdm

# Constant variable pre-loaded otherwise it'll take long time to load for individual molecule

params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
catalog = FilterCatalog.FilterCatalog(params)
params1 = FilterCatalog.FilterCatalogParams()
params1.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
catalog1 = FilterCatalog.FilterCatalog(params1)

# Drug likeness check with different functions

def lipinski_drug_like_ness(H_bond_acceptors,H_bond_doner,Molecular_Weight,LogP):
	if H_bond_acceptors > 10 or H_bond_doner > 5 or Molecular_Weight > 500 or LogP > 5 :
		return 1
	else:
		return 0

def ghose_drug_like_ness_prefer(Molecular_Weight,LogP,molecular_refractivity,number_of_atoms):
	if (1.3 > LogP or LogP > 4.1) and (230 < Molecular_Weight or  Molecular_Weight > 390) and (70 > molecular_refractivity or molecular_refractivity > 110) and (30 > number_of_atoms or number_of_atoms > 55):
		return 1
	else:
		return 0

def ghose_drug_like_ness(Molecular_Weight,LogP,molecular_refractivity,number_of_atoms):
	if (LogP > 5.6 or LogP < -0.4) or (Molecular_Weight < 160 or Molecular_Weight > 480) or (molecular_refractivity < 40 or molecular_refractivity > 130) or (number_of_atoms < 20 or number_of_atoms > 70):
		return 1
	else:
		return 0

def egan_drug_like_ness(tPSA,LogP):
	if tPSA > 131.6 or LogP > 5.88 :
		return 1
	else:
		return 0
def mugge_drug_like_ness(Molecular_Weight, H_bond_acceptors, H_bond_doner, tPSA, LogP, number_of_rings,number_of_hetro_atoms,number_of_carbon,rotatable_bonds):
	if 200 > Molecular_Weight < 600  or H_bond_acceptors > 10 or H_bond_doner > 5 or tPSA > 150 and -2 > LogP > 5 or number_of_rings > 7 or number_of_hetro_atoms < 2 or number_of_carbon < 5 or rotatable_bonds > 15:
		return 1
	else :
		return 0

def veber_drug_like_ness(tPSA,rotatable_bonds):
	if tPSA > 140 or rotatable_bonds > 10 :
		return 1
	else :
		return 0

def brenk(mol1):
	if catalog.HasMatch(mol1) :
		return 1
	else:
		return 0

def pains(mol2):
	if catalog1.HasMatch(mol2):
		return 1
	else:
		return 0

df = pd.DataFrame(columns=["SMILES", "Formula", "Total Molecular Weight","Number of Atoms","Number of Aromatic Atoms","Fraction csp3","Rotatable Bonds","NumHAcceptors", "NumHDonors","Total Molar Refractivity","tPSA","LogP","Hetro Atoms","Number of Rings","Total Carbon","Lipinski Check","Ghose Check Prefer","Ghose Drug Likeness check", "Egan Druglikeness", "Mueggie Druglikeness Check","Veber Druglikeness Check","Brenk", "PAINS","Synthetic Accessibility"])

i = 1
with open('temp.smi','r') as csvfile:
	csvreader = csv.reader(csvfile)
	for each_row in tqdm(csvreader):
		chem = each_row[0]
		mol = Chem.MolFromSmiles(chem)
		formula = CalcMolFormula(mol)
		tPSA = Descriptors.TPSA(mol)
		LogP = Descriptors.MolLogP(mol)
		H_bond_doner = Chem.Lipinski.NumHDonors(mol)
		H_bond_acceptors = Chem.Lipinski.NumHAcceptors(mol)
		Molecular_Weight = Descriptors.ExactMolWt(mol)
		molecular_refractivity = Chem.Crippen.MolMR(mol)
		number_of_atoms = mol.GetNumAtoms()
		number_of_rings = Descriptors.rdMolDescriptors.CalcNumRings(mol)
		number_of_hetro_atoms = Descriptors.rdMolDescriptors.CalcNumHeteroatoms(mol)
		rotatable_bonds = Chem.Lipinski.NumRotatableBonds(mol)
		carbon = Chem.rdqueries.AtomNumEqualsQueryAtom(6)
		number_of_carbon = len(mol.GetAtomsMatchingQuery(carbon))
		lipinski_check = lipinski_drug_like_ness(H_bond_acceptors,H_bond_doner,Molecular_Weight,LogP)
		ghose_drug_check_prefer = ghose_drug_like_ness_prefer(Molecular_Weight,LogP,molecular_refractivity,number_of_atoms)
		ghose_drug_like_ness_check = ghose_drug_like_ness(Molecular_Weight,LogP,molecular_refractivity,number_of_atoms)
		egan_drug_check = egan_drug_like_ness(tPSA,LogP)
		mugge_drug_check = mugge_drug_like_ness(Molecular_Weight, H_bond_acceptors, H_bond_doner, tPSA, LogP, number_of_rings,number_of_hetro_atoms,number_of_carbon,rotatable_bonds)
		veber_drug_check = veber_drug_like_ness(tPSA,rotatable_bonds)
		df.loc[i] = [chem, formula,Molecular_Weight, number_of_atoms, len(list(mol.GetAromaticAtoms())), CalcFractionCSP3(mol), rotatable_bonds, H_bond_acceptors, H_bond_doner, molecular_refractivity, tPSA, LogP, number_of_hetro_atoms,number_of_rings,number_of_carbon,lipinski_check,ghose_drug_check_prefer,ghose_drug_like_ness_check,egan_drug_check,mugge_drug_check,veber_drug_check,brenk(mol),pains(mol),sascorer.calculateScore(mol)]
		i = i + 1

df.to_csv('output1.csv')
