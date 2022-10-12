import csv
import os
# import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import FilterCatalog
from rdkit.Chem import rdqueries
from rdkit.Chem.rdMolDescriptors import CalcMolFormula 
from tqdm import tqdm

# Constant section
total_id = []
total_smiles = []
total_formula = []
total_num_of_hydrogen_doners = []
total_num_of_hydrogen_acceptors = []
total_molecular_weight = []
total_number_of_atoms = []
total_molar_refractivity = []
total_rotatable_bonds = []
total_hetro_atoms = []
total_number_of_rings = []
total_logP = []
total_tPSA = []
total_carbon = []
total_lipinski_check = []
total_ghose_drug_check_prefer = []
total_ghose_drug_check = []
total_egan_drug_likeness_check = []
total_mugge_drug_like_ness = []
total_veber_drug_check = []
total_brenk = []
total_pains = []

# Drug likeness check with different functions

def lipinski_drug_like_ness(H_bond_acceptors,H_bond_doner,Molecular_Weight,LogP):
	if H_bond_acceptors > 10 or H_bond_doner > 5 or Molecular_Weight > 500 or LogP > 5 :
		return "Lipinski rule error"
	else:
		return 0

def ghose_drug_like_ness_prefer(Molecular_Weight,LogP,molecular_refractivity,number_of_atoms):
	if (1.3 > LogP or LogP > 4.1) and (230 < Molecular_Weight or  Molecular_Weight > 390) and (70 > molecular_refractivity or molecular_refractivity > 110) and (30 > number_of_atoms or number_of_atoms > 55):
		return "Ghose Preferred rule error"
	else:
		return 0

def ghose_drug_like_ness(Molecular_Weight,LogP,molecular_refractivity,number_of_atoms):
	if (LogP > 5.6 or LogP < -0.4) or (Molecular_Weight < 160 or Molecular_Weight > 480) or (molecular_refractivity < 40 or molecular_refractivity > 130) or (number_of_atoms < 20 or number_of_atoms > 70):
		return "Ghose Rule error"
	else:
		return 0

def egan_drug_like_ness(tPSA,LogP):
	if tPSA > 131.6 or LogP > 5.88 :
		return "Egan Druglikeness Check Error"
	else:
		return 0
def mugge_drug_like_ness(Molecular_Weight, H_bond_acceptors, H_bond_doner, tPSA, LogP, number_of_rings,number_of_hetro_atoms,number_of_carbon,rotatable_bonds):
	if 200 > Molecular_Weight < 600  or H_bond_acceptors > 10 or H_bond_doner > 5 or tPSA > 150 and -2 > LogP > 5 or number_of_rings > 7 or number_of_hetro_atoms < 2 or number_of_carbon < 5 or rotatable_bonds > 15:
		return "Mugge Druglikeness check error"
	else :
		return 0

def veber_drug_like_ness(tPSA,rotatable_bonds):
	if tPSA > 140 or rotatable_bonds > 10 :
		return "Veber Druglikeness check error"
	else :
		return 0

def brenk(mol1):
	params = FilterCatalog.FilterCatalogParams()
	params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
	catalog = FilterCatalog.FilterCatalog(params)
	if catalog.HasMatch(mol1) :
		return 1
	else:
		return 0

def pains(mol2):
	params = FilterCatalog.FilterCatalogParams()
	params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
	catalog = FilterCatalog.FilterCatalog(params)
	if catalog.HasMatch(mol2):
		return 1
	else:
		return 0

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
		total_lipinski_check.append(lipinski_check)
		ghose_drug_check_prefer = ghose_drug_like_ness_prefer(Molecular_Weight,LogP,molecular_refractivity,number_of_atoms)
		total_ghose_drug_check_prefer.append(ghose_drug_check_prefer)
		ghose_drug_like_ness_check = ghose_drug_like_ness(Molecular_Weight,LogP,molecular_refractivity,number_of_atoms)
		total_ghose_drug_check.append(ghose_drug_like_ness_check)
		egan_drug_check = egan_drug_like_ness(tPSA,LogP)
		total_egan_drug_likeness_check.append(egan_drug_check)
		mugge_drug_check = mugge_drug_like_ness(Molecular_Weight, H_bond_acceptors, H_bond_doner, tPSA, LogP, number_of_rings,number_of_hetro_atoms,number_of_carbon,rotatable_bonds)
		total_mugge_drug_like_ness.append(mugge_drug_check)
		veber_drug_check = veber_drug_like_ness(tPSA,rotatable_bonds)
		total_veber_drug_check.append(veber_drug_check)
		total_smiles.append(chem)
		total_formula.append(formula)
		total_num_of_hydrogen_doners.append(H_bond_doner)
		total_num_of_hydrogen_acceptors.append(H_bond_acceptors)
		total_molecular_weight.append(Molecular_Weight)
		total_number_of_atoms.append(number_of_atoms)
		total_molar_refractivity.append(molecular_refractivity)
		total_rotatable_bonds.append(rotatable_bonds)
		total_hetro_atoms.append(number_of_hetro_atoms)
		total_number_of_rings.append(number_of_rings)
		total_logP.append(LogP)
		total_tPSA.append(tPSA)
		total_carbon.append(number_of_carbon)
		total_brenk.append(brenk(mol))
		total_pains.append(pains(mol))
		total_id.append(i)
		i = i + 1


zip_data = zip (total_id,total_smiles,total_formula,total_molecular_weight,total_number_of_atoms,total_rotatable_bonds,total_num_of_hydrogen_acceptors,total_num_of_hydrogen_doners,total_molar_refractivity,total_tPSA,total_logP,total_hetro_atoms,total_number_of_rings,total_carbon,total_lipinski_check,total_ghose_drug_check_prefer,total_ghose_drug_check,total_egan_drug_likeness_check,total_mugge_drug_like_ness,total_veber_drug_check,total_brenk,total_pains)

with open('output.csv', 'w') as writeFile:
	writer = csv.writer(writeFile, delimiter=",")
	writer.writerow(["ID","SMILES", "Formula", "Total Molecular Weight","Number of Atoms","Rotatable Bonds","NumHAcceptors", "NumHDonors","Total Molar Refractivity","tPSA","LogP","Hetro Atoms","Number of Rings","Total Carbon","Lipinski Check","Ghose Check Prefer","Ghose Drug Likeness check", "Egan Druglikeness", "Mueggie Druglikeness Check","Veber Druglikeness Check","Brenk", "PAINS"])
	writer.writerows(zip_data)
writeFile.close()