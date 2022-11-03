import os
import sys
import csv
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import QED
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors,FilterCatalog,rdqueries,RDConfig
from rdkit.Chem.rdMolDescriptors import CalcMolFormula, CalcFractionCSP3
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer
sys.path.append(os.path.join(RDConfig.RDContribDir, 'NP_Score'))
import npscorer

import pandas as pd

def NP_score(mol, fscore=None):
	if fscore is None:
		fscore =npscorer.readNPModel()
	if mol is None:
		raise ValueError('invalid molecule')
	fp = rdMolDescriptors.GetMorganFingerprint(mol, 2)
	bits = fp.GetNonzeroElements()

# calculating the score
	score = 0.
	for bit in bits:
		score += fscore.get(bit, 0)
	score /= float(mol.GetNumAtoms())

# preventing score explosion for exotic molecules
	if score > 4:
		score = 4. + math.log10(score - 4. + 1.)
	if score < -4:
		score = -4. - math.log10(-4. - score + 1.)
	return score 

def lipinski_drug_like_ness(H_bond_acceptors,H_bond_doner,Molecular_Weight,LogP):
	if H_bond_acceptors > 10 or H_bond_doner > 5 or Molecular_Weight > 500 or LogP > 5 :
		return "Rejected"
	else:
		return "Accepted"
def pfizer_violations_check():
	return "Need to be implement"

def gsk_violations_check():
	return "Need to be implement"

def golden_triangle_check():
	return "Need to be implement"

def ghose_drug_like_ness_prefer(Molecular_Weight,LogP,molecular_refractivity,number_of_atoms):
	if (1.3 > LogP or LogP > 4.1) and (230 < Molecular_Weight or  Molecular_Weight > 390) and (70 > molecular_refractivity or molecular_refractivity > 110) and (30 > number_of_atoms or number_of_atoms > 55):
		return "Rejected"
	else:
		return "Accepted"

def veber_drug_like_ness(tPSA,rotatable_bonds):
	if tPSA > 140 or rotatable_bonds > 10 :
		return "Rejected"
	else :
		return "Accepted"

def muegge_drug_like_ness(Molecular_Weight, H_bond_acceptors, H_bond_doner, tPSA, LogP, number_of_rings,number_of_hetro_atoms,number_of_carbon,rotatable_bonds):
	if 200 > Molecular_Weight < 600  or H_bond_acceptors > 10 or H_bond_doner > 5 or tPSA > 150 and -2 > LogP > 5 or number_of_rings > 7 or number_of_hetro_atoms < 2 or number_of_carbon < 5 or rotatable_bonds > 15:
		return "Rejected"
	else :
		return "Accepted"

def egan_drug_like_ness(tPSA,LogP):
	if tPSA > 131.6 or LogP > 5.88 :
		return 1
	else:
		return 0

params1 = FilterCatalog.FilterCatalogParams()
params1.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
catalog1 = FilterCatalog.FilterCatalog(params1)
def pains(mol2):
	if catalog1.HasMatch(mol2):
		return 1
	else:
		return 0

params = FilterCatalog.FilterCatalogParams()
params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
catalog = FilterCatalog.FilterCatalog(params)
def brenk(mol1):
	if catalog.HasMatch(mol1) :
		return 1
	else:
		return 0


i = 1
with open('../Physicochemical_property/temp.smi','r') as csvfile:
	csvreader = csv.reader(csvfile)
	for each_row in csvreader:
		chem = each_row[0]
		mol = Chem.MolFromSmiles(chem)
		MW, ALogP, HBA, HBD, PSA, ROTB, AROM, ALERTS = QED.properties( mol ) 
		SAscore = sascorer.calculateScore(mol)
		NPScore = NP_score(mol)
		fsp3 = CalcFractionCSP3(mol)
		Lipinski_Violations=lipinski_drug_like_ness(HBA,HBD,MW,ALogP)
		Pfizer_Violations = pfizer_violations_check()
		GSK_Violations = gsk_violations_check()
		Golden_Triangle = golden_triangle_check()
		Ghose_Violations = ghose_drug_like_ness_prefer(MW,ALogP,Chem.Crippen.MolMR(mol),mol.GetNumAtoms())
		Veber_Violations = veber_drug_like_ness(PSA,ROTB)
		carbon = Chem.rdqueries.AtomNumEqualsQueryAtom(6)
		Muegge_Violations = muegge_drug_like_ness(MW,HBA,HBD, PSA, ALogP, Descriptors.rdMolDescriptors.CalcNumRings(mol), Descriptors.rdMolDescriptors.CalcNumHeteroatoms(mol), len(mol.GetAtomsMatchingQuery(carbon)), ROTB )
		Egan_Violations = egan_drug_like_ness(PSA,ALogP)
		PAINS_Alerts = pains(mol)
		BRENK_Alerts = brenk(mol)
