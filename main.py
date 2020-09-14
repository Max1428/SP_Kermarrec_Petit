import argparse
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP


class proteine(object):
	def __init__(self, file):
		self.file = file
		self.cafile = file[:-4]+"_ca.pdb" #output file creation
		self.Calpha=np.array([0,0,0,0])
		self.AA_name = [] #Amino Acid name
		self.ASA = [] #Accessibility solvent area relative between 0 and 1
		self.hydrophobicityList = []

		self.trouve_calpha(self.file)
			#self.Calpha = np.concatenate((self.Calpha, np.concatenate(self.AA_name)[:,None]), axis=1)
		self.center()
		self.calc_ASA(file)

	def trouve_calpha(self, file):
		with open(file, 'r') as fillin, open(self.cafile, 'w') as fillout:
			lines = fillin.readlines()
			for line in lines:
				if not line.startswith("ATOM"):
					continue
				if line[12:16].strip() == "CA":
					fillout.write(line)
					self.Calpha = np.vstack((self.Calpha, np.float_([line[22:26],line[30:38], line[38:46], line[46:54]]))) #get the number and position of all Calpha
					self.AA_name.append(line[17:20])
			self.Calpha = np.delete(self.Calpha, (0), axis = 0) #Delete the first row of the array fill with 0.
			self.AA_name = np.array([self.AA_name])

	def center(self):
		self.centre=list(np.mean(self.Calpha, axis=0))[1:4]

	def calc_ASA(self, file):
		output = file[:-4]+".dssp"
		p = PDBParser()
		structure = p.get_structure(file[:-4].upper(), file)
		model = structure[0]
		dssp = DSSP(model, file)
		for AA in dssp:
			self.ASA.append(AA[3])


	def __str__(self):
		return("Gravity center : {}\n Accessibility solvent area : {}".format(self.centre, self.ASA))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("PDBfile", help = "PDB file", type=str)
	args = parser.parse_args()



	if args.PDBfile[-4:] != '.pdb': #Vérification de l'extension pdb
		print("Please select file with pdb extension")
	else: 
		pro = proteine(args.PDBfile)
		print(pro)



