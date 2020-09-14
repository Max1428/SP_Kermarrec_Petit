import argparse
import numpy as np
import pandas as pd
import csv
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

#Our own file
import amino_acid as aa

class sphere(object):
	def __init__(self, samples=20):
		self.points = np.array([0,0,0])
		phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians
		for i in range(samples):
			y = 1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
			radius = math.sqrt(1 - y * y)  # radius at y

			theta = phi * i  # golden angle increment

			x = math.cos(theta) * radius
			z = math.sin(theta) * radius
			self.points = np.vstack((self.points, np.array([x,y,z])))
		self.points = np.delete(self.points, (0), axis = 0)

	def __str__(self):
		return(self.points)

	def translation(self, center):
		centre=np.array(center)
		for i in range(0,len(self.points)):
			self.points_center = np.vstack((self.points_center, self.points[i]+centre))


class proteine(object):
	def __init__(self, file, aacDF):
		self.file = file
		self.cafile = file[:-4]+"_ca.pdb" #output file creation
		self.Calpha=np.array([0,0,0,0,0])
		self.AA_name = [] #Amino Acid name
		self.ASA = [] #Accessibility solvent area relative between 0 and 1
		self.ca_hydrophobe = [] #Alpha carbon solvent exposed list

		self.trouve_calpha(self.file)
			#self.Calpha = np.concatenate((self.Calpha, np.concatenate(self.AA_name)[:,None]), axis=1)
		self.center()
		self.calc_ASA(file)
		#print(str(self.AA_name))
		self.hydrophobicity(aacDF)

	def trouve_calpha(self, file): #Read pdbfile and extracing in memorr and in file all alpha carbon
		with open(file, 'r') as fillin, open(self.cafile, 'w') as fillout:
			lines = fillin.readlines()
			for line in lines:
				if not line.startswith("ATOM"):
					continue
				if line[12:16].strip() == "CA":
					fillout.write(line)
					self.Calpha = np.vstack((self.Calpha, np.float_([line[22:26],line[30:38], line[38:46], line[46:54], 0]))) #get the number and position of all Calpha
					self.AA_name.append(line[17:20].strip())
			self.Calpha = np.delete(self.Calpha, (0), axis = 0) #Delete the first row of the array fill with 0.

	def center(self): #Calculate the center of the proteine
		self.centre=list(np.mean(self.Calpha, axis=0))[1:4]

	def calc_ASA(self, file): #Calculate for each amino acid the relative accessibility solvent area, stocking it in a list
		output = file[:-4]+".dssp"
		p = PDBParser()
		structure = p.get_structure(file[:-4].upper(), file)
		model = structure[0]
		dssp = DSSP(model, file)
		for AA in dssp:
			self.ASA.append(AA[3])

	def hydrophobicity(self, aacDF):
		for i in range(0,len(self.AA_name)):
			if self.ASA[i]<0.30:
				continue
			self.Calpha[i][4]=int(aacDF.loc[[str(self.AA_name[i])],["hydrophobic"]].values)
			self.ca_hydrophobe.append(i)
			#print("AA name : {}\t hydrophobicity: {}".format(self.AA_name[i], aacDF.loc[[self.AA_name[i]],["hydrophobic"]]))
		#print(self.Calpha[:,[4]])

	def __str__(self):
		return("PDB file : {}\nGravity center : {}\n".format(self.file, self.centre))


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("PDBfile", help = "PDB file", type=str)
	parser.add_argument("--Number of point", help="Number of point used to create the sphere (default = 20point)"
		, type=int, default=20)
	args = parser.parse_args()

	aacDF=aa.amino_acid_caracteristics()
	if args.PDBfile[-4:] != '.pdb': #Vérification de l'extension pdb
		print("Please select file with pdb extension")
	else: 
		pro = proteine(args.PDBfile, aacDF)
		print(pro)



