import argparse
import numpy as np
import pandas as pd
import csv
import math
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

#Our own file
import amino_acid as aa
import plan
import entre_deux_plans

class sphere(object):
	def __init__(self, samples=20):
		self.points = np.array([0,0,0])
		self.points_center = np.array([0,0,0])

		phi = math.pi * (3. - math.sqrt(5.))  # golden angle in radians
		for i in range(0,samples):
			y = 1 - (i / float(samples*2 - 1)) * 2  # y goes from 1 to -1
			radius = math.sqrt(1 - y * y)  # radius at y

			theta = phi * i  # golden angle increment

			x = math.cos(theta) * radius
			z = math.sin(theta) * radius
			self.points = np.vstack((self.points, np.array([x,y,z])))
		self.points = np.delete(self.points, (0), axis = 0)

	def __str__(self):
		return("{}".format(self.points))

	def translation(self, center):
		centre=np.array(center)
		for i in range(0,len(self.points)):
			self.points_center = np.vstack((self.points_center, self.points[i]+centre))
		self.points_center = np.delete(self.points_center, (0), axis = 0)


class proteine(object):
	def __init__(self, file, aacDF, N):
		self.file = file
		self.cafile = file[:-4]+"_ca.pdb" #output file creation
		self.Calpha=np.array([0,0,0,0,0])
		self.AA_name = [] #Amino Acid name
		self.ASA = [] #Accessibility solvent area relative between 0 and 1
		self.ca_hydrophobe = [] #Alpha carbon solvent exposed list

		self.ca_finder(self.file)
		self.center()
		self.calc_ASA(file)
		self.hydrophobicity(aacDF)
		self.small_sphere=sphere(N)
		self.small_sphere.translation(self.centre)
		#print(self.small_sphere.points_center)
		self.maxi()
		self.scoring()


	def scoring(self):
		self.plan_score = np.array([0,0,0,0,0,0])
		for i,vector in enumerate(self.small_sphere.points_center):
			plan_1, plan_2 = plan.plan_construction(self.centre, vector.tolist(), self.dmax)
			val = entre_deux_plans.max_score(plan_1, plan_2, self.Calpha, self.ca_hydrophobe, self.dmax)
			self.plan_score = np.vstack((self.plan_score, np.float_([plan_1[0], plan_1[1], plan_1[2], plan_1[3], plan_2[3], val])))
		self.plan_score = np.delete(self.plan_score, (0), axis = 0)
		#print(self.plan_score)


	def ca_finder(self, file): #Read pdbfile and extracing in memorr and in file all alpha carbon
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

	def maxi(self):
		self.dmax=0

		for i in self.Calpha:
			d = math.sqrt((i[1]-self.centre[0])**2 +  (i[2]-self.centre[1])**2 + (i[3]-self.centre[2])**2)
			if d > self.dmax:
				self.dmax=float(d)
		print(self.dmax)

		#maxi=max(self.Calpha.max(axis=0)[1:4].tolist())
		#mini=min(self.Calpha.min(axis=0)[1:4].tolist())
		#print("Maxi ",maxi)
		#print(mini)
		#d_max=[]
		#x=maxi
		#for i in self.centre:
	#		if  (x<0 and i<0): 
		#		d_max.append(abs(x+i))
		#	else:
		#		d_max.append(abs(x-i))
		#x=mini
		#for i in self.centre:
		#	if (x<0 and i<0): 
		#		d_max.append(abs(x+i))
		#	else:
		#		d_max.append(abs(x-i))
		#self.d_maxi=int(max(d_max))
		#print(self.d_maxi)



	def calc_ASA(self, file): #Calculate for each amino acid the relative accessibility solvent area, stocking it in a list
		#output = file[:-4]+".dssp"
		
		p = PDBParser()
		structure = p.get_structure(file[:-4].upper(), file)
		model = structure[0]
		dssp = DSSP(model, file)
		#print(dssp)
		for AA in dssp:
			self.ASA.append(AA[3])



	def hydrophobicity(self, aacDF):
		for i in range(0,len(self.AA_name)):
			if self.ASA[i]>0.30:
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
	parser.add_argument("--N", help="Number of point used to create the sphere (default = 20point)"
		, type=int, default=20)
	args = parser.parse_args()

	aacDF=aa.amino_acid_caracteristics()
	if args.PDBfile[-4:] != '.pdb': #Vérification de l'extension pdb
		print("Please select file with pdb extension")
	else: 
		pro = proteine(args.PDBfile, aacDF, args.N)
		print(pro)
		



