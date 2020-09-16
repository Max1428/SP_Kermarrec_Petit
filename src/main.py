import argparse
import numpy as np
import pandas as pd
import csv
import math
import os
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
	def __init__(self, file, aacDF, N, ASAT, MBW):
		self.MBW = MBW
		self.N = N
		self.file = file
		self.cafile = file[8:-4]+"_ca.pdb" #output file creation
		self.Calpha=np.array([0,0,0,0,0])
		self.AA_name = [] #Amino Acid name
		self.ASAT = ASAT #Accessibility solvent area threshold
		self.ASA = [] #Accessibility solvent area relative between 0 and 1
		self.ca_hydrophobe = [] #Alpha carbon solvent exposed list
		self.best_tilt = [] #Tilt with the higest score
		self.bilan = "" #Bilan of AA inside the membrane

		self.ca_finder(self.file)
		self.center()
		self.calc_ASA(file)
		self.hydrophobicity(aacDF)
		self.small_sphere=sphere(N)
		self.small_sphere.translation(self.centre)
		self.maxi()
		self.scoring() #Function calculatin best score for every points created
		self.output()

	def output(self):
		plan_1 = self.best_tilt[0:4]
		plan_2=self.best_tilt[0:3]
		plan_2.append(self.best_tilt[4])
		id_list = entre_deux_plans.id_atom_plan(plan_1, plan_2, self.Calpha, self.ca_hydrophobe)
		block_list=[]
		n=0
		for x in range(0, len(id_list)):
			if(x >= n and x+1< len(id_list) and id_list[x+1]==id_list[x]+1):
				n=1
				while(1):
					if(x+n <len(id_list) and id_list[x+n] != id_list[x]+n):
						block_list.append("{}({}-{})".format(len(block_list), id_list[x], id_list[x+n-1]))
						n=x+n
						break
					else:
						n+=1
			elif x>=n:
				block_list.append("{}({})".format(len(block_list), id_list[x]))
				
		self.bilan = ",".join(block_list)

		filename=(self.file[8:-4]+".txt")
		#print(filename)
		with open(filename, "w") as fillout:
			fillout.write("PDB file : {}\n\nPoints:{}\tASSAT:{}\tMembrane width:{}\n\nGravity center : {}\n\nBest tilt : score={} a={}, b={}, c={}, d1={}, d2={}\n\nPartie intramembranaire: {}"
			.format(self.file, self.N, self.ASAT, self.MBW, self.centre, self.best_tilt[5], self.best_tilt[0], self.best_tilt[1],
			 self.best_tilt[2], self.best_tilt[3], self.best_tilt[4], self.bilan))
		os.rename(filename, "../results/"+filename)



	def scoring(self):
		score_max = 0
		for i,vector in enumerate(self.small_sphere.points_center):
			plan_1, plan_2 = plan.plan_construction(self.centre, vector.tolist(), self.dmax, self.MBW)
			score, d, dd = entre_deux_plans.max_score(plan_1, plan_2, self.Calpha, self.ca_hydrophobe, self.dmax)
			
			if score>=score_max:
				#print(d, dd)
				score_max=score
				self.best_tilt=list(plan_1[0:3])
				self.best_tilt.append(d)
				self.best_tilt.append(dd)
				self.best_tilt.append(score_max)


	def ca_finder(self, file): #Read pdbfile and extracing in memorry and in file all alpha carbon
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
		os.rename(self.cafile, "../results/"+self.cafile)

	def center(self): #Calculate the center of the proteine
		self.centre=list(np.mean(self.Calpha, axis=0))[1:4]

	def maxi(self):
		self.dmax=0
		for i in self.Calpha:
			d = math.sqrt((i[1]-self.centre[0])**2 +  (i[2]-self.centre[1])**2 + (i[3]-self.centre[2])**2)
			if d > self.dmax:
				self.dmax=float(d)


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
			if self.ASA[i]>self.ASAT:
				continue
			self.Calpha[i][4]=int(aacDF.loc[[str(self.AA_name[i])],["hydrophobic"]].values)
			self.ca_hydrophobe.append(i)


	def __str__(self):
		return("PDB file : {}\n\nPoints:{}\tASSAT:{}\tMembrane width:{}\n\nGravity center : {}\n\nBest tilt : score={} a={}, b={}, c={}, d1={}, d2={}\n\nPartie intramembranaire: {}"
			.format(self.file, self.N, self.ASAT, self.MBW, self.centre, self.best_tilt[5], self.best_tilt[0], self.best_tilt[1],
			 self.best_tilt[2], self.best_tilt[3], self.best_tilt[4], self.bilan))


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("PDBfile", help = "PDB file", type=str)
	parser.add_argument("--AAfile", help = "Mandatory file auto", type=str, default="../data/amino_acid.csv")
	parser.add_argument("--N", help="Number of point used to create the sphere (default = 20point)"
		, type=int, default=20)
	parser.add_argument("--ASAT", help="Accessibility solvent area threshold between 0 and 1 (default=0.3)",
		type=float, default=0.3)
	parser.add_argument("--MBW", help="Membrane width in Angström (default=15)",
		type=int, default=15)
	args = parser.parse_args()

	aacDF=aa.amino_acid_caracteristics(args.AAfile)
	if args.PDBfile[-4:] != '.pdb': #Vérification de l'extension pdb
		print("Please select file with pdb extension")
	else: 
		pro = proteine(args.PDBfile, aacDF, args.N, args.ASAT, args.MBW)
		print(pro)
		



