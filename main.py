import argparse
import numpy as np


class proteine(object):
	def __init__(self, file):
		self.Calpha=np.array([0,0,0,0])
		self.AA_name = []
		with open(file, 'r') as fillin:
			lines=fillin.readlines()
			for line in lines:
				#self.Calpha = np.vstack((self.Calpha, np.float_(line.split()[5:9]))) 
				self.Calpha = np.vstack((self.Calpha, np.float_([line[22:26],line[30:38], line[38:46], line[46:54]]))) #get the number and position of all Calpha
				self.AA_name.append(line[17:20])
			self.Calpha = np.delete(self.Calpha, (0), axis = 0) #Delete the first row of the array fill with 0.
			self.AA_name = np.array([self.AA_name])
			
			#self.Calpha = np.concatenate((self.Calpha, np.concatenate(self.AA_name)[:,None]), axis=1)
		self.center()

	def center(self):
		self.centre=list(np.mean(self.Calpha, axis=0))[1:4]

	def __str__(self):
		return("Gravity center : {}".format(self.centre))

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("PDBfile", help = "PDB file with only alpha carbons extracted, no header, can be done with calpha.py file", type=str)
	args = parser.parse_args()



	if args.PDBfile[-4:] != '.pdb': #Vérification de l'extension pdb
		print("Please select file with pdb extension")
	else: 
		pro = proteine(args.PDBfile)
		print(pro)



