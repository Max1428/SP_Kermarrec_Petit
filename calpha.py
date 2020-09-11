###Function reading a pdb file, extracting all alpha carbon and creating a new output file

import argparse

def trouve_calpha(pdbfile, out):

	if out=='auto':
		output = pdbfile[:-4]+"_ca.pdb" #output file creation
	else :
		output=out

	with open(pdbfile, 'r') as fillin, open(output, 'w') as fillout:
		lines = fillin.readlines()
		for line in lines:
			if not line.startswith("ATOM"):
				continue
			if line[12:16].strip() == "CA":
				fillout.write(line)



if __name__ == "__main__":

	parser = argparse.ArgumentParser()
	parser.add_argument("PDBfile", help = "PDB file from which alpha carbons are extracted", type=str)
	parser.add_argument("--OutFile", help="Select the output file name.", type=str, default='auto')
	args = parser.parse_args()


	if args.PDBfile[-4:] != '.pdb': #Vérification de l'extension pdb
		print("Veuillez recommencer en sélectionnant un fichier avec l'estension '.pdb'")
	else: 
		trouve_calpha(args.PDBfile, args.OutFile)

