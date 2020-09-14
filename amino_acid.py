import csv
import pandas as pd

def amino_acid_caracteristics():
    with open('amino_acid.csv', newline='\n') as csvfile:
        amino = csv.reader(csvfile, delimiter=',')
        labels = []
        for row in amino:
            head = list(row)
            labels.append(str(list(row)[0]))
            break
        for row in amino:
            labels.append(str(list(row)[0]))
    csvfile.close()
    with open('amino_acid.csv', newline='\n') as csvfile:
        amino = csv.reader(csvfile, delimiter=',')
        df = pd.DataFrame(amino, index=labels, columns=head)
        df.drop(columns = ['amino_acid'], inplace = True)
        df.drop(index = ['amino_acid'], inplace = True)
	return(df)

        #labels = ['amino_acid', 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
        #head = ['amino_acid', 'aromatic', 'polar', 'aliphatic', 'charged', 'negative', 'positive', 'hydrophobic', 'small', 'tiny']

if __name__ == '__main__':
	bidule()
