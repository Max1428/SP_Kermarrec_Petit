# SP_Kermarrec_Petit
Projet rapide réaliser sur une semaine dans le cadre de notre formation de Master en Biologie-informatique à l'Université de Paris.
Ce logiciel non copyrigthé calcul l'orientation de la membrane d'une protéine transmembranaire en se basant sur les propriété hydrophobes du milieu intracellulaire.

##Input :
###Obligatoire : 
Fichier pdb d'une molécule a une seule chaine, stockée dans le répertoire data
###Optionnel :
--N nombre (par défault 20) indiquant le nombre de points générés et répartis de facon homogène sur une sphère. Plus le nombre est élevé plus la précision augmente.
--ASAT nombre (compris entre 0 et 1, par défault 0.3) indiquant le critère de sélection des acides aminés exposés au solvant.
A noté qu'il semble y avoir un bug a partir de 0.6 (leproramme ne smble pas s'arreter.
--MBW nombre (en angstöm, par défault 15) indiquant la taille de la membrane.

##Output :
Un fichier text nomproteine.txt contenant :
	-un rappel des différents arguments utilisés
	-les coordonnées du centre de masse de la proteine (x,y,z), 
	-le score ainsi que les facteur (a,b,c,d1,d2) des coordonées cartésiennes du plan tilt pour résoudre l'équation ax+by+cz+d=0
	-La séquence des acides aminés présant dans la membrane
Un fichier pdb nomproteine_ca.pdb contenant tous les carbonnes alpha de la proteine.
Toutes les fichiers de sortis se trouvent dans le répertoire results

##Installation :
Creation d'un environnement conda avec toutes les dépendances.
conda env create -f projet_Membrane_MK_FP.yml
conda activate projet_MK_FP

##Utilisation :
Le programme se lance a partir du répertoire src.
Pour executer le programme placer vous dans le répertoire conda creer précédement via la commande:
conda activate nom_du_répertoire

Lancer ensuite le logiciel via la commande :
python main.py fichierpdb --option
exemple : python3 main.py ../data/5llu.pdb

