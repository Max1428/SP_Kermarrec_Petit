# Considérons deux plans respectant l'équation générale suivante : 

#   A*x + B*y + C*z + D = 0

# Ces deux plans sont parallèles et non confondues ssi leurs coefficients directeurs sont, respectivements, égaux ou proportionnels et que leur constante D est différente.

# Par conséquent, connaissant l'équation d'un premier plan (P1), il est possible de construire un plan parallèle et non confondu (P2) au plan (P1) en incrémentant la valeur de la constante D.

# Cette fonction retourne les équations des deux plans parallèles et non confudues, contruits à partir d'un vecteur normal

def plan_construction(centre, vecteur, maxi):

	vect_normal = [vecteur[0]-centre[0],vecteur[1]-centre[1],vecteur[2]-centre[2]]
    
	plan_eq_1 = [vect_normal[0], vect_normal[1], vect_normal[2], float(-1)*((vect_normal[0]*centre[0]) + (vect_normal[1]*centre[1]) + (vect_normal[2]*centre[2]))]
	plan_eq_1[3]=plan_eq_1[3]
	plan_eq_2 = list(plan_eq_1)
	plan_eq_2[3] = plan_eq_2[3] + 15
    #print(plan_eq_1[3], plan_eq_2[3])
    
	return plan_eq_1, plan_eq_2
