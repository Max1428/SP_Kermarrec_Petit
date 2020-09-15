# Considérons deux plans parallèles et non confudus, P1 et P2, d'équation :

# (P1) : (A)*x + (B)*y + (C)*z + D1
# (P2) : (A)*x + (B)*y + (C)*z + D2

# Un point, de coordonnée [Xa,Ya,Za] appartient à la portion d'espace situé entre les deux plans ssi :

#    Cas n°1 : D1 > D2   et   D2 < A*Xa + B*Ya + C*Za < D1
#
#    Cas n°2 : D2 > D1   et   D1 < A*Xa + B*Ya + C*Za < D2

# Cette fonction prend en argument les quatre constantes des plans concernés ainsi que les coordonnées du point dont on veut savoir s'il appartient à l'espace situé entre les deux plans.

import numpy

def entre_deux_plans(plan_1, plan_2, point):
	#print("Dans entre deux plans")
	resolution = float( (plan_1[0]*point[0]) + (plan_1[1]*point[1]) + (plan_1[2]*point[2]) )
	resolution = resolution*-1
	#(resolution, plan_1[3], plan_2[3])
	if resolution >= plan_1[3] and resolution <= plan_2[3] :
		#print("ouiiiii")
		return True


	else :
		#print("false")
		return False

#Fonction déterminant le nombre de carbonealpha avec hydrophobe entre les deux plans
def score_plan(plan_1, plan_2, ca_list, ca_hydro):
	#print("Dans score_plan")
	score=0

	for i in ca_hydro:
		#print(i)
		#print(ca_list[i][1:4])
		if entre_deux_plans(plan_1, plan_2, ca_list[i][1:4].tolist()):
			score += 1
	#print(score)
	return score


#fonction déterminant les planes avec le plus grand nombre de carbone alpha hydrophobes dans l'ensemble des plans possible avec 
def max_score(plan_1, plan_2, ca_list, ca_hydro, maxi):
	#print("Dans max_score")
	score_max=score_plan(plan_1, plan_2, ca_list, ca_hydro)
	tmp_plan_1 = list(plan_1)
	tmp_plan_2 = list(plan_2)

	#print(ca_hydro)
	#print(len(ca_hydro))

	for i in range(1, int(maxi)*2-14):
		tmp_plan_1[3]+=1
		tmp_plan_2[3]+=1
		#print(tmp_plan_1[3])
		#print("dvalue : {}\t{}".format(tmp_plan_1[3], tmp_plan_2[3]))
		tmp_score=score_plan(tmp_plan_1, tmp_plan_2, ca_list, ca_hydro)
		#print(tmp_score)
		#print("tmp_score = {}".format(tmp_score))
		if(tmp_score>score_max):
			score_max=tmp_score
	#print("\n")


	return score_max




# Finalement, cette fonction rapporte une valeur True ou False.
